#include "mesh.h"

Mesh::Mesh() {}

// Constructor that initializes the surface with size n x m
Mesh::Mesh(int rows, int cols, double a, bool usingPBC)
    : nodes(rows, cols), a(a), rows(rows), cols(cols),
      usingPBC(usingPBC), currentDeformation()
{
    // Calculate nrElements based on whether usingPBC is true
    if (usingPBC)
    {
        nrElements = 2 * rows * cols;
    }
    else
    {
        nrElements = 2 * (rows - 1) * (cols - 1);
    }
    nrNodes = rows * cols;

    // Now initialize elements with the calculated size
    elements.resize(nrElements);

    m_createNodes();
    m_updateFixedAndFreeNodeIds();
    m_fillNeighbours();

    // we create some elements here, but note that these will need to be replaced
    // if changes are made to the fixed and free status of the nodes.
    createElements();

    // Set ground state energy
    groundStateEnergy = TElement::calculateEnergy(1, 1, 0);
}

Mesh::Mesh(int rows, int cols, bool usingPBC) : Mesh(rows, cols, 1, usingPBC) {}

bool Mesh::isFixedNode(NodeId nodeId)
{
    return (*this)[nodeId]->fixedNode;
}

void Mesh::applyTransformation(Matrix2x2<double> transformation)
{
    // We get all the nodes in the mesh.
    for (Node &node : nodes.data)
    {
        transformInPlace(transformation, node);
    }
    // We also assume we want to transform the current deformation
    applyTransformationToSystemDeformation(transformation);
}

void Mesh::applyTransformationToFixedNodes(Matrix2x2<double> transformation)
{
    // We get the id of each node in the border
    for (NodeId &nodeId : fixedNodeIds)
    {
        transformInPlace(transformation, *(*this)[nodeId]);
    }
}

void Mesh::applyTransformationToSystemDeformation(Matrix2x2<double> transformation)
{
    currentDeformation = currentDeformation * transformation;
}

void Mesh::applyTranslation(Vector2d displacement)
{
    for (Node &node : nodes.data)
    {
        translateInPlace(node, displacement);
    }
}

void Mesh::setInitPos()
{
    for (Node &n : nodes.data)
    {
        n.setInitPos(n.pos());
    }
}

void Mesh::resetForceOnNodes()
{
    for (Node &n : nodes.data)
    {
        n.resetForce();
    }
}

Node Mesh::m_getNeighbourNode(Node node, int direction)
{
    return *(*this)[(*this)[node.id]->neighbours[direction]];
}

double Mesh::averageResolvedShearStress()
{

    double avgRSS = 0;

    for (size_t i = 0; i < elements.size(); i++)
    {
        avgRSS += elements[i].resolvedShearStress;
    }
    return avgRSS / elements.size();
}

int Mesh::nrPlasticEvents() const
{
    // Note, this is effectively the number of plastic events relative to last
    // time this function was called. We rely on the past_m3Nr in the element
    // to be updated in order to find the change since last loading step. If
    // this function is called every 100 loading steps (for example), it will
    // be the number of plasticEvents that have occured in the last 100 steps.
    // (assuming that the mrNr only increases during this period)

    // We also only update this function if the load has changed
    static int nrPlasticEvents = 0;
    static double previousLoad = load;
    if (previousLoad != load)
    {
        nrPlasticEvents = 0;
        for (size_t i = 0; i < elements.size(); i++)
        {
            if (elements[i].plasticEvent())
            {
                nrPlasticEvents += 1;
            }
        }
        previousLoad = load;
    }
    return nrPlasticEvents;
}

// Function to fix the elements of the border vector
void Mesh::fixBorderNodes()
{
    fixNodesInRow(0);
    fixNodesInColumn(0);
    fixNodesInRow(rows - 1);
    fixNodesInColumn(cols - 1);
}

void Mesh::fixNodesInRow(int row)
{
    for (int col = 0; col < cols; ++col)
    {
        nodes[row][col].fixedNode = true;
    }
    m_updateFixedAndFreeNodeIds();
}

void Mesh::fixNodesInColumn(int column)
{
    for (int row = 0; row < rows; ++row)
    {
        nodes[row][column].fixedNode = true;
    }
    m_updateFixedAndFreeNodeIds();
}

void Mesh::m_updateFixedAndFreeNodeIds()
{
    fixedNodeIds.clear();
    freeNodeIds.clear();
    for (int i = 0; i < nodes.data.size(); i++)
    {
        NodeId nodeId(i, cols);
        if (isFixedNode(nodeId))
        {
            fixedNodeIds.push_back(nodeId);
        }
        else
        {
            freeNodeIds.push_back(nodeId);
        }
    }
}

void Mesh::m_createNodes()
{
    int n = rows;
    int m = cols;

    for (int row = 0; row < n; ++row)
    {
        for (int col = 0; col < m; ++col)
        {
            // Set the x and y positions based on the surface indices and spacing "a"
            nodes[row][col] = Node(a, row, col, m);
        }
    }
}

void Mesh::m_fillNeighbours()
{
    int n = rows; // Number of rows
    int m = cols; // Number of columns

    for (int row = 0; row < n; ++row) // Iterate over rows
    {
        for (int col = 0; col < m; ++col) // Iterate over columns
        {
            // Define neighbor indices using periodic boundary conditions
            int left = (col == 0) ? m - 1 : col - 1;
            int right = (col == m - 1) ? 0 : col + 1;
            int up = (row == n - 1) ? 0 : row + 1;
            int down = (row == 0) ? n - 1 : row - 1;

            // Fill in the neighbors
            nodes[row][col].neighbours[LEFT_N] = m_makeNId(row, left);   // left
            nodes[row][col].neighbours[RIGHT_N] = m_makeNId(row, right); // right
            nodes[row][col].neighbours[UP_N] = m_makeNId(up, col);       // up
            nodes[row][col].neighbours[DOWN_N] = m_makeNId(down, col);   // down
        }
    }
}

void Mesh::createElements()
{
    // Note that neighbours must be filled before using this function.

    int rows = this->rows;
    int cols = this->cols;

    // If we are not using PBC, we skipp creating the last elements
    if (!usingPBC)
    {
        rows -= 1;
        cols -= 1;
    }

    // We construct the elements by finding the four nodes that create two
    // opposing cells, but stopping before we get to the last row and columns if
    // we don't use periodic boundary conditions.
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            // We now find the 4 nodes in the current square
            Node n1 = *(*this)[m_makeNId(row, col)];
            // If we are using PBC, we need to use the neighbours to find the
            // adjacent nodes instead of just incrementing col and row.
            // We want to move right and up. (We start the bottom left corner when indexing)
            Node n2 = m_getNeighbourNode(n1, RIGHT_N);
            Node n3 = m_getNeighbourNode(n1, UP_N);
            // n4 is now up AND right of n1
            Node n4 = m_getNeighbourNode(n2, UP_N);

            if (usingPBC)
            {
                if (row == rows - 1 && col == cols - 1)
                {
                    // If we are in the corner, we need to move n2, n3 and n4
                    m_makeGN(n2, n2.id.row, cols);
                    m_makeGN(n3, rows, n3.id.col);
                    m_makeGN(n4, rows, cols);
                }
                else if (col == cols - 1)
                {
                    // If we are in the last column, we need to move n2 and n4
                    m_makeGN(n2, n2.id.row, cols);
                    m_makeGN(n4, n4.id.row, cols);
                }
                else if (row == rows - 1)
                {
                    // If we are in the last row, we need to move n3 and n4
                    m_makeGN(n3, rows, n3.id.col);
                    m_makeGN(n4, rows, n4.id.col);
                }
            }

            // Picture these as top and bottom triangles. The bottom
            // triangle is element 1 with index e1i.
            int e1i = 2 * (row * cols + col); // Triangle 1 index
            int e2i = e1i + 1;
            elements[e1i] = TElement(n1, n2, n3);
            elements[e2i] = TElement(n2, n3, n4);
        }
    }
}

// This is just a function to avoid having to write cols
NodeId Mesh::m_makeNId(int row, int col)
{
    return NodeId(row, col, cols);
}

Vector2d Mesh::makeGhostPos(Vector2d pos, Vector2d shift)
{
    return pos + currentDeformation * shift;
}

void Mesh::m_makeGN(Node &n, int newRow, int newCol)
{
    // Think carefully about where the initial position of the ghost node should
    // be. Is it shifted by the current deformation or not? This might change
    // depending on what you want to simulate. Currently, I think the initial
    // position should not be affected by the current deformation.

    n.isGhostNode = true;
    n.ghostId = NodeId(newRow, newCol, cols + 1);

    // We now need to shift the position.
    // We subtract the old pos from the new pos to get the shift we want to apply
    n.ghostShift = {
        (newCol - n.id.col) * a,
        (newRow - n.id.row) * a,
    };

    n.setInitPos(n.init_pos() + n.ghostShift);
    n.setPos(makeGhostPos(n.pos(), n.ghostShift));
}

void Mesh::printConnectivity(bool realId)
{
    std::string sep;
    std::string end;
    if (nrNodes <= 9)
    {
        sep = "";
        end = " ";
    }
    else
    {
        sep = ",";
        end = "\n";
    }
    for (size_t i = 0; i < nrElements; i++)
    {
        TElement &e = elements[i];
        for (size_t j = 0; j < e.nodes.size(); j++)
        {
            if (realId)
            {
                std::cout << e.nodes[j].id.i << sep;
            }
            else
            {
                std::cout << e.nodes[j].ghostId.i << sep;
            }
        }
        std::cout << end;
    }
    std::cout << '\n';
}

void Mesh::updateElements()
{
#pragma omp parallel
#pragma omp for
    for (size_t i = 0; i < nrElements; i++)
    {
        elements[i].update(*this);
    }
}

void Mesh::applyForceFromElementsToNodes()
{
    for (size_t i = 0; i < nrElements; i++)
    {
        elements[i].applyForcesOnNodes((*this));
    }
}

double Mesh::calculateTotalEnergy()
{
    // we reset the maxEnergy
    maxEnergy = 0;
    // This is the total energy from all the triangles
    double totalEnergy = 0;
    for (size_t i = 0; i < nrElements; i++)
    {
        // We subtract the groundStateEnergy so that the energy is relative to that
        // (so when the system is in it's ground state, the energy is 0)
        totalEnergy += elements[i].energy - groundStateEnergy;
        // We also keep track of the highest energy value.
        if (elements[i].energy > maxEnergy)
        {
            maxEnergy = elements[i].energy;
        }
    }
    averageEnergy = totalEnergy / nrElements;
    return totalEnergy;
}

void transform(const Matrix2x2<double> &matrix, Mesh &mesh, std::vector<NodeId> nodesToTransform)
{
    // We get the adress of each node
    for (NodeId &nodeId : nodesToTransform)
    {
        transformInPlace(matrix, *mesh[nodeId]);
    }
}

std::ostream &operator<<(std::ostream &os, const Mesh &mesh)
{
    for (int i = mesh.cols - 1; i >= 0; --i)
    {
        for (int j = 0; j < mesh.rows; ++j)
        {
            Node n = mesh.nodes[i][j];
            os << n.f << "\t";
        }
        os << "\n";
    }
    return os;
}

void transform(const Matrix2x2<double> &matrix, Mesh &mesh)
{
    // Transform all nodes
    transform(matrix, mesh, mesh.fixedNodeIds);
    transform(matrix, mesh, mesh.freeNodeIds);
}

void translate(Mesh &mesh, std::vector<NodeId> nodesToTranslate, double x, double y)
{
    // We get the adress of each node
    for (NodeId &nodeId : nodesToTranslate)
    {
        translateInPlace(*mesh[nodeId], x, y);
    }
}
void translate(Mesh &mesh, double x, double y)
{
    translate(mesh, mesh.fixedNodeIds, x, y);
    translate(mesh, mesh.freeNodeIds, x, y);
}
