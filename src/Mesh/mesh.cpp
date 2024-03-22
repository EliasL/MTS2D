#include "mesh.h"

Mesh::Mesh() {}

// Constructor that initializes the surface with size n x m
Mesh::Mesh(int rows, int cols, double a, bool usingPBC)
    : nodes(rows, cols), a(a), usingPBC(usingPBC), currentDeformation()
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
    if (!usingPBC)
    {
        fixBorderNodes();
    }
    else
    {
        m_updateFixedAndFreeNodeIds();
    }
    m_fillNeighbours();
    m_createElements();

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

void Mesh::applyTransformationToGhostNodes(Matrix2x2<double> transformation)
{
    // We loop over all the elements, and all the nodes in each element and apply
    // the transformation to the position
}

void Mesh::applyTransformationToSystemDeformation(Matrix2x2<double> transformation)
{
    currentDeformation = currentDeformation * transformation;
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

int Mesh::nrPlasticEvents()
{
    // Note, this is effectively the number of plastic events relative to last
    // time this function was called. We rely on the past_m3Nr in the element
    // to be updated in order to find the change since last loading step. If
    // this function is called every 100 loading steps (for example), it will
    // be the number of plasticEvents that have occured in the last 100 steps.
    // (assuming that the mrNr only increases during this period)
    int nrPlasticEvents = 0;
    for (size_t i = 0; i < elements.size(); i++)
    {
        if (elements[i].plasticEvent())
        {
            nrPlasticEvents += 1;
        }
    }
    return nrPlasticEvents;
}

/*
The idea here is to use the periodic nodes stored in the elements instead of
using the nodes. We reuse all the code we already have for creating fixed
boundary meshes, and only update the positions and forces from the periodic
nodes in the elements. This implementation relies heavily on a specific
ordering of the elements. The elements should be created (and stored in the
elements vector) two and two, such that bottom triangular and top triangular
elements are grouped together. If we then only access the first node of the
bottom elements, we have almost covered all the nodes in the mesh already, but
we miss the nodes of the extra row and column, plus the last node in the
top left corner, since this node is only accessable via a top triangular element.
*/
Mesh Mesh::duplicateAsFixedBoundary()
{
    // We make a new mesh without pbc, but with one extra row and column to make
    // space for the elements that would wrap around in the pbc mesh.
    Mesh newMesh(nodes.rows + 1, nodes.cols + 1, false);

    // In order to

    // We need to get the possitional data from the element instead of the nodes
    // (to make use of the PeriodicNode class), but elements reference the same
    // nodes multiple times. To avoid unneccecary work, we keep track of which
    // nodes have already been copied.
    std::vector<bool> alreadyCopied(newMesh.nrNodes);

    // Iterate over each element in the mesh
    for (size_t i = 0; i < nrElements; ++i)
    {
        TElement &e = elements[i];
        // Iterate over each node in the element
        for (size_t j = 0; j < e.nodes.size(); ++j)
        {
            int newIndex = e.nodes[j].ghostId.i;
            if (!alreadyCopied[newIndex])
            {
                newMesh.nodes.data[newIndex].copyForceAndPos(e.nodes[j]);
                alreadyCopied[newIndex] = true;
            }
        }
    }

    // We also need to update the elements
    for (size_t i = 0; i < nrElements; i++)
    {
        newMesh.elements[i].copyValues(elements[i]);
    }

    // Finally, we need to transfer some other info
    newMesh.a = a;
    newMesh.averageEnergy = averageEnergy;
    newMesh.groundStateEnergy = groundStateEnergy;
    newMesh.load = load;

    return newMesh;
}

// Function to fix the elements of the border vector
void Mesh::fixBorderNodes()
{
    int n = nodes.rows;
    int m = nodes.cols;

    // Loop over the border elements only
    for (int i = 0; i < n; ++i)
    {
        nodes[i][0].fixedNode = true;
        nodes[i][m - 1].fixedNode = true;
    }
    for (int j = 0; j < m; ++j)
    {
        nodes[0][j].fixedNode = true;
        nodes[n - 1][j].fixedNode = true;
    }

    m_updateFixedAndFreeNodeIds();
}

void Mesh::m_updateFixedAndFreeNodeIds()
{
    fixedNodeIds.clear();
    freeNodeIds.clear();
    for (int i = 0; i < nodes.data.size(); i++)
    {
        NodeId nodeId(i, nodes.cols);
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
    int n = nodes.rows;
    int m = nodes.cols;

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
    int n = nodes.rows; // Number of rows
    int m = nodes.cols; // Number of columns

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

void Mesh::m_createElements()
{
    // Note that neighbours must be filled before using this function.

    int rows = nodes.rows;
    int cols = nodes.cols;

    // If we are not using PBC, we skipp creating the last elements
    if (!usingPBC)
    {
        rows -= 1;
        cols -= 1;
    }

    // We construct the elements by finding the four nodes that create two
    // opposing cells, but stopping before we get to the last row and columns if
    // we don't use periodic boundaryconditions.
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
                // std::cout << row << ' ' << col << '\n';
                if (row == rows - 1 && col == cols - 1)
                {
                    // If we are in the corner, we need to move n2, n3 and n4
                    // std::cout << "moving corner\n";
                    m_makeGN(n2, n2.id.row, cols); // Swapped cols and n2.id.row
                    m_makeGN(n3, rows, n3.id.col); // Swapped rows and n3.id.col
                    m_makeGN(n4, rows, cols);      // Swapped rows and cols
                }
                else if (col == cols - 1)
                {
                    // std::cout << "moving first column to the extra column\n";
                    // If we are in the last column, we need to move n2 and n4
                    m_makeGN(n2, n2.id.row, cols); // Swapped cols and n2.id.row
                    m_makeGN(n4, n4.id.row, cols); // Swapped cols and n4.id.row
                }
                else if (row == rows - 1)
                {
                    // std::cout << "moving bottom row up to the extra row\n";
                    // If we are in the last row, we need to move n3 and n4
                    m_makeGN(n3, rows, n3.id.col); // Swapped rows and n3.id.col
                    m_makeGN(n4, rows, n4.id.col); // Swapped rows and n3.id.col
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

// This is just a function to avoid having to write nodes.cols
NodeId Mesh::m_makeNId(int row, int col)
{
    return NodeId(row, col, nodes.cols);
}

VArray Mesh::makeGhostPos(VArray pos, VArray shift)
{
    VArray test = pos + currentDeformation * shift;
    return test;
}

void Mesh::m_makeGN(Node &n, int newRow, int newCol)
{
    n.ghostNode = true;
    n.ghostId = NodeId(newRow, newCol, nodes.cols + 1);

    // We now need to shift the position
    n.ghostShift = {
        (newCol - n.id.col) * a,
        (newRow - n.id.row) * a,
    };

    // We now use the shift to create the new position of the node
    // (Usually, there would be no deformation of the mesh, so this
    // is a bit redundant, but just in case)
    n.setInitPos(makeGhostPos(n.pos(), n.ghostShift));
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
    // #pragma omp parallel
    // #pragma omp for
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
    // This is the total energy from all the triangles
    double totalEnergy = 0;
    for (size_t i = 0; i < nrElements; i++)
    {
        // We subtract the groundStateEnergy so that the energy is relative to that
        // (so when the system is in it's ground state, the energy is 0)
        totalEnergy += elements[i].energy - groundStateEnergy;
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
    for (int i = mesh.nodes.cols - 1; i >= 0; --i)
    {
        for (int j = 0; j < mesh.nodes.rows; ++j)
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
