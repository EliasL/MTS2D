#include "mesh.h"

Mesh::Mesh() {}

// Constructor that initializes the surface with size n x m
Mesh::Mesh(int rows, int cols, double a, bool usingPBC)
    : nodes(rows, cols), a(a), usingPBC(usingPBC)
{

    // Calculate nrElements based on whether usingPBC is true
    if (usingPBC)
    {
        nrElements = 2 * rows * cols;
    }
    else
    {
        // Ensure rows and cols are greater than 1 to avoid underflow
        if (rows > 1 && cols > 1)
        {
            nrElements = 2 * (rows - 1) * (cols - 1);
        }
        else
        {
            nrElements = 0; // Handle the edge case where mesh is too small
        }
    }
    // Now initialize elements with the calculated size
    elements = std::vector<TElement>(nrElements);

    if (!usingPBC)
    {
        m_fixBorderNodes();
    }
    m_fillFixedAndFreeNodeIds();
    m_setNodePositions();
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
}

void Mesh::applyTransformationToFixedNodes(Matrix2x2<double> transformation)
{
    // We get the id of each node in the border
    for (NodeId &nodeId : fixedNodeIds)
    {
        transformInPlace(transformation, *(*this)[nodeId]);
    }
}

void Mesh::resetForceOnNodes()
{
    for (Node &n : nodes.data)
    {
        n.resetForce();
    }
}

// This is just a function to avoid having to write nodes.cols
NodeId Mesh::makeNId(int row, int col)
{
    return NodeId(row, col, nodes.cols);
}

NodeId Mesh::getNeighbourNodeId(NodeId nodeId, int direction)
{
    return (*this)[nodeId]->neighbours[direction];
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
        TElement &e = elements[i];
        if (elements[i].plasticEvent())
        {
            nrPlasticEvents += 1;
        }
    }
    return nrPlasticEvents;
}

// Function to fix the elements of the border vector
void Mesh::m_fixBorderNodes()
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
}

void Mesh::m_fillFixedAndFreeNodeIds()
{
    for (int i = 0; i < nodes.data.size(); i++)
    {
        NodeId nodeId = NodeId{i, nodes.cols};
        // If nodeId is a border node,
        if (isFixedNode(nodeId))
        {
            // we add it to the vector.
            fixedNodeIds.push_back(NodeId(i, nodes.cols));
        }
        else
        {
            freeNodeIds.push_back(NodeId(i, nodes.cols));
        }
    }
}

void Mesh::m_setNodePositions()
{
    int n = nodes.rows;
    int m = nodes.cols;

    for (int row = 0; row < n; ++row)
    {
        for (int col = 0; col < m; ++col)
        {
            // Set the x and y positions based on the surface indices and spacing "a"
            nodes[row][col].setPos(col * a, row * a);
            nodes[row][col].setInitPos(col * a, row * a);
            nodes[row][col].id = makeNId(row, col);
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
            nodes[row][col].neighbours[LEFT_N] = makeNId(row, left);   // left
            nodes[row][col].neighbours[RIGHT_N] = makeNId(row, right); // right
            nodes[row][col].neighbours[UP_N] = makeNId(up, col);       // up
            nodes[row][col].neighbours[DOWN_N] = makeNId(down, col);   // down
        }
    }
}

void Mesh::m_createElements()
{
    // Note that neighbours must be filled before using this function.

    int n = nodes.rows;
    int m = nodes.cols;

    // If we are not using PBC, we skipp creating the last elements
    if (!usingPBC)
    {
        n -= 1;
        m -= 1;
    }

    // We construct the elements by finding the four nodes that create two
    // opposing cells, but stopping before we get to the last row and columns if
    // we don't use periodic boundaryconditions.
    for (int row = 0; row < n; ++row)
    {
        for (int col = 0; col < m; ++col)
        {
            // We now find the 4 nodes in the current square
            NodeId n1 = makeNId(row, col);
            // If we are using PBC, we need to use the neighbours to find the
            // adjacent nodes instead of just incrementing col and row.
            // We want to move right and down. (We start the top left corner when indexing)
            NodeId n2 = getNeighbourNodeId(n1, RIGHT_N);
            NodeId n3 = getNeighbourNodeId(n1, UP_N);
            // n4 is now down AND right of n1
            NodeId n4 = getNeighbourNodeId(n2, UP_N);

            // Picture these as top and bottom right angle triangles. The bottom
            // triangle is element 1 with index e1i.
            int e1i = 2 * (row * m + col); // Triangle 1 index
            int e2i = e1i + 1;

            elements[e1i] = TElement{(*this)[n1], (*this)[n2], (*this)[n3]};
            elements[e2i] = TElement{(*this)[n2], (*this)[n3], (*this)[n4]};

            /*
                Helping visualization

                2       6   7   8
                1       3   4   5
                0       0   1   2

                n/m     0   1   2
            */

            // Now we deal with periodic boundary conditions. We can set offsets
            // to the positions in the triangular elements.
            // We only adjust the extra elements that are created on the last
            // row and column TODO THIS IS WRONG
            if (row == n - 1 && usingPBC)
            {
                elements[e1i].moveN[2] = true;
                elements[e1i].yOffset = a * nodes.rows;

                elements[e2i].moveN[0] = true;
                elements[e2i].yOffset = -1 * a * nodes.rows;
            }
            if (col == m - 1 && usingPBC)
            {
                elements[e1i].moveN[1] = true;
                elements[e1i].xOffset = a * nodes.cols;

                elements[e2i].moveN[1] = true;
                elements[e2i].xOffset = -1 * a * nodes.cols;
            }
        }
    }
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
        for (int j = 0; j < 2; ++j)
        {
            Node n = mesh.nodes[i][j];
            os << "(" << n.X() << ", " << n.Y() << ")\t";
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
