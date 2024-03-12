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
        nrElements = 2 * (rows - 1) * (cols - 1);
    }

    // Now initialize elements with the calculated size
    elements = std::vector<TElement>(nrElements);

    if (!usingPBC)
    {
        fixBorderNodes();
    }
    m_updateFixedAndFreeNodeIds();
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
NodeId Mesh::m_makeNId(int row, int col)
{
    return NodeId(row, col, nodes.cols);
}
// Make a periodic node
PeriodicNode Mesh::m_makePN(NodeId id, double x, double y)
{
    return PeriodicNode((*this)[id], x, y, nodes.rows, nodes.cols);
}

NodeId Mesh::m_getNeighbourNodeId(NodeId nodeId, int direction)
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

/*
The idea here is to use the periodic nodes stored in the elements instead of
using the nodes. We reuse all the code we already have for creating fixed
boundary meshes, and only update the positions and forces from the periodic
nodes in the elements. This implementation relies heavily on a spesific
ordering of the elements. The elements should be created (and stored in the
elements vector) two and two, such that bottom triangular and top triangular
elements are grouped together. If we then only access the first node of the
bottom elements, we have almost covered all the nodes in the mesh already, but
we miss the nodes of the extra row and column, pluss the last node in the
top left corner, since this node is only accessable via a top triangular element.
*/
Mesh Mesh::duplicateAsFixedBoundary() const
{
    Mesh newMesh(nodes.rows + 1, nodes.cols + 1, false);
    // we only loop through the lower triangles
    for (size_t i = 0; i < nrElements; i += 2)
    {
        // Current element e
        TElement e = elements[i];
        // We update our new mesh with the data from the old mesh
        int nodeIndex = e.nodes[0].periodicId.i;
        newMesh.nodes.data[nodeIndex].update(e.nodes[0]);
        // If on the edge, also update nodes 1 or 2
        if (e.row == nodes.rows - 1)
        {
            int nodeIndex = e.nodes[2].periodicId.i;
            newMesh.nodes.data[nodeIndex].update(e.nodes[2]);
        }
        if (e.col == nodes.cols - 1)
        {
            int nodeIndex = e.nodes[1].periodicId.i;
            newMesh.nodes.data[nodeIndex].update(e.nodes[1]);
        }
    }
    // We have now dealt with all the normal nodes, and all the periodic
    // nodes except for the last node in the top right corner. This node
    // can be accessed by choosing the very last element (assuming that
    // the top element is added last, which in time of writing, it is)
    TElement e = elements[nrElements - 1];
    if (e.row == nodes.rows - 1 && e.col == nodes.cols - 1)
    {
        int nodeIndex = e.nodes[2].periodicId.i;
        newMesh.nodes.data[nodeIndex].update(e.nodes[2]);
    }
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
}

void Mesh::m_updateFixedAndFreeNodeIds()
{
    fixedNodeIds.clear();
    freeNodeIds.clear();
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
            nodes[row][col].id = m_makeNId(row, col);
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
            NodeId n1 = m_makeNId(row, col);
            // If we are using PBC, we need to use the neighbours to find the
            // adjacent nodes instead of just incrementing col and row.
            // We want to move right and down. (We start the top left corner when indexing)
            NodeId n2 = m_getNeighbourNodeId(n1, RIGHT_N);
            NodeId n3 = m_getNeighbourNodeId(n1, UP_N);
            // n4 is now down AND right of n1
            NodeId n4 = m_getNeighbourNodeId(n2, UP_N);

            /*
            This part is a bit complicated.
            If we are using periodic boundary conditions, and we are at the
            boundary of our system, we want to move our nodes into a new
            position. Here is an example: Let's say we have a 3x2 system:

            3 4 5
            0 1 2

            node 2 will be connected to node 0 and 5 (twice). To avoid various
            problems, we add a new column and row, which changes the indexes.

            8  9  10 11
            4  5  6  7
            0  1  2  3

            The following code ensures that we know the original and periodic
            index of the nodes in the element, while also adjusting the position,
            ie. moving node 7 so that it mirrors node 4.
            */
            PeriodicNode pn1 = m_makePN(n1, 0, 0);
            PeriodicNode pn2 = m_makePN(n2, 0, 0);
            PeriodicNode pn3 = m_makePN(n3, 0, 0);
            PeriodicNode pn4 = m_makePN(n4, 0, 0);

            if (usingPBC)
            {
                if (row == rows - 1 && col == cols - 1)
                {
                    // If we are in the corner, we need to move n2, n3 and n4
                    pn2 = m_makePN(n2, a * cols, 0);
                    pn3 = m_makePN(n3, 0, a * rows);
                    pn4 = m_makePN(n4, a * cols, a * rows);
                }
                else if (col == cols - 1)
                {
                    // If we are in the last column, we need to move n2 and n4
                    pn2 = m_makePN(n2, a * cols, 0);
                    pn4 = m_makePN(n4, a * cols, 0);
                }
                else if (row == rows - 1)
                {
                    // If we are in the last column, we need to move n3 and n4
                    pn3 = m_makePN(n3, 0, a * rows);
                    pn4 = m_makePN(n4, 0, a * rows);
                }

                // And now we also apply a loading transformation and adjust
                // the length between the nodes across the border.
                // TODO
            }

            // Picture these as top and bottom triangles. The bottom
            // triangle is element 1 with index e1i.
            int e1i = 2 * (row * cols + col); // Triangle 1 index
            int e2i = e1i + 1;
            elements[e1i] = TElement{pn1, pn2, pn3, row, col};
            elements[e2i] = TElement{pn2, pn3, pn4, row, col};
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
            os << "(" << n.x() << ", " << n.y() << ")\t";
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
