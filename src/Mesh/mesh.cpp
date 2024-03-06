#include "mesh.h"

Mesh::Mesh() {}

// Constructor that initializes the surface with size n x m
Mesh::Mesh(int rows, int cols, double a) : nodes(rows, cols), a(a),
                                           nrElements(2 * (rows - 1) * (cols - 1)),
                                           elements(2 * (rows - 1) * (cols - 1))
{
    m_fixBorderNodes();
    m_fillFixedAndFreeNodeIds();
    m_setNodePositions();
    m_fillNeighbours();
    m_createElements();

    // Set ground state energy
    groundStateEnergy = TElement::calculateEnergy(1, 1, 0);
}

Mesh::Mesh(int rows, int cols) : Mesh(rows, cols, 1) {}

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
NodeId Mesh::getNodeId(int row, int col)
{
    return NodeId(row, col, nodes.cols);
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
            nodes[row][col].id = getNodeId(row, col);
        }
    }
}

// Function to fill neighbours using periodic boundary conditions
void Mesh::m_fillNeighbours()
{
    int n = nodes.rows;
    int m = nodes.cols;

    for (int row = 0; row < m; ++row)
    {
        for (int col = 0; col < n; ++col)
        {
            // Define neighbor indices using periodic boundary conditions
            int left = (col == 0) ? m - 1 : col - 1;
            int right = (col == m - 1) ? 0 : col + 1;
            int up = (row == 0) ? n - 1 : row - 1;
            int down = (row == n - 1) ? 0 : row + 1;

            // Fill in the neighbors
            nodes[row][col].neighbours[0] = getNodeId(left, row);  // left
            nodes[row][col].neighbours[1] = getNodeId(right, row); // right
            nodes[row][col].neighbours[2] = getNodeId(col, up);    // up
            nodes[row][col].neighbours[3] = getNodeId(col, down);  // down
        }
    }
}

void Mesh::m_createElements()
{
    int n = nodes.rows;
    int m = nodes.cols;

    for (int row = 0; row < m - 1; ++row)
    {
        for (int col = 0; col < n - 1; ++col)
        {
            // We now find the 4 nodes in the current square
            NodeId n1 = getNodeId(col, row);
            NodeId n2 = getNodeId(col, row + 1);
            NodeId n3 = getNodeId(col + 1, row);
            NodeId a4 = getNodeId(col + 1, row + 1);

            int e1i = 2 * (col * (m - 1) + row); // Triangle 1 index
            int e2i = e1i + 1;

            elements[e1i] = TElement{(*this)[n1], (*this)[n2], (*this)[n3]};
            elements[e2i] = TElement{(*this)[n2], (*this)[n3], (*this)[a4]};
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
