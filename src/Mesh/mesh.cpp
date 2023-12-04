#include "mesh.h"





// BOUNDARY_CONDITIONS -------

BoundaryConditions::BoundaryConditions(double load, double theta, BCF bc)
{
    this->load = load;
    this->theta = theta;
    m_calculateGradiant(bc);
}

void BoundaryConditions::m_calculateGradiant(BCF bc)
{
    switch (bc)
    {
    case BCF::macroShear:
        m_macroShear();
        break;

    default:
        throw std::invalid_argument("Invalid boundary condition function");
        break;
    }
}
/*
 */
void BoundaryConditions::m_macroShear()
{
    double perturb = 0;

    F[0][0] = (1. - load * cos(theta + perturb) * sin(theta + perturb));
    F[1][1] = (1. + load * cos(theta - perturb) * sin(theta - perturb));
    F[0][1] = load * pow(cos(theta), 2.);
    F[1][0] = -load * pow(sin(theta - perturb), 2.);
}

// GRID ----
Mesh::Mesh() {}

// Constructor that initializes the surface with size n x m
Mesh::Mesh(int rows, int cols, double a) : nodes(rows, cols), a(a),
                                           nrTriangles(2 * (rows - 1) * (cols - 1)),
                                           triangles(2 * (rows - 1) * (cols - 1)),
                                           cells(2 * (rows - 1) * (cols - 1))
{
    m_setBorderElements();
    m_fillNonBorderNodeIds();
    m_setNodePositions();
    m_fillNeighbours();
    m_createTriangles();
    // TODO We often loop over either borderNodes or nonBorderNodes, therefore
    // it might be faster if we sort the nodes vector such that the nodes come
    // in order in memory when itterating over them. This would probably save a
    // neglidgable amount of time...
}

Mesh::Mesh(int rows, int cols) : Mesh(rows, cols, 1) {}

bool Mesh::isBorder(NodeId nodeId)
{
    return (*this)[nodeId]->borderNode;
}

void Mesh::applyBoundaryConditions(BoundaryConditions bc)
{
    // We get the id of each node in the border
    for (NodeId nodeId : borderNodeIds)
    {
        transformInPlace(bc.F, *(*this)[nodeId]);
    }
}

void Mesh::resetForceOnNodes()
{
    for (Node n : nodes.data)
    {
        n.f_x = n.f_y = 0;
    }
}

// This is just a function to avoid having to write nodes.cols
NodeId Mesh::getNodeId(int row, int col)
{
    return NodeId(row, col, nodes.cols);
}

// Function to set border elements of the border vector to true
void Mesh::m_setBorderElements()
{
    int n = nodes.rows;
    int m = nodes.cols;

    // Loop over the border elements only
    for (int i = 0; i < n; ++i)
    {
        nodes[i][0].borderNode = true;
        nodes[i][m - 1].borderNode = true;
    }
    for (int j = 0; j < m; ++j)
    {
        nodes[0][j].borderNode = true;
        nodes[n - 1][j].borderNode = true;
    }
}

void Mesh::m_fillNonBorderNodeIds()
{
    for (int i = 0; i < nodes.data.size(); i++)
    {
        NodeId nodeId = NodeId{i, nodes.cols};
        // If nodeId is a border node,
        // LOG(DEBUG) << nodeId << " is border? " << isBorder(nodeId) << std::endl;
        if (isBorder(nodeId))
        {
            // we add it to the vector.
            borderNodeIds.push_back(NodeId(i, nodes.cols));
        }
        else
        {
            nonBorderNodeIds.push_back(NodeId(i, nodes.cols));
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
            nodes[row][col].x = col * a;
            nodes[row][col].y = row * a;
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

// See the bottom of the doc for explination
void Mesh::m_createTriangles()
{
    int n = nodes.rows;
    int m = nodes.cols;

    // CHECK THAT CALCULATIONS CAN BE DONE ON UPSIDE DOWN TRIANGLES WITHOUT MINUS
    for (int row = 0; row < m - 1; ++row)
    {
        for (int col = 0; col < n - 1; ++col)
        {
            // We now find the 4 nodes in the current square
            NodeId a1 = getNodeId(col, row);
            NodeId a2 = getNodeId(col, row + 1);
            NodeId a3 = getNodeId(col + 1, row);
            NodeId a4 = getNodeId(col + 1, row + 1);

            int t1i = 2 * (col * (m - 1) + row); // Triangle 1 index
            int t2i = t1i + 1;

            triangles[t1i] = Triangle{(*this)[a1], (*this)[a2], (*this)[a3]};
            triangles[t2i] = Triangle{(*this)[a2], (*this)[a3], (*this)[a4]};
        }
    }
}

void transform(const Matrix2x2<double> &matrix, Mesh &mesh, std::vector<NodeId> nodesToTransform)
{
    // We get the adress of each node
    for (NodeId nodeId : nodesToTransform)
    {
        transformInPlace(matrix, *mesh[nodeId]);
    }
}

void transform(const Matrix2x2<double> &matrix, Mesh &mesh)
{
    // Transform all nodes
    transform(matrix, mesh, mesh.borderNodeIds);
    transform(matrix, mesh, mesh.nonBorderNodeIds);
}

void translate(Mesh &mesh, std::vector<NodeId> nodesToTranslate, double x, double y)
{
    // We get the adress of each node
    for (NodeId nodeId : nodesToTranslate)
    {
        translateInPlace(*mesh[nodeId], x, y);
    }
}
void translate(Mesh &mesh, double x, double y)
{
    translate(mesh, mesh.borderNodeIds, x, y);
    translate(mesh, mesh.nonBorderNodeIds, x, y);
}