#include "mesh.h"

// NODE AND NODE_ID ------

NodeId::NodeId() : row(0), col(0), i(0) {}
NodeId::NodeId(int row_, int col_, int cols) : row(row_), col(col_), i(row_ * cols + col_) {}
NodeId::NodeId(int i_, int cols) : row(i_ / cols), col(i_ % cols), i(i_) {}

std::ostream &operator<<(std::ostream &os, const NodeId &nodeId)
{
    // This implementation is confusing because (3,4) resembles vector notation
    // where x=3 and y=4.
    // os << "Node " << nodeId.i << "(" << nodeId.col << ", " << nodeId.row << ")";

    // This implementation, while less compact, is clearer.
    os << "Node " << nodeId.i << ", row: " << nodeId.row << ", col: " << nodeId.col;
    return os;
}

Node::Node(double x_, double y_)
{
    x = x_;
    y = y_;
    f_x = f_y = 0;
    borderNode = false;
}

Node::Node() : Node(0, 0) {}

void transformInPlace(const Matrix2x2<double> &matrix, Node &n)
{
    double newX = matrix[0][0] * n.x + matrix[0][1] * n.y;
    double newY = matrix[1][0] * n.x + matrix[1][1] * n.y;
    n.x = newX;
    n.y = newY;
}

Node transform(const Matrix2x2<double> &matrix, const Node &n)
{
    Node result = n;
    transformInPlace(matrix, result);
    return result;
}

void translateInPlace(Node &n, double x, double y)
{
    n.x += x;
    n.y += y;
}

void translateInPlace(Node &n, const Node &delta, double multiplier)
{
    translateInPlace(n, multiplier * delta.x, multiplier * delta.y);
}

Node translate(const Node &n, const Node &delta, double multiplier)
{
    Node result = n;
    translateInPlace(result, delta, multiplier);
    return result;
}

// TRIANGLE -------

// e1 and e2 form basis vectors for the triangle
std::array<double, 2> Triangle::e1() const
{ // a is the lattice spacing of the gird
    return {
        a2->x - a1->x,
        a2->y - a1->y,
    };
}

std::array<double, 2> Triangle::e2() const
{ // a is the lattice spacing of the gird
    return {
        a3->x - a1->x,
        a3->y - a1->y,
    };
}

// Provices a metric tensor for the triangle
Matrix2x2<double> Triangle::metric(MetricFunction f) const
{
    // Symetric matricies would be faster, but only slightly for 2x2 matrix
    Matrix2x2<double> m;
    auto e1_ = e1();
    auto e2_ = e2();
    // There are many ways to calculate a metric. The user can specify which
    // to use.
    switch (f)
    {
    case MetricFunction::faicella:
        m[0][0] = e1_[0] * e1_[0] + e1_[1] * e1_[1];
        m[1][1] = e2_[0] * e2_[0] + e2_[1] * e2_[1];
        m[1][0] = m[0][1] = e1_[0] * e2_[0] + e1_[1] * e2_[1];
        return m;
    case MetricFunction::epsilon_lineaire: // This function looks super strange
        m[0][0] = e1_[0] - 1;
        m[1][1] = e2_[1] - 1;
        m[1][0] = m[0][1] = e2_[0];
    default:
        throw std::invalid_argument("Invalid metric function");
        break;
    }
}

std::ostream &operator<<(std::ostream &os, const Triangle &triangle)
{
    os << "a1: (" << triangle.a1->x << ", " << triangle.a1->y << "), "
       << "a2: (" << triangle.a2->x << ", " << triangle.a2->y << "), "
       << "a3: (" << triangle.a3->x << ", " << triangle.a3->y << ")";
    return os;
}

std::ostream &operator<<(std::ostream &os, const Triangle *trianglePtr)
{
    if (trianglePtr)
    {
        // Use the existing Triangle operator<< overload
        os << *trianglePtr;
    }
    else
    {
        // Handle the nullptr case
        os << "nullptr";
    }
    return os;
}

// CELL --------------

Cell::Cell(const Triangle &triangle)
{
    // Calculates D
    m_getDeformationGradiant(triangle);

    // Calculates C
    C = triangle.metric(METRICFUNCTION);

    // Calculate C_ and m
    m_lagrangeReduction();
};

Cell::Cell()
{
}

// A basis vector for the cell
double Cell::e1(int index)
{
    return F[0][index];
}

// A basis vector for the cell
double Cell::e2(int index)
{
    return F[1][index];
}

// Calculate Piola stress tensor and force on each node from current cell
// Note that each node is part of multiple cells. Therefore, the force must
// be reset after each itteration.
void Cell::setForcesOnNodes(Triangle &triangle)
{

    if (!hasComputedReducedStress)
        throw std::runtime_error("Reduced stress has not yet been calculated");

    // extended stress is not quite the "real" stress, but it is a component
    // in calculating the piola stress, which is the real stress on the cell,
    // and we can then find the force on each individual node.
    // The name extended_stress does not have much meaning.
    // TODO consider storing this variable in the cell, such that it does
    // not have to be allocated every time the function is called.
    Matrix2x2<double> extended_stress = r_s.sym_orth_conjugate(m);
    // TODO THIS IS WRONG
    P[0][0] = 2 * extended_stress[0][0] * F[0][0] + extended_stress[0][1] * F[1][0];
    P[1][0] = 2 * extended_stress[0][0] * F[0][1] + extended_stress[0][1] * F[1][1];
    P[0][1] = 2 * extended_stress[1][1] * F[1][0] + extended_stress[0][1] * F[0][0];
    P[1][1] = 2 * extended_stress[1][1] * F[1][1] + extended_stress[0][1] * F[0][1];
    // The assignment here is dependant on the shape of the cell.
    // For triangular shapes, the forces on the nodes is applied as shown
    // below. For a general shape, see Gael-notes page 2, partial N^i / partial x_j
    // on how to calculate.

    // DO GENERAL DISTRIBUTION of Piola stress

    // FORCES SHOULD BE RESET AFTER EACH ITERATION
    // MUST SUM (+=) BECAUSE THEY ARE THE SAME NODES THAT ARE IN A FLIPPED TRIANGLE
    triangle.a1->f_x += -P[0][0] - P[0][1];
    triangle.a1->f_y += -P[1][0] - P[1][1];

    triangle.a2->f_x += P[0][0];
    triangle.a2->f_y += P[1][0];

    triangle.a3->f_x += P[0][1];
    triangle.a3->f_y += P[1][1];
}

void Cell::m_getDeformationGradiant(const Triangle &triangle)
{
    auto e1_ = triangle.e1();
    auto e2_ = triangle.e2();

    F[0][0] = e1_[0];
    F[1][0] = e1_[1];
    F[0][1] = e2_[0];
    F[1][1] = e2_[1];
}

void Cell::m_lagrangeReduction()
{
    // We start by copying the values from C to the reduced matrix
    C_ = C;

    if (LINEARITY)
    {
        // If we assume linearity, we are done. m is already identity.
        return;
    }
    // And then we follow an algorithm generate both m and C_
    while (C_[0][1] < 0 || C_[1][1] < C_[0][0] || 2 * C_[0][1] > C_[0][0])
    {

        if (C_[0][1] < 0)
        {
            C_.flip(0, 1);
            m.lag_m1();
        }

        if (C_[1][1] < C_[0][0])
        {
            C_.swap(0, 0, 1, 1);
            m.lag_m2();
        }

        if (2 * C_[0][1] > C_[0][0])
        {
            // The order here matters, don't modify C_[0][1] before using it
            // to calculate C_[1][1].
            C_[1][1] += C_[0][0] - 2 * C_[0][1];
            C_[0][1] -= C_[0][0];
            m.lag_m3();
        }
    }
}

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