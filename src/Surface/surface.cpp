#include "surface.h"

// NODE AND NODE_ID ------

NodeId::NodeId() : xi(0), yi(0), i(0) {}
NodeId::NodeId(int xi_, int yi_, int cols) : xi(xi_), yi(yi_), i(xi_*cols + yi_) {}
NodeId::NodeId(int i_, int cols) : xi(i_ / cols), yi(i_ % cols), i(i_) {}

Node::Node(){}

void transformInPlace(const Matrix2x2<double> &matrix, Node &n) {
    double newX = matrix[0][0] * n.x + matrix[1][0] * n.y;
    double newY = matrix[0][1] * n.x + matrix[1][1] * n.y;
    n.x = newX;
    n.y = newY;
}

Node transform(const Matrix2x2<double> &matrix, const Node &n) {
    Node result = n;
    transformInPlace(matrix, result);
    return result;
}

void translateInPlace(Node &n, const Node &delta, double multiplier) {
    n.x += multiplier * delta.x;
    n.y += multiplier * delta.y;
}

Node translate(const Node &n, const Node &delta, double multiplier) {
    Node result = n;
    translateInPlace(result, delta, multiplier);
    return result;
}

// TRIANGLE -------

// e1 and e2 form basis vectors for the triangle
std::array<double, 2> Triangle::e1() const { // a is the lattice spacing of the gird
    return {
        a2->x - a1->x,
        a2->y - a1->y,
        };
}

std::array<double, 2> Triangle::e2() const { // a is the lattice spacing of the gird
    return {
        a3->x - a1->x,
        a3->y - a1->y,
        };
}

// Provices a metric tensor for the triangle
Matrix2x2<double> Triangle::metric(MetricFunction f) const {
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
        m[0][0] = e1_[0] -1;
        m[1][1] = e2_[1] -1;
        m[1][0] = m[0][1] = e2_[0];
    default:
        throw std::invalid_argument("Invalid metric function");
        break;
    }

}

// CELL --------------

Cell::Cell(Triangle t) {
    // Calculates D
    m_getDeformationGradiant(t);

    // Calculates C
    C  = t.metric(METRICFUNCTION);

    // Calculate C_ and m
    m_lagrangeReduction();

};

Cell::Cell(){};

// A basis vector for the cell
double Cell::e1(int index){
    return F[0][index];
}

// A basis vector for the cell
double Cell::e2(int index){
    return F[1][index];
}

// Calculate Piola stress tensor and force on each node from current cell
// Note that each node is part of multiple cells. Therefore, the force must
// be reset after each itteration.
void Cell::setForcesOnNodes(Triangle t){
    
    // extended stress is not quite the "real" stress, but it is a component
    // in calculating the piola stress, which is the real stress on the cell,
    // and we can then find the force on each individual node.
    // The name extended_stress does not have much meaning.
    // TODO consider storing this variable in the cell, such that it does
    // not have to be allocated every time the function is called.
    Matrix2x2<double> extended_stress = r_s.sym_orth_conjugate(m);
    
    P[0][0] = 2* extended_stress[0][0] * F[0][0] + extended_stress[0][1] * F[1][0];
    P[1][0] = 2* extended_stress[0][0] * F[0][1] + extended_stress[0][1] * F[1][1];
    P[0][1] = 2* extended_stress[1][1] * F[1][0] + extended_stress[0][1] * F[0][0];
    P[1][1] = 2* extended_stress[1][1] * F[1][1] + extended_stress[0][1] * F[0][1];
    // The assignment here is dependant on the shape of the cell.
    // For triangular shapes, the forces on the nodes is applied as shown
    // below. For a general shape, see Gael-notes page 2, partial N^i / partial x_j
    // on how to calculate. 

    // DO GENERAL DISTRIBUTION of Piola stress

    // FORCES SHOULD BE RESET AFTER EACH ITERATION
    // MUST SUM (+=) BECAUSE THEY ARE THE SAME NODES THAT ARE IN A FLIPPED TRIANGLE
    t.a1->f_x += -P[0][0] - P[0][1];
    t.a1->f_y += -P[1][0] - P[1][1];

    t.a2->f_x += P[0][0];
    t.a2->f_y += P[1][0];

    t.a3->f_x += P[0][1];
    t.a3->f_y += P[1][1];
}

void Cell::m_getDeformationGradiant(Triangle t){
    auto e1_ = t.e1();
    auto e2_ = t.e2();

    F[0][0] = e1_[0];
    F[0][1] = e1_[1];
    F[1][0] = e2_[0];
    F[1][1] = e2_[1];
}

void Cell::m_lagrangeReduction(){
    // We start by copying the values from C to the reduced matrix
    C_ = C;

    if(LINEARITY){
        // If we assume linearity, we are done. m is already identity.
        return;
    }
    // And then we follow an algorithm generate both m and C_
    while(C_[0][1]<0 || C_[1][1]<C_[0][0] || 2 * C_[0][1]>C_[0][0] ){
        
        if (C_[0][1]<0){
            C_.flip(0,1);
            m.lag_m1();
        }
        
        if (C_[1][1]<C_[0][0]){
            C_.swap(0,0,1,1);
            m.lag_m2();
        }

        if (2 * C_[0][1]>C_[0][0]){
            // The order here matters, don't modify C_[0][1] before using it
            // to calculate C_[1][1].
            C_[1][1] += C_[0][0] - 2 * C_[0][1];
            C_[0][1] -= C_[0][0];
            m.lag_m3();
        }
    } 
}

// BOUNDARY_CONDITIONS -------

BoundaryConditions::BoundaryConditions(double load, double theta, BCF bc) {
    this->load = load;
    this->theta = theta;
    m_calculateGradiant(bc); 
}

void BoundaryConditions::m_calculateGradiant(BCF bc){
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
void BoundaryConditions::m_macroShear(){
    double perturb = 0 ;

    F[0][0] = (1. - load*cos(theta + perturb)*sin(theta + perturb));
    F[1][1]  = (1. + load*cos(theta - perturb)*sin(theta - perturb));
    F[0][1]  = load* pow(cos(theta), 2.);
    F[1][0]  = -load* pow(sin(theta - perturb), 2.);
}

// GRID ----
Mesh::Mesh(){}

// Constructor that initializes the surface with size n x m
Mesh::Mesh(int n, int m, double a): nodes(n, m), a(a),
    nrTriangles(2*(n-1)*(m-1)), triangles(2*(n-1)*(m-1)), cells(2*(n-1)*(m-1)){
    // These functions loop over the same elements, and we
    // could be slightly more optimized by combining everything 
    // into one loop, but this is more readable, and there is no
    // need to optimize the constructor since we only construct one surface.
    m_setBorderElements();
    m_fillNonBorderNodeIds();
    m_setNodePositions();
    m_fillNeighbours();
    m_createTriangles();
}

Mesh::Mesh(int n, int m) : Mesh(n, m, 1){}

bool Mesh::isBorder(NodeId n_id){
    return (*this)[n_id]->borderNode;
}


void Mesh::applyBoundaryConditions(BoundaryConditions bc){
    // We get the id of each node in the border
    for(NodeId n_id : borderNodeIds){
        transformInPlace(bc.F, *(*this)[n_id]);
    }
}

void Mesh::resetForceOnNodes(){
    for(Node n : nodes.data){
        n.f_x=n.f_y=0;
    }
}

// This is just a function to avoid having to write nodes.cols
NodeId Mesh::getNodeId(int xi, int yi) {
    return NodeId(xi, yi, nodes.cols);
}

// Function to set border elements of the border vector to true
void Mesh::m_setBorderElements() {
    int n = nodes.rows;
    int m = nodes.cols;

    // Loop over the border elements only
    for (int i = 0; i < n; ++i) {
        nodes[i][0].borderNode = true;
        nodes[i][m - 1].borderNode = true;
    }
    for (int j = 0; j < m; ++j) {
        nodes[0][j].borderNode = true;
        nodes[n - 1][j].borderNode = true;
    }
}

void Mesh::m_fillNonBorderNodeIds(){
    for (int i = 0; i < nodes.data.size(); i++)
    {
        NodeId n_id = NodeId{i, nodes.cols};
        // If n_id is a border node,
        if (isBorder(n_id)){
            // we add it to the vector.
            borderNodeIds.push_back(NodeId(i, nodes.cols));
        } else {
            nonBorderNodeIds.push_back(NodeId(i, nodes.cols));
        }
    }
}

void Mesh::m_setNodePositions() {
    int n = nodes.rows;
    int m = nodes.cols;

    for (int xi = 0; xi < n; ++xi) {
        for (int yi = 0; yi < m; ++yi) {
            // Set the x and y positions based on the surface indices and spacing "a"
            nodes[xi][yi].x = xi * a;
            nodes[xi][yi].y = yi * a;
            nodes[xi][yi].id = getNodeId(xi, yi);
        }
    }
}

// Function to fill neighbours using periodic boundary conditions
void Mesh::m_fillNeighbours() {
    int n = nodes.rows;
    int m = nodes.cols;

    for (int xi = 0; xi < n; ++xi) {
        for (int yi = 0; yi < m; ++yi) {
            // Define neighbor indices using periodic boundary conditions
            int left = (yi == 0) ? m - 1 : yi - 1;
            int right = (yi == m - 1) ? 0 : yi + 1;
            int up = (xi == 0) ? n - 1 : xi - 1;
            int down = (xi == n - 1) ? 0 : xi + 1;

            // Fill in the neighbors
            nodes[xi][yi].neighbours[0] = getNodeId(xi, left);
            nodes[xi][yi].neighbours[1] = getNodeId(xi, right);
            nodes[xi][yi].neighbours[3] = getNodeId(down, yi);
            nodes[xi][yi].neighbours[2] = getNodeId(up, yi);
        }
    }
}

// See the bottom of the doc for explination
void Mesh::m_createTriangles(){
    int n = nodes.rows;
    int m = nodes.cols;

    // CHECK THAT CALCULATIONS CAN BE DONE ON UPSIDE DOWN TRIANGLES WITHOUT MINUS
    for (int xi = 0; xi < n-1; ++xi) {
        for (int yi = 0; yi < m-1; ++yi) {
            // We now find the 4 nodes in the current square
            NodeId a1 = getNodeId(xi, yi);
            NodeId a2 = getNodeId(xi, yi+1);
            NodeId a3 = getNodeId(xi+1, yi);
            NodeId a4 = getNodeId(xi+1, yi+1);

            int t1i = 2*(xi*(m-1)+yi); // Triangle 1 index
            int t2i = t1i+1;

            triangles[t1i] = Triangle{(*this)[a1], (*this)[a2], (*this)[a3]};
            triangles[t2i] = Triangle{(*this)[a2], (*this)[a3], (*this)[a4]};
        }
    } 
}


