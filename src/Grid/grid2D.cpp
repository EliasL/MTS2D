#include "grid2D.h"

// NODE AND NODE_ID ------

node_id::node_id() : xi(0), yi(0), i(0) {}
node_id::node_id(int xi_, int yi_, int cols) : xi(xi_), yi(yi_), i(xi_*cols + yi_) {}
node_id::node_id(int i_, int cols) : xi(i_ / cols), yi(i_ % cols), i(i_) {}

node::node(){}

node transform(const Matrix2x2<double>& matrix, const node& n) {
    node result;
    result.x = matrix[0][0] * n.x + matrix[1][0] * n.y;
    result.y = matrix[0][1] * n.x + matrix[1][1] * n.y;
    return result;
}

node translate(const node& n, const node& delta, double multiplier) {
    node result;
    result.x = n.x + multiplier * delta.x;
    result.y = n.y + multiplier * delta.y;
    return result;
}

// TRIANGLE -------

// e1 and e2 form basis vectors for the triangle
std::array<double, 2> triangle::e1() const { // a is the lattice spacing of the gird
    return {
        a2->x - a1->x,
        a2->y - a1->y,
        };
}

std::array<double, 2> triangle::e2() const { // a is the lattice spacing of the gird
    return {
        a3->x - a1->x,
        a3->y - a1->y,
        };
}

// Provices a metric tensor for the triangle
Matrix2x2<double> triangle::metric(MetricFunction f) const {
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

Cell::Cell(triangle t) {
    // Calculates D
    get_deformation_gradiant(t);

    // Calculates C
    C  = t.metric(METRICFUNCTION);

    // Calculate C_ and m
    lagrange_reduction();

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
void Cell::set_forces_on_nodes(triangle t){
    
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

void Cell::get_deformation_gradiant(triangle t){
    auto e1_ = t.e1();
    auto e2_ = t.e1();

    F[0][0] = e1_[0];
    F[0][1] = e1_[1];
    F[1][0] = e2_[0];
    F[1][1] = e2_[1];
}

void Cell::lagrange_reduction(){
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

boundary_conditions::boundary_conditions(double load, double theta) {
    // Your implementation here
    this->load = load;
    this->theta = theta;
    // ...Initialize other members as needed...
}

void boundary_conditions::calculate_gradiant(){
    switch (bcFun)
    {
    case BoundaryConditionFunction::macro_shear:
        macro_shear();
        break;
    
    default:
        throw std::invalid_argument("Invalid boundary condition function");
        break;
    }
}
/*
*/
void boundary_conditions::macro_shear(){
    double perturb = 0 ;

    F[0][0] = (1. - load*cos(theta + perturb)*sin(theta + perturb));
    F[1][1]  = (1. + load*cos(theta - perturb)*sin(theta - perturb));
    F[0][1]  = load* pow(cos(theta), 2.);
    F[1][0]  = -load* pow(sin(theta - perturb), 2.);
}

// GRID ----
Grid::Grid(){}

// Constructor that initializes the grid with size n x m
Grid::Grid(int n, int m, double a): nodes(n, m), a(a),
    nr_triangles(2*(n-1)*(m-1)), triangles(2*(n-1)*(m-1)), cells(2*(n-1)*(m-1)){
    // These functions loop over the same elements, and we
    // could be slightly more optimized by combining everything 
    // into one loop, but this is more readable, and there is no
    // need to optimize the constructor since we only construct one grid.
    setBorderElements();
    fill_non_border_node_ids();
    setNodePositions();
    fillNeighbours();
    createTriangles();
}

Grid::Grid(int n, int m) : Grid(n, m, 1){}

bool Grid::isBorder(node_id n_id){
    return (*this)[n_id]->border_node;
}


void Grid::apply_boundary_conditions(boundary_conditions bc){
    // We get the id of each node in the border
    for(node_id n_id : border_node_ids){
        double x = (*this)[n_id]->x; 
        double y = (*this)[n_id]->y; 
        double temp_x = x*bc.F[0][0] + y*bc.F[0][1];
        double temp_y = x*bc.F[1][0] + y*bc.F[1][1];
        (*this)[n_id]->x = temp_x; 
        (*this)[n_id]->y = temp_y; 
    }
}

void Grid::reset_force_on_nodes(){
    for(node n : nodes.data){
        n.f_x=n.f_y=0;
    }
}

// This is just a function to avoid having to write nodes.cols
node_id Grid::node_id_(int xi, int yi) {
    return node_id(xi, yi, nodes.cols);
}

// Function to set border elements of the border vector to true
void Grid::setBorderElements() {
    int n = nodes.rows;
    int m = nodes.cols;

    // Loop over the border elements only
    for (int i = 0; i < n; ++i) {
        nodes[i][0].border_node = true;
        nodes[i][m - 1].border_node = true;
    }
    for (int j = 0; j < m; ++j) {
        nodes[0][j].border_node = true;
        nodes[n - 1][j].border_node = true;
    }
}

void Grid::fill_non_border_node_ids(){
    for (int i = 0; i < nodes.data.size(); i++)
    {
        node_id n_id = node_id{i, nodes.cols};
        // If n_id is a border node,
        if (isBorder(n_id)){
            // we add it to the vector.
            border_node_ids.push_back(node_id(i, nodes.cols));
        } else {
            non_border_node_ids.push_back(node_id(i, nodes.cols));
        }
    }
}

void Grid::setNodePositions() {
    int n = nodes.rows;
    int m = nodes.cols;

    for (int xi = 0; xi < n; ++xi) {
        for (int yi = 0; yi < m; ++yi) {
            // Set the x and y positions based on the grid indices and spacing "a"
            nodes[xi][yi].x = xi * a;
            nodes[xi][yi].y = yi * a;
            nodes[xi][yi].id = node_id_(xi, yi);
        }
    }
}

// Function to fill neighbours using periodic boundary conditions
void Grid::fillNeighbours() {
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
            nodes[xi][yi].neighbours[0] = node_id_(xi, left);
            nodes[xi][yi].neighbours[1] = node_id_(xi, right);
            nodes[xi][yi].neighbours[3] = node_id_(down, yi);
            nodes[xi][yi].neighbours[2] = node_id_(up, yi);
        }
    }
}

// See the bottom of the doc for explination
void Grid::createTriangles(){
    int n = nodes.rows;
    int m = nodes.cols;

    // CHECK THAT CALCULATIONS CAN BE DONE ON UPSIDE DOWN TRIANGLES WITHOUT MINUS
    for (int xi = 0; xi < n-1; ++xi) {
        for (int yi = 0; yi < m-1; ++yi) {
            // We now find the 4 nodes in the current square
            node_id a1 = node_id_(xi, yi);
            node_id a2 = node_id_(xi, yi+1);
            node_id a3 = node_id_(xi+1, yi);
            node_id a4 = node_id_(xi+1, yi+1);

            int t1i = 2*(xi*(m-1)+yi); // Triangle 1 index
            int t2i = t1i+1;

            triangles[t1i] = triangle{(*this)[a1], (*this)[a2], (*this)[a3]};
            triangles[t2i] = triangle{(*this)[a2], (*this)[a3], (*this)[a4]};
        }
    } 
}


