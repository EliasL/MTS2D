#ifndef GRID2D_H
#define GRID2D_H
#pragma once

/*

A grid consists of many nodes. The grid keeps track of what
nodes are at the border and what nodes are neighbours

The grid will also create triangles in the grid and keep track
of these.
*/

#include "settings.h"
#include "matrix.h" // Include matrix.h here after the class definition
#include "matrix2x2.h"

// Node_id references the index of a specific node,
// it should not be confused with the real position
// of the node. That is contained in node.x and node.y.
struct node_id
{
    // x-index range (0, n-1)
    int xi;
    // y-index range (0, m-1)
    int yi;

    // index range (0, n*m-1)
    int i;
    
    node_id() : xi(0), yi(0), i(0) {}
    node_id(int xi, int yi, int cols) : xi(xi), yi(yi), i(xi*cols + yi) {}
    node_id(int i, int cols) : xi(i / cols), yi(i % cols), i(i) {}
};

// Properties of an node.
struct node
{
    // x position
    double x;
    // y position
    double y;

    // x force
    double f_x;
    // y force
    double f_y;

    // Unused variables
    // double v_x; // x velocity
    // double v_y; // y velocity
    // double size;
    // double mass;

    // This determines whether or not the node is a border node
    double border_node=false;

    // id to itself
    node_id id;
    // The four nearest neighbours around the node
    std::array<node_id, 4> neighbours;

    node(){}
    node(double x, double y): x(x), y(y) {}
};

node transform(const Matrix2x2<double>& matrix, const node& n) {
    node result;
    result.x = matrix[0][0] * n.x + matrix[1][0] * n.y;
    result.y = matrix[0][1] * n.x + matrix[1][1] * n.y;
    return result;
}

node translate(const node& n, const node& delta, double multiplier = 1) {
    node result;
    result.x = n.x + multiplier * delta.x;
    result.y = n.y + multiplier * delta.y;
    return result;
}







// References to three nodes that form a triangle in the grid.
struct triangle
{
    node* a1;
    node* a2;
    node* a3;

    // e1 and e2 form basis vectors for the triangle
    std::array<double, 2> e1() const { // a is the lattice spacing of the gird
        return {
            a2->x - a1->x,
            a2->y - a1->y,
            };
    }

    std::array<double, 2> e2() const { // a is the lattice spacing of the gird
        return {
            a3->x - a1->x,
            a3->y - a1->y,
            };
    }

    // Provices a metric tensor for the triangle
    Matrix2x2<double> metric(MetricFunction f = MetricFunction::faicella) const {
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
};

/**
 * A finite element cell constructed from a triangle of three nodes
 */

class Cell{
public:

    Matrix2x2<double> F; // deformation gradiant / basis vectors (e1, e2)
    Matrix2x2<double> C; // real_metrics (T(F)F)
    Matrix2x2<double> C_; // reduced_metrics
    Matrix2x2<double> m; // reduction transformation ( T(m)Cm = C_ )
    Matrix2x2<double> r_s; // reduces stress 
    // Piola-Kirchhoff stress: https://en.wikipedia.org/wiki/Piola%E2%80%93Kirchhoff_stress_tensors
    Matrix2x2<double> P; // Piola stress tensor (Piola stress is just normal stress, but calculated in a special way)
    double energy;
    bool plasticity; // not used yet

    Cell(triangle t) {
        // Calculates D
        get_deformation_gradiant(t);

        // Calculates C
        C  = t.metric(METRICFUNCTION);

        // Calculate C_ and m
        lagrange_reduction();

    };

    // Don't use
    Cell(){};

    // A basis vector for the cell
    double e1(int index){
        return F[0][index];
    }

    // A basis vector for the cell
    double e2(int index){
        return F[1][index];
    }

    // Calculate Piola stress tensor and force on each node from current cell
    // Note that each node is part of multiple cells. Therefore, the force must
    // be reset after each itteration.
    void set_forces_on_nodes(triangle t){
        
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

private:
    void get_deformation_gradiant(triangle t){
        auto e1_ = t.e1();
        auto e2_ = t.e1();

        F[0][0] = e1_[0];
        F[0][1] = e1_[1];
        F[1][0] = e2_[0];
        F[1][1] = e2_[1];
    }

   void lagrange_reduction(){
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
};


struct boundary_conditions {
    double load;
	double theta; // rotation of shear
    Matrix2x2<double> F; // Deformation gradiant.
    BoundaryConditionFunction bcFun = BOUNDARYCONDITIONFUNCTION;

    boundary_conditions(double load, double theta) : load(load), theta(theta){
        calculate_gradiant();
    }

    void calculate_gradiant(){
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

    void macro_shear(){
        double perturb = 0 ;

        F[0][0] = (1. - load*cos(theta + perturb)*sin(theta + perturb));
        F[1][1]  = (1. + load*cos(theta - perturb)*sin(theta - perturb));
        F[0][1]  = load* pow(cos(theta), 2.);
        F[1][0]  = -load* pow(sin(theta - perturb), 2.);
    }

};

// A 2D grid of nodes
class Grid {
public:
    Matrix<node> nodes;

    // Triangles (See bottom for explination)
    //std::vector<std::array<node_id, 3>> triangle_ids;
    std::vector<triangle> triangles;

    // Cells are triangles with extra numbers.
    // More precicely, they are the finite elements constructed
    // from a given triangle
    std::vector<Cell> cells;
    
    // All node_ids coresponding to nodes on the border of the system
    std::vector<node_id> border_node_ids;
    // Non-border-nodes (Used in energy minimization solver)
    std::vector<node_id> non_border_node_ids;


    // Initial distance between nearest neighbours
    double a;

    // The load on the boundary conditions
    double load;

    int nr_triangles;

    // Default Constructor
    Grid(){}

    // Constructor that initializes the grid with size n x m
    Grid(int n, int m, double a): nodes(n, m), a(a),
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

    Grid(int n, int m) : Grid(n, m, 1){}

    // With this overload, we can turn this: 
    // grid.nodes.data[id.i].x
    // into this:
    // grid[id]->x
    node* operator[](node_id id) { return &nodes.data[id.i]; }
    const node* operator[](node_id id) const { return &nodes.data[id.i]; }

    bool isBorder(node_id n_id){
        return (*this)[n_id]->border_node;
    }


    void apply_boundary_conditions(boundary_conditions bc){
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

    void reset_force_on_nodes(){
        for(node n : nodes.data){
            n.f_x=n.f_y=0;
        }
    }

private:

    // This is just a function to avoid having to write nodes.cols
    node_id node_id_(int xi, int yi) {
        return node_id(xi, yi, nodes.cols);
    }

    // Function to set border elements of the border vector to true
    void setBorderElements() {
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

    void fill_non_border_node_ids(){
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

    void setNodePositions() {
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
    void fillNeighbours() {
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
    void createTriangles(){
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

};

#endif






/*
Triangles explination.

We have a grid of nodes with a unique id, as an example, consider this 2x3 grid:

    1    2

    3    4

    5    6

This grid will have 4 triangles, and we will now see how we can construct them.
Notice first that we can consider a grid of squares made of nodes instead of
looking at the grid of nodes themselves. We have a 1x2 grid of squares. The two
squares in question are (1,2,3,4) and (3,4,5,6). We will now fit two triangles
in each of the squares. In the (1,2,3,4) square we have the (1,2,3) triangle, we
refer to this as an up-triangle, and (2,3,4) as the down triangle. Similarly, 
the bottom square (3,4,5,6) contains the triangles (3,4,5) and (4,5,6). 

Instead of looping over each node, we loop over all the (n-1)x(m-1) squares and
create a up and down triangle in each itteration of the loop. 
*/