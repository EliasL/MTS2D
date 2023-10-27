#include "matrix.h"

// Define a 2x2 matrix
template <typename T>
class Matrix2x2 {
public:
    double data[2][2];

    Matrix2x2() {
        data[0][0] = 1.0;
        data[0][1] = 0.0;
        data[1][0] = 0.0;
        data[1][1] = 1.0;
    }

    // Define the operator[] to access matrix elements
    T* operator[](int row) { return &data[row*2]; }

    // Const version of operator[]
    const T* operator[](int row) const { return &data[row * 2]; }

};

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

struct node {
    double x;
    double y;

    // Constructor
    node(double xVal, double yVal) : x(xVal), y(yVal) {}
    node(){}

    // Overload multiplication operator for node and Matrix2x2
    node operator*(const Matrix2x2<double>& matrix) const {
        node result;
        result.x = matrix[0][0] * x + matrix[0][1] * y;
        result.y = matrix[1][0] * x + matrix[1][1] * y;
        // Multiply other properties as needed
        return result;
    }

    // Overload subtraction operator for node objects
    node operator-(const node& other) const {
        node result;
        result.x = x - other.x;
        result.y = y - other.y;
        // Subtract other properties as needed
        return result;
    }

    // Overload subtraction equals operator for node objects
    node& operator-=(const node& other){
        node result;
        x -= other.x;
        y -= other.y;
        // Subtract other properties as needed
        return *this;
    }
};

struct boundary_conditions {
    double load;
	double theta; // rotation of shear
    Matrix2x2<double> F; // Deformation gradiant.

    boundary_conditions(double load, double theta) : load(load), theta(theta){
    }

    void macro_shear(){
        double perturb = 0 ;

        F[0][0] = (1. - load*cos(theta + perturb)*sin(theta + perturb));
        F[1][1]  = (1. + load*cos(theta - perturb)*sin(theta - perturb));
        F[0][1]  = load* pow(cos(theta), 2.);
        F[1][0]  = -load* pow(sin(theta - perturb), 2.);
    }

};

class Grid {
public:
    Matrix<node> nodes;

    // Non-border-nodes (Used in energy minimization solver)
    std::vector<node_id> non_border_node_ids;

    // Initial distance between nearest neighbours
    double a;

    // Default Constructor
    Grid(){}

    // Constructor that initializes the grid with size n x m
    Grid(int n, int m, double a): nodes(n, m), a(a){

    }

    Grid(int n, int m) : Grid(n, m, 1){}

    // With this overload, we can turn this: 
    // grid.nodes.data[id.i].x
    // into this:
    // grid[id]->x
    node* operator[](node_id id) { return &nodes.data[id.i]; }
    const node* operator[](node_id id) const { return &nodes.data[id.i]; }


private:

    // This is just a function to avoid having to write nodes.cols
    node_id node_id_(int xi, int yi) {
        return node_id(xi, yi, nodes.cols);
    }

};

void initial_guess(const Grid& g, const boundary_conditions& bc){

	int nr_x_elements = 5;
    node transformed_node;
    for (size_t i = 0; i < g.non_border_node_ids.size(); i++)
    {
        
        transformed_node = bc.F * g[g.non_border_node_ids[i]];
        transformed_node -= g[g.non_border_node_ids[i]];
    }
}



