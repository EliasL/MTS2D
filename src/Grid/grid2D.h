#ifndef GRID2D_H
#define GRID2D_H
#pragma once

#include "../settings.h"
#include "../Matrix/matrix.h"
#include "../Matrix/matrix2x2.h"
#include <array>
#include <vector>
#include <stdexcept>

struct node_id
{
    int xi, yi, i;
    node_id();
    node_id(int xi, int yi, int cols);
    node_id(int i, int cols);
};

struct node
{
    double x, y, f_x, f_y, border_node = false;
    node_id id;
    std::array<node_id, 4> neighbours;
    node();
    node(double x, double y);
};

node transform(const Matrix2x2<double> &matrix, const node &n);
node translate(const node &n, const node &delta, double multiplier = 1);

struct triangle
{
    node *a1;
    node *a2;
    node *a3;
    std::array<double, 2> e1() const;
    std::array<double, 2> e2() const;
    Matrix2x2<double> metric(MetricFunction f = MetricFunction::faicella) const;
};

class Cell
{
public:
    Matrix2x2<double> F, C, C_, m, r_s, P;
    double energy;
    bool plasticity;
    Cell(triangle t);
    Cell();
    double e1(int index);
    double e2(int index);
    void set_forces_on_nodes(triangle t);

private:
    void get_deformation_gradiant(triangle t);
    void lagrange_reduction();
};

struct boundary_conditions
{
    double load, theta;
    Matrix2x2<double> F;
    BoundaryConditionFunction bcFun;
    boundary_conditions(double load, double theta);
    void calculate_gradiant();
    void macro_shear();
};

class Grid
{
public:
    Matrix<node> nodes;
    std::vector<triangle> triangles;
    std::vector<Cell> cells;
    std::vector<node_id> border_node_ids, non_border_node_ids;
    double a, load;
    int nr_triangles;
    Grid();
    Grid(int n, int m, double a);
    Grid(int n, int m);

    // With this overload, we can turn this:
    // grid.nodes.data[id.i].x
    // into this:
    // grid[id]->x
    node *operator[](node_id id) { return &nodes.data[id.i]; }
    const node *operator[](node_id id) const { return &nodes.data[id.i]; }

    bool isBorder(node_id n_id);
    void apply_boundary_conditions(boundary_conditions bc);
    void reset_force_on_nodes();
    void setBorderElements();
    void fill_non_border_node_ids();
    void setNodePositions();
    void fillNeighbours();
    void createTriangles();
    node_id node_id_(int xi, int yi);
};
#endif
