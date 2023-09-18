#ifndef GRID2D_H
#define GRID2D_H
#pragma once

/*

A grid consists of many atoms. The grid keeps track of what
atoms are at the border and what atoms are neighbours

The grid will also create triangles in the grid and keep track
of these.
*/

#include "matrix.h" // Include matrix.h here after the class definition
#include "matrix2x2.h"

enum MetricFunction { faicella, epsilon_lineaire };


struct atom_id
{
    // atom_id references the index of a specific atom 
    // it should not be confused with the real position
    // of the atom. That is contained in atom.x and atom.y.
    int xi; // x-index range (0, n-1)
    int yi; // y-index range (0, m-1)
    int i; // index range (0, n*m-1)
    
    // It is quite annoying to have to pass cols 
    atom_id() : xi(0), yi(0), i(0) {}
    atom_id(int xi, int yi, int cols) : xi(xi), yi(yi), i(xi*cols + yi) {}
    atom_id(int i, int cols) : xi(i / cols), yi(i % cols), i(i) {}
};

struct atom
{
    // The possition of an atom
    // Can perhaps also contain size and mass
    double x;
    double y;
    // double size;
    // double mass;

    // id to itself
    atom_id id;
    // The four nearest neighbours around the atom
    std::array<atom_id, 4> neighbours;
};

struct triangle
{
    atom* a1;
    atom* a2;
    atom* a3;

    // e1 and e2 form basis vectors for the triangle
    std::array<double, 2> e1(double a) const { // a is the lattice spacing of the gird
        return {
            a*( a2->x - a1->x),
            a*( a2->y - a1->y),
            };
    }

    std::array<double, 2> e2(double a) const { // a is the lattice spacing of the gird
        return {
            a*( a3->x - a1->x),
            a*( a3->y - a1->y),
            };
    }

    // Provices a metric tensor for the triangle
    Matrix2x2<double> metric(double a, MetricFunction f = MetricFunction::faicella) const {
        // Symetric matricies would be faster, but only slightly for 2x2 matrix
        Matrix2x2<double> m;
        auto e1_ = e1(a);
        auto e2_ = e2(a);
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

class Grid {
public:
    // A 2D grid of atoms
    Matrix<atom> atoms;

    // Triangles (See bottom for explination)
    //std::vector<std::array<atom_id, 3>> triangle_ids;
    std::vector<triangle> triangles;

    // A table choosing whether or not an atom is a border atom
    // We use uint8_t 1 and 0 instead of a bool
    Matrix<uint8_t> border;
    
    // Initial (average) distance between atoms
    double a = 0;

    // Default Constructor
    Grid(): atoms(0, 0), border(0, 0), triangles(0){}//, triangle_ids(0){}

    // Constructor that initializes the grid with size n x m
    Grid(int n, int m, double a): atoms(n, m), border(n, m), a(a),
        triangles(2*(n-1)*(m-1)){//, triangle_ids(2*(n-1)*(m-1)) {
        
        // These functions loop over the same elements, and we
        // could be slightly more optimized by combining everything 
        // into one loop, but this is more readable, and there is no
        // need to optimize the constructor since we only construct one grid.
        setBorderElements();
        setAtomPositions();
        fillNeighbours();
        createTriangles();
    }

    Grid(int n, int m) : Grid(n, m, 1){}

    // With this overload, we can turn this: 
    // grid.atoms.data[id.i].x
    // into this:
    // grid[id]->x
    atom* operator[](atom_id id) { return &atoms.data[id.i]; }

    bool isBorder(atom_id id){
        return border.data[id.i];
    }

private:
    atom_id atom_id_(int xi, int yi) {
        return atom_id(xi, yi, atoms.cols);
    }

    // Function to set border elements of the border vector to true
    void setBorderElements() {
        int n = border.rows;
        int m = border.cols;

        // Loop over the border elements only
        for (int i = 0; i < n; ++i) {
            border[i][0] = true;
            border[i][m - 1] = true;
        }
        for (int j = 0; j < m; ++j) {
            border[0][j] = true;
            border[n - 1][j] = true;
        }
    }

    void setAtomPositions() {
        int n = atoms.rows;
        int m = atoms.cols;

        for (int xi = 0; xi < n; ++xi) {
            for (int yi = 0; yi < m; ++yi) {
                // Set the x and y positions based on the grid indices and spacing "a"
                atoms[xi][yi].x = xi * a;
                atoms[xi][yi].y = yi * a;
                atoms[xi][yi].id = atom_id_(xi, yi);
            }
        }
    }

    // Function to fill neighbours using periodic boundary conditions
    void fillNeighbours() {
        int n = atoms.rows;
        int m = atoms.cols;

        for (int xi = 0; xi < n; ++xi) {
            for (int yi = 0; yi < m; ++yi) {
                // Define neighbor indices using periodic boundary conditions
                int left = (yi == 0) ? m - 1 : yi - 1;
                int right = (yi == m - 1) ? 0 : yi + 1;
                int up = (xi == 0) ? n - 1 : xi - 1;
                int down = (xi == n - 1) ? 0 : xi + 1;

                // Fill in the neighbors
                atoms[xi][yi].neighbours[0] = atom_id_(xi, left);
                atoms[xi][yi].neighbours[1] = atom_id_(xi, right);
                atoms[xi][yi].neighbours[3] = atom_id_(down, yi);
                atoms[xi][yi].neighbours[2] = atom_id_(up, yi);
            }
        }
    }

    void createTriangles(){
        // See the bottom of the doc for explination
        int n = atoms.rows;
        int m = atoms.cols;


        for (int xi = 0; xi < n-1; ++xi) {
            for (int yi = 0; yi < m-1; ++yi) {
                // We now find the 4 atoms in the current square
                atom_id a1 = atom_id_(xi, yi);
                atom_id a2 = atom_id_(xi, yi+1);
                atom_id a3 = atom_id_(xi+1, yi);
                atom_id a4 = atom_id_(xi+1, yi+1);

                int t1i = 2*(xi*(m-1)+yi); // Triangle 1 index
                int t2i = t1i+1;
                //triangle_ids[t1i] = {a1, a2, a3};
                //triangle_ids[t2i] = {a2, a3, a4};

                triangles[t1i] = {(*this)[a1], (*this)[a2], (*this)[a3]};
                triangles[t2i] = {(*this)[a2], (*this)[a3], (*this)[a4]};
            }
        } 
    }

};

#endif






/*
Triangles explination.

We have a grid of atoms with a unique id, as an example, consider this 2x3 grid:

    1    2

    3    4

    5    6

This grid will have 4 triangles, and we will now see how we can construct them.
Notice first that we can consider a grid of squares made of atoms instead of
looking at the grid of atoms themselves. We have a 1x2 grid of squares. The two
squares in question are (1,2,3,4) and (3,4,5,6). We will now fit two triangles
in each of the squares. In the (1,2,3,4) square we have the (1,2,3) triangle, we
refer to this as an up-triangle, and (2,3,4) as the down triangle. Similarly, 
the bottom square (3,4,5,6) contains the triangles (3,4,5) and (4,5,6). 

Instead of looping over each atom, we loop over all the (n-1)x(m-1) squares and
create a up and down triangle in each itteration of the loop. 
*/