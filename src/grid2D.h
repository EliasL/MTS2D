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



struct atom_id
{
    // atom_id references the index of a specific atom 
    // it should not be confused with the real position
    // of the atom. That is contained in atom.x and atom.y.
    int xi; // x-index range (0, n-1)
    int yi; // y-index range (0, m-1)
    int i; // index range (0, n*m-1)
    
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
    // The four nearest neighbours around the atom
    std::array<atom_id, 4> neighbours;
};

struct grid {
    // A 2D grid of atoms
    Matrix<atom> atoms;

    // A table choosing whether or not an atom is a border atom
    // We use uint8_t 1 and 0 instead of a bool
    Matrix<uint8_t> border;


    // Constructor that initializes the grid with size n x m
    grid(int n, int m): atoms(n, m), border(n, m) {
        setBorderElements();  // Call the function to set border elements to true
        fillNeighbours();
        // You can add additional initialization code here if needed
    }

    // With this overload, we can turn this:
    // grid.atoms.data[id.i].x
    // into this:
    // grid[id].x
    atom operator[](atom_id id) { return atoms.data[id.i]; }


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
                atoms[xi][yi].neighbours[0] = {xi, left, m};
                atoms[xi][yi].neighbours[1] = {xi, right, m};
                atoms[xi][yi].neighbours[3] = {down, yi, m};
                atoms[xi][yi].neighbours[2] = {up, yi, m};
            }
        }
    }

    bool isBorder(atom_id id){
        return border.data[id.i];
    }

};

void setAtomPositions(grid& g, double a) {
    int n = g.atoms.rows;
    int m = g.atoms.cols;
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            // Set the x and y positions based on the grid indices and spacing "a"
            g.atoms[i][j].x = i * a;
            g.atoms[i][j].y = j * a;
        }
    }
}

#endif