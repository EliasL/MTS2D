/*

A grid consists of many atoms. The grid keeps track of what
atoms are at the border and what atoms are neighbours

The grid will also create triangles in the grid and keep track
of these.
*/
#pragma once
#include <vector>


struct atom_id
{
    // atom_id references the index of a specific atom 
    // it should not be confused with the real position
    // of the atom. That is contained in atom.x and atom.y.
    int i;
    int j;
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
    std::vector<std::vector<atom>> atoms;

    // A table choosing whether or not an atom is a border atom
    std::vector<std::vector<bool>> border;

    // Constructor that initializes the grid with size n x m
    grid(int n, int m)
        : atoms(n, std::vector<atom>(m)),
          border(n, std::vector<bool>(m)) {

        setBorderElements();  // Call the function to set border elements to true
        fillNeighbours();
        // You can add additional initialization code here if needed
    }

    // Function to set border elements of the border vector to true
    void setBorderElements() {
        int n = border.size();
        int m = border[0].size();

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
        int n = atoms.size();
        int m = atoms[0].size();

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                // Define neighbor indices using periodic boundary conditions
                int left = (j == 0) ? m - 1 : j - 1;
                int right = (j == m - 1) ? 0 : j + 1;
                int up = (i == 0) ? n - 1 : i - 1;
                int down = (i == n - 1) ? 0 : i + 1;

                // Fill in the neighbors
                atoms[i][j].neighbours[0] = {i, left};
                atoms[i][j].neighbours[1] = {i, right};
                atoms[i][j].neighbours[2] = {up, j};
                atoms[i][j].neighbours[3] = {down, j};
            }
        }
    }
};

void setAtomPositions(grid& g, double a) {
    int n = g.atoms.size();
    int m = g.atoms[0].size();
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            // Set the x and y positions based on the grid indices and spacing "a"
            g.atoms[i][j].x = i * a;
            g.atoms[i][j].y = j * a;
        }
    }
}