#include "run/doctest.h"
#include "../src/grid2D.h"  // Include the header for your grid struct

TEST_CASE("Grid Initialization") {
    // Create a grid with known dimensions for testing
    grid g(4, 4);

    // Verify that the grid dimensions are correct
    // Verify that the grid dimensions are correct
    REQUIRE(g.atoms.size() == 4);
    REQUIRE(g.atoms[0].size() == 4);
    REQUIRE(g.border.size() == 4);
    REQUIRE(g.border[0].size() == 4);
}

TEST_CASE("Accessing Grid Elements") {
    grid g(3, 3);

    // Modify an atom
    g.atoms[1][1].x = 5.0;
    g.atoms[1][1].y = 10.0;

    // Verify that the modification is reflected in the grid
    REQUIRE(g.atoms[1][1].x == 5.0);
    REQUIRE(g.atoms[1][1].y == 10.0);

    // Verify some other elements
    REQUIRE(g.atoms[0][0].x == 0.0);
    REQUIRE(g.atoms[2][2].y == 0.0);
}

// Test case for checking neighbors with periodic boundary conditions
TEST_CASE("Neighbors with Periodic Boundary Conditions") {
    grid g(3, 3);

    // Central atom
    std::array<atom_id, 4> neighbors_center = g.atoms[1][1].neighbours;
    REQUIRE(neighbors_center.size() == 4);
    REQUIRE(neighbors_center[0].i == 1);
    REQUIRE(neighbors_center[0].j == 0); // Left
    REQUIRE(neighbors_center[1].i == 1);
    REQUIRE(neighbors_center[1].j == 2); // Right
    REQUIRE(neighbors_center[2].i == 0);
    REQUIRE(neighbors_center[2].j == 1); // Up
    REQUIRE(neighbors_center[3].i == 2);
    REQUIRE(neighbors_center[3].j == 1); // Down

    // Corner atom [0][0]
    std::array<atom_id, 4> neighbors_corner = g.atoms[0][0].neighbours;
    REQUIRE(neighbors_corner.size() == 4);
    REQUIRE(neighbors_corner[0].i == 0);
    REQUIRE(neighbors_corner[0].j == 2); // Left (wrap around)
    REQUIRE(neighbors_corner[1].i == 0);
    REQUIRE(neighbors_corner[1].j == 1); // Right
    REQUIRE(neighbors_corner[2].i == 2);
    REQUIRE(neighbors_corner[2].j == 0); // Up (wrap around)
    REQUIRE(neighbors_corner[3].i == 1);
    REQUIRE(neighbors_corner[3].j == 0); // Down
}

// Test case for setting atom positions in a regular square grid
TEST_CASE("Setting Atom Positions in a Regular Grid") {
    grid g(4, 4);
    double spacing = 1.0; // Spacing between atoms

    // Call the function to set atom positions
    setAtomPositions(g, spacing);

    // Verify that atom positions are correctly set
    REQUIRE(g.atoms[0][0].x == 0.0);
    REQUIRE(g.atoms[0][0].y == 0.0);


    REQUIRE(g.atoms[2][2].x == 2 * spacing);
    REQUIRE(g.atoms[2][2].y == 2 * spacing);

    REQUIRE(g.atoms[2][3].x == 2 * spacing);
    REQUIRE(g.atoms[2][3].y == 3 * spacing);

    // You can add more checks as needed
}