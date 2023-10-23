#include "run/doctest.h"
#include "../src/grid2D.h"  // Include the header for your grid struct

TEST_CASE("Grid Initialization") {
    // Create a grid with known dimensions for testing
    Grid g(4, 4);

    // Verify that the grid dimensions are correct
    // Verify that the grid dimensions are correct
    REQUIRE(g.nodes.rows == 4);
    REQUIRE(g.nodes.cols == 4);
}

// Test case for the node_id struct
TEST_CASE("node_id Struct Test") {
    /*
    Visualization of xi, yi and i in node_id for 3x3 matrix

    x/y     0   1   2

    0       0   1   2
    1       3   4   5
    2       6   7   8
    */
   node_id id1 = {1,1,3};
   REQUIRE(id1.i == 4);
   node_id id2 = {4, 3};
   REQUIRE(id2.xi == 1);
   REQUIRE(id2.yi == 1);
}

// Test case for the node_id struct
TEST_CASE("node_id Matrix Interface Test") {
    Grid g(3,3);
    node_id id(4,3);
    g.nodes[1][1].x = 2;

    REQUIRE(g[id]->x == 2);
    REQUIRE(g.nodes.data[id.i].x == 2);
}


TEST_CASE("Accessing Grid Elements") {
    Grid g(3, 3, 0);

    // Modify an node
    g.nodes[1][1].x = 5.0;
    g.nodes[1][1].y = 10.0;

    // Verify that the modification is reflected in the grid
    REQUIRE(g.nodes[1][1].x == 5.0);
    REQUIRE(g.nodes[1][1].y == 10.0);

    // Verify some other elements
    REQUIRE(g.nodes[0][0].x == 0.0);
    REQUIRE(g.nodes[2][2].y == 0.0);
}

// Test case for checking neighbors with periodic boundary conditions
TEST_CASE("Neighbors with Periodic Boundary Conditions") {
    Grid g(3, 3);

    /*
    Visualization of xi, yi and i in node_id for 3x3 matrix

    x/y     0   1   2

    0       0   1   2
    1       3   4   5
    2       6   7   8
    */

    // Corner node [0][0]
    std::array<node_id, 4> neighbors_corner = g.nodes[0][0].neighbours;
    REQUIRE(neighbors_corner.size() == 4);
    REQUIRE(neighbors_corner[0].xi == 0);
    REQUIRE(neighbors_corner[0].yi == 2); // Left (wrap around)
    REQUIRE(neighbors_corner[1].xi == 0);
    REQUIRE(neighbors_corner[1].yi == 1); // Right
    REQUIRE(neighbors_corner[2].xi == 2);
    REQUIRE(neighbors_corner[2].yi == 0); // Up (wrap around)
    REQUIRE(neighbors_corner[3].xi == 1);
    REQUIRE(neighbors_corner[3].yi == 0); // Down

        // Element at [2][2]
    std::array<node_id, 4> neighbors_22 = g.nodes[2][2].neighbours;
    REQUIRE(neighbors_22.size() == 4);
    REQUIRE(neighbors_22[0].xi == 2);
    REQUIRE(neighbors_22[0].yi == 1); // Left
    REQUIRE(neighbors_22[1].xi == 2);
    REQUIRE(neighbors_22[1].yi == 0); // Right (wrap around)
    REQUIRE(neighbors_22[2].xi == 1);
    REQUIRE(neighbors_22[2].yi == 2); // Up
    REQUIRE(neighbors_22[3].xi == 0);
    REQUIRE(neighbors_22[3].yi == 2); // Down (wrap around)
}

// Test case for setting node positions in a regular square grid
TEST_CASE("Setting Node Positions in a Regular Grid") {
    double spacing = 1.0; // Spacing between nodes
    Grid g(4, 4, spacing);


    // Verify that node positions are correctly set
    REQUIRE(g.nodes[0][0].x == 0.0);
    REQUIRE(g.nodes[0][0].y == 0.0);


    REQUIRE(g.nodes[2][2].x == 2 * spacing);
    REQUIRE(g.nodes[2][2].y == 2 * spacing);

    REQUIRE(g.nodes[2][3].x == 2 * spacing);
    REQUIRE(g.nodes[2][3].y == 3 * spacing);

    // You can add more checks as needed
}


TEST_CASE("Create Triangles Test") {
    Grid g(3, 3); // Create a grid with 3x3 dimensions

    /*
    0   1   2
    3   4   5
    6   7   8
    */
   // The triangles should be 013 134 124 245 346 467 457 and 578

    // Ensure the number of triangles created matches the expected count
    CHECK(g.triangles.size() == 2 * (g.nodes.rows - 1) * (g.nodes.cols - 1));

    // Check some specific triangles to ensure they were correctly created
    // Replace these with actual checks based on your grid layout
    CHECK(g.triangles[0].a1->id.i == 0); // Check the first triangle's first node
    CHECK(g.triangles[0].a2->id.i == 1); // Check the first triangle's second node
    CHECK(g.triangles[0].a3->id.i == 3); // Check the first triangle's third node

    CHECK(g.triangles[1].a1->id.i == 1); // Check the second triangle's first node
    CHECK(g.triangles[1].a2->id.i == 3); // Check the second triangle's second node
    CHECK(g.triangles[1].a3->id.i == 4); // Check the second triangle's third node
    
    CHECK(g.triangles[7].a1->id.i == 5); // Check the second triangle's first node
    CHECK(g.triangles[7].a2->id.i == 7); // Check the second triangle's second node
    CHECK(g.triangles[7].a3->id.i == 8); // Check the second triangle's third node

    // Add more checks as needed for your specific grid layout
}