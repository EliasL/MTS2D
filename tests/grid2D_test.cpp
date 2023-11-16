#include "run/doctest.h"
#include "../src/Mesh/mesh.h"  // Include the header for your surface struct

TEST_CASE("Mesh Initialization") {
    // Create a surface with known dimensions for testing
    Mesh mesh(4, 4);

    // Verify that the surface dimensions are correct
    // Verify that the surface dimensions are correct
    REQUIRE(mesh.nodes.rows == 4);
    REQUIRE(mesh.nodes.cols == 4);
}

// Test case for the NodeId struct
TEST_CASE("NodeId Struct Test") {
    /*
    Visualization of xi, yi and i in NodeId for 3x3 matrix

    x/y     0   1   2

    0       0   1   2
    1       3   4   5
    2       6   7   8
    */
   NodeId id1 = {1,1,3};
   REQUIRE(id1.i == 4);
   NodeId id2 = {4, 3};
   REQUIRE(id2.xi == 1);
   REQUIRE(id2.yi == 1);
}

// Test case for the NodeId struct
TEST_CASE("NodeId Matrix Interface Test") {
    Mesh mesh(3,3);
    NodeId id(4,3);
    mesh.nodes[1][1].x = 2;
    
    REQUIRE(mesh[id]->x == 2);
    REQUIRE(mesh.nodes.data[id.i].x == 2);
}


TEST_CASE("Accessing Mesh Elements") {
    Mesh mesh(3, 3, 0);

    // Modify an Node
    mesh.nodes[1][1].x = 5.0;
    mesh.nodes[1][1].y = 10.0;

    // Verify that the modification is reflected in the surface
    REQUIRE(mesh.nodes[1][1].x == 5.0);
    REQUIRE(mesh.nodes[1][1].y == 10.0);

    // Verify some other elements
    REQUIRE(mesh.nodes[0][0].x == 0.0);
    REQUIRE(mesh.nodes[2][2].y == 0.0);
}

// Test case for checking neighbors with periodic boundary conditions
TEST_CASE("Neighbors with Periodic Boundary Conditions") {
    Mesh mesh(3, 3);

    /*
    Visualization of xi, yi and i in NodeId for 3x3 matrix

    x/y     0   1   2

    0       0   1   2
    1       3   4   5
    2       6   7   8
    */

    // Corner Node [0][0]
    std::array<NodeId, 4> neighbors_corner = mesh.nodes[0][0].neighbours;
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
    std::array<NodeId, 4> neighbors_22 = mesh.nodes[2][2].neighbours;
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

// Test case for setting Node positions in a regular square surface
TEST_CASE("Setting Node Positions in a Regular Mesh") {
    double spacing = 1.0; // Spacing between nodes
    Mesh mesh(4, 4, spacing);


    // Verify that Node positions are correctly set
    REQUIRE(mesh.nodes[0][0].x == 0.0);
    REQUIRE(mesh.nodes[0][0].y == 0.0);


    REQUIRE(mesh.nodes[2][2].x == 2 * spacing);
    REQUIRE(mesh.nodes[2][2].y == 2 * spacing);

    REQUIRE(mesh.nodes[2][3].x == 2 * spacing);
    REQUIRE(mesh.nodes[2][3].y == 3 * spacing);

    // You can add more checks as needed
}


TEST_CASE("Create Triangles Test") {
    Mesh mesh(3, 3); // Create a surface with 3x3 dimensions

    /*
    0   1   2
    3   4   5
    6   7   8
    */
   // The triangles should be 013 134 124 245 346 467 457 and 578

    // Ensure the number of triangles created matches the expected count
    CHECK(mesh.triangles.size() == 2 * (mesh.nodes.rows - 1) * (mesh.nodes.cols - 1));

    // Check some specific triangles to ensure they were correctly created
    // Replace these with actual checks based on your surface layout
    CHECK(mesh.triangles[0].a1->id.i == 0); // Check the first Triangle's first Node
    CHECK(mesh.triangles[0].a2->id.i == 1); // Check the first Triangle's second Node
    CHECK(mesh.triangles[0].a3->id.i == 3); // Check the first Triangle's third Node

    CHECK(mesh.triangles[1].a1->id.i == 1); // Check the second Triangle's first Node
    CHECK(mesh.triangles[1].a2->id.i == 3); // Check the second Triangle's second Node
    CHECK(mesh.triangles[1].a3->id.i == 4); // Check the second Triangle's third Node
    
    CHECK(mesh.triangles[7].a1->id.i == 5); // Check the second Triangle's first Node
    CHECK(mesh.triangles[7].a2->id.i == 7); // Check the second Triangle's second Node
    CHECK(mesh.triangles[7].a3->id.i == 8); // Check the second Triangle's third Node

    // Add more checks as needed for your specific surface layout
}