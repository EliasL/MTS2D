#include "run/doctest.h"
#include "../src/Mesh/mesh.h" // Include the header for your surface struct

TEST_CASE("Mesh Initialization")
{
    // Create a surface with known dimensions for testing
    Mesh mesh(4, 4);

    // Verify that the surface dimensions are correct
    REQUIRE(mesh.nodes.rows == 4);
    REQUIRE(mesh.nodes.cols == 4);
}

// Test case for the NodeId struct
TEST_CASE("NodeId Struct Test")
{
    /*
    Visualization of col, row and i in NodeId for 3x3 matrix

    x/y     0   1   2

    0       0   1   2
    1       3   4   5
    2       6   7   8
    */
    NodeId id1 = {1, 1, 3};
    REQUIRE(id1.i == 4);
    NodeId id2 = {4, 3};
    REQUIRE(id2.col == 1);
    REQUIRE(id2.row == 1);
}

// Test case for the NodeId struct
TEST_CASE("NodeId Matrix Interface Test")
{
    Mesh mesh(3, 3);
    NodeId id(4, 3);
    mesh.nodes[1][1].setPos(2.0, 0.0);

    REQUIRE(mesh[id]->X() == 2);
    REQUIRE(mesh.nodes.data[id.i].X() == 2);
}

// Test case for the fixedNode bool
TEST_CASE("FixedNode Bool Test")
{
    Mesh mesh(3, 3);
    REQUIRE(mesh.nodes[0][0].fixedNode == true);
    REQUIRE(mesh.nodes[1][1].fixedNode == false);
}

TEST_CASE("Accessing Mesh Nodes")
{
    Mesh mesh(3, 3, 1.0);

    // Modify an Node
    mesh.nodes[1][1].setPos(5.0, 10.0);

    // Verify that the modification is reflected in the surface
    REQUIRE(mesh.nodes[1][1].X() == 5.0);
    REQUIRE(mesh.nodes[1][1].Y() == 10.0);

    // Verify some other elements
    REQUIRE(mesh.nodes[0][0].X() == 0.0);
    REQUIRE(mesh.nodes[2][2].Y() == 2.0);
}

// Test case for checking neighbors with periodic boundary conditions
TEST_CASE("Neighbors with Periodic Boundary Conditions")
{
    Mesh mesh(3, 3);

    /*
    Visualization of row, col and i in NodeId for 3x3 matrix

    y\x     0   1   2

    0       0   1   2
    1       3   4   5
    2       6   7   8
    */

    // Corner Node [0][0]
    std::array<NodeId, 4> neighbors_corner = mesh.nodes[0][0].neighbours;
    REQUIRE(neighbors_corner.size() == 4);
    REQUIRE(neighbors_corner[0].row == 2);
    REQUIRE(neighbors_corner[0].col == 0); // Left (wrap around)
    REQUIRE(neighbors_corner[1].row == 1);
    REQUIRE(neighbors_corner[1].col == 0); // Right
    REQUIRE(neighbors_corner[2].row == 0);
    REQUIRE(neighbors_corner[2].col == 2); // Up (wrap around)
    REQUIRE(neighbors_corner[3].row == 0);
    REQUIRE(neighbors_corner[3].col == 1); // Down

    // Element at [2][2]
    std::array<NodeId, 4> neighbors_22 = mesh.nodes[2][2].neighbours;
    REQUIRE(neighbors_22.size() == 4);
    REQUIRE(neighbors_22[0].row == 1);
    REQUIRE(neighbors_22[0].col == 2); // Left
    REQUIRE(neighbors_22[1].row == 0);
    REQUIRE(neighbors_22[1].col == 2); // Right (wrap around)
    REQUIRE(neighbors_22[2].row == 2);
    REQUIRE(neighbors_22[2].col == 1); // Up
    REQUIRE(neighbors_22[3].row == 2);
    REQUIRE(neighbors_22[3].col == 0); // Down (wrap around)
}

// Test case for setting Node positions in a regular square surface
TEST_CASE("Setting Node Positions in a Regular Mesh")
{
    double spacing = 1.0; // Spacing between nodes
    Mesh mesh(4, 4, spacing);

    // Verify that Node positions are correctly set
    REQUIRE(mesh.nodes[0][0].X() == 0.0);
    REQUIRE(mesh.nodes[0][0].Y() == 0.0);

    REQUIRE(mesh.nodes[2][2].X() == 2 * spacing);
    REQUIRE(mesh.nodes[2][2].Y() == 2 * spacing);

    REQUIRE(mesh.nodes[2][3].X() == 3 * spacing);
    REQUIRE(mesh.nodes[2][3].Y() == 2 * spacing);

    // You can add more checks as needed
}

TEST_CASE("Create Elements Test")
{
    Mesh mesh(3, 3); // Create a surface with 3x3 dimensions

    /*
    0   1   2
    3   4   5
    6   7   8
    */
    // The elements should be 013 134 124 245 346 467 457 and 578

    // Ensure the number of elements created matches the expected count
    CHECK(mesh.elements.size() == 2 * (mesh.nodes.rows - 1) * (mesh.nodes.cols - 1));

    // Check some specific elements to ensure they were correctly created
    // Replace these with actual checks based on your surface layout
    CHECK(mesh.elements[0].n1->id.i == 0); // Check the first Element's first Node
    CHECK(mesh.elements[0].n2->id.i == 1); // Check the first Element's second Node
    CHECK(mesh.elements[0].n3->id.i == 3); // Check the first Element's third Node

    CHECK(mesh.elements[1].n1->id.i == 1); // Check the second Element's first Node
    CHECK(mesh.elements[1].n2->id.i == 3); // Check the second Element's second Node
    CHECK(mesh.elements[1].n3->id.i == 4); // Check the second Element's third Node

    CHECK(mesh.elements[7].n1->id.i == 5); // Check the second Element's first Node
    CHECK(mesh.elements[7].n2->id.i == 7); // Check the second Element's second Node
    CHECK(mesh.elements[7].n3->id.i == 8); // Check the second Element's third Node
}

TEST_CASE("Node transformation using transform function")
{
    // Create a node with an initial position
    Node originalNode;
    originalNode.setPos(1.0, 2.0);

    // Define a transformation matrix (e.g., a rotation matrix)
    Matrix2x2<double> matrix = {{0, -1}, {1, 0}};

    // Apply the transformation
    Node transformedNode = transform(matrix, originalNode);

    // Check the results
    REQUIRE(transformedNode.X() == -2.0);
    REQUIRE(transformedNode.Y() == 1.0);
}

TEST_CASE("In-place Node transformation using transformInPlace function")
{
    // Create a node with an initial position
    Node node;
    node.setPos(1.0, 2.0);

    // Define a transformation matrix (e.g., a scale matrix)
    Matrix2x2<double> matrix = {{2, 0}, {0, 2}};

    // Apply the transformation in-place
    transformInPlace(matrix, node);

    // Check the results
    REQUIRE(node.X() == 2.0);
    REQUIRE(node.Y() == 4.0);
}