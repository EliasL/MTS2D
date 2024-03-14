#include "run/doctest.h"
#include "../src/Mesh/mesh.h" // Include the header for your surface struct

TEST_CASE("Mesh Initialization")
{
    // Create a surface with known dimensions for testing
    Mesh mesh(4, 4);

    // Verify that the surface dimensions are correct
    CHECK(mesh.nodes.rows == 4);
    CHECK(mesh.nodes.cols == 4);
}

// Test case for the NodeId struct
TEST_CASE("NodeId Struct Test")
{
    /*
    Visualization of col, row and i in NodeId for 3x3 matrix

    2       6   7   8
    1       3   4   5
    0       0   1   2

    y/x     0   1   2
    */
    NodeId id1 = {1, 1, 3};
    CHECK(id1.i == 4);
    NodeId id2 = {4, 3};
    CHECK(id2.col == 1);
    CHECK(id2.row == 1);
}

// Test case for the NodeId struct
TEST_CASE("NodeId Matrix Interface Test")
{
    Mesh mesh(3, 3);
    NodeId id(1, 1, 3);
    mesh.nodes[1][1].setPos({2.0, 0.0});

    CHECK(mesh[id]->x() == 2);
    CHECK(mesh.nodes.data[id.i].x() == 2);
}

// Test case for the fixedNode bool
TEST_CASE("FixedNode Bool Test")
{
    Mesh mesh(3, 3, false);
    CHECK(mesh.nodes[0][0].fixedNode == true);
    CHECK(mesh.nodes[1][1].fixedNode == false);
}

TEST_CASE("Accessing Mesh Nodes")
{
    Mesh mesh(3, 3, 1.0);

    // Modify an Node
    mesh.nodes[1][1].setPos({5.0, 10.0});

    // Verify that the modification is reflected in the surface
    CHECK(mesh.nodes[1][1].x() == 5.0);
    CHECK(mesh.nodes[1][1].y() == 10.0);

    // Verify some other elements
    CHECK(mesh.nodes[0][0].x() == 0.0);
    CHECK(mesh.nodes[2][2].y() == 2.0);
}

// Test case for checking neighbors with periodic boundary conditions
TEST_CASE("Neighbors with Periodic Boundary Conditions")
{
    Mesh mesh(3, 3);

    /*
    Visualization of row, col and i in NodeId for 3x3 matrix

    2       6   7   8
    1       3   4   5
    0       0   1   2

    y/x     0   1   2
    */

    // Corner Node [0][0]
    std::array<NodeId, 4> neighbors_corner = mesh.nodes[0][0].neighbours;
    CHECK(neighbors_corner.size() == 4);
    CHECK(neighbors_corner[LEFT_N].col == 2);
    CHECK(neighbors_corner[LEFT_N].row == 0);
    CHECK(neighbors_corner[RIGHT_N].col == 1);
    CHECK(neighbors_corner[RIGHT_N].row == 0);
    CHECK(neighbors_corner[UP_N].col == 0);
    CHECK(neighbors_corner[UP_N].row == 1);
    CHECK(neighbors_corner[DOWN_N].col == 0);
    CHECK(neighbors_corner[DOWN_N].row == 2);

    // Element at [2][2]
    std::array<NodeId, 4> neighbors_22 = mesh.nodes[2][2].neighbours;
    CHECK(neighbors_22.size() == 4);
    CHECK(neighbors_22[LEFT_N].col == 1);
    CHECK(neighbors_22[LEFT_N].row == 2);
    CHECK(neighbors_22[RIGHT_N].col == 0);
    CHECK(neighbors_22[RIGHT_N].row == 2);
    CHECK(neighbors_22[UP_N].col == 2);
    CHECK(neighbors_22[UP_N].row == 0);
    CHECK(neighbors_22[DOWN_N].col == 2);
    CHECK(neighbors_22[DOWN_N].row == 1);
}

// Test case for setting Node positions in a regular square surface
TEST_CASE("Setting Node Positions in a Regular Mesh")
{
    double spacing = 1.0; // Spacing between nodes
    Mesh mesh(4, 4, spacing);

    // Verify that Node positions are correctly set
    CHECK(mesh.nodes[0][0].x() == 0.0);
    CHECK(mesh.nodes[0][0].y() == 0.0);

    CHECK(mesh.nodes[2][2].x() == 2 * spacing);
    CHECK(mesh.nodes[2][2].y() == 2 * spacing);

    CHECK(mesh.nodes[2][3].x() == 3 * spacing);
    CHECK(mesh.nodes[2][3].y() == 2 * spacing);

    // You can add more checks as needed
}

TEST_CASE("Create Elements Test")
{
    Mesh mesh(3, 3, false); // Create a surface with 3x3 dimensions

    /*
    6   7   8
    3   4   5
    0   1   2
    */
    // The elements should be 013 134 124 245 346 467 457 and 578

    // Ensure the number of elements created matches the expected count
    CHECK(mesh.elements.size() == 2 * (mesh.nodes.rows - 1) * (mesh.nodes.cols - 1));

    // Check some specific elements to ensure they were correctly created
    // Replace these with actual checks based on your surface layout
    CHECK(mesh.elements[0]->nodes[0].realId.i == 0); // Check the first Element's first Node
    CHECK(mesh.elements[0]->nodes[1].realId.i == 1); // Check the first Element's second Node
    CHECK(mesh.elements[0]->nodes[2].realId.i == 3); // Check the first Element's third Node

    CHECK(mesh.elements[1]->nodes[0].realId.i == 1); // Check the second Element's first Node
    CHECK(mesh.elements[1]->nodes[1].realId.i == 3); // Check the second Element's second Node
    CHECK(mesh.elements[1]->nodes[2].realId.i == 4); // Check the second Element's third Node

    CHECK(mesh.elements[7]->nodes[0].realId.i == 5); // Check the eigth Element's first Node
    CHECK(mesh.elements[7]->nodes[1].realId.i == 7); // Check the eigth Element's second Node
    CHECK(mesh.elements[7]->nodes[2].realId.i == 8); // Check the eigth Element's third Node

    mesh = Mesh(3, 3, true); // Create a surface with 3x3 dimensions and PBC

    /*
    6   7   8
    3   4   5
    0   1   2
    */

    // Ensure the number of elements created matches the expected count
    CHECK(mesh.elements.size() == 2 * (mesh.nodes.rows) * (mesh.nodes.cols));

    CHECK(mesh.elements[4]->nodes[0].realId.i == 2); // Check the fifth Element's first Node
    CHECK(mesh.elements[4]->nodes[1].realId.i == 0); // Check the fifth Element's second Node
    CHECK(mesh.elements[4]->nodes[2].realId.i == 5); // Check the fifth Element's third Node
}

TEST_CASE("Node transformation using transform function")
{
    // Create a node with an initial position
    Node originalNode;
    originalNode.setPos({1.0, 2.0});

    // Define a transformation matrix (e.g., a rotation matrix)
    Matrix2x2<double> matrix = {{0, -1}, {1, 0}};

    // Apply the transformation
    Node transformedNode = transform(matrix, originalNode);

    // Check the results
    CHECK(transformedNode.x() == -2.0);
    CHECK(transformedNode.y() == 1.0);
}

TEST_CASE("In-place Node transformation using transformInPlace function")
{
    // Create a node with an initial position
    Node node;
    node.setPos({1.0, 2.0});

    // Define a transformation matrix (e.g., a scale matrix)
    Matrix2x2<double> matrix = {{2, 0}, {0, 2}};

    // Apply the transformation in-place
    transformInPlace(matrix, node);

    // Check the results
    CHECK(node.x() == 2.0);
    CHECK(node.y() == 4.0);
}