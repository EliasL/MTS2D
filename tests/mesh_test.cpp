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
    CHECK(mesh.elements[0].id.realNodes[0].i == 0); // Check the first Element's first Node
    CHECK(mesh.elements[0].id.realNodes[1].i == 1); // Check the first Element's second Node
    CHECK(mesh.elements[0].id.realNodes[2].i == 3); // Check the first Element's third Node

    CHECK(mesh.elements[1].id.realNodes[0].i == 1); // Check the second Element's first Node
    CHECK(mesh.elements[1].id.realNodes[1].i == 3); // Check the second Element's second Node
    CHECK(mesh.elements[1].id.realNodes[2].i == 4); // Check the second Element's third Node

    CHECK(mesh.elements[7].id.realNodes[0].i == 5); // Check the eigth Element's first Node
    CHECK(mesh.elements[7].id.realNodes[1].i == 7); // Check the eigth Element's second Node
    CHECK(mesh.elements[7].id.realNodes[2].i == 8); // Check the eigth Element's third Node

    mesh = Mesh(3, 3, true); // Create a surface with 3x3 dimensions and PBC

    /*
    6   7   8
    3   4   5
    0   1   2
    */

    // Ensure the number of elements created matches the expected count
    CHECK(mesh.elements.size() == 2 * (mesh.nodes.rows) * (mesh.nodes.cols));

    CHECK(mesh.elements[4].id.realNodes[0].i == 2); // Check the fifth Element's first Node
    CHECK(mesh.elements[4].id.realNodes[1].i == 0); // Check the fifth Element's second Node
    CHECK(mesh.elements[4].id.realNodes[2].i == 5); // Check the fifth Element's third Node
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

TEST_CASE("Periodic node indexing")
{
    Mesh mesh(2, 2, true);
    /*
    Here are the real nodes
    2 3
    0 1
    with elements 012 123 103 032 230 301 321 210

    In the periodic mesh, we should have the nodes
    6 7 8
    3 4 5
    0 1 2
    with elements 013 134 124 245 346 467 457 578
    */
    int ans[8][3] = {
        {0, 1, 3}, // element 0: nodes 013
        {1, 3, 4}, // element 1: nodes 134
        {1, 2, 4}, // element 2: nodes 124
        {2, 4, 5}, // element 3: nodes 245
        {3, 4, 6}, // element 4: nodes 346
        {4, 6, 7}, // element 5: nodes 467
        {4, 5, 7}, // element 6: nodes 457
        {5, 7, 8}  // element 7: nodes 578
    };

    mesh.printConnectivity(false);
    for (size_t i = 0; i < mesh.nrElements; i++)
    {
        TElement &e = mesh.elements[i];
        for (size_t j = 0; j < e.id.realNodes.size(); j++)
        {
            CHECK(e.id.periodicNodes[j].i == ans[i][j]);
        }
    }
}

TEST_CASE("Periodic Mesh to Non-periodic mesh duplication")
{
    Mesh p_mesh(2, 2, true);

    p_mesh.nodes[0][0].setPos({0.1, 0.1});
    p_mesh.nodes[1][0].setPos({1.2, 0});
    p_mesh.nodes[0][1].setPos({0, 1.3});

    Mesh n_mesh = p_mesh.duplicateAsFixedBoundary();

    CHECK(n_mesh.nodes[2][2].x() == 2.1);
    CHECK(n_mesh.nodes[1][2].x() == 1.2);
    CHECK(n_mesh.nodes[1][2].y() == 2);
    CHECK(n_mesh.nodes[2][1].x() == 2);
    CHECK(n_mesh.nodes[2][1].y() == 1.3);

    /*
    Here are the real nodes
    2 3
    0 1
    with elements 012 123 103 032 230 301 321 210

    In the periodic mesh, we should have the nodes
    6 7 8
    3 4 5
    0 1 2
    with elements 013 134 124 245 346 467 457 578
    */
    int ans[8][3] = {
        {0, 1, 3}, // element 0: nodes 013
        {1, 3, 4}, // element 1: nodes 134
        {1, 2, 4}, // element 2: nodes 124
        {2, 4, 5}, // element 3: nodes 245
        {3, 4, 6}, // element 4: nodes 346
        {4, 6, 7}, // element 5: nodes 467
        {4, 5, 7}, // element 6: nodes 457
        {5, 7, 8}  // element 7: nodes 578
    };
    for (size_t i = 0; i < p_mesh.nrElements; i++)
    {
        TElement &e = p_mesh.elements[i];
        for (size_t j = 0; j < e.id.realNodes.size(); j++)
        {
            CHECK(n_mesh.elements[i].id.realNodes[j].i == ans[i][j]);
        }
    }
}

TEST_CASE("Simple periodic boundary transformation")
{
    /*
    We check that when accessing the position of nodes of elements that wrap
    around the system, the distances are correctly calculated using the
    currentDeformation. The effectively means that in the periodic mesh
    of a 2x2 system where we have the nodes
    6 7 8
    3 4 5
    0 1 2
    with elements 013 134 124 245 346 467 457 578,
    The nodes 2, 5, 6, 7 and 8 should effectively appear transformed.
    */

    Mesh mesh(2, 2, true);
    Matrix2x2<double> shear = {{1, 0.5},
                               {0, 1}};

    mesh.applyTransformation(shear);

    Mesh n_mesh = mesh.duplicateAsFixedBoundary();

    // The distance between node 4 and 5 should 1.
    // The position of node 5 should be (2.5,1)
    Node n4 = n_mesh.nodes.data[4];
    Node n5 = n_mesh.nodes.data[5];
    double xpos = n5.pos()[0];
    double distance = n5.pos()[0] - n4.pos()[0];
    CHECK(xpos == 2.5);
    CHECK(distance == 1);

    // The x position of node 7 and 8 should be 2 and 3 respectively
    Node n7 = n_mesh.nodes.data[7];
    Node n8 = n_mesh.nodes.data[8];
    CHECK(n7.pos()[0] == 2);
    CHECK(n8.pos()[0] == 3);
}

TEST_CASE("Simple periodic boundary transformation with load")
{
    /*
    We check that when accessing the position of nodes of elements that wrap
    around the system, the distances are correctly calculated using the
    currentDeformation. The effectively means that in the periodic mesh
    of a 2x2 system where we have the nodes
    6 7 8
    3 4 5
    0 1 2
    with elements 013 134 124 245 346 467 457 578,
    The nodes 2, 5, 6, 7 and 8 should effectively appear transformed.
    */

    Mesh mesh(2, 2, true);
    Matrix2x2<double> shear = {{1, 0.5},
                               {0, 1}};

    mesh.applyTransformationToSystemDeformation(shear);
    Mesh n_mesh = mesh.duplicateAsFixedBoundary();
    // The distance between node 4 and 5 should now be 1.5 instead of 1.
    // The position of node 5 should be (2.5,1)
    Node n4 = n_mesh.nodes.data[4];
    Node n5 = n_mesh.nodes.data[5];

    double xpos = n5.pos()[0];
    double distance = n5.pos()[0] - n4.pos()[0];
    CHECK(xpos == 2.5);
    CHECK(distance == 1.5);

    // Similar for element 346
    Node n3 = n_mesh.nodes.data[3];
    Node n6 = n_mesh.nodes.data[6];
    xpos = n6.pos()[0];
    distance = n6.pos()[0] - n3.pos()[0];
    CHECK(xpos == 1);
    CHECK(distance == 1);

    // The x position of node 7 and 8 should be 2 and 3 respectively
    Node n7 = n_mesh.nodes.data[7];
    Node n8 = n_mesh.nodes.data[8];
    CHECK(n7.pos()[0] == 2);
    CHECK(n8.pos()[0] == 3);
}

TEST_CASE("Compound periodic boundary transformation")
{
    /*
    Similar to the simple test, but now we also have an
    existing transformation applied to our system.

    6 7 8
    3 4 5
    0 1 2
    with elements 013 134 124 245 346 467 457 578.
    */

    Mesh mesh(2, 2, true);
    Matrix2x2<double> shear = {{1, 0.5},
                               {0, 1}};
    mesh.applyTransformation(shear);
    mesh.applyTransformationToSystemDeformation(shear);
    Mesh n_mesh = mesh.duplicateAsFixedBoundary();

    // The distance between node 4 and 5 should now be 1.5 instead of 1.
    // The position of node 5 should be (2.5,1)
    Node n4 = n_mesh.nodes.data[4];
    Node n5 = n_mesh.nodes.data[5];
    double xpos = n5.pos()[0];
    double distance = n5.pos()[0] - n4.pos()[0];
    CHECK(xpos == 3);
    CHECK(distance == 1.5);

    // The x position of node 7 and 8 should be 3 and 4 respectively
    Node n7 = n_mesh.nodes.data[7];
    Node n8 = n_mesh.nodes.data[8];

    CHECK(n7.pos()[0] == 3);
    CHECK(n8.pos()[0] == 4);
}