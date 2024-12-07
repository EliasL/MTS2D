#include "../src/Mesh/mesh.h" // Include the header for your surface struct
#include "run/doctest.h"

TEST_CASE("Mesh Initialization") {
  // Create a surface with known dimensions for testing
  Mesh mesh(4, 4);

  // Verify that the surface dimensions are correct
  CHECK(mesh.rows == 4);
  CHECK(mesh.cols == 4);
}

// Test case for the NodeId struct
TEST_CASE("NodeId Matrix Interface Test") {
  Mesh mesh(3, 3);
  NodeId id(1, 1, 3);
  mesh.nodes(1, 1).setPos({2.0, 0.0});

  CHECK(mesh[id]->pos()[0] == 2);
  CHECK(mesh.nodes(id.i).pos()[0] == 2);
}

// Test case for the fixedNode bool
TEST_CASE("FixedNode Bool Test") {
  Mesh mesh(3, 3, false);
  mesh.fixBorderNodes();
  CHECK(mesh.nodes(0, 0).fixedNode == true);
  CHECK(mesh.nodes(1, 1).fixedNode == false);
}

TEST_CASE("Accessing Mesh Nodes") {
  Mesh mesh(3, 3, 1.0);

  // Modify an Node
  mesh.nodes(1, 1).setPos({5.0, 10.0});

  // Verify that the modification is reflected in the surface
  CHECK(mesh.nodes(1, 1).pos()[0] == 5.0);
  CHECK(mesh.nodes(1, 1).pos()[1] == 10.0);

  // Verify some other elements
  CHECK(mesh.nodes(0, 0).pos()[0] == 0.0);
  CHECK(mesh.nodes(2, 2).pos()[1] == 2.0);
}

// Test case for checking neighbors with periodic boundary conditions
TEST_CASE("Neighbors with Periodic Boundary Conditions") {
  Mesh mesh(3, 3);

  /*
  Visualization of row, col and i in NodeId for 3x3 matrix

  2       6   7   8
  1       3   4   5
  0       0   1   2

  y/x     0   1   2
  */

  // Corner Node (0,0)
  std::array<NodeId, 4> neighbors_corner = mesh.nodes(0, 0).neighbours;
  CHECK(neighbors_corner.size() == 4);
  CHECK(neighbors_corner[LEFT_N].col == 2);
  CHECK(neighbors_corner[LEFT_N].row == 0);
  CHECK(neighbors_corner[RIGHT_N].col == 1);
  CHECK(neighbors_corner[RIGHT_N].row == 0);
  CHECK(neighbors_corner[UP_N].col == 0);
  CHECK(neighbors_corner[UP_N].row == 1);
  CHECK(neighbors_corner[DOWN_N].col == 0);
  CHECK(neighbors_corner[DOWN_N].row == 2);

  // Element at (2,2)
  std::array<NodeId, 4> neighbors_22 = mesh.nodes(2, 2).neighbours;
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
TEST_CASE("Setting Node Positions in a Regular Mesh") {
  double spacing = 1.0; // Spacing between nodes
  Mesh mesh(4, 4, spacing);

  // Verify that Node positions are correctly set
  CHECK(mesh.nodes(0, 0).pos()[0] == 0.0);
  CHECK(mesh.nodes(0, 0).pos()[1] == 0.0);

  CHECK(mesh.nodes(2, 2).pos()[0] == 2 * spacing);
  CHECK(mesh.nodes(2, 2).pos()[1] == 2 * spacing);

  CHECK(mesh.nodes(2, 3).pos()[0] == 3 * spacing);
  CHECK(mesh.nodes(2, 3).pos()[1] == 2 * spacing);

  // You can add more checks as needed
}

TEST_CASE("Create Elements Test") {
  Mesh mesh(3, 3, false); // Create a surface with 3x3 dimensions
  CHECK(mesh.nrElements == 8);

  /*
  6   7   8
  3   4   5
  0   1   2
  */
  // The elements should be 013 134 124 245 346 467 457 and 578
  // Ensure the number of elements created matches the expected count
  CHECK(mesh.elements.size() == 2 * (mesh.rows - 1) * (mesh.cols - 1));

  // Check some specific elements to ensure they were correctly created
  // Replace these with actual checks based on your surface layout
  // Check the first Element's first Node
  CHECK(mesh.elements[0].nodes[0].id.i == 0);
  // Check the first Element's second Node
  CHECK(mesh.elements[0].nodes[1].id.i == 1);
  // Check the first Element's third Node
  CHECK(mesh.elements[0].nodes[2].id.i == 3);

  // Check the second Element's first Node
  CHECK(mesh.elements[1].nodes[0].id.i == 1);
  // Check the second Element's second Node
  CHECK(mesh.elements[1].nodes[1].id.i == 3);
  // Check the second Element's third Node
  CHECK(mesh.elements[1].nodes[2].id.i == 4);

  // Check the eigth Element's first Node
  CHECK(mesh.elements[7].nodes[0].id.i == 5);
  // Check the eigth Element's second Node
  CHECK(mesh.elements[7].nodes[1].id.i == 7);
  // Check the eigth Element's third Node
  CHECK(mesh.elements[7].nodes[2].id.i == 8);

  mesh = Mesh(3, 3, true); // Create a surface with 3x3 dimensions and PBC
  CHECK(mesh.nrElements == 18);
  /*
  6   7   8
  3   4   5
  0   1   2
  */
  // Ensure the number of elements created matches the expected count
  CHECK(mesh.elements.size() == 2 * (mesh.rows) * (mesh.cols));

  CHECK(mesh.elements[4].nodes[0].id.i ==
        2); // Check the fifth Element's first Node
  CHECK(mesh.elements[4].nodes[1].id.i ==
        0); // Check the fifth Element's second Node
  CHECK(mesh.elements[4].nodes[2].id.i ==
        5); // Check the fifth Element's third Node
}

TEST_CASE("Node transformation using transform function") {
  // Create a node with an initial position
  Node originalNode;
  originalNode.setPos({1.0, 2.0});

  // Define a transformation matrix (e.g., a rotation matrix)
  Eigen::Matrix2d matrix;
  matrix << 0, -1, 1, 0;

  // Apply the transformation
  Node transformedNode = transform(matrix, originalNode);

  // Check the results
  CHECK(transformedNode.pos()[0] == -2.0);
  CHECK(transformedNode.pos()[1] == 1.0);
}

TEST_CASE("In-place Node transformation using transformInPlace function") {
  // Create a node with an initial position
  Node node;
  node.setPos({1.0, 2.0});

  // Define a transformation matrix (e.g., a scale matrix)
  Eigen::Matrix2d matrix;
  matrix << 2, 0, 0, 2;

  // Apply the transformation in-place
  transformInPlace(matrix, node);

  // Check the results
  CHECK(node.pos()[0] == 2.0);
  CHECK(node.pos()[1] == 4.0);
}

TEST_CASE("Periodic node indexing") {
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

  for (int i = 0; i < mesh.nrElements; i++) {
    TElement &e = mesh.elements[i];
    for (size_t j = 0; j < e.nodes.size(); j++) {
      CHECK(e.nodes[j].ghostId.i == ans[i][j]);
    }
  }
}

TEST_CASE("Periodic node positions") {
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

  for (int i = 0; i < mesh.nrElements; i++) {
    TElement &e = mesh.elements[i];
    for (size_t j = 0; j < e.nodes.size(); j++) {
      CHECK(e.nodes[j].ghostId.i == ans[i][j]);
    }
  }
}

TEST_CASE("Periodic node positions after transformation") {
  // Create a node with an initial position
  Node node;
  node.setPos({1.0, 2.0});

  // Define a transformation matrix (e.g., a scale matrix)
  Eigen::Matrix2d matrix;
  matrix << 2, 0, 0, 2;

  // Apply the transformation in-place
  transformInPlace(matrix, node);

  // Check the results
  CHECK(node.pos()[0] == 2.0);
  CHECK(node.pos()[1] == 4.0);
}

TEST_CASE("Each node has six elements") {

  Mesh mesh(2, 2, true);
  // Method 1: Using size() and data()
  for (int i = 0; i < mesh.nodes.size(); ++i) {
    Node n = mesh.nodes.data()[i];
    CHECK(n.elementIndices.size() == 6);
  }
}