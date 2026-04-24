#include "../src/Mesh/mesh.h" // Include the header for your surface struct
#include "run/doctest.h"
#include <cmath>
#include <iostream>

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
  Node corner = mesh.nodes(0, 0);
  CHECK(mesh.getNeighbourNode(corner, LEFT_N)->id.col() == 2);
  CHECK(mesh.getNeighbourNode(corner, LEFT_N)->id.row() == 0);
  CHECK(mesh.getNeighbourNode(corner, RIGHT_N)->id.col() == 1);
  CHECK(mesh.getNeighbourNode(corner, RIGHT_N)->id.row() == 0);
  CHECK(mesh.getNeighbourNode(corner, UP_N)->id.col() == 0);
  CHECK(mesh.getNeighbourNode(corner, UP_N)->id.row() == 1);
  CHECK(mesh.getNeighbourNode(corner, DOWN_N)->id.col() == 0);
  CHECK(mesh.getNeighbourNode(corner, DOWN_N)->id.row() == 2);

  // Element at (2,2)
  corner = mesh.nodes(2, 2);
  CHECK(mesh.getNeighbourNode(corner, LEFT_N)->id.col() == 1);
  CHECK(mesh.getNeighbourNode(corner, LEFT_N)->id.row() == 2);
  CHECK(mesh.getNeighbourNode(corner, RIGHT_N)->id.col() == 0);
  CHECK(mesh.getNeighbourNode(corner, RIGHT_N)->id.row() == 2);
  CHECK(mesh.getNeighbourNode(corner, UP_N)->id.col() == 2);
  CHECK(mesh.getNeighbourNode(corner, UP_N)->id.row() == 0);
  CHECK(mesh.getNeighbourNode(corner, DOWN_N)->id.col() == 2);
  CHECK(mesh.getNeighbourNode(corner, DOWN_N)->id.row() == 1);
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

TEST_CASE("Setting Node Positions in a Triangular Mesh") {
  const double spacing = 2.0;
  Mesh mesh(3, 3, spacing, 0, true, "major", "contiTriangular");

  const double h = std::sqrt(3.0) * 0.5 * spacing;

  CHECK(mesh.latticeBasis(0, 0) == doctest::Approx(spacing));
  CHECK(mesh.latticeBasis(0, 1) == doctest::Approx(0.5 * spacing));
  CHECK(mesh.latticeBasis(1, 0) == doctest::Approx(0.0));
  CHECK(mesh.latticeBasis(1, 1) == doctest::Approx(h));

  // Row 0
  CHECK(mesh.nodes(0, 0).pos()[0] == doctest::Approx(0.0));
  CHECK(mesh.nodes(0, 0).pos()[1] == doctest::Approx(0.0));
  CHECK(mesh.nodes(0, 1).pos()[0] == doctest::Approx(spacing));
  CHECK(mesh.nodes(0, 1).pos()[1] == doctest::Approx(0.0));

  // Row 1 (shifted by spacing/2)
  CHECK(mesh.nodes(1, 0).pos()[0] == doctest::Approx(0.5 * spacing));
  CHECK(mesh.nodes(1, 0).pos()[1] == doctest::Approx(h));
  CHECK(mesh.nodes(1, 1).pos()[0] == doctest::Approx(1.5 * spacing));
  CHECK(mesh.nodes(1, 1).pos()[1] == doctest::Approx(h));
}

TEST_CASE("Create Elements Test") {
  Mesh mesh(3, 3, false); // Create a surface with 3x3 dimensions
  CHECK(mesh.nrElements == 8);

  /*
  6   7   8
  3   4   5
  0   1   2
  */
  // The elements should be 013 413 124 524 346 746 457 and 857
  // Always starting with the corner node (such that the angle is 90), and
  // then folloing increasing node index number
  // Ensure the number of elements created matches the expected count
  CHECK(mesh.elements.size() == 2 * (mesh.rows - 1) * (mesh.cols - 1));

  // Check some specific elements to ensure they were correctly created
  // Replace these with actual checks based on your surface layout
  // Check the first Element's first Node
  CHECK(mesh.elements[0].ghostNodes[0].referenceId.i == 0);
  // Check the first Element's second Node
  CHECK(mesh.elements[0].ghostNodes[1].referenceId.i == 1);
  // Check the first Element's third Node
  CHECK(mesh.elements[0].ghostNodes[2].referenceId.i == 3);

  // Check the second Element's first Node
  CHECK(mesh.elements[1].ghostNodes[0].referenceId.i == 4);
  // Check the second Element's second Node
  CHECK(mesh.elements[1].ghostNodes[1].referenceId.i == 1);
  // Check the second Element's third Node
  CHECK(mesh.elements[1].ghostNodes[2].referenceId.i == 3);

  // Check the eigth Element's first Node
  CHECK(mesh.elements[7].ghostNodes[0].referenceId.i == 8);
  // Check the eigth Element's second Node
  CHECK(mesh.elements[7].ghostNodes[1].referenceId.i == 5);
  // Check the eigth Element's third Node
  CHECK(mesh.elements[7].ghostNodes[2].referenceId.i == 7);

  mesh = Mesh(3, 3, true); // Create a surface with 3x3 dimensions and PBC
  CHECK(mesh.nrElements == 18);
  /*
  6   7   8
  3   4   5
  0   1   2
  */
  // Ensure the number of elements created matches the expected count
  CHECK(mesh.elements.size() == 2 * (mesh.rows) * (mesh.cols));

  CHECK(mesh.elements[4].ghostNodes[0].referenceId.i ==
        2); // Check the fifth Element's first Node
  CHECK(mesh.elements[4].ghostNodes[1].referenceId.i ==
        0); // Check the fifth Element's second Node
  CHECK(mesh.elements[4].ghostNodes[2].referenceId.i ==
        5); // Check the fifth Element's third Node
}

TEST_CASE("Create Elements Test PBC") {
  Mesh mesh(3, 3, true); // Create a surface with 3x3 dimensions
  CHECK(mesh.nrElements == 18);

  // Check that each node is connected to six elements
  for (int i = 0; i < mesh.nodes.size(); ++i) {
    CHECK(mesh.nodes(i).elementCount == 6);
  }
  Mesh mesh2(3, 3, 1, 0, true, "minor"); // Create a surface with 3x3 dimensions
  CHECK(mesh2.nrElements == 18);

  // Check that each node is connected to six elements
  for (int i = 0; i < mesh2.nodes.size(); ++i) {
    CHECK(mesh2.nodes(i).elementCount == 6);
  }
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

/*
Here we used a different ordering of the nodes inside the elements
In the new version, we always make sure to have the corner node in the corner of
the element pair. In otherwords, in the reference state, we choose vectors such
that all the representative angles in each element is the one with 90 degrees.
TEST_CASE("Old Periodic node indexing") {
  Mesh mesh(2, 2, 1, 0, true, "major");

  // Here are the real nodes
  // 2 3
  // 0 1
  // with elements 012 123 103 032 230 301 321 210

  // In the periodic mesh, we should have the nodes
  // 6 7 8
  // 3 4 5
  // 0 1 2
  // with elements 013 134 124 245 346 467 457 578

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
    for (size_t j = 0; j < e.ghostNodes.size(); j++) {
      CHECK(e.ghostNodes[j].ghostId.i == ans[i][j]);
    }
  }
}
*/

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

  Mesh mesh(2, 2, 1, 0, true, "major");
  // Method 1: Using size() and data()
  for (int i = 0; i < mesh.nodes.size(); ++i) {
    Node n = mesh.nodes.data()[i];
    CHECK(n.elementCount == 6);
  }
}

TEST_CASE("Force check") {
  // TODO
  Mesh mesh(3, 3, 1, 0, false, "major");

  // Offset the center node
  // std::cout << "Starting force check\n" << std::endl;
  omp_set_num_threads(1);

  mesh.nodes(1, 1).setDisplacement({0.2, 0});
  mesh.markDirty();
  mesh.updateElements();
  mesh.applyForceFromElementsToNodes();
  // std::cout << mesh << std::endl;
  mesh.setSimNameAndDataPath("defaultName", "test_data");
  mesh.writeToVtu("reconnectTest");
  // std::cout << "After force check\n" << std::endl;
}
