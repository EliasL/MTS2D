#include "../src/Mesh/node.h" // Include the header for your surface struct
#include "run/doctest.h"

// Test case for the NodeId struct
TEST_CASE("NodeId Struct Test") {
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

TEST_CASE("Node test") {
  /*
  This is a node at the corner of a 2x2 mesh
  2 3
  0 1

  so it should have index i=3
  */
  Node n(1, 1, 1, 2);
  CHECK(n.id.i == 3);
  CHECK(n.id.col == 1);
  CHECK(n.id.row == 1);
  CHECK(n.pos()[0] == 1);
  CHECK(n.pos()[1] == 1);
  CHECK(n.init_pos()[0] == 1);
  CHECK(n.init_pos()[1] == 1);
}

TEST_CASE("Ghost node construction and properties") {
  // Create a reference node at position (1,1) in a 2x2 mesh
  Node referenceNode(1, 1, 1, 2);

  SUBCASE("Basic ghost node construction") {
    // Create a ghost node one unit to the right (col+1)
    double a = 1.0;                              // Unit spacing
    Matrix2d deformation = Matrix2d::Identity(); // No deformation
    GhostNode ghost(referenceNode, 1, 2, 2, a, deformation);

    // Check IDs
    CHECK(ghost.referenceId.i == referenceNode.id.i);
    CHECK(ghost.ghostId.row == 1);
    CHECK(ghost.ghostId.col == 2);
    CHECK(ghost.ghostId.cols == 3); // original cols (2) + 1

    // Check position calculations
    CHECK(ghost.periodShift[0] == 1.0); // (2-1)*1.0
    CHECK(ghost.periodShift[1] == 0.0); // (1-1)*1.0

    // Check positions
    CHECK(ghost.pos[0] == 2.0); // 1.0 + 1.0
    CHECK(ghost.pos[1] == 1.0); // 1.0 + 0.0
    CHECK(ghost.init_pos[0] == 2.0);
    CHECK(ghost.init_pos[1] == 1.0);

    // Check displacement matches reference node
    CHECK(ghost.u[0] == referenceNode.u()[0]);
    CHECK(ghost.u[1] == referenceNode.u()[1]);
  }

  SUBCASE("Ghost node with deformation") {
    double a = 1.0;
    Matrix2d deformation = Matrix2d::Zero();
    deformation(0, 0) = 2.0; // Stretch in x direction
    deformation(1, 1) = 0.5; // Compress in y direction

    // Create a ghost node one unit diagonally (row+1, col+1)
    GhostNode ghost(referenceNode, 2, 2, 2, a, deformation);

    // Check periodShift
    CHECK(ghost.periodShift[0] == 1.0); // (2-1)*1.0
    CHECK(ghost.periodShift[1] == 1.0); // (2-1)*1.0

    // Check deformed position
    // pos = referencePos + deformation * periodShift
    CHECK(ghost.pos[0] == 3.0); // 1.0 + (2.0*1.0 + 0.0*1.0)
    CHECK(ghost.pos[1] == 1.5); // 1.0 + (0.0*1.0 + 0.5*1.0)

    // Check initial position (undeformed)
    CHECK(ghost.init_pos[0] == 2.0); // 1.0 + 1.0
    CHECK(ghost.init_pos[1] == 2.0); // 1.0 + 1.0
  }

  SUBCASE("Ghost node with negative offset") {
    double a = 1.0;
    Matrix2d deformation = Matrix2d::Identity();

    // Create a ghost node one unit to the left (col-1)
    GhostNode ghost(referenceNode, 1, 0, 2, a, deformation);

    // Check periodShift
    CHECK(ghost.periodShift[0] == -1.0); // (0-1)*1.0
    CHECK(ghost.periodShift[1] == 0.0);  // (1-1)*1.0

    // Check positions
    CHECK(ghost.pos[0] == 0.0); // 1.0 + (-1.0)
    CHECK(ghost.pos[1] == 1.0); // 1.0 + 0.0
  }

  SUBCASE("Convenience constructor test") {
    double a = 1.0;
    Matrix2d deformation = Matrix2d::Identity();

    // This should create a ghost node at the same position as the reference
    // node
    GhostNode ghost(referenceNode, a, deformation);

    // Check IDs
    CHECK(ghost.ghostId.row == referenceNode.id.row);
    CHECK(ghost.ghostId.col == referenceNode.id.col);

    // Check periodShift is zero
    CHECK(ghost.periodShift[0] == 0.0);
    CHECK(ghost.periodShift[1] == 0.0);

    // Check positions match reference node
    CHECK(ghost.pos[0] == referenceNode.pos()[0]);
    CHECK(ghost.pos[1] == referenceNode.pos()[1]);
    CHECK(ghost.init_pos[0] == referenceNode.init_pos()[0]);
    CHECK(ghost.init_pos[1] == referenceNode.init_pos()[1]);
  }
}

TEST_CASE("Periodic Ghost Node with Simple Shear") {
  // Create a reference node at position (0,1) in a 2x2 mesh
  // This node is at the left edge (x=0), one unit up (y=1)
  Node referenceNode(1, 1, 0, 2);

  // Initial setup with identity deformation (no deformation)
  double a = 1.0;
  Matrix2d deformation = Matrix2d::Identity();

  // Create a ghost node that is "across the system"
  // We place it 2 units to the right of the reference
  // So the periodic shift should be (2,0)
  GhostNode ghost(referenceNode, 1, 2, 2, a, deformation);

  // Check the periodic shift is correct (2 units in x, 0 in y)
  CHECK(ghost.periodShift[0] == doctest::Approx(2.0)); // (2-0)*1.0
  CHECK(ghost.periodShift[1] == doctest::Approx(0.0)); // (1-1)*1.0

  // Check initial positions without deformation
  CHECK(ghost.pos[0] == doctest::Approx(2.0)); // 0.0 + 2.0
  CHECK(ghost.pos[1] == doctest::Approx(1.0)); // 1.0 + 0.0

  // Now let's create a vertical simple shear deformation
  // In this matrix, the top-right element (0,1) controls vertical shear
  // Setting it to 0.5 means: for each x unit, shift up by 0.5 units
  Matrix2d shearDeformation = Matrix2d::Identity();
  shearDeformation(1, 0) = 0.0; // No horizontal shift due to y
  shearDeformation(0, 1) = 0.0; // No horizontal shift due to y
  shearDeformation(1, 0) =
      0.5; // Vertical shift due to x: for each x unit, shift up by 0.5 units

  // Update the ghost node with the shear deformation
  ghost.updatePosition(referenceNode, shearDeformation);

  // Check positions after applying shear deformation
  // The reference node at (0,1) stays fixed because x=0 means no vertical shift
  // But the ghost node shifts upward because it has x=2, so 2*0.5=1 unit up
  CHECK(ghost.pos[0] == doctest::Approx(2.0)); // x position unchanged
  CHECK(ghost.pos[1] ==
        doctest::Approx(2.0)); // 1.0 + (0.5*2.0) = 1.0 + 1.0 = 2.0

  // Test with a different reference position to verify calculations
  Vector2d newPos = {0.0, 3.0}; // Still at x=0, but now at y=3
  referenceNode.setPos(newPos);

  // Update ghost with same shear deformation
  ghost.updatePosition(referenceNode, shearDeformation);

  // Check new positions - reference node moved up, ghost follows but still has
  // the shear effect
  CHECK(ghost.pos[0] == doctest::Approx(2.0)); // x position unchanged
  CHECK(ghost.pos[1] ==
        doctest::Approx(4.0)); // 3.0 + (0.5*2.0) = 3.0 + 1.0 = 4.0

  // Test a more complex deformation matrix
  Matrix2d complexDeformation = Matrix2d::Zero();
  complexDeformation(0, 0) = 1.5; // Stretch in x
  complexDeformation(1, 1) = 0.8; // Compress in y
  complexDeformation(1, 0) = 0.7; // Vertical shift due to x

  // Update with complex deformation
  ghost.updatePosition(referenceNode, complexDeformation);

  // Check deformed position:
  // pos = referencePos + deformation * periodShift
  // [x] = [0.0] + [1.5 0.0] * [2.0] = [0.0] + [3.0] = [3.0]
  // [y]   [3.0]   [0.7 0.8]   [0.0]   [3.0]   [1.4]   [4.4]
  CHECK(ghost.pos[0] == doctest::Approx(3.0)); // 0.0 + (1.5*2.0)
  CHECK(ghost.pos[1] == doctest::Approx(4.4)); // 3.0 + (0.7*2.0 + 0.8*0.0)
}