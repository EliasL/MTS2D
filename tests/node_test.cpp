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
  Node n(1, 1, 2);
  CHECK(n.id.i == 3);
  CHECK(n.id.col == 1);
  CHECK(n.id.row == 1);
  CHECK(n.pos()[0] == 1);
  CHECK(n.pos()[1] == 1);
  CHECK(n.init_pos()[0] == 1);
  CHECK(n.init_pos()[1] == 1);
}
TEST_CASE("Ghost node initialization and displacement") {
  // Reference node at (1,0) in a 3x3 mesh (on the edge)
  Node referenceNode(1, 0, 3);
  Matrix2d deformation = Matrix2d::Identity();

  SUBCASE("Ghost node displacement with zero deformation") {
    GhostNode ghost(&referenceNode, 1, 3, 3, deformation);

    // Check initial positions
    CHECK(ghost.init_pos[0] == doctest::Approx(3.0));
    CHECK(ghost.init_pos[1] == doctest::Approx(1.0));

    // Check displacement (should be zero initially)
    CHECK(ghost.u[0] == doctest::Approx(0.0));
    CHECK(ghost.u[1] == doctest::Approx(0.0));
  }

  SUBCASE("Ghost node displacement with uniform stretch") {
    deformation(0, 0) = 2.0;
    deformation(1, 1) = 2.0;

    Node deformedReference = referenceNode;
    deformedReference.applyDeformation(deformation);
    GhostNode ghost(&deformedReference, 1, 3, 3, deformation);

    // Expected deformed position
    CHECK(ghost.pos[0] == doctest::Approx(6.0)); // (2.0 * 3.0)
    CHECK(ghost.pos[1] == doctest::Approx(2.0)); // (2.0 * 1.0)

    // Expected displacement
    CHECK(ghost.u[0] == doctest::Approx(3));
    CHECK(ghost.u[1] == doctest::Approx(1));
  }

  SUBCASE("Ghost node displacement with shear deformation") {
    deformation = Matrix2d::Identity();
    deformation(0, 1) = 1.0; // Horizontal shear

    Node deformedReference = referenceNode;
    deformedReference.applyDeformation(deformation);
    GhostNode ghost(&deformedReference, 1, 3, 3, deformation);

    // Expected position shift due to shear
    CHECK(ghost.pos[0] == doctest::Approx(1.0 + 3.0)); // 1.0 + (1.0 * 3.0)
    CHECK(ghost.pos[1] == doctest::Approx(1));         // 3.0 + (-0.5 * 3.0)

    // Expected displacement
    CHECK(ghost.u[0] == doctest::Approx((1.0)));
    CHECK(ghost.u[1] == doctest::Approx((0)));
  }
}
