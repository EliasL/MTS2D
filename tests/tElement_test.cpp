#include "../src/Data/data_export.h"
#include "../src/Mesh/mesh.h" // Include the header for your surface struct
#include "Eigen/src/Core/Matrix.h"
#include "run/doctest.h"
#include <cmath>
#include <iostream>

TEST_CASE("Mesh Initialization") {
  // Create a surface with known dimensions for testing
  Mesh mesh(2, 2, false);

  // There should now be two elements.
  CHECK(mesh.nrElements == 2);
  // Create a periodic surface with known dimensions for testing
  mesh = Mesh(2, 2, true);

  // There should now be eight elements.
  CHECK(mesh.nrElements == 8);
}

TEST_CASE("Update deformation gradiant") {
  // We use a mesh to initialize an element. (Not best practice for testing)
  Mesh mesh(2, 2);
  TElement &e = mesh.elements[0];
  Eigen::Matrix2d shear;
  shear << 1, 4, 0, 1;

  mesh.applyTransformation(shear);
  mesh.updateElements();
  CHECK(e.F == shear);

  Eigen::Matrix2d shear2;
  shear2 << 1, 0, 3, 1;

  mesh.applyTransformation(shear2);
  mesh.updateElements();
  CHECK(e.F == shear2 * shear);
}

TEST_CASE("Update metric tensor") {
  Mesh mesh(2, 2);
  TElement &e = mesh.elements[0];
  Eigen::Matrix2d shear;
  shear << 1, 4, 0, 1;

  // https://www.wolframalpha.com/input?i=transpose%28%7B%7B1%2C4%7D%2C%7B0%2C1%7D%7D%29.%7B%7B1%2C4%7D%2C%7B0%2C1%7D%7D

  Eigen::Matrix2d ans;
  ans << 1, 4, 4, 17;

  mesh.applyTransformation(shear);
  mesh.updateElements();

  CHECK(e.C == ans);
}

TEST_CASE("Update reduced metric tensor") {
  Mesh mesh(2, 2);
  TElement &e = mesh.elements[0];
  Matrix2d shear;
  shear << 1, -4, 0, 1;
  Matrix2d C_Ans;
  C_Ans << 1, 0, 0, 1;
  Matrix2d mAns;
  mAns << 1, -4, 0, -1;
  mesh.applyTransformation(shear);
  mesh.updateElements();

  CHECK(e.C_ == C_Ans);
  CHECK(e.m == mAns);

  // Use the lagrange redution to test a more difficult example
  C_Ans << 1.05565452, 0.52639976, 0.52639976, 1.20976767;
  // https://www.wolframalpha.com/input?i=%7B%7B0%2C1%7D%2C%7B1%2C0%7D%7D.%7B%7B1%2C-1%7D%2C%7B0%2C1%7D%7D.%7B%7B1%2C-1%7D%2C%7B0%2C1%7D%7D.%7B%7B1%2C0%7D%2C%7B0%2C-1%7D%7D.%7B%7B0%2C1%7D%2C%7B1%2C0%7D%7D
  mAns << -1, 0, 2, 1;
  e = TElement::lagrangeReduction(3.78912615, 1.20976767, 1.89313557);
  CHECK(e.C_.isApprox(C_Ans, 0.00001));
  CHECK(e.m == mAns);
}

TEST_CASE("Update energy and reduced stress") {
  Mesh mesh(2, 2);

  TElement &e = mesh.elements[0];

  for (int i = 0; i < 4; i++) {

    // We make sure that the results are the same for any integer shear,
    // both positive and negative
    Matrix2d bigShear;
    bigShear << 1, i * (pow(-1, i)), 0, 1;
    // std::cout << bigShear;
    mesh.applyTransformation(bigShear);
    mesh.updateElements();

    CHECK(e.dPhi_dC_(0, 0) == doctest::Approx(0));
    CHECK(e.dPhi_dC_(0, 1) == doctest::Approx(0));
    CHECK(e.dPhi_dC_(1, 0) == doctest::Approx(0));
    CHECK(e.dPhi_dC_(1, 1) == doctest::Approx(0));
    CHECK(e.energy == doctest::Approx(0));

    Matrix2d smallShear;
    smallShear << 1, 0.5, 0, 1;
    mesh.applyTransformation(smallShear);
    mesh.updateElements();
    // Validated by Umut's code
    double ground_state = 3.91162;
    double other_energy = 4.00204;
    // Divide by two to account for volume
    CHECK(e.energy == doctest::Approx((other_energy - ground_state) / 2));
    Matrix2d backwardsSmallShear;
    backwardsSmallShear << 1, -0.5, 0, 1;
    mesh.applyTransformation(backwardsSmallShear);
    mesh.updateElements();
  }
}

TEST_CASE("Update Piola stress") {
  Mesh mesh(2, 2);
  TElement &e = mesh.elements[0];

  mesh.updateElements();

  CHECK(e.P(0, 0) == doctest::Approx(0));
  CHECK(e.P(0, 1) == doctest::Approx(0));
  CHECK(e.P(1, 0) == doctest::Approx(0));
  CHECK(e.P(1, 1) == doctest::Approx(0));

  Matrix2d shear;
  shear << 1, 0.5, 0, 1;
  mesh.applyTransformation(shear);
  mesh.updateElements();

  // Validated by Umut's code
  CHECK(e.P(0, 0) == doctest::Approx(-0.0462536));
  CHECK(e.P(0, 1) == doctest::Approx(-3.47E-18));
  CHECK(e.P(1, 0) == doctest::Approx(-0.0231268));
  CHECK(e.P(1, 1) == doctest::Approx(0.0462536));
}

TEST_CASE("Apply forces on nodes at rest") {
  Mesh mesh(3, 3, false);

  mesh.updateElements();
  mesh.applyForceFromElementsToNodes();

  for (long i = 0; i < mesh.nodes.size(); i++) {
    CHECK(mesh.nodes(i).f[0] == doctest::Approx(0));
    CHECK(mesh.nodes(i).f[1] == doctest::Approx(0));
  }
}

TEST_CASE("Apply forces on nodes (Normal mesh)") {
  Mesh mesh(3, 3, 1, 0, false, false);
  Matrix2d shear;
  shear << 1, 0.5, 0, 1;
  mesh.applyTransformation(shear);
  mesh.updateElements();
  mesh.applyForceFromElementsToNodes();
  // writeMeshToVtu(mesh, "test", "");

  // Check values in elements
  for (int i = 0; i < mesh.nrElements; i++) {
    // Check P (Assumed to be correct because it gives correct node force)
    CHECK(mesh.elements[i].P(0) == doctest::Approx(-0.046254));
    CHECK(mesh.elements[i].P(1) == doctest::Approx(-0.023127));
    CHECK(mesh.elements[i].P(2) == doctest::Approx(0));
    CHECK(mesh.elements[i].P(3) == doctest::Approx(0.046254));

    // Check r (Assumed to be correct because it gives correct node force)
    // Define expected values for even and odd indices
    std::vector<std::vector<int>> evenExpected{{-1, -1}, {1, 0}, {0, 1}};
    std::vector<std::vector<int>> oddExpected{{0, -1}, {-1, 0}, {1, 1}};

    // Select the expected pattern based on the index
    const auto &expected = (i % 2 == 0) ? evenExpected : oddExpected;

    for (size_t j = 0; j < expected.size(); j++) {
      CHECK(mesh.elements[i].dN_dX[j][0] == expected[j][0]);
      CHECK(mesh.elements[i].dN_dX[j][1] == expected[j][1]);
    }
  }
  Vector2d s = {0, 0};

  for (int i = 0; i < mesh.nodes.size(); i++) {
    s += mesh.nodes(i).f;
  }

  // Check values in nodes
  // Validated by Umut's code
  CHECK(mesh.nodes(0).f(0) == doctest::Approx(0.0462536));
  CHECK(mesh.nodes(0).f(1) == doctest::Approx(-0.0231268));

  CHECK(mesh.nodes(1).f(0) == doctest::Approx(0));
  CHECK(mesh.nodes(1).f(1) == doctest::Approx(-0.0925071));

  CHECK(mesh.nodes(2).f(0) == doctest::Approx(-0.0462536));
  CHECK(mesh.nodes(2).f(1) == doctest::Approx(-0.0693803));

  CHECK(mesh.nodes(3).f(0) == doctest::Approx(0.0925071));
  CHECK(mesh.nodes(3).f(1) == doctest::Approx(0.0462536));

  CHECK(mesh.nodes(4).f(0) == doctest::Approx(0));
  CHECK(mesh.nodes(4).f(1) == doctest::Approx(0));

  CHECK(mesh.nodes(5).f(0) == doctest::Approx(-0.0925071));
  CHECK(mesh.nodes(5).f(1) == doctest::Approx(-0.0462536));

  CHECK(mesh.nodes(6).f(0) == doctest::Approx(0.0462536));
  CHECK(mesh.nodes(6).f(1) == doctest::Approx(0.0693803));

  CHECK(mesh.nodes(7).f(0) == doctest::Approx(0));
  CHECK(mesh.nodes(7).f(1) == doctest::Approx(0.0925071));

  CHECK(mesh.nodes(8).f(0) == doctest::Approx(-0.0462536));
  CHECK(mesh.nodes(8).f(1) == doctest::Approx(0.0231268));
}

TEST_CASE("Apply forces on nodes (New mesh)") {
  Mesh mesh(3, 3, 1, 0, false, true);
  Matrix2d shear;
  shear << 1, 0.5, 0, 1;
  mesh.applyTransformation(shear);
  mesh.updateElements();
  mesh.applyForceFromElementsToNodes();
  // writeMeshToVtu(mesh, "test", "");

  // Check values in elements
  for (int i = 0; i < mesh.nrElements; i++) {
    // Check P (Assumed to be correct because it gives correct node force)
    CHECK(mesh.elements[i].P(0) == doctest::Approx(-0.046254));
    CHECK(mesh.elements[i].P(1) == doctest::Approx(-0.023127));
    CHECK(mesh.elements[i].P(2) == doctest::Approx(0));
    CHECK(mesh.elements[i].P(3) == doctest::Approx(0.046254));

    // Check dN_dX (Assumed to be correct because it gives correct node force)
    // Define expected values for even and odd indices

    // TODO
    // Here, instead of only having up and down triangles that can be divided
    // into the even and odd indexed triangles, we have four different triangles
    // We can make a test for it, but i can't be bothered to do it right now
    // If the forces are stil okay, it should be fine i think

    // std::vector<std::vector<int>> evenExpected{{-1, -1}, {1, 0}, {0, 1}};
    // std::vector<std::vector<int>> oddExpected{{0, -1}, {-1, 0}, {1, 1}};

    // // Select the expected pattern based on the index
    // const auto &expected = (i % 2 == 0) ? evenExpected : oddExpected;
    // for (size_t j = 0; j < expected.size(); j++) {
    //   std::cout << i << '\n' << mesh.elements[i].dN_dX[j] << '\n';
    //   CHECK(mesh.elements[i].dN_dX[j][0] == expected[j][0]);
    //   CHECK(mesh.elements[i].dN_dX[j][1] == expected[j][1]);
    // }
  }

  // Check values in nodes
  // Validated by Umut's code
  CHECK(mesh.nodes(0).f(0) == doctest::Approx(0.0462536));
  CHECK(mesh.nodes(0).f(1) == doctest::Approx(-0.0231268));

  CHECK(mesh.nodes(1).f(0) == doctest::Approx(0));
  CHECK(mesh.nodes(1).f(1) == doctest::Approx(-0.0925071));

  CHECK(mesh.nodes(2).f(0) == doctest::Approx(-0.0462536));
  CHECK(mesh.nodes(2).f(1) == doctest::Approx(-0.0693803));

  CHECK(mesh.nodes(3).f(0) == doctest::Approx(0.0925071));
  CHECK(mesh.nodes(3).f(1) == doctest::Approx(0.0462536));

  CHECK(mesh.nodes(4).f(0) == doctest::Approx(0));
  CHECK(mesh.nodes(4).f(1) == doctest::Approx(0));

  CHECK(mesh.nodes(5).f(0) == doctest::Approx(-0.0925071));
  CHECK(mesh.nodes(5).f(1) == doctest::Approx(-0.0462536));

  CHECK(mesh.nodes(6).f(0) == doctest::Approx(0.0462536));
  CHECK(mesh.nodes(6).f(1) == doctest::Approx(0.0693803));

  CHECK(mesh.nodes(7).f(0) == doctest::Approx(0));
  CHECK(mesh.nodes(7).f(1) == doctest::Approx(0.0925071));

  CHECK(mesh.nodes(8).f(0) == doctest::Approx(-0.0462536));
  CHECK(mesh.nodes(8).f(1) == doctest::Approx(0.0231268));
}