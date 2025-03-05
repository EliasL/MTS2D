#include "../src/Data/data_export.h"
#include "../src/Mesh/mesh.h"
#include "Eigen/src/Core/Matrix.h"
#include "Mesh/node.h"
#include "Mesh/tElement.h"
#include "run/doctest.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

template <typename DerivedA, typename DerivedB>
void printMatrixSideBySide(const Eigen::MatrixBase<DerivedA> &actual,
                           const Eigen::MatrixBase<DerivedB> &expected) {
  std::cout << "  Actual\t\t\t\tExpected\n";
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      std::cout << std::setw(12) << actual(i, j) << " ";
    }
    std::cout << "\t";
    for (int j = 0; j < 2; ++j) {
      std::cout << std::setw(12) << expected(i, j) << " ";
    }
    std::cout << '\n';
  }
}

// Custom function to check matrix equality with doctest::Approx
template <typename DerivedA, typename DerivedB>
bool checkMatrixApprox(const Eigen::MatrixBase<DerivedA> &actual,
                       const Eigen::MatrixBase<DerivedB> &expected,
                       double epsilon = 1e-15) {
  if (actual.rows() != expected.rows() || actual.cols() != expected.cols()) {
    return false;
  }

  for (int i = 0; i < actual.rows(); i++) {
    for (int j = 0; j < actual.cols(); j++) {
      if (!(actual(i, j) == doctest::Approx(expected(i, j)).epsilon(epsilon))) {
        std::cout << "Error: actual is not the same as expected.\n";
        printMatrixSideBySide(actual, expected);
        return false;
      }
    }
  }

  return true;
}

TEST_CASE("Mesh Initialization") {
  // Test non-periodic mesh
  Mesh mesh(2, 2, false);
  CHECK(mesh.nrElements == 2);

  // Test periodic mesh
  mesh = Mesh(2, 2, true);
  CHECK(mesh.nrElements == 8);
}

TEST_CASE("Element Property Updates") {
  SUBCASE("Deformation gradient") {
    Mesh mesh(3, 3);

    // First transformation
    Eigen::Matrix2d shear;
    shear << 1, 1, 0, 1;
    mesh.applyTransformation(shear);
    mesh.updateElements();
    for (const auto &e : mesh.elements) {
      CHECK(e.F == shear);
    }
    // Composite transformation
    Eigen::Matrix2d shear2;
    shear2 << 1, 1, 0, 1;
    mesh.applyTransformation(shear2);
    mesh.updateElements();
    for (const auto &e : mesh.elements) {
      CHECK(e.F == shear2 * shear);
    }

    // Composite transformation
    Eigen::Matrix2d shear3;
    shear3 << 1, 0, 1, 1;
    mesh.applyTransformation(shear3);
    mesh.updateElements();
    for (const auto &e : mesh.elements) {
      CHECK(e.F == shear3 * shear2 * shear);
    }
  }

  SUBCASE("Metric tensor") {
    Mesh mesh(2, 2);
    TElement &e = mesh.elements[0];

    Eigen::Matrix2d shear;
    shear << 1, 4, 0, 1;
    Eigen::Matrix2d expected;
    expected << 1, 4, 4, 17;

    mesh.applyTransformation(shear);
    mesh.updateElements();
    CHECK(e.C == expected);
  }

  SUBCASE("Reduced metric tensor") {
    Mesh mesh(2, 2);
    TElement &e = mesh.elements[0];

    // Simple case
    Matrix2d shear;
    shear << 1, -4, 0, 1;
    Matrix2d C_expected;
    C_expected << 1, 0, 0, 1;
    Matrix2d m_expected;
    m_expected << 1, -4, 0, -1;

    mesh.applyTransformation(shear);
    mesh.updateElements();
    CHECK(e.C_ == C_expected);
    CHECK(e.m == m_expected);

    // Complex case using Lagrange reduction
    C_expected << 1.05565452, 0.52639976, 0.52639976, 1.20976767;
    m_expected << -1, 0, 2, 1;
    e = TElement::lagrangeReduction(3.78912615, 1.20976767, 1.89313557);
    CHECK(e.C_.isApprox(C_expected, 0.00001));
    CHECK(e.m == m_expected);
  }

  SUBCASE("Energy and reduced stress") {
    Mesh mesh(2, 2);
    TElement &e = mesh.elements[0];

    // Test various shear values
    for (int i = 0; i < 4; i++) {
      // Integer shear test
      Matrix2d bigShear;
      bigShear << 1, i * (pow(-1, i)), 0, 1;
      mesh.applyTransformation(bigShear);
      mesh.updateElements();

      // Check zero energy state
      CHECK(checkMatrixApprox(e.dPhi_dC_, Matrix2d::Zero()));
      CHECK(e.energy == doctest::Approx(0));

      // Small shear test
      Matrix2d smallShear;
      smallShear << 1, 0.5, 0, 1;
      mesh.applyTransformation(smallShear);
      mesh.updateElements();

      // Validated by Umut's code
      double ground_state = 3.91162;
      double other_energy = 4.00204;
      CHECK(e.energy == doctest::Approx((other_energy - ground_state) / 2));

      // Reverse shear
      Matrix2d backwardsSmallShear;
      backwardsSmallShear << 1, -0.5, 0, 1;
      mesh.applyTransformation(backwardsSmallShear);
      mesh.updateElements();
    }
  }

  SUBCASE("Piola stress") {
    Mesh mesh(2, 2);
    TElement &e = mesh.elements[0];

    // Initial state
    mesh.updateElements();
    CHECK(checkMatrixApprox(e.P, Matrix2d::Zero()));

    // After shear
    Matrix2d shear;
    shear << 1, 0.5, 0, 1;
    mesh.applyTransformation(shear);
    mesh.updateElements();

    // Validated values
    Matrix2d expectedP;
    expectedP << -0.046253551136363619, 0, -0.023126775568181813,
        0.046253551136363619;
    CHECK(checkMatrixApprox(e.P, expectedP));
  }
}

// Helper function to check expected dN_dX values for different triangle types
void checkTriangleDN_dX(const TElement &element, bool useMajorDiagonal,
                        bool isFirstTriangle) {
  // old triangulation
  // std::vector<std::vector<int>> majorLeft{{0, -1}, {-1, 0}, {1, 1}};
  // std::vector<std::vector<int>> majorRight{{-1, -1}, {1, 0}, {0, 1}};
  // std::vector<std::vector<int>> minorLeft{{-1, 0}, {1, -1}, {0, 1}};
  // std::vector<std::vector<int>> minorRight{{0, -1}, {-1, 1}, {1, 0}};
  std::vector<std::vector<int>> majorLeft{{1, 1}, {0, -1}, {-1, 0}};
  std::vector<std::vector<int>> majorRight{{-1, -1}, {1, 0}, {0, 1}};
  std::vector<std::vector<int>> minorLeft{{1, -1}, {-1, 0}, {0, 1}};
  std::vector<std::vector<int>> minorRight{{-1, 1}, {0, -1}, {1, 0}};

  const auto &expected = useMajorDiagonal
                             ? (isFirstTriangle ? majorLeft : majorRight)
                             : (isFirstTriangle ? minorLeft : minorRight);

  for (size_t j = 0; j < expected.size(); j++) {
    // std::cout << (useMajorDiagonal ? "Major" : "Minor")
    //           << (isFirstTriangle ? " left" : " right") << '\n';
    // std::cout << element.dN_dX[0] << ", " << element.dN_dX[1] << ", "
    //           << element.dN_dX[2] << '\n';
    CHECK(element.dN_dX[j][0] == expected[j][0]);
    CHECK(element.dN_dX[j][1] == expected[j][1]);
  }
}

// Helper function to check expected node forces after shear
void checkNodeForces(const Mesh &mesh, bool isPeriodic) {
  if (isPeriodic) {
    // For periodic boundary conditions, forces should sum to zero at each node
    for (int i = 0; i < mesh.nodes.size(); i++) {
      CHECK(mesh.nodes(i).f(0) == doctest::Approx(0.0).epsilon(1e-10));
      CHECK(mesh.nodes(i).f(1) == doctest::Approx(0.0).epsilon(1e-10));
    }
  } else {
    // Expected values from Umut's code validation
    std::vector<Vector2d> expectedForces = {
        {0.0462536, -0.0231268},  // Node 0
        {0, -0.0925071},          // Node 1
        {-0.0462536, -0.0693803}, // Node 2
        {0.0925071, 0.0462536},   // Node 3
        {0, 0},                   // Node 4
        {-0.0925071, -0.0462536}, // Node 5
        {0.0462536, 0.0693803},   // Node 6
        {0, 0.0925071},           // Node 7
        {-0.0462536, 0.0231268}   // Node 8
    };

    for (int i = 0; i < mesh.nodes.size(); i++) {
      CHECK(mesh.nodes(i).f(0) == doctest::Approx(expectedForces[i][0]));
      CHECK(mesh.nodes(i).f(1) == doctest::Approx(expectedForces[i][1]));
    }
  }
}

TEST_CASE("Forces on nodes") {
  SUBCASE("Nodes at rest") {
    Mesh mesh(3, 3, false);
    mesh.updateElements();
    mesh.applyForceFromElementsToNodes();

    for (long i = 0; i < mesh.nodes.size(); i++) {
      CHECK(mesh.nodes(i).f[0] == doctest::Approx(0));
      CHECK(mesh.nodes(i).f[1] == doctest::Approx(0));
    }
  }

  const doctest::String normalMesh = "Normal mesh";
  const doctest::String newMesh = "New mesh";
  const doctest::String newPBCMesh = "New PBC mesh";
  std::vector<std::tuple<std::tuple<bool, std::string>, doctest::String>>
      meshConfigs = {
          {{false, "major"}, normalMesh},   // Normal mesh
          {{false, "alternate"}, newMesh},  // New mesh
          {{true, "alternate"}, newPBCMesh} // New PBC mesh
      };

  // Fix: Proper structured binding syntax with parentheses
  for (const auto &[config, caseName] : meshConfigs) {
    const auto &[isPeriodic, useDiagonalFlipping] = config;
    SUBCASE(caseName) {

      Mesh mesh(3, 3, 1, 0, isPeriodic, useDiagonalFlipping);
      Matrix2d shear;
      shear << 1, 0.5, 0, 1;
      mesh.applyTransformation(shear);
      mesh.updateElements();
      mesh.applyForceFromElementsToNodes();

      // Check element properties
      for (int i = 0; i < mesh.nrElements; i++) {
        // Check P values
        Matrix2d expectedP;
        expectedP << -0.046253551136363619, 0, -0.023126775568181813,
            0.046253551136363619;
        CHECK(checkMatrixApprox(mesh.elements[i].P, expectedP));

        // Check dN_dX values - for new mesh type only
        if (useDiagonalFlipping == "alternate") {
          int e1i = i / 2; // Triangle 1 base index
          int row = e1i / mesh.ePairCols();
          int col = e1i % mesh.ePairCols();
          bool useMajorDiagonal = (row + col) % 2;
          bool isFirstTriangle = (i % 2 == useMajorDiagonal);

          checkTriangleDN_dX(mesh.elements[i], useMajorDiagonal,
                             isFirstTriangle);
        }
      }

      // Check node forces
      checkNodeForces(mesh, isPeriodic);
    }
  }
}