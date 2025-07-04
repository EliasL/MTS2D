#include "../src/Data/data_export.h"
#include "../src/Mesh/mesh.h"
#include "Eigen/Core"
#include "Mesh/node.h"
#include "Mesh/tElement.h"
#include "Simulation/simulation.h"
#include "run/doctest.h"
#include <algorithm>
#include <string>

/**
 * Save mesh state to a VTU file
 * @param mesh The mesh to save
 * @param name The name suffix for the output file
 */
void save(Mesh &mesh, std::string name) {
  static int fileNr = 0;
  mesh.loadSteps = fileNr; // loadSteps is used to name the files
  mesh.updateMesh();
  writeMeshToVtu(mesh, "reconnecting", "", name);
  fileNr++;
  // createCollection(getDataPath(name, ""), getOutputPath(name, ""));
  createCollection("reconnecting/data", "reconnecting");
}

/**
 * Perform common setup for mesh reconnecting tests
 * @param meshSize Size of the mesh (number of nodes per dimension)
 * @param isPeriodic Whether to use periodic boundary conditions
 * @param displacement Amount to displace upper nodes by
 * @param prefix Prefix for output filenames
 * @return Configured mesh with displacement applied
 */
Mesh setupMeshForReconnectingTest(int meshSize, bool isPeriodic,
                                  const double shear,
                                  const std::string &prefix) {
  // Create a mesh
  Mesh mesh(meshSize, meshSize, isPeriodic);
  save(mesh, prefix + "_InitialState");

  Matrix2d trans;
  trans << 1, shear, 0, 1;
  mesh.applyTransformation(trans);

  save(mesh, prefix + "_displaced");

  return mesh;
}

/**
 * Calculate the sum of forces across all nodes
 * @param mesh The mesh to calculate forces for
 * @return Sum of all nodal forces
 */
Vector2d calculateTotalForce(const Mesh &mesh) {
  Vector2d sumForce = {0, 0};
  for (int i = 0; i < mesh.nrNodes; i++) {
    sumForce += mesh.nodes(i).f;
  }
  return sumForce;
}

/**
 * Run reconnecting tests and verify conservation properties
 * @param mesh The mesh to test
 * @param row Row index for diagonal change
 * @param col Column index for diagonal change
 * @param useLeftDiagonal Whether to use left diagonal
 * @param prefix Prefix for output filenames
 */
void testReconnectingConservation(Mesh &mesh, int row, int col,
                                  bool useLeftDiagonal,
                                  const std::string &prefix) {
  // Make a copy of the original mesh
  Mesh oldMesh = mesh;

  // Calculate initial total force
  Vector2d initialForce = calculateTotalForce(oldMesh);

  // Perform reconnecting
  mesh.setDiagonal(row, col, useLeftDiagonal);
  save(mesh, prefix + "_reconnected");

  // Calculate new total force
  Vector2d newForce = calculateTotalForce(mesh);

  // Veryfy force is balanced
  CHECK(initialForce[0] == doctest::Approx(0));
  CHECK(initialForce[1] == doctest::Approx(0));
  CHECK(newForce[0] == doctest::Approx(0));
  CHECK(newForce[1] == doctest::Approx(0));

  // Verify energy conservation
  CHECK(mesh.totalEnergy == oldMesh.totalEnergy);

  // Verify element-wise energy conservation
  for (size_t e = 0;
       e < std::min(mesh.elements.size(), oldMesh.elements.size()); e++) {
    CHECK(mesh.elements[e].energy ==
          doctest::Approx(oldMesh.elements[e].energy));
  }
  // Verify element-wise forces;
  for (int i = 0; i < mesh.nrElements; i++) {
    TElement &e = mesh.elements[i];
    Vector2d sumForce = {0, 0};

    for (int j = 0; j < e.ghostNodes.size(); j++) {
      sumForce += e.ghostNodes[j].f;
    }

    CHECK(sumForce[0] == doctest::Approx(0));
    CHECK(sumForce[1] == doctest::Approx(0));
  }

  // Compare the energy and forces on the nodes
  for (int i = 0; i < mesh.nrNodes; i++) {
    CHECK(mesh.nodes(i).f[0] == doctest::Approx(oldMesh.nodes(i).f[0]));
    CHECK(mesh.nodes(i).f[1] == doctest::Approx(oldMesh.nodes(i).f[1]));
    if (mesh.usingPBC) {
      CHECK(mesh.nodes(i).f[0] == doctest::Approx(0));
      CHECK(mesh.nodes(i).f[1] == doctest::Approx(0));
    }
  }
}

// TEST_CASE("Simple mesh reconnecting") {
//   Mesh mesh = setupMeshForReconnectingTest(2, false, -0.3, "simple");
//   testReconnectingConservation(mesh, 0, 0, false, "simple");
// }

// TEST_CASE("Simple periodic mesh reconnecting") {
//   Mesh mesh = setupMeshForReconnectingTest(2, true, -0.3, "periodic");
//   testReconnectingConservation(mesh, 1, 1, false, "periodic");
// }

// TEST_CASE("Simple periodic large deformation mesh reconnecting") {
//   Mesh mesh = setupMeshForReconnectingTest(2, true, -2.3,
//   "large_deform_periodic"); testReconnectingConservation(mesh, 1, 1, false,
//   "large_deform_periodic");
// }

// TEST_CASE("Larger periodic mesh reconnecting") {
//   Mesh mesh = setupMeshForReconnectingTest(4, true, 0.1, "large_periodic");
//   testReconnectingConservation(mesh, 1, 1, false, "large_periodic");
// }

// TEST_CASE("Complex periodic mesh reconnecting") {
//   Mesh mesh =
//       setupMeshForReconnectingTest(4, false, 0.1, "complex_large_periodic");
//   // Now we also deform the middle part a bit
//   mesh.nodes(2, 2).setDisplacement({.4, -0.0});
//   save(mesh, "complex_large_periodic_deformed");

//   testReconnectingConservation(mesh, 1, 1, false, "complex_large_periodic");
// }

TEST_CASE("Remove elements from nodes") {
  int testRow = 1;
  int testCol = 2;
  Mesh mesh(3, 3);

  SUBCASE("Remove single element from nodes") {
    std::vector<int> elementsToRemove = {10, 0};

    auto nodes = mesh.getSquareNodes(testRow, testCol);

    // Store original element connections
    std::vector<std::vector<int>> originalElementIndices;
    std::vector<std::vector<int>> originalNodeIndices;
    for (auto *node : nodes) {
      std::vector<int> e, n;
      for (int i = 0; i < node->elementCount; ++i) {
        e.push_back(node->elementIndices[i]);
        n.push_back(node->nodeIndexInElement[i]);
      }
      originalElementIndices.push_back(e);
      originalNodeIndices.push_back(n);
    }

    mesh.removeElementsFromNodes(testRow, testCol, elementsToRemove);

    auto updatedNodes = mesh.getSquareNodes(testRow, testCol);

    for (size_t i = 0; i < updatedNodes.size(); ++i) {
      auto *node = updatedNodes[i];
      std::vector<int> expectedElementIndices;
      std::vector<int> expectedNodeIndices;

      for (size_t j = 0; j < originalElementIndices[i].size(); ++j) {
        if (std::find(elementsToRemove.begin(), elementsToRemove.end(),
                      originalElementIndices[i][j]) == elementsToRemove.end()) {
          expectedElementIndices.push_back(originalElementIndices[i][j]);
          expectedNodeIndices.push_back(originalNodeIndices[i][j]);
        }
      }

      CHECK(node->elementCount ==
            static_cast<int>(expectedElementIndices.size()));
      for (size_t j = 0; j < expectedElementIndices.size(); ++j) {
        CHECK(node->elementIndices[j] == expectedElementIndices[j]);
        CHECK(node->nodeIndexInElement[j] == expectedNodeIndices[j]);
      }
    }
  }
}

TEST_CASE("Check angle after reconnecting") {

  Mesh mesh(2, 2, false, "minor");

  mesh.applyTransformation(getShear(1));
  mesh.updateElements();
  // C_12 is not really an angle, but close enough
  double oldAngle1 = mesh.elements[0].largestAngle;
  double oldAngle2 = mesh.elements[1].largestAngle;
  // std::cout << mesh.elements[0] << '\n' << mesh.elements[1] << '\n';
  CHECK(oldAngle1 == doctest::Approx(135));
  CHECK(oldAngle2 == doctest::Approx(135));
  save(mesh, "AngleCheckBeforeReconnect");
  // mesh.setDiagonal(0, 0, false);

  // save(mesh, "AngleCheckAfterSetDiagonal");
  // double newAngle1 = mesh.elements[0].largestAngle;
  // double newAngle2 = mesh.elements[1].largestAngle;
  // CHECK(newAngle1 == doctest::Approx(135));
  // CHECK(newAngle2 == doctest::Approx(135));
  // // std::cout << mesh.elements[0] << '\n' << mesh.elements[1] << '\n';

  mesh.reconnect();
  double reconnectAngle1 = mesh.elements[0].largestAngle;
  double reconnectAngle2 = mesh.elements[1].largestAngle;
  // std::cout << mesh.elements[0] << '\n' << mesh.elements[1] << '\n';
  CHECK(reconnectAngle1 == doctest::Approx(90));
  CHECK(reconnectAngle2 == doctest::Approx(90));
  save(mesh, "AngleCheckAfterReconnect");
}

TEST_CASE("Check reconnecting with PBC") {

  Mesh mesh(2, 2, true, "major");

  mesh.applyTransformation(getShear(1));
  save(mesh, "PBCBeforeReconnect0");
  mesh.nodes(0, 1).addDisplacement({0, 0.3});
  mesh.nodes(1, 0).addDisplacement({0, 0.3});
  save(mesh, "PBCBeforeReconnect1");
  mesh.nodes(0, 1).addDisplacement({0, 0.7});
  mesh.nodes(1, 0).addDisplacement({0, 0.7});
  mesh.calculateAverages();
  save(mesh, "PBCBeforeReconnect2");
  mesh.reconnect();
  // The angle node of the first element should now be moved.
  save(mesh, "PBCAfterReconnect");
  CHECK(mesh.elements[0].getAngleNode()->pos == Vector2d{2, 1});

  // Check node-element connections
}

// Helper function to perform a mesh operation sequence.
void performMeshOperation(Mesh &mesh, double firstParam, double secondParam,
                          const Vector2d &direction, const std::string &label) {
  mesh.moveMeshSection(firstParam, secondParam, direction, true, true);
  mesh.calculateAverages();
  save(mesh, label);
  mesh.reconnect();
  save(mesh, label + "AfterReconnect");
}

TEST_CASE("Check multiple reconnecting") {
  Mesh mesh(5, 5, false, "minor");

  save(mesh, "MultiReconnect0");
  mesh.applyTransformation(getShear(1));
  mesh.updateMesh();
  mesh.reconnect();
  save(mesh, "MultiReconnect1");

  // Forward operations
  for (int i = 2; i < 5; i += 2) {
    // Move horizontally
    performMeshOperation(mesh, 0, i + 0.5 - 1, {1, 0},
                         "MutiReconnectSide" + std::to_string(i));
    // Move vertically
    performMeshOperation(mesh, i + 0.5, 0, {0, 1},
                         "MutiReconnectUp" + std::to_string(i));
  }

  // Backward operations
  for (int i = 4; i > 0; i -= 2) {
    // Move horizontally in the opposite direction
    performMeshOperation(mesh, i + 0.5, 0, {0, -1},
                         "backwardsMutiUpReconnect" + std::to_string(i));
    performMeshOperation(mesh, 0, i + 0.5 - 1, {-1, 0},
                         "backwardsMutiSideReconnect" + std::to_string(i));
    // Move vertically in the opposite direction
  }
}