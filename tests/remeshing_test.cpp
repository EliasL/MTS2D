#include "../src/Data/data_export.h"
#include "../src/Mesh/mesh.h"
#include "Eigen/Core"
#include "Mesh/node.h"
#include "Mesh/tElement.h"
#include "Simulation/simulation.h"
#include "run/doctest.h"
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
  writeMeshToVtu(mesh, "remeshing", "", name);
  fileNr++;
  // createCollection(getDataPath(name, ""), getOutputPath(name, ""));
  createCollection("remeshing/data", "remeshing");
}

/**
 * Perform common setup for mesh remeshing tests
 * @param meshSize Size of the mesh (number of nodes per dimension)
 * @param isPeriodic Whether to use periodic boundary conditions
 * @param displacement Amount to displace upper nodes by
 * @param prefix Prefix for output filenames
 * @return Configured mesh with displacement applied
 */
Mesh setupMeshForRemeshingTest(int meshSize, bool isPeriodic,
                               const double shear, const std::string &prefix) {
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
 * Run remeshing tests and verify conservation properties
 * @param mesh The mesh to test
 * @param row Row index for diagonal change
 * @param col Column index for diagonal change
 * @param useLeftDiagonal Whether to use left diagonal
 * @param prefix Prefix for output filenames
 */
void testRemeshingConservation(Mesh &mesh, int row, int col,
                               bool useLeftDiagonal,
                               const std::string &prefix) {
  // Make a copy of the original mesh
  Mesh oldMesh = mesh;

  // Calculate initial total force
  Vector2d initialForce = calculateTotalForce(oldMesh);

  // Perform remeshing
  mesh.setDiagonal(row, col, useLeftDiagonal);
  save(mesh, prefix + "_remeshed");

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

// TEST_CASE("Simple mesh remeshing") {
//   Mesh mesh = setupMeshForRemeshingTest(2, false, -0.3, "simple");
//   testRemeshingConservation(mesh, 0, 0, false, "simple");
// }

// TEST_CASE("Simple periodic mesh remeshing") {
//   Mesh mesh = setupMeshForRemeshingTest(2, true, -0.3, "periodic");
//   testRemeshingConservation(mesh, 1, 1, false, "periodic");
// }

// TEST_CASE("Simple periodic large deformation mesh remeshing") {
//   Mesh mesh = setupMeshForRemeshingTest(2, true, -2.3,
//   "large_deform_periodic"); testRemeshingConservation(mesh, 1, 1, false,
//   "large_deform_periodic");
// }

// TEST_CASE("Larger periodic mesh remeshing") {
//   Mesh mesh = setupMeshForRemeshingTest(4, true, 0.1, "large_periodic");
//   testRemeshingConservation(mesh, 1, 1, false, "large_periodic");
// }

// TEST_CASE("Complex periodic mesh remeshing") {
//   Mesh mesh =
//       setupMeshForRemeshingTest(4, false, 0.1, "complex_large_periodic");
//   // Now we also deform the middle part a bit
//   mesh.nodes(2, 2).setDisplacement({.4, -0.0});
//   save(mesh, "complex_large_periodic_deformed");

//   testRemeshingConservation(mesh, 1, 1, false, "complex_large_periodic");
// }

TEST_CASE("Remove elements from nodes") {
  // Create a test square with nodes
  int testRow = 1;
  int testCol = 2;
  Mesh mesh(3, 3);

  SUBCASE("Remove single element from nodes") {

    // Elements to remove
    std::vector<int> elementsToRemove = {10};

    // Execute

    auto updatedNodes = mesh.getSquareNodes(testRow, testCol);
    mesh.removeElementsFromNodes(testRow, testCol, elementsToRemove);

    // Verify: Get the updated nodes
    updatedNodes = mesh.getSquareNodes(testRow, testCol);

    // Check Node 0:
    // TODO
    CHECK(updatedNodes[0]->elementCount == 5);
    // CHECK(updatedNodes[0]->elementIndices[0] == 11);
    // CHECK(updatedNodes[0]->elementIndices[1] == 12);
    // CHECK(updatedNodes[0]->nodeIndexInElement[0] == 1);
    // CHECK(updatedNodes[0]->nodeIndexInElement[1] == 2);

    // // Check Node 1: Should have elements 13, 14 left
    // CHECK(updatedNodes[1]->elementCount == 2);
    // CHECK(updatedNodes[1]->elementIndices[0] == 13);
    // CHECK(updatedNodes[1]->elementIndices[1] == 14);
    // CHECK(updatedNodes[1]->nodeIndexInElement[0] == 0);
    // CHECK(updatedNodes[1]->nodeIndexInElement[1] == 1);
  }
}

TEST_CASE("Check angle after remeshing") {

  Mesh mesh(2, 2, false, "minor");

  mesh.applyTransformation(getShear(1));
  mesh.updateElements();
  // C_12 is not really an angle, but close enough
  double oldAngle1 = mesh.elements[0].largestAngle;
  double oldAngle2 = mesh.elements[1].largestAngle;
  // std::cout << mesh.elements[0] << '\n' << mesh.elements[1] << '\n';
  CHECK(oldAngle1 == doctest::Approx(135));
  CHECK(oldAngle2 == doctest::Approx(135));
  save(mesh, "AngleCheckBeforeRemesh");
  // mesh.setDiagonal(0, 0, false);

  // save(mesh, "AngleCheckAfterSetDiagonal");
  // double newAngle1 = mesh.elements[0].largestAngle;
  // double newAngle2 = mesh.elements[1].largestAngle;
  // CHECK(newAngle1 == doctest::Approx(135));
  // CHECK(newAngle2 == doctest::Approx(135));
  // // std::cout << mesh.elements[0] << '\n' << mesh.elements[1] << '\n';

  mesh.remesh();
  double remeshAngle1 = mesh.elements[0].largestAngle;
  double remeshAngle2 = mesh.elements[1].largestAngle;
  // std::cout << mesh.elements[0] << '\n' << mesh.elements[1] << '\n';
  CHECK(remeshAngle1 == doctest::Approx(90));
  CHECK(remeshAngle2 == doctest::Approx(90));
  save(mesh, "AngleCheckAfterRemesh");
}

TEST_CASE("Check remeshing with PBC") {

  Mesh mesh(2, 2, true, "major");

  mesh.applyTransformation(getShear(1));
  save(mesh, "PBCBeforeRemesh0");
  mesh.nodes(0, 1).addDisplacement({0, 0.3});
  mesh.nodes(1, 0).addDisplacement({0, 0.3});
  save(mesh, "PBCBeforeRemesh1");
  mesh.nodes(0, 1).addDisplacement({0, 0.7});
  mesh.nodes(1, 0).addDisplacement({0, 0.7});
  mesh.calculateAverages();
  save(mesh, "PBCBeforeRemesh2");
  mesh.remesh();
  // The angle node of the first element should now be moved.
  CHECK(mesh.elements[0].getAngleNode()->pos == Vector2d{3, 1});
  save(mesh, "PBCAfterRemesh");
}