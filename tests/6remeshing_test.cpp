#include "../src/Data/data_export.h"
#include "../src/Mesh/mesh.h"
#include "Eigen/Core"
#include "Mesh/node.h"
#include "Mesh/tElement.h"
#include "Simulation/scenarios.h"
#include "Simulation/simulation.h"
#include "element_history_export.h"
#include "run/doctest.h"
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(IDE_LIGHTWEIGHT)
#define CGAL_TEST_SKIP *doctest::skip(true)
#else
#define CGAL_TEST_SKIP
#endif

namespace {

const std::string kRemeshingTestDataPath = "test_data";

std::string sanitizeFolderName(const std::string &name) {
  std::string out;
  out.reserve(name.size());
  bool lastWasSeparator = false;
  for (unsigned char c : name) {
    if (std::isalnum(c)) {
      out.push_back(static_cast<char>(c));
      lastWasSeparator = false;
    } else if (!lastWasSeparator) {
      out.push_back('_');
      lastWasSeparator = true;
    }
  }
  while (!out.empty() && out.back() == '_') {
    out.pop_back();
  }
  return out.empty() ? "unnamed_test" : out;
}

std::string currentRemeshingTestFolder(const std::string &groupName) {
  const doctest::ContextOptions *options = doctest::getContextOptions();
  const char *testName = (options != nullptr && options->currentTest != nullptr)
                             ? options->currentTest->m_name
                             : "manual_debug";
  return groupName + "/" + sanitizeFolderName(testName);
}

int &saveCounterForFolder(const std::string &folderName) {
  static std::map<std::string, int> counters;
  return counters[folderName];
}

void clearRemeshingFolderOnFirstUse(const std::string &folderName) {
  static std::set<std::string> clearedFolders;
  if (clearedFolders.insert(folderName).second) {
    clearOutputFolder(folderName, kRemeshingTestDataPath);
  }
}

} // namespace

/**
 * Save mesh state to a VTU file
 * @param mesh The mesh to save
 * @param name The name suffix for the output file
 */
void save(Mesh &mesh, std::string name) {
  const std::string folderName = currentRemeshingTestFolder("reconnecting");
  clearRemeshingFolderOnFirstUse(folderName);
  int &fileNr = saveCounterForFolder(folderName);
  mesh.loadSteps = fileNr; // loadSteps is used to name the files
  mesh.ensureFull();
  writeMeshToVtu(mesh, folderName, kRemeshingTestDataPath, name);
  fileNr++;
  createCollection(getDataPath(folderName, kRemeshingTestDataPath),
                   getOutputPath(folderName, kRemeshingTestDataPath));
}

void saveCurrentAndReference(Mesh &mesh, const std::string &name) {
  const std::string folderName =
      currentRemeshingTestFolder("reconnectingReferenceTest");
  clearRemeshingFolderOnFirstUse(folderName);
  int &fileNr = saveCounterForFolder(folderName);
  mesh.loadSteps = fileNr; // loadSteps is used to name the files
  mesh.ensureFull();
  writeMeshToVtu(mesh, folderName, kRemeshingTestDataPath, name);
  writeMeshToVtu(mesh, folderName, kRemeshingTestDataPath, name, false,
                 VtuFieldLevel::All, "reference", "", true);
  fileNr++;
  createCollection(getDataPath(folderName, kRemeshingTestDataPath),
                   getOutputPath(folderName, kRemeshingTestDataPath), "",
                   ".vtu", {}, COLLECTIONNAME, "", "_reference");
  createCollection(getDataPath(folderName, kRemeshingTestDataPath),
                   getOutputPath(folderName, kRemeshingTestDataPath), "",
                   ".vtu", {}, "reference_collection", "_reference", "");
}

void debugReconnect(Mesh &mesh, const std::string &name) {
  saveCurrentAndReference(mesh, name + "BeforeReconnect");
  mesh.reconnect();
  saveCurrentAndReference(mesh, name + "AfterReconnect");
  // std::cout << "F:\n" << mesh.elements[0].F << "\n";
}

void debugElement(std::array<GhostNode, 3> e, std::string name = "") {
  std::cout << name << ":\n";
  Matrix2d F = tElementF(e);
  std::cout << "F:\n" << F << "\n";
  std::cout << "ref\n";
  for (GhostNode n : e) {
    std::cout << n.ref_pos << "\n";
  }
  std::cout << "pos\n";
  for (GhostNode n : e) {
    std::cout << n.pos << "\n";
  }
}

void forceEdgeFlipOnFirstElementPair(Mesh &mesh) {
  mesh.ensureGeometry();

  TElement &e1 = mesh.elements[0];
  TElement &e2 = mesh.elements[1];
  const int e1i = e1.eIndex;
  const int e2i = e2.eIndex;
  const auto oe1 = e1.getAngleCo1Co2Nodes();
  const auto oe2 = e2.getAngleCo1Co2Nodes();
  const std::array<GhostNode, 3> n1a = {oe1[0], oe1[1], oe2[0]};
  const std::array<GhostNode, 3> n2a = {oe2[0], oe2[1], oe1[0]};

  const TElement::EdgeFlipRemeshState state1 =
      e1.evaluateEdgeFlipRemeshState(n1a, mesh.F_P_H[static_cast<size_t>(e1i)]);
  const TElement::EdgeFlipRemeshState state2 =
      e2.evaluateEdgeFlipRemeshState(n2a, mesh.F_P_H[static_cast<size_t>(e2i)]);

  mesh.removeElementFromNodes(e1);
  mesh.removeElementFromNodes(e2);
  mesh.createElementPair(state1.newGhostNodes, state2.newGhostNodes, e1i, e2i,
                         true);
  mesh.F_P_H[static_cast<size_t>(e1i)] = state1.H_new;
  mesh.F_P_H[static_cast<size_t>(e2i)] = state2.H_new;
  mesh.markDirty();
  // saveCurrentAndReference(mesh, "Edge flipped");
}

// Canonical (sorted) triple of reference node ids (linearized)
// static inline std::array<int, 3> tri_sig(const TElement &e) {
//   std::array<int, 3> s = {e.ghostNodes[0].referenceId.i,
//                           e.ghostNodes[1].referenceId.i,
//                           e.ghostNodes[2].referenceId.i};
//   std::sort(s.begin(), s.end());
//   return s;
// }

// Can be used to compare two meshes for equality
// static std::vector<std::array<int, 3>> tri_connectivity(const Mesh &m) {
//   std::vector<std::array<int, 3>> v;
//   v.reserve(m.nrElements);
//   for (const auto &e : m.elements)
//     v.push_back(tri_sig(e));
//   std::sort(v.begin(), v.end());
//   return v;
// }
// Verifies that every node's connectivity mirrors the elements array
static void check_node_connectivity_consistency(const Mesh &m) {
  // Build reverse map: for each node, which (element,localIdx) pairs reference
  // it?
  std::vector<std::vector<std::pair<int, int>>> rev(m.nrNodes);
  for (int ei = 0; ei < m.nrElements; ++ei) {
    const auto &e = m.elements[ei];
    for (int k = 0; k < 3; ++k) {
      int ni = e.ghostNodes[k].referenceId.i;
      rev[ni].push_back({ei, k});
    }
  }
  // Compare to node.connectedElements/nodeIndexInElement
  for (int ni = 0; ni < m.nrNodes; ++ni) {
    const Node &n = m.nodes(ni);
    std::vector<std::pair<int, int>> got, want = rev[ni];
    for (int t = 0; t < n.elementCount; ++t)
      got.push_back({n.connectedElements[t], n.nodeIndexInElement[t]});
    std::sort(got.begin(), got.end());
    std::sort(want.begin(), want.end());
    CHECK(got == want);
  }
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
        e.push_back(node->connectedElements[i]);
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
        CHECK(node->connectedElements[j] == expectedElementIndices[j]);
        CHECK(node->nodeIndexInElement[j] == expectedNodeIndices[j]);
      }
    }
  }
}

void setNodeState(Node &node, const Vector2d &refPos, const Vector2d &pos) {
  node.setRefPos(refPos);
  node.setPos(pos);
}

GhostNode makeGhostWithReference(const Node &node, const Vector2d &refPos) {
  GhostNode ghost(&node, Matrix2d::Identity());
  ghost.updateReferencePosition(refPos);
  return ghost;
}

void setEmpiricalEdgeCaseRealNodes(Mesh &mesh,
                                   const std::array<Vector2d, 4> &positions) {
  // Map the four real nodes in the crash dump onto a 2x2 mesh:
  // node 0 -> former refId 6
  // node 1 -> former refId 55
  // node 2 -> former refId 56
  // node 3 -> former refId 105
  const std::array<Vector2d, 4> refPositions = {
      Vector2d(6.0, 0.0), Vector2d(5.0, 1.0), Vector2d(6.0, 1.0),
      Vector2d(5.0, 2.0)};
  for (size_t i = 0; i < positions.size(); ++i) {
    setNodeState(mesh.nodes(static_cast<long>(i)), refPositions[i],
                 positions[i]);
  }
}

static inline std::array<int, 3> triSig(const TElement &e) {
  std::array<int, 3> s = {e.ghostNodes[0].referenceId.i,
                          e.ghostNodes[1].referenceId.i,
                          e.ghostNodes[2].referenceId.i};
  std::sort(s.begin(), s.end());
  return s;
}

static std::vector<std::array<int, 3>> triConnectivity(const Mesh &m) {
  std::vector<std::array<int, 3>> v;
  v.reserve(m.nrElements);
  for (const auto &e : m.elements) {
    v.push_back(triSig(e));
  }
  std::sort(v.begin(), v.end());
  return v;
}

TEST_CASE("Logged simple-shear edge flip selects finite remesh candidates") {
  Mesh mesh(2, 50, true, "major");
  Matrix2d shear;
  shear << 1.0, 0.15, 0.0, 1.0;
  mesh.currentDeformation = shear;
  mesh.load = 0.15;
  mesh.loadSteps = 1;
  mesh.nrMinItterations = 2045;
  mesh.nrMinFunctionCalls = 3674;

  setNodeState(mesh.nodes(0, 25), Vector2d(25.0, 0.0),
               Vector2d(24.8405, -0.280886));
  setNodeState(mesh.nodes(0, 26), Vector2d(26.0, 0.0),
               Vector2d(25.8398, -0.369231));
  setNodeState(mesh.nodes(1, 25), Vector2d(25.0, 1.0),
               Vector2d(24.9216, 0.714346));
  setNodeState(mesh.nodes(1, 26), Vector2d(26.0, 1.0),
               Vector2d(25.9265, 0.625194));

  mesh.removeElementFromNodes(mesh.elements[50]);
  mesh.removeElementFromNodes(mesh.elements[51]);

  const std::array<GhostNode, 3> oldE1 = {
      makeGhostWithReference(mesh.nodes(0, 25), {-0.5, -0.5}),
      makeGhostWithReference(mesh.nodes(0, 26), {0.5, -0.5}),
      makeGhostWithReference(mesh.nodes(1, 25), {-0.5, 0.5}),
  };
  const std::array<GhostNode, 3> oldE2 = {
      makeGhostWithReference(mesh.nodes(1, 26), {0.5, 0.5}),
      makeGhostWithReference(mesh.nodes(1, 25), {-0.5, 0.5}),
      makeGhostWithReference(mesh.nodes(0, 26), {0.5, -0.5}),
  };

  mesh.createElementPair(oldE1, oldE2, 50, 51, true);
  TElement &e1 = mesh.elements[50];
  TElement &e2 = mesh.elements[51];
  const auto oe1 = e1.getAngleCo1Co2Nodes();
  const auto oe2 = e2.getAngleCo1Co2Nodes();
  const std::array<GhostNode, 3> n1a = {oe1[0], oe1[1], oe2[0]};
  const std::array<GhostNode, 3> n2a = {oe2[0], oe2[1], oe1[0]};

  TElement::EdgeFlipRemeshState state1;
  TElement::EdgeFlipRemeshState state2;
  CHECK_NOTHROW(state1 = e1.evaluateEdgeFlipRemeshState(n1a, mesh.F_P_H[50]));
  CHECK_NOTHROW(state2 = e2.evaluateEdgeFlipRemeshState(n2a, mesh.F_P_H[51]));
  CHECK(state1.valid);
  CHECK(state2.valid);
  for (const GhostNode &gn : state1.newGhostNodes) {
    CHECK(gn.referenceId.i >= 0);
  }
  for (const GhostNode &gn : state2.newGhostNodes) {
    CHECK(gn.referenceId.i >= 0);
  }
}

TEST_CASE("Logged simple-shear reconnect reproduces large post-flip element") {
  Mesh mesh(2, 2, false, "minor");
  mesh.load = 0.953260;
  mesh.loadSteps = 80409;
  mesh.nrMinItterations = 162;
  mesh.nrMinFunctionCalls = 249;

  // Four-node patch extracted from:
  // LargeDetReconnect_after_reconnect_cycle1_e2992_before_*.
  //
  // Local ids:
  //   0 -> original refId 1496
  //   1 -> original refId 1448
  //   2 -> original refId 1497
  //   3 -> original refId 1644
  setNodeState(mesh.nodes(0), Vector2d(0.0, 0.0),
               Vector2d(74.5232842106733, 32.8258271758159));
  setNodeState(mesh.nodes(1), Vector2d(1.0, 0.0),
               Vector2d(75.4113500667928, 33.1159452906610));
  setNodeState(mesh.nodes(2), Vector2d(0.0, 1.0),
               Vector2d(75.4436314647575, 34.2432711760344));
  setNodeState(mesh.nodes(3), Vector2d(1.0, 1.0),
               Vector2d(74.1455855272479, 33.4291086817798));

  mesh.removeElementFromNodes(mesh.elements[0]);
  mesh.removeElementFromNodes(mesh.elements[1]);

  // Pre-reconnect cells from the captured VTU:
  //   original e2894: [1496, 1448, 1497]
  //   original e2992: [1496, 1497, 1644]
  const std::array<GhostNode, 3> e1 = {
      makeGhostWithReference(mesh.nodes(0),
                             {74.8294241033744, 33.0286644279585}),
      makeGhostWithReference(mesh.nodes(1),
                             {75.8239459987426, 33.1331928912261}),
      makeGhostWithReference(mesh.nodes(2),
                             {74.7248956401067, 34.0231863233267}),
  };
  const std::array<GhostNode, 3> e2 = {
      makeGhostWithReference(mesh.nodes(0),
                             {74.5049429084062, 33.0721647565062}),
      makeGhostWithReference(mesh.nodes(2),
                             {75.4446355291921, 33.4141848998319}),
      makeGhostWithReference(mesh.nodes(3),
                             {74.1629227650805, 34.0118573772921}),
  };

  mesh.createElementPair(e1, e2, 0, 1, true);
  mesh.markDirty();
  mesh.ensureGeometry();

  CHECK(mesh.elements[0].G.determinant() ==
        doctest::Approx(0.98362).epsilon(1e-4));
  CHECK(mesh.elements[1].G.determinant() ==
        doctest::Approx(1.1894).epsilon(1e-4));

  saveCurrentAndReference(mesh, "LargeDetReconnectBefore");
  REQUIRE(mesh.reconnect());
  saveCurrentAndReference(mesh, "LargeDetReconnectAfter");

  const std::vector<std::array<int, 3>> expectedConnectivity = {{0, 1, 3},
                                                                {1, 2, 3}};
  CHECK(triConnectivity(mesh) == expectedConnectivity);
  CHECK(mesh.elements[1].G.determinant() > 2.0);
  CHECK(mesh.elements[1].G.determinant() ==
        doctest::Approx(2.06508).epsilon(1e-4));
}

TEST_CASE("Check angle after reconnecting") {

  Mesh mesh(2, 2, false, "minor");

  mesh.applyTransformation(getShear(1));
  mesh.updateElements();
  mesh.updateAngles();
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
  mesh.updateAngles();
  double reconnectAngle1 = mesh.elements[0].largestAngle;
  double reconnectAngle2 = mesh.elements[1].largestAngle;
  // std::cout << mesh.elements[0] << '\n' << mesh.elements[1] << '\n';
  CHECK(reconnectAngle1 == doctest::Approx(90));
  CHECK(reconnectAngle2 == doctest::Approx(90));
  save(mesh, "AngleCheckAfterReconnect");
}

TEST_CASE("shear updated reference elements single flip") {
  Mesh mesh(2, 2, false, "minor");

  mesh.applyTransformation(getShear(1.5, 0.01));
  saveCurrentAndReference(mesh, "ShearUpdatedReferenceElementsBeforeReconnect");

  forceEdgeFlipOnFirstElementPair(mesh);
  saveCurrentAndReference(mesh, "ShearUpdatedReferenceElementsAfterReconnect");
}

TEST_CASE("Empirical simulation edge case: flipped pair can reproduce "
          "collapsed-current-state geometry") {
  // This setup was extracted from a simulation crash log. We intentionally
  // build a minimal two-element mesh with element-specific reference triangles
  // to reproduce an empirical remeshing edge case in isolation.
  //
  // With persistent reference-based ordering, this captured post-flip fixture
  // should be accepted and canonicalized into a valid CCW reference ordering
  // instead of being rejected based on the current collapsed geometry.
  Mesh mesh(2, 2, false, "minor");

  setEmpiricalEdgeCaseRealNodes(
      mesh, {Vector2d(5.98662, 1.14637), Vector2d(5.99003, 1.13878),
             Vector2d(7.35512, 2.04187), Vector2d(5.68484, 2.37218)});
  mesh.load = 0.65275;
  mesh.loadSteps = 50276;
  mesh.nrMinItterations = 168;
  mesh.nrMinFunctionCalls = 333;

  mesh.removeElementFromNodes(mesh.elements[0]);
  mesh.removeElementFromNodes(mesh.elements[1]);

  // These two element reference configurations come directly from the crash
  // debug output after the edge flip.
  const std::array<GhostNode, 3> e1 = {
      makeGhostWithReference(mesh.nodes(1), {-0.5, 0.5}),
      makeGhostWithReference(mesh.nodes(0), {-0.5, -0.5}),
      makeGhostWithReference(mesh.nodes(3), {0.5, -0.5}),
  };
  const std::array<GhostNode, 3> e2 = {
      makeGhostWithReference(mesh.nodes(2), {0.5, -0.5}),
      makeGhostWithReference(mesh.nodes(0), {-0.5, 0.5}),
      makeGhostWithReference(mesh.nodes(3), {0.5, 0.5}),
  };

  CHECK_NOTHROW(mesh.createElementPair(e1, e2, 0, 1, true));
  CHECK(tElementInitialArea(mesh.elements[0].ghostNodes) > 0.0);
  CHECK(tElementInitialArea(mesh.elements[1].ghostNodes) > 0.0);
}


TEST_CASE("shear updated reference elements mesh") {
  Mesh mesh(2, 3, false, "minor");

  debugReconnect(mesh, "ShearUpdatedReferenceElements");
  mesh.applyTransformation(getShear(0.4));
  debugReconnect(mesh, "ShearUpdatedReferenceElements");
  mesh.applyTransformation(getShear(0.2));
  debugReconnect(mesh, "ShearUpdatedReferenceElements");
  mesh.applyTransformation(getShear(0.4));
  debugReconnect(mesh, "ShearUpdatedReferenceElements");
  mesh.applyTransformation(getShear(0.4));
  debugReconnect(mesh, "ShearUpdatedReferenceElements");
  mesh.applyTransformation(getShear(0.1));
  debugReconnect(mesh, "ShearUpdatedReferenceElements");
}

TEST_CASE("Edge flip counting") {
  SUBCASE("Counts a single flip between steps") {
    Mesh mesh(2, 2, false, "minor");
    mesh.resetCounters(); // baseline

    mesh.applyTransformation(getShear(1));
    mesh.updateMesh();
    CHECK(mesh.reconnect());

    mesh.updateAveragesAndPlasticEvents();
    CHECK(mesh.totalEdgeFlipsInStep == 1);
    CHECK(mesh.edgeFlipsFromLastStep() == 1);
  }

  SUBCASE("Flip back within the same step does not count") {
    Mesh mesh(2, 2, false, "minor");
    mesh.resetCounters(); // baseline
    const Mesh baseline = mesh;

    mesh.applyTransformation(getShear(1));
    mesh.updateMesh();
    CHECK(mesh.reconnect());
    mesh = baseline; // back to baseline topology

    mesh.resetCounters();
    CHECK(mesh.totalEdgeFlipsInStep == 0);
    CHECK(mesh.edgeFlipsFromLastStep() == 0);
  }
}

TEST_CASE("Check reconnecting with PBC") {

  Mesh mesh(2, 2, true, "major");

  mesh.applyTransformation(getShear(1));
  save(mesh, "PBCBeforeReconnect0");
  mesh.nodes(0, 1).addDisplacement({0, 0.3});
  mesh.nodes(1, 0).addDisplacement({0, 0.3});
  mesh.markDirty();
  save(mesh, "PBCBeforeReconnect1");
  mesh.nodes(0, 1).addDisplacement({0, 0.7});
  mesh.nodes(1, 0).addDisplacement({0, 0.7});
  mesh.markDirty();
  mesh.updateAveragesAndPlasticEvents();
  save(mesh, "PBCBeforeReconnect2");
  mesh.reconnect();
  // The angle node of the first element should now be moved.
  save(mesh, "PBCAfterReconnect");
  CHECK(mesh.elements[0].getAngleNode()->pos == Vector2d{0, 1});

  // Check node-element connections
}

TEST_CASE("Adjacent elements can disagree on F_P diagnostic") {
  // This is for debugging and experimenting.
  // I conclude that it is best to not touch these elements and only
  // reconnect them when they have the same F_P value.

  return;
  const Vector2d A_ref{1.0, 1.0};
  const Vector2d B_ref{1.0, 0.0};
  const Vector2d C_ref{0.0, 1.0};
  const Vector2d D_ref{2.0, 0.0};

  const Vector2d A0{2, 1.49};
  const Vector2d B0{1.0, 0.0};
  const Vector2d C0{1.0, 1.0};
  const Vector2d D0{2.0, 0.5};

  const std::array<Vector2d, 3> ref1 = {A_ref, B_ref, C_ref};
  const std::array<Vector2d, 3> ref2 = {B_ref, D_ref, A_ref};
  const std::array<Vector2d, 3> cur1 = {A0, B0, C0};
  const std::array<Vector2d, 3> cur2 = {B0, D0, A0};
  const TElement e1 = TElement(cur1, ref1);
  const TElement e2 = TElement(cur2, ref2);
  const Matrix2d &F_P1 = e1.F_P;
  const Matrix2d &F_P2 = e2.F_P;

  std::cout << "\nF1:\n" << e1.F << "\nF2:\n" << e2.F;
  std::cout << "\nF_P1:\n" << F_P1 << "\nF_P2:\n" << F_P2;
  std::cout << "\nF_E1:\n" << e1.F_E << "\nF_E2:\n" << e2.F_E;

  if (F_P1 != F_P2) {
    INFO("\nF1:\n" << e1.F << "\nF2:\n" << e2.F);
    INFO("\nF_P1:\n" << F_P1 << "\nF_P2:\n" << F_P2);
    INFO("\nF_E1:\n" << e1.F_E << "\nF_E2:\n" << e2.F_E);
    CHECK_FALSE(F_P1 != F_P2);
    return;
  }
}

// Helper function to perform a mesh operation sequence.
void performMeshOperation(Mesh &mesh, double firstParam, double secondParam,
                          const Vector2d &direction, const std::string &label) {
  mesh.moveMeshSection(firstParam, secondParam, direction, true, true);
  mesh.updateAveragesAndPlasticEvents();
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

// This test is disabled because it is slow.
// TEST_CASE("Check multiple reconnecting with saving and loading") {
//   // Create a simple config
//   Config testConfig;
//   testConfig.setDefaultValues();
//   testConfig.rows = 20;
//   testConfig.cols = 20;
//   testConfig.loadIncrement = 0.7;
//   testConfig.scenario = "doubleDislocationTest";
//   testConfig.maxLoad = 2;
//   testConfig.reconnectingEnabled = true;
//   testConfig.name = "4x4PBCLoadingTestWithReconnectingSaveLoad";
//   testConfig.logDuringMinimization = true;

//   // Create a data path and file paths
//   std::string dataPath = "test_data/";
//   std::string dumpPath = dataPath + testConfig.name +
//   "/dumps/dump_l2.1.xml.gz";

//   // Remove old data
//   clearOutputFolder(testConfig.name, dataPath);
//   std::shared_ptr<Simulation> s =
//       std::make_shared<Simulation>(testConfig, dataPath);
//   s->mesh.fixBorderNodes();
//   s->initialize();
//   s->firstStep();

//   // Run the scenario and check CSV
//   runSimulationScenario(testConfig, dataPath, s);

//   // // Load simulation into a new object
//   // using SimPtr = std::shared_ptr<Simulation>;
//   // SimPtr loadedSim = std::make_shared<Simulation>(testConfig, dataPath);
//   // Simulation::loadSimulation(*loadedSim, dumpPath, "", dataPath, true);

//   // CHECK(loadedSim->mesh == s->mesh);
//   // if (loadedSim->mesh != s->mesh) {
//   //   std::cout << debugCompare(loadedSim->mesh, s->mesh) << std::endl;
//   // }

//   // for (TElement e : s->mesh.elements) {
//   //   if (e.C(0, 1) == 0 && e.eIndex > 0 &&
//   //       e.eIndex < s->mesh.elements.size() - 1 && e.G(0, 1) != 0) {
//   //     e.update(s->mesh);
//   //     CHECK(false);
//   //   }
//   // }

//   // Rerun
//   // runSimulationScenario(testConfig, dataPath, loadedSim);
// }

struct ElementTStepSnapshot {
  double gamma = 0.0;
  Matrix2d T47 = Matrix2d::Zero();
  Matrix2d T48 = Matrix2d::Zero();
};

static Matrix2d makeMatrix2d(double a00, double a01, double a10, double a11) {
  Matrix2d matrix;
  matrix << a00, a01, a10, a11;
  return matrix;
}

static void recordElementTStepSnapshot(Simulation &simulation, void *context) {
  auto *rows = static_cast<std::vector<ElementTStepSnapshot> *>(context);
  if (rows == nullptr) {
    throw std::runtime_error(
        "recordElementTStepSnapshot requires a valid snapshot container.");
  }

  Mesh &mesh = simulation.mesh;
  const size_t e47 = 47;
  const size_t e48 = 48;
  if (mesh.elements.size() <= e48 || mesh.F_P_H.size() <= e48) {
    throw std::runtime_error(
        "recordElementTStepSnapshot requires elements 47 and 48.");
  }

  rows->push_back({mesh.load, mesh.elements[e47].totalBranch(mesh.F_P_H[e47]),
                   mesh.elements[e48].totalBranch(mesh.F_P_H[e48])});
}

static void
checkElementTStepPattern(const std::vector<ElementTStepSnapshot> &rows) {
  static const std::array<double, 7> expectedGammas = {0.0, 0.5, 1.0, 1.5,
                                                       2.0, 2.5, 3.0};
  static const std::array<Matrix2d, 7> expectedTs = {
      makeMatrix2d(1, 0, 0, 1), makeMatrix2d(1, 0, 0, 1),
      makeMatrix2d(1, 1, 0, 1), makeMatrix2d(1, 1, 0, 1),
      makeMatrix2d(1, 1, 1, 2), makeMatrix2d(1, 1, 1, 2),
      makeMatrix2d(1, 1, 2, 3)};

  REQUIRE(rows.size() == expectedGammas.size());

  constexpr double gammaTol = 1e-12;
  constexpr double matrixTol = 1e-12;
  for (size_t i = 0; i < expectedGammas.size(); ++i) {
    INFO("step index = " << i);
    INFO("T47 = " << rows[i].T47);
    INFO("T48 = " << rows[i].T48);
    INFO("expected = " << expectedTs[i]);
    CHECK(std::abs(rows[i].gamma - expectedGammas[i]) < gammaTol);
    CHECK(rows[i].T47.isApprox(expectedTs[i], matrixTol));
    CHECK(rows[i].T48.isApprox(expectedTs[i], matrixTol));
  }
}

TEST_CASE("Generate coarse 8x8 double-dislocation inspection data") {
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 8;
  testConfig.cols = 8;
  testConfig.usingPBC = false;
  testConfig.scenario = "doubleDislocationTest";
  testConfig.reconnectionMethod = "edgeFlip";
  testConfig.reconnectEdgeLocking = true;
  testConfig.reconnectRevert = false;
  testConfig.loadIncrement = 0.5;
  testConfig.maxLoad = 3.0;
  testConfig.GP1 = 0.0;
  testConfig.GP2 = 0.0;
  testConfig.epsR = 1e-3;
  testConfig.LBFGSMaxIterations = 200;
  testConfig.logDuringMinimization = true;
  testConfig.showProgress = -1;
  testConfig.writeDumps = false;
  testConfig.forceReRun = true;
  testConfig.name = "doubleDislocation8x8Inspection";

  const std::string dataPath = "test_data";
  std::vector<ElementTStepSnapshot> rows;
  constexpr int reportElement = 48;
  std::vector<ElementStepReconnectSnapshot> reportRows;
  ElementStepReconnectLoggerContext reportLoggerContext;
  reportLoggerContext.elementIndex = reportElement;
  reportLoggerContext.rows = &reportRows;

  std::shared_ptr<Simulation> simulation =
      std::make_shared<Simulation>(testConfig, dataPath, true);
  simulation->setStepLogger(recordElementTStepSnapshot, &rows);
  simulation->setReconnectStepLogger(recordElementStepReconnectSnapshot,
                                     &reportLoggerContext);
  simulation->mesh.fixNodesInRow(0);
  simulation->mesh.fixNodesInColumn(0);
  simulation->firstStep();
  runSimulationScenario(testConfig, dataPath, simulation);

  const std::filesystem::path outputDir =
      std::filesystem::path(getOutputPath(testConfig.name, dataPath));
  const std::filesystem::path jsonPath =
      outputDir / "element_48_matrix_history.json";
  writeElementStepReconnectJson(jsonPath, reportElement,
                                "History table for element 48.", reportRows,
                                testConfig);

  checkElementTStepPattern(rows);
}

TEST_CASE("Check single reconnecting Delaunay with PBC" CGAL_TEST_SKIP) {

  Mesh mesh(3, 3, true, "major");
  mesh.applyTransformation(getShear(0.1));
  save(mesh, "SimpleDelaunayPBCBeforeReconnect");
  mesh.reconnectDelaunay();
  save(mesh, "SimpleDelaunayPBCAfterReconnect");

  // TODO make proper checks here

  CHECK(mesh.nrElements == 2 * mesh.rows * mesh.cols);
}
TEST_CASE("Check reconnecting Delaunay with PBC" CGAL_TEST_SKIP) {

  Mesh mesh(3, 3, true, "major");

  mesh.applyTransformation(getShear(0.1));
  save(mesh, "SimpleDelaunayPBCBeforeReconnect");
  for (int i = 0; i < 9; ++i) {
    mesh.applyTransformation(getShear(0.1));
    mesh.reconnectDelaunay();
    CHECK(mesh.nrElements == 2 * mesh.rows * mesh.cols);
    save(mesh, "SimpleDelaunayPBCAfterReconnect");
  }
}

TEST_CASE("reconnectDelaunay with PBC: face count, forces, "
          "non-degenerate" CGAL_TEST_SKIP) {
  Mesh m(5, 5, /*PBC=*/true, "major");

  // Make it non-trivial under PBC
  m.applyTransformation(getShear(0.1));
  m.nodes(0, 1).addDisplacement({0.0, 0.25});
  // m.nodes(4, 3).addDisplacement({0.0, -0.15});
  m.markDirty();
  m.updateMesh();

  save(m, "DelaunayBeforePBC");
  m.reconnectDelaunay();
  save(m, "DelaunayAfterPBC");
  m.applyTransformation(getShear(0.1));
  m.nodes(0, 1).addDisplacement({0.0, 0.25});
  // m.nodes(4, 3).addDisplacement({0.0, -0.15});
  m.markDirty();
  m.updateMesh();

  save(m, "DelaunayBeforePBC");
  m.reconnectDelaunay();
  save(m, "DelaunayAfterPBC");

  CHECK(m.nrElements == 2 * m.rows * m.cols);

  // total nodal force ~ 0
  Vector2d sum = {0, 0};
  for (int i = 0; i < m.nrNodes; ++i)
    sum += m.nodes(i).f;
  CHECK(sum.x() == doctest::Approx(0).epsilon(1e-10));
  CHECK(sum.y() == doctest::Approx(0).epsilon(1e-10));

  // No degenerate triangles, consistent node connectivity
  for (const TElement &e : m.elements) {
    CHECK(std::abs(e.area()) > 1e-14);
  }
  check_node_connectivity_consistency(m);
}

TEST_CASE("Compare reconnectDelaunay with edgeFlip" CGAL_TEST_SKIP) {
  Mesh m(5, 5, /*PBC=*/true, "major");

  // Make it non-trivial under PBC
  m.applyTransformation(getShear(0.5));
  m.nodes(0, 1).addDisplacement({0.0, 0.25});
  m.nodes(3, 3).addDisplacement({0.0, -0.5});
  // m.nodes(4, 3).addDisplacement({0.0, -0.15});
  m.markDirty();
  m.updateMesh();

  save(m, "beforeReconnectPBC");
  m.reconnect();
  save(m, "afterReconnectPBC_EdgeFlip");
  m.reconnectDelaunay();
  save(m, "DelaunayAfterPBC");
  m.reconnect();
  save(m, "afterReconnectPBC_EdgeFlip2");
}
