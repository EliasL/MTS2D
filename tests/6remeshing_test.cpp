#include "../src/Data/data_export.h"
#include "../src/Mesh/mesh.h"
#include "Eigen/Core"
#include "Mesh/node.h"
#include "Mesh/tElement.h"
#include "Simulation/simulation.h"
#include "run/doctest.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#if defined(IDE_LIGHTWEIGHT)
#define CGAL_TEST_SKIP *doctest::skip(true)
#else
#define CGAL_TEST_SKIP
#endif

/**
 * Save mesh state to a VTU file
 * @param mesh The mesh to save
 * @param name The name suffix for the output file
 */
void save(Mesh &mesh, std::string name) {
  const std::string dataPath = "test_data";
  static int fileNr = 0;
  mesh.loadSteps = fileNr; // loadSteps is used to name the files
  mesh.ensureFull();
  writeMeshToVtu(mesh, "reconnecting", dataPath, name);
  fileNr++;
  createCollection(getDataPath("reconnecting", dataPath),
                   getOutputPath("reconnecting", dataPath));
}

void saveCurrentAndReference(Mesh &mesh, const std::string &name) {
  const std::string dataPath = "test_data";
  static int fileNr = 0;
  mesh.loadSteps = fileNr; // loadSteps is used to name the files
  mesh.ensureFull();
  writeMeshToVtu(mesh, "reconnectingReferenceTest", dataPath, name);
  writeMeshToVtu(mesh, "reconnectingReferenceTest", dataPath, name, false,
                 VtuFieldLevel::All, "reference", "", true);
  fileNr++;
  createCollection(getDataPath("reconnectingReferenceTest", dataPath),
                   getOutputPath("reconnectingReferenceTest", dataPath), "",
                   ".vtu", {}, COLLECTIONNAME, "", "_reference");
  createCollection(getDataPath("reconnectingReferenceTest", dataPath),
                   getOutputPath("reconnectingReferenceTest", dataPath), "",
                   ".vtu", {}, "reference_collection", "_reference", "");
}

void debugReconnect(Mesh &mesh, const std::string &name) {
  saveCurrentAndReference(mesh, name + "BeforeReconnect");
  mesh.reconnect();
  saveCurrentAndReference(mesh, name + "AfterReconnect");
  std::cout << "F:\n" << mesh.elements[0].F << "\n";
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
  const std::array<GhostNode, 3> n2a = {oe2[0], oe1[0], oe2[2]};
  constexpr int nrThetaSamples = 181;
  constexpr double mu = 0.0;

  const TElement::EdgeFlipRemeshState state1 =
      e1.findBestEdgeFlipRemeshStateLinearScan(
          n1a, mesh.F_P_H[static_cast<size_t>(e1i)], nrThetaSamples, mu);
  const TElement::EdgeFlipRemeshState state2 =
      e2.findBestEdgeFlipRemeshStateLinearScan(
          n2a, mesh.F_P_H[static_cast<size_t>(e2i)], nrThetaSamples, mu);

  mesh.removeElementFromNodes(e1);
  mesh.removeElementFromNodes(e2);
  mesh.createElementPair(state1.newGhostNodes, state2.newGhostNodes, e1i, e2i,
                         true);
  mesh.F_P_H[static_cast<size_t>(e1i)] = state1.H_new;
  mesh.F_P_H[static_cast<size_t>(e2i)] = state2.H_new;
  mesh.F_P_history_list[e1i].push_back(state1.P_new);
  mesh.F_P_history_list[e2i].push_back(state2.P_new);
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
  const std::array<GhostNode, 3> n2a = {oe2[0], oe1[0], oe2[2]};

  TElement::EdgeFlipRemeshState state1;
  TElement::EdgeFlipRemeshState state2;
  CHECK_NOTHROW(state1 = e1.findBestEdgeFlipRemeshStateLinearScan(
                    n1a, mesh.F_P_H[50], 181, 0.0));
  CHECK_NOTHROW(state2 = e2.findBestEdgeFlipRemeshStateLinearScan(
                    n2a, mesh.F_P_H[51], 181, 0.0));
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

  CHECK(mesh.elements[0].G.determinant() == doctest::Approx(0.98362).epsilon(1e-4));
  CHECK(mesh.elements[1].G.determinant() == doctest::Approx(1.1894).epsilon(1e-4));

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

Matrix2d simpleHorizontalShear(double gamma) {
  Matrix2d F = Matrix2d::Identity();
  F(0, 1) = gamma;
  return F;
}

Matrix2d simpleVerticalShear(double gamma) {
  Matrix2d F = Matrix2d::Identity();
  F(1, 0) = gamma;
  return F;
}

Matrix2d areaPreservingPureShear(double amount) {
  const double lambda = 1.0 + amount;
  Matrix2d F = Matrix2d::Identity();
  F(0, 0) = lambda;
  F(1, 1) = 1.0 / lambda;
  return F;
}

std::string csvEscape(std::string text) {
  bool needsQuotes = false;
  std::string escaped;
  escaped.reserve(text.size());
  for (char c : text) {
    if (c == '"') {
      escaped += "\"\"";
      needsQuotes = true;
    } else {
      if (c == ',' || c == '\n' || c == '\r') {
        needsQuotes = true;
      }
      escaped += c;
    }
  }
  if (!needsQuotes) {
    return escaped;
  }
  return "\"" + escaped + "\"";
}

void writeMatrix2dCsvHeader(std::ostream &out, const std::string &prefix) {
  out << prefix << "_00," << prefix << "_01," << prefix << "_10," << prefix
      << "_11,";
}

void writeMatrix2dCsvValues(std::ostream &out, const Matrix2d &matrix) {
  out << matrix(0, 0) << "," << matrix(0, 1) << "," << matrix(1, 0) << ","
      << matrix(1, 1) << ",";
}

struct EdgeFlipJScanRow {
  double theta = 0.0;
  bool threw = false;
  std::string error;
  TElement::EdgeFlipRemeshState state;
};

std::vector<EdgeFlipJScanRow>
sampleEdgeFlipJ(const TElement &donor, const std::array<GhostNode, 3> &candidate,
                const Matrix2d &H_old, int thetaSamples, double mu) {
  std::vector<EdgeFlipJScanRow> rows;
  rows.reserve(static_cast<size_t>(thetaSamples));
  const double dTheta = (2.0 * M_PI) / static_cast<double>(thetaSamples - 1);
  for (int i = 0; i < thetaSamples; ++i) {
    EdgeFlipJScanRow row;
    row.theta = -M_PI + dTheta * static_cast<double>(i);
    try {
      row.state =
          donor.evaluateEdgeFlipRemeshState(candidate, H_old, row.theta, mu);
    } catch (const std::exception &ex) {
      row.threw = true;
      row.error = ex.what();
    }
    rows.push_back(row);
  }
  return rows;
}

void writeEdgeFlipJScanCsv(const std::filesystem::path &path,
                           const TElement &donor,
                           const std::array<GhostNode, 3> &candidate,
                           const Matrix2d &H_old, int thetaSamples,
                           double mu) {
  const std::vector<EdgeFlipJScanRow> rows =
      sampleEdgeFlipJ(donor, candidate, H_old, thetaSamples, mu);

  double bestJ = std::numeric_limits<double>::infinity();
  for (const EdgeFlipJScanRow &row : rows) {
    if (!row.threw && row.state.valid && row.state.J < bestJ) {
      bestJ = row.state.J;
    }
  }

  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not open edge-flip J scan CSV: " +
                             path.string());
  }
  out << std::setprecision(17);
  out << "theta,theta_degrees,theta_old,theta_old_degrees,J,elastic_jump,"
         "rotation_penalty,valid,is_best,";
  writeMatrix2dCsvHeader(out, "Q_old");
  writeMatrix2dCsvHeader(out, "Q_new");
  writeMatrix2dCsvHeader(out, "F_old");
  writeMatrix2dCsvHeader(out, "F_new");
  writeMatrix2dCsvHeader(out, "P_old");
  writeMatrix2dCsvHeader(out, "P_new");
  out << "energy_old,energy_new,";
  writeMatrix2dCsvHeader(out, "sigma_old");
  writeMatrix2dCsvHeader(out, "sigma_new");
  writeMatrix2dCsvHeader(out, "E_old");
  writeMatrix2dCsvHeader(out, "E_new");
  writeMatrix2dCsvHeader(out, "E_new_minus_E_old");
  writeMatrix2dCsvHeader(out, "H_old");
  writeMatrix2dCsvHeader(out, "H_new");
  out << "error\n";

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const Matrix2d nanMatrix = Matrix2d::Constant(nan);
  for (const EdgeFlipJScanRow &row : rows) {
    const bool evaluated = !row.threw;
    const bool valid = evaluated && row.state.valid;
    const bool isBest = valid && row.state.J == bestJ;
    const double thetaOld = evaluated ? row.state.theta_old : nan;
    const double J = evaluated ? row.state.J : nan;
    const double elasticJump = evaluated ? row.state.elastic_jump : nan;
    const double rotationPenalty = evaluated ? row.state.rotation_penalty : nan;
    const Matrix2d Q_old = evaluated ? row.state.Q_old : nanMatrix;
    const Matrix2d Q_new = evaluated ? row.state.Q_new : nanMatrix;
    const Matrix2d F_old = evaluated ? row.state.F_old : nanMatrix;
    const Matrix2d F_new = evaluated ? row.state.F_new : nanMatrix;
    const Matrix2d P_old = evaluated ? row.state.P_old : nanMatrix;
    const Matrix2d P_new = evaluated ? row.state.P_new : nanMatrix;
    const double energyOld = evaluated ? row.state.energy_old : nan;
    const double energyNew = evaluated ? row.state.energy_new : nan;
    const Matrix2d sigmaOld = evaluated ? row.state.sigma_old : nanMatrix;
    const Matrix2d sigmaNew = evaluated ? row.state.sigma_new : nanMatrix;
    const Matrix2d E_old = evaluated ? row.state.E_old : nanMatrix;
    const Matrix2d E_new = evaluated ? row.state.E_new : nanMatrix;
    const Matrix2d delta_E = evaluated ? row.state.delta_E : nanMatrix;
    const Matrix2d H_old_row = evaluated ? row.state.H_old : nanMatrix;
    const Matrix2d H_new = evaluated ? row.state.H_new : nanMatrix;
    out << row.theta << "," << row.theta * 180.0 / M_PI << "," << thetaOld
        << "," << thetaOld * 180.0 / M_PI << "," << J << "," << elasticJump
        << "," << rotationPenalty << "," << valid << "," << isBest << ",";
    writeMatrix2dCsvValues(out, Q_old);
    writeMatrix2dCsvValues(out, Q_new);
    writeMatrix2dCsvValues(out, F_old);
    writeMatrix2dCsvValues(out, F_new);
    writeMatrix2dCsvValues(out, P_old);
    writeMatrix2dCsvValues(out, P_new);
    out << energyOld << "," << energyNew << ",";
    writeMatrix2dCsvValues(out, sigmaOld);
    writeMatrix2dCsvValues(out, sigmaNew);
    writeMatrix2dCsvValues(out, E_old);
    writeMatrix2dCsvValues(out, E_new);
    writeMatrix2dCsvValues(out, delta_E);
    writeMatrix2dCsvValues(out, H_old_row);
    writeMatrix2dCsvValues(out, H_new);
    out << csvEscape(row.error) << "\n";
  }
}

void writeEdgeFlipScenarioMetadata(const std::filesystem::path &path,
                                   const Matrix2d &F, int thetaSamples,
                                   double mu) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not open edge-flip metadata CSV: " +
                             path.string());
  }
  out << std::setprecision(17);
  out << "key,value\n";
  out << "theta_samples," << thetaSamples << "\n";
  out << "mu," << mu << "\n";
  out << "deformation_00," << F(0, 0) << "\n";
  out << "deformation_01," << F(0, 1) << "\n";
  out << "deformation_10," << F(1, 0) << "\n";
  out << "deformation_11," << F(1, 1) << "\n";
}

void writeLoggedElementGeometryCsv(const std::filesystem::path &path,
                                   const std::array<Vector2d, 3> &currentNodes,
                                   const std::array<Vector2d, 3> &referenceNodes,
                                   const std::array<int, 3> &sourceRefIds) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not open logged element geometry CSV: " +
                             path.string());
  }

  out << std::setprecision(17);
  out << "node,source_ref_index,current_x,current_y,reference_x,reference_y\n";
  for (int i = 0; i < 3; ++i) {
    out << i << "," << sourceRefIds[static_cast<size_t>(i)] << ","
        << currentNodes[static_cast<size_t>(i)].x() << ","
        << currentNodes[static_cast<size_t>(i)].y() << ","
        << referenceNodes[static_cast<size_t>(i)].x() << ","
        << referenceNodes[static_cast<size_t>(i)].y() << "\n";
  }
}

TEST_CASE("Export edge flip J(theta) scans" * doctest::skip(false)) {
  struct Scenario {
    std::string name;
    Matrix2d deformation;
  };

  const std::array<Scenario, 6> scenarios = {{
      {"integer_vertical_shear", simpleVerticalShear(1.0)},
      {"integer_horizontal_shear", simpleHorizontalShear(1.0)},
      {"double_integer_horizontal_shear", simpleHorizontalShear(2.0)},
      {"half_horizontal_shear", simpleHorizontalShear(0.5)},
      {"half_pure_shear", areaPreservingPureShear(0.5)},
      {"one_pure_shear", areaPreservingPureShear(1.0)},
  }};

  constexpr int thetaSamples = 361;
  constexpr double mu = 0.0;
  const std::filesystem::path root =
      std::filesystem::path("test_data") / "edge_flip_j_scan";
  std::filesystem::create_directories(root);

  for (const Scenario &scenario : scenarios) {
    Mesh mesh(2, 2, false, "minor");
    mesh.applyTransformation(scenario.deformation);
    mesh.ensureGeometry();

    TElement &e1 = mesh.elements[0];
    TElement &e2 = mesh.elements[1];
    const auto oe1 = e1.getAngleCo1Co2Nodes();
    const auto oe2 = e2.getAngleCo1Co2Nodes();

    const std::array<GhostNode, 3> n1a = {oe1[0], oe1[1], oe2[0]};
    const std::array<GhostNode, 3> n2a = {oe2[0], oe1[0], oe2[2]};
    const std::array<GhostNode, 3> n1b = {oe1[0], oe2[0], oe1[2]};
    const std::array<GhostNode, 3> n2b = {oe2[0], oe2[1], oe1[0]};

    const std::filesystem::path scenarioDir = root / scenario.name;
    std::filesystem::create_directories(scenarioDir);
    writeEdgeFlipScenarioMetadata(scenarioDir / "metadata.csv",
                                  scenario.deformation, thetaSamples, mu);
    writeEdgeFlipJScanCsv(scenarioDir / "option_a_element_1.csv", e1, n1a,
                          mesh.F_P_H[0], thetaSamples, mu);
    writeEdgeFlipJScanCsv(scenarioDir / "option_a_element_2.csv", e2, n2a,
                          mesh.F_P_H[1], thetaSamples, mu);
    writeEdgeFlipJScanCsv(scenarioDir / "option_b_element_1.csv", e1, n1b,
                          mesh.F_P_H[0], thetaSamples, mu);
    writeEdgeFlipJScanCsv(scenarioDir / "option_b_element_2.csv", e2, n2b,
                          mesh.F_P_H[1], thetaSamples, mu);
  }

  {
    const std::array<int, 3> sourceRefIds = {2011, 2012, 1962};
    const std::array<Vector2d, 3> currentNodes = {
        Vector2d(20.2799, 39.9861),
        Vector2d(21.2521, 39.8640),
        Vector2d(20.1652, 38.9713),
    };
    Matrix2d loggedF;
    loggedF << 0.89403, 1.1521, -0.1841, 0.88195;

    Matrix2d D_current;
    D_current.col(0) = currentNodes[1] - currentNodes[0];
    D_current.col(1) = currentNodes[2] - currentNodes[0];
    const Matrix2d D_reference = loggedF.inverse() * D_current;

    const std::array<Vector2d, 3> referenceNodes = {
        Vector2d::Zero(),
        D_reference.col(0),
        D_reference.col(1),
    };

    TElement element3923(currentNodes, referenceNodes);
    element3923.eIndex = 3923;

    CHECK(element3923.F(0, 0) == doctest::Approx(loggedF(0, 0)).epsilon(1e-4));
    CHECK(element3923.F(0, 1) == doctest::Approx(loggedF(0, 1)).epsilon(1e-4));
    CHECK(element3923.F(1, 0) == doctest::Approx(loggedF(1, 0)).epsilon(1e-4));
    CHECK(element3923.F(1, 1) == doctest::Approx(loggedF(1, 1)).epsilon(1e-4));

    const std::filesystem::path scenarioDir = root / "logged_element_3923";
    std::filesystem::create_directories(scenarioDir);
    writeEdgeFlipScenarioMetadata(scenarioDir / "metadata.csv", loggedF,
                                  thetaSamples, mu);
    writeLoggedElementGeometryCsv(scenarioDir / "geometry.csv", currentNodes,
                                  referenceNodes, sourceRefIds);
    writeEdgeFlipJScanCsv(scenarioDir / "self_element.csv", element3923,
                          element3923.ghostNodes, Matrix2d::Identity(),
                          thetaSamples, mu);
  }
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
  // With the current CCW invariant for newly created elements, this particular
  // collapsed fixture is no longer accepted: representing the logged negative-F
  // state would require opposite current/reference orientation at creation
  // time, so construction should fail loudly instead of silently reordering it.
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

  CHECK_THROWS_WITH(mesh.createElementPair(e1, e2, 0, 1, true),
                    doctest::Contains("counterclockwise order"));
}

TEST_CASE("Empirical simulation edge case: reconnect reproduces logged flip" *
          doctest::skip(true)) {
  // This setup was extracted from the same simulation crash log, but here we
  // recreate the state just before the edge flip and then call the real
  // reconnect logic. This is useful for stepping through empirical remeshing
  // edge cases in the debugger.
  Mesh mesh(2, 2, false, "minor");

  setEmpiricalEdgeCaseRealNodes(
      mesh, {Vector2d(6.49119, 1.04006), Vector2d(5.60101, 1.21699),
             Vector2d(6.85385, 2.06289), Vector2d(6.07056, 2.37926)});
  mesh.load = 0.65275;
  mesh.loadSteps = 50276;
  mesh.nrMinItterations = 168;
  mesh.nrMinFunctionCalls = 332;

  mesh.removeElementFromNodes(mesh.elements[0]);
  mesh.removeElementFromNodes(mesh.elements[1]);

  // These reference triangles come directly from the pre-flip element debug
  // output in the simulation log.
  const std::array<GhostNode, 3> e1 = {
      makeGhostWithReference(mesh.nodes(0), {-0.5, -0.5}),
      makeGhostWithReference(mesh.nodes(1), {-0.5, 0.5}),
      makeGhostWithReference(mesh.nodes(2), {0.5, -0.5}),
  };
  const std::array<GhostNode, 3> e2 = {
      makeGhostWithReference(mesh.nodes(3), {0.5, 0.5}),
      makeGhostWithReference(mesh.nodes(1), {-0.5, 0.5}),
      makeGhostWithReference(mesh.nodes(2), {0.5, -0.5}),
  };

  mesh.createElementPair(e1, e2, 0, 1, true);
  mesh.markDirty();

  mesh.reconnect();

  const std::vector<std::array<int, 3>> expectedConnectivity = {{0, 1, 3},
                                                                {0, 2, 3}};
  CHECK(triConnectivity(mesh) == expectedConnectivity);

  // Replace the node positions with the logged state from the following
  // minimizer evaluation. This reproduces the same post-flip collapse path
  // seen in the simulation.
  setEmpiricalEdgeCaseRealNodes(
      mesh, {Vector2d(5.98662, 1.14637), Vector2d(5.99003, 1.13878),
             Vector2d(7.35512, 2.04187), Vector2d(5.68484, 2.37218)});
  mesh.nrMinFunctionCalls = 333;
  mesh.markDirty();

  CHECK_THROWS_WITH(mesh.updateElements(),
                    doctest::Contains("Reduction exploded in "
                                      "Mesh::updateElementsForces"));
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
