#include "mesh.h"
#include "Data/data_export.h"
#include "Mesh/node.h"
#include "Mesh/tElement.h"
#include "Simulation/randomUtils.h"
#include <Eigen/Core>
#include <Eigen/LU>
#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <omp.h>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
using Eigen::Vector2i;
// CGAL is a heavy library and slows down clangd and other tools a lot.
// We use a seperate build without CGAL.
#if !defined(IDE_LIGHTWEIGHT)
#include <CGAL/Delaunay_triangulation_2.h>
// I had some issues with inexact constructions giving different triangulations
// where they should be the same.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/centroid.h>

using K = CGAL::Exact_predicates_exact_constructions_kernel; // exact predicates
                                                             // & constructions
using Point = K::Point_2;
using Delaunay = CGAL::Delaunay_triangulation_2<K>;
// Add this just above the CGAL aliases
struct VInfo {
  int refNodeIndex; // index of the reference node in the MTS mesh
                    // (0..nrNodes-1)
};
using Vb = CGAL::Triangulation_vertex_base_with_info_2<VInfo, K>;
using Fb = CGAL::Triangulation_face_base_2<K>;
using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using DelaunayInfo = CGAL::Delaunay_triangulation_2<K, Tds>;
#endif

Mesh::Mesh() {}

// Constructor that initializes the surface with size n x m
Mesh::Mesh(int rows, int cols, double a, double QDSD, bool usingPBC,
           std::string diagonal, std::string energyFunction, double bulkModulus)
    : nodes(rows, cols), a(a), rows(rows), cols(cols), loadSteps(0),
      currentDeformation(Eigen::Matrix2d::Identity()), QDSD(QDSD),
      energyFunction(energyFunction), bulkModulus(bulkModulus),
      usingPBC(usingPBC), diagonal(diagonal) {
  // Calculate nrElements based on whether usingPBC is true
  if (usingPBC) {
    nrElements = 2 * rows * cols;
  } else {
    nrElements = 2 * (rows - 1) * (cols - 1);
  }
  nrNodes = rows * cols;

  // Now initialize elements with the calculated size
  elements.resize(nrElements);
  F_P_history_list.resize(nrElements);
  F_P_H.resize(nrElements, Matrix2d::Identity());

  updateLatticeBasis();
  m_createNodes();
  m_updateFixedAndFreeNodeIds();

  // we create some elements here, but note that these will need to be replaced
  // if changes are made to the fixed and free status of the nodes.
  createElements();
  updateAveragesAndPlasticEvents();
  resetCounters();
}

Mesh::Mesh(int rows, int cols, bool usingPBC, std::string diagonal)
    : Mesh(rows, cols, 1, 0, usingPBC, diagonal) {}
Mesh::Mesh(int rows, int cols, bool usingPBC)
    : Mesh(rows, cols, usingPBC, "major") {}

void Mesh::updateLatticeBasis() {
  if (energyFunction == "conti_triangular" ||
      energyFunction == "contiTriangular") {
    const double h = std::sqrt(3.0) * 0.5 * a;
    latticeBasis << a, 0.5 * a, 0.0, h;
  } else {
    latticeBasis = Matrix2d::Identity() * a;
  }
  markDirty();
}

bool Mesh::isFixedNode(const NodeId &nodeId) const {
  return (*this)[nodeId]->fixedNode;
}

void Mesh::addLoad(double loadChange) {
  load += loadChange;
  loadSteps++;
}

void Mesh::applyTransformation(const Matrix2d &transformation) {
  // We get all the nodes in the mesh.
  for (long i = 0; i < nodes.size(); i++) {
    transformInPlace(transformation, nodes(i));
  }
  // We also assume we want to transform the current deformation
  applyTransformationToSystemDeformation(transformation);
  markDirty();
}

void Mesh::applyTransformationToFixedNodes(const Matrix2d &transformation) {
  // We get the id of each node in the border
  for (NodeId &nodeId : fixedNodeIds) {
    transformInPlace(transformation, *(*this)[nodeId]);
  }
  markDirty();
}

void Mesh::applyTransformationToSystemDeformation(
    const Matrix2d &transformation) {
  currentDeformation = transformation * currentDeformation;
  markDirty();
}

void Mesh::applyTranslation(const Vector2d &displacement) {
  for (long i = 0; i < nodes.size(); i++) {
    translateInPlace(nodes(i), displacement);
  }
  markDirty();
}

void Mesh::setRefConfiguration() {
  for (long i = 0; i < nodes.size(); i++) {
    nodes(i).setRefPos(nodes(i).pos());
  }
  for (TElement &e : elements) {
    e.setReferenceElementFromCurrentState(*this);
  }
  referenceDeformation = currentDeformation;
  for (auto &history : F_P_history_list) {
    history.clear();
  }
  for (Matrix2d &H : F_P_H) {
    H = Matrix2d::Identity();
  }
  markDirty();
}

static inline int wrap_index(int i, int n) {
  int r = i % n;
  return r < 0 ? r + n : r;
}

Node *Mesh::getNode(const Vector2i &rowCol) {
  int x = rowCol.x(); // col
  int y = rowCol.y(); // row
  // Always assume PBC
  x = wrap_index(x, cols);
  y = wrap_index(y, rows);

  // nodes(row, col) -> (y, x)
  return &nodes(y, x);
}

GhostNode Mesh::getGhostNode(const Vector2i &rowCol) {
  int x = rowCol.x(); // col
  int y = rowCol.y(); // row
  // Always assume PBC
  int ref_x = wrap_index(x, cols);
  int ref_y = wrap_index(y, rows);

  // nodes(row, col) -> (y, x)
  Node *n = &nodes(ref_y, ref_x);
  return m_gn(n, y, x);
}

Node *Mesh::getNeighbourNode(const Node &node, const Vector2i &direction) {
  int x = node.id.col() + direction.x(); // col
  int y = node.id.row() + direction.y(); // row
  if (!usingPBC) {
    // Check bounds
    if (x < 0 || x >= cols || y < 0 || y >= rows) {
      throw std::out_of_range(
          "Attempted to access neighbour node out of mesh bounds.");
    }
  }
  return getNode(Vector2i{x, y});
}

// Function to fix the elements of the border vector
void Mesh::fixBorderNodes() {
  fixNodesInRow(0);
  fixNodesInColumn(0);
  fixNodesInRow(rows - 1);
  fixNodesInColumn(cols - 1);
  m_updateFixedAndFreeNodeIds();
}

void Mesh::fixBottomLeftCorner() {
  // This is useful in PBC to avoid translation (and maybe rotation?)

  // Fix the bottom left corner node
  nodes(0, 0).fixedNode = true;
  // Update the fixed and free node IDs
  m_updateFixedAndFreeNodeIds();
}

void Mesh::fixNodesInRow(int row) {
  // Allow for negative indexing
  if (row < 0) {
    row = cols + row;
  }
  for (int col = 0; col < cols; ++col) {
    nodes(row, col).fixedNode = true;
  }
  m_updateFixedAndFreeNodeIds();
}

void Mesh::fixNodesInColumn(int col) {
  // Allow for negative indexing
  if (col < 0) {
    col = cols + col;
  }
  for (int row = 0; row < rows; ++row) {
    nodes(row, col).fixedNode = true;
  }
  m_updateFixedAndFreeNodeIds();
}

void Mesh::m_updateFixedAndFreeNodeIds() {
  fixedNodeIds.clear();
  freeNodeIds.clear();
  for (long i = 0; i < nodes.size(); i++) {
    NodeId nodeId(i, cols);
    if (isFixedNode(nodeId)) {
      fixedNodeIds.push_back(nodeId);
    } else {
      freeNodeIds.push_back(nodeId);
    }
  }
}

void Mesh::m_createNodes() {
  int n = rows;
  int m = cols;

  for (int row = 0; row < n; ++row) {
    for (int col = 0; col < m; ++col) {
      Vector2d pos = latticeBasis * Vector2d(col, row);
      Node node(pos.x(), pos.y());
      node.id = NodeId(row, col, m);
      nodes(row, col) = node;
    }
  }
}

// Gets the four nodes and their ghost versions for a given row and column in
// the reference state. NOT in the current state.
std::array<Node *, 4> Mesh::getSquareNodes(int row, int col) {
  // We find the 4 nodes in the current square
  Node *n1 = (*this)[m_makeNId(row, col)];
  Node *n2 = getNeighbourNode(*n1, RIGHT_N);
  Node *n3 = getNeighbourNode(*n1, UP_N);
  // n4 is now up AND right of n1
  Node *n4 = getNeighbourNode(*n2, UP_N);

  // The nodes should be in this square configuration in the reference state
  // n3  n4
  // n1  n2

  // This function will find nodes across the system if necessary due to
  // periodic boundaries
  return {n1, n2, n3, n4};
}

// Gets the four nodes and their ghost versions for a given row and column
std::array<GhostNode, 4> Mesh::getSquareGhostNodes(int row, int col) {
  auto nodes = getSquareNodes(row, col);
  return m_makeGhostNodes(nodes, row, col);
}

// Calculates element indices for a given row and column
std::pair<int, int> Mesh::getElementIndices(int row, int col) {
  int e1i = 2 * (row * ePairCols() + col); // Triangle 1 index
  int e2i = e1i + 1;                       // Triangle 2 index
  return {e1i, e2i};
}

/**
 * Creates or updates two triangular elements based on the specified diagonal
 * direction
 * @param ghosts array of four ghost nodes {g1, g2, g3, g4}
 * @param e1i Index of first element
 * @param e2i Index of second element
 * @param useLeftDiagonal Whether to use left diagonal splitting
 * @param preserveNoise Whether to preserve existing noise values (true for
 * updates)
 */
void Mesh::createElementPair(const std::array<GhostNode, 4> &ghosts, int e1i,
                             int e2i, bool majorDiagonalOrder,
                             bool preserveNoise) {
  if (ghosts.size() != 4) {
    throw std::runtime_error("Expected 4 ghost nodes!");
  }
  const std::array<const GhostNode *, 4> g = {&ghosts[0], &ghosts[1],
                                              &ghosts[2], &ghosts[3]};
  createElementPair(g, e1i, e2i, majorDiagonalOrder, preserveNoise);
}

void Mesh::createElementPair(const std::array<const GhostNode *, 4> &ghostsPtr,
                             int e1i, int e2i, bool majorDiagonalOrder,
                             bool preserveNoise) {
  const std::array<GhostNode, 4> g = {*ghostsPtr[0], *ghostsPtr[1],
                                      *ghostsPtr[2], *ghostsPtr[3]};
  if (majorDiagonalOrder) {
    const std::array<GhostNode, 3> e1 = {g[0], g[1], g[2]};
    const std::array<GhostNode, 3> e2 = {g[3], g[1], g[2]};
    createElementPair(e1, e2, e1i, e2i, preserveNoise);
  } else {
    const std::array<GhostNode, 3> e1 = {g[1], g[0], g[3]};
    const std::array<GhostNode, 3> e2 = {g[2], g[0], g[3]};
    createElementPair(e1, e2, e1i, e2i, preserveNoise);
  }
}

void Mesh::createElementPair(const std::array<GhostNode, 3> &e1,
                             const std::array<GhostNode, 3> &e2, int e1i,
                             int e2i, bool preserveNoise) {
  double noise1, noise2;
  if (preserveNoise) {
    noise1 = elements[e1i].noise;
    noise2 = elements[e2i].noise;
  } else {
    noise1 = sampleNormal(1, QDSD);
    noise2 = sampleNormal(1, QDSD);
  }

  const std::array<GhostNode, 3> e1Copy = e1;
  const std::array<GhostNode, 3> e2Copy = e2;
  elements[e1i] = TElement((*this), e1Copy[0], e1Copy[1], e1Copy[2], e1i,
                           noise1, energyFunction, bulkModulus);
  elements[e2i] = TElement((*this), e2Copy[0], e2Copy[1], e2Copy[2], e2i,
                           noise2, energyFunction, bulkModulus);
}

void Mesh::createElementPair(const std::array<const GhostNode *, 3> &e1,
                             const std::array<const GhostNode *, 3> &e2,
                             int e1i, int e2i, bool preserveNoise) {
  const std::array<GhostNode, 3> e1Copy = {*e1[0], *e1[1], *e1[2]};
  const std::array<GhostNode, 3> e2Copy = {*e2[0], *e2[1], *e2[2]};
  createElementPair(e1Copy, e2Copy, e1i, e2i, preserveNoise);
}
void centerReferencePairAtOrigin(std::array<GhostNode, 4> &pairGhosts) {
  // Reference translations do not affect F/C/forces, so we keep each
  // created element pair in one common local frame by translating the 4-node
  // pair so the midpoint of the longest reference diagonal lies at the
  // origin.
  double maxDist2 = -1.0;
  int iMax = 0;
  int jMax = 1;
  for (int i = 0; i < 4; ++i) {
    for (int j = i + 1; j < 4; ++j) {
      const double dist2 =
          (pairGhosts[i].ref_pos - pairGhosts[j].ref_pos).squaredNorm();
      if (dist2 > maxDist2) {
        maxDist2 = dist2;
        iMax = i;
        jMax = j;
      }
    }
  }
  const Vector2d pairCenter =
      0.5 * (pairGhosts[iMax].ref_pos + pairGhosts[jMax].ref_pos);
  for (GhostNode &gn : pairGhosts) {
    gn.updateReferencePosition(gn.ref_pos - pairCenter);
  }
}

void Mesh::createElements() {
  // Note that neighbours must be filled before using this function.

  // We construct the elements by finding the four nodes that create two
  // opposing cells, but stopping before we get to the last row and columns if
  // we don't use periodic boundary conditions.
  for (int row = 0; row < ePairRows(); ++row) {
    for (int col = 0; col < ePairCols(); ++col) {
      auto ghosts = getSquareGhostNodes(row, col);
      auto [e1i, e2i] = getElementIndices(row, col);

      // Determine diagonal direction based on alternating pattern

      bool majorDiagonalOrder;
      if (diagonal == "major") {
        majorDiagonalOrder = true;
      } else if (diagonal == "minor") {
        majorDiagonalOrder = false;
      } else if (diagonal == "alternate") {
        majorDiagonalOrder = (row + col) % 2;
      } else {
        throw std::invalid_argument("Unkown meshing: " + diagonal);
      }
      centerReferencePairAtOrigin(ghosts);
      createElementPair(ghosts, e1i, e2i, majorDiagonalOrder, false);
    }
  }
}

// This is just a function to avoid having to write cols
NodeId Mesh::m_makeNId(int row, int col) { return NodeId(row, col, cols); }

namespace {
constexpr double kMatrixCompareTol = 1e-9;

std::string formatGhostNodeDebug(const GhostNode &gn) {
  std::ostringstream oss;
  oss << "refId=" << gn.referenceId.i << " id=(" << gn.id.x() << ","
      << gn.id.y() << ") pShift=(" << gn.periodicShift.x() << ","
      << gn.periodicShift.y() << ") pos=(" << gn.pos.x() << ", " << gn.pos.y()
      << ") ref=(" << gn.ref_pos.x() << ", " << gn.ref_pos.y() << ") u=("
      << gn.u.x() << ", " << gn.u.y() << ")";
  return oss.str();
}

std::string
formatGhostNodeGroupDebug(const std::string &groupName,
                          const std::array<GhostNode, 3> &ghostNodes) {
  std::ostringstream oss;
  oss << groupName << ":\n";
  for (int i = 0; i < 3; ++i) {
    oss << groupName << "[" << i << "]: " << formatGhostNodeDebug(ghostNodes[i])
        << "\n";
  }
  return oss.str();
}

std::string formatVector2dDebug(const Vector2d &v) {
  std::ostringstream oss;
  oss << "(" << v.x() << ", " << v.y() << ")";
  return oss.str();
}

std::string sanitizeDebugFileComponent(const std::string &text) {
  std::string result = text;
  for (char &c : result) {
    const unsigned char uc = static_cast<unsigned char>(c);
    if (!(std::isalnum(uc) || c == '_' || c == '-')) {
      c = '_';
    }
  }
  return result;
}

std::string formatElementReductionDebug(const TElement &e,
                                        bool includeGhostNodes = true) {
  std::ostringstream oss;
  oss << "eIndex: " << e.eIndex << "\n"
      << "m3Nr: " << e.m3Nr << " red_quadrant: " << e.red_quadrant << "\n"
      << "F:\n"
      << e.F << "\n"
      << "F_P:\n"
      << e.F_P << "\n"
      << "F_E:\n"
      << e.F_E << "\n"
      << "C:\n"
      << e.C << "\n"
      << "C_R:\n"
      << e.C_R << "\n"
      << "G:\n"
      << e.G << "\n"
      << "M_e:\n"
      << e.M_e << "\n"
      << "M_l:\n"
      << e.M_l << "\n"
      << "P:\n"
      << e.P << "\n"
      << "sigma:\n"
      << e.sigma << "\n";
  if (includeGhostNodes) {
    for (int i = 0; i < 3; ++i) {
      oss << "ghost[" << i << "]: " << formatGhostNodeDebug(e.ghostNodes[i])
          << "\n";
    }
  }
  return oss.str();
}

std::string formatEdgeFlipCandidateFailure(
    const Mesh &mesh, const TElement &e1, const TElement &e2,
    const std::array<GhostNode, 3> &oldGhosts1,
    const std::array<GhostNode, 3> &oldGhosts2,
    const std::array<GhostNode, 3> &candidate1,
    const std::array<GhostNode, 3> &candidate2, const Matrix2d &history1,
    const Matrix2d &history2, const std::string &phase,
    const std::string &cause) {
  std::ostringstream oss;
  oss << "Mesh::flipEdge failed.\n\n"
      << "phase: " << phase << "\n\n"
      << "cause:\n"
      << cause << "\n\n"
      << "mesh state:\n"
      << "load: " << mesh.load << "\n"
      << "loadSteps: " << mesh.loadSteps << "\n"
      << "nrMinItterations: " << mesh.nrMinItterations << "\n"
      << "nrMinFunctionCalls: " << mesh.nrMinFunctionCalls << "\n"
      << "usingPBC: " << mesh.usingPBC << "\n"
      << "currentDeformation:\n"
      << mesh.currentDeformation << "\n"
      << "referenceDeformation:\n"
      << mesh.referenceDeformation << "\n\n"
      << "element 1 before flip:\n"
      << formatElementReductionDebug(e1, false)
      << formatGhostNodeGroupDebug("oldGhosts1", oldGhosts1) << "\n"
      << "element 2 before flip:\n"
      << formatElementReductionDebug(e2, false)
      << formatGhostNodeGroupDebug("oldGhosts2", oldGhosts2) << "\n"
      << "candidate1 inherits element " << e1.eIndex << ":\n"
      << formatGhostNodeGroupDebug("candidate1", candidate1) << "\n"
      << "candidate2 inherits element " << e2.eIndex << ":\n"
      << formatGhostNodeGroupDebug("candidate2", candidate2) << "\n"
      << "F_P_H[element1]:\n"
      << history1 << "\n"
      << "F_P_H[element2]:\n"
      << history2 << "\n";
  return oss.str();
}
} // namespace

void Mesh::resetCounters() {
  nrMinItterations = 0;
  nrMinFunctionCalls = 0;
  totalEdgeFlipsInStep = 0;
  edgeFlipDeltaSinceLastStep.clear();
  resetPastPlasticCount();
}

void Mesh::resetPastPlasticCount(bool endOfStep) {
  nr_elements_with_m3_changeInStep = 0;
  for (size_t i = 0; i < elements.size(); i++) {
    elements[i].pastStepM3Nr = elements[i].m3Nr;
  }

  if (endOfStep) {
    nr_elements_with_m3_change = 0;
    for (size_t i = 0; i < elements.size(); i++) {
      elements[i].pastM3Nr = elements[i].m3Nr;
    }
  }
}

void Mesh::setSimNameAndDataPath(std::string name, std::string path) {
  simName = name;
  dataPath = path;
}

// Helper function to make ghost node

// This function automatically sets the periodic shift based on row and col
GhostNode Mesh::m_gn(const Node *n, int row, int col) {
  return GhostNode(n, row, col, cols, latticeBasis, currentDeformation,
                   referenceDeformation);
}
// This function automatically sets the periodic shift based on targetPos
GhostNode Mesh::m_gn(const Node *n, const Vector2d targetPos) {
  return GhostNode(n, targetPos, rows, cols, latticeBasis, currentDeformation,
                   referenceDeformation);
}
GhostNode Mesh::m_gn(const Node *n) {
  return GhostNode(n, n->id.row(), n->id.col(), cols, latticeBasis,
                   currentDeformation, referenceDeformation);
}

// The idea here is that we have taken our grid, and made rows and columns of
// squares of 4 and 4 nodes. This function converts from reference nodes to
// ghost nodes and makes sure that the ghost nodes are appropreately shifted
// when that is requried.
// We never want our square to span the entire system which would happen at the
// boundary if we use the "real" position of the nodes
std::array<GhostNode, 4>
Mesh::m_makeGhostNodes(const std::array<Node *, 4> refNodes, int row, int col) {
  assert(refNodes.size() == 4);

  const Node *n1 = refNodes[0];
  const Node *n2 = refNodes[1];
  const Node *n3 = refNodes[2];
  const Node *n4 = refNodes[3];

  GhostNode gn1 = m_gn(n1);
  GhostNode gn2 = m_gn(n2);
  GhostNode gn3 = m_gn(n3);
  GhostNode gn4 = m_gn(n4);

  if (usingPBC) {

    if (row == rows - 1 && col == cols - 1) {
      // If we are in the corner, we need to move n2, n3 and n4
      gn2 = m_gn(n2, n2->id.row(), cols);
      gn3 = m_gn(n3, rows, n3->id.col());
      gn4 = m_gn(n4, rows, cols);
    } else if (col == cols - 1) {
      // If we are in the last column, we need to move n2 and n4
      gn2 = m_gn(n2, n2->id.row(), cols);
      gn4 = m_gn(n4, n4->id.row(), cols);
    } else if (row == rows - 1) {
      // If we are in the last row, we need to move n3 and n4
      gn3 = m_gn(n3, rows, n3->id.col());
      gn4 = m_gn(n4, rows, n4->id.col());
    }
  }

  return {gn1, gn2, gn3, gn4};
}

void Mesh::printConnectivity(bool realId) {
  std::string sep;
  std::string end;
  if (nrNodes <= 9) {
    sep = "";
    end = " ";
  } else {
    sep = ",";
    end = "\n";
  }
  for (int i = 0; i < nrElements; i++) {
    TElement &e = elements[i];
    for (size_t j = 0; j < e.ghostNodes.size(); j++) {
      if (realId) {
        std::cout << e.ghostNodes[j].referenceId.i << sep;
      } else {
        std::cout << e.ghostNodes[j].id << sep;
      }
    }
    std::cout << end;
  }
  std::cout << '\n';
}

void Mesh::throwIfReductionExploded(const TElement &element,
                                    const std::string &context) const {
  if (element.m3Nr <= 100) {
    return;
  }

  const std::string debugBaseName =
      "ReductionExploded_" + sanitizeDebugFileComponent(context);
  try {
    Mesh debugMesh = *this;
    debugMesh.refreshCurrentGhostGeometryForDebug();
    writeMeshToVtu(debugMesh, simName, dataPath, debugBaseName, true,
                   VtuFieldLevel::All, "current");
  } catch (...) {
  }
  try {
    Mesh debugMesh = *this;
    debugMesh.refreshCurrentGhostGeometryForDebug();
    writeMeshToVtu(debugMesh, simName, dataPath, debugBaseName, true,
                   VtuFieldLevel::All, "reference", "", true);
  } catch (...) {
  }

  std::ostringstream oss;
  oss << "Reduction exploded in " << context << ".\n\n"
      << "minimization:\n"
      << "nrMinItterations: " << nrMinItterations << "\n"
      << "nrMinFunctionCalls: " << nrMinFunctionCalls << "\n"
      << "load: " << load << "\n"
      << "loadSteps: " << loadSteps << "\n\n"
      << "element:\n"
      << formatElementReductionDebug(element);

  oss << "realNodeRefs:\n";
  for (int i = 0; i < 3; ++i) {
    const Node *realNode = (*this)[element.ghostNodes[i].referenceId];
    oss << "realNodeRef[" << i << "]: refId=" << realNode->id.i << " pos=("
        << realNode->pos().x() << ", " << realNode->pos().y() << ") ref=("
        << realNode->ref_pos().x() << ", " << realNode->ref_pos().y() << ") u=("
        << realNode->u().x() << ", " << realNode->u().y() << ")\n";
  }

  if (element.eIndex >= 0 &&
      element.eIndex < static_cast<int>(lastFlipDebugStates.size())) {
    const LastFlipDebugState &debugState = lastFlipDebugStates[element.eIndex];
    if (debugState.valid) {
      oss << "\nlastFlipDebug:\n"
          << "element=" << element.eIndex << " with partner "
          << debugState.partner << "\n"
          << "minIterationsAtFlip: " << debugState.minIterationsAtFlip << "\n"
          << "minFunctionCallsAtFlip: " << debugState.minFunctionCallsAtFlip
          << "\n"
          << "deltaMinIterationsSinceFlip: "
          << (nrMinItterations - debugState.minIterationsAtFlip) << "\n"
          << "deltaMinFunctionCallsSinceFlip: "
          << (nrMinFunctionCalls - debugState.minFunctionCallsAtFlip) << "\n"
          << "applied_F_P:\n"
          << debugState.applied_F_P << "\n"
          << "oldAnchor: " << formatVector2dDebug(debugState.oldAnchor) << "\n"
          << "newAnchor: " << formatVector2dDebug(debugState.newAnchor) << "\n";
      const std::array<GhostNode, 3> postFlipSelfGhosts = element.ghostNodes;
      oss << formatGhostNodeGroupDebug("postFlipSelfGhost", postFlipSelfGhosts);
      if (debugState.partner >= 0 &&
          debugState.partner < static_cast<int>(elements.size())) {
        oss << "\npartnerElement:\n"
            << formatElementReductionDebug(elements[debugState.partner], false);
        oss << formatGhostNodeGroupDebug(
            "postFlipPartnerGhost", elements[debugState.partner].ghostNodes);
      }
      oss << formatGhostNodeGroupDebug("preFlipSelfGhost",
                                       debugState.oldSelfGhostNodes);
      oss << formatGhostNodeGroupDebug("preFlipPartnerGhost",
                                       debugState.oldPartnerGhostNodes);
    }
  }

  throw std::runtime_error(oss.str());
}

// Updates the forces on the nodes in the surface and returns the total
// energy from all the elements in the surface.
void Mesh::updateMesh() { ensureForces(); }

void Mesh::updateElements() { ensureForces(); }

void Mesh::markDirty() { updateState = UpdateState::Dirty; }

void Mesh::refreshCurrentGhostGeometryForDebug() {
  for (TElement &element : elements) {
    element.refreshCurrentGhostGeometryForDebug(*this);
  }
}

std::size_t Mesh::edgeFlipsFromLastStep() const {
  if (edgeFlipDeltaSinceLastStep.size() % 2 != 0) {
    throw std::runtime_error(
        "Mesh::edgeFlipsFromLastStep: edge delta set has odd size.");
  }
  return edgeFlipDeltaSinceLastStep.size() / 2;
}

void Mesh::ensureForces() {
  if (updateState != UpdateState::Dirty) {
    return;
  }
  updateElementsForces();
  updateState = UpdateState::Forces;
}

void Mesh::ensureGeometry() {
  if (updateState == UpdateState::Geometry ||
      updateState == UpdateState::Full) {
    return;
  }
  ensureForces();
  updateElementsGeometry();
  updateState = UpdateState::Geometry;
}

void Mesh::ensureFull() {
  if (updateState == UpdateState::Full) {
    return;
  }
  ensureGeometry();
  updateElementsFull();
  updateState = UpdateState::Full;
}

/*
This function has given me a lot of headache. When multithreading, the program
will sometimes hang, with all threads waiting for each other. I have checked
that all threads do actually reach the barrier, but the function never returns.
Using a manual barrier, i can add a timeout. This seems to work.
*/
void Mesh::updateElementsForces() {
  omp_set_dynamic(0);

  const int nNodes = nodes.size();
  const int nThreads = omp_get_max_threads();
  const size_t scratchSize = static_cast<size_t>(nThreads) * nNodes;

  if (forceScratch.size() != scratchSize || forceScratchThreads != nThreads) {
    forceScratch.resize(scratchSize);
    forceScratchThreads = nThreads;
  }

  double energy_sum = 0.0;

#pragma omp parallel reduction(+ : energy_sum)
  {
    const int tid = omp_get_thread_num();
    Vector2d *local = forceScratch.data() + static_cast<size_t>(tid) * nNodes;

    // Parallel first-touch / zeroing of this thread's stripe
    for (int i = 0; i < nNodes; ++i) {
      local[i].setZero();
    }

#pragma omp for schedule(static) nowait
    for (int i = 0; i < nrElements; ++i) {
      TElement &e = elements[i];
      e.updateForces(*this);
      energy_sum += e.energy;

      const GhostNode &g0 = e.ghostNodes[0];
      const GhostNode &g1 = e.ghostNodes[1];
      const GhostNode &g2 = e.ghostNodes[2];

      local[g0.referenceId.i] += g0.f;
      local[g1.referenceId.i] += g1.f;
      local[g2.referenceId.i] += g2.f;
    }
  }

  for (int i = 0; i < nrElements; ++i) {
    throwIfReductionExploded(elements[i], "Mesh::updateElementsForces");
  }

  totalEnergy = energy_sum;

  double maxForce = 0.0;
#pragma omp parallel for reduction(max : maxForce)
  for (int i = 0; i < nNodes; ++i) {
    Vector2d sum = Vector2d::Zero();
    for (int t = 0; t < nThreads; ++t) {
      sum += forceScratch[static_cast<size_t>(t) * nNodes + i];
    }
    nodes(i).f = sum;
    if (!nodes(i).fixedNode) {
      maxForce = std::max({maxForce, std::abs(sum[0]), std::abs(sum[1])});
    }
  }

  this->maxForce = maxForce;
}

void Mesh::updateElementsGeometry() {
#pragma omp parallel for schedule(static, 1024)
  for (int i = 0; i < nrElements; ++i) {
    elements[i].updateGeometry();
  }
  for (int i = 0; i < nrElements; ++i) {
    throwIfReductionExploded(elements[i], "Mesh::updateElementsGeometry");
  }
}

void Mesh::updateElementsFull() {
#pragma omp parallel for schedule(static, 1024)
  for (int i = 0; i < nrElements; ++i) {
    elements[i].updateFull();
  }
}

void Mesh::updateNodeForce(Node &node) {
  node.resetForce();
  for (size_t e = 0; e < node.elementCount; ++e) {
    const int elIdx = node.connectedElements[e];
    const int nodeIdx = node.nodeIndexInElement[e];
    if (elIdx < 0 || elIdx >= static_cast<int>(elements.size()) ||
        nodeIdx < 0 || nodeIdx >= 3) {
      throw std::runtime_error("updateNodeForce: invalid node connectivity.");
    }
    node.f += elements[elIdx].ghostNodes[nodeIdx].f;
  }
}

void Mesh::applyForceFromElementsToNodes() {
  // Loop over all the nodes
  // (Looping over the elements would create a problem since two threads might
  // write to the same node at the same time. By looping over the nodes, each
  // thread only writes to one node at a time. There is no problem if two
  // threads read from the same element at the same time.)
#pragma omp parallel for
  for (int i = 0; i < (int)nodes.size(); ++i) {
    Node &n = nodes(i);
    updateNodeForce(n);
  }
}

void Mesh::rebuildConnectivity() {
  for (int i = 0; i < nodes.size(); ++i) {
    Node &n = nodes(i);
    n.elementCount = 0;
    n.connectedElements.fill(-1);
    n.nodeIndexInElement.fill(-1);
  }

  for (int eIdx = 0; eIdx < elements.size(); ++eIdx) {
    const TElement &e = elements[eIdx];
    for (int i = 0; i < 3; ++i) {
      const GhostNode &gn = e.ghostNodes[i];
      Node *node = (*this)[gn.referenceId];
      int &count = node->elementCount;
      if (count < MAX_ELEMENTS_PER_NODE) {
        node->connectedElements[count] = eIdx;
        node->nodeIndexInElement[count] = i;
        ++count;
      } else {
        throw std::overflow_error("Element index overflow for node " +
                                  std::to_string(gn.referenceId.i));
      }
    }
  }
}

void Mesh::checkForces(std::vector<Node *> nodes) {
  for (Node *n : nodes) {
    updateNodeForce(*n);
    if (n->f.norm() > 1e-6) {
      std::cout << "Invalid: " << n->id << ": \t" << n->f << '\n';
      std::cout << n->f.norm() << '\n';
    } else {
      std::cout << "Valid: " << n->id << ": \t" << n->f << '\n';
      std::cout << n->f.norm() << '\n';
    }
  }
}

void Mesh::checkForces(const std::vector<GhostNode> nodes) {

  checkForces({(*this)[nodes[0].referenceId], (*this)[nodes[1].referenceId],
               (*this)[nodes[2].referenceId], (*this)[nodes[3].referenceId]});
}

inline bool inRegion(const Eigen::Matrix2d &G) {
  const double a = G(0, 0);
  const double b = G(0, 1);
  assert(G(0, 1) == G(1, 0));
  const double c = G(1, 1);
  if (b < 0) {
    return false;
  }
  if (b > std::min(a, c)) {
    return false;
  }
  return true;
}

inline bool inExtendedRegion(const Eigen::Matrix2d &G) {
  const double a = G(0, 0);
  const double b = 0.5 * (G(0, 1) + G(1, 0)); // robust off-diagonal
  const double c = G(1, 1);

  if (a <= c) {
    const double b_lo = -0.5 * a;
    const double b_hi = std::min(1.5 * a, (3.0 * a + c) / 4.0);
    return (b >= b_lo) && (b <= b_hi);
  } else {
    const double b_lo = -0.5 * c;
    const double b_hi = std::min(1.5 * c, (a + 3.0 * c) / 4.0);
    return (b >= b_lo) && (b <= b_hi);
  }
}

inline double maxCosineForMinAngle(const TElement &e) {
  // We want to compare the smalest angle, but
  // avoid expensive trigonometric functions.
  double maxCos = -1.0;
  for (int i = 0; i < 3; ++i) {
    const int next = (i + 1) % 3;
    const int prev = (i + 2) % 3;
    const Vector2d v1 = e.ghostNodes[next].pos - e.ghostNodes[i].pos;
    const Vector2d v2 = e.ghostNodes[prev].pos - e.ghostNodes[i].pos;
    const double denom = v1.norm() * v2.norm();
    if (denom <= 1e-12) {
      return 1.0;
    }
    double cosA = v1.dot(v2) / denom;
    cosA = std::clamp(cosA, -1.0, 1.0);
    if (cosA > maxCos) {
      maxCos = cosA;
    }
  }
  return maxCos;
}

static Mesh::EdgeKey getSharedReferenceEdgeKey(const TElement &e1,
                                               const TElement &e2) {
  std::array<int, 2> shared = {-1, -1};
  int sharedCount = 0;
  for (const GhostNode &g1 : e1.ghostNodes) {
    for (const GhostNode &g2 : e2.ghostNodes) {
      if (g1.referenceId.i != g2.referenceId.i) {
        continue;
      }
      if (sharedCount >= 2) {
        throw std::runtime_error(
            "getSharedReferenceEdgeKey: element pair shares more than two "
            "reference nodes.");
      }
      shared[sharedCount++] = g1.referenceId.i;
      break;
    }
  }

  if (sharedCount != 2) {
    throw std::runtime_error(
        "getSharedReferenceEdgeKey: element pair does not share a common "
        "reference edge.");
  }
  return Mesh::EdgeKey(shared[0], shared[1]);
}

static void toggleEdgeKey(Mesh::EdgeSet &edgeSet, const Mesh::EdgeKey &edge) {
  auto [it, inserted] = edgeSet.insert(edge);
  if (!inserted) {
    edgeSet.erase(it);
  }
}

bool Mesh::reconnect(bool onlyCheck, EdgeSet *lockedEdges) {
  ensureGeometry();
  meshReconnected = false;
  bool changedInAnySweep = false;

  while (true) {
    bool changedThisSweep = false;
    for (int i = 0; i < elements.size(); i++) {
      TElement &e = elements[i];

      // Is the element deformed enough to be reconnected?
      if (inRegion(e.G)) {
        continue;
      }
      // Does the element have a matching twin to reconnect with?
      const int twinIndex = e.getElementTwin(*(this));
      if (twinIndex == -1) {
        continue;
      }
      // Is the twin also deformed enough to be reconnected?
      TElement &twin = elements[twinIndex];
      if (inRegion(twin.G)) {
        continue;
      }

      // if (e.F_P != twin.F_P) {
      //   std::cout << "skipping different F_P\n";
      //   continue;
      // }

      // Is the edge being flipped locked?
      const EdgeKey sharedEdge = getSharedReferenceEdgeKey(e, twin);
      if (lockedEdges != nullptr &&
          lockedEdges->find(sharedEdge) != lockedEdges->end()) {
        continue;
      }

      if (onlyCheck) {
        return true;
      }
      try {
        checkAndFixPeriodicElementPair(e, twin);
      } catch (const std::exception &ex) {
        throw std::runtime_error(formatEdgeFlipCandidateFailure(
            *this, e, twin, e.ghostNodes, twin.ghostNodes, e.ghostNodes,
            twin.ghostNodes, F_P_H[static_cast<size_t>(i)],
            F_P_H[static_cast<size_t>(twinIndex)],
            "moving periodic element pair before edge flip", ex.what()));
      } catch (...) {
        throw std::runtime_error(formatEdgeFlipCandidateFailure(
            *this, e, twin, e.ghostNodes, twin.ghostNodes, e.ghostNodes,
            twin.ghostNodes, F_P_H[static_cast<size_t>(i)],
            F_P_H[static_cast<size_t>(twinIndex)],
            "moving periodic element pair before edge flip",
            "unknown non-std exception"));
      }

      flipEdge(e, twin);
      const EdgeKey newSharedEdge =
          getSharedReferenceEdgeKey(elements[i], elements[twinIndex]);
      toggleEdgeKey(edgeFlipDeltaSinceLastStep, sharedEdge);
      toggleEdgeKey(edgeFlipDeltaSinceLastStep, newSharedEdge);
      if (lockedEdges != nullptr) {
        lockedEdges->insert(newSharedEdge);
      }

      changedThisSweep = true;
      changedInAnySweep = true;
      meshReconnected = true;
    }

    if (onlyCheck || !changedThisSweep) {
      break;
    }
  }

  if (!onlyCheck && changedInAnySweep) {
    markDirty();
  }
  return changedInAnySweep;
}

// Delauney Helpers

#if !defined(IDE_LIGHTWEIGHT)
inline Eigen::Vector2d toEigen(const Point &p) {
  return {CGAL::to_double(p.x()), CGAL::to_double(p.y())};
}
// Unique key for a triangle based on its three reference vertex indices
struct TriKey {
  uint64_t a, b, c;
  bool operator==(const TriKey &o) const noexcept {
    return a == o.a && b == o.b && c == o.c;
  }
  // Define less-than operator for use in std::set
  bool operator<(const TriKey &o) const noexcept {
    if (a != o.a)
      return a < o.a;
    if (b != o.b)
      return b < o.b;
    return c < o.c;
  }
};

static inline TriKey makeTriKey(const DelaunayInfo::Face_handle &f) {
  std::array<uint64_t, 3> v{0, 0, 0};
  for (int k = 0; k < 3; ++k) {
    v[k] = f->vertex(k)->info().refNodeIndex;
  }
  std::sort(v.begin(), v.end());
  return TriKey{v[0], v[1], v[2]};
}

void Mesh::reconnectDelaunay() {
  // We will refer to the native MTS mesh as the "MTSMesh", and the Delaunay
  // triangulation as the "Dmesh".
  // In order to deal with periodic boundaries, we will extend the mesh
  // by three node layers in each direction.
  int extension = 3;
  int extendedRows = rows + extension * 2;
  int extendedCols = cols + extension * 2;
  int nrNodesInDMesh = extendedRows * extendedCols;
  // Note that the same NodeId will appear multiple times
  std::vector<NodeId> DNodeToMTSNode; // index in Dmesh -> NodeId in MTSMesh
  DNodeToMTSNode.reserve(nrNodesInDMesh);

  // 2) Build Delaunay with vertex info = VInfo (baseId, dr, dc)
  DelaunayInfo dt;
  for (int dr = -extension; dr < rows + extension; ++dr) {
    for (int dc = -extension; dc < cols + extension; ++dc) {
      const GhostNode gn = getGhostNode(Vector2i{dc, dr});
      const Point p(gn.pos.x(), gn.pos.y());
      auto vh = dt.insert(p);
      VInfo vi;
      vi.refNodeIndex = gn.referenceId.i;
      vh->info() = vi;
    }
  }

  // 3) Remove excess elements
  // We only keep faces with two or more vertices inside the window
  // And we also make sure to only keep unique triangles (no duplicates)
  std::set<TriKey> seenTriangles;

  auto keep_element = [&](const DelaunayInfo::Face_handle &f) {
    // Deduplicate by canonical triangle key based on reference vertex indices
    TriKey key = makeTriKey(f);
    if (seenTriangles.find(key) != seenTriangles.end()) {
      return false; // duplicate triangle
    }

    const Point &p0 = f->vertex(0)->point();
    const Point &p1 = f->vertex(1)->point();
    const Point &p2 = f->vertex(2)->point();

    // For some reason, CGAL makes flat triagnles sometimes.
    const K::FT area = CGAL::area(p0, p1, p2); // non-negative area in K::FT
    if (area < K::FT(1e-5)) {
      return false; // near-flat triangle
    }

    // Half-open base window test via exact comparisons on K::FT
    const Point c = CGAL::centroid(p0, p1, p2);
    const double cx = CGAL::to_double(c.x());
    const double cy = CGAL::to_double(c.y());

    // We keep an extra row and column if using PBC. These are the elements
    // that wrap around.
    int extra = usingPBC ? 1 : 0;
    // Keep iff 0 <= cx < Lx and 0 <= cy < Ly (half-open [0,L))
    if (!(0 <= cx && cx < (cols - 1 + extra) * a && //
          0 <= cy && cy < (rows - 1 + extra) * a)) {
      return false;
    }

    seenTriangles.insert(key);
    return true;
  };

  // Collect faces to keep
  std::vector<DelaunayInfo::Face_handle> kept_faces;
  kept_faces.reserve(
      std::distance(dt.finite_faces_begin(), dt.finite_faces_end()));
  for (auto f = dt.finite_faces_begin(); f != dt.finite_faces_end(); ++f) {
    if (keep_element(f)) {
      kept_faces.push_back(f);
    }
  }

  // 4) Sanity: number of faces should match our preallocated elements
  const int nFaces = static_cast<int>(kept_faces.size());
  if (nFaces != nrElements) {
    elements.resize(nFaces);
    F_P_history_list.resize(nFaces);
    F_P_H.resize(nFaces, Matrix2d::Identity());
    std::cerr
        << "Warning: reconnectDelaunay(): triangle count mismatch (expected "
        << nrElements << ", got " << nFaces << ")." << std::endl;
    nrElements = nFaces;
  }

  // 5) Reset per-node connectivity
  for (long i = 0; i < nodes.size(); ++i) {
    nodes(i).elementCount = 0;
  }

  // 6) Refill elements from CGAL faces (each face -> one TElement)
  int eIdx = 0;
  for (const auto &f : kept_faces) {
    std::array<GhostNode, 3> g;
    for (int k = 0; k < 3; ++k) {
      auto v = f->vertex(k);
      const VInfo &vi = v->info();
      const Node *n = &nodes(vi.refNodeIndex); // reference node in MTS mesh
      // Build ghost from DT vertex coordinate. Use the targetPos to set the
      // periodic shift.
      g[k] = m_gn(n, toEigen(v->point()));
    }

    double noise = (eIdx < (int)elements.size()) ? elements[eIdx].noise
                                                 : sampleNormal(1, QDSD);
    elements[eIdx] = TElement((*this), g[0], g[1], g[2], eIdx, noise,
                              energyFunction, bulkModulus);
    ++eIdx;
  }

  markDirty();
}
#else
void Mesh::reconnectDelaunay() {
  throw std::runtime_error(
      "reconnectDelaunay requires CGAL; build without IDE_LIGHTWEIGHT.");
}
#endif

std::array<GhostNode, 4> Mesh::getElementPairNodes(const TElement &e1,
                                                   const TElement &e2) {
  auto e1Co = e1.getCoAngleNodes();
  auto e2Co = e2.getCoAngleNodes();

  const GhostNode *el1AngleNode = e1.getAngleNode();
  const GhostNode *el2AngleNode = e2.getAngleNode();

  return {*el1AngleNode, *e1Co[0], *e2Co[1], *el2AngleNode};
}

std::vector<GhostNode>
Mesh::getUniqueNodes(const std::vector<TElement *> elements) {
  auto compById = [](const GhostNode &a, const GhostNode &b) {
    if (a.id.x() != b.id.x())
      return a.id.x() < b.id.x();
    return a.id.y() < b.id.y();
  };
  std::set<GhostNode, decltype(compById)> uniqueNodes(compById);

  for (const auto &e : elements) {
    for (const auto &gn : e->ghostNodes) {
      uniqueNodes.insert(gn);
    }
  }

  return {uniqueNodes.begin(), uniqueNodes.end()};
}

static bool nodeReferencesElement(const Node &node, int elementIndex) {
  for (int i = 0; i < node.elementCount; ++i) {
    if (node.connectedElements[i] == elementIndex) {
      return true;
    }
  }
  return false;
}

static void assertNodesDoNotReferenceElements(const std::vector<Node *> &nodes,
                                              int e1i, int e2i,
                                              const std::string &context) {
  for (const Node *node : nodes) {
    if (node == nullptr) {
      throw std::runtime_error(context + ": null node pointer.");
    }
    if (nodeReferencesElement(*node, e1i) ||
        nodeReferencesElement(*node, e2i)) {
      throw std::runtime_error(context +
                               ": removed element still present in node "
                               "connectivity.");
    }
  }
}

static void assertDistinctTriangleGhosts(const std::array<GhostNode, 3> &ghosts,
                                         const std::string &context) {
  for (int i = 0; i < 3; ++i) {
    for (int j = i + 1; j < 3; ++j) {
      if (ghosts[i].id == ghosts[j].id) {
        throw std::runtime_error(context +
                                 ": triangle contains duplicate ghost-node "
                                 "ids.");
      }
      if (ghosts[i].referenceId == ghosts[j].referenceId) {
        throw std::runtime_error(
            context +
            ": triangle contains duplicate reference-node ids.");
      }
    }
  }
}

static void assertEdgeFlipGhostPairConsistency(
    const std::array<GhostNode, 3> &e1Ghosts,
    const std::array<GhostNode, 3> &e2Ghosts, const std::string &context,
    double tol = 1e-10) {
  assertDistinctTriangleGhosts(e1Ghosts, context + " element 1");
  assertDistinctTriangleGhosts(e2Ghosts, context + " element 2");

  std::vector<std::pair<const GhostNode *, const GhostNode *>> sharedGhosts;
  sharedGhosts.reserve(2);
  for (const auto &gn1 : e1Ghosts) {
    for (const auto &gn2 : e2Ghosts) {
      if (gn1.id == gn2.id) {
        sharedGhosts.push_back({&gn1, &gn2});
      }
    }
  }

  if (sharedGhosts.size() != 2) {
    throw std::runtime_error(context +
                             ": edge-flip pair does not share exactly two "
                             "ghost nodes.");
  }

  for (const auto &shared : sharedGhosts) {
    if (shared.first->referenceId != shared.second->referenceId) {
      throw std::runtime_error(
          context + ": shared ghost nodes disagree on reference id.");
    }
    if ((shared.first->pos - shared.second->pos).norm() > tol) {
      throw std::runtime_error(
          context + ": shared ghost nodes disagree on current position.");
    }
  }
}

void Mesh::removeElementsFromNodes(
    const std::array<const GhostNode *, 4> &gNodes, int e1i, int e2i) {
  const std::array<Node *, 4> nodes = {
      (*this)[gNodes[0]->referenceId], (*this)[gNodes[1]->referenceId],
      (*this)[gNodes[2]->referenceId], (*this)[gNodes[3]->referenceId]};
  removeElementsFromNodes(nodes, e1i, e2i);
}

void Mesh::removeElementsFromNodes(const std::array<Node *, 4> &nodes, int e1i,
                                   int e2i) {
  for (Node *n : nodes) {
    int tempElementIndices[MAX_ELEMENTS_PER_NODE];
    int tempNodeIndexInElement[MAX_ELEMENTS_PER_NODE];
    int newCount = 0;

    for (int i = 0; i < n->elementCount; i++) {
      const int el = n->connectedElements[i];
      if (el == e1i || el == e2i) {
        continue;
      }
      tempElementIndices[newCount] = n->connectedElements[i];
      tempNodeIndexInElement[newCount] = n->nodeIndexInElement[i];
      newCount++;
    }

    for (int i = 0; i < newCount; i++) {
      n->connectedElements[i] = tempElementIndices[i];
      n->nodeIndexInElement[i] = tempNodeIndexInElement[i];
    }
    for (int i = newCount; i < n->elementCount; i++) {
      n->connectedElements[i] = -1;
      n->nodeIndexInElement[i] = -1;
    }
    n->elementCount = newCount;
  }
}

void Mesh::flipEdge(TElement &e1, TElement &e2) {
  // This function takes two elements that should both have large angles, and
  // reconfigures the 4 nodes into two new elements that have smaller angles.

  const int e1i = e1.eIndex;
  const int e2i = e2.eIndex;
  if (e1i < 0 || e2i < 0 || e1i >= nrElements || e2i >= nrElements) {
    throw std::runtime_error("Mesh::flipEdge: invalid element indices.");
  }
  if (F_P_H.size() != elements.size()) {
    throw std::runtime_error(
        "Mesh::flipEdge: branch history size does not match elements.");
  }
  if (lastFlipDebugStates.size() != elements.size()) {
    lastFlipDebugStates.clear();
    lastFlipDebugStates.resize(elements.size());
  }

  const std::array<GhostNode, 3> oldGhosts1 = e1.ghostNodes;
  const std::array<GhostNode, 3> oldGhosts2 = e2.ghostNodes;
  const auto oe1 = e1.getAngleCo1Co2Nodes();
  const auto oe2 = e2.getAngleCo1Co2Nodes();

  // For now we follow the PDF with donor inheritance on option a only:
  // new element 1 inherits from old e1, new element 2 inherits from old e2.
  std::array<GhostNode, 3> n1a = {oe1[0], oe1[1], oe2[0]};
  n1a[2].updateReferencePosition(oe1[2].ref_pos);
  std::array<GhostNode, 3> n2a = {oe2[0], oe1[0], oe2[2]};
  n2a[1].updateReferencePosition(oe2[1].ref_pos);
  assertEdgeFlipGhostPairConsistency(n1a, n2a,
                                     "Mesh::flipEdge raw option-A candidates");

  // Note that evaluateEdgeFlipRemeshState currently may rotate the reference
  // state of the element. This needs to be remembered if/when we try to
  // evaluate option A vs B.
  const TElement::EdgeFlipRemeshState state1 =
      e1.evaluateEdgeFlipRemeshState(n1a, F_P_H[static_cast<size_t>(e1i)]);
  const TElement::EdgeFlipRemeshState state2 =
      e2.evaluateEdgeFlipRemeshState(n2a, F_P_H[static_cast<size_t>(e2i)]);
  if (!state1.valid || !state2.valid) {
    throw std::runtime_error(formatEdgeFlipCandidateFailure(
        *this, e1, e2, oldGhosts1, oldGhosts2, n1a, n2a,
        F_P_H[static_cast<size_t>(e1i)], F_P_H[static_cast<size_t>(e2i)],
        "evaluating fixed option-A remesh state",
        "encountered invalid edge-flip state"));
  }
  assertEdgeFlipGhostPairConsistency(
      state1.newGhostNodes, state2.newGhostNodes,
      "Mesh::flipEdge evaluated option-A states");

  auto assertDonorReferencePreserved = [&](const std::array<GhostNode, 3>
                                               &oldGhosts,
                                           const std::array<GhostNode, 3>
                                               &newGhosts,
                                           const std::string &context) {
    for (int i = 0; i < 3; ++i) {
      if (!oldGhosts[i].ref_pos.isApprox(newGhosts[i].ref_pos,
                                         kMatrixCompareTol)) {
        throw std::runtime_error(
            "Mesh::flipEdge: option-A donor reference element changed for " +
            context + ".");
      }
    }
  };

  // assertDonorReferencePreserved(oldGhosts1, state1.newGhostNodes, "element
  // 1"); assertDonorReferencePreserved(oldGhosts2, state2.newGhostNodes,
  // "element 2");

  removeElementFromNodes(elements[e1i]);
  removeElementFromNodes(elements[e2i]);

  createElementPair(state1.newGhostNodes, state2.newGhostNodes, e1i, e2i, true);
  F_P_H[static_cast<size_t>(e1i)] = state1.H_new;
  F_P_H[static_cast<size_t>(e2i)] = state2.H_new;
  elements[e1i].updateForces(*this);
  elements[e1i].updateGeometry();
  elements[e1i].updateFull();
  elements[e2i].updateForces(*this);
  elements[e2i].updateGeometry();
  elements[e2i].updateFull();

  auto assertAppliedStateMatches =
      [&](int elementIndex, const TElement::EdgeFlipRemeshState &state,
          const std::string &context) {
        const TElement &element = elements[static_cast<size_t>(elementIndex)];
        for (int i = 0; i < 3; ++i) {
          if (element.ghostNodes[i].referenceId !=
              state.newGhostNodes[i].referenceId) {
            throw std::runtime_error(
                "Mesh::flipEdge: rebuilt element changed node ordering for " +
                context + ".");
          }
          if (!element.ghostNodes[i].ref_pos.isApprox(
                  state.newGhostNodes[i].ref_pos, kMatrixCompareTol)) {
            throw std::runtime_error(
                "Mesh::flipEdge: rebuilt element changed reference nodes for " +
                context + ".");
          }
        }
        if (!element.F_P.isApprox(state.P_new, kMatrixCompareTol)) {
          throw std::runtime_error(
              "Mesh::flipEdge: rebuilt element F_P does not match evaluated "
              "remesh state for " +
              context + ".");
        }
      };

  assertAppliedStateMatches(e1i, state1, "element 1");
  assertAppliedStateMatches(e2i, state2, "element 2");
  assertEdgeFlipGhostPairConsistency(elements[e1i].ghostNodes,
                                     elements[e2i].ghostNodes,
                                     "Mesh::flipEdge rebuilt pair");

  F_P_history_list[e1i].push_back(state1.P_new);
  F_P_history_list[e2i].push_back(state2.P_new);

  LastFlipDebugState debug1;
  debug1.valid = true;
  debug1.partner = e2i;
  debug1.minIterationsAtFlip = nrMinItterations;
  debug1.minFunctionCallsAtFlip = nrMinFunctionCalls;
  debug1.applied_F_P = state1.P_new;
  debug1.oldSelfGhostNodes = oldGhosts1;
  debug1.oldPartnerGhostNodes = oldGhosts2;

  LastFlipDebugState debug2;
  debug2.valid = true;
  debug2.partner = e1i;
  debug2.minIterationsAtFlip = nrMinItterations;
  debug2.minFunctionCallsAtFlip = nrMinFunctionCalls;
  debug2.applied_F_P = state2.P_new;
  debug2.oldSelfGhostNodes = oldGhosts2;
  debug2.oldPartnerGhostNodes = oldGhosts1;

  lastFlipDebugStates[e1i] = debug1;
  lastFlipDebugStates[e2i] = debug2;

  totalEdgeFlipsInStep++;
  throwIfReductionExploded(elements[e1i], "Mesh::fixElementPair");
  throwIfReductionExploded(elements[e2i], "Mesh::fixElementPair");
}

bool Mesh::checkAndFixPeriodicElementPair(TElement &e1, TElement &e2) {
  // Check if the coAngleNodes are actually the same. If they are in different
  // periodic images, we need to move them together.
  auto e1Co = e1.getCoAngleNodes();
  auto e2Co = e2.getCoAngleNodes();
  if (e1Co[0]->id == e2Co[0]->id && e1Co[1]->id == e2Co[1]->id) {
    // The element are already together, no need to move elements
    return false;
  }
  // The situation here is like this: We would usually fix an element pair
  // that looks like this:
  // c2____a3
  //  |  \  |
  // a0____c1
  // But now, this element pair is accros the periodic boundary, like this:
  //  c2____a3           c2
  //      \  |    ...    |  \
  //        c1          a0____c1
  // What we need to do now, is to move one of the elements to the other
  // element so they can be reconnected properly. We will do this by considering
  // which element is closest to the center of the system. We will then move
  // the other one.

  // We move the one that is furthest away from the center of the mesh

  double distE1 = (e1.getCom() - com).norm();
  double distE2 = (e2.getCom() - com).norm();
  if (distE1 > distE2) {
    moveElementToTwin(e1, e2);
  } else {
    moveElementToTwin(e2, e1);
  }
  return true;
}

void Mesh::moveElementToTwin(TElement &elementToMove,
                             const TElement &fixedElement) {
  auto fixedCoAngleNodes = fixedElement.getCoAngleNodes();
  auto movingCoAngleNodes = elementToMove.getCoAngleNodes();

  const Vector2i deltaShift = fixedCoAngleNodes[0]->periodicShift -
                              movingCoAngleNodes[0]->periodicShift;
  const Vector2i secondDeltaShift = fixedCoAngleNodes[1]->periodicShift -
                                    movingCoAngleNodes[1]->periodicShift;
  if ((secondDeltaShift.array() != deltaShift.array()).any()) {
    throw std::runtime_error(
        "moveElementToTwin: shared edge nodes disagree on periodic shift.");
  }

  std::array<GhostNode, 3> shiftedGhostNodes = elementToMove.ghostNodes;
  for (GhostNode &gn : shiftedGhostNodes) {
    gn.applyPeriodicShift(deltaShift, latticeBasis, currentDeformation);
  }

  // We remove the element from all the nodes it is connected to.
  const std::vector<Node *> oldNodes = {
      (*this)[elementToMove.ghostNodes[0].referenceId],
      (*this)[elementToMove.ghostNodes[1].referenceId],
      (*this)[elementToMove.ghostNodes[2].referenceId]};
  removeElementFromNodes(elementToMove);
  assertNodesDoNotReferenceElements(oldNodes, elementToMove.eIndex,
                                    elementToMove.eIndex, "moveElementToTwin");

  // overwrite the element
  elementToMove = TElement{(*this),
                           shiftedGhostNodes[0],
                           shiftedGhostNodes[1],
                           shiftedGhostNodes[2],
                           elementToMove.eIndex,
                           elementToMove.noise,
                           energyFunction,
                           bulkModulus};
}

int Mesh::countConnectionsInGhostNode(const GhostNode &gn) {
  int connections = 0;

  // First we find the reference node
  Node *n = (*this)[gn.referenceId];

  // Then we need to go through each element it is connected to, and compare
  // with the ghost node in those elements to see if they are the same ghost
  // node as the one we are considering. We check if they are the same by
  // comparing their row, col AND periodicShift.

  for (int i = 0; i < n->elementCount; i++) {
    // Get one of the elements connected to the reference node.
    TElement &e = elements[n->connectedElements[i]];
    // Find the ghost node in that element that represents the reference node.
    GhostNode &gnInElement = e.ghostNodes[n->nodeIndexInElement[i]];
    // Check if this node is the same node as our input ghost node.
    // If it is, it means that our input ghost node is connected to this
    // element. (That means that the minimum number of connections will almost
    // always be 1, since we count the element that the input ghost node comes
    // from.)
    if (gnInElement.id == gn.id) {
      connections += 1;
    }
  }

  return connections;
}

void Mesh::setDiagonal(int row, int col, bool majorDiagonalOrder) {
  // get the 4 ghost nodes of the selected section of the mesh
  const std::array<GhostNode, 4> ghosts = getSquareGhostNodes(row, col);
  // Get the indexes of the elements
  auto [e1i, e2i] = getElementIndices(row, col);

  // Now we need to remove e1i and e2i from the 4 nodes, since they will be
  // added back in the tElement constructor
  removeElementsFromNodes(row, col, {e1i, e2i});

  // Update elements, preserving existing noise values
  createElementPair(ghosts, e1i, e2i, majorDiagonalOrder, true);
}

void Mesh::removeElementFromNodes(const TElement &element) {

  std::vector<const GhostNode *> ghostNodesVector;
  ghostNodesVector.reserve(3);

  for (const auto &gn : element.ghostNodes) {
    ghostNodesVector.push_back(&gn);
  }

  removeElementsFromNodes(ghostNodesVector, {element.eIndex});
}

void Mesh::removeElementsFromNodes(int row, int col,
                                   const std::vector<int> elIndexToRemove) {
  std::array<Node *, 4> nodes = getSquareNodes(row, col);
  removeElementsFromNodes({nodes[0], nodes[1], nodes[2], nodes[3]},
                          elIndexToRemove);
}

void Mesh::removeElementsFromNodes(const std::vector<const GhostNode *> gNodes,
                                   const std::vector<int> elIndexToRemove) {
  std::vector<Node *> nodes;
  nodes.reserve(gNodes.size()); // Reserve space for efficiency

  std::transform(
      gNodes.begin(), gNodes.end(), std::back_inserter(nodes),
      [this](const GhostNode *g) { return (*this)[g->referenceId]; });

  removeElementsFromNodes(nodes, elIndexToRemove);
}

void Mesh::removeElementsFromNodes(std::vector<Node *> nodes,
                                   const std::vector<int> elIndexToRemove) {

  // We check each node
  for (Node *n : nodes) {
    std::vector<int> indexesToRemove;

    // Find indices of elements to remove
    for (int i = 0; i < n->elementCount; i++) {
      for (int j = 0; j < elIndexToRemove.size(); j++) {
        if (n->connectedElements[i] == elIndexToRemove[j]) {
          indexesToRemove.push_back(i);
          break; // Found a match, no need to check other elements
        }
      }
    }

    // If there's nothing to remove, continue to next node
    if (indexesToRemove.empty()) {
      continue;
    }

    // Create temporary arrays to store elements we want to keep
    int tempElementIndices[MAX_ELEMENTS_PER_NODE];
    int tempNodeIndexInElement[MAX_ELEMENTS_PER_NODE];
    int newCount = 0;

    // Copy only the elements we want to keep
    for (int i = 0; i < n->elementCount; i++) {
      // Check if this index should be removed
      bool shouldRemove = false;
      for (int removeIdx : indexesToRemove) {
        if (i == removeIdx) {
          shouldRemove = true;
          break;
        }
      }

      // If not marked for removal, keep it
      if (!shouldRemove) {
        tempElementIndices[newCount] = n->connectedElements[i];
        tempNodeIndexInElement[newCount] = n->nodeIndexInElement[i];
        newCount++;
      }
    }

    // Update the node with the new arrays
    for (int i = 0; i < newCount; i++) {
      n->connectedElements[i] = tempElementIndices[i];
      n->nodeIndexInElement[i] = tempNodeIndexInElement[i];
    }
    for (int i = newCount; i < n->elementCount; i++) {
      n->connectedElements[i] = -1;
      n->nodeIndexInElement[i] = -1;
    }

    // Update the count
    n->elementCount = newCount;
  }
}

// Helper function to update positions using a generic buffer and its size
void Mesh::updateNodePositions(const double *data, size_t length) {
  const size_t nr_x_values = length / 2;
  const double *xData = data;
  const double *yData = data + nr_x_values;
  const size_t freeCount = freeNodeIds.size();
  for (size_t i = 0; i < freeCount; i++) {
    Node *n = (*this)[freeNodeIds[i]];
    n->setDisplacement({xData[i], yData[i]});
  }
  markDirty();
}

void Mesh::saveNodeDisplacements(std::vector<double> &displacements) const {
  const size_t nodeCount = static_cast<size_t>(rows * cols);
  displacements.resize(2 * nodeCount);
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      const Node &n = nodes(row, col);
      const Vector2d &disp = n.u();
      const size_t idx = static_cast<size_t>(row * cols + col);
      displacements[idx] = disp[0];
      displacements[idx + nodeCount] = disp[1];
    }
  }
}

void Mesh::captureDisplacementSnapshot(DisplacementSnapshot &snapshot) const {
  saveNodeDisplacements(snapshot.displacements);
}

Mesh::DisplacementSnapshot Mesh::displacementSnapshot() const {
  DisplacementSnapshot snapshot;
  captureDisplacementSnapshot(snapshot);
  return snapshot;
}

double Mesh::rmsDistanceToMesh(const Mesh &other, bool subtractAvgU) const {
  const size_t freeCount = freeNodeIds.size();
  if (freeCount == 0) {
    return 0.0;
  }
  if (rows != other.rows || cols != other.cols ||
      freeNodeIds.size() != other.freeNodeIds.size()) {
    throw std::runtime_error(
        "rmsDistanceToMesh: mesh dimensions or free-node counts differ.");
  }

  Vector2d avgUCurrent = Vector2d::Zero();
  Vector2d avgUOther = Vector2d::Zero();
  if (subtractAvgU) {
    for (size_t i = 0; i < freeCount; ++i) {
      avgUCurrent += (*this)[freeNodeIds[i]]->u();
      avgUOther += other[other.freeNodeIds[i]]->u();
    }
    avgUCurrent /= static_cast<double>(freeCount);
    avgUOther /= static_cast<double>(freeCount);
  }

  double sum = 0.0;
  for (size_t i = 0; i < freeCount; ++i) {
    Vector2d u = (*this)[freeNodeIds[i]]->u();
    Vector2d otherU = other[other.freeNodeIds[i]]->u();
    if (subtractAvgU) {
      u -= avgUCurrent;
      otherU -= avgUOther;
    }
    const Vector2d diff = u - otherU;
    sum += diff.squaredNorm();
  }
  return std::sqrt(sum / static_cast<double>(freeCount));
}

double
Mesh::rmsDistanceToDisplacementSnapshot(const DisplacementSnapshot &snapshot,
                                        bool subtractAvgU) const {
  const size_t freeCount = freeNodeIds.size();
  if (freeCount == 0) {
    return 0.0;
  }
  const size_t nodeCount = static_cast<size_t>(rows * cols);
  if (snapshot.displacements.size() != 2 * nodeCount) {
    throw std::runtime_error(
        "Displacement snapshot size does not match node count.");
  }

  const double *xSnapU = snapshot.displacements.data();
  const double *ySnapU = snapshot.displacements.data() + nodeCount;
  Vector2d avgUCurrent = Vector2d::Zero();
  Vector2d avgUSnap = Vector2d::Zero();
  if (subtractAvgU) {
    for (size_t i = 0; i < freeCount; i++) {
      const Node *n = (*this)[freeNodeIds[i]];
      const Vector2d u = n->u();
      const size_t idx = static_cast<size_t>(freeNodeIds[i].i);
      const Vector2d snapU = Vector2d{xSnapU[idx], ySnapU[idx]};
      avgUCurrent += u;
      avgUSnap += snapU;
    }
    avgUCurrent /= static_cast<double>(freeCount);
    avgUSnap /= static_cast<double>(freeCount);
  }

  double sum = 0.0;
  for (size_t i = 0; i < freeCount; i++) {
    const Node *n = (*this)[freeNodeIds[i]];
    Vector2d u = n->u();
    const size_t idx = static_cast<size_t>(freeNodeIds[i].i);
    Vector2d snapU = Vector2d{xSnapU[idx], ySnapU[idx]};
    if (subtractAvgU) {
      u -= avgUCurrent;
      snapU -= avgUSnap;
    }
    const Vector2d diff = u - snapU;
    sum += diff.squaredNorm();
  }
  return std::sqrt(sum / static_cast<double>(freeCount));
}

void Mesh::updateBoundingBox() {
  ensureForces();
  // Reset the bounding box
  bounds[0] = -100000; // max x
  bounds[1] = 100000;  // min x
  bounds[2] = -100000; // max y
  bounds[3] = 100000;  // min y

  // Update bounding box
  for (int i = 0; i < nrElements; i++) {
    for (int j = 0; j < 3; j++) {
      // Update bounds for x-coordinate
      if (elements[i].ghostNodes[j].pos[0] > bounds[0])
        bounds[0] = elements[i].ghostNodes[j].pos[0]; // max x
      if (elements[i].ghostNodes[j].pos[0] < bounds[1])
        bounds[1] = elements[i].ghostNodes[j].pos[0]; // min x

      // Update bounds for y-coordinate
      if (elements[i].ghostNodes[j].pos[1] > bounds[2])
        bounds[2] = elements[i].ghostNodes[j].pos[1]; // max y
      if (elements[i].ghostNodes[j].pos[1] < bounds[3])
        bounds[3] = elements[i].ghostNodes[j].pos[1]; // min y
    }
  }
}

void Mesh::updateAngles() {
  // It may be strange that updateAngles is not part of updateElements, but
  // angles are not used in the simulation, and relatively expensive to
  // calculate. They could be useful for analysis, so we do save them in the vtu
  // files.
  ensureForces();
  for (auto &el : elements) {
    el.updateAngles();
  }
}

void Mesh::updateAveragesAndPlasticEvents() {
  ensureFull();

  // Note that totalEnergy has already been calculated since we use it in
  // the energy minimization

  // we reset the maxEnergy
  maxEnergy = 0;

  // We calculate total force for debugging (should be zero)
  // Vector2d totalForce = Vector2d::Zero();
  // for (int i = 0; i < nodes.size(); i++) {
  //   totalForce += nodes(i).f;
  // }
  // // Mathematically, the total force should always be 0 (even out of
  // // equilibrium), but due to some rounding errors (i think), we need to
  // allow
  // // some freedom before we declare that something is wrong.
  // if (totalForce.norm() / nrNodes > 1e-13) {
  //   std::cout << "Min step: " << nrMinItterations << '\n';
  //   std::cout << "Total force: " << totalForce << '\n';
  //   std::cout << "Big force: " << totalForce.norm() << '\n';
  //   if (endOfStep) {

  //     std::cerr << "Total force is not zero. Something is wrong.";
  //   }
  // }

  // This is the total energy from all the triangles
  nr_elements_with_m3_change = 0;
  nr_elements_with_m3_changeInStep = 0;
  redQuadrantCounts = {0, 0, 0, 0};
  redQuadrantFixedCounts = {0, 0, 0, 0};
  double totalP11 = 0;
  double totalP12 = 0;
  double totalP21 = 0;
  double totalP22 = 0;
  double totalSigma11 = 0;
  double totalSigma12 = 0;
  double totalSigma22 = 0;
  double totalSigmaTrace = 0;
  double totalReferenceRotationTheta = 0;
  sumM3Nr = 0;
  for (int i = 0; i < nrElements; i++) {
    const TElement &e = elements[i];
    totalP11 += e.P(0, 0);
    totalP12 += e.P(0, 1);
    totalP21 += e.P(1, 0);
    totalP22 += e.P(1, 1);
    totalSigma11 += e.sigma(0, 0);
    totalSigma12 += e.sigma(0, 1);
    // std::cout << e.sigma(0, 1) << "\n";
    totalSigma22 += e.sigma(1, 1);
    totalSigmaTrace += e.sigma.trace();
    totalReferenceRotationTheta += e.referenceRotationTheta();
    if (e.red_quadrant >= 1 && e.red_quadrant <= 4) {
      redQuadrantCounts[static_cast<size_t>(e.red_quadrant - 1)] += 1;
    }

    // We also keep track of the highest energy and some other things
    if (e.energy > maxEnergy) {
      maxEnergy = e.energy;
    }

    // Plastic event counters are based on pastM3Nr_fix/pastStepM3Nr_fix, which
    // are reset in resetPastPlasticCount.
    if (e.m3Nr > maxM3Nr) {
      maxM3Nr = e.m3Nr;
    }
    sumM3Nr += e.m3Nr;
    int plasticChange = e.m3Nr - e.pastM3Nr;
    if (plasticChange > maxPlasticJump) {
      maxPlasticJump = plasticChange;
    } else if (plasticChange < minPlasticJump) {
      minPlasticJump = plasticChange;
    }
    if (e.pastM3Nr != e.m3Nr) {
      nr_elements_with_m3_change += 1;
    }
    if (e.pastStepM3Nr != e.m3Nr) {
      nr_elements_with_m3_changeInStep += 1;
    }
  }

  averageEnergy = totalEnergy / nrElements;
  averageP11 = totalP11 / nrElements;
  averageP12 = totalP12 / nrElements;
  averageP21 = totalP21 / nrElements;
  averageP22 = totalP22 / nrElements;
  averageSigma11 = totalSigma11 / nrElements;
  averageSigma12 = totalSigma12 / nrElements;
  averageSigma22 = totalSigma22 / nrElements;
  averageSigmaTrace = totalSigmaTrace / nrElements;
  averageReferenceRotationTheta = totalReferenceRotationTheta / nrElements;
}

void Mesh::updateCom() {
  ensureForces();
  com = Vector2d::Zero();
  for (int i = 0; i < nrElements; i++) {
    com += elements[i].getCom();
  }
  com /= nrElements;
}

// Moves a section of the mesh based on spatial coordinates (x, y).
void Mesh::moveMeshSection(double minX, double minY, Vector2d disp,
                           bool moveFixed, bool moveFree, double maxX,
                           double maxY) {
  auto isInBounds = [&](const Node &n) {
    const Vector2d &p = n.pos();
    return (p[0] >= minX && p[0] <= maxX && p[1] >= minY && p[1] <= maxY);
  };

  if (moveFixed) {
    for (const NodeId &nId : fixedNodeIds) {
      Node *n = (*this)[nId];
      if (isInBounds(*n)) {
        n->addDisplacement(disp);
      }
    }
  }

  if (moveFree) {
    for (const NodeId &nId : freeNodeIds) {
      Node *n = (*this)[nId];
      if (isInBounds(*n)) {
        n->addDisplacement(disp);
      }
    }
  }
  markDirty();
}

void Mesh::writeToVtu(std::string filename, bool minimizationStep,
                      bool useReferenceElements) {
  writeToVtu(std::move(filename), minimizationStep, VtuFieldLevel::All, "",
             useReferenceElements);
}

void Mesh::writeToVtu(std::string filename, bool minimizationStep,
                      VtuFieldLevel level, std::string nameSuffix,
                      bool useReferenceElements) {
  ensureFull();
  // We keep one exporter and optionally emit a paired reference mesh for
  // minimization-step logging.
  writeMeshToVtu((*this), simName, dataPath, filename, minimizationStep, level,
                 nameSuffix, "", useReferenceElements);

  if (minimizationStep && !useReferenceElements) {
    const std::string referenceSuffix =
        nameSuffix.empty() ? "reference" : nameSuffix + "_reference";
    writeMeshToVtu((*this), simName, dataPath, filename, minimizationStep,
                   level, referenceSuffix, "", true);
  }
}

void transform(const Matrix2d &matrix, Mesh &mesh,
               std::vector<NodeId> nodesToTransform) {
  // We get the adress of each node
  for (NodeId &nodeId : nodesToTransform) {
    transformInPlace(matrix, *mesh[nodeId]);
  }
  mesh.markDirty();
}

std::ostream &operator<<(std::ostream &os, const Mesh &mesh) {
  os << "Mesh: " << mesh.rows << " x " << mesh.cols << " nodes\nIds:\n";
  for (int i = mesh.cols - 1; i >= 0; --i) {
    for (int j = 0; j < mesh.rows; ++j) {
      Node n = mesh.nodes(i, j);
      os << n.id.i << "\t";
    }
    os << "\n";
  }
  os << "Positions:\n";
  for (int i = mesh.cols - 1; i >= 0; --i) {
    for (int j = 0; j < mesh.rows; ++j) {
      Node n = mesh.nodes(i, j);
      // set precision to 1
      os << std::fixed << std::setprecision(1) << n.pos().transpose() << "    ";
    }
    os << "\n";
  }

  return os;
}

void transform(const Matrix2d &matrix, Mesh &mesh) {
  // Transform all nodes
  transform(matrix, mesh, mesh.fixedNodeIds);
  transform(matrix, mesh, mesh.freeNodeIds);
}

void translate(Mesh &mesh, std::vector<NodeId> nodesToTranslate, double x,
               double y) {
  // We get the adress of each node
  for (NodeId &nodeId : nodesToTranslate) {
    translateInPlace(*mesh[nodeId], x, y);
  }
  mesh.markDirty();
}
void translate(Mesh &mesh, double x, double y) {
  translate(mesh, mesh.fixedNodeIds, x, y);
  translate(mesh, mesh.freeNodeIds, x, y);
}
