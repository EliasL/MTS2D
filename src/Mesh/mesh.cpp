#include "mesh.h"
#include "Data/data_export.h"
#include "Data/logging.h"
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
#include <omp.h>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
using Eigen::Vector2i;
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
    row = rows + row;
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

void Mesh::resetCounters() {
  nrMinItterations = 0;
  nrMinFunctionCalls = 0;
  nrMinItterationsSinceLastReconnect = 0;
  nrMinFunctionCallsSinceLastReconnect = 0;
  totalEdgeFlipsInStep = 0;
  edgeFlipChosenMinusOtherEnergyInStep = 0.0;
  edgeFlipAlwaysChoseLowerEnergyInStep = true;
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
                                    std::string_view context) const {
  if (element.m3Nr <= 200) {
    return;
  }
  throw std::runtime_error(
      DebugLog::formatReductionExplosion(element, *this, context));
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
  double maxForce = 0.0;
  int maxM3InForceUpdate = 0;

  // Note that floating-point addition is not associative:
  //   (a + b) + c != a + (b + c)
  // The order of an OpenMP reduction can therefore vary, so using multiple
  // threads can sometimes lead to irreproducible simulation results.
#pragma omp parallel reduction(+ : energy_sum)                                  \
    reduction(max : maxForce, maxM3InForceUpdate)
  {
    const int tid = omp_get_thread_num();
    Vector2d *local = forceScratch.data() + static_cast<size_t>(tid) * nNodes;

    // Parallel first-touch / zeroing of this thread's stripe
    for (int i = 0; i < nNodes; ++i) {
      local[i].setZero();
    }

#pragma omp for schedule(static)
    for (int i = 0; i < nrElements; ++i) {
      TElement &e = elements[i];
      e.updateForces(*this);
      energy_sum += e.energy;
      maxM3InForceUpdate = std::max(maxM3InForceUpdate, e.m3Nr);

      const GhostNode &g0 = e.ghostNodes[0];
      const GhostNode &g1 = e.ghostNodes[1];
      const GhostNode &g2 = e.ghostNodes[2];

      local[g0.referenceId.i] += g0.f;
      local[g1.referenceId.i] += g1.f;
      local[g2.referenceId.i] += g2.f;
    }

#pragma omp for schedule(static)
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
  }

  totalEnergy = energy_sum;
  this->maxForce = maxForce;

  // Sometimes the initial guess of the lgbfs algorithm is unlucky with it's
  // first (large) guess. That sometimes leads to highly deformed elements, but
  // should not be counted as a "reduction explosion".
  if (maxM3InForceUpdate > 200 && nrMinItterationsSinceLastReconnect>10) {
    for (int i = 0; i < nrElements; ++i) {
      throwIfReductionExploded(elements[i], "Mesh::updateElementsForces");
    }
  }
}

void Mesh::updateElementsGeometry() {
#pragma omp parallel for schedule(static, 1024)
  for (int i = 0; i < nrElements; ++i) {
    elements[i].updateGeometry();
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

// Updates the forces on the nodes in the surface and returns the total
// energy from all the elements in the surface.
void Mesh::updateMesh() { ensureForces(); }

void Mesh::updateElements() { ensureForces(); }

void Mesh::markDirty() { updateState = UpdateState::Dirty; }

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
    throw std::runtime_error(DebugLog::formatDisplacementSnapshotSizeError(
        "Mesh::rmsDistanceToDisplacementSnapshot",
        snapshot.displacements.size(), 2 * nodeCount, rows, cols));
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

namespace {

struct MeshAverageTotals {
  double p11 = 0.0;
  double p12 = 0.0;
  double p21 = 0.0;
  double p22 = 0.0;
  double sigma11 = 0.0;
  double sigma12 = 0.0;
  double sigma22 = 0.0;
  double sigmaTrace = 0.0;
};

void resetAverageFields(Mesh &mesh) {
  mesh.maxEnergy = 0;
  mesh.redQuadrantCounts = {0, 0, 0, 0};
  mesh.redQuadrantFixedCounts = {0, 0, 0, 0};
  mesh.sumM3Nr = 0;
}

void accumulateElementAverages(Mesh &mesh, MeshAverageTotals &totals,
                               const TElement &e, const Matrix2d &sigma) {
  totals.p11 += e.P(0, 0);
  totals.p12 += e.P(0, 1);
  totals.p21 += e.P(1, 0);
  totals.p22 += e.P(1, 1);
  totals.sigma11 += sigma(0, 0);
  totals.sigma12 += sigma(0, 1);
  totals.sigma22 += sigma(1, 1);
  totals.sigmaTrace += sigma.trace();

  if (e.red_quadrant >= 1 && e.red_quadrant <= 4) {
    mesh.redQuadrantCounts[static_cast<size_t>(e.red_quadrant - 1)] += 1;
  }
  mesh.maxEnergy = std::max(mesh.maxEnergy, e.energy);
  mesh.maxM3Nr = std::max(mesh.maxM3Nr, e.m3Nr);
  mesh.sumM3Nr += e.m3Nr;
}

void publishElementAverages(Mesh &mesh, const MeshAverageTotals &totals) {
  if (mesh.nrElements <= 0 ||
      mesh.elements.size() != static_cast<size_t>(mesh.nrElements)) {
    throw std::runtime_error(
        "Mesh averages require a positive, consistent element count.");
  }

  const double n = static_cast<double>(mesh.nrElements);
  mesh.averageEnergy = mesh.totalEnergy / n;
  mesh.averageP11 = totals.p11 / n;
  mesh.averageP12 = totals.p12 / n;
  mesh.averageP21 = totals.p21 / n;
  mesh.averageP22 = totals.p22 / n;
  mesh.averageSigma11 = totals.sigma11 / n;
  mesh.averageSigma12 = totals.sigma12 / n;
  mesh.averageSigma22 = totals.sigma22 / n;
  mesh.averageSigmaTrace = totals.sigmaTrace / n;
}

void updateAverageFields(Mesh &mesh, bool recomputeSigmaFromForceState) {
  resetAverageFields(mesh);
  MeshAverageTotals totals;
  for (const TElement &e : mesh.elements) {
    Matrix2d sigma = e.sigma;
    if (recomputeSigmaFromForceState) {
      sigma = (1.0 / e.F.determinant()) * e.P * e.F.transpose();
    }
    accumulateElementAverages(mesh, totals, e, sigma);
  }
  publishElementAverages(mesh, totals);
}

void moveMeshSectionImpl(Mesh &mesh, double minX, double minY, Vector2d disp,
                         bool moveFixed, bool moveFree, bool hasMaxX,
                         double maxX, bool hasMaxY, double maxY) {
  auto isInBounds = [&](const Node &n) {
    const Vector2d &p = n.pos();
    if (p[0] < minX || p[1] < minY) {
      return false;
    }
    if (hasMaxX && p[0] > maxX) {
      return false;
    }
    if (hasMaxY && p[1] > maxY) {
      return false;
    }
    return true;
  };

  auto moveNodes = [&](const std::vector<NodeId> &nodeIds) {
    for (const NodeId &nId : nodeIds) {
      Node *n = mesh[nId];
      if (isInBounds(*n)) {
        n->addDisplacement(disp);
      }
    }
  };

  if (moveFixed) {
    moveNodes(mesh.fixedNodeIds);
  }
  if (moveFree) {
    moveNodes(mesh.freeNodeIds);
  }
  mesh.markDirty();
}

} // namespace

void Mesh::updateAveragesAndPlasticEvents() {
  ensureFull();
  updateAverageFields(*this, false);
  updatePlasticEventCounts();
}

void Mesh::updateForceStateAveragesAndPlasticEvents() {
  ensureForces();
  updateAverageFields(*this, true);
  updatePlasticEventCounts();
}

void Mesh::updatePlasticEventCounts() {
  int m3Change = 0;
  int m3ChangeInStep = 0;
  for (const TElement &e : elements) {
    const int plasticChange = e.m3Nr - e.pastM3Nr;
    if (plasticChange > maxPlasticJump) {
      maxPlasticJump = plasticChange;
    } else if (plasticChange < minPlasticJump) {
      minPlasticJump = plasticChange;
    }
    m3Change += (e.pastM3Nr != e.m3Nr);
    m3ChangeInStep += (e.pastStepM3Nr != e.m3Nr);
  }
  nr_elements_with_m3_change = m3Change;
  nr_elements_with_m3_changeInStep = m3ChangeInStep;
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
                           bool moveFixed, bool moveFree) {
  moveMeshSectionImpl(*this, minX, minY, disp, moveFixed, moveFree, false, 0.0,
                      false, 0.0);
}

void Mesh::moveMeshSection(double minX, double minY, Vector2d disp,
                           bool moveFixed, bool moveFree, double maxX) {
  moveMeshSectionImpl(*this, minX, minY, disp, moveFixed, moveFree, true, maxX,
                      false, 0.0);
}

void Mesh::moveMeshSection(double minX, double minY, Vector2d disp,
                           bool moveFixed, bool moveFree, double maxX,
                           double maxY) {
  moveMeshSectionImpl(*this, minX, minY, disp, moveFixed, moveFree, true, maxX,
                      true, maxY);
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
