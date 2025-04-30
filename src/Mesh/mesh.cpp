#include "mesh.h"
#include "Data/data_export.h"
#include "Eigen/src/Core/Matrix.h"
#include "Mesh/node.h"
#include "Mesh/tElement.h"
#include "Simulation/randomUtils.h"
#include <array>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

Mesh::Mesh() {}

// Constructor that initializes the surface with size n x m
Mesh::Mesh(int rows, int cols, double a, double QDSD, bool usingPBC,
           std::string diagonal)
    : nodes(rows, cols), a(a), rows(rows), cols(cols), loadSteps(0),
      currentDeformation(Eigen::Matrix2d::Identity()), QDSD(QDSD),
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
  remeshedElements.resize(nrElements);

  m_createNodes();
  m_updateFixedAndFreeNodeIds();
  m_fillNeighbours();

  // we create some elements here, but note that these will need to be replaced
  // if changes are made to the fixed and free status of the nodes.
  createElements();
  calculateAverages();
  resetCounters();
}

Mesh::Mesh(int rows, int cols, bool usingPBC, std::string diagonal)
    : Mesh(rows, cols, 1, 0, usingPBC, diagonal) {}
Mesh::Mesh(int rows, int cols, bool usingPBC)
    : Mesh(rows, cols, usingPBC, "major") {}

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
}

void Mesh::applyTransformationToFixedNodes(const Matrix2d &transformation) {
  // We get the id of each node in the border
  for (NodeId &nodeId : fixedNodeIds) {
    transformInPlace(transformation, *(*this)[nodeId]);
  }
}

void Mesh::applyTransformationToSystemDeformation(
    const Matrix2d &transformation) {
  currentDeformation = transformation * currentDeformation;
}

void Mesh::applyTranslation(const Vector2d &displacement) {
  for (long i = 0; i < nodes.size(); i++) {
    translateInPlace(nodes(i), displacement);
  }
}

void Mesh::setInitPos() {
  for (long i = 0; i < nodes.size(); i++) {
    nodes(i).setInitPos(nodes(i).pos());
  }
}

Node *Mesh::m_getNeighbourNode(const Node &node, int direction) {
  return (*this)[(*this)[node.id]->refNeighbours[direction]];
}

// Function to fix the elements of the border vector
void Mesh::fixBorderNodes() {
  fixNodesInRow(0);
  fixNodesInColumn(0);
  fixNodesInRow(rows - 1);
  fixNodesInColumn(cols - 1);
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
      // Set the x and y positions based on the surface indices and spacing "a"
      nodes(row, col) = Node(a, row, col, m);
    }
  }
}

void Mesh::m_fillNeighbours() {
  int n = rows; // Number of rows
  int m = cols; // Number of columns

  for (int row = 0; row < n; ++row) // Iterate over rows
  {
    for (int col = 0; col < m; ++col) // Iterate over columns
    {
      // Define neighbor indices using periodic boundary conditions
      int left = (col == 0) ? m - 1 : col - 1;
      int right = (col == m - 1) ? 0 : col + 1;
      int up = (row == n - 1) ? 0 : row + 1;
      int down = (row == 0) ? n - 1 : row - 1;

      // Fill in the neighbors
      nodes(row, col).refNeighbours[LEFT_N] = m_makeNId(row, left);   // left
      nodes(row, col).refNeighbours[RIGHT_N] = m_makeNId(row, right); // right
      nodes(row, col).refNeighbours[UP_N] = m_makeNId(up, col);       // up
      nodes(row, col).refNeighbours[DOWN_N] = m_makeNId(down, col);   // down
    }
  }
}

// Gets the four nodes and their ghost versions for a given row and column in
// the reference state. NOT in the current state.
std::vector<Node *> Mesh::getSquareNodes(int row, int col) {
  // We find the 4 nodes in the current square
  Node *n1 = (*this)[m_makeNId(row, col)];
  Node *n2 = m_getNeighbourNode(*n1, RIGHT_N);
  Node *n3 = m_getNeighbourNode(*n1, UP_N);
  // n4 is now up AND right of n1
  Node *n4 = m_getNeighbourNode(*n2, UP_N);

  // The nodes should be in this square configuration in the reference state
  // n3  n4
  // n1  n2

  // This function will find nodes across the system if necessary due to
  // periodic boundaries
  return {n1, n2, n3, n4};
}

// Gets the four nodes and their ghost versions for a given row and column
std::vector<GhostNode> Mesh::getSquareGhostNodes(int row, int col) {
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
void Mesh::createElementPair(const std::vector<GhostNode> &ghosts, int e1i,
                             int e2i, bool useMajorDiagonal,
                             bool preserveNoise) {

  // The nodes should be in this square configuration in the current state
  // g3  g4
  // g1  g2
  assert(ghosts.size() == 4);

  GhostNode g1 = ghosts[0];
  GhostNode g2 = ghosts[1];
  GhostNode g3 = ghosts[2];
  GhostNode g4 = ghosts[3];

  double noise1, noise2;

  if (preserveNoise) {
    noise1 = elements[e1i].noise;
    noise2 = elements[e2i].noise;
  } else {
    noise1 = sampleNormal(1, QDSD);
    noise2 = sampleNormal(1, QDSD);
  }

  // When choosing what order to give the nodes to the element, we carefully
  // choose the first node to be the corner node, such that, in the reference
  // frame, all elements have an angle of 90 degrees.

  if (useMajorDiagonal) {
    // Split using major-diagonal from top-left to bottom-right (↘)
    elements[e1i] = TElement((*this), g1, g2, g3, e1i, noise1);
    elements[e2i] = TElement((*this), g4, g2, g3, e2i, noise2);
  } else {
    // Split using minor-diagonal from top-right to bottom-left (↙)
    elements[e1i] = TElement((*this), g2, g1, g4, e1i, noise1);
    elements[e2i] = TElement((*this), g3, g1, g4, e2i, noise2);
  }
  elements[e1i].update((*this));
  elements[e2i].update((*this));
}

void Mesh::createElementPair(const std::vector<const GhostNode *> &ghostsPtr,
                             int e1i, int e2i, bool useMajorDiagonal,
                             bool preserveNoise) {

  const std::vector<GhostNode> ghosts = {*ghostsPtr[0], *ghostsPtr[1],
                                         *ghostsPtr[2], *ghostsPtr[3]};
  createElementPair(ghosts, e1i, e2i, useMajorDiagonal, preserveNoise);
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

      bool useMajorDiagonal;
      if (diagonal == "major") {
        useMajorDiagonal = true;
      } else if (diagonal == "minor") {
        useMajorDiagonal = false;
      } else if (diagonal == "alternate") {
        useMajorDiagonal = (row + col) % 2;
      } else {
        throw std::invalid_argument("Unkown meshing: " + diagonal);
      }

      createElementPair(ghosts, e1i, e2i, useMajorDiagonal, false);
    }
  }
}

// This is just a function to avoid having to write cols
NodeId Mesh::m_makeNId(int row, int col) { return NodeId(row, col, cols); }

void Mesh::resetCounters() {
  nrMinItterations = 0;
  nrMinFunctionCalls = 0;
  resetPastPlasticCount();
  // Set all elements to not remeshed
  std::fill(remeshedElements.begin(), remeshedElements.end(), false);
}

void Mesh::resetPastPlasticCount(bool endOfStep) {
  nrPlasticChangesInStep = 0;
  for (size_t i = 0; i < elements.size(); i++) {
    elements[i].pastStepM3Nr = elements[i].m3Nr;
  }

  if (endOfStep) {
    nrPlasticChanges = 0;
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
GhostNode Mesh::m_gn(const Node *n, int row, int col) {
  return GhostNode(n, row, col, cols, a, currentDeformation);
}
GhostNode Mesh::m_gn(const Node *n) {
  return GhostNode(n, n->id.row(), n->id.col(), cols, a, currentDeformation);
}

// The idea here is that we have taken our grid, and made rows and columns of
// squares of 4 and 4 nodes. This function converts from reference nodes to
// ghost nodes and makes sure that the ghost nodes are appropreately shifted
// when that is requried.
// We never want our square to span the entire system which would happen at the
// boundary if we use the "real" position of the nodes
std::vector<GhostNode>
Mesh::m_makeGhostNodes(const std::vector<Node *> refNodes, int row, int col) {
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

// Updates the forces on the nodes in the surface and returns the total
// energy from all the elements in the surface.
void Mesh::updateMesh() {
  // Now we update all the elements using the current positions of the nodes
  updateElements();

  // We then add the force from the elements back to the nodes
  applyForceFromElementsToNodes();
}

void Mesh::updateElements() {
  // It seems like some elements are faster to update than others. This is
  // a bit strange. The only thing i can think of is the lagrange reduction,
  // but that one should be qutie fast?
#pragma omp parallel for schedule(guided)
  // dynamic, nrElements / (10 * omp_get_max_threads()))
  for (int i = 0; i < nrElements; i++) {
    elements[i].update(*this);
  }

  totalEnergy = 0;
  for (int i = 0; i < nrElements; i++) {
    totalEnergy += elements[i].energy;
  }
}

void Mesh::updateNodeForce(Node &node) {
  node.resetForce();
  for (size_t e = 0; e < node.elementCount; ++e) {
    const int elementNr = node.elementIndices[e];
    const int nodeNrInElement = node.nodeIndexInElement[e];
    assert(nodeNrInElement != -1); // This should never happen
    const TElement &element = elements[elementNr];
    const GhostNode &elementNode = element.ghostNodes[nodeNrInElement];
    node.f += elementNode.f;
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

// This function goes through each pair of elements and checks if the pair
// should flip their diagonal
bool Mesh::remesh(bool lockElements) {
  hasRemeshed = false;
  for (int i = 0; i < elements.size(); i++) {
    continue;
    // If the elements are locked and we have already remeshed this element,
    // we continue
    if (lockElements && remeshedElements[i]) {
      continue;
    }

    TElement &e = elements[i];
    // Checking if the largest angle is larger than 90 degrees
    double eps = 0.1;
    // if (e.C(0, 1) < -eps || e.C(0, 1) > 1 + eps) {
    if (e.largestAngle > 90 + eps) {
      double badness = (e.largestAngle - 90) / 5;
      int twinIndex = e.getElementTwin(*(this));

      // If we found a twin
      if (twinIndex != -1) {
        TElement &twin = elements[twinIndex];

        // We check that it is similar to our current element
        if (abs(e.largestAngle - twin.largestAngle) < eps / 2 + badness) {
          // writeMeshToVtu((*this), simName, dataPath,
          //                std::to_string(twin.eIndex) + "Before remesh");
          // std::cout << "Angles: " << e.largestAngle << ", " <<
          // twin.largestAngle
          //           << ", " << abs(e.largestAngle - twin.largestAngle) <<
          //           '\n';
          // std::cout << "Pre remesh" << '\n';
          // checkForces(getElementPairNodes(e, twin));
          fixElementPair(e, twin);
          remeshedElements[i] = true;
          remeshedElements[twinIndex] = true;
          hasRemeshed = true;
          // checkForces(getUniqueNodes({&e, &twin}));

          // writeMeshToVtu((*this), simName, dataPath,
          //                std::to_string(twin.eIndex) + "After OK remesh");
        }
      }
    }
  }
  return hasRemeshed;
}

std::vector<GhostNode> Mesh::getElementPairNodes(const TElement &e1,
                                                 const TElement &e2) {

  // Note that this function only works for elements where the angle nodes are
  // "seeing" each other.
  //  Extract the nodes with large angles from each element
  const GhostNode &el1AngleNode = e1.ghostNodes[e1.angleNode];
  const GhostNode &el2AngleNode = e2.ghostNodes[e2.angleNode];

  // Extract the other nodes that are common between the elements
  auto coAngleNodes = e1.getCoAngleNodes();

  // Arrange the nodes in the required order for createElementPair:
  // {angle0, coNode1, coNode2, angle3}:
  // c2____a3
  //  |  \  |
  // a0____c1
  // With angle nodes in positions 0 and 3
  const std::vector<GhostNode> assumedOrder = {el1AngleNode, *coAngleNodes[0],
                                               *coAngleNodes[1], el2AngleNode};

  // Confirm that we have the same coAngleNodes
  assert(coAngleNodes[0]->id == e2.getCoAngleNodes()[0]->id);
  assert(coAngleNodes[1]->id == e2.getCoAngleNodes()[1]->id);

  return assumedOrder;
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

void Mesh::fixElementPair(TElement &e1, TElement &e2) {
  // This function takes two elements that should both have large angles, and
  // reconfigures the 4 nodes into two new elements that have smaller angles.

  // PBC test
  // Due to the periodic horse problem (TODO provide reference), we need to
  // check if the coAngleNodes are the same. They should have the same
  // ghostNode index. If they don't, we need to move a real node to the other
  // side of the periodic boundary. (We only need to check one. )
  const GhostNode *e1gn = e1.getCoAngleNodes()[0];
  const GhostNode *e2gn = e2.getCoAngleNodes()[0];
  if (e1gn->id != e2gn->id) {

    fixPeriodicElementPair(e1, e2);
  }
  const std::vector<GhostNode> standardOrder = getElementPairNodes(e1, e2);

  // When we give these nodes to the createElementPair function, it is
  // important To consider the order in which we give them. The function will
  // interpret the list of nodes like this {angle0, coNode1, coNode2, angle3}:
  // c2____a3
  //  |  \  |
  // a0____c1
  // Assuming that we use major diagonal (we could use either so long as we
  // change the order we give the nodes in), we will make g2 and g3 be the new
  // corner nodes. That means we should make the angle nodes go in the first
  // and last possition. The order of the two others nodes does not matter,
  // but by convention, we want to make the node with the smaller index come
  // first

  // This new order will swich what nodes are angle nodes and coAngle nodes
  const std::vector<const GhostNode *> newPairOrder = {
      &standardOrder[1], &standardOrder[0], &standardOrder[3],
      &standardOrder[2]};

  // Now we disconnect the elements from the nodes, and create new elements
  removeElementsFromNodes(newPairOrder, {e1.eIndex, e2.eIndex});

  createElementPair(newPairOrder, e1.eIndex, e2.eIndex, true);
}

void Mesh::fixPeriodicElementPair(TElement &e1, TElement &e2) {
  // The situation here is like this: We would usually fix an element pair
  // that looks like this: c2____a3
  //  |  \  |
  // a0____c1
  // But now, this element pair is accros the periodic boundary, like this:
  //  c2____a3           c2
  //      \  |    ...    |  \
  //        c1          a0____c1
  // What we need to do now, is to move one of the elements to the other
  // element so they can be remeshed properly. We will do this by considering
  // which element is closest to the center of the system. We will then move
  // the other one.

  // We simply move the one that is furthest away from the center of the mesh

  double distE1 = (e1.getCom() - com).norm();
  double distE2 = (e2.getCom() - com).norm();
  if (distE1 > distE2) {
    moveElementToTwin(e1, e2);
  } else {
    moveElementToTwin(e2, e1);
  }
}

void Mesh::moveElementToTwin(TElement &elementToMove,
                             const TElement &fixedElement) {
  // In order to move the element, we just copy the coAngleNodes from the
  // fixed element and find the appropreate periodic shift for the final angle
  // node.

  // Get the coAngleNodes of the fixed element, as these are the ghost
  // nodes we should be using to create the new element.
  auto fixedCoAngleNodes = fixedElement.getCoAngleNodes();
  // This is the node that we should move by giving it a new periodic shift
  // (We are not actually moving the reference node, but only the ghost node).
  Node *refNode = (*this)[elementToMove.ghostNodes[elementToMove.angleNode]];

  // There are only 8 possible shifts (excluding 0,0) the last element can
  // have. We try all possible shift and choose the one that minimizes the
  // distance to the coAngleNodes
  Vector2i minPShift = {};
  double minDistance = std::numeric_limits<double>::infinity();

  for (int i = -1; i < 2; i++) {
    for (int j = -1; j < 2; j++) {
      Vector2i shift{i * cols, j * rows};
      GhostNode testNode(refNode, shift, cols, a, currentDeformation);
      Vector2d p1 = fixedCoAngleNodes[0]->pos;
      Vector2d p2 = fixedCoAngleNodes[1]->pos;
      double distance = (testNode.pos - p1).norm() + (testNode.pos - p2).norm();
      if (distance < minDistance) {
        minDistance = distance;
        minPShift = shift;
      }
    }
  }

  const GhostNode cornerNode =
      GhostNode(refNode, minPShift, cols, a, currentDeformation);

  // We remove the element from all the nodes it is connected to.
  removeElementFromNodes(elementToMove);

  // overwrite the element
  elementToMove = TElement{(*this),
                           cornerNode,
                           *fixedCoAngleNodes[0],
                           *fixedCoAngleNodes[1],
                           elementToMove.eIndex,
                           elementToMove.noise};
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
    TElement &e = elements[n->elementIndices[i]];
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

void Mesh::setDiagonal(int row, int col, bool useMajorDiagonal) {
  // get the 4 ghost nodes of the selected section of the mesh
  const std::vector<GhostNode> ghosts = getSquareGhostNodes(row, col);
  // Get the indexes of the elements
  auto [e1i, e2i] = getElementIndices(row, col);

  // Now we need to remove e1i and e2i from the 4 nodes, since they will be
  // added back in the tElement constructor
  removeElementsFromNodes(row, col, {e1i, e2i});

  // Update elements, preserving existing noise values
  createElementPair(ghosts, e1i, e2i, useMajorDiagonal, true);
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
  std::vector<Node *> nodes = getSquareNodes(row, col);
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
        if (n->elementIndices[i] == elIndexToRemove[j]) {
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
        tempElementIndices[newCount] = n->elementIndices[i];
        tempNodeIndexInElement[newCount] = n->nodeIndexInElement[i];
        newCount++;
      }
    }

    // Update the node with the new arrays
    for (int i = 0; i < newCount; i++) {
      n->elementIndices[i] = tempElementIndices[i];
      n->nodeIndexInElement[i] = tempNodeIndexInElement[i];
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
}

void Mesh::updateBoundingBox() {
  // Reset the bounding box
  bounds[0] = -INFINITY; // max x
  bounds[1] = INFINITY;  // min x
  bounds[2] = -INFINITY; // max y
  bounds[3] = INFINITY;  // min y

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

void Mesh::updateNrPlasticEvents() {
  // Note, this is effectively the number of plastic events relative to last
  // time this function was called. We rely on the pastM3Nr in the element
  // to be updated in order to find the change since last loading step. For
  // example, if this function is called every 100 loading steps, it will be
  // the number of plasticEvents that have occured in the last 100 steps.
  // (assuming that the mrNr only increases during this period)

  // We also only update this function if the load has changed
  nrPlasticChanges = 0;
  nrPlasticChangesInStep = 0;
  for (size_t i = 0; i < elements.size(); i++) {
    if (elements[i].pastM3Nr != elements[i].m3Nr) {
      nrPlasticChanges += 1;
    }
    if (elements[i].pastStepM3Nr != elements[i].m3Nr) {
      nrPlasticChangesInStep += 1;
    }
  }
}

void Mesh::calculateAverages(bool endOfStep) {

  // Note that totalEnergy has already been calculated since we use it in
  // the energy minimization

  // we reset the maxEnergy
  maxEnergy = 0;
  if (endOfStep) {
    // we update the previous energy
    previousAverageEnergy = averageEnergy;
  }

  // We calculate the center of mass
  com = Vector2d::Zero();

  // We calculate total force for debugging (should be zero)
  Vector2d totalForce = Vector2d::Zero();
  for (int i = 0; i < nodes.size(); i++) {
    totalForce += nodes(i).f;
  }
  // Mathematically, the total force should always be 0, but due to some
  // rounding errors (i think), we need to allow some freedom before we declare
  // that something is wrong.
  if (totalForce.norm() / nrNodes > 1e-13) {
    std::cout << "Min step: " << nrMinItterations << '\n';
    std::cout << "Total force: " << totalForce << '\n';
    std::cout << "Big force: " << totalForce.norm() << '\n';
    if (endOfStep) {

      std::cerr << "Total force is not zero. Something is wrong.";
    }
  }

  // This is the total energy from all the triangles
  double totalRSS = 0;
  for (int i = 0; i < nrElements; i++) {
    TElement e = elements[i];
    totalRSS += e.resolvedShearStress;

    // We also keep track of the highest energy and some other things
    if (e.energy > maxEnergy) {
      maxEnergy = e.energy;
    }
    if (e.m3Nr > maxM3Nr) {
      maxM3Nr = e.m3Nr;
    }
    int plasticChange = e.m3Nr - e.pastM3Nr;
    if (plasticChange > maxPlasticJump) {
      maxPlasticJump = plasticChange;
    } else if (plasticChange < minPlasticJump) {
      minPlasticJump = plasticChange;
    }
    com += e.getCom();
  }

  com /= nrElements;

  // We subtract the ground state energy to make the values a bit nicer
  // (totalEnergy is calculated several times during minimization. This
  // function is called after minimization.)
  averageEnergy = totalEnergy / nrElements;
  delAvgEnergy = (averageEnergy - previousAverageEnergy);
  if (loadSteps == 1) {
    // On the first step, we don't have a previous energy to compare with
    delAvgEnergy = 0;
  }
  averageRSS = totalRSS / nrElements;

  // Update number of plastic events
  updateNrPlasticEvents();
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
}

void Mesh::writeToVtu(std::string filename, bool minimizationStep) {
  // For video making, we need to have accurate bounding boxes.
  updateBoundingBox();

  writeMeshToVtu((*this), simName, dataPath, filename, minimizationStep);
}

void transform(const Matrix2d &matrix, Mesh &mesh,
               std::vector<NodeId> nodesToTransform) {
  // We get the adress of each node
  for (NodeId &nodeId : nodesToTransform) {
    transformInPlace(matrix, *mesh[nodeId]);
  }
}

std::ostream &operator<<(std::ostream &os, const Mesh &mesh) {
  for (int i = mesh.cols - 1; i >= 0; --i) {
    for (int j = 0; j < mesh.rows; ++j) {
      Node n = mesh.nodes(i, j);
      os << n.id.i << "\t";
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
}
void translate(Mesh &mesh, double x, double y) {
  translate(mesh, mesh.fixedNodeIds, x, y);
  translate(mesh, mesh.freeNodeIds, x, y);
}
