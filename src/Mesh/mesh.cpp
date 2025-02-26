#include "mesh.h"
#include "Mesh/node.h"
#include "Mesh/tElement.h"
#include "Simulation/randomUtils.h"
#include <cassert>
#include <iostream>
#include <ostream>
#include <vector>

Mesh::Mesh() {}

// Constructor that initializes the surface with size n x m
Mesh::Mesh(int rows, int cols, double a, double QDSD, bool usingPBC,
           bool useDiagonalFlipping)
    : nodes(rows, cols), a(a), rows(rows), cols(cols), loadSteps(0),
      currentDeformation(Eigen::Matrix2d::Identity()), QDSD(QDSD),
      usingPBC(usingPBC), useDiagonalFlipping(useDiagonalFlipping) {
  // Calculate nrElements based on whether usingPBC is true
  if (usingPBC) {
    nrElements = 2 * rows * cols;
  } else {
    nrElements = 2 * (rows - 1) * (cols - 1);
  }
  nrNodes = rows * cols;

  // Now initialize elements with the calculated size
  elements.resize(nrElements);

  m_createNodes();
  m_updateFixedAndFreeNodeIds();
  m_fillNeighbours();

  // we create some elements here, but note that these will need to be replaced
  // if changes are made to the fixed and free status of the nodes.
  createElements();
}

Mesh::Mesh(int rows, int cols, bool usingPBC)
    : Mesh(rows, cols, 1, 0, usingPBC) {}

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
  currentDeformation = currentDeformation * transformation;
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

Node Mesh::m_getNeighbourNode(const Node &node, int direction) {
  return *(*this)[(*this)[node.id]->neighbours[direction]];
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
  for (size_t i = 0; i < elements.size(); i++) {
    if (elements[i].plasticChange) {
      nrPlasticChanges += 1;
    }
  }
}

// Function to fix the elements of the border vector
void Mesh::fixBorderNodes() {
  fixNodesInRow(0);
  fixNodesInColumn(0);
  fixNodesInRow(rows - 1);
  fixNodesInColumn(cols - 1);
}

void Mesh::fixNodesInRow(int row) {
  for (int col = 0; col < cols; ++col) {
    nodes(row, col).fixedNode = true;
  }
  m_updateFixedAndFreeNodeIds();
}

void Mesh::fixNodesInColumn(int col) {
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
      nodes(row, col).neighbours[LEFT_N] = m_makeNId(row, left);   // left
      nodes(row, col).neighbours[RIGHT_N] = m_makeNId(row, right); // right
      nodes(row, col).neighbours[UP_N] = m_makeNId(up, col);       // up
      nodes(row, col).neighbours[DOWN_N] = m_makeNId(down, col);   // down
    }
  }
}

void Mesh::createElements() {
  // Note that neighbours must be filled before using this function.

  // We construct the elements by finding the four nodes that create two
  // opposing cells, but stopping before we get to the last row and columns if
  // we don't use periodic boundary conditions.
  for (int row = 0; row < ePairRows(); ++row) {
    for (int col = 0; col < ePairCols(); ++col) {
      // We now find the 4 nodes in the current square
      // If we are using PBC, we need to use the neighbours to find the
      // adjacent nodes instead of just incrementing col and row.
      // We want to move right and up. (We start the bottom left corner when
      // indexing)
      Node n1 = *(*this)[m_makeNId(row, col)];
      Node n2 = m_getNeighbourNode(n1, RIGHT_N);
      Node n3 = m_getNeighbourNode(n1, UP_N);
      // n4 is now up AND right of n1
      Node n4 = m_getNeighbourNode(n2, UP_N);

      // This function will shift some nodes across the system if neccecary to
      // create periodic boundaries
      std::vector<GhostNode> ghosts =
          m_makeGhostNodes({n1, n2, n3, n4}, row, col);

      GhostNode g1 = ghosts[0];
      GhostNode g2 = ghosts[1];
      GhostNode g3 = ghosts[2];
      GhostNode g4 = ghosts[3];

      // Determine diagonal direction based on alternating pattern
      bool flipDiagonal = false;
      if (useDiagonalFlipping) {
        flipDiagonal = (row + col) % 2;
      }
      // Picture these as top and bottom triangles. The bottom triangle
      // is element 1 with index e1i. The top triangle is element 2 with
      // index e2i.
      int e1i = 2 * (row * ePairCols() + col); // Triangle 1 index
      int e2i = e1i + 1;                       // Triangle 2 index

      // std::cout << "Nr: " << e1i << (flipDiagonal ? " Right " : " Left ")
      //           << ((e1i % 2 == flipDiagonal) ? "Up" : "Down") << '\n';
      // std::cout << "Nr: " << e2i << (flipDiagonal ? " Right " : " Left ")
      //           << ((e2i % 2 == flipDiagonal) ? "Up" : "Down") << '\n';

      if (flipDiagonal) {
        // Split using diagonal from top-left to bottom-right (↙)
        elements[e1i] = TElement(g1, g2, g4, sampleNormal(1, QDSD));
        elements[e2i] = TElement(g1, g3, g4, sampleNormal(1, QDSD));

        m_addElementIndices({n1, n2, n4}, e1i);
        m_addElementIndices({n1, n3, n4}, e2i);
      } else {
        // Split using diagonal from bottom-left to top-right (↘)
        elements[e1i] = TElement(g1, g2, g3, sampleNormal(1, QDSD));
        elements[e2i] = TElement(g2, g3, g4, sampleNormal(1, QDSD));

        m_addElementIndices({n1, n2, n3}, e1i);
        m_addElementIndices({n2, n3, n4}, e2i);
      }
    }
  }
  // std::cout << (*this);
}

// This is just a function to avoid having to write cols
NodeId Mesh::m_makeNId(int row, int col) { return NodeId(row, col, cols); }

// Helper function to add element indices into nodes
//
// Here we add element indices into nodes so that each node knows what
// elements are connected to it.

// NB! We only need to update the nodes in the "nodes" matrix. The
// duplicate nodes that are given to the elements (n1, n2, n3, n4) do not
// need to be updated. At least I hope we don't need that information in
// those nodes.

// Add indices for the two groups of nodes
// Note here that n1..n4 are the "fake" duplicate nodes beloning to the
// e1i and e2i elements. They do however store the index of the real nodes
// in the nodes matrix. That's why you see that we do nodes(n.id.i) in
// the m_addElementIndices function.
void Mesh::m_addElementIndices(const std::vector<Node> nodeList,
                               int elementIndex) {
  int i = 0;
  for (Node n : nodeList) {
    // Reference to the current count
    int &count = nodes(n.id.i).elementCount;
    // Ensure we don't exceed the array size
    if (count < MAX_ELEMENTS_PER_NODE) {
      nodes(n.id.i).elementIndices[count] = elementIndex;
      nodes(n.id.i).nodeIndexInElement[count] = i;
      ++count; // Increment the count for the node
    } else {
      // Handle overflow (e.g., log an error or take other measures)
      std::cerr << "Error: Too many elements for node " << n.id << std::endl;
    }
    i++;
  }
};

// This function adjusts the position of a node using a shift, also taking into
// acount the current deformation of the system.
Vector2d Mesh::makeGhostPos(const Vector2d &pos, const Vector2d &shift) const {
  return pos + currentDeformation * shift;
}

void Mesh::resetCounters() {
  nrMinimizationItterations = 0;
  nrUpdateFunctionCalls = 0;

  // We also only update this function if the load has changed
  nrPlasticChanges = 0;
  for (size_t i = 0; i < elements.size(); i++) {
    elements[i].pastM3Nr = elements[i].m3Nr;
  }
}

void Mesh::setSimNameAndDataPath(std::string name, std::string path) {
  simName = name;
  dataPath = path;
}

// Helper function to make ghost node
GhostNode Mesh::m_gn(Node n, int row, int col) {
  return GhostNode(n, row, col, cols, a, currentDeformation);
}
GhostNode Mesh::m_gn(Node n) {
  return GhostNode(n, n.id.row, n.id.col, cols, a, currentDeformation);
}

// The idea here is that we have taken our grid, and made rows and columns of
// squares of 4 and 4 nodes. This function converts from reference nodes to
// ghost nodes and makes sure that the ghost nodes are appropreately shifted
// when that is requried
std::vector<GhostNode> Mesh::m_makeGhostNodes(const std::vector<Node> refNodes,
                                              int row, int col) {
  assert(refNodes.size() == 4);

  const Node &n1 = refNodes[0];
  const Node &n2 = refNodes[1];
  const Node &n3 = refNodes[2];
  const Node &n4 = refNodes[3];

  GhostNode gn1 = m_gn(n1);
  GhostNode gn2 = m_gn(n2);
  GhostNode gn3 = m_gn(n3);
  GhostNode gn4 = m_gn(n4);

  if (usingPBC) {

    if (row == rows - 1 && col == cols - 1) {
      // If we are in the corner, we need to move n2, n3 and n4
      gn2 = m_gn(n2, n2.id.row, cols);
      gn3 = m_gn(n3, rows, n3.id.col);
      gn4 = m_gn(n4, rows, cols);
    } else if (col == cols - 1) {
      // If we are in the last column, we need to move n2 and n4
      gn2 = m_gn(n2, n2.id.row, cols);
      gn4 = m_gn(n4, n4.id.row, cols);
    } else if (row == rows - 1) {
      // If we are in the last row, we need to move n3 and n4
      gn3 = m_gn(n3, rows, n3.id.col);
      gn4 = m_gn(n4, rows, n4.id.col);
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
    for (size_t j = 0; j < e.tElementNodes.size(); j++) {
      if (realId) {
        std::cout << e.tElementNodes[j].referenceId.i << sep;
      } else {
        std::cout << e.tElementNodes[j].ghostId.i << sep;
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

void Mesh::applyForceFromElementsToNodes() {
  // Loop over all the nodes
  // (Looping over the elements would create a problem since two threads might
  // write to the same node at the same time. By looping over the nodes, each
  // thread only writes to one node at a time. There is no problem if two
  // threads read from the same element at the same time.)
#pragma omp parallel for
  for (int i = 0; i < (int)nodes.size(); ++i) {
    Node &n = nodes(i);
    n.resetForce();
    for (size_t e = 0; e < n.elementCount; e++) {
      // This is for debugging and readability
      // These variables should be optimized away
      int elementNr = n.elementIndices[e];
      int nodeNrInElement = n.nodeIndexInElement[e];

      // -1 is the default value and means the element has not been assigned
      if (nodeNrInElement != -1) {
        TElement &element = elements[elementNr];
        GhostNode &elementNode = element.tElementNodes[nodeNrInElement];
        n.f += elementNode.f;
      } else {
        std::cerr << "Something is wrong in the meshing!" << std::endl;
      }
    }
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
      if (elements[i].tElementNodes[j].pos[0] > bounds[0])
        bounds[0] = elements[i].tElementNodes[j].pos[0]; // max x
      if (elements[i].tElementNodes[j].pos[0] < bounds[1])
        bounds[1] = elements[i].tElementNodes[j].pos[0]; // min x

      // Update bounds for y-coordinate
      if (elements[i].tElementNodes[j].pos[1] > bounds[2])
        bounds[2] = elements[i].tElementNodes[j].pos[1]; // max y
      if (elements[i].tElementNodes[j].pos[1] < bounds[3])
        bounds[3] = elements[i].tElementNodes[j].pos[1]; // min y
    }
  }
}

void Mesh::calculateAverages() {

  // Note that totalEnergy has already been calculated since we use it in
  // the energy minimization

  // we reset the maxEnergy
  maxEnergy = 0;
  // we update the previous energy
  previousAverageEnergy = averageEnergy;

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
  }

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
