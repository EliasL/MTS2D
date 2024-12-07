#include "mesh.h"
#include "Mesh/node.h"
#include "Simulation/randomUtils.h"
#include <iostream>

Mesh::Mesh() {}

// Constructor that initializes the surface with size n x m
Mesh::Mesh(int rows, int cols, double a, double QDSD, bool usingPBC)
    : nodes(rows, cols), a(a), rows(rows), cols(cols), loadSteps(0),
      currentDeformation(Eigen::Matrix2d::Identity()), QDSD(QDSD),
      usingPBC(usingPBC) {
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

bool Mesh::isFixedNode(NodeId nodeId) { return (*this)[nodeId]->fixedNode; }

void Mesh::addLoad(double loadChange) {
  load += loadChange;
  loadSteps++;
}

void Mesh::applyTransformation(Matrix2d transformation) {
  // We get all the nodes in the mesh.
  for (long i = 0; i < nodes.size(); i++) {
    transformInPlace(transformation, nodes(i));
  }
  // We also assume we want to transform the current deformation
  applyTransformationToSystemDeformation(transformation);
}

void Mesh::applyTransformationToFixedNodes(Matrix2d transformation) {
  // We get the id of each node in the border
  for (NodeId &nodeId : fixedNodeIds) {
    transformInPlace(transformation, *(*this)[nodeId]);
  }
}

void Mesh::applyTransformationToSystemDeformation(Matrix2d transformation) {
  currentDeformation = currentDeformation * transformation;
}

void Mesh::applyTranslation(Vector2d displacement) {
  for (long i = 0; i < nodes.size(); i++) {
    translateInPlace(nodes(i), displacement);
  }
}

void Mesh::setInitPos() {
  for (long i = 0; i < nodes.size(); i++) {
    nodes(i).setInitPos(nodes(i).pos());
  }
}

void Mesh::resetForceOnNodes() {
  for (long i = 0; i < nodes.size(); i++) {
    nodes(i).resetForce();
  }
}

Node Mesh::m_getNeighbourNode(Node node, int direction) {
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
    elements[i].updatePastM3Nr();
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
  int rows = this->rows;
  int cols = this->cols;

  // If we are not using PBC, we skipp creating the last elements
  if (!usingPBC) {
    rows -= 1;
    cols -= 1;
  }

  // We construct the elements by finding the four nodes that create two
  // opposing cells, but stopping before we get to the last row and columns if
  // we don't use periodic boundary conditions.
  for (int row = 0; row < rows; ++row) {
    for (int col = 0; col < cols; ++col) {
      // We now find the 4 nodes in the current square
      Node n1 = *(*this)[m_makeNId(row, col)];
      // If we are using PBC, we need to use the neighbours to find the
      // adjacent nodes instead of just incrementing col and row.
      // We want to move right and up. (We start the bottom left corner when
      // indexing)
      Node n2 = m_getNeighbourNode(n1, RIGHT_N);
      Node n3 = m_getNeighbourNode(n1, UP_N);
      // n4 is now up AND right of n1
      Node n4 = m_getNeighbourNode(n2, UP_N);

      if (usingPBC) {
        if (row == rows - 1 && col == cols - 1) {
          // If we are in the corner, we need to move n2, n3 and n4
          m_makeGN(n2, n2.id.row, cols);
          m_makeGN(n3, rows, n3.id.col);
          m_makeGN(n4, rows, cols);
        } else if (col == cols - 1) {
          // If we are in the last column, we need to move n2 and n4
          m_makeGN(n2, n2.id.row, cols);
          m_makeGN(n4, n4.id.row, cols);
        } else if (row == rows - 1) {
          // If we are in the last row, we need to move n3 and n4
          m_makeGN(n3, rows, n3.id.col);
          m_makeGN(n4, rows, n4.id.col);
        }
      }

      // Picture these as top and bottom triangles. The bottom
      // triangle is element 1 with index e1i.
      int e1i = 2 * (row * cols + col); // Triangle 1 index
      int e2i = e1i + 1;

      elements[e1i] = TElement(n1, n2, n3, sampleNormal(1, QDSD));
      elements[e2i] = TElement(n2, n3, n4, sampleNormal(1, QDSD));

      // Add element indices into nodes so that each node knows what elements
      // it is surrounded by
      int i = 0;
      for (Node n : {n1, n2, n3}) {
        // NB! We need to use the nodes in the mesh!
        // Not these copies that we use here (n1, n2, n3, n4)
        nodes(n.id.i).elementIndices.push_back(e1i);
        nodes(n.id.i).nodeIndexInElement.push_back(i++);
      }
      i = 0;
      for (Node n : {n2, n3, n4}) {

        // NB! We need to use the nodes in the mesh!
        // Not these copies that we use here (n1, n2, n3, n4)
        nodes(n.id.i).elementIndices.push_back(e2i);
        nodes(n.id.i).nodeIndexInElement.push_back(i++);
      }
    }
  }
}

// This is just a function to avoid having to write cols
NodeId Mesh::m_makeNId(int row, int col) { return NodeId(row, col, cols); }

// This function adjusts the position of a node using a shift, also taking into
// acount the current deformation of the system.
Vector2d Mesh::makeGhostPos(Vector2d pos, Vector2d shift) const {
  return pos + currentDeformation * shift;
}

void Mesh::resetLoadingStepFunctionCounters() {
  nrMinimizationItterations = 0;
  nrUpdateFunctionCalls = 0;
}

void Mesh::setSimNameAndDataPath(std::string name, std::string path) {
  simName = name;
  dataPath = path;
}

void Mesh::m_makeGN(Node &n, int newRow, int newCol) {
  // Think carefully about where the initial position of the ghost node should
  // be. Is it shifted by the current deformation or not? This might change
  // depending on what you want to simulate. Currently, I think the initial
  // position should not be affected by the current deformation.

  n.isGhostNode = true;
  n.ghostId = NodeId(newRow, newCol, cols + 1);

  // We now need to shift the position.
  // We subtract the old pos from the new pos to get the shift we want to apply
  n.ghostShift = {
      (newCol - n.id.col) * a,
      (newRow - n.id.row) * a,
  };

  n.setInitPos(n.init_pos() + n.ghostShift);
  n.setPos(makeGhostPos(n.pos(), n.ghostShift));
  // std::cout << n << '\n';
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
    for (size_t j = 0; j < e.nodes.size(); j++) {
      if (realId) {
        std::cout << e.nodes[j].id.i << sep;
      } else {
        std::cout << e.nodes[j].ghostId.i << sep;
      }
    }
    std::cout << end;
  }
  std::cout << '\n';
}

void Mesh::updateElements() {
// Parallel loop with reduction clauses for min and max
#pragma omp parallel for
  for (int i = 0; i < nrElements; i++) {
    elements[i].update(*this);
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
      if (elements[i].nodes[j].pos()[0] > bounds[0])
        bounds[0] = elements[i].nodes[j].pos()[0]; // max x
      if (elements[i].nodes[j].pos()[0] < bounds[1])
        bounds[1] = elements[i].nodes[j].pos()[0]; // min x

      // Update bounds for y-coordinate
      if (elements[i].nodes[j].pos()[1] > bounds[2])
        bounds[2] = elements[i].nodes[j].pos()[1]; // max y
      if (elements[i].nodes[j].pos()[1] < bounds[3])
        bounds[3] = elements[i].nodes[j].pos()[1]; // min y
    }
  }
}

void Mesh::applyForceFromElementsToNodes() {
#pragma omp parallel for
  // Loop over all the nodes
  for (int i = 0; i < nodes.size(); ++i) {
    Node *n = &nodes.data()[i];
    // Loop over all elements connected to a node
    for (int e = 0; e < 6; e++) {
      int elementNr = n->elementIndices[e];
      int nodeNrInElement = n->nodeIndexInElement[e];
      n->f += elements[elementNr].nodes[nodeNrInElement].f;
    }
  }
}

void Mesh::calculateAverages() {
  // we reset the maxEnergy
  maxEnergy = 0;
  // we update the previous energy
  previousAverageEnergy = averageEnergy;

  // This is the total energy from all the triangles
  double totalEnergy = 0;
  double totalRSS = 0;
  for (int i = 0; i < nrElements; i++) {
    totalEnergy += elements[i].energy;
    totalRSS += elements[i].resolvedShearStress;

    // We also keep track of the highest energy value.
    if (elements[i].energy > maxEnergy) {
      maxEnergy = elements[i].energy;
    }
  }

  // We subtract the ground state energy to make the values a bit nicer
  averageEnergy = totalEnergy / nrElements;
  delAvgEnergy = (averageEnergy - previousAverageEnergy);
  if (loadSteps == 1) {
    // On the first step, we don't have a previous energy to compare with
    delAvgEnergy = 0;
  }
  averageRSS = totalRSS / nrElements;
}

double Mesh::calculateTotalEnergy() {

  double totalEnergy = 0;
  for (int i = 0; i < nrElements; i++) {
    totalEnergy += elements[i].energy;
  }
  return totalEnergy;
}

// Helper function to update positions using a generic buffer and its size
void Mesh::updateNodePositions(const double *data, size_t length) {
  // The displacement is structed like this: [x1,x2,x3,x4,y1,y2,y3,y4], so we
  // need to know where the "x" end and where the "y" begin.
  int nr_x_values = length / 2;
  Node *n;

  // We loop over all the free nodes
  for (size_t i = 0; i < freeNodeIds.size(); i++) {
    n = (*this)[freeNodeIds[i]];
    // This function changes the position of the node based on the given
    // displacement and the current initial position.
    n->setDisplacement({data[i], data[i + nr_x_values]});
  }
}

// Updates the forces on the nodes in the surface and returns the total
// energy from all the elements in the surface.
void Mesh::updateMesh() {
  // First of all we need to make sure that the forces on the nodes have been
  // reset
  resetForceOnNodes();

  // Now we update all the elements using the current positions of the nodes
  updateElements();

  // We then add the force from the elements back to the nodes
  applyForceFromElementsToNodes();
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
      os << n.f << "\t";
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
