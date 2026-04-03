#ifndef MESH_H
#define MESH_H
#include <limits>
#pragma once

#include "Data/cereal_help.h"
#include "Eigen/Core"
#include "compare_macros.h"
#include "node.h"
#include "tElement.h"
#include <array>
#include <cereal/types/vector.hpp>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <type_traits>

// The neighbours should be indexed using these defines for added readability
#define LEFT_N Vector2i(-1, 0)
#define RIGHT_N Vector2i(1, 0)
#define UP_N Vector2i(0, 1)
#define DOWN_N Vector2i(0, -1)

enum class VtuFieldLevel;

/**
 * @brief The Mesh class represents a 2D surface mesh of nodes.
 *
 * It tracks nodes at the surface's border, neighbors, and manages the creation
 * of triangular elements formed by neighboring nodes.
 *
 * NB: The nodes in the elements are seperate from the nodes in the mesh. This
 * means that you should modify the nodes in the mesh before creating the
 * elements. The elements will then create copies of the nodes with the
 * modifications.
 */
class Mesh {
public:
  enum class UpdateState { Dirty, Forces, Geometry, Full };
  // A matrix representing the mesh of nodes.
  // Using RowMajored makes more natural indexing:
  /*
    6   7   8
    3   4   5
    0   1   2
  */
  Matrix<Node, Dynamic, Dynamic, Eigen::RowMajor> nodes;

  // A collection of elements formed by triangles of nodes.
  // We use pointers so that the nodes in the elements can have a reference
  // to the mesh. Otherwise, we would need to set the reference at
  // initialization.
  std::vector<TElement> elements;

  struct FlipRecord {
    int e1 = -1;
    int e2 = -1;
    std::array<Vector2i, 4> nodeIds;
  };
  std::vector<FlipRecord> flipRecords;
  int flipRecordCount = 0;
  std::vector<uint8_t> flippedThisReconnect;

  struct EdgeKey {
    int a = 0;
    int b = 0;
    EdgeKey() = default;
    EdgeKey(int i, int j) {
      if (i <= j) {
        a = i;
        b = j;
      } else {
        a = j;
        b = i;
      }
    }
    bool operator==(const EdgeKey &o) const noexcept {
      return a == o.a && b == o.b;
    }
  };

  struct EdgeKeyHash {
    std::size_t operator()(const EdgeKey &k) const noexcept {
      const auto hi = static_cast<uint64_t>(static_cast<uint32_t>(k.a));
      const auto lo = static_cast<uint64_t>(static_cast<uint32_t>(k.b));
      return static_cast<std::size_t>((hi << 32) ^ lo);
    }
  };

  using EdgeSet = std::unordered_set<EdgeKey, EdgeKeyHash>;

  // Tracks internal edges (shared by two elements) for flip counting.
  EdgeSet currentEdgeSet;
  EdgeSet previousStepEdgeSet;
  std::size_t edgeFlipsSinceLastStep = 0;

  // IDs of nodes that are on the border of the mesh.
  std::vector<NodeId> fixedNodeIds;

  // IDs of nodes that are not on the border of the mesh.
  std::vector<NodeId> freeNodeIds;

  // The characteristic dimension of the mesh.
  double a;
  // Fixed reference area for each element. We use a constant reference area
  // because after reconnecting, the reference nodes may end up being aligned
  // in a straight line, resulting in a reference area=0.
  double init_element_area = 0.0;

  // Number of rows of nodes in the mesh.
  int rows;
  // Number of columns of nodes in the mesh.
  int cols;

  // This is the number rows and columns of pairs of element in the mesh
  // Consider grouping up two and two triangular elements to form squares,
  // The number of rows and columns of this unit depends on whether or not
  // we are using periodic boundary conditions. This is useful when constructing
  // the elements.
  int ePairCols() const { return usingPBC ? cols : cols - 1; }
  int ePairRows() const { return usingPBC ? rows : rows - 1; }

  // The applied load on the mesh.
  // This variable is not used for physics. The physics are solely based on
  // the position of the boundary nodes. This value is stored for logging
  // purposes.
  double load;

  // Number of load steps taken.
  int loadSteps;

  // If we want to shear the entire mesh, we use the periodic transform, but
  // if we want to create a change in the current load without moving any of
  // the nodes, we use the periodic load (which applies only to the "distances"
  // between the periodically repeated systems).

  // We need to know how to tile the system periodically. This transformation
  // is applied to the diplacement of the periodic nodes.
  Matrix2d currentDeformation;

  // The number of triangles created in the mesh.
  int nrElements;
  // Nr of nodes
  int nrNodes;

  // We calculate the total energy during the simulation, and the average
  // energy is useful to plot, so we keep this value here for easy access.
  double totalEnergy = 0;
  double averageEnergy = 0;
  // This might also be usefull
  double maxEnergy = 0;
  double maxForce = 0;          // Max force component in mesh.
  double averageP11 = 0;        // FirstPiola[0,0]
  double averageP12 = 0;        // FirstPiola[0,1]
  double averageP21 = 0;        // FirstPiola[1,0]
  double averageP22 = 0;        // FirstPiola[1,1]
  double averageSigma11 = 0;    // cauchy stress[0,0]
  double averageSigma12 = 0;    // cauchy stress[0,1]
  double averageSigma22 = 0;    // cauchy stress[1,1]
  double averageSigmaTrace = 0; // cauchy stress[0,0] + cauchy stress[1,1]
  int maxM3Nr = 0;
  int sumM3Nr = 0;
  int maxPlasticJump = 0;
  int minPlasticJump = 0;
  std::array<int, 4> redQuadrantCounts = {0, 0, 0, 0};
  std::array<int, 4> redQuadrantFixedCounts = {0, 0, 0, 0};
  Vector2d com = {0, 0}; // Center of mass

  // Number of plastic changes is last loading step.
  int nr_elements_with_m3_fix_change = 0;
  // Number of plastic changes during the minimization of a single step.
  int nr_elements_with_m3_fix_changeInStep = 0;

  // Controls the standard deviation of the quenched dissorder in the mesh.
  double QDSD = 0;

  // Energy function settings used when creating elements.
  std::string energyFunction = "contiSquare";
  double bulkModulus = 4.0;
  // Lattice basis that maps (col,row) to physical coordinates.
  Matrix2d latticeBasis = Matrix2d::Identity();

  // Flag for using periodic or fixed boundary conditions.
  bool usingPBC;

  bool meshReconnected = false;

  // Flag for diagonal meshing.
  std::string diagonal;

  // This is the number of iterations the mesh has gone through in the current
  // loading step.
  int nrMinItterations = 0;

  // This is the number of update function calls the minimuzation algorithm has
  // used in the current loading step.
  int nrMinFunctionCalls = 0;

  // These are sometimes convenient to access through the mesh instead of the
  // simulation, so they are stored here as well.
  std::string simName = "defaultName";
  std::string dataPath = "";

  // the bounding rectangle of the mesh: max x, min x, max y, min y
  std::array<double, 4> bounds = {-1000000, 1000000, -1000000, 1000000};

  // Default constructor.
  Mesh();

  // Constructor to initialize the mesh with a specified number of rows,
  // columns, and characteristic dimension.
  Mesh(int rows, int cols, double a = 1, double QDSD = 0, bool usingPBC = true,
       std::string diagonal = "major",
       std::string energyFunction = "contiSquare", double bulkModulus = 4);

  Mesh(int rows, int cols, bool usingPBC, std::string diagonal);

  Mesh(int rows, int cols, bool usingPBC);

  // Update latticeBasis based on energyFunction and a.
  void updateLatticeBasis();

  // Overloaded indexing operator to access nodes by their NodeId.
  Node *operator[](const NodeId &id) { return &nodes(id.i); }
  // Const overloaded indexing operator to access nodes by their NodeId.
  const Node *operator[](const NodeId &id) const { return &nodes(id.i); }

  // Overloaded indexing operator to access nodes by their NodeId.
  Node *operator[](GhostNode &gn) { return &nodes(gn.referenceId.i); }

  // Overloaded indexing operator to access nodes by their NodeId.
  const Node *operator[](const GhostNode &gn) const {
    return &nodes(gn.referenceId.i);
  }

  // Determines if a node is at the border of the mesh.
  bool isFixedNode(const NodeId &n_id) const;

  // Note that this load variable is ONLY for logging. It does not affect the
  // physics of the simulation.
  // This adds a load to the mesh load variable, but also increases the
  // load steps counter. Therefore, this function should always be used
  // when increasing the load during a step.
  void addLoad(double loadChange);

  // Applies a transform to all nodes in the mesh, including the PBC.
  void applyTransformation(const Matrix2d &transformation);

  // Applies a transform to the border nodes.
  void applyTransformationToFixedNodes(const Matrix2d &transformation);

  // Applies a transform to the periodic boundary tranform.
  // (see how it affects the pos function in PeriodicNode. )
  void applyTransformationToSystemDeformation(const Matrix2d &transformation);

  // Apply translation to all nodes in the mesh.
  void applyTranslation(const Vector2d &displacement);

  // This sets the current position as the initial position of the mesh.
  void setInitPos();

  // Calculates averages.
  double averageResolvedShearStress() const;

  // Fixes the border nodes in the mesh.
  void fixBorderNodes();

  // Fixes the node in the bottom left corner of the mesh.
  void fixBottomLeftCorner();

  // Fixes the nodes in a given row.
  void fixNodesInRow(int row);

  // Fixes the nodes in a given column.
  void fixNodesInColumn(int column);

  // Print element connectivity (for debugging).
  void printConnectivity(bool realId = true);

  // Retrives the node at a given position. (with PBC wrapping if applicable)
  Node *getNode(const Vector2i &rowCol);
  // Creates a ghost node at a given position. Position will not be wrapped, but
  // still use the correct reference node.
  GhostNode getGhostNode(const Vector2i &rowCol);
  // Retrives the neighbour of a node at a given direction.
  Node *getNeighbourNode(const Node &node, const Vector2i &direction);

  // Gets 4 nodes from grid
  std::vector<Node *> getSquareNodes(int row, int col);

  // Gets 4 nodes from grid and converts to ghost nodes
  std::vector<GhostNode> getSquareGhostNodes(int row, int col);

  // Calculates element indices for a given row and column
  std::pair<int, int> getElementIndices(int row, int col);

  // Similar to getSquareGhostNodes, but using elements instead of row and col
  std::array<GhostNode, 4> getElementPairNodes(const TElement &e1,
                                               const TElement &e2);

  std::vector<GhostNode> getUniqueNodes(const std::vector<TElement *> elements);

  // Creates or updates two triangular elements based on the specified diagonal
  // direction
  void createElementPair(const std::vector<GhostNode> &ghosts, int e1i, int e2i,
                         bool useMajorDiagonal, bool preserveNoise = false);
  void createElementPair(const std::vector<const GhostNode *> &ghostsPtr,
                         int e1i, int e2i, bool useMajorDiagonal,
                         bool preserveNoise = false);
  void createElementPair(const std::array<const GhostNode *, 4> &ghostsPtr,
                         int e1i, int e2i, bool useMajorDiagonal,
                         bool preserveNoise = false);

  // Creates triangles from neighboring nodes to form the elements of the mesh.
  void createElements();

  // When we create tElements, it's index is added to the nodes so the nodes
  // know what elements they are connected to. When we reconnect, we need to
  // remove these connections.
  void removeElementsFromNodes(int row, int col,
                               const std::vector<int> elementIndexes);

  void removeElementsFromNodes(std::vector<Node *> nodes,
                               const std::vector<int> elIndexToRemove);

  void removeElementsFromNodes(const std::vector<const GhostNode *> gNodes,
                               const std::vector<int> elIndexToRemove);

  void removeElementFromNodes(const TElement &element);

  // This reconnectes one pair of triangles (4 nodes) to have their diagonal in
  // the specified direction.
  void setDiagonal(int row, int col, bool useMajorDiagonal);

  // This function goes through each pair of elements and checks if the pair
  // should flip their diagonal
  MTS_NOINLINE bool reconnect(bool onlyCheck = false);
  void resetFlipTracking();
  void undoReconnections();

  MTS_NOINLINE void reconnectDelaunay();

  // This function takes two elements that should both have large angles, and
  // reconfigures the 4 nodes into two new elements that have smaller angles.
  void fixElementPair(TElement &e1, TElement &e2);

  // Sometimes, the element pair will be accross the periodic boundary. This
  // case needs special care
  void fixPeriodicElementPair(TElement &e1, TElement &e2);

  // counts the number of elements connected to a specific ghost node. This is
  // different from the number of elements connected to the reference node of
  // that ghost node
  int countConnectionsInGhostNode(const GhostNode &gn);

  // This function moves one element so that it is appropreately next to another
  // element
  void moveElementToTwin(TElement &elementToMove, const TElement &fixedElement);

  // Uses information from the nodes to recreate a mesh structure.
  void recreateElements();

  // Loops over all elements, updates them, and accumulates forces/energy.
  void updateElements();
  // Marks cached element data as stale after node changes.
  void markDirty();
  // Ensure element forces/energy are up to date.
  void ensureForces();
  // Ensure element geometry (G/angleNode) is up to date.
  void ensureGeometry();
  // Ensure full element data (cauchy stress/angles) is up to date.
  void ensureFull();
  UpdateState getUpdateState() const { return updateState; }

  // Updates the node positions using the data array.
  void updateNodePositions(const double *data, size_t length);
  // Saves node displacements into the provided buffer (x block then y block).
  void saveNodeDisplacements(std::vector<double> &displacements) const;
  // Restores node displacements from the provided buffer.
  void restoreNodeDisplacements(const std::vector<double> &displacements);

  struct Snapshot {
    std::vector<double> displacements;
    double load = 0.0;
    int loadSteps = 0;
    Matrix2d currentDeformation = Matrix2d::Identity();
  };

  void saveSnapshot(Snapshot &snapshot) const;
  Snapshot snapshotState() const;
  void restoreState(const Snapshot &snapshot);
  double rmsDistanceToSnapshot(const Snapshot &snapshot,
                               bool subtractCom = true) const;

  // Loops through all elements connected to node and updates the force
  void updateNodeForce(Node &node);

  // Uses the ids in the elements to update the force on the nodes.
  void applyForceFromElementsToNodes();

  // Calculates averages and updates plastic event counters.
  // Should only be used AFTER minimization.
  MTS_NOINLINE void updateAveragesAndPlasticEvents();

  void updateCom();

  // Reset forces, update elements, calculate forces and energy.
  void updateMesh();

  // This function should be called at the end of each loading step to reset
  // the counters keeping track of how many times things have been called.
  // And some other stuff. Important to call after writing also.
  void resetCounters();

  void resetPastPlasticCount(bool endOfStep = true);

  void setSimNameAndDataPath(std::string name, std::string path);

  void updateBoundingBox();
  void updateAngles();

  void moveMeshSection(double minX, double minY, Vector2d disp,
                       bool moveFixed = true, bool moveFree = false,
                       double maxX = std::numeric_limits<double>().max(),
                       double maxY = std::numeric_limits<double>().max());

  void writeToVtu(std::string filename = "", bool minimizationStep = false);
  MTS_NOINLINE void writeToVtu(std::string filename, bool minimizationStep,
                               VtuFieldLevel level,
                               std::string nameSuffix = "");

private:
  // Fills in the IDs of nodes that are not at the border.
  void m_updateFixedAndFreeNodeIds();

  // Sets the positions of nodes in the mesh based on mesh dimensions and
  // spacing.
  void m_createNodes();

  // Rebuild per-node element connectivity from the current elements.
  void rebuildConnectivity();
  // Rebuild cached ghost-node pointers after load/reconnect.
  void rebuildConnectedGhostNodes();
  // Update cached ghost-node pointers for a single element.
  void updateConnectedGhostNodesForElement(int elementIndex);
  // Update elements with force/energy only (fast path).
  MTS_NOINLINE void updateElementsForces();
  // Update elements with geometry needed for reconnect (G + angleNode).
  MTS_NOINLINE void updateElementsGeometry();
  // Update elements with derived fields (F/C/m/sigma/angles).
  MTS_NOINLINE void updateElementsFull();
  void rebuildCurrentEdgeSet();
  void initializeEdgeSets();
  void updateEdgeSetForFlip(const std::array<const GhostNode *, 4> &nodeOrder,
                            bool useMajorDiagonal);
  std::size_t countEdgeSetDelta(const EdgeSet &lhs, const EdgeSet &rhs) const;
  // Flip a pair of elements using the provided node order and diagonal choice.
  void flipEdge(int e1i, int e2i,
                const std::array<const GhostNode *, 4> &nodeOrder,
                bool useMajorDiagonal, bool preserveNoise = false);
  void recordFlipPair(int e1i, int e2i,
                      const std::array<const GhostNode *, 4> &nodeOrder);
  const GhostNode *findGhostInPairById(int e1i, int e2i,
                                       const Vector2i &id) const;
  void removeElementsFromNodes(const std::array<const GhostNode *, 4> &gNodes,
                               int e1i, int e2i);
  void removeElementsFromNodes(const std::array<Node *, 4> &nodes, int e1i,
                               int e2i);

  // Creates the NodeId of a node at a given position.
  NodeId m_makeNId(int row, int col);

  // Helper function to make ghost node

  // This function automatically sets the periodic shift based on row and col
  GhostNode m_gn(const Node *n, int row, int col);
  // This function automatically sets the periodic shift based on targetPos
  GhostNode m_gn(const Node *n, const Vector2d targetPos);
  // This function uses no shift and mirrors the reference nodes position
  GhostNode m_gn(const Node *n);

  // Creates ghost nodes from reference nodes.
  std::vector<GhostNode>
  m_makeGhostNodes(const std::vector<Node *> referenceNodes, int row, int col);

  // Debugging function to confirm that forces are low after reconnecting
  void checkForces(std::vector<Node *> nodes);
  void checkForces(const std::vector<GhostNode> nodes);

  // Scratch buffer for per-thread force accumulation.
  std::vector<Vector2d> forceScratch;
  int forceScratchThreads = 0;
  UpdateState updateState = UpdateState::Dirty;

  friend class cereal::access; // Necessary to serialize private members.
  template <class Archive> void serialize(Archive &ar);
};

std::ostream &operator<<(std::ostream &os, const Mesh &mesh);

/**
 * @brief Transforms all nodes in a mesh by applying a transformation matrix.
 *
 * This function applies a linear transformation defined by a matrix to the
 * node's position.
 *
 * @param matrix The transformation matrix to apply.
 * @param mesh The mesh to transform.
 */
void transform(const Matrix2d &matrix, Mesh &mesh);
// Only transform nodes in the provided list.
void transform(const Matrix2d &matrix, Mesh &mesh,
               std::vector<NodeId> nodesToTransform);

/**
 * @brief Translates all nodes in a mesh by a given displacement.
 *
 * This function adds a displacement to the node's position, with an optional
 * multiplier to scale the displacement.
 *
 * @param mesh The original node to be translated.
 * @param x The displacement in the x direction.
 * @param y The displacement in the y direction.
 */
void translate(Mesh &mesh, double x, double y);
// Only translate nodes in the provided list.
void translate(Mesh &mesh, std::vector<NodeId> nodesToTranslate, double x,
               double y);

template <class Archive> void Mesh::serialize(Archive &ar) {
  // Serialize fields using the MAKE_NVP macro.
  ar(MAKE_NVP(nodes), MAKE_NVP(elements), MAKE_NVP(fixedNodeIds),
     MAKE_NVP(freeNodeIds), MAKE_NVP(a), MAKE_NVP(rows), MAKE_NVP(cols),
     MAKE_NVP(load), MAKE_NVP(loadSteps), MAKE_NVP(currentDeformation),
     MAKE_NVP(nrElements), MAKE_NVP(nrNodes), MAKE_NVP(totalEnergy),
     MAKE_NVP(averageEnergy), MAKE_NVP(maxEnergy), MAKE_NVP(QDSD));

  using ArchiveT = std::decay_t<Archive>;
  if constexpr (ArchiveT::is_loading::value) {
    int tmpChange = -1;
    loadWithDefault(ar, "nr_elements_with_m3_fix_change", tmpChange, -1);
    if (tmpChange == -1) {
      loadWithDefault(ar, "nrPlasticChanges", tmpChange, 0);
    }
    nr_elements_with_m3_fix_change = tmpChange;

    int tmpChangeStep = -1;
    loadWithDefault(ar, "nr_elements_with_m3_fix_changeInStep", tmpChangeStep,
                    -1);
    if (tmpChangeStep == -1) {
      loadWithDefault(ar, "nrPlasticChangesInStep", tmpChangeStep, 0);
    }
    nr_elements_with_m3_fix_changeInStep = tmpChangeStep;
  } else {
    ar(MAKE_NVP(nr_elements_with_m3_fix_change),
       MAKE_NVP(nr_elements_with_m3_fix_changeInStep));
  }

  ar(MAKE_NVP(usingPBC), MAKE_NVP(nrMinItterations),
     MAKE_NVP(nrMinFunctionCalls), MAKE_NVP(simName), MAKE_NVP(dataPath),
     MAKE_NVP(bounds));

  // Load fields with default values if they are missing from the archive.
  // This is important for backwards compatibility when loading older dumps.
  // In theory, all the fields could be loaded with default values, but this
  // is perhaps more controlled.
  LOAD_WITH_DEFAULT(ar, energyFunction, std::string("contiSquare"));
  LOAD_WITH_DEFAULT(ar, bulkModulus, 4.0);
  LOAD_WITH_DEFAULT(ar, maxM3Nr, 0);
  LOAD_WITH_DEFAULT(ar, sumM3Nr, 0);
  LOAD_WITH_DEFAULT(ar, maxPlasticJump, 0);
  LOAD_WITH_DEFAULT(ar, minPlasticJump, 0);
  LOAD_WITH_DEFAULT(ar, maxForce, 0.0);
  LOAD_WITH_DEFAULT(ar, averageP11, 0.0);
  LOAD_WITH_DEFAULT(ar, averageP12, 0.0);
  LOAD_WITH_DEFAULT(ar, averageP21, 0.0);
  LOAD_WITH_DEFAULT(ar, averageP22, 0.0);
  LOAD_WITH_DEFAULT(ar, averageSigma11, 0.0);
  LOAD_WITH_DEFAULT(ar, averageSigma12, 0.0);
  LOAD_WITH_DEFAULT(ar, averageSigma22, 0.0);
  LOAD_WITH_DEFAULT(ar, averageSigmaTrace, 0.0);
  ar(MAKE_NVP(com));

  // Reset flip tracking (not serialized).
  flipRecordCount = 0;
  flipRecords.resize(nrElements / 2);
  flippedThisReconnect.resize(nrElements);
  std::fill(flippedThisReconnect.begin(), flippedThisReconnect.end(), 0);

  if constexpr (Archive::is_loading::value) {
    updateLatticeBasis();
    for (auto &e : elements) {
      e.setInitArea(init_element_area);
    }
    rebuildConnectedGhostNodes();
    initializeEdgeSets();
    updateState = UpdateState::Dirty;
  }
}

// Cereal save function for matrices
// https://stackoverflow.com/questions/22884216/serializing-eigenmatrix-using-cereal-library
namespace cereal {

template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options,
          int _MaxRows, int _MaxCols>
inline void save(Archive &ar,
                 const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows,
                                     _MaxCols> &m) {
  // Serialize matrix dimensions
  ar(cereal::make_nvp("rows", m.rows()), cereal::make_nvp("cols", m.cols()));

  // Store data as a vector for compatibility with XML/JSON
  std::vector<_Scalar> data(m.data(), m.data() + m.size());
  ar(cereal::make_nvp("data", data));
}

template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options,
          int _MaxRows, int _MaxCols>
inline void
load(Archive &ar,
     Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> &m) {
  int rows, cols;
  ar(cereal::make_nvp("rows", rows), cereal::make_nvp("cols", cols));

  // Load vector data
  std::vector<_Scalar> data;
  ar(cereal::make_nvp("data", data));

  // Resize Eigen matrix and copy data
  m.resize(rows, cols);
  std::copy(data.begin(), data.end(), m.data());
}

} // namespace cereal

/*
   Internal helper function that compares two Mesh objects field by field.
   It takes an optional std::string pointer (debugMsg) that, if not null,
   collects messages for any differences found.
   When debugMsg is nullptr, it simply returns whether the objects are equal.
*/
inline bool compareconnectesInternal(const Mesh &lhs, const Mesh &rhs,
                                     std::string *debugMsg = nullptr,
                                     int tabNumber = 0) {
  bool equal = true;

  // Compare the node matrix.
  COMPARE_FIELD(nodes);

  // Compare the vector of TElements (which calls TElement::operator==
  // internally).
  COMPARE_FIELD(elements);

  // Compare vector<NodeId> fixedNodeIds, vector<NodeId> freeNodeIds.
  COMPARE_FIELD(fixedNodeIds);
  COMPARE_FIELD(freeNodeIds);

  // Compare doubles and ints.
  COMPARE_FIELD(a);
  COMPARE_FIELD(init_element_area);
  COMPARE_FIELD(rows);
  COMPARE_FIELD(cols);
  COMPARE_FIELD(load);
  COMPARE_FIELD(loadSteps);

  // Compare currentDeformation (Eigen::Matrix2d).
  // Eigen doesn't define operator== by default, so either define it or
  // approximate.
  COMPARE_FIELD(currentDeformation);

  COMPARE_FIELD(nrElements);
  COMPARE_FIELD(nrNodes);
  COMPARE_FIELD(totalEnergy);
  COMPARE_FIELD(averageEnergy);
  COMPARE_FIELD(maxEnergy);
  COMPARE_FIELD(maxForce);
  COMPARE_FIELD(averageP11);
  COMPARE_FIELD(averageP12);
  COMPARE_FIELD(averageP21);
  COMPARE_FIELD(averageP22);
  COMPARE_FIELD(averageSigma11);
  COMPARE_FIELD(averageSigma12);
  COMPARE_FIELD(averageSigma22);
  COMPARE_FIELD(averageSigmaTrace);
  COMPARE_FIELD(maxM3Nr);
  COMPARE_FIELD(sumM3Nr);
  COMPARE_FIELD(maxPlasticJump);
  COMPARE_FIELD(minPlasticJump);
  COMPARE_FIELD(redQuadrantCounts);
  COMPARE_FIELD(redQuadrantFixedCounts);
  COMPARE_FIELD(nr_elements_with_m3_fix_change);
  COMPARE_FIELD(nr_elements_with_m3_fix_changeInStep);
  COMPARE_FIELD(QDSD);
  COMPARE_FIELD(energyFunction);
  COMPARE_FIELD(bulkModulus);
  COMPARE_FIELD(usingPBC);
  // COMPARE_FIELD(nrMinFunctionCalls);// This one is not accurate over
  // save/load dumps
  // COMPARE_FIELD(nrMinItterations); // So we don't compare
  // this one too, just because

  // Compare strings.
  COMPARE_FIELD(simName);
  COMPARE_FIELD(dataPath);

  // Compare the raw array of doubles: bounds.
  COMPARE_FIELD(bounds);

  return equal;
}

/*
   Standard equality operator for Mesh.
   Declares compareconnectesInternal with a nullptr, so no debug messages are
   produced.
*/
inline bool operator==(const Mesh &lhs, const Mesh &rhs) {
  return compareconnectesInternal(lhs, rhs, nullptr);
}
inline bool operator!=(const Mesh &lhs, const Mesh &rhs) {
  return !(lhs == rhs);
}

inline std::string limitDebugLines(const std::string &text, size_t maxLines) {
  if (maxLines == 0 || text.empty()) {
    return text;
  }

  size_t lines = 0;
  size_t pos = 0;
  while (pos < text.size()) {
    size_t lineEnd = text.find('\n', pos);
    lines++;
    if (lines == maxLines) {
      if (lineEnd == std::string::npos) {
        return text;
      }
      if (lineEnd + 1 >= text.size()) {
        return text;
      }
      std::string out = text.substr(0, lineEnd);
      out += " ... (truncated after " + std::to_string(maxLines) + " lines)";
      return out;
    }
    if (lineEnd == std::string::npos) {
      return text;
    }
    pos = lineEnd + 1;
  }

  return text;
}

// Instead of just saying that the meshes are not the same, this function
// attempts to make it easy to see exactly what the difference is.
inline std::string debugCompare(const Mesh &lhs, const Mesh &rhs,
                                size_t maxLines = 100) {
  std::string diff;

  // Handle tElements and Nodes in a special way.
  if (lhs.elements.size() != rhs.elements.size()) {
    diff += "elements size differs; ";
  } else {
    // If sizes match, compare each element.
    for (size_t i = 0; i < lhs.elements.size(); i++) {
      if (!(lhs.elements[i] == rhs.elements[i])) {
        diff += "elements[" + std::to_string(i) + "] differs -> \n";
        // Recursively call debugCompare for TElement.
        diff += debugCompare(lhs.elements[i], rhs.elements[i], 1);
      }
    }
  }
  if (lhs.nodes.size() != rhs.nodes.size()) {
    diff += "nodes size differs; ";
  } else {
    // If sizes match, compare each element.
    for (size_t i = 0; i < lhs.nodes.size(); i++) {
      if (!(lhs.nodes(i) == rhs.nodes(i))) {
        diff += "nodes[" + std::to_string(i) + "] differs -> \n";
        // Recursively call debugCompare for TElement.
        diff += debugCompare(lhs.nodes(i), rhs.nodes(i), 1);
      }
    }
  }

  compareconnectesInternal(lhs, rhs, &diff);
  return limitDebugLines(diff, maxLines);
}

#endif
