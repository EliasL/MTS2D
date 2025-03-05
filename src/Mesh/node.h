#ifndef NODE_H
#define NODE_H
#include "Eigen/Core"
#include <ostream>
#pragma once

#include "Data/cereal_help.h"
#include "compare_macros.h"
#include <Eigen/Core>
#include <array>
#include <cereal/types/array.hpp> // Cereal serialization for std::vector
#include <string>
using namespace Eigen;

#define MAX_ELEMENTS_PER_NODE 8

/**
 * @brief Identifier for a node.
 *
 * This structure holds the indices that uniquely identify a node's
 * position within a 2D surface, as well as a flattened index for a 1D
 * representation.
 */
struct NodeId {
  int i; // Flattened index for the node within a 1D array representation of the
         // surface.
  int row, col; // Mesh position indices in the x and y directions.
  int cols; // Since we take it as an argument, we might as well save it. (We
            // use it in the PeriodicNode constructor)
  // Default constructor.
  NodeId();

  // Constructor to initialize NodeId with x and y indices and total number of
  // columns in the surface.
  NodeId(int row, int col, int cols);

  // Constructor to initialize NodeId with a flattened index and total number of
  // columns in the surface.
  NodeId(int i, int cols);

  friend std::ostream &operator<<(std::ostream &os, const NodeId &nodeId);
  template <class Archive> void serialize(Archive &ar) {
    ar(MAKE_NVP(i), MAKE_NVP(row), MAKE_NVP(col), MAKE_NVP(cols));
  }
};

/**
 * @brief Represents a node.
 *
 * Nodes are used to define the geometry of a surface and its physical
 * properties, such as forces applied at the node points.
 *
 * We make pos and init_pos be private variables to avoid the user
 * forgetting to update the dispacement.
 */
struct Node {
public:
  NodeId id;      // The identifier for this node.
  Vector2d f;     // The force experienced by the node.
  bool fixedNode; // Flag indicating if the node is fixed or not.
                  // accross the periodoc boundary.

  // Fixed-size arrays for element information
  // Fixed-size array for element indices
  std::array<int, MAX_ELEMENTS_PER_NODE> elementIndices;
  // Fixed-size array for node indices
  std::array<int, MAX_ELEMENTS_PER_NODE> nodeIndexInElement;
  // Identifiers for the neighboring nodes in the reference state.
  std::array<NodeId, 4> refNeighbours;

  int elementCount = 0; // Tracks the current number of elements

private:
  // Whenever we update x/y or init x/y, we also need to update u x/y,
  // therefore, we need to make these private and access them through functions.
  Vector2d m_pos;      // Current state x
  Vector2d m_init_pos; // Reference state X
  Vector2d m_u;        // Displacement u

public:
  // Constructor to initialize the arrays with default values
  Node();

  // Constructor to initialize a Node with coordinates.
  Node(double x, double y);
  Node(int row, int col, int cols);
  Node(double a, int row, int col, int cols);

  // Set the x and y variables
  void setPos(const Vector2d &pos);
  void addPos(const Vector2d &pos);

  // Set the initial x and y variables
  void setInitPos(const Vector2d &init_pos);

  // Set the pos using current initial pos and displacement
  void setDisplacement(const Vector2d &disp);

  // Adds a vector to the displacement of the node
  void addDisplacement(const Vector2d &dispChange);

  // Add a force to the node
  void addForce(const Vector2d &f);

  // Set f_x and f_y to 0
  void resetForce();

  void applyDeformation(const Matrix2d &deformation);

  // Getters, making them read-only from outside.
  const Vector2d &pos() const { return m_pos; }
  const Vector2d &init_pos() const { return m_init_pos; }
  const Vector2d &u() const { return m_u; }

  friend std::ostream &operator<<(std::ostream &os, const Node &node);

  friend bool compareNodesInternal(const Node &lhs, const Node &rhs,
                                   std::string *debugMsg, int tabNumber);

private:
  // Function to update displacement based on the current and initial positions.
  void updateDisplacement();

  friend class cereal::access; // Necessary to serialize private members
  template <class Archive> void serialize(Archive &ar) {
    ar(MAKE_NVP(id), MAKE_NVP(f), MAKE_NVP(fixedNode), MAKE_NVP(elementIndices),
       MAKE_NVP(nodeIndexInElement), MAKE_NVP(refNeighbours),
       MAKE_NVP(elementCount), MAKE_NVP(m_pos), MAKE_NVP(m_init_pos),
       MAKE_NVP(m_u));

    // LOAD_WITH_DEFAULT(ar, elementCount, MAX_ELEMENTS_PER_NODE);
  }
};

// Elements will not be given "real" nodes, they will only use these ghost
// nodes, which are imperfect copies of reference nodes in the matrix of "real"
// nodes in the mesh. By using GhostNodes, we can create several different
// copies of the real node, but shifted in different directions as needed to
// construct the periodic boundary.

struct GhostNode {
public:
  NodeId referenceId; // The identifier for the real node.
  NodeId ghostId;     // This id points to the row, column and index of a n+1
                      // by m+1 system.

  Vector2d f; // The force experienced by this node.

  Vector2d periodShift; // Shift from reference pos across the system
  Vector2d pos;         // Current state x
  Vector2d init_pos;    // Reference state X
  Vector2d u;           // Displacement u

  GhostNode(const Node *referenceNode, int row, int col, int cols, double a,
            const Matrix2d &currentDeformation);

  GhostNode(const Node *referenceNode, double a, const Matrix2d &deformation);

  GhostNode(const Node *referenceNode, int row, int col, int cols,
            const Matrix2d &currentDeformation);

  GhostNode(const Node *referenceNode, const Matrix2d &deformation);

  GhostNode() = default;

  void updatePosition(const Node *referenceNode,
                      const Matrix2d &currentDeformation);

  template <class Archive> void serialize(Archive &ar) {
    ar(MAKE_NVP(referenceId), MAKE_NVP(ghostId), MAKE_NVP(f),
       MAKE_NVP(periodShift), MAKE_NVP(pos), MAKE_NVP(init_pos), MAKE_NVP(u));
  }
};

std::ostream &operator<<(std::ostream &os, const GhostNode &node);

// The neighbours should be indexed using these defines for added readability
#define LEFT_N 0
#define RIGHT_N 1
#define UP_N 2
#define DOWN_N 3

/**
 * @brief Transforms a node by applying a transformation matrix.
 *
 * This function applies a linear transformation defined by a matrix to the
 * node's position.
 *
 * @param matrix The transformation matrix to apply.
 * @param n The node to transform.
 * @return The transformed node.
 */
Node transform(const Matrix2d &matrix, const Node &n);
void transformInPlace(const Matrix2d &matrix, Node &n);

/**
 * @brief Translates a node by a given displacement.
 *
 * This function adds a displacement to the node's position, with an optional
 * multiplier to scale the displacement.
 *
 * @param n The original node to be translated.
 * @param delta The displacement to apply to the node.
 * @param multiplier A scaling factor for the displacement (default is 1).
 * @return The translated node.
 */
Node translate(const Node &n, const Node &delta, double multiplier = 1);
void translateInPlace(Node &n, const Vector2d &disp, double multiplier = 1);
void translateInPlace(Node &n, double x, double y, double multiplier = 1);
void translateInPlace(Node &n, const Node &delta, double multiplier = 1);

// Overload the << operator for Vector2d
std::ostream &operator<<(std::ostream &os, const Vector2d &arr);

// ********************************************************************
// Overload the equality operator (==) for the NodeId and Node classes.
// ********************************************************************

inline bool compareNodeIdsInternal(const NodeId &lhs, const NodeId &rhs,
                                   std::string *debugMsg = nullptr,
                                   int tabNumber = 0) {
  bool equal = true;
  COMPARE_FIELD(i);
  COMPARE_FIELD(row);
  COMPARE_FIELD(col);
  COMPARE_FIELD(cols);
  return equal;
}
inline bool operator==(const NodeId &lhs, const NodeId &rhs) {
  return compareNodeIdsInternal(lhs, rhs, nullptr);
}
inline bool operator!=(const NodeId &lhs, const NodeId &rhs) {
  return !(lhs == rhs);
}
/*
   Internal helper function that compares two Node objects field by field.
   It takes an optional std::string pointer (debugMsg) that, if not null,
   collects messages for any differences found.
   When debugMsg is nullptr, it simply returns whether the objects are equal.
*/
inline bool compareNodesInternal(const Node &lhs, const Node &rhs,
                                 std::string *debugMsg = nullptr,
                                 int tabNumber = 0) {
  bool equal = true;

  // Compare public member variables.
  COMPARE_FIELD(id);
  COMPARE_FIELD(f);
  COMPARE_FIELD(fixedNode);
  COMPARE_FIELD(elementIndices);
  COMPARE_FIELD(nodeIndexInElement);
  COMPARE_FIELD(elementCount);
  COMPARE_FIELD(refNeighbours);

  // Compare private member variables directly since compareNodesInternal is a
  // friend.
  COMPARE_FIELD(m_pos);
  COMPARE_FIELD(m_init_pos);
  COMPARE_FIELD(m_u);

  return equal;
}

/*
   Internal helper function that compares two GhostNode objects field by field.
   It takes an optional std::string pointer (debugMsg) that, if not null,
   collects messages for any differences found.
   When debugMsg is nullptr, it simply returns whether the objects are equal.
*/
inline bool compareNodesInternal(const GhostNode &lhs, const GhostNode &rhs,
                                 std::string *debugMsg = nullptr,
                                 int tabNumber = 0) {
  bool equal = true;
  // Compare public member variables.
  COMPARE_FIELD(referenceId);
  COMPARE_FIELD(ghostId);
  COMPARE_FIELD(f);
  COMPARE_FIELD(periodShift);
  COMPARE_FIELD(pos);
  COMPARE_FIELD(init_pos);
  COMPARE_FIELD(u);
  return equal;
}

/*
   Standard equality operator for Node.
   This function is declared as a friend inside Node so that it can access
   private members. It calls compareNodesInternal with a nullptr to avoid
   generating debug messages.
*/
inline bool operator==(const Node &lhs, const Node &rhs) {
  return compareNodesInternal(lhs, rhs, nullptr);
}
inline bool operator!=(const Node &lhs, const Node &rhs) {
  return !(lhs == rhs);
}

/*
   Standard equality operator for GhostNode.
*/
inline bool operator==(const GhostNode &lhs, const GhostNode &rhs) {
  return compareNodesInternal(lhs, rhs, nullptr);
}
inline bool operator!=(const GhostNode &lhs, const GhostNode &rhs) {
  return !(lhs == rhs);
}

/*
   Debug function that uses the same internal comparison logic.
   It returns a string describing which fields differ between the two Node
   objects.

  Example use:
  std::string diff = debugCompare(node1, node2);
  if (!diff.empty()) {
      std::cerr << "Node differences: " << diff << std::endl;
  }

*/
inline std::string debugCompare(const Node &lhs, const Node &rhs,
                                int tabNumber = 0) {
  std::string diff;
  compareNodesInternal(lhs, rhs, &diff, tabNumber);
  return diff;
}
inline std::string debugCompare(const GhostNode &lhs, const GhostNode &rhs,
                                int tabNumber = 0) {
  std::string diff;
  compareNodesInternal(lhs, rhs, &diff, tabNumber);
  return diff;
}

#endif // NODE_H