#ifndef NODE_H
#define NODE_H
#include <vector>
#pragma once

#include <Eigen/Core>

#include <array>
#include <cereal/types/array.hpp> // Cereal serialization for std::vector

using namespace Eigen;

#define MAX_ELEMENTS_PER_NODE 6

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
    ar(i, row, col, cols);
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
  NodeId id;        // The identifier for this node.
  Vector2d f;       // The force experienced by the node.
  bool fixedNode;   // Flag indicating if the node is fixed or not.
  bool isGhostNode; // Flag indicating if it is only representing another node
                    // accross the periodoc boundary.
  NodeId ghostId; // This id points to the row, column and index of a n+1 by m+1
                  // system.
  Vector2d ghostShift; // This is the displacement from the normal position to
                       // the periodic

  // Fixed-size arrays for element information
  // Fixed-size array for element indices
  std::array<int, MAX_ELEMENTS_PER_NODE> elementIndices;
  // Fixed-size array for node indices
  std::array<int, MAX_ELEMENTS_PER_NODE> nodeIndexInElement;
  int elementCount = 0; // Tracks the current number of elements

private:
  // Whenever we update x/y or init x/y, we also need to update u x/y,
  // therefore, we need to make these private and access them through functions.
  Vector2d m_pos;
  Vector2d m_init_pos;
  Vector2d m_u;

public:
  std::array<NodeId, 4> neighbours; // Identifiers for the neighboring nodes.

  // Constructor to initialize the arrays with default values
  Node();

  // Constructor to initialize a Node with coordinates.
  Node(double x, double y);
  Node(double a, int row, int col, int cols);

  // Set the x and y variables
  void setPos(const Vector2d &pos);
  void addPos(const Vector2d &pos);

  // Set the initial x and y variables
  void setInitPos(const Vector2d &init_pos);

  // Set the pos using current initial pos and displacement
  void setDisplacement(const Vector2d &disp);

  // Add a force to the node
  void addForce(const Vector2d &f);

  // Set f_x and f_y to 0
  void resetForce();

  // Getters, making them read-only from outside.
  const Vector2d &pos() const { return m_pos; }
  const Vector2d &init_pos() const { return m_init_pos; }
  const Vector2d &u() const { return m_u; }

  friend std::ostream &operator<<(std::ostream &os, const Node &node);

  friend double tElementInitialArea(const Node &A, const Node &B,
                                    const Node &C);

private:
  // Function to update displacement based on the current and initial positions.
  void updateDisplacement();

  friend class cereal::access; // Necessary to serialize private members
  template <class Archive> void serialize(Archive &ar) {
    ar(id, f, fixedNode, isGhostNode, ghostId, ghostShift, elementIndices,
       nodeIndexInElement, neighbours, m_pos, m_init_pos, m_u);
  }
};

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
#endif