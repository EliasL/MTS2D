#ifndef NODE_H
#define NODE_H
#pragma once

#include "settings.h"
#include "Matrix/matrix.h"
#include "Matrix/matrix2x2.h"
#include "spdlog/spdlog.h"
#include <array>
#include <vector>
#include <stdexcept>

// Declaration
struct PeriodicNode;

/**
 * @brief Identifier for a node.
 *
 * This structure holds the indices that uniquely identify a node's
 * position within a 2D surface, as well as a flattened index for a 1D representation.
 */
struct NodeId
{
    int row, col; // Mesh position indices in the x and y directions.
    int i;        // Flattened index for the node within a 1D array representation of the surface.

    // Default constructor.
    NodeId();

    // Constructor to initialize NodeId with x and y indices and total number of columns in the surface.
    NodeId(int row, int col, int cols);

    // Constructor to initialize NodeId with a flattened index and total number of columns in the surface.
    NodeId(int i, int cols);

    friend std::ostream &operator<<(std::ostream &os, const NodeId &nodeId);
};

/**
 * @brief Represents a node.
 *
 * Nodes are used to define the geometry of a surface and its physical properties,
 * such as forces applied at the node points.
 */
struct Node
{
    // Whenever we update x/y or init x/y, we also need to update u x/y,
    // therefore, we need to make these private and access them through functions
private:
    double m_x, m_y;           // Coordinates of the node in the surface.
    double m_init_x, m_init_y; // Coordinates of the initial position of the node.
    double m_u_x, m_u_y;       // Displacement from the initial to the current position
public:
    double f_x, f_y;                  // Force components acting on the node.
    bool fixedNode;                   // Flag indicating if the node is fixed or not
    NodeId id;                        // The identifier for this node.
    std::array<NodeId, 4> neighbours; // Identifiers for the neighboring nodes.

    // Default constructor.
    Node();

    // Constructor to initialize a Node with coordinates.
    Node(double x, double y);

    // Set the x and y variables
    void setPos(double x, double y);

    // Set the initial x and y variables
    void setInitPos(double x, double y);

    // Add a force to the node
    void addForce(std::array<double, 2> f);

    // Set f_x and f_y to 0
    void resetForce();

    // Copys the data from a periodic node
    void copyValues(PeriodicNode node);

    // Getters, making them read-only from outside.
    double x() const { return m_x; }
    double y() const { return m_y; }
    double init_x() const { return m_init_x; }
    double init_y() const { return m_init_y; }
    double u_x() const { return m_u_x; }
    double u_y() const { return m_u_y; }

private:
    // Function to update displacement based on the current and initial positions.
    void updateDisplacement();
};

/*
The mesh has many nodes, but when using periodic boundary conditions, some nodes
have to be in two places at once to avoid elements being drawn across the system.
In order to solve this issue, elements do not use nodes directly, but instead
through this wrapper class PeriodicNode. This way, the node can still be accessed
directly through references, but seemlessly appear in a different position.
*/
struct PeriodicNode
{
    double dx, dy; // Displacement of the periodic node from the real node.
    Node *realNode;
    // In the periocid mesh, there is an aditional row and column (to prevent
    // wrapping), so the indixes will be slightly different.
    NodeId periodicId; // The identifier for this node.

    // Default constructor.
    PeriodicNode();

    // Constructor to initialize a Node with coordinates.
    PeriodicNode(Node *realNode, double dx, double dy, int rows, int cols);

    // Add a force to the real node
    void addForce(std::array<double, 2> f);

    // Getters for the real node with modified position
    double x() const { return realNode->x() + dx; }
    double y() const { return realNode->y() + dy; }
    double init_x() const { return realNode->init_x() + dx; }
    double init_y() const { return realNode->init_y() + dy; }
    double u_x() const { return realNode->u_x(); }
    double u_y() const { return realNode->u_y(); }
    double f_x() const { return realNode->f_x; }
    double f_y() const { return realNode->f_y; }
};

// The neighbours should be indexed using these defines for added readability
#define LEFT_N 0
#define RIGHT_N 1
#define UP_N 2
#define DOWN_N 3

/**
 * @brief Transforms a node by applying a transformation matrix.
 *
 * This function applies a linear transformation defined by a matrix to the node's position.
 *
 * @param matrix The transformation matrix to apply.
 * @param n The node to transform.
 * @return The transformed node.
 */
Node transform(const Matrix2x2<double> &matrix, const Node &n);
void transformInPlace(const Matrix2x2<double> &matrix, Node &n);

/**
 * @brief Translates a node by a given displacement.
 *
 * This function adds a displacement to the node's position, with an optional multiplier
 * to scale the displacement.
 *
 * @param n The original node to be translated.
 * @param delta The displacement to apply to the node.
 * @param multiplier A scaling factor for the displacement (default is 1).
 * @return The translated node.
 */
Node translate(const Node &n, const Node &delta, double multiplier = 1);
void translateInPlace(Node &n, const Node &delta, double multiplier = 1);
void translateInPlace(Node &n, double x, double y, double multiplier = 1);

#endif