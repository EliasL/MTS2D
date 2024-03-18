#ifndef NODE_H
#define NODE_H
#pragma once

#include "settings.h"
#include "Matrix/matrix.h"
#include "Matrix/matrix2x2.h"
#include "spdlog/spdlog.h"
#include <array>
#include <valarray>
#include <vector>
#include <stdexcept>

// Declarations
class Mesh;
struct PeriodicNode;

// Name simplification
using VArray = std::valarray<double>;

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
    int cols;     // Since we take it as an argument, we might as well save it. (We use it in the PeriodicNode constructor)
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
    VArray m_pos;
    VArray m_init_pos;
    VArray m_u;

public:
    VArray f;
    bool fixedNode;                   // Flag indicating if the node is fixed or not
    NodeId id;                        // The identifier for this node.
    std::array<NodeId, 4> neighbours; // Identifiers for the neighboring nodes.

    // Default constructor.
    Node();

    // Constructor to initialize a Node with coordinates.
    Node(double x, double y);

    // Set the x and y variables
    void setPos(VArray pos);

    // Set the initial x and y variables
    void setInitPos(VArray init_pos);

    // Set the pos using current initial pos and displacement
    void setDisplacement(VArray disp);

    // Add a force to the node
    void addForce(VArray f);

    // Set f_x and f_y to 0
    void resetForce();

    // Copys the data from a periodic node
    void copyValues(Mesh &mesh, PeriodicNode node);

    // Getters, making them read-only from outside.
    double x() const { return m_pos[0]; }
    double y() const { return m_pos[1]; }
    double init_x() const { return m_init_pos[0]; }
    double init_y() const { return m_init_pos[1]; }
    double u_x() const { return m_u[0]; }
    double u_y() const { return m_u[1]; }
    double f_x() const { return f[0]; }
    double f_y() const { return f[1]; }
    VArray pos() const { return m_pos; }
    VArray init_pos() const { return m_init_pos; }
    VArray u() const { return m_u; }

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
    // This is the translation required for the node to be placed on the other
    // side of the system
    VArray periodicShift;
    // Instead of using poiner to the real node, we use the id and a shared
    // pointer to the mesh. (The pointer to the mesh is needed anyway)
    NodeId realId;
    // In the non-periocid mesh, there is an aditional row and column (to prevent
    // wrapping), so the indixes will be slightly different.
    NodeId periodicId;

    // This bool tells you if it is a periodic node or a "normal" node
    bool isPeriodic;

    PeriodicNode(NodeId nodeId);
    PeriodicNode(){};

    // Add a force to the real node
    void addForce(Mesh &mesh, VArray f);

    // Getters for the real node with modified position
    VArray pos(Mesh &mesh) const;
    VArray init_pos(Mesh &mesh) const;
    VArray u(Mesh &mesh) const;
    VArray f(Mesh &mesh) const;

    // Updates periodic shift and periodicId
    void updatePeriodicity(Mesh &mesh, bool shiftX, bool shiftY);
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
void translateInPlace(Node &n, VArray disp, double multiplier = 1);
void translateInPlace(Node &n, double x, double y, double multiplier = 1);
void translateInPlace(Node &n, const Node &delta, double multiplier = 1);

#endif