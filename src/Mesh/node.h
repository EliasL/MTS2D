#ifndef NODE_H
#define NODE_H
#pragma once

#include "settings.h"
#include "Matrix/matrix.h"
#include "Matrix/matrix2x2.h"
#include "easylogging++.h"
#include <array>
#include <vector>
#include <stdexcept>


/**
 * @brief Identifier for a node.
 *
 * This structure holds the indices that uniquely identify a node's
 * position within a 2D surface, as well as a flattened index for a 1D representation.
 */
struct NodeId
{
    int row, col; // Mesh position indices in the x and y directions.
    int i;      // Flattened index for the node within a 1D array representation of the surface.
    
    // Default constructor.
    NodeId();
    
    // Constructor to initialize NodeId with x and y indices and total number of columns in the surface.
    NodeId(int row, int col, int cols);
    
    // Constructor to initialize NodeId with a flattened index and total number of columns in the surface.
    NodeId(int i, int cols);

    friend std::ostream& operator<<(std::ostream &os, const NodeId &nodeId);
};

/**
 * @brief Represents a node.
 *
 * Nodes are used to define the geometry of a surface and its physical properties,
 * such as forces applied at the node points.
 */
struct Node
{
    double x, y;       // Coordinates of the node in the surface.
    double f_x, f_y;   // Force components acting on the node.
    bool borderNode;   // Flag indicating if the node is at the border of the surface.
    NodeId id;         // The identifier for this node.
    std::array<NodeId, 4> neighbours; // Identifiers for the neighboring nodes.

    // Default constructor.
    Node();
    
    // Constructor to initialize a Node with coordinates.
    Node(double x, double y);
};

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
void translateInPlace(Node &n, double x, double y);

#endif