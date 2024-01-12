#ifndef MESH_H
#define MESH_H
#pragma once

#include "settings.h"
#include "Matrix/matrix.h"
#include "Matrix/matrix2x2.h"
#include "node.h"
#include "tElement.h"
#include "spdlog/spdlog.h"
#include <array>
#include <vector>
#include <stdexcept>

/**
 * @brief The Mesh class represents a 2D surface mesh of nodes.
 *
 * It tracks nodes at the surface's border, neighbors, and manages the creation
 * of triangular elements formed by neighboring nodes. These triangles are used
 * to represent elements with specific properties within the surface structure.
 */
class Mesh
{
public:
    // A matrix representing the surface of nodes.
    Matrix<Node> nodes;

    // A collection of elements formed by triangles of nodes.
    std::vector<TElement> elements;

    // IDs of nodes that are on the border of the surface.
    std::vector<NodeId> borderNodeIds;

    // IDs of nodes that are not on the border of the surface.
    std::vector<NodeId> nonBorderNodeIds;

    // The characteristic dimension of the surface.
    double a;

    // The applied load on the surface.
    // This variable is not used for physics. The physics are solely based on
    // the position of the boundary nodes. This value is stored for logging
    // purposes. 
    double load;

    // The number of triangles created in the surface.
    int nrElements;

    // Default constructor.
    Mesh();

    // Constructor to initialize the surface with a specified number of rows, columns, and characteristic dimension.
    Mesh(int rows, int cols, double a);

    // Constructor to initialize the surface with a specified number of rows and columns with the characteristic dimesion set to one.
    Mesh(int rows, int cols);

    // Overloaded indexing operator to access nodes by their NodeId.
    Node *operator[](NodeId id) { return &nodes.data[id.i]; }

    // Const overloaded indexing operator to access nodes by their NodeId.
    const Node *operator[](NodeId id) const { return &nodes.data[id.i]; }

    // Determines if a node is at the border of the surface.
    bool isBorder(NodeId n_id);

    // Applies a transform to all nodes in the mesh.
    void applyTransformation(Matrix2x2<double> transformation);

    // Applies a transform to the border nodes.
    void applyTransformationToBoundary(Matrix2x2<double> transformation);
    
    // Resets the forces acting on all nodes in the surface.
    void resetForceOnNodes();

    // Retrieves the NodeId for a node at a given surface position.
    NodeId getNodeId(int row, int col);



private:
    // Identifies and marks the border elements in the surface.
    void m_setBorderElements();

    // Fills in the IDs of nodes that are not at the border.
    void m_fillNonBorderNodeIds();

    // Sets the positions of nodes in the surface based on surface dimensions and spacing.
    void m_setNodePositions();

    // Fills the neighbor relationships between nodes in the surface.
    void m_fillNeighbours();

    // Creates triangles from neighboring nodes to form the elements of the surface.
    void m_createElements();

};

std::ostream &operator<<(std::ostream &os, const Mesh &mesh);

/**
 * @brief Transforms all nodes in a mesh by applying a transformation matrix.
 *
 * This function applies a linear transformation defined by a matrix to the node's position.
 * 
 * @param matrix The transformation matrix to apply.
 * @param mesh The mesh to transform.
 */
void transform(const Matrix2x2<double> &matrix, Mesh &mesh);
// Only transform nodes in the provided list
void transform(const Matrix2x2<double> &matrix, Mesh &mesh, std::vector<NodeId> nodesToTransform);

/**
 * @brief Translates all nodes in a mesh by a given displacement.
 *
 * This function adds a displacement to the node's position, with an optional multiplier
 * to scale the displacement.
 * 
 * @param mesh The original node to be translated.
 * @param x The displacement in the x direction.
 * @param y The displacement in the y direction.
 */
void translate(Mesh &mesh, double x, double y);
// Only translate nodes in the provided list
void translate(Mesh &mesh, std::vector<NodeId> nodesToTranslate, double x, double y);


#endif