#ifndef MESH_H
#define MESH_H
#pragma once

#include "settings.h"
#include "Matrix/matrix.h"
#include "Matrix/matrix2x2.h"
#include "node.h"
#include "triangle.h"
#include "cell.h"
#include "easylogging++.h"
#include <array>
#include <vector>
#include <stdexcept>





/**
 * @brief The Mesh class represents a 2D surface mesh of nodes.
 *
 * It tracks nodes at the surface's border, neighbors, and manages the creation
 * of triangular elements formed by neighboring nodes. These triangles are used
 * to represent cells with specific properties within the surface structure.
 */
class Mesh
{
public:
    // A matrix representing the surface of nodes.
    Matrix<Node> nodes;

    // A collection of triangular elements within the surface.
    std::vector<Triangle> triangles;

    // A collection of cells formed by triangles of nodes.
    std::vector<Cell> cells;

    // IDs of nodes that are on the border of the surface.
    std::vector<NodeId> borderNodeIds;

    // IDs of nodes that are not on the border of the surface.
    std::vector<NodeId> nonBorderNodeIds;

    // The characteristic dimension of the surface.
    double a;

    // The applied load on the surface.
    double load;

    // The number of triangles created in the surface.
    int nrTriangles;

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

    // Applies boundary conditions to the surface.
    void applyBoundaryConditions(BoundaryConditions bc);

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

    // Creates triangles from neighboring nodes to form the cells of the surface.
    void m_createTriangles();

};

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


//TODO This whole way of dealing with boundary conditions is overengineered. It should be simpler.
// This should be redone when you have a bit more exmerience using the BC and after talking with Umut.

/**
 * @brief Represents the boundary conditions for a Mesh.
 */
class BoundaryConditions
{
public:
    // The load applied to the boundary nodes.
    double load;

    // The angle of the applied load.
    double theta;

    // The deformation gradient F that is applied to the boundary nodes.
    Matrix2x2<double> F;

    /**
     * @brief Constructor to set up the boundary conditions with a load and an angle.
     * 
     * @param load The magnitude of the load applied to the boundary.
     * @param theta The angle at which the load is applied.
     */
    BoundaryConditions(double load, double theta, BCF bc=BCF::macroShear);

private:
    // Calculates the gradient of a field across the boundary.
    void m_calculateGradiant(BCF bc);

    /**
     * @brief Applies a macroscopic shear deformation to the boundary.
     * 
     * This function defines how the boundary will deform under a shearing
     * load, which is important for simulations involving material strain.
     */
    void m_macroShear();
};


#endif