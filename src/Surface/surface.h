#ifndef SURFACE_H
#define SURFACE_H
#pragma once

#include "settings.h"
#include "Matrix/matrix.h"
#include "Matrix/matrix2x2.h"
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
    int xi, yi; // Mesh position indices in the x and y directions.
    int i;      // Flattened index for the node within a 1D array representation of the surface.
    
    // Default constructor.
    NodeId();
    
    // Constructor to initialize NodeId with x and y indices and total number of columns in the surface.
    NodeId(int xi, int yi, int cols);
    
    // Constructor to initialize NodeId with a flattened index and total number of columns in the surface.
    NodeId(int i, int cols);
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



/**
 * @brief Represents a triangular element defined by three nodes.
 *
 * This structure is used to model a discrete element within a mesh or surface,
 * where each corner of the triangle is a node.
 */
struct Triangle
{
    // Pointers to the nodes that form the vertices of the triangle.
    Node *a1;
    Node *a2;
    Node *a3;

    /**
     * @brief Calculates the first edge vector of the triangle.
     *
     * This vector represents one side of the triangle.
     * @return Returns the first edge vector as a 2D array.
     */
    std::array<double, 2> e1() const;

    /**
     * @brief Calculates the second edge vector of the triangle.
     *
     * Similar to e1, this vector represents another side of the triangle.
     * @return Returns the second edge vector as a 2D array.
     */
    std::array<double, 2> e2() const;

    /**
     * @brief Computes the metric tensor for the triangle.
     *
     * The metric tensor provides information about the geometric properties
     * of the triangle, such as lengths of sides, angles, and area.
     * @param f A function that defines how the metric is calculated.
     * @return Returns the metric tensor as a 2x2 matrix.
     */
    Matrix2x2<double> metric(MetricFunction f = MetricFunction::faicella) const;
};


/**
 * @brief Represents a cell in a material surface, characterized by its physical properties.
 *
 * A cell is formed by a triangle of nodes and contains information about the
 * deformation gradient F, the metric tensor C, the reduced metric tensor C_,
 * the reduction transformation matrix m, the reduces stress tensor r_s and the
 * Piola stress tensor P
 */
class Cell
{
public:
    // Deformation gradient / Basis vectors
    Matrix2x2<double> F;

    // Metric tensor (C = FF^T)
    Matrix2x2<double> C;

    // Reduced metric tensor
    Matrix2x2<double> C_;

    // Reduction transformation matrix (m^TCm = C_)
    Matrix2x2<double> m;

    // Reduced stress
    Matrix2x2<double> r_s;
    
    // TODO Is the comment below accurate? Or is it just normal stress, but calculated in a special way?
    // Second Piola-Kirchhoff stress tensor, representing the stress relative to the undeformed configuration.
    Matrix2x2<double> P;

    // TODO The comment below seems plausible, but I am not quite sure. 
    // Strain energy of the cell, representing the potential energy stored due to deformation.
    double energy;

    // UNUSED TODO When implemented, rewrite comment
    // Flag indicating if the cell can undergo plastic (permanent) deformation.
    bool plasticity;

    /**
     * @brief Initializes Cell and calculates several values:
     * 
     *  the deformation gradient D, 
     *  the metric tension C, 
     *  the transformation matrix m, 
     *  and the reduced metric tension C_.
     * 
     */
    Cell(Triangle t);

    // Default constructor for the Cell.
    Cell();

    // Accesses the first basis vector of the cell (F[0])
    double e1(int index);

    // Accesses the second basis vector of the cell (F[1])
    double e2(int index);

    // Sets the forces on the nodes that form the cell's triangle.
    void setForcesOnNodes(Triangle t);

private:
    // Computes the deformation gradient for the cell based on the triangle's vertices.
    void m_getDeformationGradiant(Triangle t);

    // Performs a Lagrange reduction on C to calculate C_.
    void m_lagrangeReduction();
};


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
    Mesh(int n, int m, double a);

    // Constructor to initialize the surface with a specified number of rows and columns with the characteristic dimesion set to one.
    Mesh(int n, int m);

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
    NodeId getNodeId(int xi, int yi);

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


#endif