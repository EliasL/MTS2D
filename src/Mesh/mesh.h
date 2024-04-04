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
 * of triangular elements formed by neighboring nodes.
 *
 * NB: The nodes in the elements are seperate from the nodes in the mesh. This
 * means that you should modify the nodes in the mesh before creating the elements.
 * The elements will then create copies of the nodes with the modifications.
 */
class Mesh
{
public:
    // A matrix representing the mesh of nodes.
    Matrix<Node> nodes;

    // A collection of elements formed by triangles of nodes.
    // We use pointers so that the nodes in the elements can have a reference
    // to the mesh. Otherwise, we would need to set the reference at initialization.
    std::vector<TElement> elements;

    // IDs of nodes that are on the border of the mesh.
    std::vector<NodeId> fixedNodeIds;

    // IDs of nodes that are not on the border of the mesh.
    std::vector<NodeId> freeNodeIds;

    // The characteristic dimension of the mesh.
    double a;

    // Number of rows of nodes in the mesh.
    int rows;
    // Number of columns of nodes in the mesh.
    int cols;

    // The applied load on the mesh.
    // This variable is not used for physics. The physics are solely based on
    // the position of the boundary nodes. This value is stored for logging
    // purposes.
    double load;

    // If we want to shear the entire mesh, we use the periodic transform, but
    // if we want to create a change in the current load without moving any of
    // the nodes, we use the periodic load (which applies only to the "distances"
    // between the periodically repeated systems).

    // We need to know how to tile the system periodically. This transformation
    // is applied to the diplacement of the periodic nodes.
    Matrix2x2<double> currentDeformation;

    // The number of triangles created in the mesh.
    int nrElements;
    // Nr of nodes
    int nrNodes;

    // We calculate the total energy during the simulation, and the average
    // energy is useful to plot, so we keep this value here for easy access.
    double averageEnergy;
    // This might also be usefull
    double maxEnergy;

    // Used to make it seem like the ground state has an energy of 0
    double groundStateEnergy;

    // Flag for using periodic or fixed boundary conditions
    bool usingPBC;

    std::shared_ptr<Mesh> meshPtr;

    // Default constructor.
    Mesh();

    // Constructor to initialize the mesh with a specified number of rows, columns, and characteristic dimension.
    Mesh(int rows, int cols, double a, bool usingPBC = true);

    // Constructor to initialize the mesh with a specified number of rows and columns with the characteristic dimesion set to one.
    Mesh(int rows, int cols, bool usingPBC = true);

    // Overloaded indexing operator to access nodes by their NodeId.
    Node *operator[](NodeId id) { return &nodes.data[id.i]; }

    // Const overloaded indexing operator to access nodes by their NodeId.
    const Node *operator[](NodeId id) const { return &nodes.data[id.i]; }

    // Determines if a node is at the border of the mesh.
    bool isFixedNode(NodeId n_id);

    // Applies a transform to all nodes in the mesh, including the PBC.
    void applyTransformation(Matrix2x2<double> transformation);

    // Applies a transform to the border nodes.
    void applyTransformationToFixedNodes(Matrix2x2<double> transformation);

    // Applies a transform to the periodic boundary tranform.
    // (see how it affects the pos function in PeriodicNode )
    void applyTransformationToSystemDeformation(Matrix2x2<double> transformation);

    // Apply translation to all nodes in the mesh.
    void applyTranslation(Vector2d displacement);

    // This sets the current position as the initial position of the mesh.
    void setInitPos();

    // Resets the forces acting on all nodes in the mesh.
    void resetForceOnNodes();

    // Calculates averages
    double averageResolvedShearStress();

    // Checks for a change in the m matrixes of the elements and assumes that means plasticity TODO
    int nrPlasticEvents();

    // Fixes the border nodes in the mesh.
    void fixBorderNodes();

    // Fixes the nodes in a given row
    void fixNodesInRow(int row);

    // Fixes the nodes in a given column
    void fixNodesInColumn(int column);

    // Print element connectivity (for debugging)
    void printConnectivity(bool realId = true);

    // Creates triangles from neighboring nodes to form the elements of the mesh.
    void createElements();

    // Loops over all elements and updates them
    void updateElements();

    // Uses the ids in the elements to update the force on the nodes
    void applyForceFromElementsToNodes();

    // Calculates total and average energy. Returns the total.
    double calculateTotalEnergy();

    Vector2d makeGhostPos(Vector2d pos, Vector2d shift);

private:
    // Fills in the IDs of nodes that are not at the border.
    void m_updateFixedAndFreeNodeIds();

    // Sets the positions of nodes in the mesh based on mesh dimensions and spacing.
    void m_createNodes();

    // Fills the neighbor relationships between nodes in the mesh.
    void m_fillNeighbours();

    // Creates the NodeId of a node at a given position.
    NodeId m_makeNId(int row, int col);

    // Changes a node to a ghost node with a new row and column pos
    void m_makeGN(Node &n, int newRow, int newCol);

    // Retrives the NodeId of the neighbour of a node at a given position.
    Node m_getNeighbourNode(Node node, int direction);
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