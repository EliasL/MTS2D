#ifndef TELEMENT_H
#define TELEMENT_H
#pragma once

#include "settings.h"
#include "Matrix/matrix.h"
#include "Matrix/matrix2x2.h"
#include "Simulation/energyFunctions.h"
#include "node.h"
#include "spdlog/spdlog.h"
#include <array>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip> // Include this for std::fixed and std::setprecision

/**
 * @brief Represents a triangular element in a material surface, characterized
 * by its physical properties.
 *
 * A triangular element is formed by a triangle of nodes and contains
 * information about the deformation gradient F,
 * the metric tensor C, the reduced metric tensor C_, the reduction
 * transformation matrix m, the reduces stress tensor r_s and the
 * Piola stress tensor P.
 *
 *  The shape functions used are:
 *  N1 = 1 - ξ1 - ξ2
 *  N2 = ξ1
 *  N3 = ξ2
 *
 *  u referes to displacement.
 *  X referes to the reference state.
 *  x refers to the current state.
 *
 */
class TElement
{
public:
    // Periodic nodes that form the vertices of the element.
    // Think of periodic nodes as normal nodes. The only differece is that two
    // periodic nodes might be getting their values from the same real node.
    std::array<PeriodicNode, 3> nodes;

    // Deformation gradient
    Matrix2x2<double> F;

    // Metric tensor (C = F^TF)
    Matrix2x2<double> C;

    // Reduced metric tensor
    Matrix2x2<double> C_;

    // Reduction transformation matrix (m^TCm = C_)
    Matrix2x2<double> m;

    // Reduced stress
    Matrix2x2<double> r_s;

    // TODO Is the comment below accurate? Or is it just normal stress, but
    // calculated in a special way?
    // First Piola-Kirchhoff stress tensor, representing the stress relative
    // to the undeformed configuration.
    Matrix2x2<double> P;

    // Strain energy of the cell, representing the potential energy stored due
    // to deformation.
    double energy = 0;

    // A representation of stress that is unaffected by the directionality of
    // loading. Discontinuous yielding of pristine micro-crystals (page 216/17)
    double resolvedShearStress = 0;

    // The jacobian of the shapefunction with respect to the reference state.
    //  We use this to calculate the deformation gradiant F
    // ∂ξ/∂X
    Matrix2x2<double> dxi_dX;

    // This is the row, column position of where the element was initialized.
    // The exact relationship between this possition and the nodes of the element
    // depends on whether it is a top or bottom triangle.
    double row;
    double col;

private:
    /*
    Shape functions:
    N1 = 1 - ξ1 - ξ2
    N2 = ξ1
    N3 = ξ2

    Derivative of shape functions:
    For simplicity of implementation, these derivatives are assumed to be
    constant! If you change b1, b2 or b3, you will also need to manually change
    the implementation of the jacobian calculation. See du_dxi.
    */
    static constexpr std::array<std::array<double, 2>, 3> b = {{{-1, -1},
                                                                {1, 0},
                                                                {0, 1}}};

    // These are adjustment vectors that we multiply together with the piola
    // tensor to correctly extract the force corresponding to each node.
    // Similarly to dxi_dX, these only update once, during initialization.
    std::array<VArray, 3> r;

    // We only save data when plasticity occurs, so we keep a reference of
    // how many times m3 is applied in the lagrange reduction. If this number
    // changes from one reduction to another, we know that a plastic event has
    // occured. (ie. the energy potential suddenly has a gradient in a new
    // direction, ie. the node has fallen into a different energy well.)
    int m3Nr = 0;
    int past_m3Nr = 0;

    // Various numbers used in energy and reduced stress calculation. TODO understand and comment
    // Coresponds (somehow) to square lattice. beta=4 gives triangular lattice.
    static constexpr double beta = -0.25;
    // Bulk modulus. Controlls the contribution of the volumetric energy function. (or something)
    // called K in Umut's code
    static constexpr double mu = 4.;

public:
    // Constructor for the triangular element. Initializes the 3 defining nodes
    // and calculates the inverse state A_inv, to later be used in calculating F.
    TElement(Mesh &mesh, PeriodicNode n1, PeriodicNode n2, PeriodicNode n3, int row, int col);
    TElement(){};

    /**
     * @brief Initializes TElement and calculates several values:
     *
     *  the deformation gradient D,
     *  the metric tension C,
     *  the transformation matrix m,
     *  and the reduced metric tension C_.
     *
     */
    void update(Mesh &mesh);

    // Sets the forces on the nodes that form the cell's triangle.
    void applyForcesOnNodes(Mesh &mesh);

    // // Usefull if you only care about the energy given the C matrix.
    static double calculateEnergy(double c11, double c22, double c12);

    // // Used for testing the lagrange reuction functions
    static TElement lagrangeReduction(double c11, double c22, double c12);

    // Check if m had changed NB Can only be called once per frame!
    bool plasticEvent();

    void copyValues(const TElement &other);

private:
    // Calculate the Jacobian with respect to the displacement of the nodes
    Matrix2x2<double> du_dxi(Mesh &mesh);
    // Calculate the Jacobian with respect to the initial position of the nodes
    Matrix2x2<double> dX_dxi(Mesh &mesh);

    // Computes the deformation gradient for the cell based on the triangle's vertices.
    void m_updateDeformationGradiant(Mesh &mesh);

    // Computes the metric tensor for the triangle.
    void m_updateMetricTensor();

    // Performs a Lagrange reduction on C to calculate C_.
    void m_lagrangeReduction();

    // Calculates energy
    void m_updateEnergy();

    // Calculate reduced stress
    void m_updateReducedStress();

    // Calculate Piola stress P
    void m_updatePiolaStress();

    // Calculate the resolved-shear stress
    void m_updateResolvedShearStress();

    // Calculates the difference in displacement between two nodes
    VArray du(Mesh &mesh, int idx1, int idx2) const;

    // Calculates the difference in position between two nodes
    VArray dx(Mesh &mesh, int idx1, int idx2) const;

    // Calculates the difference in initial position between two nodes
    VArray dX(Mesh &mesh, int idx1, int idx2) const;
};

std::ostream &operator<<(std::ostream &os, const TElement &element);

#endif