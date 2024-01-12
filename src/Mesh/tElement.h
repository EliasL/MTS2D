#ifndef TELEMENT_H
#define TELEMENT_H
#pragma once

#include "settings.h"
#include "Matrix/matrix.h"
#include "Matrix/matrix2x2.h"
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
 *  The shape functions used are
 *  N1 = 1 - ξ1 - ξ2  
 *  N2 = ξ1  
 *  N3 = ξ2  
 */
class TElement
{
public:
    // Pointers to the nodes that form the vertices of the element.
    Node *n1;
    Node *n2;
    Node *n3;

    // Deformation gradient
    Matrix2x2<double> F;

    // Metric tensor (C = FF^T)
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

    // TODO The comment below seems plausible, but I am not quite sure.
    // Strain energy of the cell, representing the potential energy stored due
    // to deformation.
    double energy = 0;

private:
    /*
    Shape functions:
    N1 = 1 - ξ1 - ξ2
    N2 = ξ1
    N3 = ξ2

    Derivative of shape functions:
    For simplicity of implementation, these derivatives are assumed to be
    constant! If you change b1, b2 or b3, you will also need to manually change
    the implementation of the jacobian calculation. See m_calculateJacobian.
    */
    std::array<double, 2> b1 = {-1, -1}; // ∂N1/∂ξi (i=1,2)
    std::array<double, 2> b2 = {1, 0};   // ∂N2/∂ξi (i=1,2)
    std::array<double, 2> b3 = {0, 1};   // ∂N3/∂ξi (i=1,2)

    // The jacobian is given by ∂xi/∂ξk where xi is a sum of the three shape
    // functions N1, N2 and N3.
    Matrix2x2<double> jacobian;

    // Inverse of the jacobian, but only updated once at the beginning of the
    // simulation, hence a Ref(rence)
    Matrix2x2<double> invJacobianRef;

    // These are ajustment vectors that we multiply together with the piola
    // tensor to correctly extract the force coresponding to each node. 
    // Similarly to invJacobianRef, these only update once, during initialization.
    std::array<double, 2> r1;
    std::array<double, 2> r2;
    std::array<double, 2> r3;

    // Various numbers used in energy and reduced stress calculation. TODO understand and comment
    double burgers = 1.;
    double beta = -0.25;
    double K = 4.;

public:
    // Constructor for the triangular element. Initializes the 3 defining nodes
    // and calculates the inverse state A_inv, to later be used in calculating F.
    TElement(Node *n1, Node *n2, Node *n3);
    TElement();

    // The vector from node 1 to node 2
    std::array<double, 2> e12() const;
    // The vector from node 1 to node 3
    std::array<double, 2> e13() const;
    // The vector from node 2 to node 3
    std::array<double, 2> e23() const;

    /**
     * @brief Initializes TElement and calculates several values:
     *
     *  the deformation gradient D,
     *  the metric tension C,
     *  the transformation matrix m,
     *  and the reduced metric tension C_.
     *
     */
    void update();

    // Sets the forces on the nodes that form the cell's triangle.
    void applyForcesOnNodes();

    // Usefull if you only care about the energy given the C matrix.
    static double calculateEnergy(double c11, double c22, double c12);

    // Used for testing the lagrange reuction functions
    static TElement lagrangeReduction(double c11, double c22, double c12);

private:
    // updates current state using two vectors from the triangle
    void m_calculateJacobian();

    // Computes the deformation gradient for the cell based on the triangle's vertices.
    void m_updateDeformationGradiant();

    // Computes the metric tensor for the triangle.
    void m_updateMetricTensor();

    // Performs a Lagrange reduction on C to calculate C_.
    void m_lagrangeReduction();

    // Under certain conditions, the normal LR is very slow and can be simplified
    // This funcion usually calls the normal LR, but uses a faster method when
    // possible.
    void m_fastLagrangeReduction();

    // Calculates energy
    void m_calculateEnergy();

    // Calculate reduced stress
    void m_calculateReducedStress();

    // Calculate Piola stress P
    void m_updatePiolaStress();
};

std::ostream &operator<<(std::ostream &os, const TElement &element);

#endif