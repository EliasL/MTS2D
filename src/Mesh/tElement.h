#ifndef TELEMENT_H
#define TELEMENT_H
#pragma once

#include "settings.h"
#include "Matrix/matrix.h"
#include "Matrix/matrix2x2.h"
#include "node.h"
#include "easylogging++.h"
#include <array>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip> // Include this for std::fixed and std::setprecision

/**
 * How do we calculate F? First, we need a representation of our initial state.
 * To do this, we take any two vectors in the triangular element. 
 * We use this reference state of the element to calculate the inverse of the
 * reference state: invRefState. Using this together with the current state, 
 * we can extract any deformations that have occured since initializing the 
 * element. Example: A is reference state, C is the current state, S is deformation:
 *      In general, we have that:
 *          S=CA^-1
 *      If there have been no changes made to the element, then C=A, so
 *          S=AA^-1=I
 *      as expected.
 */


/**
 * @brief Represents a triangular element in a material surface, characterized
 * by its physical properties.
 *
 * A triangular element is formed by a triangle of nodes and contains
 * information about the inverse of the reference state, deformation gradient F,
 * the metric tensor C, the reduced metric tensor C_, the reduction
 * transformation matrix m, the reduces stress tensor r_s and the
 * Piola stress tensor P
 */
class TElement
{
public:
    // Pointers to the nodes that form the vertices of the element.
    Node *n1;
    Node *n2;
    Node *n3;

    // Two vectors in the triangle
    Matrix2x2<double> currentState;
    // Inverse of the referece state
    Matrix2x2<double> invRefState;

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
    double energy;

    // UNUSED TODO When implemented, rewrite comment
    // Flag indicating if the cell can undergo plastic (permanent) deformation.
    bool plasticity;

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

private:
    // updates current state using two vectors from the triangle
    void m_updateCurrentState();

    // Computes the deformation gradient for the cell based on the triangle's vertices.
    void m_updateDeformationGradiant();

    // Computes the metric tensor for the triangle.
    void m_updateMetricTensor();

    // Performs a Lagrange reduction on C to calculate C_.
    void m_lagrangeReduction();

    // Calculates energy and reduced stress
    void m_calculateEnergyAndReducedStress();
    void m_UNUSED_calculate_energy_and_reduced_stress();

    // Calculate Piola stress P
    void m_updatePiolaStress();

};

std::ostream &operator<<(std::ostream &os, const TElement &element);

#endif