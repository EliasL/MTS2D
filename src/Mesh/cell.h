#ifndef CELL_H
#define CELL_H
#pragma once

#include "settings.h"
#include "Matrix/matrix.h"
#include "Matrix/matrix2x2.h"
#include "node.h"
#include "triangle.h"
#include "easylogging++.h"
#include <array>
#include <vector>
#include <stdexcept>


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

    bool hasComputedReducedStress = false;

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
    Cell(const Triangle &t);

    // Default constructor for the Cell.
    Cell();

    // Sets the forces on the nodes that form the cell's triangle.
    void setForcesOnNodes(Triangle &t);

private:
    // Computes the deformation gradient for the cell based on the triangle's vertices.
    void m_getDeformationGradiant(const Triangle &triangle);

    // Performs a Lagrange reduction on C to calculate C_.
    void m_lagrangeReduction();
};

#endif