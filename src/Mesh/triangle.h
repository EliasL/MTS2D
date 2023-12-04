#ifndef TRIANGLE_H
#define TRIANGLE_H
#pragma once

#include "settings.h"
#include "Matrix/matrix.h"
#include "Matrix/matrix2x2.h"
#include "node.h"
#include "easylogging++.h"
#include <array>
#include <vector>
#include <stdexcept>



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
    
    friend std::ostream& operator<<(std::ostream &os, const Triangle &triangle);
    friend std::ostream& operator<<(std::ostream &os, const Triangle *trianglePtr);
};

#endif