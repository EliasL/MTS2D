#include "cell.h"

Cell::Cell(const Triangle &triangle)
{
    // Calculates D
    m_getDeformationGradiant(triangle);

    // Calculates C
    C = triangle.metric(METRICFUNCTION);

    // Calculate C_ and m
    m_lagrangeReduction();
};

Cell::Cell()
{
}

// Calculate Piola stress tensor and force on each node from current cell
// Note that each node is part of multiple cells. Therefore, the force must
// be reset after each itteration.
void Cell::setForcesOnNodes(Triangle &triangle)
{

    if (!hasComputedReducedStress)
        throw std::runtime_error("Reduced stress has not yet been calculated");

    // extended stress is not quite the "real" stress, but it is a component
    // in calculating the piola stress, which is the real stress on the cell,
    // and we can then find the force on each individual node.
    // The name extended_stress does not have much meaning.
    // TODO consider storing this variable in the cell, such that it does
    // not have to be allocated every time the function is called.
    Matrix2x2<double> extended_stress = r_s.sym_orth_conjugate(m);
    // TODO THIS IS WRONG
    P[0][0] = 2 * extended_stress[0][0] * F[0][0] + extended_stress[0][1] * F[1][0];
    P[1][0] = 2 * extended_stress[0][0] * F[0][1] + extended_stress[0][1] * F[1][1];
    P[0][1] = 2 * extended_stress[1][1] * F[1][0] + extended_stress[0][1] * F[0][0];
    P[1][1] = 2 * extended_stress[1][1] * F[1][1] + extended_stress[0][1] * F[0][1];
    // The assignment here is dependant on the shape of the cell.
    // For triangular shapes, the forces on the nodes is applied as shown
    // below. For a general shape, see Gael-notes page 2, partial N^i / partial x_j
    // on how to calculate.

    // DO GENERAL DISTRIBUTION of Piola stress

    // FORCES SHOULD BE RESET AFTER EACH ITERATION
    // MUST SUM (+=) BECAUSE THEY ARE THE SAME NODES THAT ARE IN A FLIPPED TRIANGLE
    triangle.a1->f_x += -P[0][0] - P[0][1];
    triangle.a1->f_y += -P[1][0] - P[1][1];

    triangle.a2->f_x += P[0][0];
    triangle.a2->f_y += P[1][0];

    triangle.a3->f_x += P[0][1];
    triangle.a3->f_y += P[1][1];
}

void Cell::m_getDeformationGradiant(const Triangle &triangle)
{
    auto e1_ = triangle.e12();
    auto e2_ = triangle.e13();

    F[0][0] = e1_[0];
    F[1][0] = e1_[1];
    F[0][1] = e2_[0];
    F[1][1] = e2_[1];
}

void Cell::m_lagrangeReduction()
{
    // We start by copying the values from C to the reduced matrix
    C_ = C;

    if (LINEARITY)
    {
        // If we assume linearity, we are done. m is already identity.
        return;
    }
    // And then we follow an algorithm generate both m and C_
    while (C_[0][1] < 0 || C_[1][1] < C_[0][0] || 2 * C_[0][1] > C_[0][0])
    {

        if (C_[0][1] < 0)
        {
            C_.flip(0, 1);
            m.lag_m1();
        }

        if (C_[1][1] < C_[0][0])
        {
            C_.swap(0, 0, 1, 1);
            m.lag_m2();
        }

        if (2 * C_[0][1] > C_[0][0])
        {
            // The order here matters, don't modify C_[0][1] before using it
            // to calculate C_[1][1].
            C_[1][1] += C_[0][0] - 2 * C_[0][1];
            C_[0][1] -= C_[0][0];
            m.lag_m3();
        }
    }
}