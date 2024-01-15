#include "tElement.h"

TElement::TElement(Node *n1, Node *n2, Node *n3) : n1(n1), n2(n2), n3(n3)
{
    // In order to calculate F later, we save the inverse of the jacobian.
    m_calculateJacobian();
    invJacobianRef = J.inverse();
    // We also need some adjustment vectors for calculating the force on each
    // node from the piola stress tensor. See applyForcesOnNodes for details.
    r1 = invJacobianRef.transpose() * b1;
    r2 = invJacobianRef.transpose() * b2;
    r3 = invJacobianRef.transpose() * b3;
}

TElement::TElement() {}

void TElement::update()
{

    // We use this public update function and hide the individual update
    // functions because the order here is very important.

    // Updates current state
    m_calculateJacobian();

    // Calculates F
    m_updateDeformationGradiant();

    // Calculates C
    m_updateMetricTensor();

    // Calculate C_ and m
    m_fastLagrangeReduction();

    // Calculate energy
    m_calculateEnergy();

    // Calculate reduced stress
    m_calculateReducedStress();

    // Calculate Piola stress P
    m_updatePiolaStress();

};

// vectors between the three nodes of the element
std::array<double, 2> TElement::e12() const
{
    return {
        n2->x - n1->x,
        n2->y - n1->y,
    };
}

std::array<double, 2> TElement::e13() const
{
    return {
        n3->x - n1->x,
        n3->y - n1->y,
    };
}

std::array<double, 2> TElement::e23() const
{
    return {
        n3->x - n2->x,
        n3->y - n2->y,
    };
}

/**
 * Jacobian Matrix Calculation for Finite Element Analysis
 *
 * Given shape functions:
 * N1 = 1 - ξ1 - ξ2
 * N2 = ξ1
 * N3 = ξ2
 *
 * x_1 = N1*x1 + N2*x2 + N3*x3 = (1 - ξ1 - ξ2)*x1 + ξ1*x2 + ξ2*x3
 * x_2 = N1*y1 + N2*y2 + N3*y3 = (1 - ξ1 - ξ2)*y1 + ξ1*y2 + ξ2*y3
 *
 * where xi and yi are the values at nodes i = 1,2,3, and x_1 and x_2 are the
 * two components of vector x. In our case, x_1 and x_2 are respectively the x
 * and y values of the nodes.
 *
 * Jacobian Matrix:
 * A = [ [∂x_1/∂e1, ∂x_1/∂e2],
 *       [∂x_2/∂e1, ∂x_2/∂e2] ],
 *
 * ∂x_1/∂e1 = -x1 + x2
 * ∂x_1/∂e2 = -x1 + x3
 * ∂x_2/∂e1 = -y1 + y2
 * ∂x_2/∂e2 = -y1 + y3
 *
 * giving us
 *
 * A = [ [-x1 + x2, -x1 + x3],
 *       [-y1 + y2, -y1 + y3] ]
 *
 * It just so happens that this can be expressed by simply using e12 and e13
 */
void TElement::m_calculateJacobian()
{
    J.setCols(e12(), e13());
}

/** TODO: Now that we are using jacobians instead of reference states. This explination needs to be updated!
 * How do we calculate F? First, we need a representation of our initial state.
 * To do this, we take any two vectors in the triangular element.
 * We use this reference state of the element to calculate the inverse of the
 * reference state: invJacobianRef. Using this together with the current state,
 * we can extract any deformations that have occured since initializing the
 * element. Example: A is reference state, C is the current state, S is deformation:
 *      If there have been no changes made to the element, then C=A, so
 *          S=AA^-1=I
 *      as expected.
 *      In general, we have that:
 *          S=CA^-1
 */
void TElement::m_updateDeformationGradiant()
{
    F = J * invJacobianRef;
}

// Provices a metric tensor for the triangle
void TElement::m_updateMetricTensor()
{
    // Discontinuous yielding of pristine micro-crystals - page 8/207
    C = F.transpose() * F;
}

void TElement::m_lagrangeReduction()
{
    // Homogeneous nucleation of dislocations as a pattern formation phenomenon - page 5
    // We start by copying the values from C to the reduced matrix
    C_ = C;
    // We should also reset m
    m = m.identity();

    // Note that we only modify C_[0][1]. At the end of the algorithm,
    // we copy C_[0][1] to C_[1][0]
    bool changed = true;
    while (changed)
    {
        changed = false;

        if (C_[0][1] < 0)
        {
            C_.flip(0, 1);
            m.lag_m1();
            changed = true;
        }

        if (C_[1][1] < C_[0][0])
        {
            C_.swap(0, 0, 1, 1);
            m.lag_m2();
            changed = true;
        }

        if (2 * C_[0][1] > C_[0][0])
        {
            // The order here matters, don't modify C_[0][1] before using it
            // to calculate C_[1][1].
            C_[1][1] += C_[0][0] - 2 * C_[0][1];
            C_[0][1] -= C_[0][0];
            m.lag_m3();
            changed = true;
        }
    }
    C_[1][0] = C_[0][1];
}

void TElement::m_fastLagrangeReduction()
{
    // If the ratio between c11 and c22 is very large, the lagrange reduction
    // algorithm takes a very long time to complete. The effective behavioure of
    // the reduction algorithm usually looks something like
    // m1 m2 m3 m3 m3 ... m3 m3 m3 m1
    // where m(1,2,3) denotes what condition is reached at each step in the
    // lagrange reduction algorithm. This long chain of m3 opperations are very
    // inefficient to do one by one, but using math, we can find how many we need
    // and do them all at once.

    // This algorithm only works when we only need to do either one or
    // zero m2 transformations.
    // TODO: I have not yet proven that this will always be the case when
    // c11 or c22 is smaller than 1, but i belive it is true.
    // Turns out that it is not true for <1, but maybe 0.5? :/
    if ((C[0][0] < 1 || C[1][1] < 1) && false)
    {
        // We start by copying the values from C to the reduced matrix
        C_ = C;
        // Reset the transformation matrix m to identity
        m = m.identity();

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

        double a = C_[0][0];
        double b = C_[0][1];
        double d = C_[1][1];
        // This is the number of times we can apply ... TODO
        int N = static_cast<int>(std::ceil(b / a - 0.5));
        C_[1][1] = -N * (b - a * N) - b * N + d;
        C_[0][1] = b - N * a;

        m.lag_m3(N);

        if (C_[0][1] < 0)
        {
            C_.flip(0, 1);
            m.lag_m1();
        }

        C_[1][0] = C_[0][1];

        Matrix2x2<double> testM = m;
        Matrix2x2<double> testC_ = C_;
        // m_lagrangeReduction();
        if(testM != m || testC_ != C_){
            // throw std::runtime_error("Error in fast lagrange reduction");
        }


    }
    else
    {
        // Use the normal Lagrange reduction process
        m_lagrangeReduction(); // Assuming this function handles v1 and v2 internally
    }
}

void TElement::m_calculateEnergy()
{
    // Uses the reduced metrics
    double c11 = C_[0][0];
    double c22 = C_[1][1];
    double c12 = C_[0][1];

    energy = -K * (log((c11 * c22 - c12 * c12) / burgers) - (c11 * c22 - c12 * c12) / burgers) + beta * (pow(1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22), 2.0) * 9.46969696969697E-4 - pow((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12), 3.0) * (4.1E1 / 9.9E1) + ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 + c22, 4.0) * (1.0 / 8.1E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * (7.0 / 1.98E2)) + pow(1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22), 2.0) * (1.7E1 / 5.28E2) + pow((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12), 3.0) * (4.0 / 1.1E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 + c22, 3.0) * (1.0 / 2.7E1) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * (8.0 / 3.3E1);
    // TODO zeroing energy
    // energy -= s.zeroing_energy;
}
void TElement::m_calculateReducedStress()
{
    // Uses the reduced metrics
    double c11 = C_[0][0];
    double c22 = C_[1][1];
    double c12 = C_[0][1];

    r_s[0][0] = -beta * (pow((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12), 2.0) * ((c11 * (1.0 / 6.0) - c12 * (2.0 / 3.0) + c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) - c22 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (4.1E1 / 3.3E1) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) + c22 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (1.0 / 5.28E2) - ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 + c22, 3.0) * (4.0 / 8.1E1) - 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 + c22, 4.0) * ((c11 * (1.0 / 6.0) - c12 * (2.0 / 3.0) + c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) - c22 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (1.0 / 8.1E1) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (7.0 / 1.98E2) + c22 * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0) * pow(c11 - c12 + c22, 4.0) * (2.0 / 8.1E1) - ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) + c22 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (7.0 / 1.98E2) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * ((c11 * (1.0 / 6.0) - c12 * (2.0 / 3.0) + c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) - c22 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (7.0 / 1.98E2) - c22 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 + c22) * (7.0 / 3.96E2)) + K * (c22 / burgers - c22 / (c11 * c22 - c12 * c12)) + pow((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12), 2.0) * ((c11 * (1.0 / 6.0) - c12 * (2.0 / 3.0) + c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) - c22 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (1.2E1 / 1.1E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) + c22 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (1.7E1 / 2.64E2) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 + c22, 3.0) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) + c22 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (1.0 / 2.7E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 + c22, 2.0) * (1.0 / 9.0) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (8.0 / 3.3E1) + c22 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 + c22, 3.0) * (1.0 / 1.8E1) - ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) + c22 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (8.0 / 3.3E1) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * ((c11 * (1.0 / 6.0) - c12 * (2.0 / 3.0) + c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) - c22 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (8.0 / 3.3E1) - c22 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 + c22) * (4.0 / 3.3E1);

    r_s[1][1] = beta * (pow((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12), 2.0) * ((c11 * (-1.0 / 6.0) + c12 * (2.0 / 3.0) - c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) + c11 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (4.1E1 / 3.3E1) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * (-pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) - c11 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (1.0 / 5.28E2) + ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 + c22, 3.0) * (4.0 / 8.1E1) - 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 + c22, 4.0) * ((c11 * (-1.0 / 6.0) + c12 * (2.0 / 3.0) - c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) + c11 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (1.0 / 8.1E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (7.0 / 1.98E2) - c11 * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0) * pow(c11 - c12 + c22, 4.0) * (2.0 / 8.1E1) - ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * (-pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) - c11 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (7.0 / 1.98E2) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * ((c11 * (-1.0 / 6.0) + c12 * (2.0 / 3.0) - c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) + c11 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (7.0 / 1.98E2) + c11 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 + c22) * (7.0 / 3.96E2)) + K * (c11 / burgers - c11 / (c11 * c22 - c12 * c12)) - pow((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12), 2.0) * ((c11 * (-1.0 / 6.0) + c12 * (2.0 / 3.0) - c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) + c11 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (1.2E1 / 1.1E1) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * (-pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) - c11 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (1.7E1 / 2.64E2) - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 + c22, 3.0) * (-pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) - c11 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (1.0 / 2.7E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 + c22, 2.0) * (1.0 / 9.0) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (8.0 / 3.3E1) + c11 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 + c22, 3.0) * (1.0 / 1.8E1) + ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * (-pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) - c11 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (8.0 / 3.3E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * ((c11 * (-1.0 / 6.0) + c12 * (2.0 / 3.0) - c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) + c11 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (8.0 / 3.3E1) - c11 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 + c22) * (4.0 / 3.3E1);

    r_s[1][0] = r_s[0][1] = pow((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12), 2.0) * ((c11 * (-2.0 / 3.0) + c12 * (8.0 / 3.0) - c22 * (2.0 / 3.0)) / (c11 * c22 - c12 * c12) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 6.0) + c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 2.0)) * (1.2E1 / 1.1E1) - K * ((c12 * 2.0) / burgers - (c12 * 2.0) / (c11 * c22 - c12 * c12)) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * 4.0 - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (4.0 / 3.0) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 3.0) - c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * 3.0) * (1.7E1 / 2.64E2) - beta * (pow((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12), 2.0) * ((c11 * (-2.0 / 3.0) + c12 * (8.0 / 3.0) - c22 * (2.0 / 3.0)) / (c11 * c22 - c12 * c12) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 6.0) + c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 2.0)) * (4.1E1 / 3.3E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * 4.0 - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (4.0 / 3.0) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 3.0) - c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * 3.0) * (1.0 / 5.28E2) + ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 + c22, 3.0) * (4.0 / 8.1E1) - 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * ((c11 * (-2.0 / 3.0) + c12 * (8.0 / 3.0) - c22 * (2.0 / 3.0)) / (c11 * c22 - c12 * c12) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 6.0) + c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 2.0)) * pow(c11 - c12 + c22, 4.0) * (1.0 / 8.1E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (7.0 / 1.98E2) - c12 * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0) * pow(c11 - c12 + c22, 4.0) * (4.0 / 8.1E1) + ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * 4.0 - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (4.0 / 3.0) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 3.0) - c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * 3.0) * (7.0 / 1.98E2) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * ((c11 * (-2.0 / 3.0) + c12 * (8.0 / 3.0) - c22 * (2.0 / 3.0)) / (c11 * c22 - c12 * c12) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 6.0) + c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 2.0)) * (c11 - c12 + c22) * (7.0 / 1.98E2) + c12 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 + c22) * (7.0 / 1.98E2)) - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 + c22, 3.0) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * 4.0 - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (4.0 / 3.0) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 3.0) - c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * 3.0) * (1.0 / 2.7E1) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 + c22, 2.0) * (1.0 / 9.0) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (8.0 / 3.3E1) - c12 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 + c22, 3.0) * (1.0 / 9.0) + ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * 4.0 - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (4.0 / 3.0) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 3.0) - c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * 3.0) * (8.0 / 3.3E1) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * ((c11 * (-2.0 / 3.0) + c12 * (8.0 / 3.0) - c22 * (2.0 / 3.0)) / (c11 * c22 - c12 * c12) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 6.0) + c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 2.0)) * (c11 - c12 + c22) * (8.0 / 3.3E1) + c12 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 + c22) * (8.0 / 3.3E1);
}

// Calculate Piola stress tensor and force on each node from current cell
void TElement::m_updatePiolaStress()
{
    //  Discontinuous yielding of pristine micro-crystals, page 16/215
    P = 2.0 * F * m * r_s * m.transpose();
}

// Note that each node is part of multiple elements. Therefore, the force must
// be reset after each itteration.
void TElement::applyForcesOnNodes()
{
    n1->addForce(P * r1);
    n2->addForce(P * r2);
    n3->addForce(P * r3);
}



// The functions below are not used in the simulation


double TElement::calculateEnergy(double c11, double c22, double c12)
{
    TElement element = TElement();
    element.C = {{c11, c12}, {c12, c22}};
    element.m_fastLagrangeReduction();
    element.m_calculateEnergy();
    return element.energy;
}

TElement TElement::lagrangeReduction(double c11, double c22, double c12)
{
    TElement element = TElement();
    element.C = {{c11, c12}, {c12, c22}};
    element.m_fastLagrangeReduction();
    return element;
}

std::ostream &operator<<(std::ostream &os, const TElement &element)
{
    os << std::fixed << std::setprecision(2); // Set precision to 2 decimal places
    os << "Energy: " << element.energy << "\t|";
    os << "n1: (" << element.n1->x << ", " << element.n1->y << "),\t"
       << "n2: (" << element.n2->x << ", " << element.n2->y << "),\t"
       << "n3: (" << element.n3->x << ", " << element.n3->y << ")";
    return os;
}