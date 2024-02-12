#include "tElement.h"

TElement::TElement(Node *n1, Node *n2, Node *n3) : n1(n1), n2(n2), n3(n3)
{
    // In order to calculate F later, we save the inverse of the jacobian of the
    // initial state of the element.
    dxi_dX = dX_dxi().inverse();
    // We also need some adjustment vectors for calculating the force on each
    // node from the piola stress tensor. See applyForcesOnNodes for details.
    r1 = dxi_dX.transpose() * b1;
    r2 = dxi_dX.transpose() * b2;
    r3 = dxi_dX.transpose() * b3;
}

TElement::TElement() {}

void TElement::update()
{
    // The order here is very important.

    // Calculates F
    m_updateDeformationGradiant();

    // Calculates C
    m_updateMetricTensor();

    // Calculate C_ and m
    m_fastLagrangeReduction();

    // If C_ is found in the rlu cache, we can load the answers for the
    // energy and the reduced stress

    // Calculate energy
    m_updateEnergy();

    // Calculate reduced stress
    m_updateReducedStress();

    // Calculate Piola stress P
    m_updatePiolaStress();
};

// vectors between the three nodes of the element

// Current Possition
std::array<double, 2> TElement::x12() const
{
    return {
        n2->x - n1->x,
        n2->y - n1->y,
    };
}
std::array<double, 2> TElement::x13() const
{
    return {
        n3->x - n1->x,
        n3->y - n1->y,
    };
}

// Displacement
std::array<double, 2> TElement::u12() const
{
    return {
        n2->u_x() - n1->u_x(),
        n2->u_y() - n1->u_y(),
    };
}
std::array<double, 2> TElement::u13() const
{
    return {
        n3->u_x() - n1->u_x(),
        n3->u_y() - n1->u_y(),
    };
}

// Initial Position
std::array<double, 2> TElement::X12() const
{
    return {
        n2->init_x - n1->init_x,
        n2->init_y - n1->init_y,
    };
}
std::array<double, 2> TElement::X13() const
{
    return {
        n3->init_x - n1->init_x,
        n3->init_y - n1->init_y,
    };
}

/**
 * Jacobian of the displacements ∂u/∂ξ
 *
 * Given shape functions:
 * N1 = 1 - ξ1 - ξ2
 * N2 = ξ1
 * N3 = ξ2
 *
 * u_1 = N1*u_x1 + N2*u_x2 + N3*u_x3 = (1 - ξ1 - ξ2)*u_x1 + ξ1*u_x2 + ξ2*u_x3
 * u_2 = N1*u_y1 + N2*u_y2 + N3*u_y3 = (1 - ξ1 - ξ2)*u_y1 + ξ1*u_y2 + ξ2*u_y3
 *
 * where u_xi and u_yi are the displacement values at nodes i = 1,2,3, 
 * and u_1 and u_2 are the two components of the displacement vector u.
 *
 * Jacobian Matrix:
 * J = [ [∂u_1/∂ξ1, ∂u_1/∂ξ2],
 *       [∂u_2/∂ξ1, ∂u_2/∂ξ2] ],
 *
 * ∂u_1/∂ξ1 = -u_x1 + u_x2
 * ∂u_1/∂ξ2 = -u_x1 + u_x3
 * ∂u_2/∂ξ1 = -u_y1 + u_y2
 * ∂u_2/∂ξ2 = -u_y1 + u_y3
 *
 * giving us
 *
 * J = [ [-u_x1 + u_x2, -u_x1 + u_x3],
 *       [-u_y1 + u_y2, -u_y1 + u_y3] ]
 * 
 * It just so happens that this can be expressed by simply using u12 and u13
 */
Matrix2x2<double> TElement::du_dxi()
{    
    // ∂u/∂ξ
    Matrix2x2<double> du_dxi;
    du_dxi.setCols(u12(), u13());
    return du_dxi;
}
/**
 * Jacobian with respect to the initial possition of the nodes ∂x/∂ξ
 *
 * Given shape functions:
 * N1 = 1 - ξ1 - ξ2
 * N2 = ξ1
 * N3 = ξ2
 *
 * X_1 = N1*x1 + N2*x2 + N3*x3 = (1 - ξ1 - ξ2)*x1 + ξ1*x2 + ξ2*x3
 * X_2 = N1*y1 + N2*y2 + N3*y3 = (1 - ξ1 - ξ2)*y1 + ξ1*y2 + ξ2*y3
 *
 * where xi and yi are the reference possitions of nodes i = 1,2,3, and
 * X_1 and X_2 are the two components of X.
 *
 * Jacobian Matrix:
 * J = [ [∂X_1/∂ξ1, ∂X_1/∂ξ2],
 *       [∂X_2/∂ξ1, ∂X_2/∂ξ2] ],
 *
 * ∂X_1/∂ξ1 = -x1 + x2
 * ∂X_1/∂ξ2 = -x1 + x3
 * ∂X_2/∂ξ1 = -y1 + y2
 * ∂X_2/∂ξ2 = -y1 + y3
 *
 * giving us
 *
 * J = [ [-x1 + x2, -x1 + x3],
 *       [-y1 + y2, -y1 + y3] ]
 *
 * It just so happens that this can be expressed by simply using X12 and X13
 */
Matrix2x2<double> TElement::dX_dxi()
{
    // ∂X/∂ξ
    Matrix2x2<double> dX_dxi;
    dX_dxi.setCols(X12(), X13());
    return dX_dxi;
}


void TElement::m_updateDeformationGradiant()
{
    // See FEMNotes pdf from Umut
    
    //                ∂u/∂x = ∂u/∂ξ * ∂ξ/∂X
    Matrix2x2<double> du_dx = du_dxi() * dxi_dX;
    F = Matrix2x2<double>::identity() + du_dx;
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
    // We should also reset m and m3Nr
    m = m.identity();
    m3Nr = 0;

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
            m3Nr += 1;
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
    // For now, we don't use this code, we always go to the else statement
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
        if (testM != m || testC_ != C_)
        {
            // throw std::runtime_error("Error in fast lagrange reduction");
        }
    }
    else
    {
        // Use the normal Lagrange reduction process
        m_lagrangeReduction(); // Assuming this function handles v1 and v2 internally
    }
}

void TElement::m_updateEnergy()
{
    // Uses the reduced metrics
    double c11 = C_[0][0];
    double c22 = C_[1][1];
    double c12 = C_[0][1];

    energy = -K * (log((c11 * c22 - c12 * c12) / burgers) - (c11 * c22 - c12 * c12) / burgers) + beta * (pow(1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22), 2.0) * 9.46969696969697E-4 - pow((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12), 3.0) * (4.1E1 / 9.9E1) + ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 + c22, 4.0) * (1.0 / 8.1E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * (7.0 / 1.98E2)) + pow(1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22), 2.0) * (1.7E1 / 5.28E2) + pow((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12), 3.0) * (4.0 / 1.1E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 + c22, 3.0) * (1.0 / 2.7E1) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * (8.0 / 3.3E1);
    // TODO zeroing energy
    // energy -= s.zeroing_energy;
}

void TElement::m_updateReducedStress()   
{
    // Uses the reduced metrics
    double c11 = C_[0][0];
    double c22 = C_[1][1];
    double c12 = C_[0][1];
 
    // eq 7 Discontinuous yielding of pristine micro-crystals, page 15/214
    // C is the lagrange reduced metric and
    // Φ is the frame indifferent elastic energy density
    // ∂Φ/∂C
    double dPhi_dC11 = -beta * (pow((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12), 2.0) * ((c11 * (1.0 / 6.0) - c12 * (2.0 / 3.0) + c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) - c22 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (4.1E1 / 3.3E1) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) + c22 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (1.0 / 5.28E2) - ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 + c22, 3.0) * (4.0 / 8.1E1) - 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 + c22, 4.0) * ((c11 * (1.0 / 6.0) - c12 * (2.0 / 3.0) + c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) - c22 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (1.0 / 8.1E1) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (7.0 / 1.98E2) + c22 * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0) * pow(c11 - c12 + c22, 4.0) * (2.0 / 8.1E1) - ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) + c22 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (7.0 / 1.98E2) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * ((c11 * (1.0 / 6.0) - c12 * (2.0 / 3.0) + c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) - c22 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (7.0 / 1.98E2) - c22 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 + c22) * (7.0 / 3.96E2)) + K * (c22 / burgers - c22 / (c11 * c22 - c12 * c12)) + pow((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12), 2.0) * ((c11 * (1.0 / 6.0) - c12 * (2.0 / 3.0) + c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) - c22 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (1.2E1 / 1.1E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) + c22 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (1.7E1 / 2.64E2) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 + c22, 3.0) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) + c22 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (1.0 / 2.7E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 + c22, 2.0) * (1.0 / 9.0) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (8.0 / 3.3E1) + c22 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 + c22, 3.0) * (1.0 / 1.8E1) - ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) + c22 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (8.0 / 3.3E1) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * ((c11 * (1.0 / 6.0) - c12 * (2.0 / 3.0) + c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) - c22 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) - c22 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (8.0 / 3.3E1) - c22 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 + c22) * (4.0 / 3.3E1);
    
    double dPhi_dC22 = beta * (pow((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12), 2.0) * ((c11 * (-1.0 / 6.0) + c12 * (2.0 / 3.0) - c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) + c11 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (4.1E1 / 3.3E1) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * (-pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) - c11 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (1.0 / 5.28E2) + ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 + c22, 3.0) * (4.0 / 8.1E1) - 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 + c22, 4.0) * ((c11 * (-1.0 / 6.0) + c12 * (2.0 / 3.0) - c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) + c11 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (1.0 / 8.1E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (7.0 / 1.98E2) - c11 * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0) * pow(c11 - c12 + c22, 4.0) * (2.0 / 8.1E1) - ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * (-pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) - c11 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (7.0 / 1.98E2) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * ((c11 * (-1.0 / 6.0) + c12 * (2.0 / 3.0) - c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) + c11 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (7.0 / 1.98E2) + c11 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 + c22) * (7.0 / 3.96E2)) + K * (c11 / burgers - c11 / (c11 * c22 - c12 * c12)) - pow((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12), 2.0) * ((c11 * (-1.0 / 6.0) + c12 * (2.0 / 3.0) - c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) + c11 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (1.2E1 / 1.1E1) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * (-pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) - c11 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (1.7E1 / 2.64E2) - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 + c22, 3.0) * (-pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) - c11 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (1.0 / 2.7E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 + c22, 2.0) * (1.0 / 9.0) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (8.0 / 3.3E1) + c11 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 + c22, 3.0) * (1.0 / 1.8E1) + ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * (-pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 3.0) - c11 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 6.0) + 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 * 2.0 - c22 * 2.0) * (c11 - c12 * 4.0 + c22) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * (3.0 / 2.0)) * (8.0 / 3.3E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * ((c11 * (-1.0 / 6.0) + c12 * (2.0 / 3.0) - c22 * (1.0 / 6.0)) / (c11 * c22 - c12 * c12) + (c11 * (1.0 / 2.0) - c22 * (1.0 / 2.0)) / (c11 * c22 - c12 * c12) + c11 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1) + c11 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 4.0)) * (8.0 / 3.3E1) - c11 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 + c22) * (4.0 / 3.3E1);

    double dPhi_dC12 = pow((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12), 2.0) * ((c11 * (-2.0 / 3.0) + c12 * (8.0 / 3.0) - c22 * (2.0 / 3.0)) / (c11 * c22 - c12 * c12) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 6.0) + c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 2.0)) * (1.2E1 / 1.1E1) - K * ((c12 * 2.0) / burgers - (c12 * 2.0) / (c11 * c22 - c12 * c12)) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * 4.0 - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (4.0 / 3.0) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 3.0) - c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * 3.0) * (1.7E1 / 2.64E2) - beta * (pow((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12), 2.0) * ((c11 * (-2.0 / 3.0) + c12 * (8.0 / 3.0) - c22 * (2.0 / 3.0)) / (c11 * c22 - c12 * c12) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 6.0) + c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 2.0)) * (4.1E1 / 3.3E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * 4.0 - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (4.0 / 3.0) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 3.0) - c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * 3.0) * (1.0 / 5.28E2) + ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 + c22, 3.0) * (4.0 / 8.1E1) - 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * ((c11 * (-2.0 / 3.0) + c12 * (8.0 / 3.0) - c22 * (2.0 / 3.0)) / (c11 * c22 - c12 * c12) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 6.0) + c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 2.0)) * pow(c11 - c12 + c22, 4.0) * (1.0 / 8.1E1) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (7.0 / 1.98E2) - c12 * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0) * pow(c11 - c12 + c22, 4.0) * (4.0 / 8.1E1) + ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * 4.0 - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (4.0 / 3.0) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 3.0) - c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * 3.0) * (7.0 / 1.98E2) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * ((c11 * (-2.0 / 3.0) + c12 * (8.0 / 3.0) - c22 * (2.0 / 3.0)) / (c11 * c22 - c12 * c12) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 6.0) + c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 2.0)) * (c11 - c12 + c22) * (7.0 / 1.98E2) + c12 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 + c22) * (7.0 / 1.98E2)) - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 + c22, 3.0) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * 4.0 - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (4.0 / 3.0) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 3.0) - c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * 3.0) * (1.0 / 2.7E1) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 + c22, 2.0) * (1.0 / 9.0) - (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (8.0 / 3.3E1) - c12 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 + c22, 3.0) * (1.0 / 9.0) + ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * (c11 - c12 + c22) * (pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * 4.0 - 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (4.0 / 3.0) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 3.0) - c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 5.0 / 2.0) * (c11 - c12 * 4.0 + c22) * 3.0) * (8.0 / 3.3E1) + (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * 1.0 / sqrt(c11 * c22 - c12 * c12) * ((c11 * (-2.0 / 3.0) + c12 * (8.0 / 3.0) - c22 * (2.0 / 3.0)) / (c11 * c22 - c12 * c12) + c12 * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 6.0) + c12 * pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 2.0) * (1.0 / 2.0)) * (c11 - c12 + c22) * (8.0 / 3.3E1) + c12 * (1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * pow(c11 - c12 * 4.0 + c22, 3.0) * (1.0 / 9.0) - pow(c11 - c22, 2.0) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 * 4.0 + c22)) * ((pow(c11 - c22, 2.0) * (1.0 / 4.0)) / (c11 * c22 - c12 * c12) + (pow(c11 - c12 * 4.0 + c22, 2.0) * (1.0 / 1.2E1)) / (c11 * c22 - c12 * c12)) * 1.0 / pow(c11 * c22 - c12 * c12, 3.0 / 2.0) * (c11 - c12 + c22) * (8.0 / 3.3E1);
    
    r_s[0][0] = dPhi_dC11;
    r_s[1][1] = dPhi_dC22;
    r_s[1][0] = r_s[0][1] = dPhi_dC12/2;
}

// Calculate Piola stress tensor and force on each node from current cell
void TElement::m_updatePiolaStress()
{
    //  Discontinuous yielding of pristine micro-crystals, page 16/215
    P = 2.0 * F * m * r_s * m.transpose();
}

double TElement::m_calculateResolvedShearStress()
{
    /**  Discontinuous yielding of pristine micro-crystals (page 216/17)
     *  resolved-shear stress = ∂W/∂α = ∫_Ω P:(∂F/∂α)dx
     *
     * We assume that the shear ∂F/∂α is always
     * ∂F/∂α = [ [0, 1],
     *           [0, 0] ]
     */
    return P[0][1];
}

// Note that each node is part of multiple elements. Therefore, the force must
// be reset after each iteration.
void TElement::applyForcesOnNodes()
{
    // TODO explain what is going on here
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
    element.m_updateEnergy();
    return element.energy;
}

TElement TElement::lagrangeReduction(double c11, double c22, double c12)
{
    TElement element = TElement();
    element.C = {{c11, c12}, {c12, c22}};
    element.m_fastLagrangeReduction();
    return element;
}

bool TElement::plasticEvent()
{
    // Note that we never initialize past_m before here, so the first call
    // will always return true.
    if (m3Nr != past_m3Nr)
    {
        past_m3Nr = m3Nr;
        return true;
    }
    else
    {
        return false;
    }
}

std::ostream &operator<<(std::ostream &os, const TElement &element)
{

    // Save the current format state of the stream
    std::ios_base::fmtflags f(os.flags());

    // Save the current format state of the stream
    std::streamsize prec = os.precision();

    os << std::fixed << std::setprecision(2); // Set precision to 2 decimal places
    os << "Energy: " << element.energy << "\t|";
    os << "n1: (" << element.n1->x << ", " << element.n1->y << "),\t"
       << "n2: (" << element.n2->x << ", " << element.n2->y << "),\t"
       << "n3: (" << element.n3->x << ", " << element.n3->y << ")";
    // Restore the saved precision state
    os.precision(prec);
    os.flags(f);
    return os;
}