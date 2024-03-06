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
    m_lagrangeReduction();

    // Calculate energy
    m_updateEnergy();

    // Calculate reduced stress
    m_updateReducedStress();

    // Calculate Piola stress P
    m_updatePiolaStress();

    // Calculate resolved shear stress
    m_updateResolvedShearStress();
};

// vectors between the three nodes of the element

// Current Possition
std::array<double, 2> TElement::x12() const
{
    return {
        n2->X() - n1->X(),
        n2->Y() - n1->Y(),
    };
}
std::array<double, 2> TElement::x13() const
{
    return {
        n3->X() - n1->X(),
        n3->Y() - n1->Y(),
    };
}

// Displacement
std::array<double, 2> TElement::u12() const
{
    return {
        n2->U_x() - n1->U_x(),
        n2->U_y() - n1->U_y(),
    };
}
std::array<double, 2> TElement::u13() const
{
    return {
        n3->U_x() - n1->U_x(),
        n3->U_y() - n1->U_y(),
    };
}

// Initial Position
std::array<double, 2> TElement::X12() const
{
    return {
        n2->Init_x() - n1->Init_x(),
        n2->Init_y() - n1->Init_y(),
    };
}
std::array<double, 2> TElement::X13() const
{
    return {
        n3->Init_x() - n1->Init_x(),
        n3->Init_y() - n1->Init_y(),
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

    //                ∂u/∂X = ∂u/∂ξ * ∂ξ/∂X
    Matrix2x2<double> du_dX = du_dxi() * dxi_dX;
    F = Matrix2x2<double>::identity() + du_dX;
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

void TElement::m_updateEnergy()
{
    energy = nameOfEnergyFunction::polynomialEnergy(C_[0][0], C_[1][1], C_[0][1], beta, mu);
}

void TElement::m_updateReducedStress()
{
    r_s = nameOfEnergyFunction::polynomialStress(C_[0][0], C_[1][1], C_[0][1], beta, mu);
}

// Calculate Piola stress tensor and force on each node from current cell
void TElement::m_updatePiolaStress()
{
    //  Discontinuous yielding of pristine micro-crystals, page 16/215
    P = 2.0 * F * m * r_s * m.transpose();
}

void TElement::m_updateResolvedShearStress()
{
    /**  Discontinuous yielding of pristine micro-crystals (page 216/17)
     *  resolved-shear stress = ∂W/∂α = ∫_Ω P:(∂F/∂α)dx
     *
     * We assume that the shear ∂F/∂α is always
     * ∂F/∂α = [ [0, 1],
     *           [0, 0] ]
     */
    resolvedShearStress = P[0][1];
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

bool TElement::plasticEvent()
{
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

// The functions below are not used in the simulation

double TElement::calculateEnergy(double c11, double c22, double c12)
{
    TElement element = TElement();
    element.C = {{c11, c12}, {c12, c22}};
    element.m_lagrangeReduction();
    element.m_updateEnergy();
    return element.energy;
}

TElement TElement::lagrangeReduction(double c11, double c22, double c12)
{
    TElement element = TElement();
    element.C = {{c11, c12}, {c12, c22}};
    element.m_lagrangeReduction();
    return element;
}

std::ostream &operator<<(std::ostream &os, const TElement &element)
{
    // Save the current format state of the stream
    std::ios_base::fmtflags f(os.flags());

    // Save the current format state of the stream
    std::streamsize prec = os.precision();

    os << std::fixed << std::setprecision(2); // Set precision to 2 decimal places
    os << "Energy: " << element.energy << "\t|";
    os << "n1: (" << element.n1->X() << ", " << element.n1->Y() << "),\t"
       << "n2: (" << element.n2->X() << ", " << element.n2->Y() << "),\t"
       << "n3: (" << element.n3->X() << ", " << element.n3->Y() << ")";
    // Restore the saved precision state
    os.precision(prec);
    os.flags(f);
    return os;
}