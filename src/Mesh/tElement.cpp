#include "tElement.h"

TElement::TElement(Node *n1, Node *n2, Node *n3) : nodes{n1, n2, n3}
{
    // Precompute this constant expression
    dxi_dX = dX_dxi().inverse();

    // Initialize the adjustment vectors
    for (size_t i = 0; i < r.size(); ++i)
    {
        r[i] = dxi_dX.transpose() * b[i];
    }
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
// Generalized function
std::array<double, 2> TElement::vectorBetweenNodes(
    std::function<double(Node *)> getX,
    std::function<double(Node *)> getY,
    int idx1,
    int idx2) const
{
    double x1 = getX(nodes[idx1]);
    double y1 = getY(nodes[idx1]);
    double x2 = getX(nodes[idx2]);
    double y2 = getY(nodes[idx2]);

    // Applying offset if necessary
    if (moveN[idx1])
    {
        x1 += xOffset;
        y1 += yOffset;
    }
    if (moveN[idx2])
    {
        x2 += xOffset;
        y2 += yOffset;
    }

    return {x2 - x1, y2 - y1};
}

// These functions look scary because of the lambda functions, but just
// remember that they are only the vectors between nodes

// Displacement difference
std::array<double, 2> TElement::u(int idx1, int idx2) const
{
    return vectorBetweenNodes(
        [](const Node *n)
        { return n->U_x(); },
        [](const Node *n)
        { return n->U_y(); },
        idx1, idx2);
}

// Position difference
std::array<double, 2> TElement::x(int idx1, int idx2) const
{
    return vectorBetweenNodes(
        [](const Node *n)
        { return n->X(); },
        [](const Node *n)
        { return n->Y(); },
        idx1, idx2);
}

// Initial-position difference
std::array<double, 2> TElement::X(int idx1, int idx2) const
{
    return vectorBetweenNodes(
        [](const Node *n)
        { return n->Init_x(); },
        [](const Node *n)
        { return n->Init_y(); },
        idx1, idx2);
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
    du_dxi.setCols(u(0, 1), u(0, 2));
    return du_dxi;
}
/**
 * Jacobian with respect to the initial position of the nodes ∂x/∂ξ
 *
 * Given shape functions:
 * N1 = 1 - ξ1 - ξ2
 * N2 = ξ1
 * N3 = ξ2
 *
 * X_1 = N1*x1 + N2*x2 + N3*x3 = (1 - ξ1 - ξ2)*x1 + ξ1*x2 + ξ2*x3
 * X_2 = N1*y1 + N2*y2 + N3*y3 = (1 - ξ1 - ξ2)*y1 + ξ1*y2 + ξ2*y3
 *
 * where xi and yi are the reference positions of nodes i = 1,2,3, and
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
    dX_dxi.setCols(X(0, 1), X(0, 2));
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
    energy = ContiPotential::energy(C_[0][0], C_[1][1], C_[0][1], beta, mu);
}

void TElement::m_updateReducedStress()
{
    r_s = ContiPotential::stress(C_[0][0], C_[1][1], C_[0][1], beta, mu);
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
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        // TODO explain what is going on here
        nodes[i]->addForce(P * r[i]);
    }
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
    for (size_t i = 0; i < element.nodes.size(); ++i)
    {
        os << "n" << (i + 1) << ": ("
           << element.nodes[i]->X() << ", "
           << element.nodes[i]->Y() << ")";
        if (i < element.nodes.size() - 1)
        {
            os << ",\t";
        }
    }
    // Restore the saved precision state
    os.precision(prec);
    os.flags(f);
    return os;
}