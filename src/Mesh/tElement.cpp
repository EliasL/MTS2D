#include "tElement.h"
#include "mesh.h"

TElement::TElement(Node n1, Node n2, Node n3)
    : nodes{n1, n2, n3}
{
    // Precompute this constant expression
    dxi_dX = dX_dxi().inverse();

    // Initialize the adjustment vectors
    for (size_t i = 0; i < r.size(); ++i)
    {
        r[i] = dxi_dX.transpose() * b[i];
    }
}

void TElement::update(Mesh &mesh)
{
    // The order here is very important.
    m_updatePosition(mesh);
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

// Position subtraction (The vector from node 1 to node 2)
VArray TElement::dx(Node &n1, Node &n2) const
{
    return n2.pos() - n1.pos();
}

// Initial-position subtraction
VArray TElement::dX(Node &n1, Node &n2) const
{
    return n2.init_pos() - n1.init_pos();
}

// Displacement subtraction
VArray TElement::du(Node &n1, Node &n2) const
{
    return n2.u() - n1.u();
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
    du_dxi.setCols(du(nodes[0], nodes[1]), du(nodes[0], nodes[2]));
    return du_dxi;
}

void TElement::copyValues(const TElement &other)
{
    // Copy each field from 'other' to 'this'
    this->F = other.F;
    this->C = other.C;
    this->C_ = other.C_;
    this->m = other.m;
    this->r_s = other.r_s;
    this->P = other.P;
    this->energy = other.energy;
    this->resolvedShearStress = other.resolvedShearStress;
    this->dxi_dX = other.dxi_dX;

    this->r = other.r;
    this->m3Nr = other.m3Nr;
    this->past_m3Nr = other.past_m3Nr;
}

/**
 * Jacobian with respect to the initial position of the nodes ∂X/∂ξ
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
    dX_dxi.setCols(dX(nodes[0], nodes[1]), dX(nodes[0], nodes[2]));
    return dX_dxi;
}

void TElement::m_updatePosition(Mesh &mesh)
{
    for (size_t i = 0; i < nodes.size(); i++)
    {
        Node *n = mesh[nodes[i].id];
        if (nodes[i].ghostNode)
        {
            nodes[i].setPos(mesh.makeGhostPos(n->pos(), nodes[i].ghostShift));
        }
        else
        {
            nodes[i].setPos(n->pos());
        }
    }
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
void TElement::applyForcesOnNodes(Mesh &mesh)
{
    // TODO explain what is going on here
    for (int i = 0; i < 3; i++)
    {
        mesh[nodes[i].id]->addForce(P * r[i]);
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

// std::ostream &operator<<(std::ostream &os, const TElement &element)
// {
//     // Save the current format state of the stream
//     std::ios_base::fmtflags f(os.flags());

//     // Save the current format state of the stream
//     std::streamsize prec = os.precision();

//     os << std::fixed << std::setprecision(2); // Set precision to 2 decimal places
//     os << "Energy: " << element.energy << "\t|";
//     for (size_t i = 0; i < element.nodes.size(); ++i)
//     {
//         VArray pos = element.nodes[i].pos();
//         os << "n" << (i + 1) << ": ("
//            << pos[0] << ", " << pos[0] << ")";
//         if (i < element.nodes.size() - 1)
//         {
//             os << ",\t";
//         }
//     }
//     // Restore the saved precision state
//     os.precision(prec);
//     os.flags(f);
//     return os;
// }
