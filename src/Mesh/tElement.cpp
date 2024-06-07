#include "tElement.h"
#include "Simulation/energyFunctions.h"
#include "mesh.h"
#include <Eigen/LU>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <stdexcept>

/*
-1.0, 1.0, 0.0,
-1.0, 0.0, 1.0
*/
Matrix<double, 2, 3> TElement::b =
    (Matrix<double, 2, 3>() << -1.0, 1.0, 0.0, -1.0, 0.0, 1.0).finished();

TElement::TElement(Node n1, Node n2, Node n3, double _noise)
    : nodes{n1, n2, n3}, F(Matrix2d::Identity()), C(Matrix2d::Identity()),
      C_(Matrix2d::Identity()), m(Matrix2d::Identity()),
      r_s(Matrix2d::Identity()), P(Matrix2d::Identity()), noise(_noise) {
  // Precompute this constant expression
  dxi_dX = dX_dxi().inverse();

  // Initialize the adjustment vectors
  for (size_t i = 0; i < r.size(); ++i) {
    r[i] = dxi_dX.transpose() * b.col(i);
  }

  // Calculate initial area
  initArea = tElementArea(n1, n2, n3);
}

void TElement::update(Mesh &mesh) {
  // The order here is very important.

  // Update nodes in element
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
Vector2d TElement::dx(Node &n1, Node &n2) const { return n2.pos() - n1.pos(); }

// Initial-position subtraction
Vector2d TElement::dX(Node &n1, Node &n2) const {
  return n2.init_pos() - n1.init_pos();
}

// Displacement subtraction
Vector2d TElement::du(Node &n1, Node &n2) const { return n2.u() - n1.u(); }

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
Matrix2d TElement::du_dxi() {
  // ∂u/∂ξ
  Matrix2d du_dxi;
  du_dxi.col(0) = du(nodes[0], nodes[1]);
  du_dxi.col(1) = du(nodes[0], nodes[2]);

  return du_dxi;
}

// Jacobian with respect to the initial position of the nodes ∂X/∂ξ
// See du_dxi for a similar working out.
Matrix2d TElement::dX_dxi() {
  // ∂X/∂ξ
  Matrix2d dX_dxi;
  dX_dxi.col(0) = dX(nodes[0], nodes[1]);
  dX_dxi.col(1) = dX(nodes[0], nodes[2]);
  return dX_dxi;
}

void TElement::m_updatePosition(Mesh &mesh) {
  for (size_t i = 0; i < nodes.size(); i++) {
    Node *n = mesh[nodes[i].id];
    if (nodes[i].isGhostNode) {
      nodes[i].setPos(mesh.makeGhostPos(n->pos(), nodes[i].ghostShift));
    } else {
      nodes[i].setPos(n->pos());
    }
  }
}

void TElement::m_updateDeformationGradiant() {
  // See FEMNotes pdf from Umut

  //                ∂u/∂X = ∂u/∂ξ * ∂ξ/∂X
  Matrix2d du_dX = du_dxi() * dxi_dX;
  F = Matrix2d::Identity() + du_dX;
}

// Provices a metric tensor for the triangle
void TElement::m_updateMetricTensor() {
  // Discontinuous yielding of pristine micro-crystals - page 8/207
  C = F.transpose() * F;
}

void lag_m1(Matrix2d &mat) {
  // Multiply by 1  0
  //             0 -1
  mat(0, 1) = -mat(0, 1);
  mat(1, 1) = -mat(1, 1);
}

void lag_m2(Matrix2d &mat) {
  // Multiply by 0 1
  //             1 0
  mat.col(0).swap(mat.col(1));
}

void lag_m3(Matrix2d &mat, double n = -1) {
  // Multiply by 1 -1
  //             0  1
  mat = mat * Matrix2d{{1, n}, {0, 1}};
}

void TElement::m_lagrangeReduction() {
  // Homogeneous nucleation of dislocations as a pattern formation phenomenon -
  // page 5 We start by copying the values from C to the reduced matrix
  C_ = C;
  // We should also reset m and m3Nr
  m = m.Identity();
  m3Nr = 0;

  // Note that we only modify C_[0][1]. At the end of the algorithm,
  // we copy C_[0][1] to C_[1][0]
  bool changed = true;
  while (changed) {
    changed = false;

    if (C_(0, 1) < 0) {
      C_(0, 1) = -C_(0, 1);
      lag_m1(m);
      changed = true;
    }

    if (C_(1, 1) < C_(0, 0)) {
      std::swap(C_(0, 0), C_(1, 1));
      lag_m2(m);
      changed = true;
    }

    if (2 * C_(0, 1) > C_(0, 0)) {
      // The order here matters, don't modify C_(0,1) before using it
      // to calculate C_(1,1).
      C_(1, 1) += C_(0, 0) - 2 * C_(0, 1);
      C_(0, 1) -= C_(0, 0);
      lag_m3(m);
      m3Nr += 1;
      changed = true;
    }
    if (m3Nr >= 1e5) {
      std::cout << nodes[0] << '\n'
                << nodes[1] << '\n'
                << nodes[2] << std::endl;
      throw std::runtime_error("Stuck in lagrange reduction.");
    }
  }
  C_(1, 0) = C_(0, 1);

  // We have had a plastic change if the number of m3 shears have changed
  // Note that pastM3Nr should be update by the mesh when a new loading step
  // has begun, since the minimization algorithm will call this function
  // many times
  plasticChange = pastM3Nr != m3Nr;
}

void TElement::m_updateEnergy() {
  double energyDensity = ContiPotential::energyDensity(
      C_(0, 0), C_(1, 1), C_(0, 1), beta, K, noise);
  energy = energyDensity; // * initArea * F.det();
}

void TElement::m_updateReducedStress() {
  r_s = ContiPotential::stress(C_(0, 0), C_(1, 1), C_(0, 1), beta, K);
}

// Calculate Piola stress tensor and force on each node from current cell
void TElement::m_updatePiolaStress() {
  //  Discontinuous yielding of pristine micro-crystals, page 16/215
  P = 2.0 * F * m * r_s * m.transpose();
}

void TElement::m_updateResolvedShearStress() {
  /**  Discontinuous yielding of pristine micro-crystals (page 216/17)
   *  resolved-shear stress = ∂W/∂α = ∫_Ω P:(∂F/∂α)dx
   *
   * We assume that the shear ∂F/∂α is always
   * ∂F/∂α = [ [0, 1],
   *           [0, 0] ]
   */
  resolvedShearStress = P(0, 1);
}

// Note that each node is part of multiple elements. Therefore, the force must
// be reset after each iteration.
void TElement::applyForcesOnNodes(Mesh &mesh) {
  // TODO explain what is going on here
  for (int i = 0; i < 3; i++) {
    mesh[nodes[i].id]->addForce(P * r[i]);
  }
}

// The functions below are not used in the simulation

double TElement::calculateEnergyDensity(double c11, double c22, double c12) {
  TElement e = TElement();
  e.C = Matrix2d{{c11, c12}, {c12, c22}};
  e.m_lagrangeReduction();
  return ContiPotential::energyDensity(e.C_(0, 0), e.C_(1, 1), e.C_(0, 1), beta,
                                       K);
}

TElement TElement::lagrangeReduction(double c11, double c22, double c12) {
  TElement element = TElement();
  element.C = Matrix2d{{c11, c12}, {c12, c22}};
  element.m_lagrangeReduction();
  return element;
}

void TElement::updatePastM3Nr() { pastM3Nr = m3Nr; }

std::ostream &operator<<(std::ostream &os, const TElement &element) {
  // Save the current format state of the stream
  std::ios_base::fmtflags f(os.flags());

  // Save the current format state of the stream
  std::streamsize prec = os.precision();

  os << std::fixed << std::setprecision(2); // Set precision to 2 decimal places
  os << "Energy: " << element.energy << "\t|";
  for (size_t i = 0; i < element.nodes.size(); ++i) {
    Vector2d pos = element.nodes[i].pos();
    os << "n" << (i + 1) << ": (" << pos[0] << ", " << pos[0] << ")";
    if (i < element.nodes.size() - 1) {
      os << ",\t";
    }
  }
  // Restore the saved precision state
  os.precision(prec);
  os.flags(f);
  return os;
}
