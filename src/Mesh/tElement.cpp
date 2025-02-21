#include "tElement.h"
#include "Eigen/Core"
#include "Simulation/energyFunctions.h"
#include "mesh.h"
#include <Eigen/LU>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <stdexcept>

/*
Shape functions:
N1 = 1 - ξ1 - ξ2
N2 = ξ1
N3 = ξ2

Derivative of shape functions:
-1.0, 1.0, 0.0,
-1.0, 0.0, 1.0
*/
Matrix<double, 2, 3> TElement::dN_dxi =
    (Matrix<double, 2, 3>() << -1.0, 1.0, 0.0, -1.0, 0.0, 1.0).finished();

TElement::TElement(Node n1, Node n2, Node n3, double noise)
    : tElementNodes{n1, n2, n3}, F(Matrix2d::Identity()),
      C(Matrix2d::Identity()), C_(Matrix2d::Identity()),
      m(Matrix2d::Identity()), dPhi_dC_(Matrix2d::Identity()),
      P(Matrix2d::Identity()), noise(noise) {

  // Precompute this constant expression
  m_update_dX_dxi();
  dxi_dX = dX_dxi.inverse();

  // Initialize the adjustment vectors
  for (size_t i = 0; i < dN_dX.size(); ++i) {
    // We transpose dxi_dX because it should be on the other side of dN_dxi,
    // but due to the shapes of the arrays, we do it like this instead.
    // We could have also made dN_dxi be row vectors, and then put things
    // in the correct order. (Well... actually, that might be the explination
    // but i don't think i really understand why we have the transpose here.)
    dN_dX[i] = dxi_dX.transpose() * dN_dxi.col(i);
  }

  // Calculate initial area
  initArea = tElementInitialArea(n1, n2, n3);

  // Calculate ground state energy density
  groundStateEnergyDensity = calculateEnergyDensity(1, 1, 0);
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

  // Calculate force on each node
  m_updateForceOnEachNode();
};

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
Matrix2d TElement::m_update_du_dxi() {
  // ∂u/∂ξ
  du_dxi.col(0) = du(tElementNodes[0], tElementNodes[1]);
  du_dxi.col(1) = du(tElementNodes[0], tElementNodes[2]);
  return du_dxi;
}

// Jacobian with respect to the initial position of the nodes ∂X/∂ξ
// See du_dxi for a similar working out.
Matrix2d TElement::m_update_dX_dxi() {
  // ∂X/∂ξ
  dX_dxi.col(0) = dX(tElementNodes[0], tElementNodes[1]);
  dX_dxi.col(1) = dX(tElementNodes[0], tElementNodes[2]);
  return dX_dxi;
}

void TElement::m_updateReducedStress() {
  dPhi_dC_ =
      ContiPotential::stress(C_(0, 0), C_(1, 1), C_(0, 1), beta, K, noise);
}

// Calculate Piola stress tensor and force on each node from current cell
void TElement::m_updatePiolaStress() {

  // Transform back from lagrange-reudced to un-reduced
  Matrix2d dPhi_dC = m * dPhi_dC_ * m.transpose();
  //  Discontinuous yielding of pristine micro-crystals, page 16/215
  // Calculate piola tensor
  P = 2.0 * F * dPhi_dC;
}

void TElement::m_updateForceOnEachNode() {
  for (int i = 0; i < 3; i++) {
    // Correct up to some constants i guess
    // tElementNodes[i].f = 0.5 * P * dN_dX[i] * dxi_dX.determinant();

    tElementNodes[i].f = P * dN_dX[i];
  }
}

void TElement::m_updatePosition(const Mesh &mesh) {
  // loop through the three nodes in the elements
  for (size_t i = 0; i < 3; i++) {
    // Get the node from the mesh (seperate from the node inside this element)
    const Node *n = mesh[tElementNodes[i].id];
    if (tElementNodes[i].isGhostNode) {
      tElementNodes[i].setPos(
          mesh.makeGhostPos(n->pos(), tElementNodes[i].ghostShift));
    } else {
      tElementNodes[i].setPos(n->pos());
    }
  }
}

void TElement::m_updateDeformationGradiant() {
  // See FEMNotes pdf from Umut

  //                ∂u/∂X = ∂u/∂ξ * ∂ξ/∂X

  // dxi_dX is already computed and is constant
  m_update_du_dxi();
  // Matrix2d du_dX = du_dxi * dxi_dX;
  F = Matrix2d::Identity() + du_dxi * dxi_dX;
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
  mat(0, 1) += mat(0, 0) * n; // Update first row, second column
  mat(1, 1) += mat(1, 0) * n; // Update second row, second column
}

void TElement::m_lagrangeReduction() {
  // Homogeneous nucleation of dislocations as a pattern formation phenomenon -
  // page 5

  // // First check if the previous m still works with the current C
  // if (pastM3Nr != -1 && m_checkIfPreviousReductionWorks()) {
  //   // We only need to remember to update the past M3Nr as this
  //   pastM3Nr = m3Nr;
  //   plasticChange = pastM3Nr != m3Nr;
  //   return;
  // }

  // We start by copying the values from C to the reduced matrix
  // Note that we only modify C_[0][1]. At the end of the algorithm,
  // we copy C_[0][1] to C_[1][0]
  C_ = C;

  // Then reset some values
  m = m.Identity();
  simple_m = simple_m.Identity();
  // We should also reset m and m3Nr
  m1Nr = 0;
  m2Nr = 0;
  m3Nr = 0;

  // Now we keep repeating this loop until C_ does not change
  bool changed = true;
  while (changed) {
    changed = false;

    if (C_(0, 1) < 0) {
      C_(0, 1) = -C_(0, 1);
      lag_m1(m);
      changed = true;
      m1Nr = +1;
    }

    if (C_(1, 1) < C_(0, 0)) {
      std::swap(C_(0, 0), C_(1, 1));
      lag_m2(m);
      changed = true;
      m2Nr += 1;
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
    if (m3Nr >= 10000) {
      std::cout << tElementNodes[0] << '\n'
                << tElementNodes[1] << '\n'
                << tElementNodes[2] << std::endl;
      throw std::runtime_error("Stuck in lagrange reduction.");
    }
  }
  C_(1, 0) = C_(0, 1);

  // We have had a plastic change if the number of m3 shears have changed
  // Note that pastM3Nr should be update by the mesh when a new loading step
  // has begun, since the minimization algorithm will call this function
  // many times
  plasticChange = pastM3Nr != m3Nr;

  // Testing elastic reduction
  // if (m1Nr % 2 == 1) {
  //   lag_m1(simple_m);
  // }
  // if (m2Nr % 2 == 1) {
  //   lag_m2(simple_m);
  // }
}

bool TElement::m_checkIfPreviousReductionWorks() {
  // we can easily check if we even need to do a reduction
  C_ = m.transpose() * C * m;
  // Return true if all the negations of the conditions are true
  return (C_(0, 1) >= 0) &&        //
         (C_(1, 1) >= C_(0, 0)) && //
         (2 * C_(0, 1) <= C_(0, 0));
}

void TElement::m_updateEnergy() {
  double energyDensity = ContiPotential::energyDensity(
      C_(0, 0), C_(1, 1), C_(0, 1), beta, K, noise);
  // Here we we multipy the energy density by the REFERENCE area.
  // Because the Piola tensor is calculated in a lagrangian reference frame, we
  // use the reference area (initArea) instead of the current area (initArea *
  // F.det()).
  energy = (energyDensity - groundStateEnergyDensity) * initArea;
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

std::ostream &operator<<(std::ostream &os, const TElement &element) {
  // Save the current format state of the stream
  std::ios_base::fmtflags f(os.flags());

  // Save the current precision state of the stream
  std::streamsize prec = os.precision();

  os << std::fixed << std::setprecision(2); // Set precision to 2 decimal places
  os << "Energy: " << element.energy << "\t|";
  for (size_t i = 0; i < element.tElementNodes.size(); ++i) {
    Vector2d pos = element.tElementNodes[i].pos();
    os << "n" << (i + 1) << ": (" << pos[0] << ", " << pos[0] << ")";
    if (i < element.tElementNodes.size() - 1) {
      os << ",\t";
    }
  }
  // Restore the saved precision state
  os.precision(prec);
  os.flags(f);
  return os;
}
