#include "tElement.h"
#include "Eigen/Core"
#include "Mesh/node.h"
#include "Simulation/energyFunctions.h"
#include "mesh.h"
#include <Eigen/Dense>
#include <Eigen/LU>
#include <array>
#include <cassert>
#include <iomanip>
#include <ios>
#include <iostream>
#include <ostream>
#include <stdexcept>
using Eigen::Matrix2d;
// Define and initialize static members
double TElement::groundStateEnergyDensity =
    TElement::computeGroundStateEnergyDensity();

/*
-1.0, 1.0, 0.0,
-1.0, 0.0, 1.0
*/
Matrix<double, 2, 3> TElement::dN_dxi =
    (Matrix<double, 2, 3>() << -1.0, 1.0, 0.0, -1.0, 0.0, 1.0).finished();

Matrix<double, 3, 2> TElement::dN_dX_fixed_ref =
    (Matrix<double, 3, 2>() << -1.0, -1.0, 1.0, 0.0, 0.0, 1.0).finished();

TElement::TElement(Mesh &mesh, GhostNode an, GhostNode cn1, GhostNode cn2,
                   int elementIndex, double noise)
    : ghostNodes{an, cn1, cn2}, F(Matrix2d::Zero()),
      F_fixed_ref(Matrix2d::Zero()), C(Matrix2d::Zero()),
      C_fixed_ref(Matrix2d::Zero()), C_R(Matrix2d::Zero()),
      C_R_fixed_ref(Matrix2d::Zero()), m(Matrix2d::Zero()),
      m_fixed_ref(Matrix2d::Zero()), S(Matrix2d::Zero()), P(Matrix2d::Zero()),
      sigma(Matrix2d::Zero()), energy(0.0), P_xy(0.0), sigma_xy(0.0),
      eIndex(elementIndex), noise(noise) {

  // Add this element to the nodes it is created by
  addElementIndices(mesh, ghostNodes, elementIndex);
  postLoadInit();
}

void TElement::postLoadInit() {

  // Shape functions
  m_update_dX_dxi();

  if (dX_dxi.determinant() == 0) {
    dxi_dX = Matrix2d::Zero();
    initArea = 0.5; // Note that this is an assumption we make
    // When we use a fixed reference, we require the initial area to be 0.5
    // for correct energy values.
  } else {

    dxi_dX = dX_dxi.inverse();
    // Calculate initial area
    initArea = tElementInitialArea(ghostNodes);
  }
  dN_dX = dN_dxi.transpose() * dxi_dX;

  assert(initArea = 0.5);

  // Calculate ground state energy density
  groundStateEnergyDensity = calculateEnergyDensity(1, 1, 0);
}
void TElement::update(const Mesh &mesh) {
  // The order here is very important.

  // Update nodes in element
  m_updatePosition(mesh);

  // Find and update largest angle
  updateAngleNode();

  // Calculates F
  m_updateDeformationGradiant();

  // Calculates C (and G)
  m_updateMetricTensor();

  // Calculate C_ and m
  m_lagrangeReduction();

  // Calculate energy
  m_updateEnergy();

  // Calculate stress
  m_updateSecondPiolaStress();
  m_updateFirstPiolaStress();
  m_updateCauchyStress();

  // Calculate resolved shear stress
  m_updateShearStress();

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
  du_dxi.col(0) = du(ghostNodes[0], ghostNodes[1]);
  du_dxi.col(1) = du(ghostNodes[0], ghostNodes[2]);
  return du_dxi;
}

// Jacobian with respect to the initial position of the nodes ∂X/∂ξ
// See du_dxi for a similar working out.
Matrix2d TElement::m_update_dX_dxi() {
  // ∂X/∂ξ
  dX_dxi.col(0) = dX(ghostNodes[0], ghostNodes[1]);
  dX_dxi.col(1) = dX(ghostNodes[0], ghostNodes[2]);
  return dX_dxi;
}

void TElement::m_updateDeformationGradiant() {

  // dxi_dX is already computed and is constant
  m_update_du_dxi();
  // Matrix2d du_dX = du_dxi * dxi_dX;
  F = Matrix2d::Identity() + du_dxi * dxi_dX;

  // F might not be invertable, so to calculate the force,
  // we use a fixed reference
  Matrix<double, 2, 3> x;
  x.col(0) = ghostNodes[0].pos;
  x.col(1) = ghostNodes[1].pos;
  x.col(2) = ghostNodes[2].pos;

  F_fixed_ref = x * dN_dX_fixed_ref;

  assert(F_fixed_ref.determinant() != 0);
}

// Provices a metric tensor for the triangle
void TElement::m_updateMetricTensor() {
  // Discontinuous yielding of pristine micro-crystals - page 8/207
  C = F.transpose() * F;
  C_fixed_ref = F_fixed_ref.transpose() * F_fixed_ref;
  assert(F_fixed_ref.determinant() != 0);

  m_update_G();
}

void TElement::m_update_G() {
  // G = [dx_12.dx_12, dx_12.dx_13;
  //      dx_13.dx_12, dx_13.dx_13]
  // We first calculate the two vectors
  int index1 = (angleNode + 1) % 3;
  int index2 = (angleNode + 2) % 3;
  // We always take the vectors to be from the angleNode to the other two nodes
  Vector2d dx12 = dx(ghostNodes[angleNode], ghostNodes[index1]);
  Vector2d dx13 = dx(ghostNodes[angleNode], ghostNodes[index2]);

  G(0, 0) = dx12.dot(dx12);
  G(0, 1) = dx12.dot(dx13);
  G(1, 0) = G(0, 1);
  G(1, 1) = dx13.dot(dx13);
}

void TElement::m_updateEnergy() {
  double energyDensity =
      ContiPotential::energyDensity(C_R_fixed_ref(0, 0), C_R_fixed_ref(1, 1),
                                    C_R_fixed_ref(0, 1), beta, K, noise);
  // Here we we multipy the energy density by the REFERENCE (initial) area.
  // Because the Piola tensor is calculated in a lagrangian reference frame, we
  // use the reference area (initArea) instead of the current area (initArea *
  // F.det()).
  energy = (energyDensity - groundStateEnergyDensity) * initArea;
}

void TElement::m_updateSecondPiolaStress() {
  // Sigma = 1/2 (∂Φ/∂C_R + (∂Φ/∂C_R)^T)
  Matrix2d capital_sigma =
      ContiPotential::stress(C_R_fixed_ref(0, 0), C_R_fixed_ref(1, 1),
                             C_R_fixed_ref(0, 1), beta, K, noise);
  // Transform back from lagrange-reudced to un-reduced
  // so it's not actually quite dPhi_dC
  S = m_fixed_ref * capital_sigma * m_fixed_ref.transpose();
}

// Calculate Piola stress tensor and force on each node from current cell
void TElement::m_updateFirstPiolaStress() {
  //  Discontinuous yielding of pristine micro-crystals, page 16/215
  // Calculate piola tensor
  P = 2.0 * F_fixed_ref * S;
}
// Calculate Piola stress tensor and force on each node from current cell
void TElement::m_updateCauchyStress() {
  //  Discontinuous yielding of pristine micro-crystals, page 16/215
  // Calculate piola tensor
  double J = F_fixed_ref.determinant(); // Jacobian
  sigma = (1.0 / J) * P * F_fixed_ref.transpose();
}

void TElement::m_updateForceOnEachNode() {
  // dPhi_du = P*dN_dX is the energy density gradient
  Matrix<double, 2, 3> dPhi_du = P * dN_dX_fixed_ref.transpose();
  for (int i = 0; i < 3; i++) {
    // Force is the negative of the gradient
    // Multipy by area since it's a energy DENSITY gradient
    ghostNodes[i].f = -dPhi_du.col(i) * initArea;
    // std::cout << "E:" << eIndex << " En: " << i
    //           << " Rn: " << ghostNodes[i].referenceId.i
    //           << " f: " << ghostNodes[i].f.transpose() << "\n";
  }
  // std::cout << '\n';
}

void TElement::m_updatePosition(const Mesh &mesh) {
  // loop through the three nodes in the elements
  for (size_t i = 0; i < 3; i++) {
    // Get the node from the mesh (seperate from the node inside this element)
    const Node *n = mesh[ghostNodes[i].referenceId];
    ghostNodes[i].updatePosition(n, mesh.currentDeformation, mesh.a);
  }

  // In order to make it obvious if we forget to update the angles, we give them
  // some invalid values here.
  largestAngle = -1.0;
  smallestAngle = -1.0;
  angleNode = -1;
  // We always update the angle node, but not always the angles. So it is easy
  // to forget.
}

void TElement::m_lagrangeReduction() {
  // We need some extra code here to handle error cases. The Fail State tells
  // us where the reduction failed, so we can print the relevant F and C for
  // debugging.
  enum class FailStage { None, Fixed, Normal };
  FailStage failed = FailStage::None;

  bool reduced = lagrangeReduction(C_R_fixed_ref, C_fixed_ref, m_fixed_ref,
                                   m1Nr, m2Nr, m3Nr);
  if (!reduced) {
    failed = FailStage::Fixed;
  } else {
    // Reduce the "normal" state only after fixed succeeds

    // The order is important. The fixed reduction should always be done first.
    // The lagrange reduction overwrites m1Nr, m2Nr, m3Nr and we are
    // mostly interested in the values after the normal reduction.
    reduced = lagrangeReduction(C_R, C, m, m1Nr, m2Nr, m3Nr);
    if (!reduced)
      failed = FailStage::Normal;
  }

  if (failed != FailStage::None) {
    // Save/restore stream state
    std::ios oldState(nullptr);
    oldState.copyfmt(std::cout);
    std::cout << std::defaultfloat << std::setprecision(2);

    // Always print iteration counts
    std::cout << "Lagrange Reduction Iteration Counts:\n"
              << "  m1Nr: " << m1Nr << "\n"
              << "  m2Nr: " << m2Nr << "\n"
              << "  m3Nr: " << m3Nr << "\n\n";

    // Print only the relevant F/C depending on where it failed
    if (failed == FailStage::Fixed) {
      std::cerr << "Lagrange reduction failed for FIXED reference state.\n";
      std::cout << "Fixed_ref Deformation Gradient F_fixed:\n"
                << F_fixed_ref << "\n\n"
                << "Fixed_ref Metric Tensor C_fixed_ref:\n"
                << C_fixed_ref << "\n\n";
    } else { // FailStage::Normal
      std::cerr << "Lagrange reduction failed for NORMAL state.\n";
      std::cout << "Deformation Gradient F:\n"
                << F << "\n\n"
                << "Metric Tensor C:\n"
                << C << "\n\n";
    }

    // Always print ghost node positions
    std::cout << "i0: " << ghostNodes[0] << "\n"
              << "i1: " << ghostNodes[1] << "\n"
              << "i2: " << ghostNodes[2] << "\n";

    std::cout.copyfmt(oldState);
    throw std::runtime_error("Stuck in lagrange reduction.\n");
  }
}

void TElement::updateAngleNode() {
  // Pick the node opposite the longest edge (largest angle in Euclidean
  // triangle). This is faster than computing angles.
  double largestLength = -1.0;

  for (int i = 0; i < 3; ++i) {
    const int next = (i + 1) % 3;
    const int prev = (i + 2) % 3;

    const Vector2d edge = ghostNodes[next].pos - ghostNodes[prev].pos;
    const double len2 = edge.squaredNorm();

    // prefer strictly longer
    if (len2 > largestLength) {
      largestLength = len2;
      angleNode = i;
    }
  }
}

void TElement::updateAngles() {
  // find the largest and smallest angle in the triangle
  double maxAngle = 0.0;
  int largestAngleIndex = 0;
  double minAngle = 180.0;

  // Compute and compare all three angles
  for (int i = 0; i < 3; i++) {
    int next = (i + 1) % 3;
    int prev = (i + 2) % 3;

    // Compute vectors from the current vertex to adjacent vertices
    Vector2d v1 = ghostNodes[next].pos - ghostNodes[i].pos;
    Vector2d v2 = ghostNodes[prev].pos - ghostNodes[i].pos;

    // Compute angle using dot product and vector magnitudes
    double magnitudeProduct = v1.norm() * v2.norm();

    // Avoid division by zero
    if (magnitudeProduct > 1e-10) {
      double cosAngle = std::clamp(v1.dot(v2) / magnitudeProduct, -1.0, 1.0);
      // TODO acos is slow.
      double angle = std::acos(cosAngle);

      // Convert to degrees
      angle *= 180.0 / M_PI;

      // Track the largest angle
      if (angle > maxAngle) {
        maxAngle = angle;
        largestAngleIndex = i;
      }
      // Track the smallest angle
      if (angle < minAngle) {
        minAngle = angle;
      }
    }
  }

  // Store the results
  angleNode = largestAngleIndex;
  largestAngle = maxAngle;
  smallestAngle = minAngle;
}

std::array<const GhostNode *, 2> TElement::getCoAngleNodes() const {
  int index1 = (angleNode + 1) % 3;
  int index2 = (angleNode + 2) % 3;
  const GhostNode *g1 = &ghostNodes[index1];
  const GhostNode *g2 = &ghostNodes[index2];
  // In order to compare angle nodes, we always sort in a consitent order
  if (g1->referenceId.i < g2->referenceId.i) {
    return {g1, g2};
  } else {
    return {g2, g1};
  }
}

GhostNode *TElement::getAngleNode() {
  GhostNode *agn = &ghostNodes[angleNode];
  return agn;
}

int TElement::getElementTwin(const Mesh &mesh) const {
  // TODO Create an edge lookup table in the mesh, and use that instead
  // Note that it needs to be updated in the case of a reconnect.
  // I check for reconnecting just more seldomly, and now this function doesn't
  // affect the performance so much, so making it faster is not so important.

  // Identify the two nodes to the side of the angle node
  auto coAngleNodes = getCoAngleNodes();
  const Node *n1 = mesh[coAngleNodes[0]->referenceId];
  const Node *n2 = mesh[coAngleNodes[1]->referenceId];

  // Find all elements that are common for both nodes and not this element
  for (int elementFromNode1 : n1->connectedElements) {
    // Skip the current element or end of valid elements
    if (elementFromNode1 == eIndex || elementFromNode1 == -1) {
      continue;
    }

    for (int elementFromNode2 : n2->connectedElements) {
      // Skip the current element or end of valid elements
      if (elementFromNode2 == eIndex || elementFromNode2 == -1) {
        continue;
      }

      // If we find an element that contains both nodes (and that is not this
      // element)
      if (elementFromNode1 == elementFromNode2) {
        // We now check that the two nodes they share are coAngleNodes
        const TElement &twin = mesh.elements[elementFromNode1];
        auto tCoAngles = twin.getCoAngleNodes();
        if ((tCoAngles[0]->referenceId == coAngleNodes[0]->referenceId) &&
            (tCoAngles[1]->referenceId == coAngleNodes[1]->referenceId)) {
          return elementFromNode1;
        }
      }
    }
  }
  // No match found
  return -1;
}

void TElement::m_updateShearStress() {
  // shear component of the first Piola-Kirchhoff stress
  P_xy = P(0, 1);

  // shear component of the Cauchy stress
  sigma_xy = sigma(0, 1);
}

// The functions below are not used in the simulation
double TElement::calculateEnergyDensity(double c11, double c22,
                                        double c12) const {
  TElement e = TElement();
  e.C_fixed_ref = Matrix2d{{c11, c12}, {c12, c22}};
  e.C = Matrix2d::Identity(); // Used to avoid error in lagrangeReduction
  e.m_lagrangeReduction();
  return ContiPotential::energyDensity(e.C_R_fixed_ref(0, 0),
                                       e.C_R_fixed_ref(1, 1),
                                       e.C_R_fixed_ref(0, 1), beta, K);
}

TElement TElement::reduce_element(double c11, double c22, double c12) {
  TElement element = TElement();
  element.C = Matrix2d{{c11, c12}, {c12, c22}};
  element.m_lagrangeReduction();
  return element;
}

Vector2d TElement::getCom() {
  return (ghostNodes[0].pos + ghostNodes[1].pos + ghostNodes[2].pos) / 3;
}
double TElement::area() const {
  return tElementArea(ghostNodes[0], ghostNodes[1], ghostNodes[2]);
}

//------- Non TElement functions

// --- tiny helpers (inlined, no temporaries) ---
inline void mul_m1(Matrix2d &m) {
  // Right-multiply by diag(1,-1): col1 -> -col1
  m(0, 1) = -m(0, 1);
  m(1, 1) = -m(1, 1);
}

inline void mul_m2(Matrix2d &m) {
  // Right-multiply by [[0,1],[1,0]]: swap columns
  std::swap(m(0, 0), m(0, 1));
  std::swap(m(1, 0), m(1, 1));
}

inline void mul_m3(Matrix2d &m, int n = -1) {
  // Right-multiply by [[1,n],[0,1]]: col1 += n*col0
  m(0, 1) += n * m(0, 0);
  m(1, 1) += n * m(1, 0);
}

/**
 * Lagrange reduction for a 2x2 symmetric metric in-place.
 * All arguments are passed by reference; no dynamic allocations.
 *
 * @param C_R       Output reduced metric (also used as the working matrix).
 * @param C_in      Input metric to start from (read-only).
 * @param m         Accumulated integer-transform matrix (updated in-place).
 * @param m1Nr      Counter for op1 (+= by ref).
 * @param m2Nr      Counter for op2 (+= by ref).
 * @param m3Nr      Counter for op3 (+= by ref).
 * @param maxLoops  Safety cap on iterations.
 * @param eps       Tolerance to avoid numerical chattering.
 * @return          true if converged before maxLoops; false otherwise.
 */
bool lagrangeReduction(Matrix2d &C_R,        // work/output: [[a,b],[b,c]]
                       const Matrix2d &C_in, // input metric
                       Matrix2d &m,          // accumulated transform
                       int &m1Nr, int &m2Nr, int &m3Nr, int maxLoops) {
  C_R = C_in;
  m.setIdentity();
  m1Nr = m2Nr = m3Nr = 0;
  double &a = C_R(0, 0);
  double &b = C_R(0, 1);
  double &c = C_R(1, 1);
  assert(a > 0 && c > 0);

  for (int iter = 0; iter < maxLoops; ++iter) {
    bool changed = false;

    // 1) make off-diagonal non-negative
    if (std::signbit(b)) {
      b = -b;
      mul_m1(m);
      ++m1Nr;
      changed = true;
    }

    // 2) enforce a <= c
    if (c < a) {
      std::swap(a, c);
      mul_m2(m);
      ++m2Nr;
      changed = true;
    }

    // 3) single-shot shear: choose n so that -a/2 <= b + n a < a/2
    // Traditionally, n is always -1, but here, we avoid repeated m3 steps
    // Instead of looping several times and doing b += -a each time,
    // we calculate how many times we would need to do that, and do it all at
    // once.
    const double two_b = b + b;
    if (two_b > a) {
      // b >= 0 here, so this is a fast nearest-integer round (ties away from 0)
      const int n = -static_cast<int>(b / a + 0.5);
      if (n != 0) {
        const double n_a = n * a;
        c += n * (two_b + n_a); // uses old b
        b += n_a;
        mul_m3(m, n);
        m3Nr += (n < 0) ? -n : n;
        changed = true;
      }
    }

    if (!changed) {
      C_R(1, 0) = C_R(0, 1); // enforce symmetry once at the end
      return true;
    }
  }

  C_R(1, 0) = C_R(0, 1);
  return false; // hit loop cap (shouldn't happen in practice)
}

void addElementIndices(Mesh &mesh, const std::array<GhostNode, 3> &nodeList,
                       int elementIndex) {
  for (size_t i = 0; i < nodeList.size(); ++i) {
    const GhostNode &gn = nodeList[i];

    // Reference to the current count
    Node *node = mesh[gn.referenceId];
    int &count = node->elementCount;
    // Ensure we don't exceed the array size
    if (count < MAX_ELEMENTS_PER_NODE) {
      node->connectedElements[count] = elementIndex;
      node->nodeIndexInElement[count] = static_cast<int>(i);
      node->connectedGhostNodes[count] = nullptr;
      ++count; // Increment the count for the node
    } else {
      // Handle overflow (e.g., log an error or take other measures)
      throw std::overflow_error("Element index overflow for node " +
                                std::to_string(gn.referenceId.i));

      std::cerr << "Error: Too many elements for node " << gn.referenceId
                << std::endl;
    }
  }
};

std::ostream &operator<<(std::ostream &os, const TElement &element) {
  // Save the current format state of the stream
  std::ios_base::fmtflags f(os.flags());

  // Save the current precision state of the stream
  std::streamsize prec = os.precision();

  os << std::fixed << std::setprecision(2); // Set precision to 2 decimal places
  os << "Energy: " << element.energy << "\t|";
  for (size_t i = 0; i < element.ghostNodes.size(); ++i) {
    Vector2d pos = element.ghostNodes[i].pos;
    os << "n" << (i + 1) << ": (" << pos[0] << ", " << pos[1] << ")";
    if (i < element.ghostNodes.size() - 1) {
      os << ",\t";
    }
  }
  // Restore the saved precision state
  os.precision(prec);
  os.flags(f);
  return os;
}

double triangleArea(Vector2d posA, Vector2d posB, Vector2d posC) {
  return 0.5 * std::abs(posA[0] * (posB[1] - posC[1]) +
                        posB[0] * (posC[1] - posA[1]) +
                        posC[0] * (posA[1] - posB[1]));
}
double tElementInitialArea(const std::array<GhostNode, 3> &gn) {
  return triangleArea(gn[0].init_pos, gn[1].init_pos, gn[2].init_pos);
}

double tElementArea(const GhostNode &A, const GhostNode &B,
                    const GhostNode &C) {
  return triangleArea(A.pos, B.pos, C.pos);
}
