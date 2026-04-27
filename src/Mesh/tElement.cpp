#include "tElement.h"
#include "Eigen/Core"
#include "Mesh/node.h"
#include "Simulation/energyFunctions.h"
#include "mesh.h"
#include "reduction.h"
#include <Eigen/Dense>
#include <Eigen/LU>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <ios>
#include <iostream>
#include <ostream>
#include <sstream>
#include <stdexcept>
using Eigen::Matrix2d;

/*
-1.0, 1.0, 0.0,
-1.0, 0.0, 1.0
*/
Matrix<double, 2, 3> TElement::dN_dxi =
    (Matrix<double, 2, 3>() << -1.0, 1.0, 0.0, -1.0, 0.0, 1.0).finished();
/*
-1.0, 1.0, 0.0,
-1.0, 0.0, 1.0
But transposed.
*/
Matrix<double, 3, 2> TElement::dN_dX_fixed_ref =
    (Matrix<double, 3, 2>() << -1.0, -1.0, 1.0, 0.0, 0.0, 1.0).finished();

TElement::TElement(Mesh &mesh, GhostNode an, GhostNode cn1, GhostNode cn2,
                   int elementIndex, double noise, std::string energyFunction,
                   double bulkModulus)
    : ghostNodes{an, cn1, cn2}, F(Matrix2d::Zero()),
      F_fixed_ref(Matrix2d::Zero()), C(Matrix2d::Zero()),
      C_fixed_ref(Matrix2d::Zero()), C_R(Matrix2d::Zero()),
      C_R_fixed_ref(Matrix2d::Zero()), m(Matrix2d::Zero()),
      m_fixed_ref(Matrix2d::Zero()), S(Matrix2d::Zero()), P(Matrix2d::Zero()),
      sigma(Matrix2d::Zero()), energy(0.0), K(bulkModulus),
      eIndex(elementIndex), noise(noise) {

  if (energyFunction == "contiSquare") {
    beta = -0.25;
  } else if (energyFunction == "contiTriangular") {
    beta = 4;
  } else {
    throw std::invalid_argument("Invalid energy function: " + energyFunction);
  }
  groundStateEnergyDensity = computeGroundStateEnergyDensity();
  // Use a fixed reference area from the mesh to preserve mass even if
  // reconnecting aligns reference nodes into a straight line (area=0).
  initArea = mesh.init_element_area;
  // Add this element to the nodes it is created by
  addElementIndices(mesh, ghostNodes, elementIndex);
  postLoadInit();
}

void TElement::postLoadInit() {

  // Shape functions
  m_update_dX_dxi();

  // Calculate initial area from the reference positions.
  // NOTE: After reconnecting, the reference nodes may end up being aligned in a
  // straight line, resulting in a reference area=0. We therefore keep a fixed
  // reference area from the mesh (mass conservation) and do not recompute it
  // from the reference positions here.
  // initArea = tElementInitialArea(ghostNodes);

  const double det = dX_dxi.determinant();
  const double det_eps = 1e-12 * std::max(1.0, std::abs(2.0 * initArea));
  if (std::abs(det) < det_eps) {
    dxi_dX = Matrix2d::Zero();
  } else {
    dxi_dX = dX_dxi.inverse();
  }
  dN_dX = dN_dxi.transpose() * dxi_dX;

  // Calculate ground state energy density
  groundStateEnergyDensity = calculateEnergyDensity(1, 1, 0);
}
void TElement::updateForces(const Mesh &mesh) {
  // Only compute fixed-reference quantities needed for energy and forces.
  m_updatePosition(mesh);
  m_updateFixedRef();
  m_lagrangeReductionFixedOnly();
  m_updateEnergy();
  m_updateSecondPiolaStress();
  m_updateFirstPiolaStress();
  m_updateForceOnEachNode();
}

void TElement::updateGeometry() {
  updateAngleNode();
  m_update_G();
}

void TElement::updateFull() {
  // Assume positions and fixed-reference quantities are already up to date.
  m_updateDeformationGradientRealOnly();
  m_updateMetricTensorRealOnly();
  m_lagrangeReductionNormalOnly();
  m_updateCauchyStress();
  updateAngles();
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

void TElement::m_updateDeformationGradientRealOnly() {
  m_update_du_dxi();
  F = Matrix2d::Identity();
  F.noalias() += du_dxi * dxi_dX;
}

void TElement::m_updateFixedRef() {
  F_fixed_ref.col(0) = ghostNodes[1].pos - ghostNodes[0].pos;
  F_fixed_ref.col(1) = ghostNodes[2].pos - ghostNodes[0].pos;
  assert(F_fixed_ref.determinant() != 0);

  const double g00 = F_fixed_ref(0, 0);
  const double g01 = F_fixed_ref(0, 1);
  const double g10 = F_fixed_ref(1, 0);
  const double g11 = F_fixed_ref(1, 1);
  C_fixed_ref(0, 0) = g00 * g00 + g10 * g10;
  C_fixed_ref(0, 1) = g00 * g01 + g10 * g11;
  C_fixed_ref(1, 0) = C_fixed_ref(0, 1);
  C_fixed_ref(1, 1) = g01 * g01 + g11 * g11;
}

void TElement::m_updateMetricTensorRealOnly() {
  const double f00 = F(0, 0);
  const double f01 = F(0, 1);
  const double f10 = F(1, 0);
  const double f11 = F(1, 1);
  C(0, 0) = f00 * f00 + f10 * f10;
  C(0, 1) = f00 * f01 + f10 * f11;
  C(1, 0) = C(0, 1);
  C(1, 1) = f01 * f01 + f11 * f11;
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
  // Discontinuous yielding of pristine micro-crystals, page 16/215
  // Sigma = 1/2 (∂Φ/∂C_R + (∂Φ/∂C_R)^T)
  // so it's not actually quite dPhi_dC
  Matrix2d capital_sigma =
      ContiPotential::stress(C_R_fixed_ref(0, 0), C_R_fixed_ref(1, 1),
                             C_R_fixed_ref(0, 1), beta, K, noise);
  // Transform back from lagrange-reudced to un-reduced
  S.noalias() = m_fixed_ref * capital_sigma * m_fixed_ref.transpose();
  S *= 2.0;
}

void TElement::m_updateFirstPiolaStress() {
  // Calculate piola tensor
  P.noalias() = F_fixed_ref * S;
}

void TElement::m_updateCauchyStress() {
  // Using the fixed ref is a trick to avoid colapsed reference states.
  // It also improves element conditioning.
  double J = F_fixed_ref.determinant(); // Jacobian
  sigma.noalias() = (1.0 / J) * P * F_fixed_ref.transpose();
}

void TElement::m_updateForceOnEachNode() {
  // Shlower, more readable code:
  // dPhi_du = P*dN_dX is the energy density gradient
  // Matrix<double, 2, 3> dPhi_du = P * dN_dX_fixed_ref.transpose();
  // Force is the negative of the gradient. Multiply by area since it's a
  // energy DENSITY gradient.
  // dN_dX_fixed_ref rows: [-1,-1], [1,0], [0,1]
  // for (int i = 0; i < 3; i++) {
  //   ghostNodes[i].f = -dPhi_du.col(i) * initArea;
  // }
  // std::cout << '\n';

  // Here we assume a specific dN_dX.
  const Vector2d p0 = P.col(0);
  const Vector2d p1 = P.col(1);
  ghostNodes[0].f = (p0 + p1) * initArea;
  ghostNodes[1].f = -p0 * initArea;
  ghostNodes[2].f = -p1 * initArea;
}

void TElement::m_updatePosition(const Mesh &mesh) {
  // loop through the three nodes in the elements
  for (size_t i = 0; i < 3; i++) {
    // Get the node from the mesh (seperate from the node inside this element)
    const Node *n = mesh[ghostNodes[i].referenceId];
    ghostNodes[i].updatePosition(n, mesh.currentDeformation, mesh.latticeBasis,
                                 mesh.referenceDeformation);
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
  m_lagrangeReductionFixedOnly();
  m_lagrangeReductionNormalOnly();
}

void TElement::m_lagrangeReductionFixedOnly() {
  bool reduced = lagrangeReduction(C_R_fixed_ref, C_fixed_ref, m_fixed_ref,
                                   m3Nr_fix, red_quadrant_fixed);
  if (!reduced) {
    std::cerr << "Lagrange reduction failed for FIXED reference state.\n"
              << "eIndex: " << eIndex << "\n"
              << "F_fixed_ref:\n"
              << F_fixed_ref << "\n"
              << "C_fixed_ref:\n"
              << C_fixed_ref << "\n";
  }
}

void TElement::m_lagrangeReductionNormalOnly() {
  bool reduced = lagrangeReduction(C_R, C, m, m3Nr, red_quadrant);
  if (!reduced) {
    std::cerr << "Lagrange reduction failed for NORMAL state.\n"
              << "eIndex: " << eIndex << "\n"
              << "F:\n"
              << F << "\n"
              << "C:\n"
              << C << "\n";
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

double TElement::computeGroundStateEnergyDensity() const {
  // Assuming calculateEnergyDensity is accessible and noise doesn't matter
  // (or set to zero)
  return ContiPotential::energyDensity(1, 1, 0, beta, K);
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
  return triangleArea(gn[0].ref_pos, gn[1].ref_pos, gn[2].ref_pos);
}

double tElementArea(const GhostNode &A, const GhostNode &B,
                    const GhostNode &C) {
  return triangleArea(A.pos, B.pos, C.pos);
}
