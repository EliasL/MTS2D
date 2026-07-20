#include "tElement.h"
#include "Data/logging.h"
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
#include <stdexcept>
#include <string>
#include <vector>
using Eigen::Matrix2d;

/*
-1.0, 1.0, 0.0,
-1.0, 0.0, 1.0
*/
Matrix<double, 2, 3> TElement::dN_dxi =
    (Matrix<double, 2, 3>() << -1.0, 1.0, 0.0, -1.0, 0.0, 1.0).finished();

TElement::TElement(Mesh &mesh, GhostNode n1, GhostNode n2, GhostNode n3,
                   int elementIndex, double noise, std::string energyFunction,
                   double bulkModulus)
    : TElement(n1, n2, n3, energyFunction, bulkModulus,
               "TElement mesh constructor eIndex=" +
                   std::to_string(elementIndex)) {

  eIndex = elementIndex;
  this->noise = noise;

  // Add this element to the nodes it is created by
  addElementIndices(mesh, ghostNodes, elementIndex);

  m_updatePosition(mesh);
  m_updateDeformationGradient();
  m_updateMetricTensor();
  m_lagrangeReduction();
  updateGeometry();
}

TElement::TElement(GhostNode n1, GhostNode n2, GhostNode n3,
                   std::string energyFunction, double bulkModulus)
    : TElement(n1, n2, n3, energyFunction, bulkModulus,
               "TElement constructor") {}

TElement::TElement(GhostNode n1, GhostNode n2, GhostNode n3,
                   std::string energyFunction, double bulkModulus,
                   const std::string &orderContext)
    : F(Matrix2d::Zero()), F_P(Matrix2d::Zero()), F_E(Matrix2d::Zero()),
      C(Matrix2d::Zero()), C_R(Matrix2d::Zero()), G(Matrix2d::Zero()),
      M_l(Matrix2d::Zero()), M_e(Matrix2d::Zero()), S(Matrix2d::Zero()),
      P(Matrix2d::Zero()), sigma(Matrix2d::Zero()), energy(0.0), K(bulkModulus),
      eIndex(-1), noise(1.0) {

  ghostNodes = orderNodes({n1, n2, n3}, orderContext);

  if (energyFunction == "contiSquare") {
    beta = -0.25;
  } else if (energyFunction == "contiTriangular") {
    beta = 4;
  } else {
    throw std::invalid_argument("Invalid energy function: " + energyFunction);
  }

  postLoadInit();

  m_updateDeformationGradient();
  m_updateMetricTensor();
  m_lagrangeReduction();
  updateGeometry();
}

TElement::TElement(std::array<Vector2d, 3> currentNodes,
                   std::array<Vector2d, 3> referenceNodes,
                   std::string energyFunction, double bulkModulus)
    : TElement(GhostNode(currentNodes[0], referenceNodes[0]),
               GhostNode(currentNodes[1], referenceNodes[1]),
               GhostNode(currentNodes[2], referenceNodes[2]), energyFunction,
               bulkModulus) {}

void TElement::postLoadInit() {
  updateAngleNode();
  updateReferenceGeometry();
  groundStateEnergyDensity = computeGroundStateEnergyDensity();
}

void TElement::updateReferenceGeometry() {
  const Matrix2d D_ref =
      referenceEdgeMatrix(ghostNodes, "TElement::updateReferenceGeometry");
  initArea = 0.5 * D_ref.determinant();
  m_update_dN_dX();
}
void TElement::updateForces(const Mesh &mesh) {
  // Only compute fixed-reference quantities needed for energy and forces.
  m_updatePosition(mesh);
  m_updateDeformationGradient();
  m_updateMetricTensor();
  m_lagrangeReduction();
  m_updateEnergy();
  m_updateSecondPiolaStress();
  m_updateFirstPiolaStress();
  m_updateForceOnEachNode();
}

void TElement::updateGeometry() {
  m_update_plastic_elastic_F();
  updateAngleNode();
  m_update_G();
}

void TElement::updateFull() {
  // Assume positions and fixed-reference quantities are already up to date.
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
 */
void TElement::m_update_du_dxi() {
  // ∂u/∂ξ
  du_dxi.col(0) = ghostNodes[1].u - ghostNodes[0].u; // du
  du_dxi.col(1) = ghostNodes[2].u - ghostNodes[0].u; // du
}

void TElement::m_update_dN_dX() {
  const Matrix2d D_ref =
      referenceEdgeMatrix(ghostNodes, "TElement::m_update_dN_dX");

  const double det = D_ref.determinant();
  const double det_eps = 1e-12 * std::max(1.0, std::abs(2.0 * initArea));
  if (std::abs(det) < det_eps) {
    dxi_dX = Matrix2d::Zero();
    throw std::runtime_error("Unexpected zero determinant element!");
  } else {
    dxi_dX = D_ref.inverse();
  }
  dN_dX = dN_dxi.transpose() * dxi_dX;
}

void TElement::m_updateDeformationGradient() {
  m_update_du_dxi();
  F = Matrix2d::Identity();
  F.noalias() += du_dxi * dxi_dX;
}

void TElement::m_updateMetricTensor() {
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
  Vector2d dx12 = ghostNodes[index1].pos - ghostNodes[angleNode].pos;
  Vector2d dx13 = ghostNodes[index2].pos - ghostNodes[angleNode].pos;

  G(0, 0) = dx12.dot(dx12);
  G(0, 1) = dx12.dot(dx13);
  G(1, 0) = G(0, 1);
  G(1, 1) = dx13.dot(dx13);
}

void TElement::m_updateEnergy() {
  double energyDensity = ContiPotential::energyDensity(
      C_R(0, 0), C_R(1, 1), C_R(0, 1), beta, K, noise);
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
      ContiPotential::stress(C_R(0, 0), C_R(1, 1), C_R(0, 1), beta, K, noise);
  // Transform back from lagrange-reudced to un-reduced
  S.noalias() = M_l * capital_sigma * M_l.transpose();
  S *= 2.0;
}

void TElement::m_updateFirstPiolaStress() {
  // Calculate piola tensor
  P.noalias() = F * S;
}

void TElement::m_updateCauchyStress() {
  // Using the fixed ref is a trick to avoid colapsed reference states.
  // It also improves element conditioning.
  double J = F.determinant(); // Jacobian
  sigma.noalias() = (1.0 / J) * P * F.transpose();
}

void TElement::m_updateForceOnEachNode() {
  // Shlower, more readable code:
  // dPhi_du = P*dN_dX is the energy density gradient
  Matrix<double, 2, 3> dPhi_du = P * dN_dX.transpose();
  // Force is the negative of the gradient. Multiply by area since it's a
  // energy DENSITY gradient.
  // dN_dX_fixed_ref rows: [-1,-1], [1,0], [0,1]
  for (int i = 0; i < 3; i++) {
    ghostNodes[i].f = -dPhi_du.col(i) * initArea;
  }
}

void TElement::m_updatePosition(const Mesh &mesh) {
  // loop through the three nodes in the elements
  for (size_t i = 0; i < 3; i++) {
    // Get the node from the mesh (seperate from the node inside this element)
    const Node *n = mesh[ghostNodes[i].referenceId];
    ghostNodes[i].updateCurrentPosition(n, mesh.currentDeformation,
                                        mesh.latticeBasis,
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
  bool reduced = lagrangeReduction(C_R, C, M_e, &M_l, m3Nr, red_quadrant);
  if (!reduced) {
    std::cerr << DebugLog::formatLagrangeReductionFailure(*this) << "\n";
  }
}

void TElement::m_update_plastic_elastic_F() {
  F_E = F * M_e;
  F_P = M_e.inverse();
}

void TElement::updateAngleNode() {
  // Pick the node opposite the longest edge (largest angle in Euclidean
  // triangle). This is faster than computing angles.
  angleNode = findAngleNode(ghostNodes);
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

std::array<const GhostNode *, 2> TElement::getCoNodesByIndex() const {
  if (angleNode < 0 || angleNode >= 3) {
    throw std::runtime_error(
        DebugLog::formatInvalidAngleNode("TElement::getCoNodesByIndex", *this));
  }
  const GhostNode *first = &ghostNodes[(angleNode + 1) % 3];
  const GhostNode *second = &ghostNodes[(angleNode + 2) % 3];
  if (first->referenceId.i < second->referenceId.i) {
    return {first, second};
  }
  return {second, first};
}

std::array<const GhostNode *, 2> TElement::getCoNodesCCW() const {
  if (angleNode < 0 || angleNode >= 3) {
    throw std::runtime_error(
        DebugLog::formatInvalidAngleNode("TElement::getCoNodesCCW", *this));
  }
  int index1 = (angleNode + 1) % 3;
  int index2 = (angleNode + 2) % 3;
  return {&ghostNodes[index1], &ghostNodes[index2]};
}

const GhostNode *TElement::getAngleNode() const {
  if (angleNode < 0 || angleNode >= 3) {
    throw std::runtime_error(
        DebugLog::formatInvalidAngleNode("TElement::getAngleNode", *this));
  }
  const GhostNode *agn = &ghostNodes[angleNode];
  return agn;
}

int TElement::getElementTwin(const Mesh &mesh) const {
  // Mesh::reconnect uses a lookup table for this. This slower local search is
  // kept for analysis/export code and targeted tests.

  // Identify the two nodes to the side of the angle node
  auto coAngleNodes = getCoNodesByIndex();
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
        const TElement &twin = mesh.elements[elementFromNode1];
        auto tCoAngles = twin.getCoNodesByIndex();
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

std::array<const GhostNode, 3> TElement::getAngleCo1Co2Nodes() const {
  auto co = getCoNodesCCW();
  return {ghostNodes[angleNode], *co[0], *co[1]};
}

void TElement::setReferenceElement(const std::array<Vector2d, 3> &refNodes) {
  for (int i = 0; i < 3; ++i) {
    ghostNodes[i].updateReferencePosition(refNodes[i]);
  }
  updateReferenceGeometry();
}

void TElement::setReferenceElementFromCurrentState(const Mesh &mesh) {
  m_updatePosition(mesh);
  setReferenceElement();
}

void TElement::setReferenceElement() {
  // This function initialized the reference element to use the current
  // node positions as the reference. This should not be used during the
  // simulation. For normal reference updates after reconnecting, the old
  // reference should be transformed by F_P.
  setReferenceElement(
      {ghostNodes[0].pos, ghostNodes[1].pos, ghostNodes[2].pos});
}

Vector2d TElement::referenceCentroidShiftToCurrent() const {
  Vector2d currentCentroid = Vector2d::Zero();
  Vector2d referenceCentroid = Vector2d::Zero();
  for (const GhostNode &gn : ghostNodes) {
    currentCentroid += gn.pos;
    referenceCentroid += gn.ref_pos;
  }
  // This is only used for VTU export of the disconnected reference mesh.
  // It translates the reference triangle to the current triangle's centroid
  // without mutating the stored simulation state.
  return (currentCentroid - referenceCentroid) / 3.0;
}

Matrix2d TElement::totalBranch(const Matrix2d &history) const {
  return F_P * history;
}

TElement::EdgeFlipRemeshState TElement::evaluateEdgeFlipRemeshState(
    const std::array<GhostNode, 3> &newGhostNodes,
    const Matrix2d &H_old) const {
  EdgeFlipRemeshState state;
  state.P_old = F_P;
  state.F_old = F;
  state.E_old = state.F_old * state.P_old.inverse();
  state.H_old = H_old;

  state.energy_old = edgeFlipStateEnergy(C_R, beta, K, noise,
                                         groundStateEnergyDensity, initArea);
  const Matrix2d S_old =
      edgeFlipStateSecondPiolaStress(C_R, M_l, beta, K, noise);
  state.sigma_old = edgeFlipStateCauchyStress(state.F_old, S_old);

  auto updateCandidateState = [&](const std::array<GhostNode, 3> &candidate) {
    state.newGhostNodes = prepareEdgeFlipCandidate(
        candidate, "TElement::evaluateEdgeFlipRemeshState eIndex=" +
                       std::to_string(eIndex));
    const Matrix2d D_current = currentEdgeMatrix(state.newGhostNodes);
    const Matrix2d D_ref = referenceEdgeMatrix(
        state.newGhostNodes, "TElement::evaluateEdgeFlipRemeshState");
    state.F_new = D_current * D_ref.inverse();
    state.C_new = state.F_new.transpose() * state.F_new;

    Matrix2d M_e_new = Matrix2d::Identity();
    Matrix2d M_l_new = Matrix2d::Identity();
    int new_m3Nr = 0;
    int new_quadrant = 0;
    plasticReduction(state.C_R_new, state.C_new, M_e_new, &M_l_new, new_m3Nr,
                     new_quadrant, 0.0);
    state.P_new = M_e_new.inverse();
    state.E_new = state.F_new * state.P_new.inverse();
    state.delta_E = state.E_new - state.E_old;
    state.energy_new = edgeFlipStateEnergy(
        state.C_R_new, beta, K, noise, groundStateEnergyDensity,
        tElementInitialArea(state.newGhostNodes));
    const Matrix2d S_new =
        edgeFlipStateSecondPiolaStress(state.C_R_new, M_l_new, beta, K, noise);
    state.sigma_new = edgeFlipStateCauchyStress(state.F_new, S_new);
  };

  updateCandidateState(newGhostNodes);

  state.H_new = state.P_new.inverse() * state.P_old * state.H_old;

  state.elastic_jump = state.delta_E.squaredNorm();
  state.J = state.elastic_jump;
  state.valid = true;

  return state;
}

double TElement::calculateEnergyDensity(double c11, double c22,
                                        double c12) const {
  TElement e = TElement();
  e.C = Matrix2d{{c11, c12}, {c12, c22}};
  e.m_lagrangeReduction();
  return ContiPotential::energyDensity(e.C_R(0, 0), e.C_R(1, 1), e.C_R(0, 1),
                                       beta, K);
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
  os << "Energy: " << element.energy << "|\t";
  for (size_t i = 0; i < element.ghostNodes.size(); ++i) {
    Vector2d pos = element.ghostNodes[i].pos;
    os << "n" << (i + 1) << ": (" << pos[0] << ", " << pos[1] << ")";
    if (i < element.ghostNodes.size() - 1) {
      os << ",\t";
    }
  }
  os << "\nArea: " << element.area() << "|\t";
  for (size_t i = 0; i < element.ghostNodes.size(); ++i) {
    Vector2d refPos = element.ghostNodes[i].ref_pos;
    os << "r_n" << (i + 1) << ": (" << refPos[0] << ", " << refPos[1] << ")";
    if (i < element.ghostNodes.size() - 1) {
      os << ",\t";
    }
  }
  os << '\n';
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

double tElementArea(const std::array<GhostNode, 3> &E) {
  return tElementArea(E[0], E[1], E[2]);
}
double tElementArea(const GhostNode &A, const GhostNode &B,
                    const GhostNode &C) {
  return triangleArea(A.pos, B.pos, C.pos);
}

Matrix2d tElementF(const std::array<GhostNode, 3> &E) {
  return TElement(E[0], E[1], E[2]).F;
}

double polarRotationAngle2D(const std::array<GhostNode, 3> &E) {
  return polarRotationAngle2D(tElementF(E));
}

double polarRotationAngle2D(const Matrix2d &F) {
  return std::atan2(F(1, 0) - F(0, 1), F(0, 0) + F(1, 1));
}

Matrix2d polarRotation2D(const Matrix2d &F) {
  double theta = std::atan2(F(1, 0) - F(0, 1), F(0, 0) + F(1, 1));
  return rotationMatrix(theta);
}

double squareTraceStretch(const Matrix2d &F) {
  // Computes the square of the trace of the stretch in the polar decomposition
  // of F.
  // F=RU return tr(U)^2
  // It would look better to use tr(U-I), but it is slightly more expensive, and
  // it doesn't matter for our purposes.

  if (F.determinant() <= 0) {
    throw std::runtime_error("This function requires positive determinant F!");
  }
  double a = F(0, 0);
  double b = F(0, 1);
  double c = F(1, 0);
  double d = F(1, 1);

  return (a + d) * (a + d) + (b - c) * (b - c);
}

double distanceFromIntegerShear(const Matrix2d &F) {
  Matrix2d F_P;
  return distanceFromIntegerShear(F, F_P);
}

double distanceFromIntegerShear(const Matrix2d &F, Matrix2d &F_P_out) {
  Matrix2d C = F.transpose() * F;
  Matrix2d C_R;
  Matrix2d M_e;
  plasticReduction(C_R, C, M_e);
  F_P_out = M_e.inverse();
  Matrix2d F_E = F * M_e;

  return squareTraceStretch(F_E);
}

std::array<GhostNode, 3> orderNodes(std::array<GhostNode, 3> unorderedNodes) {
  return orderNodes(unorderedNodes, "orderNodes");
}

std::array<GhostNode, 3> orderNodes(std::array<GhostNode, 3> unorderedNodes,
                                    const std::string &context) {
  // Persistent element storage order follows the reference triangle, not the
  // current deformation. This keeps reference-node labels stable across
  // updates, while angleNode is still recomputed from current positions.
  int angleNode = findReferenceAngleNode(unorderedNodes);

  if (angleNode < 0 || angleNode >= 3) {
    throw std::runtime_error("orderNodes: invalid angleNode in " + context +
                             ".");
  }
  int index1 = (angleNode + 1) % 3;
  int index2 = (angleNode + 2) % 3;
  GhostNode a1 = unorderedNodes[angleNode];
  GhostNode g1 = unorderedNodes[index1];
  GhostNode g2 = unorderedNodes[index2];

  const Vector2d edge1 = g1.ref_pos - a1.ref_pos;
  const Vector2d edge2 = g2.ref_pos - a1.ref_pos;
  const double det = edge1.x() * edge2.y() - edge1.y() * edge2.x();
  const double det_eps =
      1e-12 * std::max({1.0, edge1.squaredNorm(), edge2.squaredNorm()});
  if (std::abs(det) < det_eps) {
    throw std::runtime_error("orderNodes: degenerate reference triangle in " +
                             context + ".");
  }

  if (det > 0.0) {
    return {a1, g1, g2};
  }
  return {a1, g2, g1};
}

std::array<GhostNode, 3>
prepareEdgeFlipCandidate(const std::array<GhostNode, 3> &nodes,
                         const std::string &context) {
  std::array<GhostNode, 3> candidate = nodes;
  updateReferencePositions(candidate, closestSquareReferenceNodes(candidate));
  return orderNodes(candidate, context);
}

int findAngleNode(const std::array<GhostNode, 3> &nodes) {
  // Pick the node opposite the longest edge (largest angle in Euclidean
  // triangle). This is faster than computing angles.
  std::array<double, 3> oppositeEdgeLen2 = {0.0, 0.0, 0.0};

  for (int i = 0; i < 3; ++i) {
    const int next = (i + 1) % 3;
    const int prev = (i + 2) % 3;

    const Vector2d edge = nodes[next].pos - nodes[prev].pos;
    oppositeEdgeLen2[i] = edge.squaredNorm();
  }

  int index = 0;
  if (oppositeEdgeLen2[1] > oppositeEdgeLen2[index]) {
    index = 1;
  }
  if (oppositeEdgeLen2[2] > oppositeEdgeLen2[index]) {
    index = 2;
  }

  return index;
}

int findReferenceAngleNode(const std::array<GhostNode, 3> &nodes) {
  double largestLength = -1.0;
  int index = -1;

  for (int i = 0; i < 3; ++i) {
    const int next = (i + 1) % 3;
    const int prev = (i + 2) % 3;

    const Vector2d edge = nodes[next].ref_pos - nodes[prev].ref_pos;
    const double len2 = edge.squaredNorm();

    if (len2 > largestLength) {
      largestLength = len2;
      index = i;
    }
  }
  return index;
}

Matrix2d currentEdgeMatrix(const std::array<GhostNode, 3> &nodes) {
  Matrix2d D_current;
  D_current.col(0) = nodes[1].pos - nodes[0].pos;
  D_current.col(1) = nodes[2].pos - nodes[0].pos;
  return D_current;
}

Matrix2d referenceEdgeMatrix(const std::array<GhostNode, 3> &nodes,
                             const std::string &context) {
  Matrix2d D_ref;
  D_ref.col(0) = nodes[1].ref_pos - nodes[0].ref_pos;
  D_ref.col(1) = nodes[2].ref_pos - nodes[0].ref_pos;
  const double det = D_ref.determinant();
  if (std::abs(det) < 1e-12) {
    throw std::runtime_error(context +
                             ": unexpected degenerate reference element.");
  }
  if (det < 0.0) {
    throw std::runtime_error(context +
                             ": reference nodes must stay in counterclockwise "
                             "order.");
  }
  return D_ref;
}

double edgeFlipStateEnergy(const Matrix2d &C_R, double beta, double K,
                           double noise, double groundStateEnergyDensity,
                           double initArea) {
  const double energyDensity = ContiPotential::energyDensity(
      C_R(0, 0), C_R(1, 1), C_R(0, 1), beta, K, noise);
  return (energyDensity - groundStateEnergyDensity) * initArea;
}

Matrix2d edgeFlipStateSecondPiolaStress(const Matrix2d &C_R,
                                        const Matrix2d &M_l, double beta,
                                        double K, double noise) {
  Matrix2d capital_sigma =
      ContiPotential::stress(C_R(0, 0), C_R(1, 1), C_R(0, 1), beta, K, noise);
  Matrix2d S = M_l * capital_sigma * M_l.transpose();
  S *= 2.0;
  return S;
}

Matrix2d edgeFlipStateCauchyStress(const Matrix2d &F, const Matrix2d &S) {
  const double J = F.determinant();
  return (1.0 / J) * (F * S) * F.transpose();
}

std::array<Vector2d, 3>
closestSquareReferenceNodes(const std::array<GhostNode, 3> &nodes) {
  const int currentAngleNode = findAngleNode(nodes);
  if (currentAngleNode < 0 || currentAngleNode >= 3) {
    throw std::runtime_error(
        "closestSquareReferenceNodes: invalid current angle node.");
  }

  const int co1Index = (currentAngleNode + 1) % 3;
  const int co2Index = (currentAngleNode + 2) % 3;
  const Vector2d u = nodes[co1Index].pos - nodes[currentAngleNode].pos;
  const Vector2d v = nodes[co2Index].pos - nodes[currentAngleNode].pos;

  struct Candidate {
    Vector2d angleCorner;
    Vector2d adjacentCorner1;
    Vector2d adjacentCorner2;
  };

  const std::array<Candidate, 4> candidates = {{
      {Vector2d(-0.5, 0.5), Vector2d(0.5, 0.5), Vector2d(-0.5, -0.5)},  // ADB
      {Vector2d(-0.5, -0.5), Vector2d(0.5, -0.5), Vector2d(-0.5, 0.5)}, // ADC
      {Vector2d(0.5, 0.5), Vector2d(-0.5, 0.5), Vector2d(0.5, -0.5)},   // BCA
      {Vector2d(0.5, -0.5), Vector2d(-0.5, -0.5), Vector2d(0.5, 0.5)},  // BCD
  }};

  std::array<Vector2d, 3> bestReferenceNodes;
  const Candidate &firstCandidate = candidates[0];
  {
    const Vector2d leg1 =
        firstCandidate.adjacentCorner1 - firstCandidate.angleCorner;
    const Vector2d leg2 =
        firstCandidate.adjacentCorner2 - firstCandidate.angleCorner;
    const double score12 = u.dot(leg1) + v.dot(leg2);
    const double score21 = u.dot(leg2) + v.dot(leg1);
    bestReferenceNodes[currentAngleNode] = firstCandidate.angleCorner;
    if (score12 >= score21) {
      bestReferenceNodes[co1Index] = firstCandidate.adjacentCorner1;
      bestReferenceNodes[co2Index] = firstCandidate.adjacentCorner2;
    } else {
      bestReferenceNodes[co1Index] = firstCandidate.adjacentCorner2;
      bestReferenceNodes[co2Index] = firstCandidate.adjacentCorner1;
    }
  }
  double bestScore = std::max(
      u.dot(firstCandidate.adjacentCorner1 - firstCandidate.angleCorner) +
          v.dot(firstCandidate.adjacentCorner2 - firstCandidate.angleCorner),
      u.dot(firstCandidate.adjacentCorner2 - firstCandidate.angleCorner) +
          v.dot(firstCandidate.adjacentCorner1 - firstCandidate.angleCorner));

  for (size_t candidateIndex = 1; candidateIndex < candidates.size();
       ++candidateIndex) {
    const Candidate &candidate = candidates[candidateIndex];
    const Vector2d leg1 = candidate.adjacentCorner1 - candidate.angleCorner;
    const Vector2d leg2 = candidate.adjacentCorner2 - candidate.angleCorner;

    const double score12 = u.dot(leg1) + v.dot(leg2);
    const double score21 = u.dot(leg2) + v.dot(leg1);

    std::array<Vector2d, 3> candidateReferenceNodes;
    candidateReferenceNodes[currentAngleNode] = candidate.angleCorner;

    const double score = std::max(score12, score21);
    if (score12 >= score21) {
      candidateReferenceNodes[co1Index] = candidate.adjacentCorner1;
      candidateReferenceNodes[co2Index] = candidate.adjacentCorner2;
    } else {
      candidateReferenceNodes[co1Index] = candidate.adjacentCorner2;
      candidateReferenceNodes[co2Index] = candidate.adjacentCorner1;
    }

    if (score > bestScore) {
      bestScore = score;
      bestReferenceNodes = candidateReferenceNodes;
    }
  }

  return bestReferenceNodes;
}

void updateReferencePositions(std::array<GhostNode, 3> &nodes,
                              const std::array<Vector2d, 3> &refNodes) {
  for (size_t i = 0; i < nodes.size(); ++i) {
    nodes[i].updateReferencePosition(refNodes[i]);
  }
}
