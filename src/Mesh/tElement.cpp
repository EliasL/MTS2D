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
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using Eigen::Matrix2d;

std::array<GhostNode, 3> orderNodes(std::array<GhostNode, 3> unorderedNodes);
std::array<GhostNode, 3> orderNodes(std::array<GhostNode, 3> unorderedNodes,
                                    const std::string &context);
int findAngleNode(std::array<GhostNode, 3> nodes);
Matrix2d rotationMatrix2D(double theta);
double wrappedAngleDifference(double theta, double thetaReference);

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
  const Matrix2d D_R = referenceEdgeMatrix();
  initArea = 0.5 * D_R.determinant();
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

// Jacobian with respect to the initial position of the nodes ∂X/∂ξ
// See du_dxi for a similar working out.
// Note that we use the special reference element positions. Not the true
// reference positions of the "real" nodes.
void TElement::m_update_dX_dxi() {
  // ∂X/∂ξ
  dX_dxi.col(0) = ghostNodes[1].ref_pos - ghostNodes[0].ref_pos; // dX
  dX_dxi.col(1) = ghostNodes[2].ref_pos - ghostNodes[0].ref_pos; // dX
}

void TElement::m_update_dN_dX() {
  // Shape functions
  m_update_dX_dxi();

  const double det = dX_dxi.determinant();
  const double det_eps = 1e-12 * std::max(1.0, std::abs(2.0 * initArea));
  if (std::abs(det) < det_eps) {
    dxi_dX = Matrix2d::Zero();
    throw std::runtime_error("Unexpected zero determinant element!");
  } else {
    dxi_dX = dX_dxi.inverse();
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
    std::cerr << "Lagrange reduction failed for FIXED reference state.\n"
              << "eIndex: " << eIndex << "\n"
              << "F_fixed_ref:\n"
              << F << "\n"
              << "C_fixed_ref:\n"
              << C << "\n";
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

std::array<const GhostNode *, 2> TElement::getCoAngleNodes() const {
  if (angleNode < 0 || angleNode >= 3) {
    throw std::runtime_error("TElement::getCoAngleNodes: invalid angleNode.");
  }
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

const GhostNode *TElement::getAngleNode() const {
  if (angleNode < 0 || angleNode >= 3) {
    throw std::runtime_error("TElement::getAngleNode: invalid angleNode.");
  }
  const GhostNode *agn = &ghostNodes[angleNode];
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

std::array<const GhostNode, 3> TElement::getAngleCo1Co2Nodes() const {
  auto co = getCoAngleNodes();
  return {ghostNodes[angleNode], *co[0], *co[1]};
}

void TElement::setReferenceElement(const std::array<Vector2d, 3> &refNodes) {
  for (int i = 0; i < 3; ++i) {
    ghostNodes[i].updateReferencePosition(refNodes[i]);
  }
  updateReferenceGeometry();
}

void TElement::setReferenceElement(
    const std::array<const GhostNode, 3> &refNodes) {
  std::array<bool, 3> refNodeUsed = {false, false, false};
  std::array<bool, 3> currentNodeMatched = {false, false, false};

  // First match identical ghost-node ids so shared vertices keep their
  // reference positions across an edge flip.
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (refNodeUsed[j]) {
        continue;
      }
      if (ghostNodes[i].id != refNodes[j].id) {
        continue;
      }
      ghostNodes[i].updateReferencePosition(refNodes[j].ref_pos);
      refNodeUsed[j] = true;
      currentNodeMatched[i] = true;
      break;
    }
  }

  // Fall back to the remaining input order for vertices that do not have a
  // matching ghost-node id in the source reference element.
  for (int i = 0; i < 3; ++i) {
    if (currentNodeMatched[i]) {
      continue;
    }
    bool assigned = false;
    for (int j = 0; j < 3; ++j) {
      if (refNodeUsed[j]) {
        continue;
      }
      ghostNodes[i].updateReferencePosition(refNodes[j].ref_pos);
      refNodeUsed[j] = true;
      assigned = true;
      break;
    }
    if (!assigned) {
      throw std::runtime_error(
          "TElement::setReferenceElement: failed to assign all reference "
          "nodes.");
    }
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

void TElement::deformReferenceElement(Matrix2d F, Vector2d oldAnchor,
                                      Vector2d newAnchor) {
  for (int i = 0; i < 3; i++) {
    ghostNodes[i].transformReferencePosition(F, oldAnchor, newAnchor);
  }
  updateReferenceGeometry();
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

Matrix2d TElement::referenceEdgeMatrix() const {
  Matrix2d D_R;
  D_R.col(0) = ghostNodes[1].ref_pos - ghostNodes[0].ref_pos;
  D_R.col(1) = ghostNodes[2].ref_pos - ghostNodes[0].ref_pos;
  const double det = D_R.determinant();
  if (std::abs(det) < 1e-12) {
    throw std::runtime_error(
        "TElement::referenceEdgeMatrix: unexpected degenerate reference "
        "element.");
  }
  if (det < 0.0) {
    throw std::runtime_error(
        "TElement::referenceEdgeMatrix: reference nodes must stay in "
        "counterclockwise order.");
  }
  return D_R;
}

Matrix2d TElement::referenceRotation() const {
  const Matrix2d D_R = referenceEdgeMatrix();
  return rotationMatrix2D(polarRotationAngle2D(D_R));
}

double TElement::referenceRotationTheta() const {
  return polarRotationAngle2D(referenceEdgeMatrix());
}

Matrix2d TElement::totalBranch(const Matrix2d &history) const {
  return F_P * history;
}

double TElement::totalBranchTheta(const Matrix2d &history) const {
  return polarRotationAngle2D(totalBranch(history));
}

TElement::EdgeFlipRemeshState TElement::evaluateEdgeFlipRemeshState(
    const std::array<GhostNode, 3> &newGhostNodes, const Matrix2d &H_star,
    double theta, double mu) const {
  if (mu < 0.0) {
    throw std::runtime_error(
        "TElement::evaluateEdgeFlipRemeshState: mu must be non-negative.");
  }

  // Prepared from "History Update Across an Edge Flip" (14.05.26) by Sylvain
  // Patinet. This follows the note's notation directly:
  // donor state (P_star, Q_star, F_star, E_star), trial reference rotation
  // Q(theta), then the resulting new state (F_new, C_new, P_new, E_new, H_new).
  const Matrix2d P_star = F_P;
  const Matrix2d Q_star = referenceRotation();
  const Matrix2d F_star = F;
  const Matrix2d E_star = F_star * P_star.inverse();

  const Matrix2d D0 = Matrix2d::Identity();

  EdgeFlipRemeshState state;
  state.theta = theta;
  state.theta_star = polarRotationAngle2D(Q_star);
  state.newGhostNodes = orderNodes(
      newGhostNodes, "TElement::evaluateEdgeFlipRemeshState eIndex=" +
                         std::to_string(eIndex));

  Matrix2d D_C_new;
  D_C_new.col(0) = state.newGhostNodes[1].pos - state.newGhostNodes[0].pos;
  D_C_new.col(1) = state.newGhostNodes[2].pos - state.newGhostNodes[0].pos;

  const Matrix2d A = D_C_new * D0.inverse();
  state.Q_new = rotationMatrix2D(theta);
  state.F_new = A * state.Q_new.transpose();
  state.C_new = state.Q_new * A.transpose() * A * state.Q_new.transpose();

  Matrix2d M_e_new = Matrix2d::Identity();
  elasticReduction(state.C_R_new, state.C_new, M_e_new);
  state.P_new = M_e_new.inverse();
  state.E_new = state.F_new * state.P_new.inverse();
  state.H_new = state.P_new.inverse() * P_star * H_star;

  const double dtheta = wrappedAngleDifference(theta, state.theta_star);
  state.elastic_jump = (state.E_new - E_star).squaredNorm();
  state.rotation_penalty = mu * dtheta * dtheta;
  state.J = state.elastic_jump + state.rotation_penalty;
  state.valid = std::isfinite(state.J) && state.F_new.allFinite() &&
                state.C_new.allFinite() && state.P_new.allFinite() &&
                state.E_new.allFinite() && state.H_new.allFinite();

  // Store the chosen FE reference element directly on the returned ghost
  // nodes so flipEdge() can reuse the evaluated state without re-deriving it.
  setReferenceElementRotation(state.newGhostNodes, theta);

  return state;
}

TElement::EdgeFlipRemeshState TElement::findBestEdgeFlipRemeshStateLinearScan(
    const std::array<GhostNode, 3> &newGhostNodes, const Matrix2d &H_star,
    int nrThetaSamples, double mu) const {
  if (nrThetaSamples < 2) {
    throw std::runtime_error(
        "TElement::findBestEdgeFlipRemeshStateLinearScan: expected at least "
        "2 theta samples.");
  }

  const double dTheta = (2.0 * M_PI) / static_cast<double>(nrThetaSamples - 1);

  EdgeFlipRemeshState bestState;
  bool foundFiniteCandidate = false;
  int rejectedCandidates = 0;
  std::ostringstream firstRejected;
  for (int i = 0; i < nrThetaSamples; ++i) {
    const double theta = -M_PI + dTheta * static_cast<double>(i);
    EdgeFlipRemeshState candidate =
        evaluateEdgeFlipRemeshState(newGhostNodes, H_star, theta, mu);
    if (!candidate.valid) {
      if (rejectedCandidates == 0) {
        firstRejected << "first rejected candidate:\n"
                      << "theta: " << theta << "\n"
                      << "J: " << candidate.J << "\n"
                      << "elastic_jump: " << candidate.elastic_jump << "\n"
                      << "rotation_penalty: " << candidate.rotation_penalty
                      << "\n"
                      << "F_new:\n"
                      << candidate.F_new << "\n"
                      << "C_new:\n"
                      << candidate.C_new << "\n"
                      << "P_new:\n"
                      << candidate.P_new << "\n"
                      << "E_new:\n"
                      << candidate.E_new << "\n"
                      << "H_new:\n"
                      << candidate.H_new << "\n";
      }
      rejectedCandidates++;
      continue;
    }
    if (!foundFiniteCandidate || candidate.J < bestState.J) {
      bestState = candidate;
      foundFiniteCandidate = true;
    }
  }

  if (!foundFiniteCandidate) {
    std::ostringstream oss;
    oss << "TElement::findBestEdgeFlipRemeshStateLinearScan: no finite "
           "edge-flip candidates.\n"
        << "eIndex: " << eIndex << "\n"
        << "nrThetaSamples: " << nrThetaSamples << "\n"
        << "mu: " << mu << "\n"
        << "rejectedCandidates: " << rejectedCandidates << "\n"
        << "F:\n"
        << F << "\n"
        << "F_P:\n"
        << F_P << "\n"
        << "H_star:\n"
        << H_star << "\n";
    for (int i = 0; i < 3; ++i) {
      const GhostNode &gn = newGhostNodes[i];
      oss << "inputGhost[" << i << "]: refId=" << gn.referenceId.i
          << " id=(" << gn.id.x() << "," << gn.id.y() << ")"
          << " pShift=(" << gn.periodicShift.x() << ","
          << gn.periodicShift.y() << ")"
          << " pos=(" << gn.pos.x() << ", " << gn.pos.y() << ")"
          << " ref=(" << gn.ref_pos.x() << ", " << gn.ref_pos.y() << ")"
          << " u=(" << gn.u.x() << ", " << gn.u.y() << ")\n";
    }
    oss << firstRejected.str();
    throw std::runtime_error(oss.str());
  }

  return bestState;
}

double TElement::edgeFlipElasticJumpObjective(
    const std::array<GhostNode, 3> &newGhostNodes, const Matrix2d &H_star,
    double theta, double mu) const {
  return evaluateEdgeFlipRemeshState(newGhostNodes, H_star, theta, mu).J;
}

void TElement::refreshCurrentGhostGeometryForDebug(const Mesh &mesh) {
  m_updatePosition(mesh);
  updateAngleNode();
  m_update_G();
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

Matrix2d rotationMatrix2D(double theta) {
  Matrix2d R;
  R(0, 0) = std::cos(theta);
  R(1, 1) = std::cos(theta);
  R(0, 1) = -std::sin(theta);
  R(1, 0) = std::sin(theta);
  return R;
}

double wrappedAngleDifference(double theta, double thetaReference) {
  constexpr double twoPi = 6.283185307179586476925286766559;
  return std::remainder(theta - thetaReference, twoPi);
}

double polarRotationAngle2DStaticReference(std::array<GhostNode, 3> E) {
  setReferenceElementRotation(E, 0);
  return polarRotationAngle2D(E);
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
  elasticReduction(C_R, C, M_e);
  F_P_out = M_e.inverse();
  Matrix2d F_E = F * M_e;

  return squareTraceStretch(F_E);
}

namespace {
std::string formatDebugVector(const Vector2d &v) {
  std::ostringstream oss;
  oss << std::setprecision(17) << "(" << v.x() << ", " << v.y() << ")";
  return oss.str();
}

std::string formatDebugVector(const Vector2i &v) {
  std::ostringstream oss;
  oss << "(" << v.x() << ", " << v.y() << ")";
  return oss.str();
}

double signedAreaTwice(const Vector2d &a, const Vector2d &b,
                       const Vector2d &c) {
  const Vector2d ab = b - a;
  const Vector2d ac = c - a;
  return ab.x() * ac.y() - ab.y() * ac.x();
}

std::string formatOrderNodesInput(const std::array<GhostNode, 3> &nodes,
                                  int angleNode, double det, double det_eps,
                                  const std::string &context) {
  std::ostringstream oss;
  oss << std::setprecision(17)
      << "orderNodes: unexpected degenerate current triangle.\n"
      << "context: " << context << "\n"
      << "selectedAngleNode: " << angleNode << "\n"
      << "currentSignedArea2(input order): "
      << signedAreaTwice(nodes[0].pos, nodes[1].pos, nodes[2].pos) << "\n"
      << "referenceSignedArea2(input order): "
      << signedAreaTwice(nodes[0].ref_pos, nodes[1].ref_pos, nodes[2].ref_pos)
      << "\n"
      << "selectedDet: " << det << "\n"
      << "detEpsilon: " << det_eps << "\n"
      << "pairwiseCurrentDistances: d01=" << (nodes[1].pos - nodes[0].pos).norm()
      << " d12=" << (nodes[2].pos - nodes[1].pos).norm()
      << " d20=" << (nodes[0].pos - nodes[2].pos).norm() << "\n"
      << "pairwiseReferenceDistances: d01="
      << (nodes[1].ref_pos - nodes[0].ref_pos).norm()
      << " d12=" << (nodes[2].ref_pos - nodes[1].ref_pos).norm()
      << " d20=" << (nodes[0].ref_pos - nodes[2].ref_pos).norm() << "\n";

  bool duplicateReferenceIds = false;
  bool duplicateGhostIds = false;
  for (int i = 0; i < 3; ++i) {
    for (int j = i + 1; j < 3; ++j) {
      duplicateReferenceIds =
          duplicateReferenceIds || nodes[i].referenceId == nodes[j].referenceId;
      duplicateGhostIds = duplicateGhostIds || nodes[i].id == nodes[j].id;
    }
  }
  oss << "duplicateReferenceIds: " << duplicateReferenceIds << "\n"
      << "duplicateGhostIds: " << duplicateGhostIds << "\n";

  for (int i = 0; i < 3; ++i) {
    const GhostNode &gn = nodes[i];
    oss << "node[" << i << "]: refId=" << gn.referenceId.i
        << " id=" << formatDebugVector(gn.id)
        << " periodicShift=" << formatDebugVector(gn.periodicShift)
        << " pos=" << formatDebugVector(gn.pos)
        << " ref=" << formatDebugVector(gn.ref_pos)
        << " u=" << formatDebugVector(gn.u) << "\n";
  }
  return oss.str();
}
} // namespace

std::array<GhostNode, 3> orderNodes(std::array<GhostNode, 3> unorderedNodes) {
  return orderNodes(unorderedNodes, "orderNodes");
}

std::array<GhostNode, 3> orderNodes(std::array<GhostNode, 3> unorderedNodes,
                                    const std::string &context) {
  // We always want a defined local FE ordering:
  // angle node, then the two co-nodes in counterclockwise order around it.
  int angleNode = findAngleNode(unorderedNodes);

  if (angleNode < 0 || angleNode >= 3) {
    throw std::runtime_error("orderNodes: invalid angleNode in " + context +
                             ".");
  }
  int index1 = (angleNode + 1) % 3;
  int index2 = (angleNode + 2) % 3;
  GhostNode a1 = unorderedNodes[angleNode];
  GhostNode g1 = unorderedNodes[index1];
  GhostNode g2 = unorderedNodes[index2];

  const Vector2d edge1 = g1.pos - a1.pos;
  const Vector2d edge2 = g2.pos - a1.pos;
  const double det = edge1.x() * edge2.y() - edge1.y() * edge2.x();
  const double det_eps =
      1e-12 * std::max({1.0, edge1.squaredNorm(), edge2.squaredNorm()});
  if (std::abs(det) < det_eps) {
    throw std::runtime_error(formatOrderNodesInput(unorderedNodes, angleNode,
                                                   det, det_eps, context));
  }

  if (det > 0.0) {
    return {a1, g1, g2};
  }
  return {a1, g2, g1};
}

int findAngleNode(std::array<GhostNode, 3> nodes) {
  // Pick the node opposite the longest edge (largest angle in Euclidean
  // triangle). This is faster than computing angles.
  double largestLength = -1.0;
  int index = -1;

  for (int i = 0; i < 3; ++i) {
    const int next = (i + 1) % 3;
    const int prev = (i + 2) % 3;

    const Vector2d edge = nodes[next].pos - nodes[prev].pos;
    const double len2 = edge.squaredNorm();

    // prefer strictly longer
    if (len2 > largestLength) {
      largestLength = len2;
      index = i;
    }
  }
  return index;
}

void setReferenceElementRotation(std::array<GhostNode, 3> &nodes,
                                 double theta) {
  // Start from the canonical counterclockwise right-isosceles reference
  // triangle and rotate it by theta. We order first so both current and
  // reference node slots use the same local FE orientation.
  nodes = orderNodes(nodes, "setReferenceElementRotation");
  nodes[0].ref_pos = {-0.5, -0.5};
  nodes[1].ref_pos = {0.5, -0.5};
  nodes[2].ref_pos = {-0.5, 0.5};
  const Matrix2d R = rotationMatrix2D(theta);
  for (int i = 0; i < 3; i++) {
    nodes[i].transformReferencePosition(R);
  }
}
