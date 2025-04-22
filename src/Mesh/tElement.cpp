#include "tElement.h"
#include "Eigen/Core"
#include "Eigen/src/Core/Matrix.h"
#include "Mesh/node.h"
#include "Simulation/energyFunctions.h"
#include "mesh.h"
#include <Eigen/LU>
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <stdexcept>

// Define and initialize static members
double TElement::groundStateEnergyDensity =
    TElement::computeGroundStateEnergyDensity();

TElement::TElement(Mesh &mesh, GhostNode an, GhostNode cn1, GhostNode cn2,
                   int elementIndex, double noise)
    : ghostNodes{an, cn1, cn2}, F(Matrix2d::Identity()),
      C(Matrix2d::Identity()), C_R(Matrix2d::Identity()),
      m(Matrix2d::Identity()), sigma(Matrix2d::Identity()),
      P(Matrix2d::Identity()), eIndex(elementIndex), noise(noise) {

  // Add this element to the nodes it is created by
  addElementIndices(mesh, {an, cn1, cn2}, elementIndex);

  // Shape functions
  // Matrix<double, 2, 3> dN_dxi;
  dN_dX << -1.0, -1.0, //
      /* */ 1.0, 0.0,  //
      /* */ 0.0, 1.0;
  // Matrix<double, 2, 3> dX_dxi;
}

void TElement::update(const Mesh &mesh) {
  // The order here is very important.

  // Update nodes in element
  m_updatePosition(mesh);

  // Find and update largest angle
  updateLargestAngle();

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

void TElement::m_updateDeformationGradiant() {

  // Deformed coordinates x
  Matrix<double, 2, 3> x;
  x.col(0) = ghostNodes[0].pos;
  x.col(1) = ghostNodes[1].pos;
  x.col(2) = ghostNodes[2].pos;

  F = x * dN_dX;

  assert(F.determinant() != 0);
}

// Provices a metric tensor for the triangle
void TElement::m_updateMetricTensor() {
  // Discontinuous yielding of pristine micro-crystals - page 8/207
  C = F.transpose() * F;
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

void TElement::m_updateReducedStress() {
  // sigma = 1/2 (∂Φ/∂C_R + (∂Φ/∂C_R)^T)
  sigma =
      ContiPotential::stress(C_R(0, 0), C_R(1, 1), C_R(0, 1), beta, K, noise);
}

// Calculate Piola stress tensor and force on each node from current cell
void TElement::m_updatePiolaStress() {
  // Transform back from lagrange-reudced to un-reduced
  // sigma = 1/2 (∂Φ/∂C_R + (∂Φ/∂C_R)^T)
  // so it's not actually quite dPhi_dC
  Matrix2d dPhi_dC = m * sigma * m.transpose();
  //  Discontinuous yielding of pristine micro-crystals, page 16/215
  // Calculate piola tensor
  P = 2.0 * F * dPhi_dC;
}

void TElement::m_updateForceOnEachNode() {
  // dPhi_du = P*dN_dX is the energy density gradient
  Matrix<double, 2, 3> dPhi_du = P * dN_dX.transpose();
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

  int maxLoops = 1e6;

  // We start by copying the values from C to the reduced matrix
  // Note that we only modify C_[0][1]. At the end of the algorithm,
  // we copy C_[0][1] to C_[1][0]
  C_R = C;

  // Then reset some values
  m = m.Identity();
  simple_m = simple_m.Identity();
  // We should also reset m and m3Nr
  m1Nr = 0;
  m2Nr = 0;
  m3Nr = 0;

  // Now we keep repeating this loop until C_ does not change
  bool changed = true;
  // To prevent an infinite loop, we use a for loop just in case
  for (int i = 0; i < maxLoops; i++) {
    changed = false;

    if (C_R(0, 1) < 0) {
      C_R(0, 1) = -C_R(0, 1);
      lag_m1(m);
      changed = true;
      m1Nr = +1;
    }

    if (C_R(1, 1) < C_R(0, 0)) {
      std::swap(C_R(0, 0), C_R(1, 1));
      lag_m2(m);
      changed = true;
      m2Nr += 1;
    }

    if (2 * C_R(0, 1) > C_R(0, 0)) {
      // The order here matters, don't modify C_(0,1) before using it
      // to calculate C_(1,1).
      C_R(1, 1) += C_R(0, 0) - 2 * C_R(0, 1);
      C_R(0, 1) -= C_R(0, 0);
      lag_m3(m);
      m3Nr += 1;
      changed = true;
    }
    // If we have not changed, we break out of the loop
    if (changed == false) {
      break;
    }
  }

  // If changed is true, it means we hit the max loop iterations
  if (changed) {
    std::cout << "Lagrange Reduction Iteration Counts:\n"
              << "  m1Nr: " << m1Nr << "\n"
              << "  m2Nr: " << m2Nr << "\n"
              << "  m3Nr: " << m3Nr << "\n\n";
    std::cout << "Deformation Gradient F:\n"
              << F << "\n\n"
              << "Metric Tensor C:\n"
              << C << "\n\n"
              << "Reduced Metric Tensor C_R:\n"
              << C_R << "\n\n";
    std::cout << "Ghost Nodes:\n"
              << "  Node 0: " << ghostNodes[0] << "\n"
              << "  Node 1: " << ghostNodes[1] << "\n"
              << "  Node 2: " << ghostNodes[2] << "\n";
    throw std::runtime_error("Stuck in lagrange reduction.\n");
  }

  C_R(1, 0) = C_R(0, 1);
}

void TElement::updateLargestAngle() {
  double maxAngle = 0.0;
  int largestAngleIndex = 0;

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
      // TODO acos is slow. Use C12
      double angle = std::acos(cosAngle);

      // Convert to degrees
      angle *= 180.0 / M_PI;

      // Track the largest angle
      if (angle > maxAngle) {
        maxAngle = angle;
        largestAngleIndex = i;
      }
    }
  }

  // Store the results
  angleNode = largestAngleIndex;
  largestAngle = maxAngle;
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

  // Identify the two nodes to the side of the angle node
  auto coAngleNodes = getCoAngleNodes();
  const Node *n1 = mesh[coAngleNodes[0]->referenceId];
  const Node *n2 = mesh[coAngleNodes[1]->referenceId];

  // Find all elements that are common for both nodes and not this element
  for (int elementFromNode1 : n1->elementIndices) {
    // Skip the current element or end of valid elements
    if (elementFromNode1 == eIndex || elementFromNode1 == -1) {
      continue;
    }

    for (int elementFromNode2 : n2->elementIndices) {
      // Skip the current element or end of valid elements
      if (elementFromNode2 == eIndex || elementFromNode2 == -1) {
        continue;
      }

      // If we find an element that contains both nodes (and that is not this
      // element)
      if (elementFromNode1 == elementFromNode2) {
        // We now check that the two nodes they share are coAngleNodes
        TElement twin = mesh.elements[elementFromNode1];
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
double TElement::calculateEnergyDensity(double c11, double c22,
                                        double c12) const {
  TElement e = TElement();
  e.C = Matrix2d{{c11, c12}, {c12, c22}};
  e.m_lagrangeReduction();
  return ContiPotential::energyDensity(e.C_R(0, 0), e.C_R(1, 1), e.C_R(0, 1),
                                       beta, K);
}

TElement TElement::lagrangeReduction(double c11, double c22, double c12) {
  TElement element = TElement();
  element.C = Matrix2d{{c11, c12}, {c12, c22}};
  element.m_lagrangeReduction();
  return element;
}

Vector2d TElement::getCom() {
  return (ghostNodes[0].pos + ghostNodes[1].pos + ghostNodes[2].pos) / 3;
}

//------- Non TElement functions

void addElementIndices(Mesh &mesh, const std::array<GhostNode, 3> nodeList,
                       int elementIndex) {
  int i = 0;
  for (GhostNode gn : nodeList) {

    // Reference to the current count
    int &count = mesh[gn]->elementCount;
    // Ensure we don't exceed the array size
    if (count < MAX_ELEMENTS_PER_NODE) {
      mesh[gn]->elementIndices[count] = elementIndex;
      mesh[gn]->nodeIndexInElement[count] = i;
      ++count; // Increment the count for the node
    } else {
      // Handle overflow (e.g., log an error or take other measures)
      throw std::overflow_error("Element index overflow for node " +
                                std::to_string(gn.referenceId.i));

      std::cerr << "Error: Too many elements for node " << gn.referenceId
                << std::endl;
    }
    i++;
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
double tElementInitialArea(const GhostNode &A, const GhostNode &B,
                           const GhostNode &C) {
  return triangleArea(A.init_pos, B.init_pos, C.init_pos);
}

double tElementArea(const GhostNode &A, const GhostNode &B,
                    const GhostNode &C) {
  return triangleArea(A.pos, B.pos, C.pos);
}
