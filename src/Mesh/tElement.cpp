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

/*
Shape functions:
N1 = 1 - ξ1 - ξ2
N2 = ξ1
N3 = ξ2

Derivative of shape functions:
-1.0, 1.0, 0.0,
-1.0, 0.0, 1.0
*/
const Matrix<double, 2, 3> TElement::dN_dxi =
    (Matrix<double, 2, 3>() << -1.0, 1.0, 0.0, //
     /*                     */ -1.0, 0.0, 1.0)
        .finished();

TElement::TElement(Mesh &mesh, GhostNode an, GhostNode cn1, GhostNode cn2,
                   int elementIndex, double noise)
    : ghostNodes{an, cn1, cn2}, F(Matrix2d::Identity()),
      C(Matrix2d::Identity()), C_(Matrix2d::Identity()),
      m(Matrix2d::Identity()), dPhi_dC_(Matrix2d::Identity()),
      P(Matrix2d::Identity()), eIndex(elementIndex), noise(noise) {

  // Add this element to the nodes it is created by
  addElementIndices(mesh, {an, cn1, cn2}, elementIndex);

  dN_dX[0] = {-1, -1};
  dN_dX[1] = {1, 0};
  dN_dX[2] = {0, 1};
  dxi_dX = Matrix2d::Identity();
  dX_dxi = Matrix2d::Identity();

  // Calculate initial area
  initArea = tElementInitialArea(an, cn1, cn2);

  // Calculate ground state energy density
  groundStateEnergyDensity = calculateEnergyDensity(1, 1, 0);

  // Calculate angle to check that it is 90 degrees. (For debugging, we don't)
  // want to accidentally create elements that are not 90 degrees
  m_updateLargestAngle();
  // assert(largestAngle == 90);
  //  Another debugging assetion test. To make sure we know exactly what we are
  //  doing, we have the convention of always putting the node with the lower
  //  index as the first vector node. This is not strictly neccecary, but can
  //  help with debugging and understanding if we can always make this
  //  assumption
  // assert(vn1.ghostId.i < vn2.ghostId.i);
}

void TElement::update(const Mesh &mesh) {
  // The order here is very important.

  // Update nodes in element
  m_updatePosition(mesh);

  // Find and update largest angle
  m_updateLargestAngle();

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
  du_dxi.col(0) = du(ghostNodes[CORNER_NODE], ghostNodes[VECTOR_NODE1]);
  du_dxi.col(1) = du(ghostNodes[CORNER_NODE], ghostNodes[VECTOR_NODE2]);
  return du_dxi;
}

// Jacobian with respect to the initial position of the nodes ∂X/∂ξ
// See du_dxi for a similar working out.
Matrix2d TElement::m_update_dX_dxi() {
  // ∂X/∂ξ
  dX_dxi.col(0) = dX(ghostNodes[CORNER_NODE], ghostNodes[VECTOR_NODE1]);
  dX_dxi.col(1) = dX(ghostNodes[CORNER_NODE], ghostNodes[VECTOR_NODE2]);
  // Actually, we will make every element always use the same reference state.
  // This works because of the lagrange reduction reducing all states to the
  // same anyway
  // dX_dxi << 1, 0, 0, 1;
  return dX_dxi;
}

void TElement::m_updateDeformationGradiant() {
  // See FEMNotes pdf from Umut

  //                ∂u/∂X = ∂u/∂ξ * ∂ξ/∂X

  // dxi_dX is already computed and is constant
  m_update_du_dxi();
  // Matrix2d du_dX = du_dxi * dxi_dX;
  // F = Matrix2d::Identity() + du_dxi * dxi_dX;

  // Deformed coordinates x
  Matrix<double, 2, 3> x;
  x.col(0) = ghostNodes[0].pos;
  x.col(1) = ghostNodes[1].pos;
  x.col(2) = ghostNodes[2].pos;

  Matrix<double, 2, 3> dN_dX;
  dN_dX << -1.0, 1.0, 0.0, //
      /**/ -1.0, 0.0, 1.0;

  F = x * dN_dX.transpose();
  assert(F.determinant() != 0);
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

  int maxLoops = 1e6;

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
  // To prevent an infinite loop, we use a for loop just in case
  for (int i = 0; i < maxLoops; i++) {
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
              << "Displacement Jacobian du_dxi:\n"
              << du_dxi << "\n\n"
              << "Inverse Jacobian dxi_dX:\n"
              << dxi_dX << "\n\n"
              << "Metric Tensor C:\n"
              << C << "\n\n"
              << "Reduced Metric Tensor C_:\n"
              << C_ << "\n\n";
    std::cout << "Ghost Nodes:\n"
              << "  Node 0: " << ghostNodes[0] << "\n"
              << "  Node 1: " << ghostNodes[1] << "\n"
              << "  Node 2: " << ghostNodes[2] << "\n";
    throw std::runtime_error("Stuck in lagrange reduction.\n");
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

    ghostNodes[i].f = P * dN_dX[i];
  }
}

void TElement::m_updatePosition(const Mesh &mesh) {
  // loop through the three nodes in the elements
  for (size_t i = 0; i < 3; i++) {
    // Get the node from the mesh (seperate from the node inside this element)
    const Node *n = mesh[ghostNodes[i].referenceId];
    ghostNodes[i].updatePosition(n, mesh.currentDeformation, mesh.a);
  }
}

void TElement::m_updateLargestAngle() {
  double maxAngle = 0.0;
  int largestAngleIndex = 0;

  // Compute and compare all three angles
  for (int i = 0; i < 3; i++) {
    int next = (i + 1) % 3;
    int prev = (i + 2) % 3;

    // Compute vectors from the current vertex to adjacent vertices
    Vector2d v1 = dx(ghostNodes[i], ghostNodes[next]);
    Vector2d v2 = dx(ghostNodes[i], ghostNodes[prev]);

    // Compute angle using dot product and vector magnitudes
    double magnitudeProduct = v1.norm() * v2.norm();

    // Avoid division by zero
    if (magnitudeProduct > 1e-10) {
      double cosAngle = std::clamp(v1.dot(v2) / magnitudeProduct, -1.0, 1.0);
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
  return {g1, g2};
}

GhostNode *TElement::getAngleNode() {
  GhostNode *agn = &ghostNodes[angleNode];
  return agn;
}

int TElement::getElementTwin(const Mesh &mesh) const {
  // We start by identifying the two nodes to the side of the angle node
  auto coAngleNodes = getCoAngleNodes();
  const Node *n1 = mesh[coAngleNodes[0]->referenceId];
  const Node *n2 = mesh[coAngleNodes[1]->referenceId];

  // We now find the element that is common for both nodes, and is not this
  // element

  for (int elementFromNode1 : n1->elementIndices) {
    // Skip the current element
    if (elementFromNode1 == eIndex) {
      continue;
    }
    // Check if we've reached the end of valid elements for node1
    else if (elementFromNode1 == -1) {
      break;
    }

    for (int elementFromNode2 : n2->elementIndices) {

      // Skip the current element
      if (elementFromNode2 == eIndex) {
        continue;
      }

      // Check if we've reached the end of valid elements for node2
      else if (elementFromNode2 == -1) {
        break;
      }
      // If we find the matching element, finally check if the element does in
      // fact share the same ghost nodes
      else if (elementFromNode1 == elementFromNode2) {
        // Due to periodic boundary conditions, it is possible to have multiple
        // paris of elements that share two coNodes
        // TODO
        return elementFromNode1;
      }
    }
  }
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
  return ContiPotential::energyDensity(e.C_(0, 0), e.C_(1, 1), e.C_(0, 1), beta,
                                       K);
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
    int &count = mesh.nodes(gn.referenceId.i).elementCount;
    // Ensure we don't exceed the array size
    if (count < MAX_ELEMENTS_PER_NODE) {
      mesh.nodes(gn.referenceId.i).elementIndices[count] = elementIndex;
      mesh.nodes(gn.referenceId.i).nodeIndexInElement[count] = i;
      ++count; // Increment the count for the node
    } else {
      // Handle overflow (e.g., log an error or take other measures)
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

double tElementInitialArea(const GhostNode &A, const GhostNode &B,
                           const GhostNode &C) {
  Vector2d posA = A.init_pos;
  Vector2d posB = B.init_pos;
  Vector2d posC = C.init_pos;

  double area = 0.5 * std::abs(posA[0] * (posB[1] - posC[1]) +
                               posB[0] * (posC[1] - posA[1]) +
                               posC[0] * (posA[1] - posB[1]));
  return area;
}