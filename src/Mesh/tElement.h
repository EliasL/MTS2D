#ifndef TELEMENT_H
#define TELEMENT_H
#include "Data/cereal_help.h"
#pragma once
#include "Eigen/Core"
#include "compare_macros.h"
#include "node.h"
#include <array>
#include <cereal/types/array.hpp> // Cereal serialization for std::vector
#include <iostream>

// Profiling helper: prevent inlining of selected hot functions.
#ifndef MTS_NOINLINE
#if defined(__GNUC__) || defined(__clang__)
#define MTS_NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define MTS_NOINLINE __declspec(noinline)
#else
#define MTS_NOINLINE
#endif
#endif

// A triangle can be described by two vectors. In this element, we compute
// the angle between these two vectors. To be clear about which vectors we
// choose, we have one corner_node, and to vector nodes. When we use a vector
// from one node to another in the triangle, it will always be from the
// corner_node to one of the two vector nodes. (Except with the angle)
#define CORNER_NODE 0
#define VECTOR_NODE1 1
#define VECTOR_NODE2 2

using namespace Eigen;

// Declaration
class Mesh;

/**
 * @brief Represents a triangular element in a material surface, characterized
 * by its physical properties.
 *
 * A triangular element is formed by a triangle of nodes and contains
 * information about the deformation gradient F,
 * the metric tensor C, the reduced metric tensor C_, the reduction
 * transformation matrix m, the reduces stress tensor r_s and the
 * Piola stress tensor P.
 *
 *  The shape functions used are:
 *  N1 = 1 - ξ1 - ξ2
 *  N2 = ξ1
 *  N3 = ξ2
 *
 *  u referes to displacement.
 *  X referes to the reference state.
 *  x refers to the current state.
 *
 */
class TElement {
public:
  // Constructor for the triangular element. Initializes the 3 defining
  // GhostNodes and calculates the inverse state A_inv, to later be used in
  // calculating F.
  // Angle node, coAngleNode1, coAngleNode2
  TElement(Mesh &mesh, GhostNode an, GhostNode cn1, GhostNode cn2,
           int elementIndex, double noise = 1,
           std::string energyFunction = "contiSquare", double bulkModulus = 4);
  TElement(GhostNode an, GhostNode cn1, GhostNode cn2,
           std::string energyFunction = "contiSquare", double bulkModulus = 4);
  TElement(std::array<Vector2d, 3> currentNodes,
           std::array<Vector2d, 3> referenceNodes,
           std::string energyFunction = "contiSquare", double bulkModulus = 4);
  TElement() {};
  // Id of nodes associated with elements
  // Don't modify this list, create a new TElement instead. This is so that
  // addElementIndices is run properly (and not forgotten about)
  // Note that the ghost nodes will have a different reference positions than
  // their real node counterparts. This is a complicated and technical detail.
  // It is related to reconnection and preserving plastic deformation history.
  std::array<GhostNode, 3> ghostNodes;

  Matrix2d F;   // Deformation gradient with respect to reference element
  Matrix2d F_P; // Plastic deformation
  Matrix2d F_E; // Elastic deformation

  // Metric tensor (C = F^TF)
  Matrix2d C;

  // Reduced metric tensor
  Matrix2d C_R;

  // Element metric (Not associated with reference state. Just a gram matrix of
  // vecotrs in the element)
  Matrix2d G;

  Matrix2d M_l; // Lagrange reduction transformation matrix (m^TCm = C_)
  Matrix2d M_e; // Elastic reduction transofmration matrix

  // Second Piola-Kirchhoff stress tensor, representing the stress in the
  // reference configuration.
  Matrix2d S;

  // First Piola-Kirchhoff stress tensor, representing the stress in the
  // reference configuration.
  Matrix2d P;

  // Cauchy Stress tensor, representing the stress in the current configuration.
  Matrix2d sigma;

  // Strain energy of the cell, representing the potential energy stored due
  // to deformation.
  double energy = 0;
  // Energy function parameters. Using the contiSquare potential, setting beta
  // to -0.25 gives us a square potential, and setting it to 4 gives us a
  // triangular potential.
  double beta = -0.25;
  // Bulk modulus. Controlls the contribution of the volumetric energy function.
  double K = 4.0;

  // Derivatives
  static Matrix<double, 2, 3> dN_dxi;

  // The jacobian of the shapefunction with respect to the reference state.
  //  We use this to calculate the deformation gradiant F
  // ∂ξ/∂X
  Matrix2d dxi_dX;
  Matrix2d du_dxi;
  Matrix2d dX_dxi;

  // These are adjustment vectors that we multiply together with the piola
  // tensor to correctly extract the force corresponding to each node.
  // Similarly to dxi_dX, these only update once, during initialization.
  Matrix<double, 3, 2> dN_dX;

  // We only save data when plasticity occurs, so we keep a reference of
  // how many times m3 is applied in the lagrange reduction. If this number
  // changes from one reduction to another, we know that a plastic event has
  // occured. (ie. the energy potential suddenly has a gradient in a new
  // direction, ie. the node has fallen into a different energy well.)
  int m3Nr = 0;
  // This keeps track of the number of m3 shears in the previous STABLE STATE
  // for the fixed-reference reduction (used to detect plastic changes without
  // requiring a full update).
  int pastM3Nr = 0;
  // This is the number of m3 shears that have occured in the last minimization
  // step (fixed-reference reduction).
  int pastStepM3Nr = 0;

  // Elastic-domain quadrant labels (1..4), or 0 if outside.
  int red_quadrant = 0;

  // Index of element. Used for debugging.
  int eIndex;

  // A noise value which will slightly distort the volumetric energy term of the
  // element.
  double noise = 1;

  double largestAngle = 0;
  double smallestAngle = 0;
  // This is the node from which the angle is largest
  int angleNode = 0;

private:
  // Initial area (reference triangle area in the undeformed configuration).
  // Used to scale energy density and forces.
  double initArea = 0.0;

  // A variable to store the ground state energy to set our ground state energy
  // to be zero
  double groundStateEnergyDensity = 0.0;

  // Function to compute the ground state energy density
  double computeGroundStateEnergyDensity() const;

public:
  void setInitArea(double area) { initArea = area; }

  void postLoadInit();
  /**
   * @brief Initializes TElement and calculates several values:
   *
   *  the deformation gradient D,
   *  the metric tension C,
   *  the transformation matrix m,
   *  and the reduced metric tension C_.
   *
   */

  // Update forces/energy using fixed reference.
  MTS_NOINLINE void updateForces(const Mesh &mesh);
  // Update angleNode + G using current ghost node positions.
  void updateGeometry();
  // Update remaining derived fields (F/C/m/sigma/angles) assuming forces and
  // geometry are already up to date.
  void updateFull();

  // Usefull if you only care about the energy given the C matrix.
  double calculateEnergyDensity(double c11, double c22, double c12) const;

  // Used for testing the lagrange reuction functions
  static TElement reduce_element(double c11, double c22, double c12);

  // Find the longest edge, and set the angleNode to the node opposite that edge
  void updateAngleNode();
  // Compute all angles in mesh, and store the largest one
  void updateAngles();

  // Two elements can be seen as forming a rombus together. This function
  // returns the index of the element that is accross from the node forming the
  // largest angle in the current element, but only if that element is similar
  // in shape. If the other element is not similar, it returns -1
  int getElementTwin(const Mesh &mesh) const;

  std::array<const GhostNode, 3> getAngleCo1Co2Nodes() const;
  std::array<const GhostNode, 3> getGhostNodes() const {
    return getAngleCo1Co2Nodes();
  }
  void setReferenceElement();
  void setReferenceElement(const std::array<Vector2d, 3> &refNodes);
  void setReferenceElement(const std::array<const GhostNode, 3> &refNodes);
  void setReferenceElementFromCurrentState(const Mesh &mesh);
  void deformReferenceElement(Matrix2d F, Vector2d oldAnchor = Vector2d::Zero(),
                              Vector2d newAnchor = Vector2d::Zero());
  Vector2d referenceCentroidShiftToCurrent() const;
  void refreshCurrentGhostGeometryForDebug(const Mesh &mesh);
  Matrix2d referenceEdgeMatrix() const;
  Matrix2d referenceRotation() const;
  double referenceRotationTheta() const;
  Matrix2d totalBranch(const Matrix2d &history) const;
  double totalBranchTheta(const Matrix2d &history) const;
  struct EdgeFlipRemeshState {
    std::array<GhostNode, 3> newGhostNodes;
    double theta = 0.0;
    double theta_star = 0.0;
    double elastic_jump = std::numeric_limits<double>::infinity();
    double rotation_penalty = 0.0;
    double J = std::numeric_limits<double>::infinity();
    Matrix2d Q_new = Matrix2d::Identity();
    Matrix2d F_new = Matrix2d::Identity();
    Matrix2d C_new = Matrix2d::Identity();
    Matrix2d C_R_new = Matrix2d::Identity();
    Matrix2d P_new = Matrix2d::Identity();
    Matrix2d E_new = Matrix2d::Identity();
    Matrix2d H_new = Matrix2d::Identity();
  };
  EdgeFlipRemeshState
  evaluateEdgeFlipRemeshState(const std::array<GhostNode, 3> &newGhostNodes,
                              const Matrix2d &H_star, double theta,
                              double mu = 0.0) const;
  EdgeFlipRemeshState findBestEdgeFlipRemeshStateLinearScan(
      const std::array<GhostNode, 3> &newGhostNodes, const Matrix2d &H_star,
      int nrThetaSamples = 1000, double mu = 0.0) const;
  double
  edgeFlipElasticJumpObjective(const std::array<GhostNode, 3> &newGhostNodes,
                               const Matrix2d &H_star, double theta,
                               double mu = 0.0) const;

  std::array<const GhostNode *, 2> getCoAngleNodes() const;

  const GhostNode *getAngleNode() const;

  // Get center of mass of the element
  Vector2d getCom();

  double area() const;

private:
  void updateReferenceGeometry();

  // Copy the displacement from the real nodes to the nodes in the element
  void m_updatePosition(const Mesh &mesh);

  void m_update_du_dxi();
  void m_update_dX_dxi();
  void m_update_dN_dX();

  // Computes the deformation gradient for the cell based on the triangle's
  // vertices.
  void m_updateDeformationGradient();

  // Computes the metric tensor for the triangle (real C only).
  void m_updateMetricTensor();

  // Calculates the element metric G
  void m_update_G();

  // Performs a Lagrange reduction on C to calculate C_.
  void m_lagrangeReduction();
  // Normal lagrange reduction only (assumes C is up to date).
  void m_elasticReduction();

  void m_update_plastic_elastic_F();

  // Calculates energy Phi
  MTS_NOINLINE void m_updateEnergy();

  // Calculate reduced stress
  // Gradient of energy function Phi with respect to reduced metric tensor C_
  MTS_NOINLINE void m_updateSecondPiolaStress();

  // Calculate Piola stress P
  MTS_NOINLINE void m_updateFirstPiolaStress();

  // Calculate Cauchy stress sigma
  void m_updateCauchyStress();

  // Calculate the resolved-shear stress
  void m_updateShearStress();

  // Calculate force on each node
  MTS_NOINLINE void m_updateForceOnEachNode();

  friend class cereal::access;
  template <class Archive> void save(Archive &ar) const {
    ar(MAKE_NVP(ghostNodes), MAKE_NVP(m3Nr), MAKE_NVP(pastM3Nr),
       MAKE_NVP(pastStepM3Nr), MAKE_NVP(eIndex), MAKE_NVP(noise),
       MAKE_NVP(dX_dxi), MAKE_NVP(beta), MAKE_NVP(K));
  }
  template <class Archive> void load(Archive &ar) {
    ar(MAKE_NVP(ghostNodes), MAKE_NVP(m3Nr), MAKE_NVP(eIndex), MAKE_NVP(noise),
       MAKE_NVP(dX_dxi));
    LOAD_WITH_DEFAULT(ar, pastM3Nr, 0);
    LOAD_WITH_DEFAULT(ar, pastStepM3Nr, 0);
    // Backward compatibility: older dumps used pastM3Nr/pastStepM3Nr.
    int oldPastM3 = 0;
    int oldPastStepM3 = 0;
    loadWithDefault(ar, "pastM3Nr", oldPastM3, 0);
    loadWithDefault(ar, "pastStepM3Nr", oldPastStepM3, 0);
    if (pastM3Nr == 0 && oldPastM3 != 0) {
      pastM3Nr = oldPastM3;
    }
    if (pastStepM3Nr == 0 && oldPastStepM3 != 0) {
      pastStepM3Nr = oldPastStepM3;
    }
    // Backward compatibility: ignore old m1Nr/m2Nr fields if present.
    int dummy_m1 = 0;
    int dummy_m2 = 0;
    loadWithDefault(ar, "m1Nr", dummy_m1, 0);
    loadWithDefault(ar, "m2Nr", dummy_m2, 0);
    LOAD_WITH_DEFAULT(ar, beta, beta);
    LOAD_WITH_DEFAULT(ar, K, K);
    // Backward compatibility: ignore old cached shorthand fields.
    double ignoredReferenceTheta = 0.0;
    double ignoredThetaOffset = 0.0;
    loadWithDefault(ar, "referenceTheta", ignoredReferenceTheta, 0.0);
    loadWithDefault(ar, "thetaOffset", ignoredThetaOffset, 0.0);
    // Backward compatibility: older dumps may store the two co-nodes in the
    // opposite order. Swap the full ghost-node entries so the new CCW
    // reference-orientation invariant holds without changing F.
    Matrix2d loadedReferenceEdges;
    loadedReferenceEdges.col(0) = ghostNodes[1].ref_pos - ghostNodes[0].ref_pos;
    loadedReferenceEdges.col(1) = ghostNodes[2].ref_pos - ghostNodes[0].ref_pos;
    if (loadedReferenceEdges.determinant() < 0.0) {
      std::swap(ghostNodes[1], ghostNodes[2]);
    }
    postLoadInit();
  }

  // Giving access to private variables
  friend bool compareTElementsInternal(const TElement &lhs, const TElement &rhs,
                                       std::string *debugMsg, int tabNumber);
};
// Non-member functions

// This function updates the list of connected elements in the real nodes
void addElementIndices(Mesh &mesh, const std::array<GhostNode, 3> &nodeList,
                       int elementIndex);

double triangleArea(Vector2d posA, Vector2d posB, Vector2d posC);
double tElementArea(const std::array<GhostNode, 3> &E);
double tElementArea(const GhostNode &A, const GhostNode &B, const GhostNode &C);
double tElementInitialArea(const std::array<GhostNode, 3> &gn);

Matrix2d tElementF(const std::array<GhostNode, 3> &E);

double polarRotationAngle2D(const Matrix2d &F);
double polarRotationAngle2D(const std::array<GhostNode, 3> &E);
double polarRotationAngle2DStaticReference(std::array<GhostNode, 3> E);

double squareTraceStretch(const Eigen::Matrix2d &F);

double distanceFromIntegerShear(const Matrix2d &F);
double distanceFromIntegerShear(const Matrix2d &F, Matrix2d &F_P_out);

void setReferenceElementRotation(std::array<GhostNode, 3> &nodes, double theta);

// Management functions

std::ostream &operator<<(std::ostream &os, const TElement &element);

inline bool compareTElementsInternal(const TElement &lhs, const TElement &rhs,
                                     std::string *debugMsg = nullptr,
                                     int tabNumber = 0) {
  bool equal = true;

  // Compare public members.
  COMPARE_FIELD(ghostNodes);
  COMPARE_FIELD(C);
  COMPARE_FIELD(C_R);
  COMPARE_FIELD(M_l);
  COMPARE_FIELD(S);
  COMPARE_FIELD(P);
  COMPARE_FIELD(energy);
  COMPARE_FIELD(dN_dX);
  COMPARE_FIELD(m3Nr);
  COMPARE_FIELD(pastM3Nr);
  COMPARE_FIELD(pastStepM3Nr);
  COMPARE_FIELD(red_quadrant);
  COMPARE_FIELD(eIndex);
  COMPARE_FIELD(noise);
  COMPARE_FIELD(largestAngle);
  COMPARE_FIELD(angleNode);

  // Compare private members.
  COMPARE_FIELD(initArea);
  COMPARE_FIELD(beta);
  COMPARE_FIELD(K);

  return equal;
}

/*
   Standard equality operator for TElement.
   Declared as a friend in TElement, it calls compareTElementsInternal without
   generating debug messages.
*/
inline bool operator==(const TElement &lhs, const TElement &rhs) {
  return compareTElementsInternal(lhs, rhs, nullptr);
}
inline bool operator!=(const TElement &lhs, const TElement &rhs) {
  return !(lhs == rhs);
}

/*
   Debug function for TElement that uses the same internal comparison logic.
   Returns a string describing which fields differ between the two objects.
*/
inline std::string debugCompare(const TElement &lhs, const TElement &rhs,
                                int tabNumber = 0) {
  std::string diff;

  // If sizes match, compare each element
  for (size_t i = 0; i < lhs.ghostNodes.size(); i++) {
    if (!(lhs.ghostNodes[i] == rhs.ghostNodes[i])) {
      diff += std::string(tabNumber, '\t') + "tElementNodes[" +
              std::to_string(i) + "] differs -> \n";
      // Recursively call debugCompare for TElement
      diff += debugCompare(lhs.ghostNodes[i], rhs.ghostNodes[i], tabNumber + 1);
    }
  }
  compareTElementsInternal(lhs, rhs, &diff, tabNumber);
  return diff;
}

#endif
