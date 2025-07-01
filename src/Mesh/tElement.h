#ifndef TELEMENT_H
#define TELEMENT_H
#include "Data/cereal_help.h"
#pragma once
#include "Eigen/Core"
#include "Simulation/energyFunctions.h"
#include "compare_macros.h"
#include "node.h"
#include <array>
#include <cereal/types/array.hpp> // Cereal serialization for std::vector
#include <iostream>

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
  // Id of nodes associated with elements
  // Don't modify this list, create a new TElement instead. This is so that
  // addElementIndices is run properly (and not forgotten about)
  std::array<GhostNode, 3> ghostNodes;

  // Deformation gradient
  Matrix2d F;
  // F_fixed_ref is the deformation gradient using a fixed reference
  // configuration. The true F can sometimes be non-invertable, due to
  // reconnecting creating a 1D reference configuration. F_fixed_ref is used to
  // calculate the force.
  Matrix2d F_fixed_ref;

  // Metric tensor (C = F^TF)
  Matrix2d C;
  Matrix2d C_fixed_ref;

  // Reduced metric tensor
  Matrix2d C_R_fixed_ref;

  // Reduction transformation matrix (m^TCm = C_)
  Matrix2d m;
  Matrix2d simple_m;

  // Reduced stress
  Matrix2d sigma;

  // First Piola-Kirchhoff stress tensor, representing the stress in the
  // reference configuration.
  Matrix2d P;

  // Strain energy of the cell, representing the potential energy stored due
  // to deformation.
  double energy = 0;

  // A representation of stress that is unaffected by the directionality of
  // loading. Discontinuous yielding of pristine micro-crystals (page 216/17)
  double resolvedShearStress = 0;

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
  static Matrix<double, 3, 2> dN_dX_fixed_ref;

  // We only save data when plasticity occurs, so we keep a reference of
  // how many times m3 is applied in the lagrange reduction. If this number
  // changes from one reduction to another, we know that a plastic event has
  // occured. (ie. the energy potential suddenly has a gradient in a new
  // direction, ie. the node has fallen into a different energy well.)
  int m3Nr = 0;
  // This keeps track of the number of m3 shears in the previous STABLE STATE
  // Not previous lagrange reduction, and not previous minimization step.
  // Stable state! Since the load changed.
  int pastM3Nr = 0;
  // This is the number of m3 shears that have occured in the last minimization
  // step.
  int pastStepM3Nr = 0;

  // For completeness, we keep track of m1 and m2 as well
  int m1Nr = 0;
  int m2Nr = 0;

  // Index of element. Used for debugging.
  int eIndex;

  // A noise value which will slightly distort the volumetric energy term of the
  // element.
  double noise = 1;

  // The dotproduct of the element vectors in the current reference. Can be used
  // as the angle of the element.
  // Note that C_12 is not always
  // representative of the angle of the element as it uses the deformation
  // gradient from the reference configuration
  double largestAngle = 0;
  // This is the node from which the angle is largest
  int angleNode = 0;

private:
  // Various numbers used in energy and reduced stress calculation. TODO
  // understand and comment Coresponds (somehow) to square lattice. beta=4 gives
  // triangular lattice.
  static constexpr double beta = -0.25;
  // Bulk modulus. Controlls the contribution of the volumetric energy function.
  static constexpr double K = 4.;

  // Initial area
  // This is used together with the determinant of the deformation gradient
  // to get the current area, and the energy density function to get the
  // energy
  // Now with a constant reference state, we fix it to 0.5
  double initArea = 0.5;

  // A variable to store the ground state energy to set our ground state energy
  // to be zero
  static double groundStateEnergyDensity;

  // Function to compute the ground state energy density
  static double computeGroundStateEnergyDensity() {
    // Assuming calculateEnergyDensity is accessible and noise doesn't matter
    // (or set to zero)
    return ContiPotential::energyDensity(1, 1, 0, -0.25, 4.0);
  }

public:
  // Constructor for the triangular element. Initializes the 3 defining
  // GhostNodes and calculates the inverse state A_inv, to later be used in
  // calculating F.
  // Angle node, coAngleNode1, coAngleNode2
  TElement(Mesh &mesh, GhostNode an, GhostNode cn1, GhostNode cn2,
           int elementIndex, double noise = 1);
  TElement() {};
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

  void update(const Mesh &mesh);

  // Usefull if you only care about the energy given the C matrix.
  double calculateEnergyDensity(double c11, double c22, double c12) const;

  // Used for testing the lagrange reuction functions
  static TElement lagrangeReduction(double c11, double c22, double c12);

  // Compute all angles in mesh, and store the largest one
  void updateLargestAngle();

  // Two elements can be seen as forming a rombus together. This function
  // returns the index of the element that is accross from the node forming the
  // largest angle in the current element, but only if that element is similar
  // in shape. If the other element is not similar, it returns -1
  int getElementTwin(const Mesh &mesh) const;

  std::array<const GhostNode *, 2> getCoAngleNodes() const;

  GhostNode *getAngleNode();

  // Get center of mass of the element
  Vector2d getCom();

private:
  // Copy the displacement from the real nodes to the nodes in the element
  void m_updatePosition(const Mesh &mesh);

  Matrix2d m_update_du_dxi();
  Matrix2d m_update_dX_dxi();

  // Computes the deformation gradient for the cell based on the triangle's
  // vertices.
  void m_updateDeformationGradiant();

  // Computes the metric tensor for the triangle.
  void m_updateMetricTensor();

  // Performs a Lagrange reduction on C to calculate C_.
  void m_lagrangeReduction();

  // Calculates energy Phi
  void m_updateEnergy();

  // Calculate reduced stress
  // Gradient of energy function Phi with respect to reduced metric tensor C_
  void m_updateReducedStress();

  // Calculate Piola stress P
  void m_updatePiolaStress();

  // Calculate the resolved-shear stress
  void m_updateResolvedShearStress();

  // Calculate force on each node
  void m_updateForceOnEachNode();

  // Position subtraction (The vector from node 1 to node 2)
  Vector2d const dx(const GhostNode &n1, const GhostNode &n2) const {
    return n2.pos - n1.pos;
  }

  // Initial-position subtraction
  Vector2d const dX(const GhostNode &n1, const GhostNode &n2) const {
    return n2.init_pos - n1.init_pos;
  }

  // Displacement subtraction
  Vector2d const du(const GhostNode &n1, const GhostNode &n2) const {
    return n2.u - n1.u;
  }

  friend class cereal::access;
  template <class Archive> void save(Archive &ar) const {
    ar(MAKE_NVP(ghostNodes), MAKE_NVP(m3Nr), MAKE_NVP(pastM3Nr),
       MAKE_NVP(pastStepM3Nr), MAKE_NVP(m1Nr), MAKE_NVP(m2Nr), MAKE_NVP(eIndex),
       MAKE_NVP(noise), MAKE_NVP(dX_dxi));
  }
  template <class Archive> void load(Archive &ar) {
    ar(MAKE_NVP(ghostNodes), MAKE_NVP(m3Nr), MAKE_NVP(pastM3Nr),
       MAKE_NVP(pastStepM3Nr), MAKE_NVP(m1Nr), MAKE_NVP(m2Nr), MAKE_NVP(eIndex),
       MAKE_NVP(noise), MAKE_NVP(dX_dxi));
    postLoadInit();
  }

  // Giving access to private variables
  friend bool compareTElementsInternal(const TElement &lhs, const TElement &rhs,
                                       std::string *debugMsg, int tabNumber);
};

// This function updates the list of connected elements in the real nodes
void addElementIndices(Mesh &mesh, const std::array<GhostNode, 3> nodeList,
                       int elementIndex);

std::ostream &operator<<(std::ostream &os, const TElement &element);

inline bool compareTElementsInternal(const TElement &lhs, const TElement &rhs,
                                     std::string *debugMsg = nullptr,
                                     int tabNumber = 0) {
  bool equal = true;

  // Compare public members.
  COMPARE_FIELD(ghostNodes);
  COMPARE_FIELD(F);
  COMPARE_FIELD(F_fixed_ref);
  COMPARE_FIELD(C);
  COMPARE_FIELD(C_R_fixed_ref);
  COMPARE_FIELD(m);
  COMPARE_FIELD(simple_m);
  COMPARE_FIELD(sigma);
  COMPARE_FIELD(P);
  COMPARE_FIELD(energy);
  COMPARE_FIELD(resolvedShearStress);
  COMPARE_FIELD(dN_dX);
  COMPARE_FIELD(m3Nr);
  COMPARE_FIELD(pastM3Nr);
  COMPARE_FIELD(pastStepM3Nr);
  COMPARE_FIELD(m1Nr);
  COMPARE_FIELD(m2Nr);
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

double triangleArea(Vector2d posA, Vector2d posB, Vector2d posC);
double tElementInitialArea(const std::array<GhostNode, 3> &gn);
double tElementArea(const GhostNode &A, const GhostNode &B, const GhostNode &C);

#endif