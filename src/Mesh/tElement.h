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
  std::array<GhostNode, 3> tElementNodes;

  // Deformation gradient
  Matrix2d F;

  // Metric tensor (C = F^TF)
  Matrix2d C;

  // Reduced metric tensor
  Matrix2d C_;

  // Reduction transformation matrix (m^TCm = C_)
  Matrix2d m;
  Matrix2d simple_m;

  // Reduced stress
  Matrix2d dPhi_dC_;

  // First Piola-Kirchhoff stress tensor, representing the stress in the
  // reference configuration.
  Matrix2d P;

  // Strain energy of the cell, representing the potential energy stored due
  // to deformation.
  double energy = 0;

  // A representation of stress that is unaffected by the directionality of
  // loading. Discontinuous yielding of pristine micro-crystals (page 216/17)
  double resolvedShearStress = 0;

  // The jacobian of the shapefunction with respect to the reference state.
  //  We use this to calculate the deformation gradiant F
  // ∂ξ/∂X
  Matrix2d dxi_dX;
  // Derivatives place holders
  Matrix2d du_dxi;
  Matrix2d dX_dxi;

  // These are adjustment vectors that we multiply together with the piola
  // tensor to correctly extract the force corresponding to each node.
  // Similarly to dxi_dX, these only update once, during initialization.
  std::array<Vector2d, 3> dN_dX;

  // A flag to indicate whether or not a plastic event has occured
  bool plasticChange = false;

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

  // For completeness, we keep track of m1 and m2 as well
  int m1Nr = 0;
  int m2Nr = 0;

  // Index of element. Used for debugging.
  int eIndex;

  // A noise value which will slightly distort the volumetric energy term of the
  // element.
  double noise = 1;

private:
  /*
  Shape functions:
  N1 = 1 - ξ1 - ξ2
  N2 = ξ1
  N3 = ξ2

  Derivative of shape functions:
  For simplicity of implementation, these derivatives are assumed to be
  constant! If you change b1, b2 or b3, you will also need to manually change
  the implementation of the jacobian calculation. See du_dxi.

  dN_dxi =
  -1.0, 1.0, 0.0,
  -1.0, 0.0, 1.0

  */
  static Matrix<double, 2, 3> dN_dxi;

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
  double initArea = 0;

  // A variable to store the ground state energy to set our ground state energy
  // to be zero
  double groundStateEnergyDensity = 0;

public:
  // Constructor for the triangular element. Initializes the 3 defining
  // GhostNodes and calculates the inverse state A_inv, to later be used in
  // calculating F.
  TElement(Mesh &mesh, GhostNode n1, GhostNode n2, GhostNode n3,
           int elementIndex, double noise = 1);
  TElement() {};

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
  static double calculateEnergyDensity(double c11, double c22, double c12);

  // Used for testing the lagrange reuction functions
  static TElement lagrangeReduction(double c11, double c22, double c12);

private:
  // Calculate the Jacobian with respect to the displacement of the nodes
  Matrix2d m_update_du_dxi();
  // Calculate the Jacobian with respect to the initial position of the nodes
  Matrix2d m_update_dX_dxi();

  // Copy the displacement from the real nodes to the nodes in the element
  void m_updatePosition(const Mesh &mesh);

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

  // Performs a check to see if the previous lagrange reduction still works
  bool m_checkIfPreviousReductionWorks();

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

  // Before, I used to serialize the elements as well, but they can be
  // reconstructed from the nodes. That will be usefull later anyway.

  friend class cereal::access;
  template <class Archive> void serialize(Archive &ar) {
    ar(MAKE_NVP(tElementNodes), MAKE_NVP(F), MAKE_NVP(C), MAKE_NVP(C_),
       MAKE_NVP(m), MAKE_NVP(dPhi_dC_), MAKE_NVP(P), MAKE_NVP(energy),
       MAKE_NVP(resolvedShearStress), MAKE_NVP(dxi_dX), MAKE_NVP(du_dxi),
       MAKE_NVP(dX_dxi), MAKE_NVP(dN_dX), MAKE_NVP(plasticChange),
       MAKE_NVP(m3Nr), MAKE_NVP(pastM3Nr), MAKE_NVP(m1Nr), MAKE_NVP(m2Nr),
       MAKE_NVP(eIndex), MAKE_NVP(simple_m), MAKE_NVP(noise),
       MAKE_NVP(initArea), MAKE_NVP(groundStateEnergyDensity));
  }

  // Giving access to private variables
  friend bool compareTElementsInternal(const TElement &lhs, const TElement &rhs,
                                       std::string *debugMsg, int tabNumber);
};

// This function updates the list of connected elements in the real nodes
void addElementIndices(Mesh &mesh, const std::vector<GhostNode> nodeList,
                       int elementIndex);

std::ostream &operator<<(std::ostream &os, const TElement &element);

inline bool compareTElementsInternal(const TElement &lhs, const TElement &rhs,
                                     std::string *debugMsg = nullptr,
                                     int tabNumber = 0) {
  bool equal = true;

  // Compare public members.
  COMPARE_FIELD(tElementNodes);
  COMPARE_FIELD(F);
  COMPARE_FIELD(C);
  COMPARE_FIELD(C_);
  COMPARE_FIELD(m);
  COMPARE_FIELD(simple_m);
  COMPARE_FIELD(dPhi_dC_);
  COMPARE_FIELD(P);
  COMPARE_FIELD(energy);
  COMPARE_FIELD(resolvedShearStress);
  COMPARE_FIELD(dxi_dX);
  COMPARE_FIELD(du_dxi);
  COMPARE_FIELD(dX_dxi);
  COMPARE_FIELD(dN_dX);
  COMPARE_FIELD(plasticChange);
  COMPARE_FIELD(m3Nr);
  COMPARE_FIELD(pastM3Nr);
  COMPARE_FIELD(m1Nr);
  COMPARE_FIELD(m2Nr);
  COMPARE_FIELD(eIndex);
  COMPARE_FIELD(noise);

  // Compare private members.
  COMPARE_FIELD(initArea);
  COMPARE_FIELD(groundStateEnergyDensity);
  COMPARE_FIELD(dN_dxi);
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
  for (size_t i = 0; i < lhs.tElementNodes.size(); i++) {
    if (!(lhs.tElementNodes[i] == rhs.tElementNodes[i])) {
      diff += std::string(tabNumber, '\t') + "tElementNodes[" +
              std::to_string(i) + "] differs -> \n";
      // Recursively call debugCompare for TElement
      diff += debugCompare(lhs.tElementNodes[i], rhs.tElementNodes[i],
                           tabNumber + 1);
    }
  }
  compareTElementsInternal(lhs, rhs, &diff, tabNumber);
  return diff;
}

double tElementInitialArea(const GhostNode &A, const GhostNode &B,
                           const GhostNode &C);

#endif