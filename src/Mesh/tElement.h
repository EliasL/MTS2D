#ifndef TELEMENT_H
#define TELEMENT_H
#include "Eigen/Core"
#pragma once

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
  std::array<Node, 3> nodes;

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
  Matrix2d r_s;

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

  // These are adjustment vectors that we multiply together with the piola
  // tensor to correctly extract the force corresponding to each node.
  // Similarly to dxi_dX, these only update once, during initialization.
  std::array<Vector2d, 3> r;

  // A flag to indicate whether or not a plastic event has occured
  bool plasticChange = false;

  // We only save data when plasticity occurs, so we keep a reference of
  // how many times m3 is applied in the lagrange reduction. If this number
  // changes from one reduction to another, we know that a plastic event has
  // occured. (ie. the energy potential suddenly has a gradient in a new
  // direction, ie. the node has fallen into a different energy well.)
  int m3Nr = 0;
  // This keeps track of the number of m3 shears in the previous lagrange
  // reduction
  int pastM3Nr = 0;

  // For completeness, we keep track of m1 and m2 as well
  int m1Nr = 0;
  int m2Nr = 0;

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

  b =
  -1.0, 1.0, 0.0,
  -1.0, 0.0, 1.0

  */
  static Matrix<double, 2, 3> b;

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

  // A noise value which will slightly distort the volumetric energy term of the
  // element.
  double noise = 1;

  // A variable to store the ground state energy to set our ground state energy
  // to be zero
  double groundStateEnergyDensity = 0;

public:
  // Constructor for the triangular element. Initializes the 3 defining nodes
  // and calculates the inverse state A_inv, to later be used in calculating F.
  TElement(Node n1, Node n2, Node n3, double noise = 1);
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

  void update(Mesh &mesh);

  // Usefull if you only care about the energy given the C matrix.
  static double calculateEnergyDensity(double c11, double c22, double c12);

  // Used for testing the lagrange reuction functions
  static TElement lagrangeReduction(double c11, double c22, double c12);

  // Updates the past number of m3 steps. This should be done in the simulation
  // as this is where we keep track of different loading steps
  void updatePastM3Nr();

private:
  Matrix2d du_dxi;
  Matrix2d dX_dxi;
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

  // Calculates energy
  void m_updateEnergy();

  // Calculate reduced stress
  void m_updateReducedStress();

  // Calculate Piola stress P
  void m_updatePiolaStress();

  // Calculate the resolved-shear stress
  void m_updateResolvedShearStress();

  // Calculate force on each node
  void m_updateForceOnEachNode();

  // Position subtraction (The vector from node 1 to node 2)
  Vector2d const dx(const Node &n1, const Node &n2) const {
    return n2.pos() - n1.pos();
  }

  // Initial-position subtraction
  Vector2d const dX(const Node &n1, const Node &n2) const {
    return n2.init_pos() - n1.init_pos();
  }

  // Displacement subtraction
  Vector2d const du(const Node &n1, const Node &n2) const {
    return n2.u() - n1.u();
  }

  friend class cereal::access;
  template <class Archive> void serialize(Archive &ar) {
    ar(nodes, F, C, C_, m, r_s, P, energy, resolvedShearStress, dxi_dX, r,
       plasticChange, m3Nr, pastM3Nr, m1Nr, m2Nr, simple_m, noise, initArea,
       groundStateEnergyDensity);
  }
};

std::ostream &operator<<(std::ostream &os, const TElement &element);

#endif