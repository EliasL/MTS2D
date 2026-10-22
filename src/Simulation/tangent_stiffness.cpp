// Tangent assembly adapted from Sylvain's tangent-preconditioning implementation.
#include "Simulation/tangent_stiffness.h"
#include "Simulation/energyFunctions.h"
#include "finite.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

using Eigen::Matrix2d;
using Eigen::Matrix3d;
using Eigen::SparseMatrix;
using Eigen::Triplet;
using Eigen::Vector3d;
using Eigen::VectorXd;

namespace {

// Tiny second-order automatic-differentiation type for the reduced metric
// variables (C11, C22, C12).  It gives the gradient and Hessian of one
// element's energy density without finite differences.
struct SecondOrderScalar {
  double value = 0.0;
  Vector3d gradient = Vector3d::Zero();
  Matrix3d hessian = Matrix3d::Zero();

  SecondOrderScalar() = default;
  SecondOrderScalar(double value_) : value(value_) {}

  static SecondOrderScalar variable(double value, int derivativeIndex) {
    SecondOrderScalar out(value);
    out.gradient[derivativeIndex] = 1.0;
    return out;
  }
};

SecondOrderScalar operator+(const SecondOrderScalar &a,
                            const SecondOrderScalar &b) {
  SecondOrderScalar out;
  out.value = a.value + b.value;
  out.gradient = a.gradient + b.gradient;
  out.hessian = a.hessian + b.hessian;
  return out;
}

SecondOrderScalar operator-(const SecondOrderScalar &a,
                            const SecondOrderScalar &b) {
  SecondOrderScalar out;
  out.value = a.value - b.value;
  out.gradient = a.gradient - b.gradient;
  out.hessian = a.hessian - b.hessian;
  return out;
}

SecondOrderScalar operator*(const SecondOrderScalar &a,
                            const SecondOrderScalar &b) {
  SecondOrderScalar out;
  out.value = a.value * b.value;
  out.gradient = a.gradient * b.value + b.gradient * a.value;
  out.hessian = a.hessian * b.value + b.hessian * a.value +
                a.gradient * b.gradient.transpose() +
                b.gradient * a.gradient.transpose();
  return out;
}

SecondOrderScalar pow(const SecondOrderScalar &x, double exponent) {
  SecondOrderScalar out;
  out.value = std::pow(x.value, exponent);
  const double first = exponent * std::pow(x.value, exponent - 1.0);
  const double second =
      exponent * (exponent - 1.0) * std::pow(x.value, exponent - 2.0);
  out.gradient = first * x.gradient;
  out.hessian = first * x.hessian + second * x.gradient * x.gradient.transpose();
  return out;
}

SecondOrderScalar log(const SecondOrderScalar &x) {
  SecondOrderScalar out;
  out.value = std::log(x.value);
  out.gradient = x.gradient / x.value;
  out.hessian =
      x.hessian / x.value -
      (x.gradient * x.gradient.transpose()) / (x.value * x.value);
  return out;
}

struct MetricEnergyDerivatives {
  Vector3d gradient = Vector3d::Zero();
  Matrix3d hessian = Matrix3d::Zero();
};

// Differentiate the shared constitutive energy in reduced-metric variables,
// then assemble the nodal tangent through F -> C -> M_l^T C M_l.
MetricEnergyDerivatives metricEnergyDerivatives(double C11, double C22,
                                                double C12, double beta,
                                                double K, double noise) {
  const SecondOrderScalar c11 = SecondOrderScalar::variable(C11, 0);
  const SecondOrderScalar c22 = SecondOrderScalar::variable(C22, 1);
  const SecondOrderScalar c12 = SecondOrderScalar::variable(C12, 2);

  const SecondOrderScalar phi =
      ContiPotential::energyDensityFormula(c11, c22, c12, beta, K, noise);

  MetricEnergyDerivatives out;
  out.gradient = phi.gradient;
  out.hessian = phi.hessian;
  return out;
}

Vector3d reducedMetricVariation(const Matrix2d &F, const Matrix2d &M_l,
                                const Matrix2d &dF) {
  const Matrix2d dC = dF.transpose() * F + F.transpose() * dF;
  const Matrix2d dC_R = M_l.transpose() * dC * M_l;
  return Vector3d(dC_R(0, 0), dC_R(1, 1), dC_R(0, 1));
}

Vector3d reducedMetricSecondVariation(const Matrix2d &M_l,
                                      const Matrix2d &dF_a,
                                      const Matrix2d &dF_b) {
  const Matrix2d d2C = dF_a.transpose() * dF_b + dF_b.transpose() * dF_a;
  const Matrix2d d2C_R = M_l.transpose() * d2C * M_l;
  return Vector3d(d2C_R(0, 0), d2C_R(1, 1), d2C_R(0, 1));
}

Matrix2d nodalDeformationGradientDerivative(const TElement &element,
                                            int localNode,
                                            int component) {
  Matrix2d dF = Matrix2d::Zero();
  dF.row(component) = element.dN_dX.row(localNode);
  return dF;
}


} // namespace

// Second derivative of phi(F^T F) along F = I + gamma e_x e_y^T.
double referenceShearModulus(double beta, double bulkModulus) {
  const auto derivatives = metricEnergyDerivatives(1, 1, 0, beta, bulkModulus, 1);
  return 2 * derivatives.gradient[1] + derivatives.hessian(2, 2);
}

TangentStiffnessSystem assembleTangentStiffness(const Mesh &inputMesh,
                                               TangentAssembly assembly) {
  TangentStiffnessSystem result;
  const bool diagonalOnly = assembly == TangentAssembly::Diagonal;

  Mesh scratch;
  const Mesh *mesh = &inputMesh;
  if (inputMesh.getUpdateState() == Mesh::UpdateState::Dirty) {
    scratch = inputMesh;
    scratch.ensureForces();
    mesh = &scratch;
  }

  result.nrFreeNodes = static_cast<int>(mesh->freeNodeIds.size());
  result.nrFullDofs = 2 * result.nrFreeNodes;
  result.fullToActive.assign(static_cast<size_t>(result.nrFullDofs), -1);

  std::vector<int> freeIndexFromNode(static_cast<size_t>(mesh->nrNodes), -1);
  for (int i = 0; i < result.nrFreeNodes; ++i) {
    freeIndexFromNode[static_cast<size_t>(mesh->freeNodeIds[i].i)] = i;
  }

  std::vector<char> removeDof(static_cast<size_t>(result.nrFullDofs), 0);
  if (!diagonalOnly && mesh->fixedNodeIds.empty() && result.nrFreeNodes > 0) {
    // Periodic/free systems have two translational zero modes.  Removing one
    // real node's x/y DOFs makes the tangent factorizable while leaving all
    // relative displacements represented.
    removeDof[0] = 1;
    removeDof[static_cast<size_t>(result.nrFreeNodes)] = 1;
  }

  for (int i = 0; i < result.nrFullDofs; ++i) {
    if (!removeDof[static_cast<size_t>(i)]) {
      result.fullToActive[static_cast<size_t>(i)] = result.nrActiveDofs++;
      result.activeToFull.push_back(i);
    }
  }

  std::vector<Triplet<double>> triplets;
  if (diagonalOnly) {
    result.diagonal = VectorXd::Zero(result.nrActiveDofs);
  } else {
    triplets.reserve(static_cast<size_t>(mesh->nrElements) * 36);
  }

  for (const TElement &element : mesh->elements) {
    const double area = element.getInitArea();
    result.referenceVolume += area;

    const MetricEnergyDerivatives metricDerivatives = metricEnergyDerivatives(
        element.C_R(0, 0), element.C_R(1, 1), element.C_R(0, 1), element.beta,
        element.K, element.noise);

    std::array<Matrix2d, 6> dF_dof;
    std::array<Vector3d, 6> dc_dof;
    std::array<int, 6> activeDofByLocal;

    for (int localNode = 0; localNode < 3; ++localNode) {
      const int freeIndex = freeIndexFromNode[static_cast<size_t>(
          element.ghostNodes[localNode].referenceId.i)];
      for (int component = 0; component < 2; ++component) {
        const int localDof = 2 * localNode + component;
        dF_dof[static_cast<size_t>(localDof)] =
            nodalDeformationGradientDerivative(element, localNode, component);
        dc_dof[static_cast<size_t>(localDof)] = reducedMetricVariation(
            element.F, element.M_l, dF_dof[static_cast<size_t>(localDof)]);

        int activeDof = -1;
        if (freeIndex >= 0) {
          const int fullDof =
              component == 0 ? freeIndex : result.nrFreeNodes + freeIndex;
          activeDof = result.fullToActive[static_cast<size_t>(fullDof)];
        }
        activeDofByLocal[static_cast<size_t>(localDof)] = activeDof;
      }
    }

    for (int p = 0; p < 6; ++p) {
      const int activeP = activeDofByLocal[static_cast<size_t>(p)];
      if (activeP < 0) {
        continue;
      }

      for (int q = 0; q < 6; ++q) {
        const int activeQ = activeDofByLocal[static_cast<size_t>(q)];
        if (activeQ < 0) {
          continue;
        }
        if (diagonalOnly && activeP != activeQ) {
          continue;
        }

        const Vector3d d2c_pq = reducedMetricSecondVariation(
            element.M_l, dF_dof[static_cast<size_t>(p)],
            dF_dof[static_cast<size_t>(q)]);
        // Element Hessian block d^2E/(du_p du_q), assembled into the active
        // DOF system after fixed nodes and translational pins have been removed.
        const double stiffness =
            area *
            (dc_dof[static_cast<size_t>(p)].dot(
                 metricDerivatives.hessian *
                 dc_dof[static_cast<size_t>(q)]) +
             metricDerivatives.gradient.dot(d2c_pq));
        if (diagonalOnly) {
          // Compare global indices: on tiny periodic meshes, different local
          // corners can refer to the same physical displacement component.
          result.diagonal[activeP] += stiffness;
        } else {
          triplets.emplace_back(activeP, activeQ, stiffness);
        }
      }
    }
  }

  if (!finiteValue(result.referenceVolume) || result.referenceVolume <= 0.0) {
    return result;
  }
  if (diagonalOnly) {
    for (int i = 0; i < result.diagonal.size(); ++i) {
      if (!finiteValue(result.diagonal[i])) {
        return result;
      }
    }
    result.valid = true;
    return result;
  }
  for (const auto &triplet : triplets) {
    if (!finiteValue(triplet.value())) {
      return result;
    }
  }

  result.stiffness =
      SparseMatrix<double>(result.nrActiveDofs, result.nrActiveDofs);
  result.stiffness.setFromTriplets(triplets.begin(), triplets.end());
  result.stiffness.makeCompressed();
  result.valid = true;
  return result;
}
