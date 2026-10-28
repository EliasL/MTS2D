#include "stiffness_preconditioner.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {
Eigen::VectorXd asVector(const alglib::real_1d_array &array) {
  return Eigen::Map<const Eigen::VectorXd>(array.getcontent(), array.length());
}
} // namespace

bool StiffnessPreconditioner::build(const Mesh &mesh,
                                  const alglib::real_1d_array &physical,
                                  PreconditionerKind requestedKind) {
  ready = false;
  kind = requestedKind;
  shift = 0;
  tangent = assembleTangentStiffness(
      mesh, kind == PreconditionerKind::Diagonal ? TangentAssembly::Diagonal
                                               : TangentAssembly::Full);
  if (!tangent.valid || tangent.nrActiveDofs == 0 ||
      physical.length() != tangent.nrFullDofs) {
    return false;
  }

  double scale = 1;
  for (int i = 0; i < tangent.nrActiveDofs; ++i) {
    const double value = kind == PreconditionerKind::Diagonal
                             ? tangent.diagonal[i]
                             : tangent.stiffness.coeff(i, i);
    scale = std::max(scale, std::abs(value));
  }
  if (kind == PreconditionerKind::Diagonal) {
    // ALGLIB requires a positive Hessian diagonal, not its inverse. Floor small
    // or negative entries so the native preconditioner remains positive.
    const double floor = 1e-8 * scale;
    diagonal.setlength(tangent.nrFullDofs);
    for (int i = 0; i < tangent.nrFullDofs; ++i) {
      diagonal[i] = std::max(tangent.diagonal[i], floor);
    }
    ready = true;
  } else {
    factor.analyzePattern(tangent.stiffness);
    // Preserve Sylvain's shift sequence. If none yields a positive definite
    // factor, the caller uses ordinary optimization for this solve.
    for (double relativeShift : {0.0, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2}) {
      shift = relativeShift * scale;
      if (factor.factorize(tangent.stiffness, shift)) {
        ready = true;
        break;
      }
    }
  }
  if (!ready) {
    return false;
  }
  freeNodeIds = mesh.freeNodeIds;
  fixedNodeIds = mesh.fixedNodeIds;
  buildOrigin = asVector(physical);
  startAt(physical);
  return true;
}

bool StiffnessPreconditioner::compatibleWith(
    const Mesh &mesh, PreconditionerKind requestedKind) const {
  return ready && kind == requestedKind && freeNodeIds == mesh.freeNodeIds &&
         fixedNodeIds == mesh.fixedNodeIds;
}

void StiffnessPreconditioner::startAt(const alglib::real_1d_array &physical) {
  if (!ready || physical.length() != tangent.nrFullDofs) {
    throw std::invalid_argument("Invalid preconditioner origin.");
  }
  origin = asVector(physical);
}

void StiffnessPreconditioner::physicalToTransformed(
    const alglib::real_1d_array &physical, alglib::real_1d_array &z) const {
  if (!ready || !usesTransformedCoordinates() ||
      physical.length() != tangent.nrFullDofs) {
    throw std::invalid_argument("Invalid physical preconditioner coordinates.");
  }
  Eigen::VectorXd active(tangent.nrActiveDofs);
  for (int i = 0; i < tangent.nrActiveDofs; ++i) {
    const int full = tangent.activeToFull[i];
    active[i] = physical[full] - origin[full];
  }
  const Eigen::VectorXd transformed = factor.toTransformed(active);
  z.setlength(tangent.nrActiveDofs);
  for (int i = 0; i < tangent.nrActiveDofs; ++i) {
    z[i] = transformed[i];
  }
}

void StiffnessPreconditioner::transformedToPhysical(
    const alglib::real_1d_array &z, alglib::real_1d_array &physical) const {
  if (!ready || !usesTransformedCoordinates() ||
      z.length() != tangent.nrActiveDofs) {
    throw std::invalid_argument("Invalid transformed preconditioner coordinates.");
  }
  const Eigen::VectorXd active = factor.toPhysical(asVector(z));
  physical.setlength(tangent.nrFullDofs);
  for (int i = 0; i < tangent.nrFullDofs; ++i) {
    physical[i] = origin[i];
  }
  for (int i = 0; i < tangent.nrActiveDofs; ++i) {
    const int full = tangent.activeToFull[i];
    physical[full] += active[i];
  }
}

void StiffnessPreconditioner::physicalGradientToTransformed(
    const alglib::real_1d_array &physicalGrad,
    alglib::real_1d_array &gradient) const {
  if (!ready || !usesTransformedCoordinates() ||
      physicalGrad.length() != tangent.nrFullDofs) {
    throw std::invalid_argument("Invalid physical preconditioner gradient.");
  }
  Eigen::VectorXd active(tangent.nrActiveDofs);
  for (int i = 0; i < tangent.nrActiveDofs; ++i) {
    active[i] = physicalGrad[tangent.activeToFull[i]];
  }
  // This is the chain rule for the changed coordinates, not an approximation
  // to the physical gradient. Together with the inverse coordinate map it
  // gives the optimizer a stiffness-aware physical search direction.
  active = factor.transformGradient(active);
  // ALGLIB passes a fixed-size proxy here; resizing that proxy is forbidden.
  if (gradient.length() != tangent.nrActiveDofs) {
    gradient.setlength(tangent.nrActiveDofs);
  }
  for (int i = 0; i < tangent.nrActiveDofs; ++i) {
    gradient[i] = active[i];
  }
}

double StiffnessPreconditioner::maxDisplacementSinceBuild(
    const alglib::real_1d_array &physical) const {
  Eigen::Vector2d mean = Eigen::Vector2d::Zero();
  if (fixedNodeIds.empty()) {
    for (int i = 0; i < tangent.nrFreeNodes; ++i) {
      mean.x() += physical[i] - buildOrigin[i];
      mean.y() += physical[tangent.nrFreeNodes + i] -
                  buildOrigin[tangent.nrFreeNodes + i];
    }
    mean /= tangent.nrFreeNodes;
  }
  double maximum = 0;
  for (int i = 0; i < tangent.nrFreeNodes; ++i) {
    const double dx = physical[i] - buildOrigin[i] - mean.x();
    const double dy = physical[tangent.nrFreeNodes + i] -
                      buildOrigin[tangent.nrFreeNodes + i] - mean.y();
    maximum = std::max(maximum, std::sqrt(dx * dx + dy * dy));
  }
  return maximum;
}
