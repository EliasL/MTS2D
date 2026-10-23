#pragma once

#include "tangent_stiffness.h"
#include <Eigen/SparseCholesky>
#include <optimization.h>
#include <cereal/cereal.hpp>

struct PreconditionerStatistics {
  size_t builds = 0;
  size_t reuses = 0;
  size_t fallbacks = 0;
  double buildMilliseconds = 0;

  template <class Archive> void serialize(Archive &ar) {
    ar(cereal::make_nvp("builds", builds), cereal::make_nvp("reuses", reuses),
       cereal::make_nvp("fallbacks", fallbacks),
       cereal::make_nvp("buildMilliseconds", buildMilliseconds));
  }
};

enum class PreconditionerKind { ExactTangent, Diagonal };

// The exact mode uses Sylvain's sparse change of coordinates. Ignoring the
// sparse permutation:
// K + shift I = L L^T, z = L^T (u - origin), grad_z = L^-1 grad_u.
// The factor stays fixed throughout one minimization; the physical energy and
// forces are still evaluated by the usual mesh code.
// The diagonal mode stores only positive stiffness entries for ALGLIB's native
// diagonal preconditioner. It leaves displacement coordinates unchanged.
class StiffnessPreconditioner {
public:
  bool build(const Mesh &mesh, const alglib::real_1d_array &physical,
             PreconditionerKind kind = PreconditionerKind::ExactTangent);
  bool compatibleWith(
      const Mesh &mesh,
      PreconditionerKind kind = PreconditionerKind::ExactTangent) const;
  bool usesTransformedCoordinates() const {
    return kind == PreconditionerKind::ExactTangent;
  }
  const alglib::real_1d_array &diagonalValues() const { return diagonal; }
  int activeDofs() const { return tangent.nrActiveDofs; }
  double diagonalShift() const { return shift; }
  double maxDisplacementSinceBuild(const alglib::real_1d_array &physical) const;

  // Recenter each solve without refactorizing. In particular, a node pinned to
  // remove translations must stay at its current guess, not an older origin.
  void startAt(const alglib::real_1d_array &physical);
  void physicalToTransformed(const alglib::real_1d_array &physical,
                             alglib::real_1d_array &z) const;
  void transformedToPhysical(const alglib::real_1d_array &z,
                             alglib::real_1d_array &physical) const;
  void physicalGradientToTransformed(const alglib::real_1d_array &physicalGrad,
                                     alglib::real_1d_array &gradient) const;

private:
  TangentStiffnessSystem tangent;
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> factor;
  Eigen::VectorXd origin;
  Eigen::VectorXd buildOrigin;
  std::vector<NodeId> freeNodeIds;
  std::vector<NodeId> fixedNodeIds;
  alglib::real_1d_array diagonal;
  PreconditionerKind kind = PreconditionerKind::ExactTangent;
  double shift = 0;
  bool ready = false;
};
