#pragma once

#include "Mesh/mesh.h"
#include <Eigen/Sparse>
#include <vector>

enum class TangentAssembly { Full, Diagonal };

struct TangentStiffnessSystem {
  Eigen::SparseMatrix<double> stiffness;
  Eigen::VectorXd diagonal;
  double referenceVolume = 0;
  int nrFreeNodes = 0;
  int nrFullDofs = 0;
  int nrActiveDofs = 0;
  std::vector<int> fullToActive;
  std::vector<int> activeToFull;
  bool valid = false;
};

// Differentiate the energy on the current smooth branch, with connectivity and
// lattice reductions held fixed. Remove fixed DOFs and, for an unanchored mesh,
// one node's x/y DOFs to remove uniform translations in the full system.
// Diagonal assembly keeps all free DOFs and never constructs a sparse matrix;
// ALGLIB applies that diagonal directly in physical coordinates.
TangentStiffnessSystem assembleTangentStiffness(
    const Mesh &mesh, TangentAssembly assembly = TangentAssembly::Full);
