#pragma once

#include <Eigen/Sparse>
#include <memory>

enum class SparseBackend { Eigen, Cholmod };

// P (K + shift I) P^T = L L^T. Both backends read the lower triangle.
// Analyze once, then refactorize matrices with that same sparsity pattern.
class SparseCholesky {
public:
  static bool available(SparseBackend backend);
  static SparseBackend defaultBackend();
  static const char *name(SparseBackend backend);
  explicit SparseCholesky(SparseBackend backend = defaultBackend());
  ~SparseCholesky();
  SparseCholesky(const SparseCholesky &) = delete;
  SparseCholesky &operator=(const SparseCholesky &) = delete;

  void analyzePattern(const Eigen::SparseMatrix<double> &matrix);
  bool factorize(const Eigen::SparseMatrix<double> &matrix, double shift);
  Eigen::VectorXd toTransformed(const Eigen::VectorXd &physical) const; // L^T P x
  Eigen::VectorXd toPhysical(const Eigen::VectorXd &z) const; // P^T L^-T z
  Eigen::VectorXd transformGradient(const Eigen::VectorXd &g) const; // L^-1 P g
  Eigen::VectorXd solve(const Eigen::VectorXd &rhs) const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl;
};
