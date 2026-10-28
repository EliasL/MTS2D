#include "Simulation/sparse_cholesky.h"
#include "Simulation/simulation.h"
#include "Simulation/tangent_stiffness.h"
#include "Simulation/finite.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using Clock = std::chrono::steady_clock;
static double milliseconds(Clock::time_point start) {
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}
static double median(std::vector<double> values) {
  std::sort(values.begin(), values.end());
  return values[values.size() / 2];
}

// Run explicitly in Release; this benchmark is never part of the test suite.
int main(int argc, char **argv) {
  try {
    const int size = argc > 1 ? std::stoi(argv[1]) : 100;
    const int repeats = argc > 2 ? std::stoi(argv[2]) : 3;
    if (size < 4 || repeats < 1) throw std::invalid_argument("Expected size >= 4 and repeats >= 1.");
    Mesh mesh(size, size, 1, 0.01, true, "major", "contiSquare", 4);
    mesh.applyTransformation(getShear(0.08));
    mesh.ensureForces();
    auto start = Clock::now();
    const auto tangent = assembleTangentStiffness(mesh);
    const double assemblyMs = milliseconds(start);
    if (!tangent.valid) throw std::runtime_error("Invalid benchmark tangent.");
    const auto &matrix = tangent.stiffness;
    Eigen::VectorXd rhs(matrix.rows());
    for (int i = 0; i < rhs.size(); ++i) rhs[i] = std::cos(0.7 * i);
    std::cout << "backend,size,dofs,nnz,assembly_ms,analysis_ms,factor_ms,two_half_solves_ms,relative_residual\n";
    Eigen::VectorXd reference;
    for (auto backend : {SparseBackend::Eigen, SparseBackend::Cholmod}) {
      if (!SparseCholesky::available(backend)) continue;
      std::vector<double> analyses, factors, solves;
      double residual = 0;
      for (int run = 0; run < repeats; ++run) {
        SparseCholesky factor(backend);
        start = Clock::now();
        factor.analyzePattern(matrix);
        analyses.push_back(milliseconds(start));
        start = Clock::now();
        if (!factor.factorize(matrix, 0)) throw std::runtime_error("Benchmark tangent is not positive definite.");
        factors.push_back(milliseconds(start));
        Eigen::VectorXd x;
        start = Clock::now();
        for (int solve = 0; solve < 20; ++solve) {
          x = factor.toPhysical(factor.transformGradient(rhs));
        }
        solves.push_back(milliseconds(start) / 20);
        residual = (matrix * x - rhs).norm() / rhs.norm();
        if (!finiteValue(residual) || residual > 1e-7) throw std::runtime_error("Benchmark solve residual is too large.");
        if (reference.size() == 0) reference = x;
        if ((x - reference).norm() > 1e-7 * reference.norm()) throw std::runtime_error("Backends disagree on the physical solution.");
      }
      std::cout << std::setprecision(9) << SparseCholesky::name(backend) << ','
                << size << ',' << matrix.rows() << ',' << matrix.nonZeros() << ','
                << assemblyMs << ',' << median(analyses) << ',' << median(factors)
                << ',' << median(solves) << ',' << residual << '\n';
    }
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
