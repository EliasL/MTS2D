#include "Simulation/sparse_cholesky.h"
#include "Simulation/simulation.h"
#include "Simulation/tangent_stiffness.h"
#include "run/doctest.h"
#include <cmath>

TEST_CASE("Sparse Cholesky backends preserve physical inverse and coordinate duality") {
  for (bool periodic : {false, true}) {
    Mesh mesh(5, 5, 1, 0.01, periodic, "major", "contiSquare", 4);
    if (!periodic) mesh.fixBorderNodes();
    mesh.applyTransformation(getShear(0.08));
    const auto tangent = assembleTangentStiffness(mesh);
    REQUIRE(tangent.valid);
    const int n = tangent.nrActiveDofs;
    Eigen::VectorXd v(n), g(n);
    for (int i = 0; i < n; ++i) {
      v[i] = std::sin(i + 0.2);
      g[i] = std::cos(0.7 * i);
    }
    Eigen::VectorXd reference;
    for (auto backend : {SparseBackend::Eigen, SparseBackend::Cholmod}) {
      if (!SparseCholesky::available(backend)) continue;
      CAPTURE(periodic);
      const std::string backendName = SparseCholesky::name(backend);
      CAPTURE(backendName);
      SparseCholesky factor(backend);
      factor.analyzePattern(tangent.stiffness);
      REQUIRE(factor.factorize(tangent.stiffness, 0));
      const Eigen::VectorXd z = factor.toTransformed(v);
      CHECK((factor.toPhysical(z) - v).norm() < 1e-10 * v.norm());
      const Eigen::VectorXd gz = factor.transformGradient(g);
      CHECK(gz.dot(z) == doctest::Approx(g.dot(v)).epsilon(1e-10));
      const Eigen::VectorXd inverse = factor.toPhysical(gz);
      CHECK((tangent.stiffness * inverse - g).norm() < 1e-10 * g.norm());
      CHECK((factor.solve(g) - inverse).norm() < 1e-10 * inverse.norm());
      if (reference.size() == 0) reference = inverse;
      else CHECK((inverse - reference).norm() < 1e-10 * reference.norm());
      // Reuse the numerical factor and CHOLMOD workspace with a different RHS.
      CHECK((factor.toPhysical(factor.toTransformed(g)) - g).norm() < 1e-10 * g.norm());
    }
  }
}

TEST_CASE("Sparse Cholesky retries an indefinite matrix using the same symbolic analysis") {
  Eigen::SparseMatrix<double> matrix(2, 2);
  std::vector<Eigen::Triplet<double>> entries{{0, 0, 1}, {1, 1, 1},
                                            {0, 1, 2}, {1, 0, 2}};
  matrix.setFromTriplets(entries.begin(), entries.end());
  const Eigen::Vector2d rhs(0.3, -0.7);
  for (auto backend : {SparseBackend::Eigen, SparseBackend::Cholmod}) {
    if (!SparseCholesky::available(backend)) continue;
    const std::string backendName = SparseCholesky::name(backend);
    CAPTURE(backendName);
    SparseCholesky factor(backend);
    factor.analyzePattern(matrix);
    CHECK_FALSE(factor.factorize(matrix, 0));
    CHECK_FALSE(factor.factorize(matrix, 0.25));
    REQUIRE(factor.factorize(matrix, 2));
    const Eigen::VectorXd x = factor.toPhysical(factor.transformGradient(rhs));
    CHECK((matrix * x + 2 * x - rhs).norm() < 1e-12);
    CHECK((factor.toPhysical(factor.toTransformed(rhs)) - rhs).norm() < 1e-12);
  }
}
