#include "../src/Mesh/reduction.h"
#include "run/doctest.h"

static void checkReductionRecovery(const Matrix2d &C, const char *label) {
  Matrix2d C_R = Matrix2d::Zero();
  Matrix2d M = Matrix2d::Identity();
  int m1Nr = 0;
  int m2Nr = 0;
  int m3Nr = 0;

  const bool ok = lagrangeReduction(C_R, C, M, m1Nr, m2Nr, m3Nr);
  REQUIRE_MESSAGE(ok, label, " lagrangeReduction did not converge");

  const Matrix2d C_R_test = M.transpose() * C * M;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      CHECK(doctest::Approx(C_R_test(i, j)).epsilon(1e-12) == C_R(i, j));
    }
  }
}

TEST_CASE("Lagrange reduction recovery via M") {
  Matrix2d F;
  F << 1, 1, 0, 1;
  const Matrix2d C = F.transpose() * F;
  checkReductionRecovery(C, "easy");

  Matrix2d A;
  A << 1, -1.4, 0, 1;
  const Matrix2d F2 = F * F * F.transpose() * A;
  const Matrix2d C2 = F2.transpose() * F2;
  checkReductionRecovery(C2, "medium");

  Matrix2d C3;
  C3 << 5252720.5, -12370489., -12370489., 29133284.;
  checkReductionRecovery(C3, "hard");
}
