#include "Simulation/energyFunctions.h"
#include "run/doctest.h"
#include <algorithm>

namespace {
struct MaterialCase {
  double c11, c22, c12, beta, bulk, noise, energy;
};

// Reference values from the ordinary energy function before sharing its formula.
constexpr MaterialCase cases[] = {
    {1.0, 1.0, 0.0, -0.25, 4.0, 1.0, 3.9116161616161609},
    {1.04, 0.97, 0.12, -0.25, 4.0, 0.83, 4.0076019388134556},
    {0.8, 1.3, -0.2, -0.25, 0.5, 1.17, 0.75542337812042026},
    {1.0, 1.0, 0.0, 4.0, 4.0, 1.0, 4.0404040404040407},
    {1.12, 0.94, 0.23, 4.0, 7.0, 0.91, 7.0670288925938971},
    {1.7, 0.75, -0.31, 4.0, 2.0, 1.09, 2.0235733309590502},
};
} // namespace

TEST_CASE("Shared constitutive energy preserves ordinary and combined evaluation") {
  for (const auto &sample : cases) {
    CAPTURE(sample.beta);
    CAPTURE(sample.noise);
    const double energy = ContiPotential::energyDensity(
        sample.c11, sample.c22, sample.c12, sample.beta, sample.bulk, sample.noise);
    CHECK(energy == doctest::Approx(sample.energy).epsilon(1e-13));

    Matrix2d combinedStress;
    const double combinedEnergy = ContiPotential::energyDensityAndStress(
        sample.c11, sample.c22, sample.c12, sample.beta, sample.bulk, sample.noise,
        &combinedStress);
    CHECK(combinedEnergy == doctest::Approx(sample.energy).epsilon(1e-13));
    const Matrix2d stress = ContiPotential::stress(
        sample.c11, sample.c22, sample.c12, sample.beta, sample.bulk, sample.noise);
    CHECK((stress - combinedStress).norm() < 1e-11 * std::max(1.0, stress.norm()));
  }
}

TEST_CASE("Generated stress differentiates the shared constitutive energy") {
  for (const auto &sample : cases) {
    CAPTURE(sample.beta);
    CAPTURE(sample.noise);
    const Matrix2d stress = ContiPotential::stress(
        sample.c11, sample.c22, sample.c12, sample.beta, sample.bulk, sample.noise);
    // C12 represents both off-diagonal entries of the symmetric metric.
    const Vector3d expected(stress(0, 0), stress(1, 1), 2 * stress(0, 1));
    const Vector3d metric(sample.c11, sample.c22, sample.c12);
    const double h = 1e-6;
    for (int i = 0; i < 3; ++i) {
      Vector3d plus = metric;
      Vector3d minus = metric;
      plus[i] += h;
      minus[i] -= h;
      const double energyPlus = ContiPotential::energyDensity(
          plus[0], plus[1], plus[2], sample.beta, sample.bulk, sample.noise);
      const double energyMinus = ContiPotential::energyDensity(
          minus[0], minus[1], minus[2], sample.beta, sample.bulk, sample.noise);
      CHECK((energyPlus - energyMinus) / (2 * h) ==
            doctest::Approx(expected[i]).epsilon(1e-7).scale(1));
    }
  }
}
