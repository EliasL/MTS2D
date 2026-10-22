#include "Simulation/simulation.h"
#include "Simulation/stiffness_preconditioner.h"
#include "run/doctest.h"
#include <cmath>

namespace {
alglib::real_1d_array displacements(const Mesh &mesh) {
  const int count = mesh.freeNodeIds.size();
  alglib::real_1d_array u;
  u.setlength(2 * count);
  for (int i = 0; i < count; ++i) {
    const Vector2d value = mesh[mesh.freeNodeIds[i]]->u();
    u[i] = value.x();
    u[count + i] = value.y();
  }
  return u;
}

Eigen::VectorXd gradient(Mesh &mesh, const alglib::real_1d_array &u) {
  mesh.updateNodePositions(u.getcontent(), u.length());
  mesh.ensureForces();
  const int count = mesh.freeNodeIds.size();
  Eigen::VectorXd result(2 * count);
  for (int i = 0; i < count; ++i) {
    const auto &node = *mesh[mesh.freeNodeIds[i]];
    result[i] = -node.f.x();
    result[count + i] = -node.f.y();
  }
  return result;
}
} // namespace

TEST_CASE("Tangent preconditioning defaults on and validates configuration") {
  CHECK(initializeConfig({}).preconditioner == "exactTangent");
  CHECK(initializeConfig({{"preconditioner", "none"}}).preconditioner == "none");
  CHECK(initializeConfig({{"preconditioner", "diagonal"}}).preconditioner ==
        "diagonal");
  CHECK_THROWS_AS(initializeConfig({{"preconditioner", "unknown"}}),
                  std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"preconditionerUpdate", "unknown"}}),
                  std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"preconditionerDisplacementThreshold", "0"}}),
                  std::invalid_argument);
}

TEST_CASE("Assembled tangent matches finite differences of physical forces") {
  for (bool periodic : {false, true}) {
    for (const auto &potential : {"contiSquare", "contiTriangular"}) {
      Mesh mesh(5, 5, 1, 0.01, periodic, "major", potential, 4);
      if (!periodic) {
        mesh.fixBorderNodes();
      }
      mesh.applyTransformation(getShear(0.08));
      mesh.ensureForces();
      const auto system = assembleTangentStiffness(mesh);
      REQUIRE(system.valid);
      const auto origin = displacements(mesh);
      Eigen::VectorXd direction(system.nrActiveDofs);
      auto plus = origin;
      auto minus = origin;
      const double h = 1e-6;
      for (int i = 0; i < system.nrActiveDofs; ++i) {
        direction[i] = std::sin(i + 0.3);
        const int full = system.activeToFull[i];
        plus[full] += h * direction[i];
        minus[full] -= h * direction[i];
      }
      const Eigen::VectorXd gPlus = gradient(mesh, plus);
      const Eigen::VectorXd gMinus = gradient(mesh, minus);
      Eigen::VectorXd difference(system.nrActiveDofs);
      for (int i = 0; i < system.nrActiveDofs; ++i) {
        const int full = system.activeToFull[i];
        difference[i] = (gPlus[full] - gMinus[full]) / (2 * h);
      }
      const Eigen::VectorXd predicted = system.stiffness * direction;
      CHECK((predicted - difference).norm() / predicted.norm() < 1e-5);
    }
  }
}

TEST_CASE("Diagonal assembly matches the full tangent without its translation pins") {
  for (bool periodic : {false, true}) {
    for (const auto &potential : {"contiSquare", "contiTriangular"}) {
      Mesh mesh(5, 5, 1, 0.01, periodic, "major", potential, 4);
      if (!periodic) {
        mesh.fixBorderNodes();
      }
      mesh.applyTransformation(getShear(0.08));
      mesh.ensureForces();
      const auto full = assembleTangentStiffness(mesh);
      const auto diagonal =
          assembleTangentStiffness(mesh, TangentAssembly::Diagonal);
      REQUIRE(full.valid);
      REQUIRE(diagonal.valid);
      CHECK(diagonal.stiffness.rows() == 0);
      CHECK(diagonal.nrActiveDofs == diagonal.nrFullDofs);
      for (int i = 0; i < full.nrActiveDofs; ++i) {
        CHECK(diagonal.diagonal[full.activeToFull[i]] ==
              doctest::Approx(full.stiffness.coeff(i, i)).epsilon(1e-12));
      }
    }
  }
}

TEST_CASE("Diagonal assembly includes coincident periodic element corners") {
  Mesh mesh(1, 2, 1, 0, true, "major", "contiSquare", 4);
  mesh.applyTransformation(getShear(0.08));
  mesh.ensureForces();
  const auto diagonal =
      assembleTangentStiffness(mesh, TangentAssembly::Diagonal);
  REQUIRE(diagonal.valid);
  const auto origin = displacements(mesh);
  const double h = 1e-6;
  for (int i = 0; i < origin.length(); ++i) {
    auto plus = origin;
    auto minus = origin;
    plus[i] += h;
    minus[i] -= h;
    const auto gPlus = gradient(mesh, plus);
    const auto gMinus = gradient(mesh, minus);
    CHECK(diagonal.diagonal[i] ==
          doctest::Approx((gPlus[i] - gMinus[i]) / (2 * h)).epsilon(1e-5));
  }
}

TEST_CASE("Native diagonal preconditioner floors nonpositive curvature") {
  Mesh mesh(3, 3, 1, 0, false, "major", "contiSquare", 4);
  mesh.fixBorderNodes();
  for (auto &element : mesh.elements) {
    element.K = -100;
  }
  mesh.markDirty();
  mesh.ensureForces();
  const auto diagonal =
      assembleTangentStiffness(mesh, TangentAssembly::Diagonal);
  REQUIRE(diagonal.valid);
  REQUIRE(diagonal.diagonal.minCoeff() < 0);
  StiffnessPreconditioner preconditioner;
  REQUIRE(preconditioner.build(mesh, displacements(mesh),
                               PreconditionerKind::Diagonal));
  CHECK_FALSE(preconditioner.usesTransformedCoordinates());
  for (int i = 0; i < diagonal.nrFullDofs; ++i) {
    CHECK(preconditioner.diagonalValues()[i] > 0);
  }
}

TEST_CASE("Transformed gradient differentiates the unchanged physical energy") {
  Mesh mesh(5, 5, 1, 0.01, true, "major", "contiSquare", 4);
  mesh.applyTransformation(getShear(0.08));
  mesh.ensureForces();
  auto origin = displacements(mesh);
  StiffnessPreconditioner preconditioner;
  REQUIRE(preconditioner.build(mesh, origin));
  REQUIRE(preconditioner.diagonalShift() == 0);
  alglib::real_1d_array z;
  preconditioner.physicalToTransformed(origin, z);
  CHECK(z.length() == origin.length() - 2);
  for (int i = 0; i < z.length(); ++i) {
    z[i] = 1e-3 * std::sin(i + 0.2);
  }
  alglib::real_1d_array physical;
  preconditioner.transformedToPhysical(z, physical);
  alglib::real_1d_array roundTrip;
  preconditioner.physicalToTransformed(physical, roundTrip);
  for (int i = 0; i < z.length(); ++i) {
    CHECK(roundTrip[i] == doctest::Approx(z[i]).epsilon(1e-9));
  }
  const auto g = gradient(mesh, physical);
  alglib::real_1d_array physicalGradient;
  physicalGradient.setcontent(g.size(), g.data());
  alglib::real_1d_array transformedGradient;
  preconditioner.physicalGradientToTransformed(physicalGradient,
                                              transformedGradient);
  const double h = 1e-5;
  auto plus = z;
  auto minus = z;
  double predicted = 0;
  for (int i = 0; i < z.length(); ++i) {
    const double direction = std::cos(i + 0.1);
    plus[i] += h * direction;
    minus[i] -= h * direction;
    predicted += transformedGradient[i] * direction;
  }
  preconditioner.transformedToPhysical(plus, physical);
  gradient(mesh, physical);
  const double energyPlus = mesh.totalEnergy;
  preconditioner.transformedToPhysical(minus, physical);
  gradient(mesh, physical);
  const double energyMinus = mesh.totalEnergy;
  CHECK((energyPlus - energyMinus) / (2 * h) ==
        doctest::Approx(predicted).epsilon(1e-5).scale(1e-7));

  // Reusing the factor must not move the periodic translation anchor back to
  // the displacement at which the factor was first built.
  origin[0] += 0.2;
  origin[origin.length() / 2] -= 0.1;
  preconditioner.startAt(origin);
  preconditioner.physicalToTransformed(origin, z);
  preconditioner.transformedToPhysical(z, physical);
  for (int i = 0; i < origin.length(); ++i) {
    CHECK(physical[i] == origin[i]);
  }
}

TEST_CASE("Preconditioned and ordinary minimizers find the same elastic state") {
  for (bool periodic : {false, true}) {
    for (const auto &minimizer : {"LBFGS", "CG"}) {
      const std::string modes[] = {"none", "diagonal", "exactTangent"};
      double energies[3];
      for (int mode = 0; mode < 3; ++mode) {
        INFO("periodic=", periodic, ", minimizer=", std::string(minimizer),
             ", mode=", modes[mode]);
        Config config;
        config.setDefaultValues();
        config.rows = 7;
        config.cols = 7;
        config.boundary = periodic ? "periodic" : "fixed";
        config.minimizer = minimizer;
        config.preconditioner = modes[mode];
        config.initialGuessNoise = 1e-3;
        config.epsR = 1e-6;
        config.showProgress = -1;
        config.forceReRun = true;
        Simulation simulation(config, "", false);
        if (!periodic) {
          simulation.mesh.fixBorderNodes();
        }
        simulation.mesh.applyTransformation(getShear(0.08));
        simulation.initSolver();
        simulation.applyLoadStepToGuess();
        simulation.addNoiseToGuess();
        simulation.minimize(false);
        energies[mode] = simulation.mesh.totalEnergy;
        CHECK(simulation.mesh.maxForce < 2e-6);
        CHECK(simulation.preconditionerStatistics().fallbacks == 0);
        CHECK(simulation.preconditionerStatistics().builds == size_t(mode != 0));
        if (mode != 0) {
          simulation.applyLoadStepToGuess();
          simulation.minimize(false);
          CHECK(simulation.preconditionerStatistics().reuses == 1);
          simulation.config.preconditioner = "none";
          simulation.applyLoadStepToGuess();
          simulation.minimize(false);
          CHECK(simulation.mesh.maxForce < 2e-6);
        }
      }
      CHECK(energies[0] == doctest::Approx(energies[1]).epsilon(1e-8));
      CHECK(energies[0] == doctest::Approx(energies[2]).epsilon(1e-8));
    }
  }
}

TEST_CASE("Changing preconditioner modes clears native solver state") {
  for (bool periodic : {false, true}) {
    for (const std::string minimizer : {"LBFGS", "CG"}) {
      Config config;
      config.setDefaultValues();
      config.rows = config.cols = 5;
      config.boundary = periodic ? "periodic" : "fixed";
      config.minimizer = minimizer;
      config.showProgress = -1;
      config.forceReRun = true;
      Simulation simulation(config, "", false);
      if (!periodic) {
        simulation.mesh.fixBorderNodes();
      }
      simulation.mesh.applyTransformation(getShear(0.08));
      simulation.initSolver();
      for (const std::string mode : {"diagonal", "exactTangent", "diagonal",
                                     "none", "diagonal", "none"}) {
        CAPTURE(mode);
        CAPTURE(minimizer);
        CAPTURE(periodic);
        simulation.config.preconditioner = mode;
        simulation.applyLoadStepToGuess();
        simulation.addNoiseToGuess(1e-3);
        simulation.minimize(false);
        const int expected = mode == "diagonal" ? 2 : 0;
        if (minimizer == "LBFGS") {
          CHECK(simulation.dataLink.LBFGS_state->c_ptr()->prectype == expected);
        } else {
          CHECK(simulation.dataLink.CG_state->c_ptr()->prectype == expected);
        }
        CHECK(simulation.preconditionerStatistics().fallbacks == 0);
        CHECK(simulation.mesh.maxForce < config.epsR);
      }
    }
  }
}
