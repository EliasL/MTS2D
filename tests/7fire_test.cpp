#include "Data/cereal_help.h"
#include "FIRE.h"
#include "Param.h"
#include "Simulation/simulation.h"
#include "run/doctest.h"
#include <Eigen/Core>
#include <cmath>

TEST_CASE("Removed FIRE step-size tolerance is rejected") {
  CHECK_THROWS_AS(initializeConfig({{"FIREEpsx", "1e-6"}}), std::invalid_argument);
}

TEST_CASE("FIRE velocity resets do not imply force convergence") {
  FIREpp::FIREParam<double> parameters(1, 1);
  parameters.epsilon_R = 1e-5;
  parameters.epsilon = parameters.epsilon_rel = parameters.delta = 0;
  parameters.dt_start = parameters.dt_max = 1.5;
  parameters.max_iterations = 2;
  Eigen::VectorXd x = Eigen::VectorXd::Ones(1);
  double energy;
  int terminationType = 0;
  auto quadratic = [](Eigen::VectorXd &position, Eigen::VectorXd &gradient, void *) {
    gradient = position;
    return 0.5 * position.squaredNorm();
  };
  FIREpp::FIRESolver<double> solver(parameters);
  CHECK(solver.minimize(quadratic, x, energy, nullptr, terminationType) == 2);
  CHECK(terminationType == 4);
  CHECK(x.cwiseAbs().maxCoeff() > parameters.epsilon_R);
}

TEST_CASE("FIRE residual tolerance bounds components rather than the L2 norm") {
  FIREpp::FIREParam<double> parameters(1, 2);
  parameters.epsilon_R = 1e-3;
  parameters.epsilon = parameters.epsilon_rel = parameters.delta = 0;
  Eigen::VectorXd x = Eigen::VectorXd::Constant(2, 0.75e-3);
  auto quadratic = [](Eigen::VectorXd &position, Eigen::VectorXd &gradient, void *) {
    gradient = position;
    return 0.5 * position.squaredNorm();
  };
  double energy;
  int terminationType = 0;
  FIREpp::FIRESolver<double> solver(parameters);
  solver.minimize(quadratic, x, energy, nullptr, terminationType);
  CHECK(x.norm() > parameters.epsilon_R);
  CHECK(terminationType == 3);

  parameters.epsilon_R = 0.99999;
  x.setOnes();
  solver.minimize(quadratic, x, energy, nullptr, terminationType);
  CHECK(terminationType == 1);
  CHECK(x.cwiseAbs().maxCoeff() < parameters.epsilon_R);
  parameters.epsilon_R = -1;
  CHECK_THROWS_AS(parameters.check_param(), std::invalid_argument);
}

TEST_CASE("FIRE refreshes the configured residual tolerance before each solve") {
  Config config;
  config.setDefaultValues();
  config.rows = config.cols = 4;
  config.boundary = "periodic";
  config.minimizer = "FIRE";
  config.nrThreads = 1;
  config.showProgress = -1;
  config.epsR = 0;
  Simulation simulation(config, "", false);
  simulation.initializeLoadedStateWithoutOutput();
  simulation.mesh.applyTransformation(getShear(0.08));
  simulation.applyLoadStepToGuess();
  simulation.addNoiseToGuess(1e-3);
  simulation.config.epsR = 1e-4;
  simulation.minimize(false);
  CHECK(simulation.mesh.maxForce < simulation.config.epsR);
  CHECK(simulation.FIRERep.termType == 1);
}
