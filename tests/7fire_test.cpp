#include "Data/cereal_help.h"
#include "FIRE.h"
#include "Param.h"
#include "Simulation/simulation.h"
#include "run/doctest.h"
#include <Eigen/Core>
#include <cmath>

TEST_CASE("FIRE checks EpsX on the full step which triggers a reset") {
  FIREpp::FIREParam<double> parameters(1, 1);
  parameters.epsilon = 0;
  parameters.epsilon_x = 1.2;
  parameters.epsilon_rel = 0;
  parameters.delta = 0;
  parameters.dt_start = 1.5;
  parameters.dt_max = 1.5;
  parameters.max_iterations = 2;

  Eigen::VectorXd x(1);
  x[0] = 1.0;
  double energy = 0;
  int terminationType = 0;
  auto quadratic = [](Eigen::VectorXd &position, Eigen::VectorXd &gradient,
                      void *) {
    gradient = position;
    return 0.5 * position.squaredNorm();
  };

  FIREpp::FIRESolver<double> solver(parameters);
  const int iterations =
      solver.minimize(quadratic, x, energy, nullptr, terminationType);

  CHECK(iterations == 1);
  CHECK(terminationType == 2);
  CHECK(x[0] == doctest::Approx(0.26875));
  CHECK(energy == doctest::Approx(0.5 * x.squaredNorm()));
}

TEST_CASE("FIRE does not treat an uphill reset as EpsX convergence") {
  FIREpp::FIREParam<double> parameters(1, 1);
  parameters.epsilon = 0;
  parameters.epsilon_x = 0.7;
  parameters.epsilon_rel = 0;
  parameters.delta = 0;
  parameters.dt_start = 1.5;
  parameters.dt_max = 1.5;
  parameters.max_iterations = 2;

  Eigen::VectorXd x(1);
  x[0] = 1.0;
  double energy = 0;
  int terminationType = 0;
  auto quadratic = [](Eigen::VectorXd &position, Eigen::VectorXd &gradient,
                      void *) {
    gradient = position;
    return 0.5 * position.squaredNorm();
  };

  FIREpp::FIRESolver<double> solver(parameters);
  const int iterations =
      solver.minimize(quadratic, x, energy, nullptr, terminationType);

  CHECK(iterations == 2);
  CHECK(terminationType == 4);
}

TEST_CASE("FIRE does not check EpsX on ordinary downhill steps") {
  FIREpp::FIREParam<double> parameters(1, 1);
  parameters.epsilon = 0;
  parameters.epsilon_x = 1e-2;
  parameters.epsilon_rel = 0;
  parameters.delta = 0;
  parameters.max_iterations = 2;

  Eigen::VectorXd x(1);
  x[0] = 1.0;
  double energy = 0;
  int terminationType = 0;
  auto quadratic = [](Eigen::VectorXd &position, Eigen::VectorXd &gradient,
                      void *) {
    gradient = position;
    return 0.5 * position.squaredNorm();
  };

  FIREpp::FIRESolver<double> solver(parameters);
  const int iterations =
      solver.minimize(quadratic, x, energy, nullptr, terminationType);

  CHECK(iterations == 2);
  CHECK(terminationType == 4);
}

TEST_CASE("Legacy dumps default FIRE EpsX to disabled") {
  Simulation simulation;
  REQUIRE_NOTHROW(loadFromFile("tests/oldDumps/dump_l0.20.xml.gz", simulation));
  CHECK(simulation.config.FIREEpsx == 0.0);
}
