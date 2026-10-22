#include "Simulation/step_termination.h"
#include "Simulation/tangent_stiffness.h"
#include "Simulation/finite.h"
#include "run/doctest.h"
#include <cfloat>
#include <cstdlib>

TEST_CASE("Step targets ignore zero values and combine by all or any") {
  Config config;
  config.setDefaultValues();
  StepTermination state;
  for (const std::string mode : {"all", "any"}) {
    config.stopWhen = mode;
    CHECK_FALSE(state.reached(config));
  }
  config.stopAfterSteps = 3;
  config.stopAfterPlasticSteps = 2;
  state.steps = 3;
  CHECK(state.reached(config));
  config.stopWhen = "all";
  CHECK_FALSE(state.reached(config));
  state.plasticSteps = 2;
  CHECK(state.reached(config));
  const auto xml = serializeToXml(state);
  StepTermination restored;
  deserializeFromXml(xml, restored);
  CHECK(restored.steps == 3);
  CHECK(restored.plasticSteps == 2);
}

TEST_CASE("Step classification uses relaxation energy, actual increment and reference area") {
  Config config;
  config.setDefaultValues();
  config.stopAfterPlasticSteps = 2;
  StepTermination state;
  state.referenceVolume = 2;
  state.shearModulus = 4;
  // Equality is plastic: release/area = mu * increment^2 / 2 = 0.03125.
  state.record(config, 10, 9, 9.0625, 0.125);
  CHECK(state.plasticSteps == 1);
  CHECK(state.energyDrops == 1);
  state.record(config, 9, 10, 10, 0.125);
  CHECK(state.elasticSteps == 1);
  CHECK(state.energyIncreases == 1);
  state.record(config, 10, 10, 10.03125, 0.125);
  CHECK(state.elasticSteps == 2);
  // The same release is plastic for the shorter actual step.
  state.record(config, 10, 10, 10.03125, 0.0625);
  CHECK(state.plasticSteps == 2);
  CHECK(state.targetsMet);
  CHECK(state.stopReason == "sample_targets");
  CHECK(state.steps == 4);
  CHECK(referenceShearModulus(-0.25, 4) == doctest::Approx(401.0 / 66.0));

  StepTermination larger;
  larger.referenceVolume = 8;
  larger.shearModulus = 4;
  larger.record(config, 40, 36, 36.25, 0.125);
  CHECK(larger.plasticSteps == 1);
  CHECK(larger.energyDrops == 1);
}

TEST_CASE("Energy tolerance defines neutral changes and rejects invalid relaxation") {
  Config config;
  config.setDefaultValues();
  config.stepEnergyTolerance = 0.125;
  StepTermination state;
  state.referenceVolume = 2;
  state.record(config, 1, 1.25, 1.25, 0.1);
  state.record(config, 1.25, 1, 1, 0.1);
  CHECK(state.energyDrops == 0);
  CHECK(state.energyIncreases == 0);
  config.stopAfterElasticSteps = 5;
  state.shearModulus = 4;
  state.record(config, 1, 1, 0.875, 0.1);
  CHECK(state.elasticSteps == 1);
  CHECK_THROWS_AS(state.record(config, 1, 1, 0.5, 0.1), std::runtime_error);
  CHECK(state.steps == 3);
  CHECK(finiteValue(DBL_MAX));
  for (const char *invalid : {"nan", "inf", "-inf"}) {
    const double value = std::strtod(invalid, nullptr);
    CHECK_FALSE(finiteValue(value));
    CHECK_THROWS_AS(state.record(config, value, 1, 1, 0.1), std::runtime_error);
    CHECK_THROWS_AS(state.record(config, 1, value, 1, 0.1), std::runtime_error);
    CHECK_THROWS_AS(state.record(config, 1, 1, value, 0.1), std::runtime_error);
  }
}

TEST_CASE("Step termination configuration rejects undefined choices") {
  CHECK_THROWS_AS(initializeConfig({{"stopOnEnergyDrop", "true"}}), std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"stopWhen", "sometimes"}}), std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"stopAfterSteps", "-1"}}), std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"stopAfterSteps", "1.5"}}), std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"stepEnergyTolerance", "nan"}}), std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"stopAfterPlasticSteps", "1"},
                                   {"loading", "noRelaxation"}}), std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"stopAfterElasticSteps", "1"},
                                   {"energyFunction", "contiTriangular"}}), std::invalid_argument);
}
