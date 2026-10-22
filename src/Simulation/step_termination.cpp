#include "step_termination.h"
#include "finite.h"
#include "tangent_stiffness.h"
#include <algorithm>
#include <array>
#include <stdexcept>

void StepTermination::initialize(const Mesh &mesh, const Config &config) {
  if (referenceVolume == 0) {
    for (const auto &element : mesh.elements) {
      const double area = element.getInitArea();
      if (!finiteValue(area) || area <= 0) {
        throw std::runtime_error("Step counting requires positive reference areas.");
      }
      referenceVolume += area;
    }
  }
  if (!finiteValue(referenceVolume) || referenceVolume <= 0) {
    throw std::runtime_error("Invalid step-counting reference volume.");
  }
  if (config.stopAfterPlasticSteps || config.stopAfterElasticSteps) {
    shearModulus = referenceShearModulus(-0.25, config.bulkModulus);
    if (!finiteValue(shearModulus) || shearModulus <= 0) {
      throw std::runtime_error("Plastic classification requires a positive reference modulus.");
    }
  }
  targetsMet = reached(config);
  stopReason = targetsMet ? "sample_targets" : "running";
}

bool StepTermination::reached(const Config &config) const {
  const std::array<std::pair<int, std::uint64_t>, 5> targets{{
      {config.stopAfterSteps, steps}, {config.stopAfterPlasticSteps, plasticSteps},
      {config.stopAfterElasticSteps, elasticSteps},
      {config.stopAfterEnergyDrops, energyDrops},
      {config.stopAfterEnergyIncreases, energyIncreases}}};
  const bool all = config.stopWhen == "all";
  bool enabled = false, satisfied = all;
  for (const auto &[target, count] : targets) {
    if (target > 0) {
      enabled = true;
      satisfied = all ? satisfied && count >= static_cast<std::uint64_t>(target)
                      : satisfied || count >= static_cast<std::uint64_t>(target);
    }
  }
  return enabled && satisfied;
}

void StepTermination::record(const Config &config, double beforeEnergy,
                             double finalEnergy, double guessEnergy,
                             double increment) {
  if (!finiteValue(beforeEnergy) || !finiteValue(finalEnergy) ||
      !finiteValue(guessEnergy) || !finiteValue(increment) || increment <= 0) {
    throw std::runtime_error("Nonfinite energy or invalid completed step.");
  }
  const double change = (finalEnergy - beforeEnergy) / referenceVolume;
  const double release = (guessEnergy - finalEnergy) / referenceVolume;
  if (!finiteValue(change) || !finiteValue(release)) {
    throw std::runtime_error("Nonfinite step energy density.");
  }
  const double tolerance = config.stepEnergyTolerance;
  const bool classify = config.stopAfterPlasticSteps || config.stopAfterElasticSteps;
  if (classify) {
    if (release < -tolerance) {
      throw std::runtime_error("Relaxation increased energy beyond stepEnergyTolerance.");
    }
    // rho cancels from kappa >= mu/(2*rho). Use the actual increment.
    const double threshold = 0.5 * shearModulus * increment * increment;
    if (!finiteValue(threshold) || threshold <= 0) {
      throw std::runtime_error("Invalid plastic-step detection threshold.");
    }
    stepClass = std::max(0.0, release) >= threshold ? "plastic" : "elastic";
    plasticSteps += stepClass == "plastic";
    elasticSteps += stepClass == "elastic";
  } else {
    stepClass = "unclassified";
  }
  ++steps;
  energyDrops += change < -tolerance;
  energyIncreases += change > tolerance;
  targetsMet = reached(config);
  stopReason = targetsMet ? "sample_targets" : "running";
}
