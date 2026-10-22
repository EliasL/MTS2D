#pragma once
#include "Data/param_parser.h"
#include "Mesh/mesh.h"
#include <cstdint>

struct StepTermination {
  std::uint64_t steps = 0, plasticSteps = 0, elasticSteps = 0;
  std::uint64_t energyDrops = 0, energyIncreases = 0;
  double referenceVolume = 0, shearModulus = 0;
  std::string stepClass = "initialization", stopReason = "running";
  bool targetsMet = false;

  void initialize(const Mesh &mesh, const Config &config);
  bool reached(const Config &config) const;
  void record(const Config &config, double beforeEnergy, double finalEnergy,
              double guessEnergy, double increment);

  template <class Archive> void serialize(Archive &ar) {
    ar(MAKE_NVP(steps), MAKE_NVP(plasticSteps), MAKE_NVP(elasticSteps),
       MAKE_NVP(energyDrops), MAKE_NVP(energyIncreases),
       MAKE_NVP(referenceVolume), MAKE_NVP(shearModulus), MAKE_NVP(stepClass),
       MAKE_NVP(stopReason), MAKE_NVP(targetsMet));
  }
};
