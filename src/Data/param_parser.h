#ifndef param_parser_H
#define param_parser_H
#pragma once

#include "cereal_help.h"
#include <cereal/cereal.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/string.hpp>
#include <iostream>
#include <string>

#include <Param.h>

// I am sorry about writing a custom parser.
// I tried to add a Yaml parser, but they were problamatic

struct Config {
  // Simulation settings
  std::string name;
  int rows;
  int cols;
  bool usingPBC;
  std::string scenario;
  int nrThreads;
  int seed;
  double QDSD;
  double initialGuessNoise;
  std::string meshDiagonal;

  // Loading settings
  double startLoad;
  double loadIncrement;
  double maxLoad;

  // Minimizer settings
  std::string minimizer; // FIRE / LBFGS

  // Max residual force (using the same value across all algorithms)
  double epsR;

  // LBFGS settings
  int LBFGSNrCorrections;
  double LBFGSScale;
  double LBFGSEpsg;
  double LBFGSEpsf;
  double LBFGSEpsx;
  int LBFGSMaxIterations;
  // Conjugate Gradient settings
  double CGScale;
  double CGEpsg;
  double CGEpsf;
  double CGEpsx;
  int CGMaxIterations;
  // FIRE settings
  double finc;
  double fdec;
  double alphaStart;
  double falpha;
  double dtStart;
  double dtMax; // 10 times dtStart
  double dtMin; // 0.002 times dtStart
  double maxCompS;
  double eps;
  double epsRel;
  double delta;
  int maxIt;

  // Logging settings
  bool logDuringMinimization;
  // If a certain percentage of elements go through a m3 transformation, we log
  double plasticityEventThreshold;
  // If an energy drop is above this is threshold, we log
  double energyDropThreshold;
  int showProgress;

  // Other
  std::string configPath;

  // The program will check if the folder already contains a completed
  // simulation If the simulation is complete, the program will terminate and
  // not rerun the simulation unless forceReRun is true
  // Currently not really used
  bool forceReRun;

  void setDefaultValues();

  friend std::ostream &operator<<(std::ostream &os, const Config &config);
  std::string str() const;

  void updateParam(FIREpp::FIREParam<double> &param);

  template <class Archive> void serialize(Archive &ar) {
    // General simulation settings
    ar(MAKE_NVP(name), MAKE_NVP(rows), MAKE_NVP(cols), MAKE_NVP(usingPBC),
       MAKE_NVP(scenario), MAKE_NVP(nrThreads), MAKE_NVP(seed), MAKE_NVP(QDSD),
       MAKE_NVP(initialGuessNoise));

    LOAD_WITH_DEFAULT(ar, meshDiagonal, std::string("major"));

    // Load settings
    ar(MAKE_NVP(startLoad), MAKE_NVP(loadIncrement), MAKE_NVP(maxLoad));

    // Max force allowed
    LOAD_WITH_DEFAULT(ar, epsR, 0.0);

    // LBFGS minimizer settings
    ar(MAKE_NVP(minimizer), MAKE_NVP(LBFGSNrCorrections), MAKE_NVP(LBFGSScale),
       MAKE_NVP(LBFGSEpsg), MAKE_NVP(LBFGSEpsf), MAKE_NVP(LBFGSEpsx),
       MAKE_NVP(LBFGSMaxIterations));

    // CG minimizer settings
    ar(MAKE_NVP(CGScale), MAKE_NVP(CGEpsg), MAKE_NVP(CGEpsf), MAKE_NVP(CGEpsx),
       MAKE_NVP(CGMaxIterations));

    // Simulation step settings
    ar(MAKE_NVP(finc), MAKE_NVP(fdec), MAKE_NVP(alphaStart), MAKE_NVP(falpha),
       MAKE_NVP(dtStart), MAKE_NVP(dtMax), MAKE_NVP(dtMin), MAKE_NVP(maxCompS),
       MAKE_NVP(eps), MAKE_NVP(epsRel), MAKE_NVP(delta), MAKE_NVP(maxIt));

    // Stopping conditions and progress display
    ar(MAKE_NVP(plasticityEventThreshold), MAKE_NVP(energyDropThreshold),
       MAKE_NVP(showProgress));

    LOAD_WITH_DEFAULT(ar, logDuringMinimization, false);

    // File paths and execution options
    ar(MAKE_NVP(configPath), MAKE_NVP(forceReRun));
  }
};

std::map<std::string, std::string> parseParams(const std::string &filename);

// Function to initialize Config from a map
Config initializeConfig(const std::map<std::string, std::string> &configMap);

Config parseConfigFile(std::string configFile);

#endif
