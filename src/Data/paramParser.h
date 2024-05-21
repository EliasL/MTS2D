#ifndef PARAMPARSER_H
#define PARAMPARSER_H
#pragma once

#include <cereal/types/string.hpp>
#include <iostream>
#include <map>
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
  double plasticityEventThreshold;

  // Loading settings
  double startLoad;
  double loadIncrement;
  double maxLoad;
  double noise;

  // Minimizer settings
  std::string minimizer; // FIRE / LBFGS
  // LBFGS settings
  int LBFGSNrCorrections;
  double LBFGSScale;
  double LBFGSEpsg;
  double LBFGSEpsf;
  double LBFGSEpsx;
  int LBFGSMaxIterations;
  // FIRE settings
  double finc;
  double fdec;
  double alphaStart;
  double falpha;
  double dtStart;
  double dtStartMax;
  double dtMax; // 10 times dtStart
  double dtMin; // 0.002 times dtStart
  double eps;
  double epsRel;
  double delta;
  int maxIt;

  // Logging settings
  int showProgress;

  // Other
  std::string configPath;

  friend std::ostream &operator<<(std::ostream &os, const Config &config);
  std::string str() const;

  friend class cereal::access;
  template <class Archive> void serialize(Archive &ar) {
    ar(name, rows, cols, usingPBC, scenario, nrThreads, seed,
       plasticityEventThreshold, startLoad, loadIncrement, maxLoad, noise,
       minimizer, LBFGSNrCorrections, LBFGSScale, LBFGSEpsg, LBFGSEpsf,
       LBFGSEpsx, LBFGSMaxIterations, finc, fdec, alphaStart, falpha, dtStart,
       dtStartMax, dtMax, dtMin, eps, epsRel, delta, maxIt, showProgress,
       configPath);
  }

  void updateParam(FIREpp::FIREParam<double> &param);
};

std::map<std::string, std::string> parseParams(const std::string &filename);

// Function to initialize Config from a map
Config initializeConfig(const std::map<std::string, std::string> &configMap);

Config parseConfigFile(std::string configFile);

#endif
