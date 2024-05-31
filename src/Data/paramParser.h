#ifndef PARAMPARSER_H
#define PARAMPARSER_H
#pragma once

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
  double quenchedDisorderStrength;
  double initialGuessNoise;

  // Loading settings
  double startLoad;
  double loadIncrement;
  double maxLoad;

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
  double dtMax; // 10 times dtStart
  double dtMin; // 0.002 times dtStart
  double maxCompS;
  double eps;
  double epsRel;
  double delta;
  int maxIt;

  // Logging settings
  double plasticityEventThreshold;
  int showProgress;

  // Other
  std::string configPath;

  // The program will check if the folder already contains a completed
  // simulation If the simulation is complete, the program will terminate and
  // not rerun the simulation unless forceReRun is true
  bool forceReRun;

  friend std::ostream &operator<<(std::ostream &os, const Config &config);
  std::string str() const;

  void updateParam(FIREpp::FIREParam<double> &param);
};

std::map<std::string, std::string> parseParams(const std::string &filename);

// Function to initialize Config from a map
Config initializeConfig(const std::map<std::string, std::string> &configMap);

Config parseConfigFile(std::string configFile);

#endif
