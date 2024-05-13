#ifndef PARAMPARSER_H
#define PARAMPARSER_H
#pragma once

#include <cereal/types/string.hpp>
#include <iostream>
#include <map>
#include <string>
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
  std::string minimizer;
  int nrCorrections;
  double scale;
  double epsg;
  double epsf;
  double epsx;
  int maxIterations;
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
       minimizer, nrCorrections, scale, epsg, epsf, epsx, maxIterations,
       showProgress);
  }
};

std::map<std::string, std::string> parseParams(const std::string &filename);

// Function to initialize Config from a map
Config initializeConfig(const std::map<std::string, std::string> &configMap);

Config parseConfigFile(std::string configFile);

#endif
