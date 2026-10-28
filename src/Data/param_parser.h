#ifndef param_parser_H
#define param_parser_H
#pragma once

#include "cereal_help.h"
#include "Simulation/sparse_cholesky.h"
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
  std::string loading = "relaxation";
  std::string deformation = "simpleShear";
  std::string boundary = "fixed";
  std::string directions = "horizontal";
  double segmentLoad = 0;
  bool repeat = false;
  double stepNoise = 0;
  double reversibilityTolerance = 1e-4;
  int checkpointInterval = 0;
  int stopAfterSteps = 0;
  int stopAfterPlasticSteps = 0;
  int stopAfterElasticSteps = 0;
  int stopAfterEnergyDrops = 0;
  int stopAfterEnergyIncreases = 0;
  std::string stopWhen = "all";
  double stepEnergyTolerance = 1e-12; // Energy per reference area.
  int nrThreads;
  int seed;
  double QDSD;
  double initialGuessNoise;
  // Affinely predict free-node positions after a loading step.
  bool affineInitialGuess = true;
  std::string meshDiagonal;
  std::string reconnectionMethod; // "none", "edgeFlip", "delaunay"
  bool reconnectRevert;
  bool reconnectEdgeLocking;
  bool minimum_transport_distance = true; // false selects maximum flip-distance
  std::string energyFunction;     // "contiSquare", "contiTriangular"
  double bulkModulus;

  // Loading settings
  double startLoad;
  double loadIncrement;
  double totalLoad;

  // Minimizer settings
  std::string minimizer; // FIRE / LBFGS / CG
  std::string preconditioner; // exactTangent / diagonal / none (LBFGS and CG)
  std::string tangentBackend = SparseCholesky::name(SparseCholesky::defaultBackend());
  std::string preconditionerUpdate; // adaptive / every / frozen
  double preconditionerDisplacementThreshold;

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
  bool fullMinimizationLogging;
  bool writeDebugVTUs;
  bool writeDumps;
  bool writeMeshVTUs;
  // Targeted replay option: save and test forward steps without m3 changes.
  // This is disabled for ordinary simulations.
  bool saveElasticReversibilityStates;
  int maximumSavedElasticReversibilityStates;
  int nrVTUFrames;
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
  bool forceReRun = false;

  void setDefaultValues();
  void validatePreconditioner() const;
  SparseBackend selectedTangentBackend() const;
  bool usesPeriodicBoundary() const { return boundary == "periodic" || boundary == "periodicAnchored"; }

  friend std::ostream &operator<<(std::ostream &os, const Config &config);
  std::string str() const;

  void updateParam(FIREpp::FIREParam<double> &param);

  template <class Archive> void serialize(Archive &ar) {
    // General simulation settings
    ar(MAKE_NVP(name), MAKE_NVP(rows), MAKE_NVP(cols),
       MAKE_NVP(reconnectionMethod));
    ar(MAKE_NVP(boundary), MAKE_NVP(loading), MAKE_NVP(deformation),
       MAKE_NVP(directions), MAKE_NVP(segmentLoad), MAKE_NVP(repeat),
       MAKE_NVP(stepNoise), MAKE_NVP(reversibilityTolerance),
       MAKE_NVP(checkpointInterval));
    ar(MAKE_NVP(stopAfterSteps), MAKE_NVP(stopAfterPlasticSteps),
       MAKE_NVP(stopAfterElasticSteps), MAKE_NVP(stopAfterEnergyDrops),
       MAKE_NVP(stopAfterEnergyIncreases), MAKE_NVP(stopWhen),
       MAKE_NVP(stepEnergyTolerance));
    ar(MAKE_NVP(nrThreads), MAKE_NVP(seed), MAKE_NVP(QDSD),
       MAKE_NVP(initialGuessNoise));

    ar(MAKE_NVP(affineInitialGuess), MAKE_NVP(meshDiagonal),
       MAKE_NVP(reconnectRevert), MAKE_NVP(reconnectEdgeLocking),
       MAKE_NVP(minimum_transport_distance), MAKE_NVP(energyFunction),
       MAKE_NVP(bulkModulus));

    // Load settings
    ar(MAKE_NVP(startLoad), MAKE_NVP(loadIncrement), MAKE_NVP(totalLoad));

    // Max force allowed
    ar(MAKE_NVP(epsR));

    // LBFGS minimizer settings
    ar(MAKE_NVP(minimizer), MAKE_NVP(LBFGSNrCorrections), MAKE_NVP(LBFGSScale),
       MAKE_NVP(LBFGSEpsg), MAKE_NVP(LBFGSEpsf), MAKE_NVP(LBFGSEpsx),
       MAKE_NVP(LBFGSMaxIterations));
    ar(MAKE_NVP(preconditioner), MAKE_NVP(preconditionerUpdate),
       MAKE_NVP(preconditionerDisplacementThreshold));

    // CG minimizer settings
    ar(MAKE_NVP(CGScale), MAKE_NVP(CGEpsg), MAKE_NVP(CGEpsf), MAKE_NVP(CGEpsx),
       MAKE_NVP(CGMaxIterations));

    // Simulation step settings
    ar(MAKE_NVP(finc), MAKE_NVP(fdec), MAKE_NVP(alphaStart), MAKE_NVP(falpha),
       MAKE_NVP(dtStart), MAKE_NVP(dtMax), MAKE_NVP(dtMin), MAKE_NVP(maxCompS),
       MAKE_NVP(eps));
    ar(MAKE_NVP(epsRel), MAKE_NVP(delta), MAKE_NVP(maxIt));

    // Stopping conditions and progress display
    ar(MAKE_NVP(plasticityEventThreshold), MAKE_NVP(energyDropThreshold),
       MAKE_NVP(showProgress));

    ar(MAKE_NVP(logDuringMinimization), MAKE_NVP(fullMinimizationLogging),
       MAKE_NVP(writeDebugVTUs), MAKE_NVP(writeDumps), MAKE_NVP(writeMeshVTUs),
       MAKE_NVP(saveElasticReversibilityStates),
       MAKE_NVP(maximumSavedElasticReversibilityStates),
       MAKE_NVP(nrVTUFrames));

    // File paths and execution options
    ar(MAKE_NVP(configPath), MAKE_NVP(forceReRun));
    // This trailing field is absent only in legacy checkpoints.
    if constexpr (Archive::is_loading::value) {
      tangentBackend = "eigen";
      if (ar.getNodeName() != nullptr) ar(MAKE_NVP(tangentBackend));
      selectedTangentBackend();
    } else {
      ar(MAKE_NVP(tangentBackend));
    }
  }
};

std::map<std::string, std::string> parseParams(const std::string &filename);

// Function to initialize Config from a map
Config initializeConfig(const std::map<std::string, std::string> &configMap);

Config parseConfigFile(std::string configFile,
                       const std::string &tangentBackendDefault = "");

// Global quiet flag for suppressing console output in tests or batch runs.
inline bool g_quiet = false;
inline void setQuiet(bool quiet) { g_quiet = quiet; }
inline bool isQuiet() { return g_quiet; }

#endif
