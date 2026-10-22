#include "param_parser.h"
#include "Simulation/protocols.h"
#include <cfloat>
#include <set>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
// Define a macro to get the variable name as a string
#define GET_VALUE(var) getValue(configMap, #var, var, recognized)

void Config::setDefaultValues() {
  // General settings
  name = "defaultName";
  rows = 3;
  cols = 3;
  reconnectionMethod = "none"; // "none", "edgeFlip", "delaunay"
  loading = "relaxation";
  deformation = "simpleShear";
  boundary = "fixed";
  directions = "horizontal";
  segmentLoad = 0;
  repeat = false;
  stepNoise = 0;
  reversibilityTolerance = 1e-4;
  checkpointInterval = 0;
  stopAfterSteps = 0;
  stopAfterPlasticSteps = 0;
  stopAfterElasticSteps = 0;
  stopAfterEnergyDrops = 0;
  stopAfterEnergyIncreases = 0;
  stopWhen = "all";
  stepEnergyTolerance = 1e-12;
  nrThreads = 1;
  seed = 0;
  QDSD = 0.0;
  initialGuessNoise = 0.05;
  affineInitialGuess = true;
  meshDiagonal = "major";
  energyFunction = "contiSquare"; // "contiSquare", "contiTriangular"
  bulkModulus = 4;                //  BulkModulus
  reconnectRevert = true;
  reconnectEdgeLocking = false;
  minimum_transport_distance = true;

  // Loading settings
  startLoad = 0.0;
  loadIncrement = 0.1;
  totalLoad = 1.0;

  // Minimizer settings
  minimizer = "LBFGS";
  preconditioner = "exactTangent";
  preconditionerUpdate = "adaptive";
  preconditionerDisplacementThreshold = 0.05;

  // Max residual force allowed. Minimize until all forces are smaller than this
  epsR = 1e-5;

  // LBFGS-specific settings
  LBFGSNrCorrections = 10;
  LBFGSScale = 1.0;
  LBFGSEpsg = 1e-20;
  LBFGSEpsf = 1e-20;
  LBFGSEpsx = 1e-20;
  LBFGSMaxIterations = 0;

  // CG-specific settings
  CGScale = 1.0;
  CGEpsg = 1e-20;
  CGEpsf = 1e-20;
  CGEpsx = 1e-20;
  CGMaxIterations = 0;

  // FIRE-specific settings
  finc = 1.1;
  fdec = 0.5;
  alphaStart = 0.1;
  falpha = 0.99;
  dtStart = 0.1;
  dtMax = dtStart * 3;
  dtMin = dtStart * 0.000001;
  maxCompS = 0.01;
  eps = 1e-20;
  epsRel = 0;
  delta = 0;
  maxIt = 100000;

  // Additional settings
  logDuringMinimization = false;
  fullMinimizationLogging = false;
  writeDebugVTUs = false;
  writeDumps = true;
  writeMeshVTUs = true;
  saveElasticReversibilityStates = false;
  maximumSavedElasticReversibilityStates = 0;
  nrVTUFrames = 200;
  plasticityEventThreshold = 0.01;
  energyDropThreshold = 0.001;
  showProgress = 1;

  // Config path default
  configPath = "";
  forceReRun = false;
}

void Config::validatePreconditioner() const {
  if (preconditioner != "none" && preconditioner != "exactTangent" &&
      preconditioner != "diagonal") {
    throw std::invalid_argument(
        "preconditioner must be exactTangent, diagonal, or none.");
  }
  if (preconditionerUpdate != "adaptive" && preconditionerUpdate != "every" &&
      preconditionerUpdate != "frozen") {
    throw std::invalid_argument(
        "preconditionerUpdate must be adaptive, every, or frozen.");
  }
  if (!(preconditionerDisplacementThreshold > 0 &&
        preconditionerDisplacementThreshold <= DBL_MAX)) {
    throw std::invalid_argument(
        "preconditionerDisplacementThreshold must be positive.");
  }
}

void Config::updateParam(FIREpp::FIREParam<double> &param) {
  param.epsilon_R = epsR;
  param.finc = finc;
  param.fdec = fdec;
  param.alpha_start = alphaStart;
  param.falpha = falpha;
  param.dt_start = dtStart;
  param.dt_max = dtMax;
  param.dt_min = dtMin;
  param.max_component_step = maxCompS;
  param.epsilon = eps;
  param.epsilon_rel = epsRel;
  param.delta = delta;
  param.max_iterations = maxIt;
}

std::string Config::str() const {
  std::stringstream ss;
  ss << *this;
  return ss.str();
}
std::ostream &operator<<(std::ostream &os, const Config &config) {
  os << "Name: " << config.name << "\n"
     << "Rows, Cols: " << config.rows << ", " << config.cols << "\n"
     << "Boundary Conditions: " << (config.usesPeriodicBoundary() ? "PBC" : "NPBC") << "\n"
     << "Reconnection method: " << config.reconnectionMethod << "\n"
     << "Loading: " << config.loading << "\n"
     << "Deformation: " << config.deformation << "\n"
     << "Directions: " << config.directions << "\n"
     << "Segment load: " << config.segmentLoad << "\n"
     << "Repeat: " << config.repeat << "\n"
     << "Step noise: " << config.stepNoise << "\n"
     << "Reversibility tolerance: " << config.reversibilityTolerance << "\n"
     << "Checkpoint interval: " << config.checkpointInterval << "\n"
     << "stopAfterSteps: " << config.stopAfterSteps << "\n"
     << "stopAfterPlasticSteps: " << config.stopAfterPlasticSteps << "\n"
     << "stopAfterElasticSteps: " << config.stopAfterElasticSteps << "\n"
     << "stopAfterEnergyDrops: " << config.stopAfterEnergyDrops << "\n"
     << "stopAfterEnergyIncreases: " << config.stopAfterEnergyIncreases << "\n"
     << "stopWhen: " << config.stopWhen << "\n"
     << "stepEnergyTolerance: " << config.stepEnergyTolerance << "\n"
     << "Number of Threads: " << config.nrThreads << "\n"
     << "Seed: " << config.seed << "\n"
     << "Quenched disorder standard deviation: " << config.QDSD << "\n"
     << "Initial guess noise: " << config.initialGuessNoise << "\n"
     << "Affine initial guess: " << config.affineInitialGuess << "\n"
     << "Mesh diagonal: " << config.meshDiagonal << "\n"
     << "Reconnect revert: " << config.reconnectRevert << "\n"
     << "Reconnect edge locking: " << config.reconnectEdgeLocking << "\n"
     << "Minimum edge flip transport distance: "
     << config.minimum_transport_distance << "\n"
     << "Loading Settings:\n"
     << "  Start Load: " << config.startLoad << "\n"
     << "  Load Increment: " << config.loadIncrement << "\n"
     << "  Total Load: " << config.totalLoad << "\n"
     << "Minimizer: " << config.minimizer << "\n"
     << "Preconditioner: " << config.preconditioner << "\n"
     << "Preconditioner update: " << config.preconditionerUpdate << "\n"
     << "Preconditioner displacement threshold: "
     << config.preconditionerDisplacementThreshold << "\n";
  if (config.minimizer == "LBFGS") {
    os << "LBFGS Settings:\n"
       << "  Number of Corrections: " << config.LBFGSNrCorrections << "\n"
       << "  Scale: " << config.LBFGSScale << "\n"
       << "  EpsR: " << config.epsR << "\n"
       << "  EpsG: " << config.LBFGSEpsg << "\n"
       << "  EpsF: " << config.LBFGSEpsf << "\n"
       << "  EpsX: " << config.LBFGSEpsx << "\n"
       << "  Max LBFGS Iterations: " << config.LBFGSMaxIterations << "\n";
  }
  if (config.minimizer == "CG") {
    os << "CG Settings:\n"
       << "  Scale: " << config.CGScale << "\n"
       << "  EpsR: " << config.epsR << "\n"
       << "  EpsG: " << config.CGEpsg << "\n"
       << "  EpsF: " << config.CGEpsf << "\n"
       << "  EpsX: " << config.CGEpsx << "\n"
       << "  Max LBFGS Iterations: " << config.CGMaxIterations << "\n";
  }

  if (config.minimizer == "FIRE") {
    os << "FIRE Settings:\n"
       << "  Time step Increment Factor (finc): " << config.finc << "\n"
       << "  Time step Decrement Factor (fdec): " << config.fdec << "\n"
       << "  Alpha Start: " << config.alphaStart << "\n"
       << "  Alpha Factor (falpha): " << config.falpha << "\n"
       << "  Time Step Start (dtStart): " << config.dtStart << "\n"
       << "  Max Time Step (dtMax): " << config.dtMax << "\n"
       << "  Min Time Step (dtMin): " << config.dtMin << "\n"
       << "  Max component Step (maxCompS): " << config.maxCompS << "\n"
       << "  EpsR: " << config.epsR << "\n"
       << "  Epsilon: " << config.eps << "\n"
       << "  Epsilon Relative (epsRel): " << config.epsRel << "\n"
       << "  Delta: " << config.delta << "\n"
       << "  Max FIRE Iterations: " << config.maxIt << "\n";
  }

  os << "Plasticity event threshold: " << config.plasticityEventThreshold
     << "\n"
     << "Energy drop threshold: " << config.energyDropThreshold << "\n"
     << "Show progress: " << config.showProgress << "\n"
     << "Log during minimization: " << config.logDuringMinimization << "\n"
     << "Full minimization logging: " << config.fullMinimizationLogging << "\n"
     << "Write debug VTUs: " << config.writeDebugVTUs << "\n"
     << "Write dumps: " << config.writeDumps << "\n"
     << "Write ordinary mesh VTUs: " << config.writeMeshVTUs << "\n"
     << "Save elastic reversibility states: "
     << config.saveElasticReversibilityStates << "\n"
     << "Maximum saved elastic reversibility states: "
     << config.maximumSavedElasticReversibilityStates << "\n"
     << "Number of VTU frames: " << config.nrVTUFrames << "\n"
     << "Config path: " << config.configPath << "\n";
  return os;
}

namespace fs = std::filesystem;

std::map<std::string, std::string> parseParams(const std::string &filename) {
  std::map<std::string, std::string> config;
  std::ifstream file(filename);
  std::string line;

  // Check if the filename is empty
  if (filename.empty()) {
    std::cerr << "Error: Filename is empty." << std::endl;
    throw std::invalid_argument("Filename cannot be empty.");
  }

  // Extract the file name from the path and assign it to "name"
  fs::path filePath(filename);
  fs::path absolutePath = fs::absolute(filePath); // Convert to absolute path

  if (!file) { // Check if the file was successfully opened
    throw std::runtime_error("File not found: " +
                             filename); // Throw an error if not
  }

  while (std::getline(file, line)) {
    // Remove comments (anything after '#')
    std::string::size_type commentPos = line.find('#');
    if (commentPos != std::string::npos) {
      line = line.substr(0, commentPos);
    }

    // Trim the line for whitespace before parsing
    line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());

    std::istringstream is_line(line);
    std::string key;
    if (std::getline(is_line, key, '=')) {
      std::string value;
      if (std::getline(is_line, value)) {
        config[key] = value;
      }
    }
  }

  return config;
}

template <typename T>
void getValue(const std::map<std::string, std::string> &configMap,
              const std::string &fullKey, T &variable,
              std::set<std::string> &recognized) {
  const std::string key = fullKey.substr(fullKey.find_last_of('.') + 1);
  recognized.insert(key);
  const auto found = configMap.find(key);
  if (found == configMap.end()) {
    return;
  }
  const std::string &value = found->second;
  std::size_t consumed = value.size();
  if constexpr (std::is_same_v<T, int>) {
    variable = std::stoi(value, &consumed);
  } else if constexpr (std::is_same_v<T, double>) {
    variable = std::stod(value, &consumed);
  } else if constexpr (std::is_same_v<T, std::string>) {
    variable = value;
  } else if constexpr (std::is_same_v<T, bool>) {
    if (value != "true" && value != "false" && value != "1" && value != "0") {
      throw std::invalid_argument(key + " must be true or false.");
    }
    variable = value == "true" || value == "1";
  }
  if (consumed != value.size()) {
    throw std::invalid_argument("Invalid numeric value for " + key + ": " + value);
  }
}

Config initializeConfig(const std::map<std::string, std::string> &configMap) {
  Config config;
  config.setDefaultValues();
  std::set<std::string> recognized;
  GET_VALUE(config.loading);
  GET_VALUE(config.deformation);
  GET_VALUE(config.boundary);
  GET_VALUE(config.directions);
  GET_VALUE(config.segmentLoad);
  GET_VALUE(config.repeat);
  GET_VALUE(config.stepNoise);
  GET_VALUE(config.reversibilityTolerance);
  GET_VALUE(config.checkpointInterval);
  GET_VALUE(config.stopAfterSteps);
  GET_VALUE(config.stopAfterPlasticSteps);
  GET_VALUE(config.stopAfterElasticSteps);
  GET_VALUE(config.stopAfterEnergyDrops);
  GET_VALUE(config.stopAfterEnergyIncreases);
  GET_VALUE(config.stopWhen);
  GET_VALUE(config.stepEnergyTolerance);
  GET_VALUE(config.name);
  GET_VALUE(config.rows);
  GET_VALUE(config.cols);
  GET_VALUE(config.reconnectionMethod);
  GET_VALUE(config.nrThreads);
  GET_VALUE(config.seed);
  GET_VALUE(config.QDSD);
  GET_VALUE(config.initialGuessNoise);
  GET_VALUE(config.affineInitialGuess);
  GET_VALUE(config.meshDiagonal);
  GET_VALUE(config.reconnectRevert);
  GET_VALUE(config.reconnectEdgeLocking);
  GET_VALUE(config.minimum_transport_distance);
  GET_VALUE(config.energyFunction);
  GET_VALUE(config.bulkModulus);
  GET_VALUE(config.startLoad);
  GET_VALUE(config.loadIncrement);
  GET_VALUE(config.totalLoad);
  GET_VALUE(config.minimizer);
  GET_VALUE(config.preconditioner);
  GET_VALUE(config.preconditionerUpdate);
  GET_VALUE(config.preconditionerDisplacementThreshold);
  GET_VALUE(config.epsR);
  GET_VALUE(config.LBFGSNrCorrections);
  GET_VALUE(config.LBFGSScale);
  GET_VALUE(config.LBFGSEpsg);
  GET_VALUE(config.LBFGSEpsf);
  GET_VALUE(config.LBFGSEpsx);
  GET_VALUE(config.LBFGSMaxIterations);
  GET_VALUE(config.CGScale);
  GET_VALUE(config.CGEpsg);
  GET_VALUE(config.CGEpsf);
  GET_VALUE(config.CGEpsx);
  GET_VALUE(config.CGMaxIterations);
  GET_VALUE(config.finc);
  GET_VALUE(config.fdec);
  GET_VALUE(config.alphaStart);
  GET_VALUE(config.falpha);
  GET_VALUE(config.dtStart);
  GET_VALUE(config.maxCompS);
  GET_VALUE(config.dtMax);
  GET_VALUE(config.dtMin);
  GET_VALUE(config.eps);
  GET_VALUE(config.epsRel);
  GET_VALUE(config.delta);
  GET_VALUE(config.maxIt);
  GET_VALUE(config.logDuringMinimization);
  GET_VALUE(config.fullMinimizationLogging);
  GET_VALUE(config.writeDebugVTUs);
  GET_VALUE(config.writeDumps);
  GET_VALUE(config.writeMeshVTUs);
  GET_VALUE(config.saveElasticReversibilityStates);
  GET_VALUE(config.maximumSavedElasticReversibilityStates);
  GET_VALUE(config.nrVTUFrames);
  GET_VALUE(config.plasticityEventThreshold);
  GET_VALUE(config.energyDropThreshold);
  GET_VALUE(config.showProgress);
  GET_VALUE(config.forceReRun);

  for (const auto &[key, value] : configMap) {
    if (recognized.count(key) == 0) {
      throw std::invalid_argument("Unknown configuration key: " + key +
          ". Protocols use loading, deformation, boundary, and directions.");
    }
  }
  validateProtocolConfig(config);
  config.validatePreconditioner();
  return config;
}

Config parseConfigFile(std::string configPath) {
  auto confMap = parseParams(configPath);
  Config conf = initializeConfig(confMap);
  conf.configPath = configPath;
  return conf;
}
