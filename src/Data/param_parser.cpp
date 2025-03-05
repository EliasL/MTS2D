#include "param_parser.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
// Define a macro to get the variable name as a string
#define GET_VALUE(configMap, var, defaultValue)                                \
  getValue(configMap, #var, var, defaultValue)

void Config::setDefaultValues() {
  // General settings
  name = "default_name";
  rows = 3;
  cols = 3;
  usingPBC = false;
  scenario = "simpleShear";
  nrThreads = 1;
  seed = 0;
  QDSD = 0.0;
  initialGuessNoise = 0.05;
  meshDiagonal = "major";

  // Loading settings
  startLoad = 0.0;
  loadIncrement = 0.1;
  maxLoad = 1.0;

  // Minimizer settings
  minimizer = "LBFGS";

  // Max residual force allowed. Minimize until all forces are smaller than this
  epsR = 1e-7;

  // LBFGS-specific settings
  LBFGSNrCorrections = 10;
  LBFGSScale = 1.0;
  LBFGSEpsg = 1e-6;
  LBFGSEpsf = 1e-6;
  LBFGSEpsx = 1e-6;
  LBFGSMaxIterations = 0;

  // CG-specific settings
  CGScale = 1.0;
  CGEpsg = 1e-6;
  CGEpsf = 1e-6;
  CGEpsx = 1e-6;
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
  eps = 1e-6;
  epsRel = 0;
  delta = 0;
  maxIt = 100000;

  // Additional settings
  logDuringMinimization = false;
  plasticityEventThreshold = 0.01;
  energyDropThreshold = 0.001;
  showProgress = true;

  // Config path default
  configPath = "";
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
     << "Boundary Conditions: " << (config.usingPBC ? "PBC" : "NPBC") << "\n"
     << "Scenario: " << config.scenario << "\n"
     << "Number of Threads: " << config.nrThreads << "\n"
     << "Seed: " << config.seed << "\n"
     << "Quenched disorder standard deviation: " << config.QDSD << "\n"
     << "Initial guess noise: " << config.initialGuessNoise << "\n"
     << "Mesh diagonal: " << config.meshDiagonal << "\n"
     << "Loading Settings:\n"
     << "  Start Load: " << config.startLoad << "\n"
     << "  Load Increment: " << config.loadIncrement << "\n"
     << "  Max Load: " << config.maxLoad << "\n"
     << "Minimizer: " << config.minimizer << "\n";
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
       << "  Max FIRE Iterations: " << config.delta << "\n";
  }

  os << "Plasticity event threshold: " << config.plasticityEventThreshold
     << "\n"
     << "Energy drop threshold: " << config.energyDropThreshold << "\n"
     << "Show progress: " << config.showProgress << "\n"
     << "Log during minimization: " << config.logDuringMinimization << "\n"
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
  std::cout << "Attempting to open file: " << absolutePath
            << std::endl; // Debug print

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

// We use a macro to automatically insert the string key variable
template <typename T>
void getValue(const std::map<std::string, std::string> &configMap,
              const std::string &fullKey, T &variable, T defaultValue) {

  // The fullKet will be on the form config.something (see initializeConfig)
  std::string key = fullKey.substr(fullKey.find_last_of('.') + 1);

  // Check if the key exists in the map before accessing it
  auto it = configMap.find(key);
  if (it == configMap.end()) {
    std::cout << "Warning: Missing key '" << key << "' in configMap.\n"
              << "\tUsing default value: " << defaultValue << std::endl;
    variable = defaultValue;
    return;
  }

  // Now, safely access the value from configMap
  const std::string &value = it->second;

  if constexpr (std::is_same_v<T, int>) {
    variable = std::stoi(value);
  } else if constexpr (std::is_same_v<T, double>) {
    variable = std::stod(value);
  } else if constexpr (std::is_same_v<T, std::string>) {
    variable = value;
  } else if constexpr (std::is_same_v<T, bool>) {
    std::string val = value;
    std::transform(val.begin(), val.end(), val.begin(), ::tolower);
    variable = (val == "true" || val == "1");
  } else {
    throw std::runtime_error("Unsupported type");
  }
}

Config initializeConfig(const std::map<std::string, std::string> &configMap) {
  Config config;
  // We use a macro to automatically copy the variable name and use it as a key
  GET_VALUE(configMap, config.name, std::string(""));
  GET_VALUE(configMap, config.rows, 0);
  GET_VALUE(configMap, config.cols, 0);
  GET_VALUE(configMap, config.usingPBC, true);
  GET_VALUE(configMap, config.scenario, std::string(""));
  GET_VALUE(configMap, config.nrThreads, 0);
  GET_VALUE(configMap, config.seed, 0);
  GET_VALUE(configMap, config.QDSD, 0.0);
  GET_VALUE(configMap, config.initialGuessNoise, 0.0);
  GET_VALUE(configMap, config.meshDiagonal, std::string("major"));

  GET_VALUE(configMap, config.startLoad, 0.0);
  GET_VALUE(configMap, config.loadIncrement, 0.0);
  GET_VALUE(configMap, config.maxLoad, 0.0);

  GET_VALUE(configMap, config.minimizer, std::string(""));
  GET_VALUE(configMap, config.epsR, 0.0);

  GET_VALUE(configMap, config.LBFGSNrCorrections, 0);
  GET_VALUE(configMap, config.LBFGSScale, 0.0);
  GET_VALUE(configMap, config.LBFGSEpsg, 0.0);
  GET_VALUE(configMap, config.LBFGSEpsf, 0.0);
  GET_VALUE(configMap, config.LBFGSEpsx, 0.0);
  GET_VALUE(configMap, config.LBFGSMaxIterations, 0);

  GET_VALUE(configMap, config.CGScale, 0.0);
  GET_VALUE(configMap, config.CGEpsg, 0.0);
  GET_VALUE(configMap, config.CGEpsf, 0.0);
  GET_VALUE(configMap, config.CGEpsx, 0.0);
  GET_VALUE(configMap, config.CGMaxIterations, 0);

  GET_VALUE(configMap, config.finc, 0.0);
  GET_VALUE(configMap, config.fdec, 0.0);
  GET_VALUE(configMap, config.alphaStart, 0.0);
  GET_VALUE(configMap, config.falpha, 0.0);
  GET_VALUE(configMap, config.dtStart, 0.0);
  GET_VALUE(configMap, config.maxCompS, 0.0);
  GET_VALUE(configMap, config.dtMax, 0.0);
  GET_VALUE(configMap, config.dtMin, 0.0);
  GET_VALUE(configMap, config.eps, 0.0);
  GET_VALUE(configMap, config.epsRel, 0.0);
  GET_VALUE(configMap, config.delta, 0.0);
  GET_VALUE(configMap, config.maxIt, 0);

  GET_VALUE(configMap, config.logDuringMinimization, false);
  GET_VALUE(configMap, config.plasticityEventThreshold, 0.0);
  GET_VALUE(configMap, config.energyDropThreshold, 0.0);
  GET_VALUE(configMap, config.showProgress, 1);

  return config;
}

Config parseConfigFile(std::string configPath) {
  auto confMap = parseParams(configPath);
  Config conf = initializeConfig(confMap);
  conf.configPath = configPath;
  return conf;
}
