#include "paramParser.h"
#include <filesystem>
#include <fstream>
#include <sstream>
// Define a macro to get the variable name as a string
#define GET_VALUE(configMap, var) getValue(configMap, #var, var)

void Config::updateParam(FIREpp::FIREParam<double> &param) {
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
     << "\n"
     << "Loading Settings:\n"
     << "  Start Load: " << config.startLoad << "\n"
     << "  Load Increment: " << config.loadIncrement << "\n"
     << "  Max Load: " << config.maxLoad << "\n"
     << "Minimizer: " << config.minimizer << "\n"
     << "LBFGS Settings:\n"
     << "  Number of Corrections: " << config.LBFGSNrCorrections << "\n"
     << "  Scale: " << config.LBFGSScale << "\n"
     << "  EpsG: " << config.LBFGSEpsg << "\n"
     << "  EpsF: " << config.LBFGSEpsf << "\n"
     << "  EpsX: " << config.LBFGSEpsx << "\n"
     << "  Max LBFGS Iterations: " << config.LBFGSMaxIterations << "\n";

  if (config.minimizer != "LBFGS") {
    os << "FIRE Settings:\n"
       << "  Time step Increment Factor (finc): " << config.finc << "\n"
       << "  Time step Decrement Factor (fdec): " << config.fdec << "\n"
       << "  Alpha Start: " << config.alphaStart << "\n"
       << "  Alpha Factor (falpha): " << config.falpha << "\n"
       << "  Time Step Start (dtStart): " << config.dtStart << "\n"
       << "  Max Time Step (dtMax): " << config.dtMax << "\n"
       << "  Min Time Step (dtMin): " << config.dtMin << "\n"
       << "  Max component Step (maxCompS): " << config.maxCompS << "\n"
       << "  Epsilon: " << config.eps << "\n"
       << "  Epsilon Relative (epsRel): " << config.epsRel << "\n"
       << "  Delta: " << config.delta << "\n"
       << "  Max FIRE Iterations: " << config.delta << "\n";
  }

  os << "Plasticity Event Threshold: " << config.plasticityEventThreshold
     << "\n"
     << "Show Progress: " << config.showProgress << "\n"
     << "Config Path: " << config.configPath << "\n";
  return os;
}

namespace fs = std::filesystem;

std::map<std::string, std::string> parseParams(const std::string &filename) {
  std::map<std::string, std::string> config;
  std::ifstream file(filename);
  std::string line;

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
              const std::string &fullKey, T &variable) {
  std::string key = fullKey.substr(fullKey.find_last_of('.') + 1);
  try {
    if constexpr (std::is_same_v<T, int>) {
      variable = std::stoi(configMap.at(key));
    } else if constexpr (std::is_same_v<T, double>) {
      variable = std::stod(configMap.at(key));
    } else if constexpr (std::is_same_v<T, std::string>) {
      variable = configMap.at(key);
    } else if constexpr (std::is_same_v<T, bool>) {
      variable = std::stoi(configMap.at(key)) == 1;
    } else {
      throw std::runtime_error("Unsupported type");
    }
  } catch (const std::out_of_range &) {
    std::cerr << "Error: Missing key '" << key << "' in configMap."
              << std::endl;
    throw;
  } catch (const std::invalid_argument &) {
    std::cerr << "Error: Invalid value for key '" << key << "' in configMap."
              << std::endl;
    throw;
  }
}

Config initializeConfig(const std::map<std::string, std::string> &configMap) {
  Config config;
  // We use a macro to automatically copy the variable name and use it as a key
  GET_VALUE(configMap, config.name);
  GET_VALUE(configMap, config.rows);
  GET_VALUE(configMap, config.cols);
  GET_VALUE(configMap, config.usingPBC);
  GET_VALUE(configMap, config.scenario);
  GET_VALUE(configMap, config.nrThreads);
  GET_VALUE(configMap, config.seed);
  GET_VALUE(configMap, config.QDSD);
  GET_VALUE(configMap, config.initialGuessNoise);

  GET_VALUE(configMap, config.startLoad);
  GET_VALUE(configMap, config.loadIncrement);
  GET_VALUE(configMap, config.maxLoad);

  GET_VALUE(configMap, config.minimizer);
  GET_VALUE(configMap, config.LBFGSNrCorrections);
  GET_VALUE(configMap, config.LBFGSScale);
  GET_VALUE(configMap, config.LBFGSEpsg);
  GET_VALUE(configMap, config.LBFGSEpsf);
  GET_VALUE(configMap, config.LBFGSEpsx);
  GET_VALUE(configMap, config.LBFGSMaxIterations);

  GET_VALUE(configMap, config.finc);
  GET_VALUE(configMap, config.fdec);
  GET_VALUE(configMap, config.alphaStart);
  GET_VALUE(configMap, config.falpha);
  GET_VALUE(configMap, config.dtStart);
  GET_VALUE(configMap, config.maxCompS);
  GET_VALUE(configMap, config.dtMax);
  GET_VALUE(configMap, config.dtMin);
  GET_VALUE(configMap, config.eps);
  GET_VALUE(configMap, config.epsRel);
  GET_VALUE(configMap, config.delta);
  GET_VALUE(configMap, config.maxIt);

  GET_VALUE(configMap, config.plasticityEventThreshold);
  GET_VALUE(configMap, config.showProgress);

  return config;
}

Config parseConfigFile(std::string configPath) {
  auto confMap = parseParams(configPath);
  Config conf = initializeConfig(confMap);
  conf.configPath = configPath;
  return conf;
}
