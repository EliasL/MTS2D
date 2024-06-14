#include "paramParser.h"
#include <filesystem>
#include <fstream>
#include <sstream>

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
     << "  Max LBFGS Iterations: " << config.LBFGSMaxIterations << "\n"
     << "FIRE Settings:\n"
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
     << "  Max FIRE Iterations: " << config.delta << "\n"
     << "Plasticity Event Threshold: " << config.plasticityEventThreshold
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

Config initializeConfig(const std::map<std::string, std::string> &configMap) {
  Config config;

  // Simulation and General Settings
  config.name = configMap.at("name");
  config.rows = std::stoi(configMap.at("rows"));
  config.cols = std::stoi(configMap.at("cols"));
  config.usingPBC = std::stoi(configMap.at("usingPBC")) == 1;
  config.scenario = configMap.at("scenario");
  config.nrThreads = std::stoi(configMap.at("nrThreads"));
  config.seed = std::stoi(configMap.at("seed"));
  config.QDSD = std::stod(configMap.at("QDSD"));
  config.initialGuessNoise = std::stod(configMap.at("initialGuessNoise"));

  // Loading Settings
  config.startLoad = std::stod(configMap.at("startLoad"));
  config.loadIncrement = std::stod(configMap.at("loadIncrement"));
  config.maxLoad = std::stod(configMap.at("maxLoad"));

  // Minimizer Settings
  config.minimizer = configMap.at("minimizer");
  // Specific settings for LBFGS
  config.LBFGSNrCorrections = std::stoi(configMap.at("LBFGSNrCorrections"));
  config.LBFGSScale = std::stod(configMap.at("LBFGSScale"));
  config.LBFGSEpsg = std::stod(configMap.at("LBFGSEpsg"));
  config.LBFGSEpsf = std::stod(configMap.at("LBFGSEpsf"));
  config.LBFGSEpsx = std::stod(configMap.at("LBFGSEpsx"));
  config.LBFGSMaxIterations = std::stoi(configMap.at("LBFGSMaxIterations"));

  // Specific settings for FIRE (if used in the project)
  config.finc = std::stod(configMap.at("finc"));
  config.fdec = std::stod(configMap.at("fdec"));
  config.alphaStart = std::stod(configMap.at("alphaStart"));
  config.falpha = std::stod(configMap.at("falpha"));
  config.dtStart = std::stod(configMap.at("dtStart"));
  config.maxCompS = std::stod(configMap.at("maxCompS"));
  config.dtMax = std::stod(configMap.at("dtMax"));
  config.dtMin = std::stod(configMap.at("dtMin"));
  config.eps = std::stod(configMap.at("eps"));
  config.epsRel = std::stod(configMap.at("epsRel"));
  config.delta = std::stod(configMap.at("delta"));
  config.maxIt = std::stoi(configMap.at("maxIt"));

  // Logging Settings
  config.plasticityEventThreshold =
      std::stod(configMap.at("plasticityEventThreshold"));
  config.showProgress = std::stoi(configMap.at("showProgress"));

  return config;
}

Config parseConfigFile(std::string configPath) {
  auto confMap = parseParams(configPath);
  Config conf = initializeConfig(confMap);
  conf.configPath = configPath;
  return conf;
}
