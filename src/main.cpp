#include "Data/data_export.h"
#include "Data/param_parser.h"
#include "Simulation/experiments.h"
#include "Simulation/simulation.h"
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unistd.h>

namespace {

class DumpLoadError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

std::string trim(const std::string &str) {
  size_t first = str.find_first_not_of(' ');
  if (std::string::npos == first) {
    return str;
  }
  size_t last = str.find_last_not_of(' ');
  return str.substr(first, (last - first + 1));
}

void handleInputArgs(int argc, char *argv[]) {
  std::string configPath;
  std::string dumpPath;
  std::string outputPath;
  // The program will check if the folder already contains a completed
  // simulation. If the simulation is complete, the program will terminate
  // and not rerun the simulation unless forceReRun is true.
  bool forceReRun = false;
  double makeDumpAt = -1.0;

  int opt;
  const option longOptions[] = {
      {"makeDumpAt", required_argument, nullptr, 1000},
      {nullptr, 0, nullptr, 0},
  };
  while ((opt = getopt_long(argc, argv, "c:d:o:r", longOptions, nullptr)) !=
         -1) {
    switch (opt) {
    case 'c':
      configPath = trim(optarg);
      break;
    case 'd':
      dumpPath = trim(optarg);
      break;
    case 'o':
      outputPath = trim(optarg);
      break;
    case 'r':
      forceReRun = true;
      break;
    case 1000:
      makeDumpAt = std::stod(trim(optarg));
      break;
    default:
      std::cerr << "Usage: " << argv[0]
                << " -c <Config File> -d <Dump File> [-o <Output Path>] [-r] "
                   "[--makeDumpAt <load>]\n";
      return;
    }
  }

  // Check if dumpPath is provided
  if (!dumpPath.empty()) {
    auto sPtr = std::make_shared<Simulation>(); // Create a shared pointer to
                                                // a new Simulation object

    if (!configPath.empty()) { // Check if configPath is provided
      std::cout << "Overwriting simulation settings using\n - " << configPath
                << '\n';
    } else {
      // Try to find configPath in the same folder as dumpPath
      configPath = searchForConfig(dumpPath);
      if (configPath.empty()) { // If no configPath is found
        // Technically, we don't need a config file, since we can use the
        // settings from the dump, but to force the user to always keep a
        // config file in the folder, we throw an error here by design.
        std::cerr << "Error! No config provided or found in the same folder as "
                     "the dump file.\n";
        return; // Exit the function
      }
    }
    // Load and resume the simulation using the provided or found configPath
    try {
      Simulation::loadSimulation(*sPtr, dumpPath, configPath, outputPath,
                                 forceReRun);
    } catch (const cereal::Exception &error) {
      throw DumpLoadError("Failed to load dump " + dumpPath + ": " +
                          error.what());
    }

    std::cout << "Resuming simulation using " << dumpPath << '\n'
              << " - Config File: " << sPtr->config.name << '\n'
              << " - Data Path: " << sPtr->dataPath << '\n'
              << sPtr->config << '\n'
              << " - Current Load: " << sPtr->mesh.load << '\n';
    if (makeDumpAt > sPtr->mesh.load) {
      std::cout << "Making dump at: " << makeDumpAt << '\n';
    }
    std::cout << std::endl;
    runSimulationExperiment(sPtr->config, sPtr->dataPath, sPtr,
                            makeDumpAt); // Run the simulation experiment

    // If dumpPath is not provided but configPath is
  } else {
    if (configPath.empty()) {
      // Try to find configPath in the same folder as dumpPath
      configPath = searchForConfig(dumpPath);
      if (configPath.empty()) { // If no configPath is found
        // Technically, we don't need a config file, since we can use the
        // settings from the dump, but to force the user to always keep a
        // config file in the folder, we throw an error here by design.
        std::cerr << "Error! No config provided or found in the same folder as "
                     "the dump file.\n";
        return; // Exit the function
      }
    }
    Config config = parseConfigFile(configPath); // Parse the configuration file
    config.forceReRun = forceReRun; // Set the forceReRun flag in the config

    if (outputPath.empty()) {        // If outputPath is not provided
      outputPath = findOutputPath(); // Find the output path
    }

    // Run the simulation with the provided configuration
    std::cout << "Running simulation with:\n"
              << " - Config File: " << configPath << '\n'
              << " - Data Path: " << outputPath << '\n'
              << config << std::endl;

    runSimulationExperiment(config, outputPath, nullptr,
                            makeDumpAt); // Run the simulation experiment
  }
}

} // namespace

int main(int argc, char *argv[]) {
  try {
    // This starts the simulation.
    handleInputArgs(argc, argv);
    // If the simulation is already complete,
    // the program will terminate:
  } catch (const SimulationAlreadyComplete &e) {
    std::cout << e.what() << std::endl;
    return EXIT_SUCCESS;
  } catch (const DumpLoadError &error) {
    std::cerr << error.what() << std::endl;
    return 2;
  }
}
