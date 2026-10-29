#include "Data/data_export.h"
#include "Data/param_parser.h"
#include "Simulation/protocols.h"
#include "Simulation/simulation.h"
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <filesystem>
#include <sys/file.h>
#include <fcntl.h>
#include <fstream>
#include <cereal/archives/json.hpp>

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

class RunLock {
  int fd = -1;
public:
  explicit RunLock(const std::filesystem::path &folder) {
    fd = open((folder / ".run.lock").c_str(), O_CREAT | O_RDWR, 0600);
    if (fd < 0 || flock(fd, LOCK_EX | LOCK_NB) != 0) {
      if (fd >= 0) close(fd);
      throw std::runtime_error("Cannot acquire simulation run lock: " + folder.string());
    }
  }
  ~RunLock() { if (fd >= 0) close(fd); }
};

void runFolder(const std::string &path, bool resume, double makeDumpAt) {
  namespace fs = std::filesystem;
  const auto folder = fs::canonical(path);
  RunLock lock(folder);
  const auto conf = folder / "config.conf";
  Config config = parseConfigFile(conf.string());
  if (config.name != folder.filename().string())
    throw std::runtime_error("Configuration name differs from run folder");
  if (fs::exists(folder / "termination.json")) {
    StepTermination termination;
    std::ifstream stream(folder / "termination.json");
    cereal::JSONInputArchive archive(stream);
    archive(cereal::make_nvp("termination", termination));
    if (termination.targetsMet || simulationAlreadyComplete(config.name, folder.parent_path().string(), config.startLoad + config.totalLoad))
      throw SimulationAlreadyComplete("Simulation already complete");
  }
  std::shared_ptr<Simulation> latest;
  fs::path selected;
  std::string selectedState;
  bool conflicting = false;
  if (resume && fs::exists(folder / "dumps")) {
    for (const auto &entry : fs::directory_iterator(folder / "dumps")) {
      const auto name = entry.path().filename().string();
      if (!(name.rfind("checkpoint_", 0) == 0 || name.rfind("dump_l", 0) == 0) ||
          name.size() < 7 || name.substr(name.size() - 7) != ".xml.gz") continue;
      auto candidate = std::make_shared<Simulation>();
      loadFromFile(entry.path().string(), *candidate);
      if (!latest || candidate->mesh.loadSteps > latest->mesh.loadSteps) {
        latest = candidate;
        selected = entry.path();
        selectedState.clear();
        conflicting = false;
      } else if (candidate->mesh.loadSteps == latest->mesh.loadSteps) {
        if (selectedState.empty())
          selectedState = serializeToXml(latest->mesh) + serializeToXml(latest->termination) + serializeToXml(latest->config);
        const auto state = serializeToXml(candidate->mesh) + serializeToXml(candidate->termination) + serializeToXml(candidate->config);
        if (state != selectedState) conflicting = true;
      }
    }
  }
  if (conflicting) throw std::runtime_error("Conflicting checkpoints at the latest loading step");
  if (latest) {
    std::cout << "Resuming checkpoint " << selected << " step=" << latest->mesh.loadSteps
              << " load=" << latest->mesh.load << std::endl;
    Simulation::loadSimulation(*latest, selected.string(), conf.string(), folder.parent_path().string(), false);
    runSimulationProtocol(latest->config, latest->dataPath, latest, makeDumpAt);
  } else {
    if (fs::exists(folder / "macroData.csv"))
      throw std::runtime_error("Existing scientific output requires --resume and a completed checkpoint");
    runSimulationProtocol(config, folder.parent_path().string(), nullptr, makeDumpAt);
  }
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
  bool resume = false;

  int opt;
  const option longOptions[] = {
      {"makeDumpAt", required_argument, nullptr, 1000},
      {"resume", no_argument, nullptr, 1001},
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
    case 1001:
      resume = true;
      break;
    case 1000:
      makeDumpAt = std::stod(trim(optarg));
      break;
    default:
      std::cerr << "Usage: " << argv[0]
                << " -c <Config File> -d <Dump File> [-o <Output Path>] [-r] "
                   "[--makeDumpAt <load>]\n";
      throw std::invalid_argument("Invalid command-line arguments");
    }
  }
  if (optind < argc) {
    if (optind + 1 != argc || !configPath.empty() || !dumpPath.empty() || !outputPath.empty())
      throw std::invalid_argument("Use one run folder or explicit -c/-d/-o arguments");
    runFolder(argv[optind], resume, makeDumpAt);
    return;
  }
  if (resume) throw std::invalid_argument("--resume requires a run folder");

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
        throw std::invalid_argument("No configuration provided or found");
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
    runSimulationProtocol(sPtr->config, sPtr->dataPath, sPtr,
                            makeDumpAt); // Run the simulation protocol

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
        throw std::invalid_argument("No configuration provided or found");
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

    runSimulationProtocol(config, outputPath, nullptr,
                            makeDumpAt); // Run the simulation protocol
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
  } catch (const std::exception &error) {
    std::cerr << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
