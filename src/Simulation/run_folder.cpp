#include "run_folder.h"

#include "Data/data_export.h"
#include "Data/param_parser.h"
#include "protocols.h"
#include "simulation.h"
#include <algorithm>
#include <cereal/archives/json.hpp>
#include <cfloat>
#include <cmath>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <regex>
#include <stdexcept>
#include <sys/file.h>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {
namespace fs = std::filesystem;

class RunLock {
  int fd = -1;
public:
  explicit RunLock(const fs::path &folder) {
    fd = open((folder / ".run.lock").c_str(), O_CREAT | O_RDWR, 0600);
    if (fd < 0 || flock(fd, LOCK_EX | LOCK_NB) != 0) {
      if (fd >= 0) close(fd);
      throw std::runtime_error("Cannot acquire simulation run lock: " + folder.string());
    }
  }
  RunLock(const RunLock &) = delete;
  RunLock &operator=(const RunLock &) = delete;
  ~RunLock() { if (fd >= 0) close(fd); }
};

struct LoadBounds {
  double lower, upper;
};

std::optional<LoadBounds> namedLoadBounds(const std::string &name) {
  static const std::regex pattern(R"(dump_l([+-]?[0-9]+(?:\.([0-9]+))?(?:[eE]([+-]?[0-9]+))?)\.xml\.gz)");
  std::smatch match;
  if (!std::regex_match(name, match, pattern)) return std::nullopt;
  try {
    const double load = std::stod(match[1].str());
    const double exponent = match[3].matched ? std::stod(match[3].str()) : 0.0;
    const double unit = std::pow(10.0, exponent - match[2].length());
    const double padding = unit / 2 + 8 * DBL_EPSILON * std::abs(load);
    if (std::isfinite(load - padding) && std::isfinite(load + padding))
      return LoadBounds{load - padding, load + padding};
  } catch (const std::out_of_range &) {
  }
  return std::nullopt;
}

struct Checkpoint {
  fs::path path;
  std::shared_ptr<Simulation> state;
};

std::string checkpointState(const Simulation &s) {
  return serializeToXml(s.mesh) + serializeToXml(s.termination) +
         serializeToXml(s.config);
}

Checkpoint latestCheckpoint(const fs::path &dumps, std::optional<double> lastLoad) {
  Checkpoint selected;
  auto &latest = selected.state;
  std::string selectedState;
  bool conflicting = false;
  struct Candidate {
    fs::path path;
    std::optional<LoadBounds> bounds;
  };
  std::vector<Candidate> candidates;
  for (const auto &entry : fs::directory_iterator(dumps)) {
    const auto name = entry.path().filename().string();
    if (!(name.rfind("checkpoint_", 0) == 0 || name.rfind("dump_l", 0) == 0) ||
        name.size() < 7 || name.substr(name.size() - 7) != ".xml.gz") continue;
    candidates.push_back({entry.path(), namedLoadBounds(name)});
  }
  if (lastLoad) {
    auto priority = [&](const Candidate &candidate) {
      if (!candidate.bounds) return std::pair{1, 0.0};
      const auto [lower, upper] = *candidate.bounds;
      return std::pair{lower <= *lastLoad && *lastLoad <= upper ? 0 : 2, -upper};
    };
    std::sort(candidates.begin(), candidates.end(), [&](const auto &left, const auto &right) {
      return priority(left) < priority(right);
    });
  }

  std::vector<std::pair<int, double>> inspected;
  std::vector<const Candidate *> skipped;
  bool namesMatch = true;
  auto inspect = [&](const Candidate &file) {
    auto candidate = std::make_shared<Simulation>();
    loadFromFile(file.path.string(), *candidate);
    const double load = candidate->mesh.load;
    inspected.emplace_back(candidate->mesh.loadSteps, load);
    if (file.bounds && !(file.bounds->lower <= load && load <= file.bounds->upper))
      namesMatch = false;
    if (!latest || candidate->mesh.loadSteps > latest->mesh.loadSteps) {
      latest = candidate;
      selected.path = file.path;
      selectedState.clear();
      conflicting = false;
    } else if (candidate->mesh.loadSteps == latest->mesh.loadSteps) {
      if (selectedState.empty()) selectedState = checkpointState(*latest);
      if (checkpointState(*candidate) != selectedState) conflicting = true;
    }
  };
  for (const auto &file : candidates) {
    if (lastLoad && namesMatch && latest && file.bounds && file.bounds->upper < latest->mesh.load) {
      skipped.push_back(&file);
    } else {
      inspect(file);
    }
  }
  // Native dump names encode the increasing load coordinate, even when the
  // physical deformation reverses. Check that contract for the opened files.
  std::sort(inspected.begin(), inspected.end());
  const auto inconsistent = [](const auto &left, const auto &right) {
    return left.first == right.first ? left.second != right.second : left.second >= right.second;
  };
  if (!namesMatch || std::adjacent_find(inspected.begin(), inspected.end(), inconsistent) != inspected.end()) {
    for (const auto *file : skipped) inspect(*file);
  }
  if (conflicting)
    throw std::runtime_error("Conflicting checkpoints at the latest loading step");
  if (latest)
    std::cout << "Checked " << inspected.size() << " of " << candidates.size()
              << " checkpoint files." << std::endl;
  return selected;
}
} // namespace

void runSimulationFolder(const std::string &path, bool resume, double makeDumpAt) {
  const auto folder = fs::canonical(path);
  RunLock lock(folder);
  const auto conf = folder / "config.conf";
  Config config = parseConfigFile(conf.string());
  if (config.name != folder.filename().string())
    throw std::runtime_error("Configuration name differs from run folder");
  const auto outputRoot = folder.parent_path().string();
  const auto macro = folder / "macroData.csv";
  const auto summary = folder / "termination.json";
  if (fs::exists(summary)) {
    StepTermination termination;
    std::ifstream stream(summary);
    cereal::JSONInputArchive archive(stream);
    archive(cereal::make_nvp("termination", termination));
    if (termination.targetsMet ||
        simulationAlreadyComplete(config.name, outputRoot, config.startLoad + config.totalLoad))
      throw SimulationAlreadyComplete("Simulation already complete");
  }

  const auto dumps = folder / "dumps";
  if (resume && fs::exists(dumps)) {
    const auto checkpoint = latestCheckpoint(dumps, lastCompleteMacroLoad(macro.string()));
    if (const auto &state = checkpoint.state) {
      std::cout << "Resuming checkpoint " << checkpoint.path
                << " step=" << state->mesh.loadSteps << " load=" << state->mesh.load << std::endl;
      Simulation::loadSimulation(*state, checkpoint.path.string(), conf.string(), outputRoot, false);
      runSimulationProtocol(state->config, state->dataPath, state, makeDumpAt);
      return;
    }
  }
  if (fs::exists(macro))
    throw std::runtime_error("Existing scientific output requires --resume and a completed checkpoint");
  runSimulationProtocol(config, outputRoot, nullptr, makeDumpAt);
}
