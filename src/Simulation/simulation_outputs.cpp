#include "simulation.h"
#include "Data/cereal_help.h"
#include "Data/data_export.h"
#include "Data/param_parser.h"
#include "settings.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
std::string defaultDumpName(double load, double startLoad, double maxLoad) {
  double range = std::abs(maxLoad - startLoad);
  int precision =
      std::max(1, static_cast<int>(std::ceil(-std::log10(range / 10.0))));
  double roundedLoad = std::round(load * std::pow(10.0, precision)) /
                       std::pow(10.0, precision);

  std::ostringstream oss;
  oss << std::fixed << std::setprecision(precision) << roundedLoad;
  return "dump_l" + oss.str();
}

std::string exactLoadDumpName(double load) {
  std::ostringstream oss;
  oss << std::setprecision(std::numeric_limits<double>::max_digits10) << load;
  return "dump_l" + oss.str();
}
} // namespace

void Simulation::addCsvColumnRaw(std::string name, CsvGetter getter) {
  csvColumns.push_back({std::move(name), std::move(getter)});
}

bool Simulation::hasCsvColumn(const std::string &name) const {
  for (const auto &col : csvColumns) {
    if (col.name == name) {
      return true;
    }
  }
  return false;
}

bool Simulation::tryRecoverCsvColumn(const std::string &name) {
  if (name == "is_reversible" || name == "rev_d") {
    addReversibilityCsvColumns();
    return true;
  }
  return false;
}

static std::vector<std::string> splitCsvHeaderLine(const std::string &line) {
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) {
    out.push_back(item);
  }
  return out;
}

void Simulation::recoverCsvColumnsFromFile(const std::string &csvPath) {
  std::ifstream file(csvPath);
  if (!file.is_open()) {
    return;
  }
  std::string header;
  if (!std::getline(file, header)) {
    return;
  }
  const auto headers = splitCsvHeaderLine(header);
  for (const auto &name : headers) {
    if (!hasCsvColumn(name)) {
      tryRecoverCsvColumn(name);
    }
  }
}

// Default CSV columns. Use (expr) to keep commas inside a single macro arg.
#define DEFAULT_CSV_COLS(X)                                                    \
  X("load_step", s.mesh.loadSteps)                                             \
  X("load", s.mesh.load)                                                       \
  X("total_energy", s.mesh.totalEnergy)                                        \
  X("total_energy_change", s.energyHistory.loadStepTotalEnergyChange)          \
  X("total_init_energy", s.energyHistory.initialGuessTotalEnergy)              \
  X("total_e_change_from_init",                                                \
    s.energyHistory.totalEnergyChangeFromInitialGuess)                         \
  X("avg_energy", s.mesh.averageEnergy)                                        \
  X("avg_energy_change", s.energyHistory.loadStepAverageEnergyChange)          \
  X("avg_init_energy", s.energyHistory.initialGuessAverageEnergy)              \
  X("avg_e_change_from_init",                                                  \
    s.energyHistory.averageEnergyChangeFromInitialGuess)                       \
  X("min_iter_total_energy_change", s.energyHistory.minIterTotalEnergyChange)  \
  X("min_iter_avg_energy_change", s.energyHistory.minIterAverageEnergyChange)  \
  X("max_energy", s.mesh.maxEnergy)                                            \
  X("max_force", s.mesh.maxForce)                                              \
  X("avg_sigma11", s.mesh.averageSigma11)                                      \
  X("avg_sigma12", s.mesh.averageSigma12)                                      \
  X("avg_sigma22", s.mesh.averageSigma22)                                      \
  X("avg_init_sigma11", s.energyHistory.initialGuessAverageSigma11)            \
  X("avg_init_sigma12", s.energyHistory.initialGuessAverageSigma12)            \
  X("avg_init_sigma22", s.energyHistory.initialGuessAverageSigma22)            \
  X("avg_sigma12_change_from_init",                                            \
    s.energyHistory.averageSigma12ChangeFromInitialGuess)                      \
  /* avg_P* are diagnostic element means. avg_P12 mixes material directions  \
     from independently oriented element reference maps and is not the       \
     energy-conjugate shear stress for the left-multiplicative affine step;   \
     prefer avg_sigma12 when J is approximately one. */                       \
  X("avg_P11", s.mesh.averageP11)                                              \
  X("avg_P12", s.mesh.averageP12)                                              \
  X("avg_P21", s.mesh.averageP21)                                              \
  X("avg_P22", s.mesh.averageP22)                                              \
  X("avg_init_P11", s.energyHistory.initialGuessAverageP11)                    \
  X("avg_init_P12", s.energyHistory.initialGuessAverageP12)                    \
  X("avg_init_P21", s.energyHistory.initialGuessAverageP21)                    \
  X("avg_init_P22", s.energyHistory.initialGuessAverageP22)                    \
  X("participationFraction", s.participationFraction)                          \
  X("m3_participationFraction", s.m3ParticipationFraction)                     \
  X("nr_elements_with_m3_change", s.mesh.nr_elements_with_m3_change)           \
  X("nr_red_q1", s.mesh.redQuadrantCounts[0])                                  \
  X("nr_red_q2", s.mesh.redQuadrantCounts[1])                                  \
  X("nr_red_q3", s.mesh.redQuadrantCounts[2])                                  \
  X("nr_red_q4", s.mesh.redQuadrantCounts[3])                                  \
  X("max_m3_nr", s.mesh.maxM3Nr)                                               \
  X("sum_m3", s.mesh.sumM3Nr)                                                  \
  X("max_positive_plastic_jump", s.mesh.maxPlasticJump)                        \
  X("max_negative_plastic_jump", s.mesh.minPlasticJump)                        \
  X("nr_iterations", s.mesh.nrMinItterations)                                  \
  X("nr_func_evals", s.mesh.nrMinFunctionCalls)                                \
  X("nr_edge_flips", s.mesh.edgeFlipsFromLastStep())                           \
  X("nr_total_edge_flips", s.mesh.totalEdgeFlipsInStep)                        \
  X("nr_reconnect_cycles", s.reconnectingCycles())                             \
  X("reconnect_stop_reason", s.reconnectStopReasonName())                      \
  X("rejected_reconnect_energy_delta",                                         \
    s.rejectedReconnectEnergyDelta())                                          \
  X("edge_flip_chosen_minus_other_energy",                                     \
    s.mesh.edgeFlipChosenMinusOtherEnergyInStep)                               \
  X("edge_flip_always_chose_lower_energy",                                     \
    s.mesh.edgeFlipAlwaysChoseLowerEnergyInStep)                               \
  X("LBFGS_Term_reason", s.LBFGSRep.termType)                                  \
  X("CG_Term_reason", s.CGRep.termType)                                        \
  X("FIRE_Term_reason", s.FIRERep.termType)                                    \
  X("run_time", (s.timer.RTString()))                                          \
  X("minimization_time", (s.timer.RTString("minimization", 7)))                \
  X("write_time", (s.timer.RTString("write", 7)))                              \
  X("est_time_remaining", s.timer.oldETRString)                                \
  X("cmX", s.mesh.com[0])                                                      \
  X("cmY", s.mesh.com[1])                                                      \
  X("maxX", s.mesh.bounds[0])                                                  \
  X("minX", s.mesh.bounds[1])                                                  \
  X("maxY", s.mesh.bounds[2])                                                  \
  X("minY", s.mesh.bounds[3])

void Simulation::addDefaultCsvColumns() {
  if (csvDefaultsAdded) {
    return;
  }
  csvDefaultsAdded = true;
// If you are trying to understand how to add columns, look at:
// addCsvColumn(name, [](const auto &s) { return (expr); })
// This can be used to add a column.
#define ADD_COL(name, expr)                                                    \
  addCsvColumn(name, [](const auto &s) { return (expr); });
  DEFAULT_CSV_COLS(ADD_COL)
#undef ADD_COL
}

#undef DEFAULT_CSV_COLS

void Simulation::m_updateProgress() {
  using namespace std::chrono;

  progress =
      std::clamp((mesh.load - startLoad) / (maxLoad - startLoad), 0.0, 1.0);

  int intProgress = static_cast<int>(progress * 100);

  std::string consoleProgressMessage =
      std::to_string(intProgress) + "%"         //
      + " RT: " + timer.RTString()              //
      + "\tETR: " + timer.ETRString(progress)   //
      + "\tLoad: " + std::to_string(mesh.load); //

  // Use static variables to track the last progress and the last update
  // time
  static int oldProgress = -1;
  static int firstProgress = -1;
  static int lastDump = -1;
  static auto lastUpdateTime = steady_clock::now();

  auto now = steady_clock::now();
  auto timeSinceLastUpdate =
      duration_cast<seconds>(now - lastUpdateTime).count();

  bool shouldUpdate =
      !isQuiet() && (oldProgress != intProgress || timeSinceLastUpdate >= 20);

  if (shouldUpdate) {
    oldProgress = intProgress;
    firstProgress = (firstProgress == -1) ? intProgress : firstProgress;
    lastUpdateTime = now;

    std::cout << consoleProgressMessage << std::endl;

    if (!csvFile.is_open()) {
      throw std::runtime_error("CSV file stream is not open.");
    }
    csvFile.flush(); // Update the file on disk

    // Check if we should make a dump
    if (intProgress % 5 == 0 && intProgress != lastDump &&
        firstProgress != intProgress) {
      lastDump = intProgress;
      forceDumpAfterStep = true;
    }
  }
}

void Simulation::writeToFile(bool forceWrite, std::string fileName) {
  timer.Start("write");
  // We write to the cvs file every time this function is called
  writeToCsv(csvFile, (*this));
  // These are writing date much less often
  m_writeMesh(forceWrite);
  if (forceWrite || !fileName.empty()) {
    m_writeDump(true, fileName);
  }
  if (config.logDuringMinimization) {
    // If we are logging minimization, we want to keep the folder only for
    // steps with enough plasticity or energy drop.
    const bool hasPlasticEvents =
        mesh.nr_elements_with_m3_change >
        mesh.nrElements * config.plasticityEventThreshold;
    const bool hasEnergyDrop =
        -energyHistory.totalEnergyChangeFromInitialGuess >
        config.energyDropThreshold;
    const bool periodicKeep = (mesh.loadSteps % 1000 == 0);
    const bool keepMinFolder =
        forceWrite || hasPlasticEvents || hasEnergyDrop || periodicKeep;

    if (minCsvSubfolder.empty()) {
      timer.Stop("write");
      return;
    }

    const std::string minFolderPath =
        getOutputPath(simName, dataPath) + "/" + minCsvSubfolder;

    if (keepMinFolder) {
      if (periodicKeep && !hasPlasticEvents && !hasEnergyDrop) {
        mesh.writeToVtu("", true, VtuFieldLevel::All, "periodicKeep");
      }
      // Create a collection after each step we keep.
      createCollection(
          minFolderPath, minFolderPath,
          "^.*_minStep=[0-9]+\\.([0-9]+)(?:_[^.]+)?\\.[0-9]+\\.vtu$", ".vtu",
          {}, COLLECTIONNAME, "", "_reference");
      createCollection(
          minFolderPath, minFolderPath,
          "^.*_minStep=[0-9]+\\.([0-9]+)(?:_[^.]+)?\\.[0-9]+\\.vtu$", ".vtu",
          {}, "reference_collection", "_reference", "");
    } else {
      namespace fs = std::filesystem;
      fs::path folderPath(minFolderPath);
      if (fs::exists(folderPath) && fs::is_directory(folderPath)) {
        if (minCsvFile.is_open()) {
          minCsvFile.close();
        }
        fs::remove_all(folderPath);
        minCsvSubfolder.clear();
      }
    }
  }
  timer.Stop("write");
}

void Simulation::discardFailedCurrentLoadStepCsvRow() {
  csvFile.flush();
  if (!csvFile) {
    throw std::runtime_error("Could not flush crash diagnostic CSV row.");
  }
  csvFile.close();
  if (csvFile.fail()) {
    throw std::runtime_error("Could not close CSV before trimming crash row.");
  }

  const std::filesystem::path csvPath =
      std::filesystem::path(getOutputPath(simName, dataPath)) /
      (std::string(MACRODATANAME) + ".csv");
  if (!trimCsvFile(csvPath.string(), *this, true)) {
    throw std::runtime_error(
        "Expected crash diagnostic CSV row was missing after retry.");
  }
  csvFile.open(csvPath, std::ios::app);
  if (!csvFile.is_open()) {
    throw std::runtime_error("Could not reopen CSV after trimming crash row.");
  }
}

void Simulation::m_writeMesh(bool forceWrite) {
  if (!config.writeMeshVTUs) {
    return;
  }
  // Only if there are lots of plastic events will we want to save the data.
  // If we save every frame, it requires too much storage.
  // (A 100x100 system loaded from 0.15 to 1 with steps of 1e-5 would take
  // up 180GB) At the same time, if there are few large avalanvhes, we might
  // go long without saving data. In order to get a good framerate for an
  // animation, we want to ensure that not too much happens between frames.
  // The default nrVTUFrames=200 preserves the old hard-coded 0.005 relative
  // load spacing, but short benchmarks can lower it in their generated config.
  if (config.nrVTUFrames <= 0) {
    throw std::runtime_error("nrVTUFrames must be positive.");
  }
  static double lastLoadWritten = 0;
  const double loadRange = maxLoad - startLoad;
  if (loadRange <= 0) {
    throw std::runtime_error("maxLoad must be larger than startLoad.");
  }
  const double relativeLoadChange =
      std::abs(mesh.load - lastLoadWritten) / loadRange;
  if ((mesh.nr_elements_with_m3_change >
       mesh.nrElements *
           config.plasticityEventThreshold) || // Lots of plastic change
      (-energyHistory.totalEnergyChangeFromInitialGuess >
       config.energyDropThreshold) ||                    // Large energy drop
      (abs(mesh.load - lastLoadWritten) > 0.005) ||      // Absolute change
      (relativeLoadChange > 1.0 / config.nrVTUFrames) || // Relative change
      forceWrite)                                        // Force write
  {
    mesh.writeToVtu();
    lastLoadWritten = mesh.load;
  }
}

void Simulation::m_writeDump(bool forceWrite, std::string name) {
  if (!config.writeDumps) {
    return;
  }
  // When do we create save states?
  // I'm thinking I want to do one halfway no matter how short the
  // simulation is, but then outside of that, i'm thinking once per hour is
  // okay.
  auto now = std::chrono::steady_clock::now();
  static auto lastSaveTime =
      now; // Since this is a static variable, this line is only run once
  static bool firstSaveDone = false;
  static int lastDefaultDumpLoadStep = -1;
  auto elapsedSinceLastSave =
      std::chrono::duration_cast<std::chrono::seconds>(now - lastSaveTime);
  double midPointLoad = (startLoad + maxLoad) / 2;

  // Save every one hours
  static const std::chrono::hours saveFrequency(1);
  const bool makeTargetDump =
      makeDumpAt >= 0.0 &&
      std::abs(mesh.load - makeDumpAt) <= 0.5 * std::abs(loadIncrement);

  bool shouldSave =
      (mesh.load >= midPointLoad &&
       !firstSaveDone) || // Check for first save at midpoint
      (elapsedSinceLastSave >= saveFrequency) || // Hourly save
      mesh.load + loadIncrement / 2 > maxLoad || // Check for final save
      makeTargetDump ||                          // User-requested load dump
      forceWrite;                                // Check if forced

  if (shouldSave) {
    const bool usingDefaultName = name.empty();
    if (usingDefaultName && lastDefaultDumpLoadStep == mesh.loadSteps) {
      return;
    }
    const std::string dumpName = usingDefaultName && makeTargetDump
                                     ? exactLoadDumpName(mesh.load)
                                     : name;
    saveSimulation(dumpName);
    if (usingDefaultName) {
      lastDefaultDumpLoadStep = mesh.loadSteps;
    }

    // Perhaps a bit strange, but this seems like a nice time to also
    // create/update the pvd file. (Sometimes it can be nice to have this
    // file before the simulation is done)
    gatherDataFiles();
    lastSaveTime = now;
    firstSaveDone = true;
  }
}

void Simulation::finishStep(bool reconnect) {
  // if (reconnect) {
  //   mesh.reconnect();
  // }
  // Calculate averages without updating full element data. Geometry still
  // matters because it refreshes F_P/F_E after reductions.
  mesh.ensureGeometry();
  mesh.updateForceStateAveragesAndPlasticEvents();
  mesh.updateCom();
  mesh.updateBoundingBox();
  updateEnergyHistory(true);
  computeParticipationFraction();
  computeM3ParticipationFraction();

  // Updates progress
  m_updateProgress();
  writeToFile();
  if (stepLogger != nullptr) {
    stepLogger(*this, stepLoggerContext);
  }

  // reset some counters
  const bool forceDump = forceDumpAfterStep;
  forceDumpAfterStep = false;
  mesh.resetCounters();
  m_writeDump(forceDump);
  loadingStepReplayCheckpoint.valid = false;

  assert(initialized); // Make sure the simulation has been initialized
}

void Simulation::finishSimulation() {
  gatherDataFiles();
  m_writeDump(true);
  csvFile.flush(); // Update the file on disk
  minCsvFile.flush();
  // timer.PrintAllRuntimes();
}

void Simulation::gatherDataFiles() {
  // This creates a pvd file that links all the utv files together.
  createCollection(getDataPath(simName, dataPath),
                   getOutputPath(simName, dataPath));
}

std::string Simulation::saveSimulation(std::string fileName_) {
  timer.Save();
  std::string fileName;

  // Generate filename dynamically if not provided
  if (fileName_.empty()) {
    fileName = defaultDumpName(mesh.load, startLoad, maxLoad) + ".xml.gz";
  } else {
    fileName = fileName_ + ".xml.gz";
  }

  std::string dumpPath =
      std::filesystem::path(getDumpPath(simName, dataPath)) / fileName;

  saveToFile(dumpPath, *this);

  if (!isQuiet()) {
    std::cout << "Dump saved to: " << dumpPath << std::endl;
  }
  return dumpPath;
}

void Simulation::saveCrashDump() {
  saveSimulation("crash_" + defaultDumpName(mesh.load, startLoad, maxLoad));
}

void Simulation::loadSimulation(Simulation &s, const std::string &dumpPath,
                                const std::string &conf, std::string outputPath,
                                std::optional<bool> forceReRun) {

  namespace fs = std::filesystem; // Alias for filesystem

  // Extract XML from .gz or .xml
  loadFromFile(dumpPath, s);
  s.initialized = false;

  // Load config file if provided
  if (!conf.empty()) {
    s.config = parseConfigFile(conf);
    if (s.rows != s.config.rows || s.cols != s.config.cols) {
      throw std::invalid_argument("Mesh size mismatch: Loaded (" +
                                  std::to_string(s.rows) + ", " +
                                  std::to_string(s.cols) + ") vs Config (" +
                                  std::to_string(s.config.rows) + ", " +
                                  std::to_string(s.config.cols) + ")");
    }
  }

  bool quiet = (isQuiet() || s.config.showProgress == -1);
  if (!quiet) {
    std::cout << "Loading simulation from " << dumpPath << std::endl;
  }

  // Update output path
  if (!outputPath.empty()) {
    s.dataPath = outputPath;
  }
  // Test of output path still exsists
  if (!fs::exists(s.dataPath)) {
    std::string oldPath = s.dataPath;
    if (!quiet) {
      std::cout << "Old output path is not found!" << std::endl;
    }

    s.dataPath = findOutputPath();
    if (!s.dataPath.empty()) {
      if (!quiet) {
        std::cout << "Old path:\n\t" << oldPath << "\nNew path:\n\t"
                  << s.dataPath << std::endl;
      }

      // Ask user for confirmation
      char choice;
      if (!quiet) {
        std::cout << "Do you want to continue with the new path? (y/n): "
                  << std::endl;
        std::cin >> choice;
      } else {
        choice = 'y';
      }

      if (choice != 'y' && choice != 'Y') {
        if (!quiet) {
          std::cout << "Aborting due to user rejection of new path."
                    << std::endl;
        }
        throw std::runtime_error("Old output path not found!");
      }
    }
  }

  quiet = (isQuiet() || s.config.showProgress == -1);
  if (!quiet) {
    std::cout << "Loading config..." << std::endl;
  }
  s.m_loadConfig(s.config);

  // Ensure mesh uses the resolved output path (e.g. after -o)
  s.mesh.setSimNameAndDataPath(s.simName, s.dataPath);

  // Use the parameter if provided, otherwise use the one from the config file
  bool effectiveForceReRun =
      forceReRun.has_value() ? *forceReRun : s.config.forceReRun;
  if (s.mesh.load + s.loadIncrement / 2 >= s.maxLoad && !effectiveForceReRun) {
    throw SimulationAlreadyComplete("Simulation already complete");
  }
  if (!quiet) {
    std::cout << "Saving config..." << std::endl;
  }
  saveConfigFile(s.config, s.dataPath);
  s.addDefaultCsvColumns();
  {
    const std::filesystem::path csvPath =
        std::filesystem::path(getOutputPath(s.simName, s.dataPath)) /
        (MACRODATANAME + std::string(".csv"));
    s.recoverCsvColumnsFromFile(csvPath.string());
  }
  s.csvFile = initCsvFile(s.simName, s.dataPath, s);

  if (!quiet) {
    std::cout << "Initializing..." << std::endl;
  }
  s.initialize();
  s.mesh.updateMesh();
  s.mesh.updateAveragesAndPlasticEvents();
  s.mesh.updateAngles();
  if (!quiet) {
    std::cout << "Simulation loaded!" << std::endl;
  }
}
