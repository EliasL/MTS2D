#include "simulation.h"
#include "protocols.h"
#include "Data/data_export.h"
#include "Data/logging.h"
#include "Data/param_parser.h"
#include "Mesh/node.h"
#include "randomUtils.h"
#include "settings.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <stdexcept>
#include <string>

namespace {

// These full-field VTUs are useful when debugging the reconnection algorithm,
// but are too costly for ordinary minimization replays. Keep this as a code
// flag so enabling minimization logging does not implicitly enable them.
constexpr bool kWriteReconnectionDebugVTUs = false;

} // namespace

Simulation::Simulation(Config config_, std::string _dataPath,
                       bool cleanDataPath) {
  dataPath = _dataPath;

  // This function initializes a lot of the variables using the config file
  m_loadConfig(config_);

  timer = Timer();

  mesh = Mesh(rows, cols, 1, config.QDSD, config.usesPeriodicBoundary(), config.meshDiagonal,
              config.energyFunction, config.bulkModulus);
  mesh.minimum_transport_distance = config.minimum_transport_distance;
  mesh.load = startLoad;
  mesh.setSimNameAndDataPath(simName, dataPath);
  addDefaultCsvColumns();

  if (simulationAlreadyComplete(simName, dataPath, endLoad()) &&
      !config.forceReRun) {
    throw SimulationAlreadyComplete("Simulation already complete");
  }
  if (cleanDataPath) {

    clearOutputFolder(simName, dataPath);
    saveConfigFile(config, dataPath);
    // Create and open file
  }
}

void Simulation::initialize() {
  if (initialized) {
    return;
  }
  timer.addKey("minimization");
  // Initialization should be done after nodes have been moved and fixed as
  // desired. The elements created by the function below are copies and do not
  // dynamically update. (the update function only updates the position,
  // energy and stress)

  // initialization of the csv file needs to be done after the correct loadStep
  // has been loaded
  // TODO: Keeping the csv file open during the entire simulation might be
  // unwise. We should perhaps save batches of lines, and open and close the
  // file as needed.
  addDefaultCsvColumns();
  csvFile = initCsvFile(simName, dataPath, *this);
  csvColumnsLocked = true;

  // We assume that the nodes already contain information about the mesh
  // structure, therefore, we recreate the elements
  // mesh.recreateElements();

  // we update the solvers
  initSolver();

  // Start simulation timer
  timer.Start();
  initialized = true;
}

void Simulation::initializeLoadedStateWithoutOutput() {
  initialized = false;
  m_loadConfig(config);
  timer.addKey("minimization");
  initSolver();
  mesh.updateMesh();
  mesh.updateAveragesAndPlasticEvents();
  mesh.updateAngles();
  initialized = true;
}

void Simulation::firstStep() {
  {
    // Keep the same physical initial state across solver and thread sweeps.
    // Restore before saving outputs, including on an exception.
    struct RestoreSettings {
      Config &target;
      Config saved;
      int threads;
      ~RestoreSettings() {
        target = std::move(saved);
        omp_set_num_threads(threads);
      }
    } restore{config, config, omp_get_max_threads()};
    config.minimizer = "LBFGS";
    config.preconditioner = "none";
    config.LBFGSNrCorrections = 10;
    config.LBFGSScale = 1.0;
    config.LBFGSEpsg = config.LBFGSEpsf = config.LBFGSEpsx = 1e-20;
    config.LBFGSMaxIterations = 0;
    config.epsR = 1e-6;
    config.nrThreads = 1;
    omp_set_num_threads(1);
    if (!initialized) {
      initialize();
    }

    if (mesh.loadSteps == 0 && startLoad != 0.0) {
      initializeProtocolDeformation(*this);
    }

    // Pretend to add load to correctly count load steps
    mesh.addLoad(0);
    // Set initial guess to the current position (starting load)
    applyLoadStepToGuess();
    // Add noise (from a seed) to trigger an avalanche
    addNoiseToGuess();
    // Minimizes the energy by moving the free nodes in the mesh
    minimize();
  }
  // Updates progress and writes to file
  finishStep();

  // Reset LBFGS report
  // Usually done before minimization, but will not be done if we are not using
  // LBFGS for the other loading steps
  LBFGS_report = alglib::minlbfgsreport();
  LBFGSRep = SimReport(LBFGS_report);
}

bool Simulation::keepLoading() {
  return mesh.load + loadIncrement <= endLoad() + 0.5 * loadIncrement;
}

bool Simulation::m_reconnect(Mesh::EdgeSet *lockedEdges) {
  // std::cout << "Reconnecting\n";
  if (reconnectionMethod == "edgeFlip") {
    const bool changed = mesh.reconnect(false, lockedEdges);
    if (changed) {
      stiffnessPreconditioner.reset();
    }
    return changed;
  } else if (reconnectionMethod == "delaunay") {
    mesh.reconnectDelaunay();
    stiffnessPreconditioner.reset();
    return true; // TODO check if Delaunay made any changes
  } else if (reconnectionMethod == "none") {
    // Do nothing
  } else {
    throw std::invalid_argument("Unknown reconnection method: " +
                                reconnectionMethod);
  }
  return true;
}

void Simulation::reconnectWithoutMinimization() {
  timer.Start("minimization");
  minCsvSubfolder.clear();
  minCsvHeaders.clear();
  nrReconnectingCycles = 0;
  reconnectStopReason = ReconnectStopReason::NotAttempted;
  rejectedReconnectEnergyDeltaValue = 0.0;

  mesh.updateMesh();
  mesh.updateForceStateAveragesAndPlasticEvents();
  updateEnergyHistory(false);

  if (reconnectionMethod == "none") {
    logMinimizationState();
    timer.Stop("minimization");
    return;
  }

  if (reconnectStepLogger != nullptr) {
    reconnectStepLogger(*this, ReconnectStepStage::BeforeReconnect,
                        reconnectStepLoggerContext);
  }

  const bool useEdgeLocking =
      config.reconnectEdgeLocking && reconnectionMethod == "edgeFlip";
  reconnectLockedEdges.clear();

  while (true) {
    nrReconnectingCycles++;
    const bool meshChanged =
        m_reconnect(useEdgeLocking ? &reconnectLockedEdges : nullptr);
    if (!meshChanged) {
      reconnectStopReason = ReconnectStopReason::NoTopologyChange;
      break;
    }
    mesh.updateMesh();
    if (reconnectionMethod == "delaunay") {
      reconnectStopReason = ReconnectStopReason::Completed;
      break;
    }
  }

  if (reconnectStepLogger != nullptr) {
    reconnectStepLogger(*this, ReconnectStepStage::AfterReconnect,
                        reconnectStepLoggerContext);
  }
  logMinimizationState();
  timer.Stop("minimization");
}

void Simulation::saveMeshCheckpoint() {
  meshCheckpoint = mesh;
  hasMeshCheckpoint = true;
}

void Simulation::restoreMeshCheckpoint() {
  if (!hasMeshCheckpoint) {
    throw std::runtime_error(
        "Simulation::restoreMeshCheckpoint: no mesh checkpoint saved.");
  }
  mesh = meshCheckpoint;
  stiffnessPreconditioner.reset();
}

const char *Simulation::reconnectStopReasonName() const {
  switch (reconnectStopReason) {
  case ReconnectStopReason::NoTopologyChange:
    return "no_topology_change";
  case ReconnectStopReason::NonImproving:
    return "non_improving";
  case ReconnectStopReason::Completed:
    return "completed";
  case ReconnectStopReason::NotAttempted:
  default:
    return "not_attempted";
  }
}

void Simulation::saveLoadingStepReplayCheckpoint(const DeformationStep &step) {
  LoadingStepReplayCheckpoint &checkpoint = loadingStepReplayCheckpoint;
  checkpoint.mesh = mesh;
  checkpoint.alglibDisplacements = alglibNodeDisplacements;
  checkpoint.lbfgsState = LBFGS_state;
  checkpoint.cgState = CG_state;
  checkpoint.fireDisplacements = FIRENodeDisplacements;
  checkpoint.energyHistory = energyHistory;
  checkpoint.fireRep = FIRERep;
  checkpoint.lbfgsRep = LBFGSRep;
  checkpoint.cgRep = CGRep;
  checkpoint.step = step;
  checkpoint.loadIncrement = loadIncrement;
  checkpoint.valid = true;
}

void Simulation::restoreLoadingStepReplayCheckpoint() {
  LoadingStepReplayCheckpoint &checkpoint = loadingStepReplayCheckpoint;
  if (!checkpoint.valid) {
    throw std::runtime_error("Simulation::restoreLoadingStepReplayCheckpoint: "
                             "no loading-step replay checkpoint saved.");
  }
  mesh = checkpoint.mesh;
  stiffnessPreconditioner.reset();
  alglibNodeDisplacements = checkpoint.alglibDisplacements;
  LBFGS_state = checkpoint.lbfgsState;
  CG_state = checkpoint.cgState;
  FIRENodeDisplacements = checkpoint.fireDisplacements;
  energyHistory = checkpoint.energyHistory;
  FIRERep = checkpoint.fireRep;
  LBFGSRep = checkpoint.lbfgsRep;
  CGRep = checkpoint.cgRep;
  loadIncrement = checkpoint.loadIncrement;
}

void Simulation::minimize(bool reconnect) {
  try {
    minimizeImpl(reconnect);
  } catch (...) {
    if (tryAdjustedMinimization(reconnect)) {
      return;
    }
    replayMinimizationAfterError(reconnect, std::current_exception());
  }
}

bool Simulation::tryAdjustedMinimization(bool reconnect) {
  if (debugReplayActive || !loadingStepReplayCheckpoint.valid ||
      config.minimizer != "LBFGS") {
    return false;
  }

  const int originalCorrections = config.LBFGSNrCorrections;
  const int problemDimension = alglibNodeDisplacements.length();
  if (originalCorrections < 1 || originalCorrections >= problemDimension) {
    return false;
  }

  const int adjustedCorrections = originalCorrections + 1;
  const DeformationStep step = loadingStepReplayCheckpoint.step;

  // The checkpoint is the state before the affine step. Keep this dump even
  // when the adjusted retry succeeds so the original failure is diagnosable.
  restoreLoadingStepReplayCheckpoint();
  saveCrashDump();

  config.LBFGSNrCorrections = adjustedCorrections;
  try {
    resetLBFGSState(adjustedCorrections);
    applyDeformation(step, false);
    minimizeImpl(reconnect);
  } catch (...) {
    const std::string adjustedMessage =
        DebugLog::exceptionMessage(std::current_exception());
    config.LBFGSNrCorrections = originalCorrections;
    restoreLoadingStepReplayCheckpoint();
    resetLBFGSState(originalCorrections);
    if (!isQuiet()) {
      std::cerr << "LBFGS M adjustment retry failed; continuing with debug "
                   "replay.\n"
                << adjustedMessage << '\n';
    }
    return false;
  }

  // m_minimize wrote a crash diagnostic row before throwing. The adjusted
  // retry has now produced the replacement state for the same load step, so
  // discard that failed row before finishStep writes the accepted one.
  discardFailedCurrentLoadStepCsvRow();

  // The accepted state came from M+1, but all subsequent steps must use the
  // configured value again. Preserve the successful report while rebuilding
  // the state for the next minimization.
  config.LBFGSNrCorrections = originalCorrections;
  resetLBFGSState(originalCorrections);
  if (!isQuiet()) {
    std::cout << "A crash was avoided by M adjustment at load=" << mesh.load
              << " step=" << mesh.loadSteps << std::endl;
  }
  return true;
}

void Simulation::replayMinimizationAfterError(
    bool reconnect, std::exception_ptr originalError) {
  if (debugReplayActive || !loadingStepReplayCheckpoint.valid) {
    std::rethrow_exception(originalError);
  }

  const std::string originalMessage = DebugLog::exceptionMessage(originalError);
  const DeformationStep step = loadingStepReplayCheckpoint.step;
  restoreLoadingStepReplayCheckpoint();
  saveCrashDump();
  debugReplayActive = true;

  try {
    applyDeformation(step);
    config.logDuringMinimization = true;
    minimizeImpl(reconnect);
  } catch (...) {
    debugReplayActive = false;
    throw std::runtime_error(DebugLog::formatDebugReplayFailure(
        originalMessage, DebugLog::exceptionMessage(std::current_exception())));
  }

  debugReplayActive = false;
  throw std::runtime_error(
      DebugLog::formatDebugReplayDidNotReproduce(originalMessage));
}

void Simulation::minimizeImpl(bool reconnect) {
  /*
  Pseudo code for minimization with reconnection:

  testEnergy = minimize()
  do:
    saveState()
    bestEnergy = testEnergy
    reconnect()
    testEnergy = minimize()
  while (testEnergy < bestEnergy)

  Revert to saved state.
  */
  nrReconnectingCycles = 0;
  reconnectStopReason = ReconnectStopReason::NotAttempted;
  rejectedReconnectEnergyDeltaValue = 0.0;

  if (mesh.freeNodeIds.size() == 0) {
    return;
  }
  if (reconnectionMethod == "none") {
    reconnect = false; // We don't need to run the minimization multiple times
                       // if we don't reconnect
  }
  timer.Start("minimization");
  minCsvSubfolder.clear();
  minCsvHeaders.clear();

  // If we log during minimization, we need a new file for each minimization
  if (config.logDuringMinimization) {
    minCsvSubfolder =
        std::string(DATAFOLDERPATH) + "/" + getMinDataSubFolder(mesh);
    minCsvFile = initCsvFile(simName, dataPath, *this, minCsvSubfolder);
    minCsvHeaders = readCsvHeaders(simName, dataPath, minCsvSubfolder);
  }

  // Save mesh before minimization for calculation of paticipation fraction
  mesh.captureDisplacementSnapshot(beforeMinimization);

  saveMeshCheckpoint();
  // First minimization before any reconnection cycles.
  m_minimize();
  mesh.updatePlasticEventCounts();
  saveMeshCheckpoint();
  if (reconnectStepLogger != nullptr) {
    reconnectStepLogger(*this, ReconnectStepStage::BeforeReconnect,
                        reconnectStepLoggerContext);
  }

  // std::cout << "Step: " << mesh.loadSteps << '\n';

  // We only reconnect if there is an energy drop since last load step
  // if (mesh.totalEnergy > energyHistory.prevLoadStepTotalEnergy) {
  //   reconnect = false;
  // }
  if (mesh.nr_elements_with_m3_changeInStep == 0) {
    reconnect = false;
  }

  if (!reconnect) { // If we don't reconnect, we are done after one
    // minimization
    if (reconnectStepLogger != nullptr) {
      reconnectStepLogger(*this, ReconnectStepStage::AfterReconnect,
                          reconnectStepLoggerContext);
    }
    logMinimizationState();
    // Save mesh after minimization for calculation of participation fraction
    mesh.captureDisplacementSnapshot(afterMinimization);
    timer.Stop("minimization");
    return;
  }

  double bestEnergy = mesh.totalEnergy;
  const bool useReconnectRevert = config.reconnectRevert;
  const bool useEdgeLocking =
      config.reconnectEdgeLocking && reconnectionMethod == "edgeFlip";
  reconnectLockedEdges.clear();
  bool meshChanged = false;
  while (true) {
    nrReconnectingCycles++;
    if (nrReconnectingCycles % 10 == 0) {
      std::cout << "Step: " << mesh.loadSteps
                << " Reconnections: " << nrReconnectingCycles << std::endl;
    }

    if (config.logDuringMinimization && kWriteReconnectionDebugVTUs) {
      mesh.writeToVtu("", true, VtuFieldLevel::All, "pre");
    }
    meshChanged = m_reconnect(useEdgeLocking ? &reconnectLockedEdges : nullptr);
    if (config.logDuringMinimization && kWriteReconnectionDebugVTUs) {
      mesh.writeToVtu("", true, VtuFieldLevel::All, "post");
    }
    if (!meshChanged) {
      reconnectStopReason = ReconnectStopReason::NoTopologyChange;
      break;
    }
    m_minimize();
    if (!useReconnectRevert) {
      saveMeshCheckpoint();
      continue;
    }
    if (mesh.totalEnergy < bestEnergy) {
      bestEnergy = mesh.totalEnergy;
      saveMeshCheckpoint();
      continue;
    }
    reconnectStopReason = ReconnectStopReason::NonImproving;
    rejectedReconnectEnergyDeltaValue = mesh.totalEnergy - bestEnergy;
    if (config.logDuringMinimization && kWriteReconnectionDebugVTUs) {
      mesh.writeToVtu("", true, VtuFieldLevel::All, "rejectedReconnect");
    }
    restoreMeshCheckpoint();
    break;
  }
  if (reconnectStepLogger != nullptr) {
    reconnectStepLogger(*this, ReconnectStepStage::AfterReconnect,
                        reconnectStepLoggerContext);
  }
  logMinimizationState();

  // Save mesh after minimization for calculation of paticipation fraction
  mesh.captureDisplacementSnapshot(afterMinimization);

  timer.Stop("minimization");
}

void Simulation::applyLoadStepToGuess(const Matrix2d &stepTransform) {
  // stepTransform is the incremental deformation since the last guess.
  // Composing it with the current displacements twice applies two increments.
  const int nr_x_values = alglibNodeDisplacements.length() / 2;
  const Matrix2d I = Matrix2d::Identity();
  const Matrix2d A = stepTransform - I; // often small

  for (size_t i = 0; i < mesh.freeNodeIds.size(); i++) {
    const Node *n = mesh[mesh.freeNodeIds[i]];

    // u = pos - init_pos (small)
    const Vector2d u = n->u();

    // stable form: (stepTransform - I)*init_pos + stepTransform*u
    // Explanation: u2 = stepTransform*x-x0, x = x0 + u1
    // Insert: u2 = stepTransform(x0+u1)-x0
    // Collect x0: u2 = (stepTransform-I)*x0 + stepTransform*u1
    const Vector2d next_displacement = A * n->ref_pos() + stepTransform * u;

    alglibNodeDisplacements[i] = next_displacement.x();
    alglibNodeDisplacements[i + nr_x_values] = next_displacement.y();
    FIRENodeDisplacements[i] = next_displacement.x();
    FIRENodeDisplacements[i + nr_x_values] = next_displacement.y();
  }

  // Not strictly necessary for minimization (the solver applies displacements
  // internally), but keeping mesh positions in sync avoids confusion. If you
  // rebuild a guess by calling setInitialGuess again, the mesh displacements
  // are used as the baseline.
  updateNodePositions(mesh, alglibNodeDisplacements);
}

// Core function to add (gausian) noise to a double array
void addNoiseToArray(double *data, size_t length, double noise) {
  double dataSum = 0;
  for (size_t i = 0; i < length; i++) {
    // Generate random noise from a normal distribution with mean 0 and
    // stddev noise
    data[i] += rng::sampleNormal(0.0, noise);
    dataSum += data[i];
  }
  // Print a noise hash to easily confirm two identical simulations
  // use as many decimal places as possible
  if (!isQuiet()) {
    std::cout << "Noise hash: " << std::scientific << std::setprecision(15)
              << dataSum << std::endl;
  }
}

// Overload for alglib::real_1d_array
void addNoise(alglib::real_1d_array &array, double noise) {
  addNoiseToArray(array.getcontent(), array.length(), noise);
}

// Overload for Eigen::VectorXd
void addNoise(Eigen::VectorXd &vector, double noise) {
  addNoiseToArray(vector.data(), vector.size(), noise);
}

// adds noise to both alglib and fire guess
void Simulation::addNoiseToGuess(double customNoise) {
  // Choose whether or not to use noise from argument or config file
  double effectiveNoise =
      customNoise == -1 ? config.initialGuessNoise : customNoise;

  if (config.minimizer == "FIRE") {
    addNoise(FIRENodeDisplacements, effectiveNoise);
  } else {
    addNoise(alglibNodeDisplacements, effectiveNoise);
  }
  // Noise only updates the displacement arrays. If you need the mesh node
  // positions to reflect the noisy guess (e.g., before saving a mesh copy),
  // call updateNodePositions(mesh,
  // alglibNodeDisplacements/FIRENodeDisplacements) explicitly.
}

void Simulation::m_loadConfig(Config config_) {
  config_.validatePreconditioner();
  // We save this for serialization
  config = config_;
  setQuiet(isQuiet() || config.showProgress == -1);
  // We fix the random seed to get reproducable results
  rng::setSeed(config.seed);
  // Set the the number of threads
  int maxThreads = omp_get_max_threads();
  int suggestedThreads = std::max(1, static_cast<int>(maxThreads * 0.75));

  if (config.nrThreads == 0) {
    config.nrThreads = suggestedThreads;
  } else if (config.nrThreads > maxThreads) {
    if (!isQuiet()) {
      std::cout << "Too many threads! Wanted " << config.nrThreads
                << ", but only " << maxThreads
                << " are available. Reducing nrThreads to: " << suggestedThreads
                << std::endl;
    }
    config.nrThreads = suggestedThreads;
  }

  omp_set_num_threads(config.nrThreads);
  // This dissables nested loops. We do not want any of these to be
  // happening.
  omp_set_max_active_levels(1);

  // Assign values from Config to Simulation members
  simName = config.name;
  mesh.simName = simName;
  mesh.energyFunction = config.energyFunction;
  mesh.minimum_transport_distance = config.minimum_transport_distance;
  mesh.bulkModulus = config.bulkModulus;
  mesh.updateLatticeBasis();
  rows = config.rows;
  cols = config.cols;

  startLoad = config.startLoad;
  loadIncrement = std::abs(config.loadIncrement);
  totalLoad = config.totalLoad;
  // This flag prevents mesh.reconnect() from being called in the
  // minimization function. If mehs.reconnect() is called somewhere else,
  // reconnection will still occur.
  reconnectionMethod = config.reconnectionMethod;
}

Matrix2d getShear(double load, double theta) {
  // perturb is currently unused. If it will be used, it should be
  // implemeted propperly.
  double perturb = 0;

  Matrix2d trans;
  trans(0, 0) = (1. - load * cos(theta + perturb) * sin(theta + perturb));
  trans(1, 1) = (1. + load * cos(theta - perturb) * sin(theta - perturb));
  trans(0, 1) = load * pow(cos(theta), 2.);
  trans(1, 0) = -load * pow(sin(theta - perturb), 2.);

  return trans;
}

void Simulation::addStepNoise() {
  applyLoadStepToGuess();
  addNoiseToGuess(config.stepNoise);
  updateNodePositions(mesh, alglibNodeDisplacements);
}
