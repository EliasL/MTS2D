#include "simulation.h"
#include "Data/data_export.h"
#include "Data/logging.h"
#include "Data/param_parser.h"
#include "Mesh/node.h"
#include "randomUtils.h"
#include "settings.h"
#include <Eigen/LU>
#include <FIRE.h>
#include <Param.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <optimization.h>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

Simulation::Simulation(Config config_, std::string _dataPath,
                       bool cleanDataPath) {
  dataPath = _dataPath;

  // This function initializes a lot of the variables using the config file
  m_loadConfig(config_);

  timer = Timer();

  mesh = Mesh(rows, cols, 1, config.QDSD, config.usingPBC, config.meshDiagonal,
              config.energyFunction, config.bulkModulus);
  mesh.load = startLoad;
  mesh.setSimNameAndDataPath(simName, dataPath);
  addDefaultCsvColumns();

  if (simulationAlreadyComplete(simName, dataPath, maxLoad) &&
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

  // We assume that the nodes already contain information about the mesh
  // structure, therefore, we recreate the elements
  // mesh.recreateElements();

  // we update the solvers
  initSolver();

  // Start simulation timer
  timer.Start();
  initialized = true;
}

void Simulation::initSolver() {

  // Alglib Initialization preparation
  int nrFreeNodes = mesh.freeNodeIds.size();
  alglibNodeDisplacements.setlength(2 * nrFreeNodes);
  // Set values to zero
  for (int i = 0; i < 2 * nrFreeNodes; i++) {
    alglibNodeDisplacements[i] = 0;
  }

  // Adjust nrCorrections based on the number of free nodes
  // alglib doesn't like that nr corrections is larger than nr of free nodes
  if (config.LBFGSNrCorrections > 2 * nrFreeNodes) {
    config.LBFGSNrCorrections = 2 * nrFreeNodes;
  }

  if (nrFreeNodes != 0) {
    // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgscreate
    // Initialize the state variable
    alglib::minlbfgscreate(config.LBFGSNrCorrections, alglibNodeDisplacements,
                           LBFGS_state);

    alglib::mincgcreate(alglibNodeDisplacements, CG_state);
  }

  // FIRE Initialization
  FIRENodeDisplacements = VectorXd::Zero(2 * nrFreeNodes);
  FIRE_param = FIREpp::FIREParam<double>(mesh.nrNodes, 2);
  config.updateParam(FIRE_param);

  // We update the datalink with the latest information about out simulation
  // data
  dataLink = DataLink(this);
  // TEMP
  // config.logDuringMinimization = true;
}

void Simulation::firstStep() {
  if (!initialized) {
    initialize();
  }

  if (mesh.loadSteps == 0 && startLoad != 0.0) {
    Matrix2d startLoadTransform = getShear(startLoad);
    mesh.applyTransformation(startLoadTransform);
  }

  // If it is the first step, we always minimize with the same algorithm to
  // ensure each seed has the same STABLE starting point
  // And we'll change the settings so that we get the same starting point for
  // all simulations
  double oldMaxForceAllowed = *dataLink.maxForceAllowed;
  double newMaxForceAllowed = 1e-6;
  *dataLink.maxForceAllowed = newMaxForceAllowed;

  std::string oldMinimizer = config.minimizer;
  std::string newMinimizer = "LBFGS";
  config.minimizer = newMinimizer;

  // Pretend to add load to correctly count load steps
  mesh.addLoad(0);
  // Set initial guess to the current position (starting load)
  applyLoadStepToGuess();
  // Add noise (from a seed) to trigger an avalanche
  addNoiseToGuess();
  // Minimizes the energy by moving the free nodes in the mesh
  minimize();
  // Updates progress and writes to file
  finishStep();

  // Reset settings
  *dataLink.maxForceAllowed = oldMaxForceAllowed; // Restore original value
  config.minimizer = oldMinimizer;

  // Reset LBFGS report
  // Usually done before minimization, but will not be done if we are not using
  // LBFGS for the other loading steps
  LBFGS_report = alglib::minlbfgsreport();
  LBFGSRep = SimReport(LBFGS_report);
}

bool Simulation::keepLoading() {
  double nextLoad = mesh.load + loadIncrement;
  double buffer = 0.5 * loadIncrement;

  if (loadIncrement > 0) {
    return nextLoad <= maxLoad + buffer;
  } else {
    return nextLoad >= startLoad + buffer;
  }
}

// Helper function
void Simulation::m_minimize(bool rough) {
  if (rough) {
    double roughMaxForceAllowed = 1e-2;
    dataLink.maxForceAllowed = &roughMaxForceAllowed;
  }

  auto fail = [&](std::string_view caughtType, const std::string &msg) -> void {
    const std::string error = DebugLog::formatMinimizationFailure(
        caughtType, msg, mesh, config.minimizer, rough);
    if (!isQuiet() && !debugReplayActive) {
      std::cerr << error << '\n';
    }
    writeToFile(true, "CrashAtLoad:" + std::to_string(mesh.load));
    throw std::runtime_error(error);
  };

  int nrIter = 0;
  try {
    if (config.minimizer == "FIRE") {
      m_minimizeWithFIRE();
      nrIter = FIRERep.nrIter;
    } else if (config.minimizer == "LBFGS") {
      m_minimizeWithLBFGS();
      nrIter = LBFGSRep.nrIter;
    } else if (config.minimizer == "CG") {
      m_minimizeWithCG();
      nrIter = CGRep.nrIter;
    } else {
      throw std::invalid_argument("Unknown minimizer: " + config.minimizer);
    }
    mesh.nrMinItterations += nrIter;
    mesh.nrMinItterationsSinceLastReconnect =
        std::max(mesh.nrMinItterationsSinceLastReconnect, nrIter);
  } catch (const alglib::ap_error &e) {
    fail("ALGLIB error", e.msg);
  } catch (const std::exception &e) {
    fail("Standard exception", e.what());
  } catch (...) {
    fail("Unknown exception", "");
  }
  // if (FIRERep.termType == -3) {
  //  writeToFile(true);
  //  throw std::runtime_error("Energy too high");
  //}
  if (LBFGSRep.termType == 1) {
    // mesh.writeToVtu("badStopStep" + std::to_string(mesh.loadSteps));
    // If maxForceAllowed is smaller than 1e-19, we assume it is zero and not
    // used
    if (mesh.maxForce > 1.5 * (*dataLink.maxForceAllowed) &&
        *dataLink.maxForceAllowed > 1e-19) {
      if (!isQuiet()) {
        std::cout << "Warning: LBFGS stopped with termType 1 but max force is "
                     "still high: "
                  << mesh.maxForce << '\n';
      }
    }
  }
  if (rough) {
    dataLink.maxForceAllowed = &config.epsR;
  }
}

bool Simulation::m_reconnect(Mesh::EdgeSet *lockedEdges) {
  // std::cout << "Reconnecting\n";
  if (reconnectionMethod == "edgeFlip") {
    return mesh.reconnect(false, lockedEdges);
  } else if (reconnectionMethod == "delaunay") {
    mesh.reconnectDelaunay();
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

  nrReconnectingCycles = 0;
  const bool useEdgeLocking =
      config.reconnectEdgeLocking && reconnectionMethod == "edgeFlip";
  reconnectLockedEdges.clear();

  while (true) {
    nrReconnectingCycles++;
    const bool meshChanged =
        m_reconnect(useEdgeLocking ? &reconnectLockedEdges : nullptr);
    if (!meshChanged) {
      break;
    }
    mesh.updateMesh();
    if (reconnectionMethod == "delaunay") {
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
}

void Simulation::saveLoadingStepReplayCheckpoint(const Matrix2d &affineStep) {
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
  checkpoint.affineStep = affineStep;
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
    replayMinimizationAfterError(reconnect, std::current_exception());
  }
}

void Simulation::replayMinimizationAfterError(
    bool reconnect, std::exception_ptr originalError) {
  if (debugReplayActive || !loadingStepReplayCheckpoint.valid) {
    std::rethrow_exception(originalError);
  }

  const std::string originalMessage = DebugLog::exceptionMessage(originalError);
  const Matrix2d affineStep = loadingStepReplayCheckpoint.affineStep;
  restoreLoadingStepReplayCheckpoint();
  debugReplayActive = true;

  try {
    applyAffineStep(affineStep);
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
    minCsvFile = initCsvFile(simName, dataPath, *this, minCsvSubfolder, false);
    minCsvHeaders = readCsvHeaders(simName, dataPath, minCsvSubfolder);
  }

  // Save mesh before minimization for calculation of paticipation fraction
  mesh.captureDisplacementSnapshot(beforeMinimization);

  saveMeshCheckpoint();
  // First minimization (If we reconnect, we also run a rough minimization)
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

  nrReconnectingCycles = 0;
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

    if (config.logDuringMinimization) {
      mesh.writeToVtu("", true, VtuFieldLevel::All, "pre");
    }
    meshChanged = m_reconnect(useEdgeLocking ? &reconnectLockedEdges : nullptr);
    if (config.logDuringMinimization) {
      mesh.writeToVtu("", true, VtuFieldLevel::All, "post");
    }
    if (!meshChanged) {
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
    if (config.writeDebugVTUs) {
      mesh.writeToVtu("", true, VtuFieldLevel::All, "deadEnd");
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

void Simulation::m_minimizeWithLBFGS() {
  timer.Start("LBFGSMinimization");
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgsrestartfrom
  // We reset and reuse the state instead of initializing it again
  // (The hessian is reset and not preserved)
  alglib::minlbfgsrestartfrom(LBFGS_state, alglibNodeDisplacements);

  // Set termination condition, ei. when is the solution good enough
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgssetcond
  alglib::minlbfgssetcond(LBFGS_state, config.LBFGSEpsg, config.LBFGSEpsf,
                          config.LBFGSEpsx, config.LBFGSMaxIterations);

  // Connect the chosen state to the minimization state
  minState = MinState(LBFGS_state);

  //  This is where the heavy calculations happen
  alglib::minlbfgsoptimize(LBFGS_state, alglibEnergyAndGradient,
                           iterationLogger, &dataLink);

  alglib::minlbfgsresults(LBFGS_state, alglibNodeDisplacements, LBFGS_report);
  LBFGSRep.nms = timer.Stop("LBFGSMinimization");
  LBFGSRep = SimReport(LBFGS_report);
}

void Simulation::m_minimizeWithCG() {
  timer.Start("CGMinimization");
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_mincgrestartfrom
  // We reset and reuse the state instead of initializing it again
  alglib::mincgrestartfrom(CG_state, alglibNodeDisplacements);

  // Set termination condition, ei. when is the solution good enough
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_mincgsetcond
  alglib::mincgsetcond(CG_state, config.CGEpsg, config.CGEpsf, config.CGEpsx,
                       config.CGMaxIterations);

  // TODO!
  // minState = MinState(CG_state);

  // This is where the heavy calculations happen
  // The null pointer can be replaced with a logging function
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_mincgoptimize
  alglib::mincgoptimize(CG_state, alglibEnergyAndGradient, iterationLogger,
                        &dataLink);

  alglib::mincgresults(CG_state, alglibNodeDisplacements, CG_report);
  CGRep.nms = timer.Stop("CGMinimization");
  CGRep = SimReport(CG_report);
}

template <typename ArrayType>
void updateMeshAndComputeForces(DataLink *dataLink, const ArrayType &disp,
                                double &energy, ArrayType &grad,
                                int nr_x_values) {
  Mesh *mesh = dataLink->mesh;

  // Update mesh position from the result of the
  // previous minimization
  updateNodePositions(*mesh, disp);

  // Calculate energy and forces
  mesh->updateMesh();
  assert(mesh->getUpdateState() == Mesh::UpdateState::Forces);

  // Total energy, only used for minimization
  energy = mesh->totalEnergy;

  // Update gradient in the minimization
  updateGradArray(mesh, grad, nr_x_values);
  double maxForce = mesh->maxForce;
  // int it = dataLink->LBFGS_state->c_ptr()->repiterationscount;

  // Determine if the minimization is done
  if (maxForce < *dataLink->maxForceAllowed) {
    // stop the minimization
    // Check if the current step is a valid step
    // TODO: Even if the max force is small right now, it might become large
    // when lbfgs selects the step with the lowest energy.
    alglib::minlbfgsrequesttermination(*dataLink->LBFGS_state);
    alglib::mincgrequesttermination(*dataLink->CG_state);
  }
  // We start to reconnect once we are 'close' to a solution
  // And only reconnect every 100 iterations
  // else if (mesh->nr_elements_with_m3_fix_changeInStep > 0) {

  //   mesh->reconnect(false, true);
  //   if (mesh->reconnectRequired) {
  //     // stop the minimization
  //     alglib::minlbfgsrequesttermination(*dataLink->LBFGS_state);
  //     alglib::mincgrequesttermination(*dataLink->CG_state);
  //   }
  // }
  // Since we don't use the x displacement argument in the iteration logger,
  // we just pass default parameter
  iterationLogger(alglib::real_1d_array(), energy, dataLink);
  // Update pastPlastic count so we can save files when it changes
  // mesh->resetPastPlasticCount(false);
  if (mesh->nrMinFunctionCalls == 0 &&
      !dataLink->s->config.logDuringMinimization) {
    // Ensure initial-guess stats are computed.
    // If we log during minimization, these will be updated in the iteration
    // logger instead
    mesh->updateForceStateAveragesAndPlasticEvents();
    dataLink->s->updateEnergyHistory(false);
  }
  mesh->nrMinFunctionCalls++;
  mesh->nrMinFunctionCallsSinceLastReconnect++;
}

void alglibEnergyAndGradient(const alglib::real_1d_array &disp, double &energy,
                             alglib::real_1d_array &grad, void *dataLinkPtr) {
  DataLink *dataLink = reinterpret_cast<DataLink *>(dataLinkPtr);
  updateMeshAndComputeForces(dataLink, disp, energy, grad, grad.length() / 2);
}

double FIREEnergyAndGradient(Eigen::VectorXd &disp, Eigen::VectorXd &grad,
                             void *dataLinkPtr) {
  double energy;
  DataLink *dataLink = reinterpret_cast<DataLink *>(dataLinkPtr);
  updateMeshAndComputeForces(dataLink, disp, energy, grad, grad.size() / 2);
  return energy;
}

void Simulation::m_minimizeWithFIRE() {
  timer.Start("FIREMinimization");
  FIREpp::FIRESolver<double> s = FIREpp::FIRESolver<double>(FIRE_param);
  double energy;
  FIRERep.nrIter = s.minimize(FIREEnergyAndGradient, FIRENodeDisplacements,
                              energy, &dataLink, FIRERep.termType);
  FIRERep.nms = timer.Stop("FIREMinimization");
  FIRERep.nfev = mesh.nrMinFunctionCalls;

  if (FIRERep.termType == -3) {
    writeToFile(true);
  }
  // We first do FIRE, but then use the result as an initial guess for LBFGS
  // setInitialGuess();
  // m_minimizeWithLBFGS();
}

// Overload for alglib::real_1d_array
void updateNodePositions(Mesh &mesh, const alglib::real_1d_array &disp) {
  mesh.updateNodePositions(disp.getcontent(), disp.length());
}

// Overload for Eigen::VectorXd
void updateNodePositions(Mesh &mesh, const Eigen::VectorXd &disp) {
  mesh.updateNodePositions(disp.data(), disp.size());
}

// Updates gradient array from current node forces.
template <typename ArrayType>
void updateGradArray(Mesh *mesh, ArrayType &grad, int nr_x_values) {
  for (int i = 0; i < nr_x_values; i++) {
    NodeId n_id = mesh->freeNodeIds[i];
    const auto &node = (*mesh)[n_id];

    // Store force values in the array
    double fx = node->f[0];
    double fy = node->f[1];

    // Gadient is the oposite sign of the force
    grad[i] = -fx;
    grad[nr_x_values + i] = -fy;
  }
}

void Simulation::applyPreviousMinimizationCorrectionToGuess() {
  const size_t nodeCount = static_cast<size_t>(mesh.rows * mesh.cols);
  const size_t expectedSize = 2 * nodeCount;
  if (beforeMinimization.displacements.size() != expectedSize ||
      afterMinimization.displacements.size() != expectedSize) {
    throw std::runtime_error(
        "Simulation::applyPreviousMinimizationCorrectionToGuess: "
        "displacement snapshot sizes do not match mesh size.");
  }

  const int nr_x_values = alglibNodeDisplacements.length() / 2;
  if (nr_x_values != static_cast<int>(mesh.freeNodeIds.size()) ||
      FIRENodeDisplacements.size() != 2 * nr_x_values) {
    throw std::runtime_error(
        "Simulation::applyPreviousMinimizationCorrectionToGuess: solver arrays "
        "do not match free-node count.");
  }

  const double *beforeX = beforeMinimization.displacements.data();
  const double *beforeY = beforeX + nodeCount;
  const double *afterX = afterMinimization.displacements.data();
  const double *afterY = afterX + nodeCount;

  for (int i = 0; i < nr_x_values; i++) {
    const NodeId nodeId = mesh.freeNodeIds[i];
    const size_t nodeIndex = static_cast<size_t>(nodeId.i);
    const Node *n = mesh[nodeId];
    const Vector2d correction(afterX[nodeIndex] - beforeX[nodeIndex],
                              afterY[nodeIndex] - beforeY[nodeIndex]);
    const Vector2d nextDisplacement = n->u() + correction;
    alglibNodeDisplacements[i] = nextDisplacement.x();
    alglibNodeDisplacements[i + nr_x_values] = nextDisplacement.y();
    FIRENodeDisplacements[i] = nextDisplacement.x();
    FIRENodeDisplacements[i + nr_x_values] = nextDisplacement.y();
  }
  updateNodePositions(mesh, alglibNodeDisplacements);
}

void Simulation::applyLoadStepToGuess(const Matrix2d &T) {
  // NOTE: T is expected to be the *incremental* deformation since the last
  // guess. This function composes T with the current displacements, so calling
  // it twice with the same T applies the deformation twice.
  const int nr_x_values = alglibNodeDisplacements.length() / 2;
  const Matrix2d I = Matrix2d::Identity();
  const Matrix2d A = T - I; // often small

  for (size_t i = 0; i < mesh.freeNodeIds.size(); i++) {
    const Node *n = mesh[mesh.freeNodeIds[i]];

    // u = pos - init_pos (small)
    const Vector2d u = n->u();

    // stable form: (T - I)*init_pos + T*u
    // Explination: u2 = Tx-x0, x = x0 + u1
    // Insert: u2 = T(x0+u1)-x0
    // Collect x0: u2 = (T-I)x0+Tu1
    const Vector2d next_displacement = A * n->ref_pos() + T * u;

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

void Simulation::applyAffineStep(const Matrix2d &T) {
  if (!debugReplayActive) {
    saveLoadingStepReplayCheckpoint(T);
  }
  mesh.addLoad(loadIncrement);
  if (mesh.usingPBC) {
    mesh.applyTransformationToSystemDeformation(T);
  } else {
    mesh.applyTransformationToFixedNodes(T);
  }
  applyLoadStepToGuess(T);

  // I tried to enhance the inital guess. It works well for low deformations,
  // but is worse further out into the simulation.
  constexpr bool enablePreviousMinimizationCorrectionPredictor = false;
  if (!enablePreviousMinimizationCorrectionPredictor) {
    return;
  }

  const bool hasPreviousCorrection =
      !beforeMinimization.displacements.empty() &&
      !afterMinimization.displacements.empty();
  if (beforeMinimization.displacements.empty() !=
      afterMinimization.displacements.empty()) {
    throw std::runtime_error(
        "Simulation::applyAffineStep: exactly one minimization displacement "
        "snapshot is available.");
  }
  if (energyHistory.loadStepTotalEnergyChange > 0.0 && hasPreviousCorrection) {
    applyPreviousMinimizationCorrectionToGuess();
  }
}

// Core function to add (gausian) noise to a double array
void addNoiseToArray(double *data, size_t length, double noise) {
  double dataSum = 0;
  for (size_t i = 0; i < length; i++) {
    // Generate random noise from a normal distribution with mean 0 and
    // stddev noise
    data[i] += sampleNormal(0.0, noise);
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
  // We save this for serialization
  config = config_;
  setQuiet(isQuiet() || config.showProgress == -1);
  // We fix the random seed to get reproducable results
  setSeed(config.seed);
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
  mesh.bulkModulus = config.bulkModulus;
  mesh.updateLatticeBasis();
  rows = config.rows;
  cols = config.cols;

  startLoad = config.startLoad;
  loadIncrement = config.loadIncrement;
  maxLoad = config.maxLoad;
  // This flag prevents mesh.reconnect() from being called in the
  // minimization function. If mehs.reconnect() is called somewhere else,
  // reconnection will still occur.
  reconnectionMethod = config.reconnectionMethod;
}

// This can be used in the LBFGS optimization
void iterationLogger(const alglib::real_1d_array &x, double energy,
                     void *dataLinkPtr) {
  (void)x; // Explicitly casting x to void to silence unused parameter warning

  DataLink *dataLink = reinterpret_cast<DataLink *>(dataLinkPtr);

  // dataPath looks like: "/Volumes/data/MTS2D_output/"

  Mesh *mesh = dataLink->mesh;

  // auto state = dataLink->minState;

  // Increment the iteration count
  // TODO Does not work for CG or FIRE, only for LBFGS
  int it = dataLink->LBFGS_state->c_ptr()->repiterationscount;
  int nrFc = mesh->nrMinFunctionCalls;
  mesh->nrMinItterationsSinceLastReconnect =
      std::max(mesh->nrMinItterationsSinceLastReconnect, it);

  // Check if iteration count is a multiple of 5000
  if (nrFc % 5000 == 0 && nrFc > 0) {
    if (!isQuiet()) {
      std::cout << "Warning (step " << mesh->loadSteps << "): " << it
                << " iterations and " << nrFc << " function calls."
                << std::endl;
    }
  }

  // Save when energy shifts by >10% vs last saved, otherwise every 1000
  // function calls.
  int saveEvery = 1000;
  static int lastSavedFc = -1;
  static bool hasLastSavedEnergy = false;
  static double lastSavedEnergy = 0.0;

  if (nrFc == 0) {
    lastSavedFc = -1;
    hasLastSavedEnergy = false;
    lastSavedEnergy = 0.0;
  }

  if (dataLink->s->config.logDuringMinimization) {
    dataLink->s->timer.Start("write");
    // Write to the CSV file
    dataLink->s->logMinimizationState();
    bool energyJump = false;
    if (!hasLastSavedEnergy) {
      energyJump = true;
    } else {
      double denom = std::max(std::abs(lastSavedEnergy), 1e-12);
      energyJump = std::abs(energy - lastSavedEnergy) / denom > 0.10;
    }

    bool stepSave = (lastSavedFc < 0) || (abs(nrFc - lastSavedFc) >= saveEvery);

    if (energyJump || stepSave) {
      // Write mesh to file
      mesh->writeToVtu("", true, VtuFieldLevel::Minimal);
      lastSavedFc = nrFc;
      lastSavedEnergy = energy;
      hasLastSavedEnergy = true;
    }
    dataLink->s->timer.Stop("write");
  }
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

DataLink::DataLink(Simulation *simulation) {
  s = simulation;
  mesh = &s->mesh;
  stopSignal = nullptr; // This is set right before minimization
  minState = &s->minState;
  LBFGS_state = &s->LBFGS_state;
  CG_state = &s->CG_state;
  maxForceAllowed = &s->config.epsR;
}
