#include "simulation.h"
#include "Data/cereal_help.h"
#include "Data/data_export.h"
#include "Data/logging.h"
#include "Data/param_parser.h"
#include "Eigen/src/Core/Matrix.h"
#include "Mesh/node.h"
#include "randomUtils.h"
#include "settings.h"
#include <FIRE.h>
#include <Param.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <omp.h>
#include <optimization.h>
#include <ostream>
#include <stdexcept>
#include <string>

Simulation::Simulation(Config config_, std::string _dataPath,
                       bool cleanDataPath) {
  dataPath = _dataPath;

  // This function initializes a lot of the variables using the config file
  m_loadConfig(config_);

  timer = Timer();

  mesh = Mesh(rows, cols, 1, config.QDSD, config.usingPBC, config.meshDiagonal);
  mesh.load = startLoad;
  mesh.setSimNameAndDataPath(simName, dataPath);

  if (simulationAlreadyComplete(simName, dataPath, maxLoad) &&
      !config.forceReRun) {
    std::cout << "Simulation already complete\n";
    exit(EXIT_SUCCESS);
  }
  if (cleanDataPath) {

    clearOutputFolder(simName, dataPath);
    saveConfigFile(config, dataPath);
    // Create and open file
  }
}

void Simulation::initialize() {
  // Initialization should be done after nodes have been moved and fixed as
  // desired. The elements created by the function below are copies and do not
  // dynamically update. (the update function only updates the position,
  // energy and stress)

  // initialization of the csv file needs to be done after the correct loadStep
  // has been loaded
  csvFile = initCsvFile(simName, dataPath, *this);

  // We assume that the nodes already contain information about the mesh
  // structure, therefore, we recreate the elements
  // mesh.recreateElements();

  // we update the solvers
  initSolver();

  // Start simulation timer
  timer.Start();
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

  // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgscreate
  // Initialize the state variable
  alglib::minlbfgscreate(config.LBFGSNrCorrections, alglibNodeDisplacements,
                         LBFGS_state);

  alglib::mincgcreate(alglibNodeDisplacements, CG_state);

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
  setInitialGuess();
  // Add noise (from a seed) to trigger the already exsisting instability
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

  if (loadIncrement > 0) {
    return nextLoad < maxLoad;
  } else {
    return nextLoad > startLoad;
  }
}

void Simulation::minimize(bool reconnect) {
  timer.Start("minimization");

  if (config.logDuringMinimization) {
    minCsvFile = initCsvFile(simName, dataPath, *this,
                             std::string(DATAFOLDERPATH) + "/" +
                                 getMinDataSubFolder(mesh));
  }
  bool repeatMinimization = false;
  std::string error_message;
  int maxReconnecting = 20000;
  int currentReconnecting = 0;
  do {
    try {
      if (config.minimizer == "FIRE") {
        m_minimizeWithFIRE();
      } else if (config.minimizer == "LBFGS") {
        m_minimizeWithLBFGS();
      } else if (config.minimizer == "CG") {
        m_minimizeWithCG();
      } else {
        throw std::invalid_argument("Unknown minimizer: " + config.minimizer);
      }
    } catch (const alglib::ap_error &e) {
      error_message = "ALGLIB error: " + std::string(e.msg);
    } catch (const std::exception &e) {
      error_message = "Standard exception: " + std::string(e.what());
    } catch (...) {
      error_message = "Unknown exception caught!";
    }

    if (!error_message.empty()) {
      std::cerr << error_message << std::endl;
      writeToFile(true, "CrashAtLoad:" + std::to_string(mesh.load));
      std::exit(EXIT_FAILURE);
    }

    if (FIRERep.termType == -3) {
      // writeToFile(true);
      // throw std::runtime_error("Energy too high");
    }
    bool badStop = false;
    if (LBFGSRep.termType == 1 && !badStop) {
      badStop = true;
      // mesh.writeToVtu("badStopStep" + std::to_string(mesh.loadSteps));
      std::cout << "LBFGS stopped with termType 1, "
                << "max force: " << mesh.maxForce << '\n';
    }

    currentReconnecting++;
    if (currentReconnecting % 20 == 0) {
      std::cout << "Step: " << mesh.loadSteps
                << " Reconnections: " << currentReconnecting << "\n";
    }
    if (currentReconnecting > maxReconnecting) {
      std::cout << "Step: " << mesh.loadSteps << ". Too many reconnections!"
                << '\n';
      mesh.reconnectRequired = false;
    }
    mesh.nrMinItterations += dataLink.LBFGS_state->c_ptr()->repiterationscount;
    if (mesh.reconnectRequired && reconnect && reconnectingEnabled) {
      // We need to reconnect, so we set the flag to true
      repeatMinimization = true;
      mesh.reconnect(true);
    } else {
      repeatMinimization = false;
    }
  } while (repeatMinimization);

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
  //  The null pointer can be replaced with a logging function
  //  https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgsoptimize
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
  // TODO give the mesh the property largest force and update this one during
  // energy calculation in a thread safe way... try to perhaps disable checking
  // the stopping with other criteria? To be a bit faster?

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

  // Total energy, only used for minimization
  energy = mesh->totalEnergy;

  // Update gradient in the minimization
  double maxForce = updateGradArray(mesh, grad, nr_x_values);
  mesh->maxForce = maxForce;
  int it = dataLink->LBFGS_state->c_ptr()->repiterationscount;

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
  // And only reconnect every 10 iterations
  else if (maxForce / 1000 < *dataLink->maxForceAllowed && it % 100 == 0) {

    mesh->reconnect(true, true);
    if (mesh->reconnectRequired) {
      // stop the minimization
      alglib::minlbfgsrequesttermination(*dataLink->LBFGS_state);
      alglib::mincgrequesttermination(*dataLink->CG_state);
    }
  }
  // Since we don't use the x displacement argument in the iteration logger,
  // we just pass default parameter
  iterationLogger(alglib::real_1d_array(), energy, dataLink);
  mesh->nrMinFunctionCalls++;
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

// Returns max force component
template <typename ArrayType>
double updateGradArray(Mesh *mesh, ArrayType &grad, int nr_x_values) {
  double maxForce = 0;

  // #pragma omp parallel for reduction(max : maxForce)
  for (int i = 0; i < nr_x_values; i++) {
    NodeId n_id = mesh->freeNodeIds[i];
    const auto &node = (*mesh)[n_id];

    // Store force values in the array
    double fx = node->f[0];
    double fy = node->f[1];

    // Gadient is the oposite sign of the force
    grad[i] = -fx;
    grad[nr_x_values + i] = -fy;

    // Update max force component
    maxForce = std::max({maxForce, std::abs(fx), std::abs(fy)});
  }

  return maxForce;
}

// This function modifies the nodeDisplacements variable used in the solver
// NB Note that they are displacements, not positions
void Simulation::setInitialGuess(Matrix2d guessTransformation) {
  // Our initial guess will be that all particles have shifted by the same
  // transformation as the border.
  // The disp is structed like this: [x1,x2,x3,x4,y1,y2,y3,y4], so we
  // need to know where the "x" end and where the "y" begin.
  int nr_x_values = alglibNodeDisplacements.length() / 2; // Shifts to y section

  const Node *n; // free node

  // We loop over all the nodes that are not on the border, ie. the innside
  // nodes.
  for (size_t i = 0; i < mesh.freeNodeIds.size(); i++) {
    n = mesh[mesh.freeNodeIds[i]];
    Vector2d next_displacement = guessTransformation * n->pos() - n->init_pos();
    alglibNodeDisplacements[i] = next_displacement[0];
    alglibNodeDisplacements[i + nr_x_values] = next_displacement[1];

    FIRENodeDisplacements[i] = next_displacement[0];
    FIRENodeDisplacements[i + nr_x_values] = next_displacement[1];
  }
}

// Core function to add (gausian) noise to a double array
void addNoiseToArray(double *data, size_t length, double noise) {
  for (size_t i = 0; i < length; i++) {
    // Generate random noise from a normal distribution with mean 0 and stddev
    // noise
    data[i] += sampleNormal(0.0, noise);
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
  addNoise(alglibNodeDisplacements, effectiveNoise);

  addNoise(FIRENodeDisplacements, effectiveNoise);
}

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

  // Use static variables to track the last progress and the last update time
  static int oldProgress = -1;
  static int firstProgress = -1;
  static int lastDump = -1;
  static auto lastUpdateTime = steady_clock::now();

  auto now = steady_clock::now();
  auto timeSinceLastUpdate =
      duration_cast<seconds>(now - lastUpdateTime).count();

  bool shouldUpdate = config.showProgress == 1 &&
                      (oldProgress != intProgress || timeSinceLastUpdate >= 20);

  if (shouldUpdate) {
    oldProgress = intProgress;
    firstProgress = (firstProgress == -1) ? intProgress : firstProgress;
    lastUpdateTime = now;

    std::cout << consoleProgressMessage << std::endl;

    // Check if we should make a dump
    if (intProgress % 5 == 0 && intProgress != lastDump &&
        firstProgress != intProgress) {
      lastDump = intProgress;
      m_writeDump(true);
    }
  }
}

void Simulation::writeToFile(bool forceWrite, std::string fileName) {
  timer.Start("write");
  // We write to the cvs file every time this function is called
  writeToCsv(csvFile, (*this));
  // These are writing date much less often
  m_writeMesh(forceWrite);
  m_writeDump(forceWrite, fileName);
  if (config.logDuringMinimization) {
    // If we are logging minimization, we want to create a collection after each
    // step

    createCollection(
        getDataPath(simName, dataPath) + "/" + getMinDataSubFolder(mesh),
        getDataPath(simName, dataPath) + "/" + getMinDataSubFolder(mesh),
        ".*_minStep=[0-9]+.([0-9]+)_.*");
  }
  timer.Stop("write");
}

void Simulation::m_writeMesh(bool forceWrite) {
  // Only if there are lots of plastic events will we want to save the data.
  // If we save every frame, it requires too much storage.
  // (A 100x100 system loaded from 0.15 to 1 with steps of 1e-5 would take up
  // 180GB)
  // At the same time, if there are few large avalanvhes, we might go long
  // without saving data. In order to get a good framerate for an animation, we
  // want to ensure that not too much happens between frames. The following
  // enures that we at least have 200 frames of states over the course of
  // loading, but also don't miss any big events
  static double lastLoadWritten = 0;
  if ((mesh.nrPlasticChanges >
       mesh.nrElements *
           config.plasticityEventThreshold) || // Lots of plastic change
      (-mesh.delAvgEnergy > config.energyDropThreshold &&
       mesh.nrPlasticChanges > 0) || // Large energy drop (with at least one pc)
      (abs(mesh.load - lastLoadWritten) > 0.005) || // Absolute change
      (abs(mesh.load - lastLoadWritten) / (maxLoad - startLoad) >
       0.005) ||  // Relative change
      forceWrite) // Force write
  {
    mesh.writeToVtu();
    lastLoadWritten = mesh.load;
  }
}

void Simulation::m_writeDump(bool forceWrite, std::string name) {
  // When do we create save states?
  // I'm thinking I want to do one halfway no matter how short the simulation
  // is, but then outside of that, i'm thinking once per hour is okay.
  auto now = std::chrono::steady_clock::now();
  static auto lastSaveTime =
      now; // Since this is a static variable, this line is only run once
  static bool firstSaveDone = false;
  auto elapsedSinceLastSave =
      std::chrono::duration_cast<std::chrono::seconds>(now - lastSaveTime);
  double midPointLoad = (startLoad + maxLoad) / 2;

  // Save every one hours
  static const std::chrono::hours saveFrequency(1);

  if ((mesh.load >= midPointLoad &&
       !firstSaveDone) || // Check for first save at midpoint
      (elapsedSinceLastSave >= saveFrequency) || // Hourly save
      forceWrite)                                // Check if forced
  {
    saveSimulation(name);

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
  //   Calculate averages
  mesh.calculateAverages();
  // Updates progress
  m_updateProgress();
  writeToFile();

  // reset some counters
  mesh.resetCounters();
}

void Simulation::m_loadConfig(Config config_) {
  // We save this for serialization
  config = config_;
  // We fix the random seed to get reproducable results
  setSeed(config.seed);
  // Set the the number of threads
  int maxThreads = omp_get_max_threads();
  int suggestedThreads = std::max(1, static_cast<int>(maxThreads * 0.75));

  if (config.nrThreads == 0) {
    config.nrThreads = suggestedThreads;
  } else if (config.nrThreads > maxThreads) {
    std::cout << "Too many threads! Wanted " << config.nrThreads
              << ", but only " << maxThreads
              << " are available. Reducing nrThreads to: " << suggestedThreads
              << std::endl;
    config.nrThreads = suggestedThreads;
  }

  omp_set_num_threads(config.nrThreads);
  // This dissables nested loops. We do not want any of these to be happening.
  omp_set_max_active_levels(1);

  // Assign values from Config to Simulation members
  simName = config.name;
  mesh.simName = simName;
  rows = config.rows;
  cols = config.cols;

  startLoad = config.startLoad;
  loadIncrement = config.loadIncrement;
  maxLoad = config.maxLoad;
  reconnectingEnabled = config.reconnectingEnabled;
}

void Simulation::finishSimulation() {
  gatherDataFiles();
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
    double range = std::abs(maxLoad - startLoad);
    int precision =
        std::max(1, static_cast<int>(std::ceil(-std::log10(range / 10.0))));
    double roundedLoad = std::round(mesh.load * std::pow(10.0, precision)) /
                         std::pow(10.0, precision);

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << roundedLoad;
    fileName = "dump_l" + oss.str() + ".xml.gz"; // Use .zip extension
  } else {
    fileName = fileName_ + ".xml.gz";
  }

  std::string dumpPath =
      std::filesystem::path(getDumpPath(simName, dataPath)) / fileName;

  saveToFile(dumpPath, *this);

  std::cout << "Dump saved to: " << dumpPath << std::endl;
  return dumpPath;
}

void Simulation::loadSimulation(Simulation &s, const std::string &dumpPath,
                                const std::string &conf, std::string outputPath,
                                const bool forceReRun) {
  std::cout << "Loading simulation from " << dumpPath << std::endl;

  // Extract XML from .gz or .xml
  loadFromFile(dumpPath, s);

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

  // Update output path
  if (!outputPath.empty()) {
    s.dataPath = outputPath;
  }
  std::cout << "Loading config..." << std::endl;
  s.m_loadConfig(s.config);

  if (s.mesh.load >= s.maxLoad && !forceReRun) {
    std::cout << "Simulation already complete\n";
    exit(EXIT_SUCCESS);
  }
  std::cout << "Saving config..." << std::endl;
  saveConfigFile(s.config, s.dataPath);
  s.csvFile = initCsvFile(s.simName, s.dataPath, s);

  std::cout << "Initializing..." << std::endl;
  s.initialize();

  std::cout << "Done!" << std::endl;
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

  // Check if iteration count is a multiple of 5000
  if (nrFc % 5000 == 0 && nrFc > 0) {
    std::cout << "Warning (step " << mesh->loadSteps << "): " << it
              << " iterations and " << nrFc << " function calls." << std::endl;
  }

  // We only save every 100 steps, unless there is a plastic change
  int saveEvery = 100;
  static int lastSavedFc = 0;

  if (dataLink->s->config.logDuringMinimization) {
    dataLink->s->timer.Start("write");

    mesh->calculateAverages(false);
    writeToCsv(dataLink->s->minCsvFile, *dataLink->s);
    if (mesh->nrPlasticChangesInStep > 0 ||
        abs(nrFc - lastSavedFc) >= saveEvery || nrFc == 0) {
      // Write to the CSV file
      // Write mesh to file
      mesh->writeToVtu("", true);
      dataLink->s->timer.Stop("write");

      lastSavedFc = nrFc;
      // Update pastPlastic count so we can save files when it changes
      mesh->resetPastPlasticCount(false);
    }
  }
}

Matrix2d getShear(double load, double theta) {
  // perturb is currently unused. If it will be used, it should be implemeted
  // propperly.
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