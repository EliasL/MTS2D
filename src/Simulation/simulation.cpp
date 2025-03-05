#include "simulation.h"
#include "Data/cereal_help.h"
#include "Data/data_export.h"
#include "Data/logging.h"
#include "Data/param_parser.h"
#include "Eigen/src/Core/Matrix.h"
#include "Mesh/node.h"
#include "randomUtils.h"
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
  mesh.setSimNameAndDataPath(name, dataPath);

  if (simulationAlreadyComplete(name, dataPath, maxLoad) &&
      !config.forceReRun) {
    std::cout << "Simulation already complete\n";
    exit(EXIT_SUCCESS);
  }
  if (cleanDataPath) {

    clearOutputFolder(name, dataPath);
    saveConfigFile(config);
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
  csvFile = initCsvFile(name, dataPath, *this);

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

  if (config.logDuringMinimization) {
    // I'm not sure if this is the intended way to enable logging, but
    // it works
    LBFGS_state.c_ptr()->xrep = true;
    CG_state.c_ptr()->xrep = true;
    // TODO FIRE
  }
}

void Simulation::firstStep() {

  // If it is the first step, we always minimize with the same algorithm to
  // ensure each seed has the same STABLE starting point
  // And we'll change the settings so that we get the same starting point for
  // all simulations
  double oldValue = *dataLink.maxForceAllowed;
  double newValue = 1e-6;
  *dataLink.maxForceAllowed = newValue;

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
  *dataLink.maxForceAllowed = oldValue; // Restore original value
  config.minimizer = oldMinimizer;

  // Reset LBFGS report
  // Usually done before minimization, but will not be done if we are not using
  // LBFGS for the other loading steps
  LBFGS_report = alglib::minlbfgsreport();
  LBFGSRep = SimReport(LBFGS_report);
}

bool Simulation::keepLoading() {
  // Determine direction based on the sign of loadIncrement
  int direction = loadIncrement > 0 ? 1 : -1;
  return mesh.load * direction < maxLoad * direction;
}

void Simulation::minimize() {
  timer.Start("minimization");

  // We need to reset some counters
  mesh.resetCounters();

  if (config.logDuringMinimization) {
    minCsvFile =
        initCsvFile(name, dataPath, *this,
                    "minimizationData/step" + std::to_string(mesh.loadSteps));
  }

  if (config.minimizer == "FIRE") {
    m_minimizeWithFIRE();
  } else if (config.minimizer == "LBFGS") {
    m_minimizeWithLBFGS();
  } else if (config.minimizer == "CG") {
    m_minimizeWithCG();
  } else {
    std::cout << config.minimizer << std::endl;
    throw std::invalid_argument("Unknown minimizer");
  }
  if (FIRERep.termType == -3) {
    // writeToFile(true);
    // throw std::runtime_error("Energy too high");
  }
  timer.Stop("minimization");
}

void Simulation::m_minimizeWithLBFGS() {
  timer.Start("LBFGSMinimization");
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgsrestartfrom
  // We reset and reuse the state instead of initializing it again
  alglib::minlbfgsrestartfrom(LBFGS_state, alglibNodeDisplacements);

  // Set termination condition, ei. when is the solution good enough
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgssetcond
  alglib::minlbfgssetcond(LBFGS_state, config.LBFGSEpsg, config.LBFGSEpsf,
                          config.LBFGSEpsx, config.LBFGSMaxIterations);

  // Connect stopSignal in mesh to state
  dataLink.stopSignal = &LBFGS_state.c_ptr()->userterminationneeded;
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

  // Connect stopSignal in mesh to state
  dataLink.stopSignal = &CG_state.c_ptr()->userterminationneeded;
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
                                double &energy, ArrayType &force,
                                int nr_x_values) {
  Mesh *mesh = dataLink->mesh;
  mesh->nrUpdateFunctionCalls++;

  // Update mesh position
  updateNodePositions(*mesh, disp);

  // Calculate energy and forces
  mesh->updateMesh();

  // Total energy, only used for minimization
  energy = mesh->totalEnergy;

  // Update forces in the minimization
  double maxForce = updateForceArray(mesh, force, nr_x_values);
  mesh->maxForce = maxForce;

  // Determine if the minimization is done
  if (maxForce < *dataLink->maxForceAllowed) {
    *dataLink->stopSignal = true;
  } else {
    // Sometimes, the minimization algorithm finds a state where max force is
    // low enough, but then moves away from it. Therefore, we want to re-set the
    // flag to false if it does so
    // This is not very well understood yet. Can be looked into.
    *dataLink->stopSignal = false;
  }
}

void alglibEnergyAndGradient(const alglib::real_1d_array &disp, double &energy,
                             alglib::real_1d_array &force, void *dataLinkPtr) {
  DataLink *dataLink = reinterpret_cast<DataLink *>(dataLinkPtr);
  updateMeshAndComputeForces(dataLink, disp, energy, force, force.length() / 2);
}

double FIREEnergyAndGradient(Eigen::VectorXd &disp, Eigen::VectorXd &force,
                             void *dataLinkPtr) {
  double energy;
  DataLink *dataLink = reinterpret_cast<DataLink *>(dataLinkPtr);
  updateMeshAndComputeForces(dataLink, disp, energy, force, force.size() / 2);
  return energy;
}

void Simulation::m_minimizeWithFIRE() {
  timer.Start("FIREMinimization");
  FIREpp::FIRESolver<double> s = FIREpp::FIRESolver<double>(FIRE_param);
  double energy;
  FIRERep.nrIter = s.minimize(FIREEnergyAndGradient, FIRENodeDisplacements,
                              energy, &dataLink, FIRERep.termType);
  FIRERep.nms = timer.Stop("FIREMinimization");
  FIRERep.nfev = mesh.nrUpdateFunctionCalls;

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
double updateForceArray(Mesh *mesh, ArrayType &force, int nr_x_values) {
  double maxForce = 0;

  // #pragma omp parallel for reduction(max : maxForce)
  for (int i = 0; i < nr_x_values; i++) {
    NodeId n_id = mesh->freeNodeIds[i];
    const auto &node = (*mesh)[n_id];

    // Store force values in the array
    double fx = node->f[0];
    double fy = node->f[1];

    force[i] = fx;
    force[nr_x_values + i] = fy;

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

void Simulation::writeToFile(bool forceWrite) {
  timer.Start("write");
  // We write to the cvs file every time this function is called
  writeToCsv(csvFile, (*this));
  // These are writing date much less often
  m_writeMesh(forceWrite);
  m_writeDump(forceWrite);
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
    writeMeshToVtu(mesh, name, dataPath);
    lastLoadWritten = mesh.load;
  }
}

void Simulation::m_writeDump(bool forceWrite) {
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
    saveSimulation();

    // Perhaps a bit strange, but this seems like a nice time to also
    // create/update the pvd file. (Sometimes it can be nice to have this
    // file before the simulation is done)
    gatherDataFiles();

    // Previously, I was updating the bounding box every step, but this is
    // overkill. Once every now and then like this is fine.
    mesh.updateBoundingBox();

    lastSaveTime = now;
    firstSaveDone = true;
  }
}

void Simulation::finishStep() {
  mesh.remesh();
  //  Calculate averages
  mesh.calculateAverages();
  // Update number of plastic events
  mesh.updateNrPlasticEvents();
  // Updates progress
  m_updateProgress();
  writeToFile();
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
  name = config.name;
  rows = config.rows;
  cols = config.cols;

  startLoad = config.startLoad;
  loadIncrement = config.loadIncrement;
  maxLoad = config.maxLoad;
  // In the simulations, we usually loop
  //   while (load < maxLoad) {...
  // Because of machine presision (i think), simulations often
  // go one step further than they need to. To prevent this,
  // we slightly lower the max load from what the user sets
  maxLoad -= loadIncrement / 100;
}

void Simulation::finishSimulation() {
  gatherDataFiles();
  // timer.PrintAllRuntimes();
}

void Simulation::gatherDataFiles() {
  // This creates a pvd file that links all the utv files together.
  createCollection(getDataPath(name, dataPath), getOutputPath(name, dataPath));
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
      std::filesystem::path(getDumpPath(name, dataPath)) / fileName;

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
  s.m_loadConfig(s.config);

  if (s.mesh.load >= s.maxLoad && !forceReRun) {
    std::cout << "Simulation already complete\n";
    exit(EXIT_SUCCESS);
  }

  saveConfigFile(s.config);
  s.csvFile = initCsvFile(s.name, s.dataPath, s);
  s.initialize();
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
  mesh->nrMinimizationItterations++;
  int it = mesh->nrMinimizationItterations;
  int nrF = mesh->nrUpdateFunctionCalls;

  // Check if iteration count is a multiple of 300
  if (mesh->nrMinimizationItterations == 0) {
    std::cout << "Warning: " << it << " iterations and " << nrF
              << " function calls.\n";
  }
  writeMeshToVtu(*mesh, mesh->simName, mesh->dataPath);
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