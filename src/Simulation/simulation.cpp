#include "simulation.h"
#include "Data/dataExport.h"
#include "Data/logging.h"
#include "Data/paramParser.h"
#include "Eigen/src/Core/Matrix.h"
#include "cereal/archives/binary.hpp"
#include "randomUtils.h"
#include <FIRE.h>
#include <Param.h>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <optimization.h>
#include <ostream>
#include <stdexcept>
#include <string>

Simulation::Simulation(Config config_, std::string _dataPath) {
  dataPath = _dataPath;

  // This function initializes a lot of the variables using the config file
  m_loadConfig(config_);

  timer = Timer();

  mesh = Mesh(rows, cols, 1, config.QDSD, config.usingPBC);
  mesh.load = config.startLoad;
  mesh.setSimNameAndDataPath(name, dataPath);

  if (simulationAlreadyComplete(name, dataPath, maxLoad) &&
      !config.forceReRun) {
    std::cout << "Simulation already complete\n";
    exit(EXIT_SUCCESS);
  }

  clearOutputFolder(name, dataPath);
  createDataFolder(name, dataPath);
  saveConfigFile(config);

  // Create and open file
  csvFile = initCsvFile(name, dataPath, *this);
}

void Simulation::initialize() {
  // Initialization should be done after nodes have been moved and fixed as
  // desired. The elements created by the function below are copies and do not
  // dynamically update. (the update function only updates the position,
  // energy and stress)
  mesh.createElements();

  // we update the solvers
  initSolver();

  // Start simulation timer
  timer.Start();
  // Give some feedback that the process has started
  m_updateProgress();
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
}

void Simulation::minimize() {
  timer.Start("minimization");

  // We explicitly reset the LBFGS_report because it won't be overwritten
  // if we use a different minimization algorithm
  if (mesh.loadSteps == 2) {
    LBFGS_report = alglib::minlbfgsreport();
    LBFGSRep = SimReport(LBFGS_report);
  }

  // If it is the first step, we always minimize with the same algorithm to
  // ensure each seed has the same STABLE starting point
  if (mesh.loadSteps == 1) {
    m_minimizeWithLBFGS();
  } else if (config.minimizer == "FIRE") {
    m_minimizeWithFIRE();
  } else if (config.minimizer == "LBFGS") {
    m_minimizeWithLBFGS();
  } else if (config.minimizer == "CG") {
    m_minimizeWithCG();
  } else {
    std::cout << config.minimizer << std::endl;
    throw std::invalid_argument("Unknown solver");
  }
  if (FIRERep.tType == -3) {
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

  // This is where the heavy calculations happen
  // The null pointer can be replaced with a logging function
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgsoptimize
  alglib::minlbfgsoptimize(LBFGS_state, alglibEnergyAndGradient, nullptr,
                           &mesh);

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

  // This is where the heavy calculations happen
  // The null pointer can be replaced with a logging function
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_mincgoptimize
  alglib::mincgoptimize(CG_state, alglibEnergyAndGradient, nullptr, &mesh);

  alglib::mincgresults(CG_state, alglibNodeDisplacements, CG_report);
  CGRep.nms = timer.Stop("CGMinimization");
  CGRep = SimReport(CG_report);
}

template <typename ArrayType>
void updateMeshAndComputeForces(Mesh *mesh, const ArrayType &disp,
                                double &energy, ArrayType &force,
                                int nr_x_values) {
  mesh->nrUpdateFunctionCalls++;

  // Update mesh position
  updateNodePositions(*mesh, disp);

  // Calculate energy and forces
  mesh->updateMesh();
  // Total energy, only used for minimization
  energy = mesh->totalEnergy;

  // Update forces
  updateForceArray(mesh, force, nr_x_values);
}

void alglibEnergyAndGradient(const alglib::real_1d_array &disp, double &energy,
                             alglib::real_1d_array &force, void *meshPtr) {
  Mesh *mesh = reinterpret_cast<Mesh *>(meshPtr);
  updateMeshAndComputeForces(mesh, disp, energy, force, force.length() / 2);
}

double FIREEnergyAndGradient(Eigen::VectorXd &disp, Eigen::VectorXd &force,
                             void *meshPtr) {
  Mesh *mesh = reinterpret_cast<Mesh *>(meshPtr);
  double energy;
  updateMeshAndComputeForces(mesh, disp, energy, force, force.size() / 2);
  return energy;
}

void Simulation::m_minimizeWithFIRE() {
  timer.Start("FIREMinimization");
  FIREpp::FIRESolver<double> s = FIREpp::FIRESolver<double>(FIRE_param);
  double energy;
  FIRERep.nrIter = s.minimize(FIREEnergyAndGradient, FIRENodeDisplacements,
                              energy, &mesh, FIRERep.tType);
  FIRERep.nms = timer.Stop("FIREMinimization");
  FIRERep.nfev = mesh.nrUpdateFunctionCalls;

  if (FIRERep.tType == -3) {
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

template <typename ArrayType>
void updateForceArray(Mesh *mesh, ArrayType &force, int nr_x_values) {
  // multithreading seems to be slower
  // #pragma omp parallel for
  for (int i = 0; i < nr_x_values; i++) {
    NodeId n_id = mesh->freeNodeIds[i];
    force[i] = (*mesh)[n_id]->f[0];
    force[nr_x_values + i] = (*mesh)[n_id]->f[1];
  }
}

// This function modifies the nodeDisplacements variable used in the solver
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

// TODO
// This can be used in the BFGS optimization
void iteration_logger(const alglib::real_1d_array &x, double func,
                      void *meshPtr) {
  (void)x; // Explicitly casting x to void to silence unused parameter warning
  (void)func;

  // Cast the void pointer back to Mesh pointer
  Mesh *mesh = reinterpret_cast<Mesh *>(meshPtr);

  // Increment the iteration count
  mesh->nrMinimizationItterations++;
  int it = mesh->nrMinimizationItterations;
  int nrF = mesh->nrUpdateFunctionCalls;

  // Check if iteration count is a multiple of 300
  if (mesh->nrMinimizationItterations % 300 == 0) {
    std::cout << "Warning: " << it << " iterations and " << nrF
              << " function calls.\n";
    writeMeshToVtu(*mesh, mesh->simName, mesh->dataPath);
  }
}

void Simulation::m_updateProgress() {
  using namespace std::chrono;

  progress = (mesh.load - startLoad) / (maxLoad - startLoad);

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
      (-mesh.delAvgEnergy > config.energyDropThreshold) || // Large energy drop
      (abs(mesh.load - lastLoadWritten) > 0.005) ||        // Absolute change
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
  // Calculate averages
  mesh.calculateAverages();
  // Update number of plastic events
  mesh.updateNrPlasticEvents();
  // Resets function counters (used for logging only)
  mesh.resetLoadingStepFunctionCounters();
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
  if (config.nrThreads == 0) {
    config.nrThreads = omp_get_max_threads();
  }
  omp_set_num_threads(config.nrThreads);

  // Assign values from Config to Simulation members
  name = config.name;
  rows = config.rows;
  cols = config.cols;

  startLoad = config.startLoad;
  loadIncrement = config.loadIncrement;
  maxLoad = config.maxLoad;
}

void Simulation::finishSimulation() {
  gatherDataFiles();
  // timer.PrintAllRuntimes();
}

void Simulation::gatherDataFiles() {
  // This creates a pvd file that links all the utv files together.
  createCollection(getDataPath(name, dataPath), getOutputPath(name, dataPath));
}

void Simulation::saveSimulation(std::string fileName_) {
  timer.Save();
  std::string fileName;
  if (fileName_ == "") {
    // Adjust epsilon proportionally to mesh.load
    double epsilon =
        std::numeric_limits<double>::epsilon() * std::abs(mesh.load) * 10.0;
    double roundedLoad = std::ceil((mesh.load - epsilon) * 10.0) / 10.0;

    // Use ostringstream to convert rounded load to string with fixed precision
    std::ostringstream oss;
    oss << std::fixed << roundedLoad;

    // Remove trailing zeros and the decimal point if necessary
    std::string loadStr = oss.str();
    loadStr.erase(loadStr.find_last_not_of('0') + 1, std::string::npos);
    if (loadStr.back() == '.') {
      loadStr.pop_back();
    }
    fileName = "dump_l" + loadStr + ".mtsb";
  } else {
    fileName = fileName_ + ".mtsb";
  }

  std::ofstream ofs(getDumpPath(name, dataPath) + fileName, std::ios::binary);
  cereal::BinaryOutputArchive oarchive(ofs);
  oarchive(*this);
  std::cout << "Dump saved to: " << getDumpPath(name, dataPath) + fileName
            << std::endl;
}

void Simulation::loadSimulation(Simulation &s, const std::string &file,
                                const std::string &conf, std::string outputPath,
                                const bool forceReRun) {
  // Open the file in binary mode
  std::ifstream ifs(file, std::ios::binary);

  // Check if the file exists and is open
  if (!ifs.is_open()) {
    // Handle the case where the file doesn't exist
    std::cerr << "Error: File '" << file << "' could not be opened."
              << std::endl;
    throw std::runtime_error("File does not exist or cannot be opened.");
  }

  // If the file is open, proceed with deserialization
  cereal::BinaryInputArchive iarchive(ifs);
  iarchive(s); // Deserialize the object from the input archive

  // Load config file
  s.config = parseConfigFile(conf);
  // Assert that mesh size has not been changed
  if (s.rows != s.config.rows || s.cols != s.config.cols) {
    std::ostringstream errorMessage;
    errorMessage << "Mesh size cannot be changed. Loaded: (" << s.rows << ", "
                 << s.cols << ") does not match config: (" << s.config.rows
                 << ", " << s.config.cols << ")";
    throw std::invalid_argument(errorMessage.str());
  }
  // Update to new dataPath if provided
  if (!outputPath.empty()) {
    s.dataPath = outputPath;
  }

  s.m_loadConfig(s.config);

  if (simulationAlreadyComplete(s.name, s.dataPath, s.maxLoad) && !forceReRun) {
    std::cout << "Simulation already complete\n";
    exit(EXIT_SUCCESS);
  }

  // If we have changed the settings, we might need to make a new folder
  createDataFolder(s.name, s.dataPath);

  saveConfigFile(s.config);
  s.csvFile = initCsvFile(s.name, s.dataPath, s);
  s.initSolver();
  s.timer.Start();
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