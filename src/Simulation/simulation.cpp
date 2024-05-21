#include "simulation.h"
#include "Data/dataExport.h"
#include "Data/paramParser.h"
#include "Eigen/src/Core/Matrix.h"
#include "cereal/archives/binary.hpp"
#include "settings.h"
#include <FIRE.h>
#include <Param.h>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>

Simulation::Simulation(Config config_, std::string _dataPath) {
  dataPath = _dataPath;

  // This function initializes a lot of the variables using the config file
  m_loadConfig(config_);

  timer = Timer();

  mesh = Mesh(rows, cols, config.usingPBC);
  mesh.load = config.startLoad;
  mesh.setSimNameAndDataPath(name, dataPath);

  clearOutputFolder(name, dataPath);
  createDataFolder(name, dataPath);
  saveConfigFile(config);

  // Create and open file
  csvFile = initCsvFile(name, dataPath);

  std::cout << config;
}

void Simulation::initialize() {
  // Initialization should be done after nodes have been moved and fixed as
  // desired. The elements created by the function below are copies and do not
  // dynamically update. (the update function only updates the position,
  // energy and stress)
  mesh.createElements();

  // we update the LBFGS solver
  initSolver();

  // Start simulation timer
  timer.Start();
  // Give some feedback that the process has started
  m_updateProgress();
}

void Simulation::initSolver() {

  // LBFGS Initialization preparation
  int nrFreeNodes = mesh.freeNodeIds.size();
  LBFGSNodeDisplacements.setlength(2 * nrFreeNodes);
  // Set values to zero
  for (int i = 0; i < 2 * nrFreeNodes; i++) {
    LBFGSNodeDisplacements[i] = 0;
  }

  // Adjust nrCorrections based on the number of free nodes
  // LBFGS doesn't like that nr corrections is larger than nr of free nodes
  if (config.LBFGSNrCorrections > 2 * nrFreeNodes) {
    config.LBFGSNrCorrections = 2 * nrFreeNodes;
  }

  // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgscreate
  // Initialize the state variable
  alglib::minlbfgscreate(config.LBFGSNrCorrections, LBFGSNodeDisplacements,
                         LBFGS_state);

  // FIRE Initialization
  FIRENodeDisplacements = VectorXd::Zero(2 * nrFreeNodes);
  FIRE_param = FIREpp::FIREParam<double>(mesh.nrNodes, 2);
  config.updateParam(FIRE_param);
}

void Simulation::minimize() {
  timer.Start("minimization");

  if (config.minimizer == "FIRE") {
    minimizeWithFIRE();
  } else if (config.minimizer == "LBFGS") {
    minimizeWithLBFGS();
  } else {
    std::cout << config.minimizer << std::endl;
    throw std::invalid_argument("Unknown solver");
  }
  if (LBFGS_report.terminationtype == -3) {
    writeToFile(true);
    throw std::runtime_error("Energy too high");
  }
  timer.Stop("minimization");
}

void Simulation::minimizeWithLBFGS() {
  timer.Start("AlglibMinimization");
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgsrestartfrom
  // We reset and reuse the state instead of initializing it again
  alglib::minlbfgsrestartfrom(LBFGS_state, LBFGSNodeDisplacements);

  // Set termination condition, ei. when is the solution good enough
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgssetcond
  alglib::minlbfgssetcond(LBFGS_state, config.LBFGSEpsg, config.LBFGSEpsf,
                          config.LBFGSEpsx, config.LBFGSMaxIterations);

  // This is where the heavy calculations happen
  // The null pointer can be replaced with a logging function
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgsoptimize
  alglib::minlbfgsoptimize(LBFGS_state, LBFGSEnergyAndGradient, nullptr, &mesh);

  // TODO Collecting and analysing these reports could be a usefull tool for
  // optimization
  alglib::minlbfgsresults(LBFGS_state, LBFGSNodeDisplacements, LBFGS_report);
  LBFGSRep.nms = timer.Stop("AlglibMinimization");
  LBFGSRep = SimReport(LBFGS_report);
}

template <typename ArrayType>
void updateMeshAndComputeForces(Mesh *mesh, const ArrayType &disp,
                                double &energy, ArrayType &force,
                                int nr_x_values) {
  mesh->nrUpdateFunctionCalls++;

  // Update mesh position
  updatePositionOfMesh(*mesh, disp);

  // Calculate energy and forces
  energy = calcEnergyAndForces(*mesh);

  // Update forces
  for (int i = 0; i < nr_x_values; i++) {
    NodeId n_id = mesh->freeNodeIds[i];
    force[i] = (*mesh)[n_id]->f[0];
    force[nr_x_values + i] = (*mesh)[n_id]->f[1];
  }
}

void LBFGSEnergyAndGradient(const alglib::real_1d_array &disp, double &energy,
                            alglib::real_1d_array &force, void *meshPtr) {
  Mesh *mesh = reinterpret_cast<Mesh *>(meshPtr);
  updateMeshAndComputeForces(mesh, disp, energy, force, force.length() / 2);
  // Additional operations specific to alglibCalcEnergyAndGradient can go here
}

double FIREEnergyAndGradient(Eigen::VectorXd &disp, Eigen::VectorXd &force,
                             void *meshPtr) {
  Mesh *mesh = reinterpret_cast<Mesh *>(meshPtr);
  double energy;
  updateMeshAndComputeForces(mesh, disp, energy, force, force.size() / 2);
  // Additional operations specific to FIREEnergyAndGradient can go here
  return energy;
}

void Simulation::minimizeWithFIRE() {
  timer.Start("FIREMinimization");
  FIREpp::FIRESolver<double> s = FIREpp::FIRESolver<double>(FIRE_param);
  double energy;
  FIRERep.nrIter = s.minimize(FIREEnergyAndGradient, FIRENodeDisplacements,
                              energy, &mesh, FIRERep.terminationType);
  FIRERep.nms = timer.Stop("FIREMinimization");
  FIRERep.nfev = mesh.nrUpdateFunctionCalls;

  // We first do FIRE, but then use the result as an initial guess for LBFGS
  setInitialGuess();
  minimizeWithLBFGS();
}

// Updates the forces on the nodes in the surface and returns the total
// energy from all the elements in the surface.
double calcEnergyAndForces(Mesh &mesh) {
  // First of all we need to make sure that the forces on the nodes have been
  // reset
  mesh.resetForceOnNodes();

  // Now we update all the elements using the current positions of the nodes
  mesh.updateElements();

  // We then add the force from the elements back to the nodes
  mesh.applyForceFromElementsToNodes();

  return mesh.calculateTotalEnergy();
}
// TODO do the same for force and guess
// Helper function to update positions using a generic buffer and its size
static void updateMeshPositions(Mesh &mesh, const double *data, size_t length) {
  // The displacement is structed like this: [x1,x2,x3,x4,y1,y2,y3,y4], so we
  // need to know where the "x" end and where the "y" begin.
  int nr_x_values = length / 2;
  Node *n;

  // We loop over all the free nodes
  for (size_t i = 0; i < mesh.freeNodeIds.size(); i++) {
    n = mesh[mesh.freeNodeIds[i]];
    // This function changes the position of the node based on the given
    // displacement and the current initial position.
    n->setDisplacement({data[i], data[i + nr_x_values]});
  }
}

// Overload for alglib::real_1d_array
void updatePositionOfMesh(Mesh &mesh, const alglib::real_1d_array &disp) {
  updateMeshPositions(mesh, disp.getcontent(), disp.length());
}

// Overload for Eigen::VectorXd
void updatePositionOfMesh(Mesh &mesh, const Eigen::VectorXd &disp) {
  updateMeshPositions(mesh, disp.data(), disp.size());
}

// This function modifies the nodeDisplacements variable used in the solver
void Simulation::setInitialGuess(Matrix2d guessTransformation) {
  // Our initial guess will be that all particles have shifted by the same
  // transformation as the border.
  // The disp is structed like this: [x1,x2,x3,x4,y1,y2,y3,y4], so we
  // need to know where the "x" end and where the "y" begin.
  int nr_x_values = LBFGSNodeDisplacements.length() / 2; // Shifts to y section

  const Node *n; // free node

  // We loop over all the nodes that are not on the border, ie. the innside
  // nodes.
  for (size_t i = 0; i < mesh.freeNodeIds.size(); i++) {
    n = mesh[mesh.freeNodeIds[i]];
    Vector2d next_displacement = guessTransformation * n->pos() - n->init_pos();
    LBFGSNodeDisplacements[i] = next_displacement[0];
    LBFGSNodeDisplacements[i + nr_x_values] = next_displacement[1];

    FIRENodeDisplacements[i] = next_displacement[0];
    FIRENodeDisplacements[i + nr_x_values] = next_displacement[1];
  }
}

// Core function to add noise to a double array
void addNoiseToArray(double *data, size_t length, double noise) {

  for (size_t i = 0; i < length; i++) {
    // Generate random noise in the range [-noise, noise]
    data[i] += ((double)rand() / RAND_MAX) * 2 * noise - noise;
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
  double effectiveNoise = customNoise == -1 ? config.noise : customNoise;
  addNoise(LBFGSNodeDisplacements, effectiveNoise);

  addNoise(FIRENodeDisplacements, effectiveNoise);
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

  progress = 100 * (mesh.load - startLoad) / (maxLoad - startLoad);

  // Always construct the progress message for logging
  int intProgress = static_cast<int>(progress);

  std::string consoleProgressMessage = std::to_string(intProgress) + "%" //
                                       + " RT: " + getRunTime()          //
                                       + " ETR: " + getEstimatedRemainingTime();

  // Construct a separate log message that includes the load and number of
  // plastic events
  std::string logProgressMessage = consoleProgressMessage                   //
                                   + " Load: " + std::to_string(mesh.load); //
  // + " dt_start: " + std::to_string(config.dtStart);

  // Use static variables to track the last progress and the last update time
  static int oldProgress = -1;
  static int firstProgress = -1;
  static int lastDump = -1;
  static auto lastUpdateTime = std::chrono::steady_clock::now();

  // Check if time since last update is more than 20 seconds or if progress has
  // changed
  auto now = std::chrono::steady_clock::now();
  int timeSinceLastUpdate =
      std::chrono::duration_cast<std::chrono::seconds>(now - lastUpdateTime)
          .count();

  if (config.showProgress == 1 &&
      (oldProgress != intProgress || timeSinceLastUpdate >= 20)) {
    // Update oldProgress and lastUpdateTime
    oldProgress = intProgress;
    if (firstProgress == -1)
      firstProgress = intProgress;

    lastUpdateTime = now; // Update the last update time

    // Output the progress message
    std::cout << logProgressMessage << std::endl;

    // We create a dump every 5 percent
    if (intProgress % 5 == 0 && intProgress != 0 && intProgress != lastDump &&
        firstProgress != intProgress) {
      lastDump = intProgress;
      writeToFile(true);

      timer.PrintAllRuntimes();
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
  // 180GB) OR If there are few large avalanvhes, we might go long without
  // saving data In order to get a good framerate for an animation, we want to
  // ensure that not too much happens between frames. This enures that we at
  // least have 200 frames of states over the course of loading
  static double lastLoadWritten = 0;
  if ((mesh.nrPlasticChanges >
       mesh.nrElements *
           config.plasticityEventThreshold) ||      // Lots of plastic change
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

    lastSaveTime = now;
    firstSaveDone = true;
  }
}

void Simulation::finishStep() {
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
  srand(config.seed);
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
  timer.PrintAllRuntimes();
}

void Simulation::gatherDataFiles() {
  // This creates a pvd file that links all the utv files together.
  createCollection(getDataPath(name, dataPath), getOutputPath(name, dataPath),
                   COLLECTIONNAME);
}

void Simulation::saveSimulation(std::string fileName_) {
  timer.Save();

  std::string fileName;
  if (fileName_ == "") {
    fileName = "Dump_l" + std::to_string(mesh.load) //
               + "_" + getCurrentDate() + ".mtsb";  // Read as MTS-Binary
  } else {
    fileName = fileName_ + ".mtsb";
  }

  std::ofstream ofs(getDumpPath(name, dataPath) + fileName,
                    std::ios::binary); // Make sure to open in binary mode
  cereal::BinaryOutputArchive oarchive(ofs);
  oarchive(*this); // Serialize the object to the output archive
  // This is also usefull to be able to see what simulation is running.
  std::cout << "Dump saved to: " << getDumpPath(name, dataPath) + fileName
            << std::endl;
}

void Simulation::loadSimulation(Simulation &s, const std::string &file,
                                const std::string &conf) {
  std::ifstream ifs(file, std::ios::binary); // Make sure to open in binary mode
  cereal::BinaryInputArchive iarchive(ifs);
  iarchive(s); // Serialize the object from the input archive

  // If we have a config file, we should load the new config
  if (conf != "") {
    s.config = parseConfigFile(conf);
    // Assert that mesh size has not been changed
    if (s.rows != s.config.rows || s.cols != s.config.cols) {
      std::invalid_argument("Mesh size cannot be changed!");
    }
  }

  s.m_loadConfig(s.config);
  // If we have changed the settings, we might need to make a new folder
  createDataFolder(s.name, s.dataPath);
  saveConfigFile(s.config);
  s.csvFile = initCsvFile(s.name, s.dataPath);
  s.initSolver();
  std::cout << s.config;
  s.timer.Start();
}

std::string Simulation::getRunTime() const { return timer.RunTimeString(); }

std::string Simulation::getEstimatedRemainingTime() const {
  return FormatDuration(calculateETR(timer.RunTime(), progress / 100));
}

// Function to calculate the Estimated Time Remaining (ETR) using progress
// fraction
std::chrono::milliseconds calculateETR(std::chrono::milliseconds elapsed,
                                       float progressFraction) {
  if (progressFraction <= 0) {
    return std::chrono::milliseconds(
        0); // Avoid division by zero if no progress
  }
  double elapsedSeconds = elapsed.count() / 1000.0;
  double rate = progressFraction / elapsedSeconds;
  if (rate == 0) {
    return std::chrono::milliseconds::min(); // Avoid infinity if rate
                                             // calculates to zero
  }
  long long etrInMilliseconds =
      static_cast<long long>(((1 - progressFraction) / rate) * 1000);

  // Ignore negative values
  if (etrInMilliseconds < 0) {
    etrInMilliseconds = 0;
  }
  return std::chrono::milliseconds(etrInMilliseconds);
}

void printReport(const alglib::minlbfgsreport &report) {
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgsresults
  std::cout << "Optimization Report:\n";
  std::cout << "\tIterations Count: " << report.iterationscount << '\n';
  std::cout << "\tNumber of Function Evaluations: " << report.nfev << '\n';
  std::cout << "\tTermination Reason: ";
  switch (report.terminationtype) {
  case -8:
    std::cout << "Infinite or NAN values in function/gradient";
    break;
  case -2:
    std::cout << "Rounding errors prevent further improvement";
    break;
  case -1:
    std::cout << "Incorrect parameters were specified";
    break;
  case 1:
    std::cout << "Relative function improvement is no more than EpsF";
    break;
  case 2:
    std::cout << "Relative step is no more than EpsX";
    break;
  case 4:
    std::cout << "Gradient norm is no more than EpsG";
    break;
  case 5:
    std::cout << "MaxIts steps was taken";
    break;
  case 7:
    std::cout << "Stopping conditions are too stringent, further improvement "
                 "is impossible";
    break;
  case 8:
    std::cout << "Terminated by user request";
    break;
  default:
    std::cout << "Unknown termination reason";
  }
  std::cout << std::endl;
}

// New method to print nodeDisplacements in (x, y) pairs
void printNodeDisplacementsGrid(alglib::real_1d_array nodeDisplacements) {
  int nr_x_values = nodeDisplacements.length() / 2;

  std::cout << "Node Displacements (x, y):" << std::endl;

  // Calculate the grid size for printing, assuming a rectangular (not
  // necessarily square) layout
  int gridSizeX = std::ceil(std::sqrt(nr_x_values)); // Width of the grid
  int gridSizeY =
      std::ceil(double(nr_x_values) /
                gridSizeX); // Height of the grid, ensuring all nodes fit

  for (int y = gridSizeY - 1; y >= 0;
       y--) { // Start from the bottom row to have (0,0) in the bottom left
    for (int x = 0; x < gridSizeX; x++) {
      int index = y * gridSizeX + x;
      if (index < nr_x_values) { // Ensure index is within the range of node
                                 // displacements
        std::cout << std::setw(10) << "(" << nodeDisplacements[index] << ", "
                  << nodeDisplacements[index + nr_x_values] << ") ";
      } else {
        // Print placeholders for grid positions without a corresponding node
        std::cout << std::setw(10) << "(--, --) ";
      }
    }
    std::cout << std::endl; // New line for each row of the grid
  }
}