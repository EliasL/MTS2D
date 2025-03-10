#ifndef SIMULATION_H
#define SIMULATION_H
#include <string>
#pragma once

#include <omp.h>
#include <sys/ioctl.h>
#include <unistd.h>

#include "Data/cereal_help.h"
#include "Data/data_export.h"
#include "Data/logging.h"
#include "Data/param_parser.h"

#include "Mesh/mesh.h"

// Alglib
#include <optimization.h>
#include <stdafx.h>

// Eigen
#include <Eigen/Core>

// Cereal
#include <cereal/types/string.hpp>

class Simulation;

/**
 * @brief A object used to provide access to data inside the minimization loop
 */
struct DataLink {
  // Gives access to the simulation variables
  Simulation *s;
  // The mesh we are minimizing the energy of
  Mesh *mesh;

  // Alglib has a state.userterminationneeded flag to stop the minimization.
  // Connect this pointer to any stop flag within a minimization algorithm and
  // stop the minimization.
  bool *stopSignal;

  // The mesh is provided to the minimization function, and sometimes, it is
  // nice to have access to the minimization state. We can get that here
  MinState *minState;
  alglib::minlbfgsstate *LBFGS_state;
  alglib::mincgstate *CG_state;

  // We minimize until all the components of the forces between all the nodes
  // is smaller than this value (Unless another stopping criteria has been
  // reached)
  double *maxForceAllowed;

  DataLink() {};
  DataLink(Simulation *s);
};

/**
 * @brief A object used to controll a loading simulation
 *
 * The object should be initialized with a yaml settings file.
 */
class Simulation {

public:
  // Initializes using a config file
  Simulation(Config config, std::string dataPath, bool cleanDataPath = false);
  Simulation() = default;

  // This should run when the mesh is properly pepared. (so we know which nodes
  // are fixed and which nodes are free.)
  void initialize();

  // Sets nrFreeNodes and prepares params and displacement vectors for the
  // minimization algorithms.
  // If some changes are made to the number of fixed nodes mid-simulation, this
  // should be used.
  void initSolver();

  // The first step is special. In order to get the same state across many
  // different settings, we always use the same settings and the same
  // minimization algorithm. That will give simulations (with the same seed) a
  // common stable starting point (as opposed to a common UNSTABLE starting
  // point)
  void firstStep();

  bool keepLoading();

  // Chooses a minimization method and keeps track of minimization time
  void minimize();

  // Our initial guess will be that all particles have shifted by the same
  // transformation as the border.
  void setInitialGuess(Matrix2d guessTransform = Eigen::Matrix2d::Identity());

  void addNoiseToGuess(double customNoise = -1);

  void finishStep();

  // Does some final touches and makes a collection of all the .vtu files in
  // the data folder.
  void finishSimulation();

  // Creates a pvd file that points to all the vtu files in the data folder.
  void gatherDataFiles();

  // Save the simulation to a XML file. Leave fileName empty for default name.
  // Returns the path of the XML file.
  std::string saveSimulation(std::string fileName_ = "");

  // Creates a vtu file of the current state of the simulation
  void writeToFile(bool forceWrite = false, std::string name = "");

  static void loadSimulation(Simulation &s, const std::string &file,
                             const std::string &conf, std::string outputPath,
                             const bool forceOverWrite = false);

  // Object used to provide access to various values inside the minimization
  // loop
  DataLink dataLink;

  // The mesh we do our simulations on.
  Mesh mesh;

  // Loading parameters
  double startLoad;
  double loadIncrement;
  double maxLoad;
  // A number from 0 to 1 of loading completion
  double progress;
  // Dimension of mesh
  int rows, cols;

  // Folder name
  std::string name;
  // Path to the output data
  std::string dataPath;

  // Config object
  Config config;

  // Timer to log simulation time
  Timer timer;

  MinState minState;

  // A report to gather information about the minimization
  SimReport FIRERep;
  SimReport LBFGSRep;
  SimReport CGRep;

private:
  // Uses minlbfgsoptimize to minimize the energy of the system.
  void m_minimizeWithLBFGS();

  // uses the FIRE algorithm to minimize the energy of the system.
  void m_minimizeWithFIRE();

  // Uses the conjugate gradient algorithm to minimize the energy of the system.
  void m_minimizeWithCG();

  // The csv file where we write meta data about each simulation step
  std::ofstream csvFile;

  // The csv file where we write meta data about the internals steps of the
  // minimization algorithm
  std::ofstream minCsvFile;

  // Variables alglib uses to give feedback on what happens in the
  // optimization function
  alglib::minlbfgsstate LBFGS_state;
  alglib::minlbfgsreport LBFGS_report;

  alglib::mincgstate CG_state;
  alglib::mincgreport CG_report;

  // FIRE parameters
  FIREpp::FIREParam<double> FIRE_param;

  // These values represents the current x and y displacements from the
  // initial position of the simulation
  alglib::real_1d_array alglibNodeDisplacements;
  VectorXd FIRENodeDisplacements;

  friend class cereal::access;
  template <class Archive> void serialize(Archive &ar);

  // Updates the progress (no physics)
  void m_updateProgress();

  // Logs the progress and writes data to disk
  void m_writeMesh(bool forceWrite = false);
  void m_writeDump(bool forceWrite = false, std::string name = "");

  // reads the config values to local variables
  void m_loadConfig(Config config);

  // Give DataLink access to private variables
  friend struct DataLink;
};

/**
 * @brief Calculates the energy and forces using the current state of the
 * mesh, pluss a displacement. The first time this function is called, this
 * displacement is the initial guess that we chose.
 *
 * First we update the positions of the mesh, ie, we nudge each node using the
 * values in displacement. Then we calculate the energy and forces for these
 * new positions. Finally, we update the forces so that the minimization
 * function knows how much, and in what direction to nudge the nodes in for
 * the next simulation step.
 *
 * @param mesh A pointer to our mesh
 * @param displacement An array of displacement values for the nodes
 * @param energy The energy of the mesh
 * @param force The force on each node
 */
template <typename ArrayType>
void updateMeshAndComputeForces(DataLink *dataLink, const ArrayType &disp,
                                double &energy, ArrayType &force,
                                int nr_x_values);

// The two following functions use updateMeshAndComputeForces
void alglibEnergyAndGradient(const alglib::real_1d_array &displacement,
                             double &energy, alglib::real_1d_array &force,
                             void *dataLink);
double FIREEnergyAndGradient(Eigen::VectorXd &disp, Eigen::VectorXd &force,
                             void *dataLink);

// Using the nodeDisplacements, we update the position of the nodes
void updateNodePositions(DataLink &dataLink,
                         const alglib::real_1d_array &displacement);
// Overload for Eigen::VectorXd
void updateNodePositions(DataLink &dataLink, const Eigen::VectorXd &disp);

// Returns the max force component found (used for stopping criteria)
template <typename ArrayType>
double updateForceArray(Mesh *mesh, ArrayType &force, int nr_x_values);

// Creates a simple shear tranformation matrix
Matrix2d getShear(double load, double theta = 0);

// Logs information about non-equilibrium states that occur during minimization.
// Creates a new folder for each loading step during the simulation.
// Warning! Can create a lot of data
void iterationLogger(const alglib::real_1d_array &x, double energy,
                     void *dataLink);

template <class Archive> void Simulation::serialize(Archive &ar) {
  ar(MAKE_NVP(rows), MAKE_NVP(cols), MAKE_NVP(mesh), MAKE_NVP(dataPath),
     MAKE_NVP(timer), MAKE_NVP(name), MAKE_NVP(config));
}

#endif