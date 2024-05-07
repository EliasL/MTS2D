#ifndef SIMULATION_H
#define SIMULATION_H
#pragma once

#include <omp.h>
#include <sys/ioctl.h>
#include <unistd.h>

#include "Data/dataExport.h"
#include "Data/logging.h"
#include "Data/paramParser.h"

#include "Mesh/mesh.h"

// Alglib
#include <optimization.h>
#include <stdafx.h>

// Eigen
#include <Eigen/Core>

/**
 * @brief A object used to controll a loading simulation
 *
 * The object should be initialized with a yaml settings file.
 */
class Simulation {

public:
  // Initializes using a config file
  Simulation(Config config, std::string dataPath);
  Simulation() = default;

  // This should run when the mesh is properly pepared. (so we know which nodes
  // are fixed and which nodes are free.)
  void initialize();

  // If some changes are made to the number of fixed nodes mid-simulation, this
  // should be used.
  void initSolver();

  // Chooses a minimization method and keeps track of minimization time
  void minimize();

  // Uses minlbfgsoptimize to minimize the energy of the system.
  void minimize_with_alglib();

  // uses the FIRE algorithm to minimize the energy of the system.
  void minimize_with_FIRE();

  // Get the report from the minimization
  const alglib::minlbfgsreport &getReport() const;

  // Our initial guess will be that all particles have shifted by the same
  // transformation as the border.
  void setInitialGuess(Matrix2d guessTransformation);

  void addNoiseToGuess(double customNoise = -1);

  void finishStep();

  // Does some final touches and makes a collection of all the .vtu files in
  // the data folder.
  void finishSimulation();

  // Creates a pvd file that points to all the vtu files in the data folder.
  void gatherDataFiles();

  // Save the simulation to a binary file
  void saveSimulation();

  // Creates a vtu file of the current state of the simulation
  void writeToFile(bool forceWrite = false);

  static void loadSimulation(Simulation &s, const std::string &file);

  // gets run time
  std::string getRunTime() const;
  std::string getEstimatedRemainingTime() const;

  // The mesh we do our simulations on.
  Mesh mesh;

  // Loading parameters
  double startLoad;
  double loadIncrement;
  double maxLoad;
  // A percentage from 0 to 100 of loading completion
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

private:
  std::ofstream csvFile;

  // Amount of noise in the first inital guess
  double noise;

  // These values represents the current x and y displacements from the
  // initial position of the simulation
  alglib::real_1d_array alglibNodeDisplacements;
  VectorXd FIRENodeDisplacements;

  // M in the documentation-
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgscreate
  int nrCorrections;
  // Sets the scale of the
  double scale;
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgssetcond
  double epsg; // epsilon gradiant. A tolerance for how small the gradiant
               // should be before termination.
  double epsf; // epsilon function. A tolerance for how small the change in
               // the value of the function between itterations should be
               // before termination.
  double epsx; // epsilon x-step. A tolerance for how small the step between
               // itterations should be before termination.
  alglib::ae_int_t maxIterations; // Maximum itterations
  // When all the values above are chosen to be 0, a small epsx is
  // automatically chosen This is unphysical, so for accademic purposes, we
  // should instead use a small epsf value to guarantee a certan level of
  // accuracy.

  // Variables alglib uses to give feedback on what happens in the
  // optimization function
  alglib::minlbfgsstate state;
  alglib::minlbfgsreport report;

  // showProgress can be 0, 1. 0 is nothing, 1 is minimal, and 2 is no longer
  // used
  int showProgress;

  double plasticityEventThreshold;

  friend class cereal::access;
  template <class Archive> void serialize(Archive &ar);

  // Updates the progress (no physics)
  void m_updateProgress();

  // Logs the progress and writes data to disk
  void m_writeMesh(bool forceWrite = false);
  void m_writeDump(bool forceWrite = false);

  // reads the config values to local variables
  void m_readConfig(Config config);

  // Alglib doesn't like that nr corrections is larger than nr of free nodes
  void m_adjustNrCorrections();
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
 * @param displacement An array of displacement values for the nodes
 * @param energy The energy of the mesh
 * @param force The force on each node
 * @param meshPtr A pointer to our mesh
 */
void alglib_calc_energy_and_gradiant(const alglib::real_1d_array &displacement,
                                     double &energy,
                                     alglib::real_1d_array &force,
                                     void *meshPtr);

// Updates the forces on the nodes in the surface and returns the total
// energy from all the elements in the surface.
double calc_energy_and_forces(Mesh &mesh);

// Using the nodeDisplacements, we update the position of the nodes
void updatePositionOfMesh(Mesh &mesh,
                          const alglib::real_1d_array &displacement);
// Overload for Eigen::VectorXd
void updatePositionOfMesh(Mesh &mesh, const Eigen::VectorXd &disp);

// Creates a simple shear tranformation matrix
Matrix2d getShear(double load, double theta = 0);

// Adds a random vector with components between +-noise
void addNoise(alglib::real_1d_array &displacement, double noise);

void printReport(const alglib::minlbfgsreport &report);

// Function to calculate the Estimated Time Remaining (ETR) using progress
// fraction
std::chrono::milliseconds calculateETR(std::chrono::milliseconds elapsed,
                                       float progressFraction);

// Debug function to see nodeDisplacements
void printNodeDisplacementsGrid(alglib::real_1d_array nodeDisplacements);

#endif

template <class Archive> inline void Simulation::serialize(Archive &ar) {
  ar(mesh, config, dataPath, timer);
}