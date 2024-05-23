#ifndef SIMULATION_H
#define SIMULATION_H
#include <cstddef>
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
  void minimizeWithLBFGS();

  // uses the FIRE algorithm to minimize the energy of the system.
  void minimizeWithFIRE();

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

  // Save the simulation to a binary file. Leave fileName empty for default name
  void saveSimulation(std::string fileName_ = "");

  // Creates a vtu file of the current state of the simulation
  void writeToFile(bool forceWrite = false);

  static void loadSimulation(Simulation &s, const std::string &file,
                             const bool forceOverWrite = false,
                             const std::string &conf = "");

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

  // A report to gather information about the minimization
  SimReport FIRERep;
  SimReport LBFGSRep;

private:
  std::ofstream csvFile;

  // These values represents the current x and y displacements from the
  // initial position of the simulation
  alglib::real_1d_array LBFGSNodeDisplacements;
  VectorXd FIRENodeDisplacements;

  // Variables alglib uses to give feedback on what happens in the
  // optimization function
  alglib::minlbfgsstate LBFGS_state;
  alglib::minlbfgsreport LBFGS_report;

  // FIRE parameters
  FIREpp::FIREParam<double> FIRE_param;

  friend class cereal::access;
  template <class Archive> void serialize(Archive &ar);

  // Updates the progress (no physics)
  void m_updateProgress();

  // Logs the progress and writes data to disk
  void m_writeMesh(bool forceWrite = false);
  void m_writeDump(bool forceWrite = false);

  // reads the config values to local variables
  void m_loadConfig(Config config);
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
void updateMeshAndComputeForces(Mesh *mesh, const ArrayType &disp,
                                double &energy, ArrayType &force,
                                int nr_x_values);

// The two following functions use updateMeshAndComputeForces
void LBFGSEnergyAndGradient(const alglib::real_1d_array &displacement,
                            double &energy, alglib::real_1d_array &force,
                            void *meshPtr);
double FIREEnergyAndGradient(Eigen::VectorXd &disp, Eigen::VectorXd &force,
                             void *meshPtr);

// Updates the forces on the nodes in the surface and returns the total
// energy from all the elements in the surface.
double calcEnergyAndForces(Mesh &mesh);

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