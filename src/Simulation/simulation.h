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

// Cereal
#include <cereal/types/string.hpp>

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

  // Save the simulation to a XML file. Leave fileName empty for default name
  // Returns the path of the XML file
  std::string saveSimulation(std::string fileName_ = "");

  // Creates a vtu file of the current state of the simulation
  void writeToFile(bool forceWrite = false);

  static void loadSimulation(Simulation &s, const std::string &file,
                             const std::string &conf, std::string outputPath,
                             const bool forceOverWrite = false);

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
void alglibEnergyAndGradient(const alglib::real_1d_array &displacement,
                             double &energy, alglib::real_1d_array &force,
                             void *meshPtr);
double FIREEnergyAndGradient(Eigen::VectorXd &disp, Eigen::VectorXd &force,
                             void *meshPtr);

// Using the nodeDisplacements, we update the position of the nodes
void updateNodePositions(Mesh &mesh, const alglib::real_1d_array &displacement);
// Overload for Eigen::VectorXd
void updateNodePositions(Mesh &mesh, const Eigen::VectorXd &disp);

template <typename ArrayType>
void updateForceArray(Mesh *mesh, ArrayType &force, int nr_x_values);

// Creates a simple shear tranformation matrix
Matrix2d getShear(double load, double theta = 0);

#endif

template <class Archive> void Simulation::serialize(Archive &ar) {
  ar(cereal::make_nvp("rows", rows), cereal::make_nvp("cols", cols),
     cereal::make_nvp("mesh", mesh), cereal::make_nvp("dataPath", dataPath),
     cereal::make_nvp("timer", timer), cereal::make_nvp("name", name),
     cereal::make_nvp("config", config));
}