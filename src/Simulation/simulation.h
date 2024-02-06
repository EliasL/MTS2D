#ifndef SIMULATION_H
#define SIMULATION_H
#pragma once

#include <vector>
#include <sys/ioctl.h>
#include <unistd.h>
#include <omp.h>
#include <optional>

#include "settings.h"
#include "Matrix/matrix2x2.h"
#include "Mesh/mesh.h"
#include "Data/dataExport.h"
#include "Data/paramParser.h"
#include "Data/logging.h"
#include "spdlog/spdlog.h"
#include <indicators/block_progress_bar.hpp>
#include <indicators/progress_bar.hpp>
#include <indicators/cursor_control.hpp>

#include "stdafx.h"
#include "interpolation.h"
#include "specialfunctions.h"
#include "linalg.h"
#include "statistics.h"
#include "alglibmisc.h"
#include "iostream"

/**
 * @brief A object used to controll a loading simulation
 *
 * The object should be initialized with a yaml settings file.
 */
class Simulation
{
public:
    // Initializes using a config file
    Simulation(std::string configFile, std::optional<std::string> dataPath);

    // Main run function
    void run_simulation();

private:
    std::string name;
    std::string dataPath;
    // nx is the number of nodes in the x direction, likewise for ny.
    int nx, ny;

    // Loading parameters
    double startLoad;
    double loadIncrement;
    double maxLoad;
    double noise;

    // Boundary conditon transformation
    Matrix2x2<double> loadStepTransform;

    // The mesh we do our simulations on
    Mesh mesh;

    // These values represents the current x and y displacements from the initial
    // position of the simulation
    alglib::real_1d_array nodeDisplacements;

    // M in the documentation- https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgscreate
    int nrCorrections;
    // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgssetcond
    double epsg;                    // epsilon gradiant. A tolerance for how small the gradiant should be before termination.
    double epsf;                    // epsilon function. A tolerance for how small the change in the value of the function between itterations should be before termination.
    double epsx;                    // epsilon x-step. A tolerance for how small the step between itterations should be before termination.
    alglib::ae_int_t maxIterations; // Maximum itterations
    // When all the values above are chosen to be 0, a small epsx is automatically chosen
    // This is unphysical, so for accademic purposes, we should instead use a small epsf value to guarantee a certan level of accuracy.

    // Variables alglib uses to give feedback on what happens in the optimization function
    alglib::minlbfgsstate state;
    alglib::minlbfgsreport report;

    // showProgress can be 0, 1 or 2. 0 is nothing, 1 is minimal, and 2 is showing a full progress bar
    int showProgress;

    // Timer to log simulation time
    Timer timer;
    double plasticityEventThreshold;

    // Uses minlbfgsoptimize to minimize the energy of the system
    void m_minimize_with_alglib();

    // Our initial guess will be that all particles have shifted by the same
    // transformation as the border.
    void m_initialGuess();

    // Updates the progress bar (Just visual, no physics)
    void m_updateProgress(double load);

    // Logs the progress and writes data to disk
    void m_writeToDisk(double load);

    // Does some final touches and makes a collection of all the .vtu files in
    // the data folder
    void m_exit();
};

/**
 * @brief Calculates the energy and forces using the current state of the mesh,
 *  pluss a displacement. The first time this function is called, this
 * displacement is the initial guess that we chose.
 *
 * First we update the positions of the mesh, ie, we nudge each node using the
 * values in displacement. Then we calculate the energy and forces for these new positions.
 * Finally, we update the forces so that the minimization function knows how
 * much, and in what direction to nudge the nodes in for the next simulation step.
 *
 * @param displacement An array of displacement values for the nodes
 * @param energy The energy of the mesh
 * @param force The force on each node
 * @param meshPtr A pointer to our mesh
 */
void alglib_calc_energy_and_gradiant(const alglib::real_1d_array &displacement,
                                     double &energy,
                                     alglib::real_1d_array &force, void *meshPtr);

// Updates the forces on the nodes in the surface and returns the total
// energy from all the elements in the surface.
double calc_energy_and_forces(Mesh &mesh);

// Using the nodeDisplacements, we update the position of the nodes
void updatePossitionOfMesh(Mesh &mesh, const alglib::real_1d_array &displacement);

// Creates a simple shear tranformation matrix
Matrix2x2<double> getShear(double load, double theta = 0);

// Adds a random vector with components between +-noise
void addNoise(alglib::real_1d_array &displacement, double noise);

void printReport(const alglib::minlbfgsreport &report);

indicators::BlockProgressBar &getBar();


// Function to calculate the Estimated Time Remaining (ETR) using progress fraction
long long calculateETR(long long elapsedMilliseconds, float progressFraction);
#endif
