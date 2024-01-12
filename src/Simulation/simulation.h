#ifndef SIMULATION_H
#define SIMULATION_H
#pragma once

#include <vector>

#include "settings.h"
#include "Matrix/matrix2x2.h"
#include "Mesh/mesh.h"
#include "Simulation/meshManipulations.h"
#include "Data/dataExport.h"
#include "spdlog/spdlog.h"

#include "stdafx.h"
#include "interpolation.h"
#include "specialfunctions.h"
#include "linalg.h"
#include "statistics.h"
#include "alglibmisc.h"

void printReport(const alglib::minlbfgsreport &report)
{
    // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgsoptimize
    std::cout << "Optimization Report:\n";
    std::cout << "\tIterations Count: " << report.iterationscount << '\n';
    std::cout << "\tNumber of Function Evaluations: " << report.nfev << '\n';
    std::cout << "\tTermination Reason: ";
    switch (report.terminationtype)
    {
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
        std::cout << "Stopping conditions are too stringent, further improvement is impossible";
        break;
    case 8:
        std::cout << "Terminated by user request";
        break;
    default:
        std::cout << "Unknown termination reason";
    }
    std::cout << std::endl;
}

void updatePossitionOfMesh(Mesh &mesh, const alglib::real_1d_array displacement,
                           double multiplier = 1)
{

    // The displacement is structed like this: [x1,x2,x3,x4,y1,y2,y3,y4], so we
    // need to know where the "x" end and where the "y" begin.
    int nr_x_elements = displacement.length() / 2; // Shifts to y section

    Node *n; // Non border, inside node

    // We loop over all the nodes that are not on the border, ie. the innside nodes.
    for (size_t i = 0; i < mesh.nonBorderNodeIds.size(); i++)
    {
        n = mesh[mesh.nonBorderNodeIds[i]];
        n->x = n->init_x + displacement[i] * multiplier;
        n->y = n->init_y + displacement[i + nr_x_elements] * multiplier;
    }
}

// Updates the forces on the nodes in the surface and returns the total
// energy from all the elements in the surface.
double calc_energy_and_forces(Mesh &mesh)
{
    // First of all we need to make sure that the forces on the nodes have been
    // reset
    mesh.resetForceOnNodes();

    // This is the total energy from all the triangles
    double total_energy = 0;

    // TODO We could check if we can make total energy a reduced variable, 
    // and make addForce a omp critical function and test for performance gains 
    #pragma omp parallel
    #pragma omp for
    for (size_t i = 0; i < mesh.nrElements; i++)
    {
        mesh.elements[i].update();
    }

    // Now that the forces on the nodes have been reset, and
    // the elements updated in parallel, we can sum up the energy
    // and forces.
    for (size_t i = 0; i < mesh.nrElements; i++)
    {
        mesh.elements[i].applyForcesOnNodes();
        total_energy += mesh.elements[i].energy;
    }
    return total_energy;
}

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
                                     alglib::real_1d_array &force, void *meshPtr)
{

    // Cast the void pointer back to Mesh pointer
    Mesh *mesh = reinterpret_cast<Mesh *>(meshPtr);

    int nr_x_values = force.length() / 2;

    // Update mesh position
    updatePossitionOfMesh(*mesh, displacement);

    // Calculate energy and forces
    energy = calc_energy_and_forces(*mesh);

    // Update forces
    for (size_t i = 0; i < nr_x_values; i++)
    {
        NodeId a_id = mesh->nonBorderNodeIds[i];
        force[i] = (*mesh)[a_id]->f_x;
        force[nr_x_values + i] = (*mesh)[a_id]->f_y;
    }
    // writeToVtu(*mesh, "minimizing");
}

// Our initial guess will be that all particles have shifted by the same
// transformation as the border.
void initialGuess(const Mesh &mesh, const Matrix2x2<double> &transformation,
                  alglib::real_1d_array &displacement)
{
    // The displacement is structed like this: [x1,x2,x3,x4,y1,y2,y3,y4], so we
    // need to know where the "x" end and where the "y" begin.
    int nr_x_elements = displacement.length() / 2; // Shifts to y section

    Node transformed_node; // These are temporary variables for readability
    const Node *n;         // Non border, inside node

    // We loop over all the nodes that are not on the border, ie. the innside nodes.
    for (size_t i = 0; i < mesh.nonBorderNodeIds.size(); i++)
    {
        n = mesh[mesh.nonBorderNodeIds[i]];
        transformed_node = transform(transformation, *n);
        // Subtract the initial possition to get the displacement
        translateInPlace(transformed_node, n->init_x, n->init_y, -1.0); // node1.position - node2.position
        displacement[i] = transformed_node.x;
        displacement[i + nr_x_elements] = transformed_node.y;
    }
}

void addNoise(alglib::real_1d_array &displacement, double noise)
{
    int nr_x_elements = displacement.length() / 2; // Shifts to y section

    Node transformed_node; // These are temporary variables for readability
    const Node *n;         // Non border, inside node

    // We loop over all the nodes that are not on the border, ie. the innside nodes.
    for (size_t i = 0; i < nr_x_elements; i++)
    {        
        // Generate random noise in the range [-noise, noise]
        double noise_x = ((double)rand() / RAND_MAX) * 2 * noise - noise;
        double noise_y = ((double)rand() / RAND_MAX) * 2 * noise - noise;

        // Add noise to the displacement
        displacement[i] += noise_x;
        displacement[i + nr_x_elements] += noise_y;
    }
}

void run_simulation()
{
    int s = 50; // Square length, Can't be smaller than 3, because with 2, all nodes would be fixed
    int nx = s;
    int ny = s;
    int n = nx * ny;

    int nrThreads = NUMEROFTHREADS;
    // If we have a small number of particles, we need fewer threads
    // (we can't have more threads than elements)
    nrThreads = (nrThreads > s) ? nrThreads : s;

    Mesh mesh = {nx, ny};

    alglib::real_1d_array nodeDisplacements;

    // We set the length of the array to be the same as the number of nodes that
    // are not on the boarder. These are the only nodes we will modify.
    int nrNonBorderNodes = mesh.nonBorderNodeIds.size();
    nodeDisplacements.setlength(2 * nrNonBorderNodes);

    // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgssetcond
    double epsg = 0;             // epsilon gradiant. A tolerance for how small the gradiant should be before termination.
    double epsf = 0;             // epsilon function. A tolerance for how small the change in the value of the function between itterations should be before termination.
    double epsx = 0;             // epsilon x-step. A tolerance for how small the step between itterations should be before termination.
    alglib::ae_int_t maxits = 0; // Maximum itterations
    // When all the values above are chosen to be 0, a small epsx is automatically chosen
    // This is unphysical, so for accademic purposes, we should instead use a small epsf value to guarantee a certan level of accuracy.
    alglib::minlbfgsstate state;
    alglib::minlbfgsreport report;

    double startLoad = 0.15;
    double loadIncrement = 0.01;
    double maxLoad = 0.7;

    // Boundary conditon transformation
    Matrix2x2<double> loadStepTransform = getShear(loadIncrement);

    clearOutputFolder();
    createDataFolder();

    // Prepare initial load condition
    mesh.applyTransformation(getShear(startLoad));

    for (double load = startLoad; load < maxLoad; load += loadIncrement)
    {
        // We shift the boundary nodes according to the loadIncrement
        mesh.applyTransformationToBoundary(loadStepTransform);

        // This is only used for logging purposes (no physics)
        mesh.load = load;

        // Modifies nodeDisplacements
        initialGuess(mesh, loadStepTransform, nodeDisplacements);
        // Give the guess some noise
        addNoise(nodeDisplacements, 0.05);

        alglib::minlbfgscreate(1, nodeDisplacements, state);

        // Set termination condition, ei. when is the solution good enough
        // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgssetcond
        alglib::minlbfgssetcond(state, epsg, epsf, epsx, maxits);

        // The null pointer can be replaced with a logging function!
        alglib::minlbfgsoptimize(state, alglib_calc_energy_and_gradiant, nullptr, &mesh);

        // TODO Collecting and analysing these reports could be a usefull tool for optimization
        alglib::minlbfgsresults(state, nodeDisplacements, report);

        // printReport(report);
        spdlog::info("{}",load);
        writeToVtu(mesh, "Relaxed");
    }

    leanvtk::createCollection(OUTPUTFOLDERPATH SUBFOLDERPATH DATAFOLDERPATH);
}

#endif
