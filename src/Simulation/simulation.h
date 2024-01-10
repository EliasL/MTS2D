#ifndef SIMULATION_H
#define SIMULATION_H
#pragma once

#include <vector>

#include "settings.h"
#include "Matrix/matrix2x2.h"
#include "Mesh/mesh.h"
#include "Simulation/meshManipulations.h"
#include "Data/dataExport.h"
#include "easylogging++.h"

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
    std::cout << "Iterations Count: " << report.iterationscount << '\n';
    std::cout << "Number of Function Evaluations: " << report.nfev << '\n';
    std::cout << "Termination Reason: ";
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

void updatePossitionOfMesh(Mesh &mesh, const alglib::real_1d_array displacement)
{

    // The displacement is structed like this: [x1,x2,x3,x4,y1,y2,y3,y4], so we
    // need to know where the "x" end and where the "y" begin.
    int nr_x_elements = displacement.length() / 2; // Shifts to y section

    Node *innside_node;

    // We loop over all the nodes that are not on the border, ie. the innside nodes.
    for (size_t i = 0; i < mesh.nonBorderNodeIds.size(); i++)
    {
        innside_node = mesh[mesh.nonBorderNodeIds[i]];
        innside_node->x += displacement[i];
        innside_node->y += displacement[i + nr_x_elements];
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

    // TODO Parllelalize this for-loop
    // #pragma omp parallel
    // #pragma omp for reduction(+:energy_thread)
    for (size_t i = 0; i < mesh.nrElements; i++)
    {
        // TODO, we still read the state of nodes from multiple elements
        // at the same time, but this seems to be okay according to ChatGPT.
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
    LOG(INFO) << displacement[0] << ", " << displacement[1];

    // Calculate energy and forces
    energy = calc_energy_and_forces(*mesh);
    LOG(INFO) << energy;

    // Update forces
    for (size_t i = 0; i < nr_x_values; i++)
    {
        NodeId a_id = mesh->nonBorderNodeIds[i];
        force[i] = (*mesh)[a_id]->f_x;
        force[nr_x_values + i] = (*mesh)[a_id]->f_y;
        LOG(INFO) << a_id.i;
    }

    LOG(INFO) << force[0] << ", " << force[1];
    LOG(INFO) << "\n";
    writeToVtu(*mesh);
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
    const Node *innside_node;

    // We loop over all the nodes that are not on the border, ie. the innside nodes.
    for (size_t i = 0; i < mesh.nonBorderNodeIds.size(); i++)
    {
        innside_node = mesh[mesh.nonBorderNodeIds[i]];
        transformed_node = transform(transformation, *innside_node); // F * node.position
        // Subtract the initial possition to get the displacement
        translateInPlace(transformed_node, *innside_node, -1); // node1.position - node2.position
        displacement[i] = transformed_node.x;
        displacement[i + nr_x_elements] = transformed_node.y;
    }
}

void run_simulation()
{
    int s = 3; // Square length, Can't be smaller than 3, because with 2, all nodes would be fixed
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

    double loadIncrement = 0.001;
    double maxLoad = 0.01;
    double theta = 0;

    // Boundary conditon transformation
    Matrix2x2<double> wrongBcTransform = Matrix2x2<double>::identity();
    Matrix2x2<double> bcTransform = getShear(loadIncrement, theta);

    clearOutputFolder();
    setLoggingOutput();
    createDataFolder();
    mesh.nodes[1][1].setPos(1.2, 1);
    writeToVtu(mesh);

    for (double load = 0; load < maxLoad; load += loadIncrement)
    {
        // mesh.applyTransformationToBoundary(bcTransform);
        mesh.load = load;
        // Modifies nodeDisplacements
        initialGuess(mesh, wrongBcTransform, nodeDisplacements);

        alglib::minlbfgscreate(1, nodeDisplacements, state);

        // Set termination condition, ei. when is the solution good enough
        // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgssetcond
        alglib::minlbfgssetcond(state, epsg, epsf, epsx, maxits);

        // The null pointer can be replaced with a logging function!
        alglib::minlbfgsoptimize(state, alglib_calc_energy_and_gradiant, nullptr, &mesh);

        // TODO Collecting and analysing these reports could be a usefull tool for optimization
        alglib::minlbfgsresults(state, nodeDisplacements, report);

        printReport(report);

        writeToVtu(mesh);
    }
    // Note that you can't / don't need to use + between two defined strings
    leanvtk::createCollection(OUTPUTFOLDERPATH SUBFOLDERPATH DATAFOLDERPATH);
}

#endif
