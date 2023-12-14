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

void updatePossitionOfMesh(Mesh &mesh, const alglib::real_1d_array displacement, double scalar = 1.0)
{

    // The displacement is structed like this: [x1,x2,x3,x4,y1,y2,y3,y4], so we
    // need to know where the "x" end and where the "y" begin.
    int nr_x_elements = displacement.length() / 2; // Shifts to y section

    Node *innside_node;

    // We loop over all the nodes that are not on the border, ie. the innside nodes.
    for (size_t i = 0; i < mesh.nonBorderNodeIds.size(); i++)
    {
        innside_node = mesh[mesh.nonBorderNodeIds[i]];
        innside_node->x += displacement[i] * scalar;
        innside_node->y += displacement[i + nr_x_elements] * scalar;
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
    double total_energy;

    // TODO Parllelalize this for-loop
    // #pragma omp parallel
    // #pragma omp for reduction(+:energy_thread)
    for (size_t i = 0; i < mesh.nrElements; i++)
    {
        // TODO, we still read the state of nodes from multiple elements
        // at the same time, but this seems to be okay according to ChatGPT.
        mesh.elements[i].update();
    }

    // Now that the forces on the nodes have been reset, amd
    // the elements updated in parallel, we can sum up the energy
    // and forces.
    for (size_t i = 0; i < mesh.nrElements; i++)
    {
        // this update function updates the forces on the nodes
        mesh.elements[i].applyForcesOnNodes();
        total_energy += mesh.elements[i].energy;
        // std::cout << "El " << i << ": " << mesh.elements[i] << '\n';
    }
    // std::cout << '\n';
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

void run_simulation()
{
    int nrThreads = 8; // Must be between 1 and nr of cpus on machine.

    int s = 3; // Can't be smaller than 3, because with 2, all nodes would be fixed
    int nx = s;
    int ny = s;
    int n = nx * ny;

    // If we have a small number of particles, we need fewer threads
    nrThreads = (nrThreads > s) ? nrThreads : s;


    Mesh mesh = {nx, ny};


    alglib::real_1d_array nodeDisplacements;

    // We set the length of the array to be the same as the number of nodes that
    // are not on the boarder. These are the only nodes we will modify.
    int nrNonBorderNodes = mesh.nonBorderNodeIds.size();
    nodeDisplacements.setlength(2 * nrNonBorderNodes);
    for (int i = 0; i < nrNonBorderNodes; ++i)
    {
        nodeDisplacements[i] = 0;
    }

    // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgssetcond
    double epsg = 0;             // epsilon gradiant. A tolerance for how small the gradiant should be before termination.
    double epsf = 0.0001;        // epsilon function. A tolerance for how small the change in the value of the function between itterations should be before termination.
    double epsx = 0;             // epsilon x-step. A tolerance for how small the step between itterations should be before termination.
    alglib::ae_int_t maxits = 0; // Maximum itterations
    // When all the values above are chosen to be 0, a small epsx is automatically chosen
    // This is unphysical, so for accademic purposes, we should instead use a small epsf value to guarantee a certan level of accuracy.
    alglib::minlbfgsstate state;
    alglib::minlbfgsreport report;

    double load = 0.05;
    double theta = 0;

    // Boundary conditon transformation
    Matrix2x2<double> bcTransform = getShear(load, theta);
    Matrix2x2<double> wrongBcTransform = getShear(load * 1.5, theta + 1);

    write_to_vtu(mesh, "Init1");
    mesh.applyTransformationToBoundary(bcTransform);

    write_to_vtu(mesh, "BC2");

    // Modifies nodeDisplacements
    initialGuess(mesh, wrongBcTransform, nodeDisplacements);

    alglib::minlbfgscreate(1, nodeDisplacements, state);

    // Set termination condition, ei. when is the solution good enough
    // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgssetcond
    alglib::minlbfgssetcond(state, epsg, epsf, epsx, maxits);

    // Temporairaly update the positions of the mesh with the guess so
    // we can see what the guess was.
    updatePossitionOfMesh(mesh, nodeDisplacements);
    write_to_vtu(mesh, "Guess3");
    // Revert back to original position.
    updatePossitionOfMesh(mesh, nodeDisplacements, -1);

    // The null pointer can be replaced with a logging function!
    alglib::minlbfgsoptimize(state, alglib_calc_energy_and_gradiant, nullptr, &mesh);

    // TODO Collecting and analysing these reports could be a usefull tool for optimization
    alglib::minlbfgsresults(state, nodeDisplacements, report);

    printReport(report);

    write_to_vtu(mesh, "Relaxed4");
}

#endif
