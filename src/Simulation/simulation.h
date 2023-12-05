#ifndef SIMULATION_H
#define SIMULATION_H
#pragma once

#include <vector>

#include "settings.h"
#include "Matrix/matrix2x2.h"
#include "Mesh/mesh.h"
#include "Data/dataExport.h"
#include "easylogging++.h"

#include "stdafx.h"
#include "interpolation.h"
#include "specialfunctions.h"
#include "linalg.h"
#include "statistics.h"
#include "alglibmisc.h"

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
// energy from all the cells in the surface.
double calc_energy_and_forces(Mesh &mesh)
{
    // This is the total energy from all the triangles
    double total_energy;

    // TODO Parllelalize this for-loop
    // #pragma omp parallel
    // #pragma omp for reduction(+:energy_thread)
    for (size_t i = 0; i < mesh.nrTriangles; i++)
    {

        // Create a new cell and store it in the cells array of the surface
        // The constructor of the Cell calculates D, C, m and C_ (See tool tip)
        mesh.cells[i] = Cell(mesh.triangles[i]);
        total_energy += mesh.cells[i].energy;
        // Set the forces on the nodes in the cell
        mesh.cells[i].setForcesOnNodes(mesh.triangles[i]);
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
    Mesh *mesh = reinterpret_cast<Mesh*>(meshPtr);

    int nr_x_values = force.length() / 2;

    // Update mesh position
    updatePossitionOfMesh(*mesh, displacement);
    // Calculate energy and forces
    energy = calc_energy_and_forces(*mesh);

    // Update forces
    // TODO TODO (Big TODO) update everything to always use real_1d_array to avoid copying and moving
    // Grad values are stored as follows:
    // given three nodes with positions (x1,y1), (x2,y2) and (x3,y3),
    // force = [x1, x2, x3, y1, y2, y3]. 
    // By finding n_a, as a halfway point, we can correctly assign the values
    for (size_t i = 0; i < nr_x_values; i++)
    {
        NodeId a_id = mesh->nonBorderNodeIds[i];
        force[i] = (*mesh)[a_id]->f_x;
        force[nr_x_values + i] = (*mesh)[a_id]->f_y;
    }
}

// Our initial guess will be that all particles have shifted by the same
// transformation as the border.
void initialGuess(const Mesh &mesh, const BoundaryConditions &bc,
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
        transformed_node = transform(bc.F, *innside_node);     // F * node.position
        // Subtract the initial possition to get the displacement
        translateInPlace(transformed_node, *innside_node, -1); // node1.position - node2.position
        displacement[i] = transformed_node.x;
        displacement[i + nr_x_elements] = transformed_node.y;
    }
}

void printReport(const alglib::minlbfgsreport &report)
{
    // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgsoptimize
    std::cout << "Optimization Report:" << std::endl;
    std::cout << "Iterations Count: " << report.iterationscount << std::endl;
    std::cout << "Number of Function Evaluations: " << report.nfev << std::endl;
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

    int nx = 4;
    int ny = 4;
    int n = nx * ny;

    Mesh mesh = {nx, ny};

    // while(load <1.){

    mesh.resetForceOnNodes();

    alglib::real_1d_array node_displacements;

    // We set the length of the array to be the same as the number of nodes that
    // are not on the boarder. These are the only nodes we will modify.
    node_displacements.setlength(mesh.nonBorderNodeIds.size());
    for (int i = 0; i < node_displacements.length(); ++i)
        node_displacements(i) = 0;

    // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgssetcond
    double epsg = 0;            // epsilon gradiant. A tolerance for how small the gradiant should be before termination.
    double epsf = 0;            // epsilon function. A tolerance for how small the change in the value of the function between itterations should be before termination.
    double epsx = 0;            // epsilon x-step. A tolerance for how small the step between itterations should be before termination.
    alglib::ae_int_t maxits = 0; // Maximum itterations
    // When all the values above are chosen to be 0, a small epsx is automatically chosen
    // This is unphysical, so for accademic purposes, we should instead use a small epsf value to guarantee a certan level of accuracy.
    alglib::minlbfgsstate state;
    alglib::minlbfgsreport report;

    double load = 0.02;
    double theta = 0;
    BoundaryConditions bc = BoundaryConditions{load, theta};
    write_to_a_ovito_file(mesh, "1Init");
    mesh.applyBoundaryConditions(bc);

    write_to_a_ovito_file(mesh, "2BC");

    // Modifies node_displacements
    initialGuess(mesh, bc, node_displacements);
    // Makes the guess slightly wrong
    translate(mesh, mesh.nonBorderNodeIds, -0.02, 0);

    alglib::minlbfgscreate(1, node_displacements, state);

    // Set termination condition, ei. when is the solution good enough
    // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgssetcond
    alglib::minlbfgssetcond(state, epsg, epsf, epsx, maxits);

    write_to_a_ovito_file(mesh, "3Guess");

    // The null pointer can be replaced with a logging function!
    alglib::minlbfgsoptimize(state, alglib_calc_energy_and_gradiant, nullptr, &mesh);

    // TODO Collecting and analysing these reports could be a usefull tool for optimization
    alglib::minlbfgsresults(state, node_displacements, report);

    printReport(report);

    write_to_a_ovito_file(mesh, "4Relaxed");
    /*
    if(singleton.c.linearity == false)
        plas= plasticity_gl2z(singleton.C_.current_metrics,current_metrics_t0);
    if(singleton.c.linearity == true)
        plas= plasticity(singleton.C_.current_metrics,current_metrics_t0);


    if(inner==adaptive.size()-1){

        //save the data just before the avalanche
        if(std::round(plas) >= nx/4){
            //save picture before the avalanche
            save_results(singleton.c,starting_point_keep);
            write_to_vtk(singleton.c,t);
            write_to_a_ovito_file(singleton.c,setnew,t++);
        }

// 					cout<<"inner"<<"--->"<<inner<<endl;

        cout<<load<<"--->"<<plas;

        prediction_failure++;

        if( adaptive.size()==1)
            cout<<"; not failure: "<<prediction_failure<<endl;

        if( adaptive.size()==1){
            if(std::round(plas) > 0 || prediction_failure  >= 10){
                prediction_failure=0;}
            else if(prediction_failure_special  >= 10){
                prediction_failure_special=0;
            }
        }


        if( adaptive.size()!=1){
            cout<<"; failure: "<<prediction_failure<<endl;

        }


        if( adaptive.size()!=1)
            if(std::round(plas) > 0)
                prediction_failure=0;


        break;
    }

// 				if(!plasticity(singleton.C_.current_metrics,current_metrics_t0)){
    if(std::round(plas) <= 0  ){

// 					if(n_load++%10==0){
            cout<<load<<"> "<<plas;
            cout<<"; succes:"<<" "<<prediction_failure<<endl;
            prediction_failure_special = 0;

// 					}


        break;
    }

    //plasticity took place
    std::copy(starting_point_keep.getcontent(),
    starting_point_keep.getcontent()+starting_point_keep.length(),
    node_displacements.getcontent());
    */
}

#endif
