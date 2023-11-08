
#include <vector>

#include "../settings.h"
#include "../Matrix/matrix2x2.h"
#include "../Utility/singelton.h"
#include "../Grid/grid2D.h"
#include "../Energy/energy_and_stress_calculations.h"
#include "../Data/dataExport.h"

#include "stdafx.h"
#include "interpolation.h"
#include "specialfunctions.h"
#include "linalg.h"
#include "statistics.h"
#include "alglibmisc.h"

// Updates the forces on the nodes in the grid and returns the total
// energy from all the cells in the grid.
double calc_energy_and_forces(Grid &g)
{

    // This is the total energy from all the triangles
    double total_energy;

    // TODO Parllelalize this for-loop
    // #pragma omp parallel
    // #pragma omp for reduction(+:energy_thread)
    for (size_t i = 0; i < g.nr_triangles; i++)
    {

        // Create a new cell and store it in the cells array of the grid
        g.cells[i] = Cell(g.triangles[i]);
        // Calculate energy and redused stress (The result is stored in the cell)
        calculate_energy_and_reduced_stress(g.cells[i]);

        total_energy += g.cells[i].energy;

        // Set the forces on the nodes in triangle t.
        g.cells[i].set_forces_on_nodes(g.triangles[i]);
    }
    return total_energy;
}

void alglib_calc_energy_and_gradiant(const alglib::real_1d_array &xy,
                                     double &energy,
                                     alglib::real_1d_array &grad, void *ptr)
{
    Singelton &s = Singelton::getInstance();

    energy = calc_energy_and_forces(s.g);

    int nr_x_values = grad.length() / 2;
    // Grad values are stored as follows:
    // given three nodes with positions (x1,y1), (x2,y2) and (x3,y3),
    // grad = [x1, x2, x3, y1, y2, y3]. (Note that grad stores forces, not possitions. xy stores possitions (I think))
    // By finding n_a, as a halfway point, we can correctly assign the values
    for (size_t i = 0; i < grad.length(); i++)
    {
        node_id a_id = s.g.non_border_node_ids[i];
        grad[i] = s.g[a_id]->f_x;
        grad[nr_x_values + i] = s.g[a_id]->f_y;
    }
}

// Our initial guess will be that all particles have shifted by the same transformation as the
// border.
void initial_guess(const Grid &g, const boundary_conditions &bc,
                   alglib::real_1d_array &displacement)
{

    // The displacement is structed like this: [x1,x2,x3,x4,y1,y2,y3,y4], so we
    // need to know where the "x" end and where the "y" begin.
    int nr_x_elements = displacement.length() / 2; // Shifts to y section

    node transformed_node; // These are temporary variables for readability
    const node *innside_node;

    // We loop over all the nodes that are not on the border, ie. the innside nodes.
    for (size_t i = 0; i < g.non_border_node_ids.size(); i++)
    {
        innside_node = g[g.non_border_node_ids[i]];
        transformed_node = transform(bc.F, *innside_node);                 // F * node.position
        transformed_node = translate(transformed_node, *innside_node, -1); // node1.position - node2.position
        displacement[i] = transformed_node.x;
        displacement[i + nr_x_elements] = transformed_node.y;
    }

}

void run_simulation()
{

    int nx, ny = 5;
    int n = nx * ny;

    Singelton &s = Singelton::getInstance();
    s.setGridSize(nx, ny);
    
    std::cout << "Running";

    // while(load <1.){

    s.g.reset_force_on_nodes();

    alglib::real_1d_array starting_point;
    starting_point.setlength(2 * (nx - 2) * (ny - 2));
    for (int i = 0; i < starting_point.length(); ++i)
        starting_point(i) = 0;

    std::vector<std::vector<Cell>> current_metrics_t0;
    current_metrics_t0.resize(n);
    for (int i = 0; i < n; ++i)
        current_metrics_t0[i].resize(4);

    for (int i = 0; i < starting_point.length(); ++i)
        starting_point(i) = 0;

    // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgssetcond
    double epsg = 0;             // epsilon gradiant. A tolerance for how small the gradiant should be before termination.
    double epsf = 0;             // epsilon function. A tolerance for how small the change in the value of the function between itterations should be before termination.
    double epsx = 0;             // epsilon x-step. A tolerance for how small the step between itterations should be before termination.
    alglib::ae_int_t maxits = 0; // Maximum itterations
    // When all the values above are chosen to be 0, a small epsx is automatically chosen
    // This is unphysical, so for accademic purposes, we should instead use a small epsf value to guarantee a certan level of accuracy.
    alglib::minlbfgsstate state;
    alglib::minlbfgsreport report;

    double load = 0.01;
    double theta = 0;
    boundary_conditions bc = boundary_conditions{load, theta};
    s.g.apply_boundary_conditions(bc);

    //initial_guess(s.g, bc, starting_point);

    alglib::minlbfgscreate(10, starting_point, state);
    // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgssetcond
    alglib::minlbfgssetcond(state, epsg, epsf, epsx, maxits);

    alglib::minlbfgsoptimize(state, alglib_calc_energy_and_gradiant);

    alglib::minlbfgsresults(state, starting_point, report);
    write_to_a_ovito_file(s.g);
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
    starting_point.getcontent());
    */
}
