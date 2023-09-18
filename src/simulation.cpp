
#include <vector>
#include "matrix2x2.h"
#include "singelton.h"
#include "grid2D.h"

#include "stdafx.h"
#include "interpolation.h"
#include "specialfunctions.h"
#include "linalg.h"
#include "statistics.h"
#include "alglibmisc.h"

struct basis
{
    // This is the basis of the grid
	double e1[2];
	double e2[2];
};

struct cell{
    Matrix2x2<double> data;
    bool plasticity; // not used yet
};



void calc_energy_forces(struct Grid& g,  alglib::real_1d_array &grad){
	
	struct cell metrics,metrics_reduced,metrics_elastic;
    
    //#pragma omp parallel
	//#pragma omp for reduction(+:energy_thread)
	for (const triangle t : g.triangles) {
		
				
//---------------------------------METRICS------------------------------------		

			//linear or not
			metrics.data  = t.metric(g.a, MetricFunction::faicella); //c.fptr3(v);
  	 	    metrics_reduced =  lagrange_reduction(metrics);
//   	 	    metrics_elastic =  elastic_reduction(metrics);
  	 	    
//   	 	    cout<<"metrics_reduced"<<endl;
//   	 	    metrics_reduced.m.print();
//   	 	    
//   	 	    cout<<"metrics_elastic"<<endl;
//   	 	    metrics_elastic.m.print();

//   	 	    c.current_metrics[i][k].m = metrics_reduced.m;
// 			if(metrics_reduced.plasticity == true)
// 			 c.plasticity = true;//metrics_reduced.plasticity;


//  			temp = coarse_grain_fast(metrics_reduced);
			temp = c.fptr(metrics_reduced)  ;
			energy_thread += temp.energy ;//-  c.disorder_x[i]*metrics_reduced.c12 ;

// 	 		stress = c.fptr2(metrics_reduced);
	 		stress = temp.r_stress;
	 		temp.r_stress.c12 = temp.r_stress.c12 ;//-  c.disorder_x[i] ;

			d =  c.fptr4(stress, v, metrics_reduced.m);
//---------------------------------FORCES ON  NODES------------------------------------		


			c.grp1[i][l].x   = a*d.gr[0].x ; 
			c.grp1[i][l].y   = a*d.gr[0].y ; 

			c.grp1[i][l+1].x = a*d.gr[1].x ; 
			c.grp1[i][l+1].y = a*d.gr[1].y ; 

			c.grp1[i][l+2].x = a*d.gr[2].x ; 	
			c.grp1[i][l+2].y = a*d.gr[2].y ; 


	 	}
 
	for (const auto& element : c.myVector) {
		 int i = element.first;
		 int k = element.second;
		 int l=3*k;
	
			struct punto_stru p=c.p[i];
			c.gr[i].x   += c.grp1[i][l].x;
			c.gr[i].y   += c.grp1[i][l].y;

			c.gr[p.nn[k]].x   += c.grp1[i][l+1].x;
			c.gr[p.nn[k]].y   += c.grp1[i][l+1].y;
			c.gr[p.nn[k+1]].x += c.grp1[i][l+2].x;
			c.gr[p.nn[k+1]].y += c.grp1[i][l+2].y;  
    	 
    }
    
	c.energy = energy_thread;
}
void function1_foralglib2D_DOF_version(const alglib::real_1d_array &x, double &func, alglib::real_1d_array &grad, void *ptr) 
{
    S& s = S::getInstance();
	calc_energy_forces(s.g, grad);

}
void run_simulation(){

    int nx, ny = 5;
    int n=nx*ny;

    S& s = S::getInstance();
    s.setGridSize(nx, ny);

    //while(load <1.){			
    
    
        //choose_bc(setnew,bc);
        //apply_BC_body(singleton.c,setnew);
        //apply_BC_CG(singleton.c);

        // initial_guess for the displacement
        //initial_guess( singleton.c, set_incr,  starting_point);

        alglib::real_1d_array starting_point;
		starting_point.setlength(2*(nx-2)*(ny-2));
		
		std::vector <std::vector<cell>> current_metrics_t0;
		current_metrics_t0.resize(n);
		for (int i = 0; i < n; ++i)
                current_metrics_t0[i].resize(4);

		for (int i=0;i<starting_point.length();++i)
			starting_point(i)=0;

        double epsg = 0;
        double epsf = 0;
        double epsx = 0;
        alglib::ae_int_t maxits = 0;
        alglib::minlbfgsstate state;
        alglib::minlbfgsreport rep;



        alglib::minlbfgscreate(10, starting_point, state);
        alglib::minlbfgssetcond(state, epsg, epsf, epsx, maxits);
 
        alglib::minlbfgsoptimize(state, function1_foralglib2D_DOF_version);

        alglib::minlbfgsresults(state, starting_point, rep);
        
        save_results(singleton.c,starting_point);
        if(singleton.c.linearity == false)
            plas= plasticity_gl2z(singleton.c.current_metrics,current_metrics_t0);
        if(singleton.c.linearity == true)
            plas= plasticity(singleton.c.current_metrics,current_metrics_t0);


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
        
// 				if(!plasticity(singleton.c.current_metrics,current_metrics_t0)){
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

    }

    
    
// 			std::cout<<"---"<<std::endl;
    //this copies the output of ALGLIB into the singletion c.p = y = x+u = c.pfix + starting_point
// 			transfer_from_alglib_to_std_v2(singleton.c,starting_point);

    std::copy(starting_point.getcontent(),
    starting_point.getcontent()+starting_point.length(), 
    starting_point_keep.getcontent());

    save_results(singleton.c,starting_point);


    save_elastic_state_gl2z(singleton.c.current_metrics,current_metrics_t0);



    // Copy the ALGLIB array to a std::vector
//     		std::vector<double> myVector(starting_point.getcontent(), starting_point.getcontent() + starting_point.length());
// 			std::cout<<"maxval displacement "<<maxval(myVector)<<std::endl;
// 			std::cout<<"minval displacement "<<minval(myVector)<<std::endl;
// 			std::cout<<"maxloc displacement "<<maxloc(myVector)<<std::endl;
// 			std::cout<<"minloc displacement "<<minloc(myVector)<<std::endl;


    if(std::round(plas) >= nx/4){
        write_to_vtk(singleton.c,t);
        write_to_a_ovito_file(singleton.c,setnew,t++);
    }
    
    if(prediction_failure<1){
        adaptive.resize(2);
        adaptive[0] = 2*load_increment;
        adaptive[1] = load_increment;

    }

    else if(prediction_failure>=1){
        adaptive.resize(1);
        adaptive[0] = load_increment;		
    }


// 	

    
    struct average_stress sum_stress = sum_2d_cella(singleton.c.stress,setnew);
    double sum_energy = sum_2d_vector(singleton.c.energy_local);

// 				cout<<"after: "<<endl;
// 
// 				cout<<"12: "<<sum_stress.stress12<<endl;
// 				cout<<"11: "<<sum_stress.stress11<<endl;
// 				cout<<"22: "<<sum_stress.stress22<<endl;
// 				cout<<"to: "<<sum_stress.stress_total<<endl;
// 
// 				cout<<"ee: "<<sum_energy<<endl;

    
    write_data_to_files(load,sum_energy,sum_stress,sum_p_energy,sum_p_stress,plas);
    


    
    //update the nodes on the boundary with the new load
// 			std::cout<<" load at "<<load<<std::endl;


}
}