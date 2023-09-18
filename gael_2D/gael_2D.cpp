#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <algorithm>    // std::min
#include <sstream>
#include <fstream>
#include <iostream>
#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <array>
#include <iomanip>

#include <vector>
#include <numeric>
#include <random>

#include <random>
//#include "Random123/philox.h"
//#include "Random123/boxmuller.hpp"

#include <iostream>
#include <sys/stat.h>
#include "optimization.h" //ALGLIB
//written by umut
// #include "common.h"
#include "common_2D.h"
#include "utility_functions.h"
#include "boundary_conditions.h"


#include "stdafx.h"
#include "interpolation.h"
#include "specialfunctions.h"
#include "linalg.h"
#include "statistics.h"
#include "alglibmisc.h"

using namespace alglib;

using namespace std;

bool createFolder(std::string path)
{
    // Check if folder already exists
    struct stat info;
    if (stat(path.c_str(), &info) == 0 && S_ISDIR(info.st_mode))
    {
        std::cout << "Folder already exists" << std::endl;
        return true;
    }

    // Create the folder
    int status = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (status != 0)
    {
        std::cerr << "Error creating folder" << std::endl;
        return false;
    }

    std::cout << "Folder created successfully" << std::endl;
    return true;
}

double maxval(std::vector<double>& m)
{
    double max;
    max  = *max_element(m.begin(), m.end());
    return max;

}

double minval(std::vector<double>& m)
{
    double min;
    min  = *min_element(m.begin(), m.end());
    return min;

}



int maxloc(std::vector<double>& m)
{
    int max;
    max  = max_element(m.begin(), m.end()) - m.begin();
    return max;

}

int minloc(std::vector<double>& m)
{
    int min;
    min  = min_element(m.begin(), m.end()) - m.begin();
    return min;

}




double sum_vector(const std::vector<double>& v) {
  double sum = std::accumulate(v.begin(), v.end(), 0.0);
  return sum;
}

double sum_2d_vector(const std::vector<std::vector<double>>& v) {
  double sum = 0.0;
    	fem_2D::calculation& singleton =  fem_2D::calculation::getInstance();

  	for (const auto& element : singleton.c.myVector){
     		int i = element.first;
			int k = element.second;

    sum += singleton.c.energy_local[i][k];
  }
  return sum;
}

double contraction (const struct cella_stru& stress, const struct boundary_conditions& setnew){

	struct base_stru n;
	struct base_stru a;
	struct base_stru co;

	co.e1[0] = 1.; co.e1[1] = 0.;
	co.e2[0] = 0.; co.e2[1] = 1.;


	// GLIDING PLANE
	a.e1[0] = co.e1[0] ;
	a.e1[1] = co.e1[1] ;
	//NORMAL TO  GLIDING PLANE 
	n.e1[0] = co.e2[0];
	n.e1[1] = co.e2[1];


	alglib::real_1d_array ad, nd; 
	ad.setlength(2);
	nd.setlength(2);

	ad(0) = a.e1[0];
	ad(1) = a.e1[1];
	nd(0) = n.e1[0];
	nd(1) = n.e1[1];


	alglib::real_2d_array shear;
	shear.setlength(2,2);

	for(int i = 0; i < ad.length(); i++) {
 	   for(int j = 0; j < nd.length(); j++) {
      	  shear(i, j) = ad(i) * nd(j);
    	}
	}	
	
	return  shear(0,0)*stress.c11+
			shear(1,1)*stress.c22+
			shear(0,1)*stress.c12+
			shear(1,0)*stress.c21;

// 	return  stress.c11*setnew.f.f11+
// 			stress.c22*setnew.f.f22+
// 			stress.c12*setnew.f.f12/2.+
// 			stress.c12*setnew.f.f21/2.;
}



struct average_stress sum_2d_cella(const std::vector<std::vector<cella_stru>>& v,const struct boundary_conditions& setnew) {
  
  	double sum = 0.0;
  	struct average_stress temp;
  
  	fem_2D::calculation& singleton =  fem_2D::calculation::getInstance();

  	for (const auto& element : singleton.c.myVector){
  	
  		int i = element.first;
		int k = element.second;

  		temp.stress12 += v[i][k].c12;
  		temp.stress21 += v[i][k].c21;

  		temp.stress11 += v[i][k].c11;
  		temp.stress22 += v[i][k].c22;
  		temp.stress_total += contraction(v[i][k],setnew) ;

    }

  return temp;
}


string IntToStr(int n){
	stringstream result;
	result << 10000 + n;
	return result.str();
}

// std::vector<double> brownian(double av, double std, int n, int& seed){
//         typedef r123::Philox2x32 RNG;
//         RNG rng;
//         RNG::ctr_type c={{}};
//         RNG::ukey_type uk={{}};
//         uk[0] = seed; // some user_supplied_seed
//         RNG::key_type k=uk;
//         std::vector<double> r;
//         r.resize(n);

//         for (int i = 0; i<n; ++i){
//                 c[0] = i; // some loop-dependent application variable 
//                 c[1] = i+n; // another loop-dependent application variable 
//                 RNG::ctr_type rnd_numb = rng(c, k);
//                 // use the random values in r for some operation related to
//                 // this iteration on objectid
//                 r[i] = std*r123::boxmuller(rnd_numb.v[0], rnd_numb.v[1]).x + av;
//         };
//         return r;
// };


std::vector<double>  wiener(double av, double std,int n, int seed)
{


	std::vector<double> r;
	r.resize(n);
	static std::random_device rd;
	static std::default_random_engine generator;
// 	if (i == 0){
// 		generator.seed(12);   //Now this is seeded same each time.
// 		i++;
// 	}
		generator.seed( seed); //Now this is seeded differently each time.

// 		generator.seed( 12 );   //Now this is seeded differently each time.



	std::normal_distribution<double> dis(av, std);
// 	std::uniform_real_distribution<double> dis(-0.001, 0.001);




	for (int i = 0; i<n; ++i)
		r[i] = dis(generator);


    

// 		cout<<"max r= "<<maxval(r)<<endl;
// 		cout<<"min r= "<<minval(r)<<endl;

	return r;

}

double square_energy(double r){

double E=2.71828182846;



double a  = 1.;
double c1 = 2.*a;
double c2 = 2.*a;
double b1 = 8.;
double b2 = 8.;

double r2=1.41421;

double r1=1.;


//Molecular-dynamics study of the melting of hexagonal and square lattices in two dimensions
// 	double energy=-2./pow(E,8*pow(-1.425 + r,2)) - 2./pow(E,8*pow(-1 + r,2)) + pow(r,-12);

double energy = -(c1/pow(E,b1*pow(r - r1,2))) - c2/pow(E,b1*pow(r - r2,2)) + a/pow(r,12);
return energy;

}

//dphi/dr
double square_energy_der(double r){

	double E=2.71828182846;
    double a  = 1.;
    double c1 = 2.*a;
    double c2 = 2.*a;
    double b1 = 8.;
    double b2 = 8.;
  
	double r2=1.41421;

	double r1=1.;


// double der = (-43.884*(-1.79 + r))/pow(E,26.5*pow(-1.79 + r,2)) - 12/pow(r,13) + 12/pow(r,7);

//Molecular-dynamics study of the melting of hexagonal and square lattices in two dimensions
//    double der = (32*(-1.425 + r))/pow(E,8*pow(-1.425 + r,2)) + 
//    (32*(-1 + r))/pow(E,8*pow(-1 + r,2)) - 12/pow(r,13);

 double der =  (-12*a)/pow(r,13) + (2*b1*c1*(r - r1))/pow(E,b1*pow(r - r1,2)) + 
   (2*b1*c2*(r - r2))/pow(E,b1*pow(r - r2,2));

	return der;

}

void atomistic_grid_square(struct conf_stru& c){
	double xx,yy;

	std::vector<double> x(0),y(0);

    double specialx =1.06619;		
	
	double specialy = specialx;
	
	for(int j=0;j<8;j++){

		for(int i=0;i<8;i++){	

			xx = specialx*(double)i;
			yy = specialx*(double)j;

			x.push_back(xx);// = xx; 
			y.push_back(yy);// = yy;
		}
	}


	double length0_x = *max_element(x.begin(), x.end()) +  specialx;
	double length0_y = *max_element(y.begin(), y.end()) +  specialy;
	
	c.length0_x=length0_x;
	c.length0_y=length0_y;


	cout<<"x_length = "<<c.length0_x<<"\n";
	cout<<"y_length = "<<c.length0_y<<"\n";

 	cout<<"VERBOSE:: NUMBER OF ATOMS in SUB-GRID"<<x.size()<<endl;
 	
 	
 	for (int n=0;n<x.size();n++){
		c.atom_x_co.push_back((1./length0_x)*x[n]) ;
		c.atom_y_co.push_back((1./length0_y)*y[n]) ;
	}
}


struct energy_stress coarse_grain_fast(const struct cella_stru& metric){

	fem_2D::calculation& singleton =  fem_2D::calculation::getInstance();

	int nx=8;
	int ny=8;
	int n = nx*ny;
	int counter=0;

	double e0=0;
	double f0=0;
	double sc11=0;
	double sc22=0;
	double sc12=0;
	double r_cutoff_sq = 6.5*6.5;
	
	struct energy_stress temp;
	double c11= metric.c11;
	double c22= metric.c22;
	double c12= metric.c12;
	double sqrtdetc=sqrt(c11*c22-c12*c12);
	

		
	for(int i=ny*nx/2+nx/2 -1;i<ny*nx/2+nx/2 ;i++){

		for(int j=0;j<n;j++){
			
			if(i==j)
				continue;
		
			double dx_tilde = (singleton.c.atom_x_co[i]  -singleton.c.atom_x_co[j]);
			double dy_tilde = (singleton.c.atom_y_co[i]  -singleton.c.atom_y_co[j]);



			if(abs(dx_tilde) > 0.5) 
				dx_tilde -=   copysign(1., dx_tilde) ;
			if(abs(dy_tilde) > 0.5) 
				dy_tilde -=   copysign(1., dy_tilde) ;

			double dx2 = singleton.c.length0_x*dx_tilde ;
			double dy2 = singleton.c.length0_y*dy_tilde ;


			double r2 =  dx2*dx2*metric.c11 + dy2*dy2*metric.c22 + 2.*dx2*dy2*metric.c12 ;
			
			if( r2 >= r_cutoff_sq )
				continue;
		
			double r = sqrt(r2);
			
			e0 +=  0.5*singleton.c.atom_energy(r);
			double tempf = 0.5*singleton.c.atom_stress(r)/r;
			
			//dphi/dC
			sc11 +=   0.5*tempf*dx2*dx2;
			sc22 +=   0.5*tempf*dy2*dy2;
			sc12 +=	  (tempf*dx2*dy2);

// 				cout<<"energy"<<e0<<endl;

		}
		
	}
		

		temp.energy  =  e0 - singleton.c.zeroing;
		temp.r_stress.c11 = sc11;
		temp.r_stress.c22 = sc22;
		temp.r_stress.c12 = sc12;

		
		return temp;
					
	

}

struct cella_stru lagrange_reduction(struct cella_stru c){

		fem_2D::calculation& singleton =  fem_2D::calculation::getInstance();
		
		if(singleton.c.linearity == true){
		
			c.m = c.m.identity();
			c.c12  = c.c12;
			return c;
		}
		
		struct matrix_stru lag_m1,lag_m2,lag_m3;
		
		lag_m1.m11 = 1.; lag_m1.m22 = -1.;  lag_m1.m12 = 0. ; lag_m1.m21=0.; 
		//np.array([[1., 0.],[0.,-1.]])
		lag_m2.m11 = 0.; lag_m2.m22 = 0. ;  lag_m2.m12 = 1. ; lag_m2.m21=1.;  
		//np.array([[0., 1.],[1., 0.]])
		lag_m3.m11 = 1.; lag_m3.m22 = 1. ;  lag_m3.m12 = -1.; lag_m3.m21=0.;
		//np.array([[1.,-1.],[0., 1.]])
		
		//set to identity the m matrix			
		c.m = c.m.identity();

	while(c.c12<0 || c.c22<c.c11 || 2 * c.c12>c.c11 ){
		
		if (c.c12<0){
			c.c12 = -c.c12;
			c.m = c.m.multiply(lag_m1);

		}
		
		if (c.c22<c.c11){
			double a;
			a = c.c11;
			c.c11 = c.c22;
			c.c22 = a;
			c.m = c.m.multiply(lag_m2);
		}

		if (2 * c.c12>c.c11){
			struct cella_stru d;
			d.c11 = c.c11;
			d.c12 = c.c12 - c.c11;
			d.c22 = c.c22 + c.c11 - 2 * c.c12;
			c.m = c.m.multiply(lag_m3);
			c.c11 = d.c11; 
			c.c12 = d.c12; 
			c.c22 = d.c22; 
		}
		

		
	}
    
	return c;
}

struct cella_stru elastic_reduction(struct cella_stru c){

		
		
		fem_2D::calculation& singleton =  fem_2D::calculation::getInstance();
		
		if(singleton.c.linearity == true){
		
			c.m = c.m.identity();
			c.m_el = c.m.identity();
			c.c11 = 1.;
			c.c22 = 1.;
			c.c12  = c.c12;
	
			return c;
		}

		
		static struct matrix_stru lag_m1,lag_m2,lag_m3, m_cap;
		
		struct matrix_stru m1,m2,m3,m_el;
	
// 	 	fem_2D::calculation& singleton =  fem_2D::calculation::getInstance();

		
		lag_m1.m11 = 1.; lag_m1.m22 = -1.;  lag_m1.m12 = 0. ; lag_m1.m21=0.; 
		//np.array([[1., 0.],[0.,-1.]])
		lag_m2.m11 = 0.; lag_m2.m22 = 0. ;  lag_m2.m12 = 1. ; lag_m2.m21=1.;  
		//np.array([[0., 1.],[1., 0.]])
		lag_m3.m11 = 1.; lag_m3.m22 = 1. ;  lag_m3.m12 = -1.; lag_m3.m21=0.;
		//np.array([[1.,-1.],[0., 1.]])
		
		//set to identity the m matrix
// 		c.m    = c.m.identity();
		c.m_el = c.m_el.identity();


		
		struct cella_stru metric_el;
		struct cella_stru metric_o;
		
		

		metric_el.m.m11 = metric_o.m.m11 = c.c11;
		metric_el.m.m22 = metric_o.m.m22 = c.c22;
		metric_el.m.m12 = metric_o.m.m12 = c.c12;
		metric_el.m.m21 = metric_o.m.m21 = c.c12;
		
		int i =0; 

	while(c.c12<0 || c.c22<c.c11 || 2 * c.c12>c.c11 ){
	
			m1  = m1.identity();
			m2  = m2.identity();
			m3  = m3.identity();
		    m_el =m_el.identity();

		
		if (c.c12<0){
			c.c12 = -c.c12;
// 			c.m = c.m.multiply(lag_m1);
			m1=lag_m1;

		}
		
		if (c.c22<c.c11){
			double a;
			a = c.c11;
			c.c11 = c.c22;
			c.c22 = a;
// 			c.m = c.m.multiply(lag_m2);
			m2=lag_m2;
		}

		if (2 * c.c12>c.c11){
			struct cella_stru d;
			d.c11 = c.c11;
			d.c12 = c.c12 - c.c11;
			d.c22 = c.c22 + c.c11 - 2 * c.c12;
// 			c.m = c.m.multiply(lag_m3);
			c.c11 = d.c11; 
			c.c12 = d.c12; 
			c.c22 = d.c22; 
			m3=lag_m3;
		}
// 	m_cap =   np.matmul(m1,m2)
// 	m_cap =   np.matmul(m_cap,m3)
// 	m_cap =   np.matmul(m_cap,m2)
// 	m_cap =   np.matmul(m_cap,m1)


			m_el =   m1.multiply(m2).multiply(m3).multiply(m2).multiply(m1)    ;// np.matmul(m1,m2)
			c.m_el =  c.m_el.multiply(m_el);
			
			
			metric_el.m=    c.m_el.transpose().multiply(metric_o.m).multiply(c.m_el);
// 			cout<<"counter: "<<i++<<endl;
// 			metric_el.m.print();
// 			metric_el.m=    metric_el.m.multiply(c.m_el);  //  np.matmul(np.matmul(np.transpose(m_cap),metrico),m_cap)
			

		
	}
	
    
			c.c11 = metric_el.m.m11;
			c.c22 = metric_el.m.m22;
			c.c12 = metric_el.m.m12;

//  		  if( singleton.c.load>= 0.235){
// 			  cout<<"t"<<endl;
// 			  cout<<c.m_el.m11<<endl;
// 			  cout<<c.m_el.m22<<endl;
// 			  cout<<c.m_el.m12<<endl;
// 			  cout<<c.m_el.m21<<endl;

// 		  }   
	return c;
}





struct energy_stress boyer_energy(const cella_stru&  metrics){

 	fem_2D::calculation& singleton =  fem_2D::calculation::getInstance();

 
    
    double phi0=0;
    double dphi=0;

    
    static double scl  = 1.0661;
	
	double c11 = metrics.c11; 
	double c22 = metrics.c22; 
	double c12 = metrics.c12;
	double sc11=0;
	double sc22=0;
	double sc12=0;
	
	struct energy_stress temp;
	
	static double scl_sq = scl*scl;
    
//for (int s = -limit; s <= limit; s++) {
//         for (int l = -limit; l <= limit; l++) {
    
    	for (const auto& element : singleton.c.idx_coarse_grain){
		  	
 		  	int s = element.first;
 		  	int l = element.second;
  
            if(s==0 && l==0)
            	continue;
            
			double r = scl * std::sqrt(s*s*c11 + 2*s*l*c12 + l*l*c22);
			double tmp=singleton.c.atom_energy(r);

			phi0 += tmp;
		 
			double der = singleton.c.atom_stress(r);
			
			dphi = 0.5*der/r;

			sc11 += 0.5*dphi*s*s *scl_sq;
			sc22 += 0.5*dphi*l*l *scl_sq;
			sc12 += dphi*s*l*scl_sq;
	
				
            
        }
    
    //}
    phi0 *= 0.5;
  
	temp.energy  =  phi0 - singleton.c.zeroing;
	temp.r_stress.c11 = sc11;
	temp.r_stress.c22 = sc22;
	temp.r_stress.c12 = sc12;

    
    return temp;

}

struct cella_stru boyer_denergy(const cella_stru&  metrics){

    
        
    double r1 = 1.;
    double r2 = 1.425;
    double b1 = 8.;
    double b2 = 8.;
    double a  = 1.;
    double c1 = 2.*a;
    double c2 = 0.;
    
//     if(lattice == "triangular") {
//         c2 = 0.;
//     }
//     if(lattice == "square") {
        c2 = c1;
//     }
    
    double scl  = 1.0661;
    double cut  = 2.5;
    double dphi = 0.;
    double E = 2.71828182846;
    double tmp_d = 0.;
    


	double c11 = metrics.c11; 
	double c22 = metrics.c22; 
	double c12 = metrics.c12;


    struct cella_stru  stresses;
stresses.c11=0;
stresses.c22=0;
stresses.c12=0;

int limit=5;
    
    for (int s = -limit; s < limit; s++) {
        for (int l = -limit; l < limit; l++) {
            
            if ((s != 0) || (l != 0)) {
                
                double r = scl * std::sqrt((s*s)*c11 + 2.*s*l*c12 + (l*l)*c22);
                

                
//                 tmp_d  = (32*(-1.425 + r))/std::pow(E,8*std::pow(-1.425 + r, 2));
//                 tmp_d += (32*(-1 + r))/std::pow(E,8*std::pow(-1 + r, 2)) - 12/std::pow(r, 13);
  			   double der = (32*(-1.425 + r))/pow(E,8*pow(-1.425 + r,2)) + 
  				 (32*(-1 + r))/pow(E,8*pow(-1 + r,2)) - 12/pow(r,13);

                dphi += 0.5*der/r;

                stresses.c11 += 0.5*dphi*s*s;
                stresses.c22 += 0.5*dphi*l*l;
                stresses.c12 += dphi*s*l;
            }
        }
    }
     
    
    return stresses;

}

struct energy_stress conti_energy(const cella_stru&  metrics){
 

	double energy = 0.;
	double burgers   = 1.;
	double beta = -0.25;
	double K=4.;
// 	static double aa= 1/pow(2.,1./6.);
// 
// 	static double aa12 = pow(aa,12);
// 	static double aa6 = pow(aa,6);

// 	double theta = 4.;
// 	double gamma = 10.;


	double c11 = metrics.c11;
	double c22 = metrics.c22;
	double c12 = metrics.c12;

    struct energy_stress temp;
	fem_2D::calculation& singleton =  fem_2D::calculation::getInstance();


	//log detC 

   energy  =-K*(log((c11*c22-c12*c12)/burgers)-(c11*c22-c12*c12)/burgers)+beta*(pow(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22),2.0)*9.46969696969697E-4-pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),3.0)*(4.1E1/9.9E1)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,4.0)*(1.0/8.1E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(7.0/1.98E2))+pow(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22),2.0)*(1.7E1/5.28E2)+pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),3.0)*(4.0/1.1E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,3.0)*(1.0/2.7E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(8.0/3.3E1);
//    energy  -= -K*(log((1.)/burgers)-(1.)/burgers);



	temp.energy  =  energy - singleton.c.zeroing;
	temp.r_stress.c11  = -beta*(pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(4.1E1/3.3E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.0/5.28E2)-((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,3.0)*(4.0/8.1E1)-1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,4.0)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(1.0/8.1E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(7.0/1.98E2)+c22*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0)*pow(c11-c12+c22,4.0)*(2.0/8.1E1)-((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(7.0/1.98E2)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(7.0/1.98E2)-c22*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(7.0/3.96E2))+K*(c22/burgers-c22/(c11*c22-c12*c12))+pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(1.2E1/1.1E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.7E1/2.64E2)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,3.0)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.0/2.7E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,2.0)*(1.0/9.0)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(8.0/3.3E1)+c22*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12+c22,3.0)*(1.0/1.8E1)-((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(8.0/3.3E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(8.0/3.3E1)-c22*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(4.0/3.3E1);

	temp.r_stress.c22  =   beta*(pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(4.1E1/3.3E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.0/5.28E2)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,3.0)*(4.0/8.1E1)-1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,4.0)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(1.0/8.1E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(7.0/1.98E2)-c11*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0)*pow(c11-c12+c22,4.0)*(2.0/8.1E1)-((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(7.0/1.98E2)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(7.0/1.98E2)+c11*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(7.0/3.96E2))+K*(c11/burgers-c11/(c11*c22-c12*c12))-pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(1.2E1/1.1E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.7E1/2.64E2)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,3.0)*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.0/2.7E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,2.0)*(1.0/9.0)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(8.0/3.3E1)+c11*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12+c22,3.0)*(1.0/1.8E1)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(8.0/3.3E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(8.0/3.3E1)-c11*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(4.0/3.3E1);

	temp.r_stress.c12  = pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*(1.2E1/1.1E1)-K*((c12*2.0)/burgers-(c12*2.0)/(c11*c22-c12*c12))+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(1.7E1/2.64E2)-beta*(pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*(4.1E1/3.3E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(1.0/5.28E2)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,3.0)*(4.0/8.1E1)-1.0/pow(c11*c22-c12*c12,2.0)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*pow(c11-c12+c22,4.0)*(1.0/8.1E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(7.0/1.98E2)-c12*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0)*pow(c11-c12+c22,4.0)*(4.0/8.1E1)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(7.0/1.98E2)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*(c11-c12+c22)*(7.0/1.98E2)+c12*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(7.0/1.98E2))-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,3.0)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(1.0/2.7E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,2.0)*(1.0/9.0)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(8.0/3.3E1)-c12*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12+c22,3.0)*(1.0/9.0)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(8.0/3.3E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*(c11-c12+c22)*(8.0/3.3E1)+c12*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(8.0/3.3E1);





   return temp;
}


struct energy_stress linear_energy(const cella_stru&  metrics){
 

	double energy = 0.;
	double K=1.6;
// 	static double aa= 1/pow(2.,1./6.);
// 
// 	static double aa12 = pow(aa,12);
// 	static double aa6 = pow(aa,6);

// 	double theta = 4.;
// 	double gamma = 10.;


	double c11 = metrics.c11;
	double c22 = metrics.c22;
	double c12 = metrics.c12;

    struct energy_stress temp;
	fem_2D::calculation& singleton =  fem_2D::calculation::getInstance();
	double Pi = alglib :: pi ();

	//log detC 
	energy  = (1.-cos(2*Pi *c12 )) / pow(2*Pi,2.)  + K*c11*c11/2.;



	temp.energy  =  energy ;
	temp.r_stress.c12 =  sin(2*Pi*c12)/(2.*Pi) ; 
	temp.r_stress.c22 =  0.; 
	temp.r_stress.c11 =  K*c11; 



   return temp;
}



struct cella_stru conti_stress(const cella_stru&  metrics){


	double energy = 0.;
	double burgers   = 1.;
	double beta = -0.25;
	double K=4.;
	double theta = 4.;
	double gamma = 10.;

	double c11 = metrics.c11;
	double c22 = metrics.c22;
	double c12 = metrics.c12;


    struct cella_stru  stresses;


  stresses.c11 = -beta*(pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(4.1E1/3.3E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.0/5.28E2)-((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,3.0)*(4.0/8.1E1)-1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,4.0)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(1.0/8.1E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(7.0/1.98E2)+c22*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0)*pow(c11-c12+c22,4.0)*(2.0/8.1E1)-((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(7.0/1.98E2)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(7.0/1.98E2)-c22*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(7.0/3.96E2))+K*(c22/burgers-c22/(c11*c22-c12*c12))+pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(1.2E1/1.1E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.7E1/2.64E2)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,3.0)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.0/2.7E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,2.0)*(1.0/9.0)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(8.0/3.3E1)+c22*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12+c22,3.0)*(1.0/1.8E1)-((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(8.0/3.3E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(8.0/3.3E1)-c22*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(4.0/3.3E1);
  stresses.c22 =  beta*(pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(4.1E1/3.3E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.0/5.28E2)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,3.0)*(4.0/8.1E1)-1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,4.0)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(1.0/8.1E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(7.0/1.98E2)-c11*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0)*pow(c11-c12+c22,4.0)*(2.0/8.1E1)-((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(7.0/1.98E2)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(7.0/1.98E2)+c11*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(7.0/3.96E2))+K*(c11/burgers-c11/(c11*c22-c12*c12))-pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(1.2E1/1.1E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.7E1/2.64E2)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,3.0)*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.0/2.7E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,2.0)*(1.0/9.0)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(8.0/3.3E1)+c11*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12+c22,3.0)*(1.0/1.8E1)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(8.0/3.3E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(8.0/3.3E1)-c11*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(4.0/3.3E1);
  stresses.c12 = pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*(1.2E1/1.1E1)-K*((c12*2.0)/burgers-(c12*2.0)/(c11*c22-c12*c12))+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(1.7E1/2.64E2)-beta*(pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*(4.1E1/3.3E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(1.0/5.28E2)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,3.0)*(4.0/8.1E1)-1.0/pow(c11*c22-c12*c12,2.0)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*pow(c11-c12+c22,4.0)*(1.0/8.1E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(7.0/1.98E2)-c12*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0)*pow(c11-c12+c22,4.0)*(4.0/8.1E1)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(7.0/1.98E2)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*(c11-c12+c22)*(7.0/1.98E2)+c12*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(7.0/1.98E2))-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,3.0)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(1.0/2.7E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,2.0)*(1.0/9.0)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(8.0/3.3E1)-c12*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12+c22,3.0)*(1.0/9.0)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(8.0/3.3E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*(c11-c12+c22)*(8.0/3.3E1)+c12*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(8.0/3.3E1);
//   stresses.c11 = -beta*(pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(4.1E1/3.3E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.0/5.28E2)-((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,3.0)*(4.0/8.1E1)-1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,4.0)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(1.0/8.1E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(7.0/1.98E2)+c22*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0)*pow(c11-c12+c22,4.0)*(2.0/8.1E1)-((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(7.0/1.98E2)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(7.0/1.98E2)-c22*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(7.0/3.96E2))+pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(1.2E1/1.1E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.7E1/2.64E2)-c22*gamma*pow(burgers-c11*c22+c12*c12,3.0)*4.0+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,3.0)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.0/2.7E1)-c22*theta*(burgers-c11*c22+c12*c12)*2.0-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,2.0)*(1.0/9.0)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(8.0/3.3E1)+c22*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12+c22,3.0)*(1.0/1.8E1)-((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(8.0/3.3E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(8.0/3.3E1)-c22*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(4.0/3.3E1);
//   stresses.c22 = beta*(pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(4.1E1/3.3E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.0/5.28E2)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,3.0)*(4.0/8.1E1)-1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,4.0)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(1.0/8.1E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(7.0/1.98E2)-c11*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0)*pow(c11-c12+c22,4.0)*(2.0/8.1E1)-((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(7.0/1.98E2)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(7.0/1.98E2)+c11*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(7.0/3.96E2))-pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(1.2E1/1.1E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.7E1/2.64E2)-c11*gamma*pow(burgers-c11*c22+c12*c12,3.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,3.0)*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.0/2.7E1)-c11*theta*(burgers-c11*c22+c12*c12)*2.0-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,2.0)*(1.0/9.0)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(8.0/3.3E1)+c11*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12+c22,3.0)*(1.0/1.8E1)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(8.0/3.3E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(8.0/3.3E1)-c11*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(4.0/3.3E1);
//   stresses.c12 = pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*(1.2E1/1.1E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(1.7E1/2.64E2)-beta*(pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*(4.1E1/3.3E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(1.0/5.28E2)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,3.0)*(4.0/8.1E1)-1.0/pow(c11*c22-c12*c12,2.0)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*pow(c11-c12+c22,4.0)*(1.0/8.1E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(7.0/1.98E2)-c12*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0)*pow(c11-c12+c22,4.0)*(4.0/8.1E1)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(7.0/1.98E2)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*(c11-c12+c22)*(7.0/1.98E2)+c12*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(7.0/1.98E2))+c12*gamma*pow(burgers-c11*c22+c12*c12,3.0)*8.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,3.0)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(1.0/2.7E1)+c12*theta*(burgers-c11*c22+c12*c12)*4.0+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,2.0)*(1.0/9.0)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(8.0/3.3E1)-c12*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12+c22,3.0)*(1.0/9.0)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(8.0/3.3E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*(c11-c12+c22)*(8.0/3.3E1)+c12*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(8.0/3.3E1);


    return stresses;
}

struct cella_stru riduci_matrix_reductions_stress(const struct cella_stru& c, const struct matrix_stru& m){

	struct cella_stru d;

	d.c11 = c.c11*m.m11*m.m11 + c.c12*m.m11*m.m12 + c.c22*m.m12*m.m12;
	d.c22 = c.c11*m.m21*m.m21 + c.c12*m.m21*m.m22 + c.c22*m.m22*m.m22;
	d.c12 = m.m21*(2 * c.c11*m.m11 + c.c12*m.m12) + m.m22*(c.c12*m.m11 + 2 * c.c22*m.m12);

	return d;

}



void write_to_vtk(const struct conf_stru& c, int t){

        int id;
        string filename;

        int nx =  c.nx;
        int ny =  c.ny;
        int n =   c.n;
        int nt=   c.myVector.size();

        filename = get_outputPaht() + "/dir_vtk/vtk_" + IntToStr(t) + ".vtk";
        std::ofstream par;
        par.open(filename.c_str()); // nom du fichier qui contient le maillage
                                    // 
        par <<"# vtk DataFile Version 1.0" << std::endl;
        par <<"2D Unstructured Grid of Linear Triangles" << std::endl;
        par <<"ASCII" << std::endl;
        par <<" " << std::endl;

        par <<"DATASET UNSTRUCTURED_GRID" << std::endl;
        //// representation du solide deformé
        par <<"POINTS " << n<< " float" << std::endl;
        for (int id=0;id<n;id++) par<< c.p[id].x <<" "<< c.p[id].y <<"  "<< 0 <<std::endl;
        par <<"CELLS " <<  nt<< " " << (4)*nt<<std::endl; // 4 car 3 dim + 1
        for (const auto& element : c.myVector) {
                int i = element.first;
                int k = element.second;
                struct punto_stru p = c.p[i];
                par<<"3 "<< i <<" "<< p.nn[k]<<" "<<p.nn[k+1]<<std::endl;
        };
        par <<" " <<std::endl;
        par <<"CELL_TYPES " << nt<<std::endl;
        for (int id=0;id<nt;id++) par<<"5"<<std::endl;
        /////////////////////Visualisation des champs sur la configuration deformée
        par <<" " <<std::endl;
        par <<"POINT_DATA " << n <<std::endl;
        /// energie potentielle
        par << "SCALARS Energy float" << std::endl;
        par << "LOOKUP_TABLE default " << std::endl;
         
        for (int id=0;id<n;id++) par << c.energy_local[id][0] << std::endl;

        par <<" " <<std::endl;
        par <<"CELL_DATA " << nt <<std::endl;
        par <<"TENSORS Stress float " <<std::endl;
//      par << "LOOKUP_TABLE default " << std::endl;
        for (const auto& element : c.myVector) {
                int i = element.first;
                int k = element.second;
         
                par << c.stress[i][k].c11 << " " << c.stress[i][k].c12 << " " << 0 << std::endl;
                par << c.stress[i][k].c12 << " " << c.stress[i][k].c22 << " " << 0 << std::endl;
                par << 0 << " " << 0 << " " << 0 << std::endl;
                par << " " << std::endl;
        };

        //// deviateur
//      par <<"SCALARS Deviator float"<<std::endl;
//      par <<"LOOKUP_TABLE default " <<std::endl;
//      for (int id=0;id<n;id++) {par<<  devVMdes[][id]<<std::endl; }
};

void write_to_a_ovito_file( struct conf_stru& c,const struct boundary_conditions& setnew, int t){


	int nx =  c.nx;
	int ny =  c.ny;
	int n =   c.p.size();

	string filename;




	filename = get_outputPaht() + "/dir_xyz/for_ovito_" + IntToStr(t) + ".xyz";

	ofstream filestr;
	filestr.open(filename.c_str());

	filestr <<c.p.size()<< endl;
	filestr << " " << endl;

	for (int i = 0; i<c.p.size(); ++i){

// 		  int i = element.first;
// 		  int k = element.second;
		int k =0;

		double sqroot_detC = sqrt(c.current_metrics[i][k].det());

		filestr << std::scientific << std::setprecision(16)
			<< c.p[i].x << " "
			<< c.p[i].y << " "
			<< (c.energy_local[i][k] ) << " "
			<< (c.current_metrics[i][k].c12 ) << " "
			<< (c.current_metrics[i][k].c11 ) << " "
			<< (c.current_metrics[i][k].c22 ) << " "
			<< (c.stress[i][k].c12 ) << " "
			<< (c.stress[i][k].c11 ) << " "
			<< (c.stress[i][k].c22 ) << " "
			<< contraction(c.stress[i][k],setnew) << " "

			<< sqroot_detC << " "
			<< c.load
			<< endl;
	}

 	filestr.close();

//****************************************************************************************



}






struct de_stru Piola1(const struct cella_stru& dtemp,const struct base_stru& v1, const struct matrix_stru& m){

	
	//M\SigmaM^T
	
	struct cella_stru temp2 = riduci_matrix_reductions_stress(dtemp, m);

	struct de_stru  d;	
	
		
	//P11 = phi,11
	d.gr[1].x =  (temp2.c11 * 2 * v1.e1[0] + temp2.c12*v1.e2[0]);
	//P21 = phi,21
	d.gr[1].y =  (temp2.c11 * 2 * v1.e1[1] + temp2.c12*v1.e2[1]);
	//P12 = phi,12
	d.gr[2].x =  (temp2.c22 * 2 * v1.e2[0] + temp2.c12*v1.e1[0]);
	//P22 = phi,22
	d.gr[2].y =  (temp2.c22 * 2 * v1.e2[1] + temp2.c12*v1.e1[1]);

	

	//forces  on the node 0     
    d.gr[0].x = -d.gr[1].x - d.gr[2].x;
 	d.gr[0].y = -d.gr[1].y - d.gr[2].y;


	return d;


}


struct de_stru cauchy(const struct cella_stru& dtemp,const struct base_stru& v1, const struct matrix_stru& m){

	struct de_stru  d;	
	
	//node east
	d.gr[1].x =  dtemp.c11;
	d.gr[1].y =  0;
	//node north
	d.gr[2].x =  dtemp.c12;
	d.gr[2].y =  0;


	//forces  on the node 0     
    d.gr[0].x = -d.gr[1].x - d.gr[2].x;
 	d.gr[0].y = -d.gr[1].y - d.gr[2].y;
	


	return d;


}





//(F^TF)

struct cella_stru faicella(const struct base_stru& v1){
	struct cella_stru metrics;

	metrics.c11 = v1.e1[0] * v1.e1[0] + v1.e1[1] * v1.e1[1];
	metrics.c22 = v1.e2[0] * v1.e2[0] + v1.e2[1] * v1.e2[1];
	metrics.c12 = v1.e1[0] * v1.e2[0] + v1.e1[1] * v1.e2[1];
	return metrics;
}


//(F^T+F-2I)/2
struct cella_stru epsilon_lineaire(const struct base_stru& v1){
	struct cella_stru metrics;

	metrics.c11 = v1.e1[0] -1;
	metrics.c22 = v1.e2[1] -1;;
	metrics.c12 = v1.e2[0];
	return metrics;
}

void init(struct conf_stru& c)
{
	int nx = c.nx;
	int ny = c.ny;
	cout<<"nx: "<<nx<<" ny: "<<ny<<" n: "<<c.p.size()<<endl;

	double dnx =  0.;//dimensions::nx;
	double dny =  0.;//c.dny;

	double dc = 1;

	int ix, iy, ntot = c.p.size();
	double h=1;
		int n_temp=0;

	for (int i = 0; i<c.p.size(); ++i){
		int iy = i / nx;
		int ix = i % nx;


		double dix = double(ix - dnx);
		double diy = double(iy - dny);

		c.p[i].x = dc*(dix - dnx)*h;
		c.p[i].y = dc*(diy - dny)*h;

	
		
		c.pfix[i].x = dc*(dix - dnx)*h;
		c.pfix[i].y = dc*(diy - dny)*h;

		c.disp[i].x = 0;
		c.disp[i].y = 0;

		//  NEIGHBORDS OF POINTS TRIANGULATION

		c.p[i].nn[0] = (iy)*nx + ix + 1; //EAST
		c.p[i].nn[1] = (iy + 1)*nx + ix; //NORTH    
		c.p[i].nn[2] = (iy)*nx + ix - 1; //WEST
		c.p[i].nn[3] = (iy - 1)*nx + ix; //SOUTH
		



		//IF PERIODIC BOUNDARY CONDITION WANTED	 OF POINTS

		c.p[i].boundary = 0;
		
		if (ix == nx - 1){
			c.p[i].nn[0] = (iy)*nx; 
			c.p[i].boundary = 1;
		}	         //EAST	    
		if (iy == ny - 1){
			c.p[i].nn[1] = ix; 
			c.p[i].boundary = 1;
		}        //NORTH			
		if (ix == 0){
			c.p[i].nn[2] = (iy)*nx + nx - 1; 
			c.p[i].boundary = 1;
		}	 //WEST	    
		if (iy == 0){
			c.p[i].nn[3] = ((ny - 1))*nx + ix;
			c.p[i].boundary = 1;
		} //SOUTH			

		c.p[i].nn[4] = c.p[i].nn[0];
// 		if (ix == nx - 1)
// 			c.p[i].nn[0] = (iy)*nx; 	         //EAST	    
// 		if (iy == ny - 1)
// 			c.p[i].nn[1] = ix;              //NORTH			
// 		if (ix == 0)
// 			c.p[i].nn[2] = (iy)*nx + nx - 1; 	 //WEST	    
// 		if (iy == 0)
// 			c.p[i].nn[3] = (iy + (ny - 1))*nx + ix; //SOUTH			

		
		
// 		cout<<"ix: "<<ix<<" iy: "<<iy<<"i: "<<i<<endl;
// 		cout<<c.p[i].nn[0]<<" "<<c.p[i].nn[1]<<" "<<c.p[i].nn[2]<<" "<<c.p[i].nn[3]<<" "<<c.p[i].nn[4]<<endl;


		if(c.p[i].boundary != 1){
	    	c.idx_alglib.push_back(std::make_pair(i, n_temp++));
	    }
		if(c.p[i].boundary == 1){
	    	c.idx_boundary.push_back(i);
	    }



	}


	// if fix bouandaries are used some of the triangles do exist; if it is the case c.bc_cond[i][k]=0; else c.bc_cond[i][k]=2
	for(int k =0;k<4;k++)
		for(int i = 0; i<c.p.size(); ++i){
			int iy = i / nx;
			int ix = i % nx;
			c.bc_cond[i][k] = 2;	
			//change the boundary condition on x
			if((ix == nx-1 && k == 0) || (ix == 0 && k == 2)  )
				c.bc_cond[i][k] = 0;
			//change the boundary condition on y, overrides on corners
			if((iy == ny-1 && k == 0) || (iy == 0 && k == 2)  )
				c.bc_cond[i][k] = 0;
			
			
			if(c.bc_cond[i][k] == 2 && (k==0 || k==2)  ){
				c.myVector.push_back(std::make_pair(i, k));
				
			}

		}
	
	cout<<"number_of_triangles: "<<c.myVector.size()<<std::endl;
		
	string neigbor_list  =  "neigbor_list.txt";

	ofstream file_neigbor_list = openFileInOutputDirectory(neigbor_list);
	for(int i=0;i<c.p.size();++i){
		file_neigbor_list << std::scientific << std::setprecision(16)<<i<<" "<<
		c.p[i].nn[0]<<" "<<
		c.p[i].nn[1]<<" "<<
		c.p[i].nn[2]<<" "<<
		c.p[i].nn[3]<<endl;
	
	}

	file_neigbor_list.close();



	string triangles  =  "triangles.txt";

	ofstream file_triangles = openFileInOutputDirectory(triangles);
	for (const auto& element : c.myVector){
		  int i = element.first;
		  int k = element.second;
		  struct punto_stru p = c.p[i];

		  file_triangles << std::scientific << std::setprecision(16)<<i<<" "<<
		  p.nn[k]<<" "<<
		  p.nn[k+1]<<" "<<endl;
	
	}

	file_triangles.close();

	int limit = 2;
// 	c.idx_coarse_grain.resize(0);
	
    for (int s = -limit; s <= limit; s++) {
        for (int l = -limit; l <= limit; l++) {
		
		

			c.idx_coarse_grain.push_back(std::make_pair(s,l));
			
			
		}
		
	}
	
	cout<<"idx_coarse_grain: "<<c.idx_coarse_grain.size()<<std::endl;


}

void apply_BC_CG(struct conf_stru& c){
						
			
	for (const auto& i : c.idx_boundary) {


		 c.p[i].x     = c.full[i].x; 
		 c.gr[i].x = 0;
	 
		 c.p[i].y     = c.full[i].y; 
		 c.gr[i].y = 0;
	}
	 
	
	
}


void calc_energy_forces(struct conf_stru& c){

	int n = c.p.size();
	
	
// 	struct punto_stru zeroStruct{0, 0};
//     std::fill(c.gr.begin(), c.gr.end(), zeroStruct);
	
	for (int i=0; i<n; ++i){
		c.gr[i].x=c.gr[i].y=0;	
		c.grp1[i][0].x = c.grp1[i][0].y=0;
		c.grp1[i][1].x = c.grp1[i][1].y=0.;
		c.grp1[i][2].x = c.grp1[i][2].y=0.;
		c.grp1[i][3].x = c.grp1[i][3].y=0.;
		c.grp1[i][4].x = c.grp1[i][4].y=0.;
		c.grp1[i][5].x = c.grp1[i][5].y=0.;
		c.grp1[i][6].x = c.grp1[i][6].y=0.;
		c.grp1[i][7].x = c.grp1[i][7].y=0.;
		c.grp1[i][8].x = c.grp1[i][8].y=0.;
		c.grp1[i][9].x = c.grp1[i][9].y=0.;
	
	}
	
	//triangle direction
	int a = 1;
	double energy_thread=0;
	static std::array<double, 4> shape;
	shape[0] = 1.;
	shape[1] = -1.;
	shape[2] = -1.;
	shape[3] = -1.;

	#pragma omp parallel
	#pragma omp for reduction(+:energy_thread)

	for (const auto& element : c.myVector) {
		  int i = element.first;
		  int k = element.second;
		  int a;
		  int l=3*k;

		  a = shape[k];

		//loop over nodes
// 		for (int i=0; i<n; i++){

			
            struct punto_stru p = c.p[i];
			struct de_stru d;
// 			double energy_thread;
			//REAL BASIS VECTORS
			struct base_stru v;
			//REAL METRICS
			struct cella_stru metrics,metrics_reduced,metrics_elastic;
			//STRESS
			struct cella_stru cauchy,stress;
			struct matrix_stru m;
			struct energy_stress temp;

			

//---------------------------------CONSTRUCT THE TRIANGLE------------------------------------				
			//positions of node origin
			d.p[0] = p;
			//positions of node east or west if k=2
			d.p[1] = c.p[p.nn[k]];
			//positions of node north  or south if k=2
			d.p[2] = c.p[p.nn[k+1]];
			

//---------------------------------Deformation gradient------------------------------------		
			v.e1[0]= a*(d.p[1].x-d.p[0].x) ; // F11
			v.e1[1]= a*(d.p[1].y-d.p[0].y) ; // F21
			v.e2[0]= a*(d.p[2].x-d.p[0].x) ; // F12
			v.e2[1]= a*(d.p[2].y-d.p[0].y) ; // F22		
				
//---------------------------------METRICS------------------------------------		

			//linear or not
			metrics  = c.fptr3(v);
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


void save_results(struct conf_stru& c,  alglib::real_1d_array& starting_point){
	
		int n = c.p.size();
		

// 	for (int i=0; i<n; ++i)
// 		c.gr[i].x=c.gr[i].y=c.energy_local[i][0]=c.energy_local[i][2]=0;	

	int n_a = starting_point.length()/2;
	for (const auto& element : c.idx_alglib) {
		 int i = element.first;
		 int k = element.second;
		 c.p[i].x = starting_point[k]     + c.pfix[i].x ;
		 c.p[i].y = starting_point[k+n_a] + c.pfix[i].y ;
	} 
	 


	//triangle direction
	int a = 1;
	double energy_thread=0;
	std::array<double, 4> shape;
	shape[0] = 1.;
	shape[1] = -1.;
	shape[2] = -1.;
	shape[3] = -1.;

	#pragma omp parallel
	#pragma omp for 

	for (const auto& element : c.myVector) {
		  int i = element.first;
		  int k = element.second;
		  int a;
		  int l=3*k;

		  a = shape[k];

		//loop over nodes
// 		for (int i=0; i<n; i++){

			
            struct punto_stru p = c.p[i];
			struct de_stru d;
// 			double energy_thread;
			//REAL BASIS VECTORS
			struct base_stru v;
			//REAL METRICS
			struct cella_stru metrics,metrics_reduced,metrics_elastic;
			//STRESS
			struct cella_stru cauchy,stress;
			struct matrix_stru m, def_grad,def_grad_inv;
			struct energy_stress temp;
			double detf;
			double P11 ;//= d.gr[1].x;
			double P12  ;//= d.gr[2].x;
			double P21  ;//= d.gr[1].y;
			double P22  ;//= d.gr[2].y;
	
			 

//---------------------------------CONSTRUCT THE TRIANGLE------------------------------------				
			//positions of node origin
			d.p[0] = p;
			//positions of node east or west if k=2
			d.p[1] = c.p[p.nn[k]];
			//positions of node north  or south if k=2
			d.p[2] = c.p[p.nn[k+1]];
			

//---------------------------------Deformation gradient------------------------------------		
			v.e1[0]= a*(d.p[1].x-d.p[0].x) ; // F11
			v.e1[1]= a*(d.p[1].y-d.p[0].y) ; // F21
			v.e2[0]= a*(d.p[2].x-d.p[0].x) ; // F12
			v.e2[1]= a*(d.p[2].y-d.p[0].y) ; // F22		
			
			
// 			def_grad.m11 = v.e1[0];
// 			def_grad.m22 = v.e2[1];
// 			def_grad.m12 = v.e2[0];
// 			def_grad.m21 = v.e1[1];
			
// 			detF = sqrt(v.e1[0]*v.e2[1] - v.e1[1]*v.e2[0]); 
// 			v.e1[0]/=detF;
// 			v.e1[1]/=detF;
// 			v.e2[0]/=detF;
// 			v.e2[1]/=detF;
				
//---------------------------------METRICS------------------------------------		

			metrics  = c.fptr3(v);
  	 	    metrics_reduced =  lagrange_reduction(metrics);
  	 	    metrics_elastic  =  elastic_reduction(metrics);

  	 	    c.current_metrics[i][k].m      = metrics_reduced.m;
  	 	    c.current_metrics[i][k].m_el   = metrics_elastic.m_el;
  	 	    
  	 	    c.current_metrics[i][k].c11    = metrics.c11;
  	 	    c.current_metrics[i][k].c22    = metrics.c22;
  	 	    c.current_metrics[i][k].c12    = metrics.c12;
  	 	    





 			temp = c.fptr(metrics_reduced);

	 		stress = temp.r_stress;
	 		d =  c.fptr4(stress, v, metrics_reduced.m);

//     		 detf =    (v.e1[0]*v.e2[1] - v.e2[0]*v.e1[1]) ;
			 
			 
			 cauchy.c11 = P11 = d.gr[1].x;
			 cauchy.c12 = P12 = d.gr[2].x;
			 cauchy.c21 = P21 = d.gr[1].y;
			 cauchy.c22 = P22 = d.gr[2].y;
// 	
// 			 cauchy.c12  = (p_f21*v.e1[0]+p_f22*v.e2[0])/detf;
// 			 cauchy.c11  = (p_f21*v.e1[0]+p_f22*v.e2[0])/detf;
// 
// 			 cauchy.c21  =  0;//(p_f11*v.e1[0]+p_f12*v.e2[0])/detf;
// 	         cauchy.c22  =  0;//(p_f21*v.e1[1]+p_f22*v.e2[1])/detf;

// 			 def_grad_inv = def_grad.inverse(def_grad);
// 			 cout<<"itself: "<<endl;
// 			 def_grad.print();
// 			 cout<<"inverse: "<<endl;
// 			 def_grad_inv.print();


//     if (std::isnan(def_grad_inv.m11)) // Checking if the variable is NaN
//     {
//         std::cout << "Error: m11 is NaN. Exiting the program." << std::endl;
//         exit(EXIT_FAILURE); // Exiting the program with an error status
//     }
//     
//         if (std::isnan(def_grad_inv.m22)) // Checking if the variable is NaN
//     {
//         std::cout << "Error: m22 is NaN. Exiting the program." << std::endl;
//         exit(EXIT_FAILURE); // Exiting the program with an error status
//     }
// 
// 
// 
//     if (std::isnan(def_grad_inv.m12)) // Checking if the variable is NaN
//     {
//         std::cout << "Error: m12 is NaN. Exiting the program." << std::endl;
//         exit(EXIT_FAILURE); // Exiting the program with an error status
//     }
// 
// 
//     if (std::isnan(def_grad_inv.m21)) // Checking if the variable is NaN
//     {
//         std::cout << "Error: m21 is NaN. Exiting the program." << std::endl;
//         exit(EXIT_FAILURE); // Exiting the program with an error status
//     }


// 	        cauchy.c11 = def_grad_inv.m21 * P11 + def_grad_inv.m22 * P21;
// 			cauchy.c12 = def_grad_inv.m21 * P12 + def_grad_inv.m22 * P22;
// 			cauchy.c21 = def_grad_inv.m11 * P11 + def_grad_inv.m12 * P21;
//  			cauchy.c22  =def_grad_inv.m11 * P12 + def_grad_inv.m12 * P22;

// 			 cauchy.c22  = (p_f21*v.e1[0]+p_f22*v.e2[0])/detf;

// 			 cauchy.c12  =  0;//(p_f11*v.e1[0]+p_f12*v.e2[0])/detf;
// 	         cauchy.c21  =  0;//(p_f21*v.e1[1]+p_f22*v.e2[1])/detf;


// 
// S[0][0] = invF[1][0] * P[0][0] + invF[1][1] * P[1][0];
// S[0][1] = invF[1][0] * P[0][1] + invF[1][1] * P[1][1];
// S[1][0] = invF[0][0] * P[0][0] + invF[0][1] * P[1][0];
// S[1][1] = invF[0][0] * P[0][1] + invF[0][1] * P[1][1];

// 			C_{12} = F_{11}P_{21} + F_{12}P_{22}	
// 
//  			C_{11} = F_{11}P_{11} + F_{12}P_{12}
// 
// 			C_{22} = F_{21}P_{21} + F_{22}P_{22}

 			c.energy_local[i][k]    = temp.energy;
 			c.stress[i][k] = cauchy;

//---------------------------------FORCES ON  NODES------------------------------------		




	 	}
 
    
	c.energy = sum_2d_vector(c.energy_local);
}

void choose_bc( boundary_conditions& set1,const string& bc){

		if(bc.compare("shear") ==0)
			set1.macro_shear();
// 		if(BLK1::macro_def.compare("shear_comb") ==0)
// 			set1.macro_multipled_loading();
// 		if(BLK1::macro_def.compare("tension") ==0)
// 			set1.macro_compression();
// 		if(BLK1::macro_def.compare("uni_tension") ==0)
// 			set1.uni_macro_compression();
// 		if(BLK1::macro_def.compare("shear_rhombic") ==0)
// 			set1.macro_shear_rhombic();
// 		if(BLK1::macro_def.compare("shear_rect") ==0)
// 			set1.macro_shear_rectangular();
// 		if(BLK1::macro_def.compare("rhombic_tension") ==0)
// 			set1.macro_rhombic_compression();
// 		if(BLK1::bc_mode.compare("frame") ==0)
// 			set1.def_frame(force_energy::calculation::getInstance().c,set1);

		
}


void apply_BC_body(struct conf_stru& c, const struct boundary_conditions& set){

	   for (const auto& i : c.idx_boundary) {
			
		c.full[i].x = c.pfix[i].x * set.f.f11 + c.pfix[i].y*set.f.f12;
		c.full[i].y = c.pfix[i].x * set.f.f21 + c.pfix[i].y*set.f.f22;
	}

}


void apply_BC_body_full(struct conf_stru& c, const struct boundary_conditions& set){

	for(int i =0;i<c.p.size();++i){	
			
		c.full[i].x = c.pfix[i].x * set.f.f11 + c.pfix[i].y*set.f.f12;
		c.full[i].y = c.pfix[i].x * set.f.f21 + c.pfix[i].y*set.f.f22;
	}

}

void deform_current_state(struct conf_stru& c, const struct boundary_conditions& set){

	for(int i =0;i<c.p.size();++i){	
			
		double x = c.p[i].x * set.f.f11 + c.p[i].y*set.f.f12;
		double y = c.p[i].x * set.f.f21 + c.p[i].y*set.f.f22;
		c.p[i].x  = x;
		c.p[i].y  = y; 
	}

}


void initial_guess(const struct conf_stru& c, struct boundary_conditions& set, alglib::real_1d_array& starting_point){

	int mm = starting_point.length()/2;

	for (const auto& element : c.idx_alglib) {
		 int i = element.first;
		 int k = element.second;
				
		starting_point[k]    = -c.pfix[i].x  + (c.p[i].x * set.f.f11 + c.p[i].y*set.f.f12);
		starting_point[k+mm] = -c.pfix[i].y  + (c.p[i].x * set.f.f21 + c.p[i].y*set.f.f22);
	}

}


double  plasticity(std::vector <std::vector<cella_stru> > current_metrics_t,std::vector <std::vector<cella_stru> > current_metrics_t0){

 	int sum = 0;
	fem_2D::calculation& singleton =  fem_2D::calculation::getInstance();

	for (const auto& element : singleton.c.myVector) {
		  int i = element.first;
		  int k = element.second;
   		 sum += 	std::round(current_metrics_t[i][k].c12) - std::round(current_metrics_t0[i][k].c12);
   }
   
   return sum;
}


double   plasticity_gl2z( std::vector <std::vector<cella_stru> >& current_metrics_t,
					      std::vector <std::vector<cella_stru> >& current_metrics_t0){

 	double sum = 0;
 	int    counter=0;
	fem_2D::calculation& singleton =  fem_2D::calculation::getInstance();

	for (const auto& element : singleton.c.myVector) {
		  int i = element.first;
		  int k = element.second;

   		 sum += 	std::abs(std::round(current_metrics_t[i][k].m_el.m11) - std::round(current_metrics_t0[i][k].m_el.m11));
   		 sum += 	std::abs(std::round(current_metrics_t[i][k].m_el.m22) - std::round(current_metrics_t0[i][k].m_el.m22));
   		 sum += 	std::abs(std::round(current_metrics_t[i][k].m_el.m12) - std::round(current_metrics_t0[i][k].m_el.m12));
   		 sum += 	std::abs(std::round(current_metrics_t[i][k].m_el.m21) - std::round(current_metrics_t0[i][k].m_el.m21));
   		 if(sum>0){
//    		 	cout<<"t-----"<<endl;
//    		 	current_metrics_t[i][k].m_el.print();
//    		 	cout<<"t0-----"<<endl;
//    		 	current_metrics_t0[i][k].m_el.print();
			sum = 0;
   		 	counter++;
   		 }

   }
   
   return (double) counter;
   
//    if (sum == 0) return false;
//    else return true;
}


void save_elastic_state_gl2z(const std::vector <std::vector<cella_stru> >& current_metrics,
					               std::vector <std::vector<cella_stru> >& current_metrics_t0){
	
	fem_2D::calculation& singleton =  fem_2D::calculation::getInstance();

	for (const auto& element : singleton.c.myVector) {
			 int i = element.first;
			 int k = element.second;
			 current_metrics_t0[i][k].m_el  = current_metrics[i][k].m_el;
			 current_metrics_t0[i][k].m     = current_metrics[i][k].m;
			 current_metrics_t0[i][k].c11     = current_metrics[i][k].c11;
			 current_metrics_t0[i][k].c22     = current_metrics[i][k].c22;
			 current_metrics_t0[i][k].c12     = current_metrics[i][k].c12;

// 	   	 	    c.current_metrics[i][k].c11    = metrics.c11;
//   	 	    c.current_metrics[i][k].c22    = metrics.c22;
//   	 	    c.current_metrics[i][k].c12    = metrics.c12;

	}
}


void transfer_from_alglib_to_std_v2( struct conf_stru& c,const alglib::real_1d_array& starting_point){

	int mm = starting_point.length()/2;



	for (const auto& element : c.idx_alglib) {
		 int i = element.first;
		 int k = element.second;
		 c.p[i].x = starting_point[k]    + c.pfix[i].x ;
		 c.p[i].y = starting_point[k+mm] + c.pfix[i].y ;

	}

}



void function1_foralglib2D_DOF_version(const alglib::real_1d_array &x, double &func, alglib::real_1d_array &grad, void *ptr) 
{

	fem_2D::calculation& singleton =  fem_2D::calculation::getInstance();
	
	int n = singleton.c.p.size();
	int n_a = x.length()/2;
	
	int ntemp = 0;

	func=0;



	for (const auto& element : singleton.c.idx_alglib) {
		 int i = element.first;
		 int k = element.second;
		 singleton.c.p[i].x = x[k]     + singleton.c.pfix[i].x ;
		 singleton.c.p[i].y = x[k+n_a] + singleton.c.pfix[i].y ;
	}

	
//     apply_BC_CG(singleton.c);
	//calculate forces and energy on each node
	calc_energy_forces(singleton.c);
	
	
	func = 	singleton.c.energy;
	


	for (const auto& element : singleton.c.idx_alglib) {
		 int i = element.first;
		 int k = element.second;
		 grad[k]      =  singleton.c.gr[i].x; ;
		 grad[k+n_a]  =  singleton.c.gr[i].y ;

	}
}

void write_data_to_files(double load, double sum_energy,
struct average_stress sum_stress,
double sum_p_energy,
struct average_stress sum_p_stress,
double plas){

	
	fem_2D::calculation& singleton =  fem_2D::calculation::getInstance();

	fstream filestr;
	string filename;
	filename=get_outputPaht() + "alpha_stress.dat";
	filestr.open (filename, fstream::in | fstream::out | fstream::app);
	filestr << std::scientific << std::setprecision(16)  <<" "<<load << " "<<sum_p_stress.stress_total/singleton.c.myVector.size()<<endl ;

	filestr << std::scientific << std::setprecision(16)  <<" "<<load << " "<<sum_stress.stress_total/singleton.c.myVector.size()<<endl ;
	filestr.close();

	filename=get_outputPaht() + "alpha_energy.dat";
	filestr.open (filename, fstream::in | fstream::out | fstream::app);
	filestr << std::scientific << std::setprecision(16)  <<" "<<load << " "<<sum_p_energy/singleton.c.myVector.size()<<endl ;
	filestr << std::scientific << std::setprecision(16)  <<" "<<load << " "<<sum_energy/singleton.c.myVector.size()<<endl ;
	filestr.close();

	if(std::round(plas) > 0){
		filename=get_outputPaht() + "alpha_stress12_drop.dat";
		filestr.open (filename, fstream::in | fstream::out | fstream::app);
		filestr << std::scientific << std::setprecision(16)  <<" "<<load << " "<<(sum_p_stress.stress12-sum_stress.stress12)<<endl ;
		filestr.close();

		filename=get_outputPaht() + "alpha_stress11_drop.dat";
		filestr.open (filename, fstream::in | fstream::out | fstream::app);
		filestr << std::scientific << std::setprecision(16)  <<" "<<load << " "<<(sum_p_stress.stress11-sum_stress.stress11)<<endl ;
		filestr.close();

		filename=get_outputPaht() + "alpha_stress22_drop.dat";
		filestr.open (filename, fstream::in | fstream::out | fstream::app);
		filestr << std::scientific << std::setprecision(16)  <<" "<<load << " "<<(sum_p_stress.stress22-sum_stress.stress22)<<endl ;
		filestr.close();

		filename=get_outputPaht() + "alpha_stresstotal_drop.dat";
		filestr.open (filename, fstream::in | fstream::out | fstream::app);
		filestr << std::scientific << std::setprecision(16)  <<" "<<load << " "<<(sum_p_stress.stress_total-sum_stress.stress_total)
		<<" "<<sum_p_stress.stress_total/singleton.c.myVector.size()<<endl ;
		filestr.close();



		filename=get_outputPaht() + "alpha_energy_drop.dat";
		filestr.open (filename, fstream::in | fstream::out | fstream::app);
		filestr << std::scientific << std::setprecision(16)  <<" "<<load << " "<<(sum_p_energy-sum_energy)
		<<" "<<sum_p_stress.stress_total/singleton.c.myVector.size()<<endl ;
		filestr.close();

		filename=get_outputPaht() + "alpha_plastic_drop.dat";
		filestr.open (filename, fstream::in | fstream::out | fstream::app);
		filestr << std::scientific << std::setprecision(16)  <<" "<<load << " "<<plas
		<<" "<<sum_p_stress.stress_total/singleton.c.myVector.size()<<endl ;
		filestr.close();

	}

}

void memory(alglib::real_1d_array &x,struct conf_stru& c, double& load_value){



		ifstream disp;
		
		string filename = "/Users/usalman/gl_avalanches/conti_potential/n100/dir1683645378/dir_xyz/for_ovito_10028.xyz";

		disp.open(filename.c_str());
		double temp;
		int nol;
		disp >> nol;
		cout<<"number of lines in the file: "<<" "<<nol<<endl;
		
		for (int i = 0; i<c.p.size(); ++i){
			disp >> std::scientific >> std::setprecision(16)>> c.p[i].x >> c.p[i].y
				 >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> load_value ;
		
		}

		cout<<"load_value: "<<" "<<load_value<<endl;

		int mm = x.length()/2;

		for (const auto& element : c.idx_alglib) {
			 int i = element.first;
			 int k = element.second;
		
			x[k]    = -c.pfix[i].x  + c.p[i].x ;
			x[k+mm] = -c.pfix[i].y  + c.p[i].y ;
		}



}



int main(int argc, char** argv) {

	int nx;
	int ny;
	double load;
	double load_increment;
   
    if (argc > 4) {
        // read nx and ny from command-line arguments
        nx =             std::atoi(argv[1]);
        ny              = std::atoi(argv[2]);
        load           = std::atof(argv[3]);
        load_increment = std::atof(argv[4]);

        std::cout << "the value of nx: "<<nx<<endl;
        std::cout << "the value of ny: "<<ny<<endl;
        std::cout << "the value of load: "<<load<<endl;
        std::cout << "the value of load_increment: "<<load_increment<<endl;


    } else {
        // prompt user for nx and ny
        //std::cout << "Enter the value of nx: ";
        //std::cin >> nx;
        //std::cout << "Enter the value of ny: ";
        //std::cin >> ny;
		nx = 5;
		ny = 5;
		load=0.14;
		load_increment=0.001;
    }
	int n=nx*ny;
	 
	int method =2; //  1 is 2D with forces set to zero at the boundaries
				   //  method 2 will be coded later	

	//find the time from bing bang
	auto now = std::chrono::system_clock::now().time_since_epoch();
	// Convert the duration to an integer seed
	int seed = static_cast<unsigned int>(now.count());//         std::cin >> seed;


	ofstream sizes = openFileInOutputDirectory("sizes.dat");
	sizes << nx << endl;
	sizes << ny << endl;
	sizes.close();
	ofstream save_seed = openFileInOutputDirectory("seed.dat");
	save_seed << seed << endl;
	save_seed.close();


	
	std::string path = get_outputPaht() + "/dir_xyz";
    createFolder(path);
	path = get_outputPaht() + "/dir_vtk";
    createFolder(path);
	
	
	fem_2D::calculation& singleton =  fem_2D::calculation::getInstance();

	singleton.c.nx=nx;
	singleton.c.ny=ny;
	singleton.c.n=n;
	singleton.setSize(n,nx,ny);
	
	//init FEM grid
	init(singleton.c);


	//for atomistic change to &coarse_grain_fast or boyer_energy linear_energy
	singleton.c.linearity = false;
	singleton.c.fptr  = &conti_energy;
	//not  used in the current code
	singleton.c.fptr2 = &boyer_energy;

// 	singleton.c.fptr  = &conti_energy;
// 	singleton.c.fptr2 = &conti_stress;

	//epsilon_lineaire, faicella
	singleton.c.fptr3 = &faicella;
	//Piola1, cauchy
	singleton.c.fptr4 = &Piola1;
	//atom_energy
	singleton.c.atom_energy = &square_energy;
	//atom_stress
	singleton.c.atom_stress = &square_energy_der;
	
	atomistic_grid_square(singleton.c);
	cella_stru metric_local;
	metric_local.c11=1;metric_local.c22=1;metric_local.c12=0;
	singleton.c.zeroing = 0;

	energy_stress temp = singleton.c.fptr(metric_local);
	singleton.c.zeroing = temp.energy;
	cout<<"the minimum of the energy: "<< singleton.c.zeroing<<endl;

	double wall0_GLOBAL = get_wall_time();
	double cpu0_GLOBAL  = get_cpu_time();

//////2D FINITE ELEMENTS	

	if(method ==2){	
	
		// singleton().c can be used anywhere with include common_2D.h
		// singleton().c  is a structure of type conf_stru 
		
		//apply Boundary conditions
		double theta = 0. ;
		struct boundary_conditions setnew;
		setnew.setParameters(theta,load);
		string bc="shear";
		choose_bc(setnew,bc);
		//it changes c.full
		apply_BC_body(singleton.c,setnew);
		//it changes c.p
		deform_current_state(singleton.c,setnew);

		//this is created to find initial displacement field
		struct boundary_conditions set_init;
		set_init.setParameters(theta,0);
		choose_bc(set_init,bc);


		int t=0;
		write_to_vtk(singleton.c,t);
		write_to_a_ovito_file(singleton.c,setnew,t++);
		
		//alglib solvers
		//starting_point is the initial guess for the displacement field
		//function1_foralglib2D is the function that calculates the energy and its derivative
		
		alglib::real_1d_array starting_point;
		starting_point.setlength(2*(nx-2)*(ny-2));
		
		std::vector <std::vector<cella_stru> > current_metrics_t0;
		current_metrics_t0.resize(n);
		for (int i = 0; i < n; ++i)
                current_metrics_t0[i].resize(4);

		for (int i=0;i<starting_point.length();++i)
			starting_point(i)=0;


		// std::vector<double> disorder  = brownian(0,0.0001,starting_point.length(),seed);
		std::vector<double> disorder  = wiener(0,0.1,starting_point.length(),seed);
// 		singleton.c.disorder_x   = wiener(0,0.,singleton.c.disorder_x.size(),seed);

// 

		initial_guess( singleton.c, set_init,  starting_point);

		alglib::real_1d_array starting_point_keep;
		starting_point_keep.setlength(2*(nx-2)*(ny-2)  );
		
		for (int i=0;i<starting_point_keep.length();++i)
			starting_point_keep[i]=starting_point[i];
		
		int mm;
		if(singleton.c.linearity == true)  mm = starting_point_keep.length()/2;
		if(singleton.c.linearity == false) mm = starting_point_keep.length();

		for (int i=0;i<mm;++i){
			starting_point_keep[i]+=disorder[i]; 
			starting_point[i]+=disorder[i];
		}

// 		// initial_guess for the displacement
// 		struct boundary_conditions set_incr;
// 		set_incr.setParameters(theta,load_increment);
// 		choose_bc(set_incr,bc);
// 		initial_guess( singleton.c, set_incr,  starting_point);

		// Copy the ALGLIB array to a std::vector
		std::vector<double> myVector2(starting_point.getcontent(), starting_point.getcontent() + starting_point.length());
		std::cout<<"maxval init "<<maxval(myVector2)<<std::endl;
		std::cout<<"minval init "<<minval(myVector2)<<std::endl;

		//adaptive load preparation
		std::vector<double> adaptive(0); // adaptive loading values

		// adaptive.push_back(1000*load_increment);
		// adaptive.push_back(100*load_increment);
		// adaptive.push_back(25*load_increment);
		// adaptive.push_back(10*load_increment);
		// adaptive.push_back(5*load_increment);
		adaptive.push_back(1*load_increment);

			
		//loading loop starting from the zero load
// 		for(int n_load = 0; n_load < 10000; n_load++){
		int n_load = 0;
		int prediction_failure= 0;
		int prediction_failure_special=0;
	

		setnew.setParameters(0,load);
		singleton.c.load=load;

		bc="shear";
		choose_bc(setnew,bc);
		apply_BC_body(singleton.c,setnew);
		apply_BC_CG(singleton.c);

		save_results(singleton.c,starting_point_keep);
		save_elastic_state_gl2z(singleton.c.current_metrics,current_metrics_t0);
		int alglib_counter=0;
	
		int memory_function=0;
		
		if(memory_function==1){
			memory(starting_point,singleton.c, load);

			std::copy(starting_point.getcontent(),
			starting_point.getcontent()+starting_point.length(), 
			starting_point_keep.getcontent());

			setnew.setParameters(0,load);
			singleton.c.load=load;

			bc="shear";
			choose_bc(setnew,bc);
			apply_BC_body(singleton.c,setnew);
			apply_BC_CG(singleton.c);

			save_results(singleton.c,starting_point_keep);
			save_elastic_state_gl2z(singleton.c.current_metrics,current_metrics_t0);
			t=28;
			write_to_vtk(singleton.c,t);
			write_to_a_ovito_file(singleton.c,setnew,t++);

		}
	
		while(load <1.){			
			
		
			double load_current = load;
			struct average_stress sum_p_stress ;
			double sum_p_energy ;
			double plas;

			for(int inner=0;inner<adaptive.size();++inner){
			
				
				load = load_current+adaptive[inner];
				setnew.setParameters(0,load);
				singleton.c.load=load;

				bc="shear";
				choose_bc(setnew,bc);
				apply_BC_body(singleton.c,setnew);
				apply_BC_CG(singleton.c);

				// initial_guess for the displacement


				struct boundary_conditions set_incr;
				double load_increment_temp = adaptive[inner];
				set_incr.setParameters(0,load_increment_temp);
				choose_bc(set_incr,bc);
				
				//find the initial positions to be deformed
				save_results(singleton.c,starting_point_keep);		
				
				
				//deform the current conf and find displacement
				initial_guess( singleton.c, set_incr,  starting_point);
 				//save  the deformed  current conf and find displacement
 				save_results(singleton.c,starting_point);
 				
 				
 				
				sum_p_stress   = sum_2d_cella(singleton.c.stress,setnew);
				sum_p_energy   = sum_2d_vector(singleton.c.energy_local);
				

				double epsg = 0;
				double epsf = 0;
				double epsx = 0;
				alglib::ae_int_t maxits = 0;
				alglib::minlbfgsstate state;
				alglib::minlbfgsreport rep;



				alglib::minlbfgscreate(10, starting_point, state);
				alglib::minlbfgssetcond(state, epsg, epsf, epsx, maxits);

// 				std::cout<<" ALGLIB OPTIMISATION STARTS "<<std::endl;

				double wall0 = get_wall_time();
				double cpu0  = get_cpu_time();
				
  					  
  								
				alglib::minlbfgsoptimize(state, function1_foralglib2D_DOF_version);

				double wall1 = get_wall_time();
				double cpu1  = get_cpu_time();

// 				cout << "Wall Time = " << wall1 - wall0 << endl;
// 				cout << "CPU Time  = " << cpu1  - cpu0  << endl;
// 				cout << "parallelisation ratio   = " << (cpu1  - cpu0)/(wall1 - wall0)  << endl;

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

	double wall1_GLOBAL = get_wall_time();
	double cpu1_GLOBAL  = get_cpu_time();


	cout << "Wall Time = " << wall1_GLOBAL  - wall0_GLOBAL  << endl;
	cout << "CPU Time  = " << cpu1_GLOBAL   - cpu0_GLOBAL   << endl;

	}

		
	


}