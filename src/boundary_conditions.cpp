
#include "boundary_conditions.h"

struct def_grad boundary_conditions::macro_shear()
{

	double perturb = 0 ;


	f.f11 = (1. - load_in*cos(theta_d + perturb)*sin(theta_d + perturb));
	f.f22 = (1. + load_in*cos(theta_d - perturb)*sin(theta_d - perturb));
	f.f12 = load_in* pow(cos(theta_d), 2.);
	f.f21 = -load_in* pow(sin(theta_d - perturb), 2.);

	return f;


}


struct def_grad boundary_conditions::macro_shear_rhombic()
{

	f.f11 = cosh(load_in);
	f.f22 = cosh(load_in);
	f.f12 = sinh(load_in);
	f.f21 = sinh(load_in);

	return f;


}

struct def_grad boundary_conditions::macro_shear_rectangular()
{


	f.f11 = cosh(load_in/2.) - sinh(load_in/2.);
	f.f22 = cosh(load_in/2.) + sinh(load_in/2.);
	f.f12 = 0.;
	f.f21 = 0.;


	return f;


}





struct def_grad boundary_conditions::macro_compression()
{

	
	f.f11 = sqrt(load_in);	
	f.f22 = sqrt(load_in);
	f.f12 = 0.;
	f.f21 = 0.;




	return f;

}






void  boundary_conditions::setParameters(double value1, double value2)
{
	theta_d = value1* alglib :: pi () / 180.0;
	load_in = value2;

	// 	return f;

}




//it deforms c.pfix


	

