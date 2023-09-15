#pragma once
#include "structures.h"
#include "optimization.h"


class boundary_conditions {
public:
	double theta_d, load_in;
	struct  def_grad f;      // Member functions declaration

	struct def_grad macro_shear();
	struct def_grad macro_shear_rhombic();
	struct def_grad macro_shear_rectangular();
	struct def_grad macro_compression();

	void setParameters(double, double);

};