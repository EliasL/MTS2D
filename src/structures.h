#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <array>


using std::string;


struct punto_stru
{
    double x;
    double y;
    int  boundary;
    int  bou_idx;

    std::array<int, 7> nn; // Indices of neighbors  with initial capacity of 7



//     punto_stru() : x(0), y(0) {} // Initialize 
};


struct de_stru
{
	struct punto_stru p[3];
	struct punto_stru gr[3];
};


struct base_stru
{
	double e1[2];
	double e2[2];
};

struct matrix_stru
{
	double m11;
	double m12;
	double m22;
	double m21;

	matrix_stru& operator+(const matrix_stru& rhs)
	{
		m11 += rhs.m11;
		m22 += rhs.m22;
		m21 += rhs.m21;
		m12 += rhs.m12;

		return *this;
	}

	matrix_stru& operator-(const matrix_stru& rhs)
	{
		m11 -= rhs.m11;
		m22 -= rhs.m22;
		m21 -= rhs.m21;
		m12 -= rhs.m12;

		return *this;
	}

	matrix_stru& operator/(const matrix_stru& rhs)
	{
		m11 /= rhs.m11;
		m22 /= rhs.m22;
		m21 /= rhs.m21;
		m12 /= rhs.m12;

		return *this;
	}



	void print(void);

	matrix_stru transpose(void);
	
	
	matrix_stru inverse(const matrix_stru&);
	matrix_stru multiply(const matrix_stru&);
	matrix_stru zero(const double&);
	matrix_stru identity();
	int check(void);
	double det(void);
// 	matrix_stru() : m11(0), m22(0),  m12(0),  m21(0) {} // Initialize 


// 	matrix_stru curl1(const  de_stru& , const  de_stru& ,int );

};

void inline matrix_stru::print()
{
	using std::cout;
	using std::endl;
	cout << " m11= " << m11 << "\n" <<
		" m22= " << m22 << "\n" <<
		" m12= " << m12 << "\n" <<
		" m21= " << m21 << endl;
		
	cout <<"------------"<<endl;

}

int inline matrix_stru::check()
{
	int i = 0;
	if (m11 == 1 && m22 == 1 && m12 == 0 && m21 == 0)
		i = 1;

	return i;

}

matrix_stru inline matrix_stru::transpose()
{
	matrix_stru mt;
	
	mt.m11=m11;
	mt.m22=m22;
	mt.m12=m21;
	mt.m21=m12;

	return mt;

}

double inline matrix_stru::det()
{
	matrix_stru mt;
	mt.m11=m11;
	mt.m22=m22;
	mt.m12=m21;
	mt.m21=m12;
	
	double det =mt.m11*mt.m22-mt.m12*mt.m21;

	return det;

}
matrix_stru inline matrix_stru::inverse(const  matrix_stru& x)
{
	matrix_stru mt;
	double det =x.m11*x.m22-x.m12*x.m21;
	
	mt.m11=m22/det;
	mt.m22=m11/det;
	mt.m12=-m12/det;
	mt.m21=-m21/det;

	return mt;

}

matrix_stru inline matrix_stru::multiply(const matrix_stru& rhs)
{
	matrix_stru k;
	k.m11 = m11*rhs.m11+ m12*rhs.m21;
	k.m12 = m11*rhs.m12+ m12*rhs.m22;
	k.m21 = m21*rhs.m11+ m22*rhs.m21;
	k.m22 = m21*rhs.m12+ m22*rhs.m22;

	return k;
	
}



matrix_stru inline matrix_stru::zero(const double& rhs)
{
	matrix_stru k;
	k.m11 = rhs;
	k.m12 = rhs;
	k.m21 = rhs;
	k.m22 = rhs;

	return k;
	
}

matrix_stru inline matrix_stru::identity()
{
	matrix_stru k;
	k.m11 = 1.;
	k.m22 = 1.;
	k.m21 = 0.;
	k.m12 = 0.;

	return k;
	
}



struct cella_stru
{
	double c11;
	double c22;
	double c12;
	double c21;

	bool plasticity;

	struct matrix_stru m;
	struct matrix_stru m_el;
	
	cella_stru() : c11(0), c22(0),  c12(0),  c21(0) {} // Initialize 

	double det(void);

};


double inline cella_stru::det()
{
	
	
	return c11*c22-c12*c12;


}


struct def_grad
{
	double f11;
	double f22;
	double f12;
	double f21;
	def_grad inverse(const def_grad&);
	def_grad multiply(const def_grad&);
};



struct energy_stress
{
	double energy;
	struct cella_stru r_stress;

};


struct average_stress
{
	double stress11;
	double stress22;
	double stress12;
	double stress21;
	double stress_total;
    average_stress() : stress11(0), stress22(0),stress21(0), stress12(0),stress_total(0) {} // Initialize 


};


struct conf_stru
{
	int n,nx,ny;
	int inc;
	double energy,load;
	double dnx, dny;
	double total_triangular_number;
	double zeroing;

	bool plasticity;
	
	bool linearity;

	std::vector<punto_stru> p, pfix,cpold,disp,gr,full;
	std::vector<double> disorder_x;
	std::vector<double> atom_x_co,atom_y_co;

	std::vector <std::vector<double> > energy_local;
	std::vector <std::vector<cella_stru> > current_metrics;
	std::vector <std::vector<cella_stru> > stress;


	std::vector <std::vector<punto_stru> > grp1;
	std::vector <std::vector<int> > bc_cond;
	std::vector<std::pair<int, int>> myVector;

    std::vector<std::pair<int, int>> idx_alglib;

    std::vector<std::pair<int, int>> idx_coarse_grain;

	std::vector<double> idx_boundary;
	double length0_x,length0_y;


    
    struct energy_stress   (*fptr)(const struct cella_stru& c); //Declare a function pointer to double with  params
    struct energy_stress   (*fptr2)(const struct cella_stru& c); //Declare a function pointer to double with  params
    struct cella_stru   (*fptr3)(const struct base_stru& v); //Declare a function pointer to double with  params
    struct de_stru      (*fptr4)(const struct cella_stru& dtemp,const struct base_stru& v1, const struct matrix_stru& m); //Declare a function pointer to double with  params
    double      (*atom_energy)(const double r); //Declare a function pointer to double with  params
    double      (*atom_stress)(const double r); //Declare a function pointer to double with  params



};