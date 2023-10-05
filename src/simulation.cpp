
#include <vector>
#include "matrix2x2.h"
#include "singelton.h"
#include "grid2D.h"
#include "settings.h"

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

class Cell{
public:
    Matrix2x2<double> D; // deformation gradiant / basis vectors (e1, e2)
    Matrix2x2<double> C; // real_metrics
    Matrix2x2<double> C_; // reduced_metrics
    Matrix2x2<double> m; // reduction transformation ( T(m)Cm = C_ )
    Matrix2x2<double> r_s; // reduces stress
    // Piola-Kirchhoff stress: https://en.wikipedia.org/wiki/Piola%E2%80%93Kirchhoff_stress_tensors
    Matrix2x2<double> P; // First Piola-Kirchhoff stress tensor 
    Matrix2x2<double> S; // Second Piola-Kirchhoff stress tensor
    double energy;
    bool plasticity; // not used yet

    Cell(triangle t, double a) {
        // Calculates D
        get_deformation_gradiant(t, a);

        // Calculates C
        C  = t.metric(a, MetricFunction::faicella);

        // Calculate C_ and m
        lagrange_reduction();

        // Updates reduced stress r_s and energy
        calculate_energy_and_reduced_stress();
    };

    // Avoid using this one
    Cell(){};

    // A basis vector for the cell
    double Cell::e1(int index){
        return D[0][index];
    }

    // A basis vector for the cell
    double Cell::e2(int index){
        return D[1][index];
    }

private:
    void Cell::get_deformation_gradiant(triangle t, double a){
        auto e1_ = t.e1(a);
        auto e2_ = t.e1(a);

        D[0][0] = e1_[0];
        D[0][1] = e1_[1];
        D[1][0] = e2_[0];
        D[1][1] = e2_[1];
    }

   void Cell::lagrange_reduction(){
        // We start by copying the values from C to the reduced matrix
        C_ = C;

        if(LINEARITY){
            // If we assume linearity, we are done. m is already identity.
            return;
        }
        // And then we follow an algorithm generate both m and C_
        while(C_[0][1]<0 || C_[1][1]<C_[0][0] || 2 * C_[0][1]>C_[0][0] ){
            
            if (C_[0][1]<0){
                C_.flip(0,1);
                m.lag_m1();
            }
            
            if (C_[1][1]<C_[0][0]){
                C_.swap(0,0,1,1);
                m.lag_m2();
            }

            if (2 * C_[0][1]>C_[0][0]){
                // The order here matters, don't modify C_[0][1] before using it
                // to calculate C_[1][1].
                C_[1][1] += C_[0][0] - 2 * C_[0][1];
                C_[0][1] -= C_[0][0];
                m.lag_m3();
            }
        } 
    }


    void Cell::calculate_energy_large(){
        double burgers   = 1.;
        double beta = -0.25;
        double K=4.;

        // Uses the reduced metrics
        double c11 = C_[0][0];
        double c22 = C_[1][1];
        double c12 = C_[0][1];

        Singelton& s = Singelton::getInstance();
        energy  = -K*(log((c11*c22-c12*c12)/burgers)-(c11*c22-c12*c12)/burgers)+beta*(pow(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22),2.0)*9.46969696969697E-4-pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),3.0)*(4.1E1/9.9E1)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,4.0)*(1.0/8.1E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(7.0/1.98E2))+pow(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22),2.0)*(1.7E1/5.28E2)+pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),3.0)*(4.0/1.1E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,3.0)*(1.0/2.7E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(8.0/3.3E1);

        energy -= s.zeroing_energy;
        r_s[0][0] = -beta*(pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(4.1E1/3.3E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.0/5.28E2)-((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,3.0)*(4.0/8.1E1)-1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,4.0)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(1.0/8.1E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(7.0/1.98E2)+c22*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0)*pow(c11-c12+c22,4.0)*(2.0/8.1E1)-((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(7.0/1.98E2)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(7.0/1.98E2)-c22*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(7.0/3.96E2))+K*(c22/burgers-c22/(c11*c22-c12*c12))+pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(1.2E1/1.1E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.7E1/2.64E2)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,3.0)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.0/2.7E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,2.0)*(1.0/9.0)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(8.0/3.3E1)+c22*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12+c22,3.0)*(1.0/1.8E1)-((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)+c22*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(8.0/3.3E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*((c11*(1.0/6.0)-c12*(2.0/3.0)+c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)-c22*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)-c22*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(8.0/3.3E1)-c22*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(4.0/3.3E1);

        r_s[1][1] = beta*(pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(4.1E1/3.3E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.0/5.28E2)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,3.0)*(4.0/8.1E1)-1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,4.0)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(1.0/8.1E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(7.0/1.98E2)-c11*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0)*pow(c11-c12+c22,4.0)*(2.0/8.1E1)-((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(7.0/1.98E2)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(7.0/1.98E2)+c11*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(7.0/3.96E2))+K*(c11/burgers-c11/(c11*c22-c12*c12))-pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(1.2E1/1.1E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.7E1/2.64E2)-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,3.0)*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(1.0/2.7E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,2.0)*(1.0/9.0)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(8.0/3.3E1)+c11*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12+c22,3.0)*(1.0/1.8E1)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/3.0)-c11*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/6.0)+1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11*2.0-c22*2.0)*(c11-c12*4.0+c22)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*(3.0/2.0))*(8.0/3.3E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*((c11*(-1.0/6.0)+c12*(2.0/3.0)-c22*(1.0/6.0))/(c11*c22-c12*c12)+(c11*(1.0/2.0)-c22*(1.0/2.0))/(c11*c22-c12*c12)+c11*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1)+c11*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/4.0))*(8.0/3.3E1)-c11*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(4.0/3.3E1);

        r_s[0][1] = pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*(1.2E1/1.1E1)-K*((c12*2.0)/burgers-(c12*2.0)/(c11*c22-c12*c12))+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(1.7E1/2.64E2)-beta*(pow((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12),2.0)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*(4.1E1/3.3E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(1.0/5.28E2)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12+c22,3.0)*(4.0/8.1E1)-1.0/pow(c11*c22-c12*c12,2.0)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*pow(c11-c12+c22,4.0)*(1.0/8.1E1)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(7.0/1.98E2)-c12*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0)*pow(c11-c12+c22,4.0)*(4.0/8.1E1)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(7.0/1.98E2)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*(c11-c12+c22)*(7.0/1.98E2)+c12*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(7.0/1.98E2))-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,3.0)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(1.0/2.7E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12+c22,2.0)*(1.0/9.0)-(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(8.0/3.3E1)-c12*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12+c22,3.0)*(1.0/9.0)+((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/sqrt(c11*c22-c12*c12)*(c11-c12+c22)*(pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*4.0-1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,2.0)*(4.0/3.0)+c12*1.0/pow(c11*c22-c12*c12,5.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/3.0)-c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,5.0/2.0)*(c11-c12*4.0+c22)*3.0)*(8.0/3.3E1)+(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*1.0/sqrt(c11*c22-c12*c12)*((c11*(-2.0/3.0)+c12*(8.0/3.0)-c22*(2.0/3.0))/(c11*c22-c12*c12)+c12*1.0/pow(c11*c22-c12*c12,2.0)*pow(c11-c12*4.0+c22,2.0)*(1.0/6.0)+c12*pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,2.0)*(1.0/2.0))*(c11-c12+c22)*(8.0/3.3E1)+c12*(1.0/pow(c11*c22-c12*c12,3.0/2.0)*pow(c11-c12*4.0+c22,3.0)*(1.0/9.0)-pow(c11-c22,2.0)*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12*4.0+c22))*((pow(c11-c22,2.0)*(1.0/4.0))/(c11*c22-c12*c12)+(pow(c11-c12*4.0+c22,2.0)*(1.0/1.2E1))/(c11*c22-c12*c12))*1.0/pow(c11*c22-c12*c12,3.0/2.0)*(c11-c12+c22)*(8.0/3.3E1);
    }

    void Cell::calculate_energy_and_reduced_stress(){
        double burgers   = 1.;
        double beta = -0.25;
        double K=4.;

        // Uses the reduced metrics
        double c11 = C_[0][0];
        double c22 = C_[1][1];
        double c12 = C_[0][1];

        /* #region Variables */

        double v1_ = c11*c22; // 316 occurances
        double v2_ = c12*c12; // 316 occurances
        double v3_ = c12+c22; // 40 occurances
        double v4_ = c11-c22; // 115 occurances
        double v5_ = c11-v3_; // 40 occurances
        double v6_ = v1_-v2_; // 316 occurances
        double v7_ = pow(v6_,3.0/2.0); // 123 occurances
        double v8_ = pow(c11-(c12*4.0)+c22,3.0); // 50 occurances
        double v9_ = pow(v4_,2.0); // 115 occurances
        double v10_ = pow(c11-(c12*4.0)+c22,2.0); // 65 occurances
        double v11_ = pow(v6_,2.0); // 37 occurances
        double v12_ = pow(v5_,4.0); // 7 occurances
        double v13_ = pow(v5_,3.0); // 10 occurances
        double v14_ = pow(v6_,5.0/2.0); // 33 occurances
        double v15_ = pow(v6_,3.0); // 3 occurances
        double v16_ = pow(v5_,2.0); // 3 occurances
        double v17_ = pow(1.0/v7_*v8_*(1.0/9.0)-(v9_*1.0/v7_)*(c11-(c12*4.0)+c22),2.0); // 2 occurances
        double v18_ = pow((v9_*(1.0/4.0))/v6_+(v10_*(1.0/1.2E1))/v6_,3.0); // 2 occurances
        double v19_ = pow((v9_*(1.0/4.0))/v6_+(v10_*(1.0/1.2E1))/v6_,2.0); // 6 occurances
        double v20_ = sqrt(v6_); // 20 occurances
        double v21_ = (1.0/4.0); // 37 occurances
        double v22_ = (1.0/1.2E1); // 37 occurances
        double v23_ = (1.0/8.1E1); // 4 occurances
        double v24_ = (1.0/9.0); // 37 occurances
        double v25_ = (v9_*1.0/v7_); // 48 occurances
        double v26_ = (c12*4.0); // 58 occurances
        double v27_ = (7.0/1.98E2); // 11 occurances
        double v28_ = (1.0/2.7E1); // 4 occurances
        double v29_ = (8.0/3.3E1); // 11 occurances
        double v30_ = (1.0/6.0); // 30 occurances
        double v31_ = (2.0/3.0); // 15 occurances
        double v32_ = (1.0/2.0); // 25 occurances
        double v33_ = (4.1E1/3.3E1); // 3 occurances
        double v34_ = (1.0/3.0); // 15 occurances
        double v35_ = (c11*2.0-c22*2.0); // 10 occurances
        double v36_ = (3.0/2.0); // 10 occurances
        double v37_ = (1.0/5.28E2); // 3 occurances
        double v38_ = (4.0/8.1E1); // 4 occurances
        double v39_ = (2.0/8.1E1); // 2 occurances
        double v40_ = (7.0/3.96E2); // 2 occurances
        double v41_ = (1.2E1/1.1E1); // 3 occurances
        double v42_ = (1.7E1/2.64E2); // 3 occurances
        double v43_ = (1.0/1.8E1); // 2 occurances
        double v44_ = (4.0/3.3E1); // 2 occurances
        double v45_ = (-1.0/6.0); // 5 occurances
        double v46_ = (-2.0/3.0); // 5 occurances
        double v47_ = (8.0/3.0); // 5 occurances
        double v48_ = (c12*2.0); // 2 occurances
        double v49_ = (4.0/3.0); // 5 occurances
        double v50_ = (v9_*v21_); // 27 occurances
        double v51_ = (v10_*v22_); // 27 occurances
        double v52_ = (c11-v26_+c22); // 58 occurances
        double v53_ = (c11*v30_-c12*v31_+c22*v30_); // 5 occurances
        double v54_ = (c11*v32_-c22*v32_); // 10 occurances
        double v55_ = (c11*v45_+c12*v31_-c22*v30_); // 5 occurances
        double v56_ = (c11*v46_+c12*v47_-c22*v31_); // 5 occurances
        double v57_ = (v50_/v6_+v51_/v6_); // 27 occurances
        double v58_ = (1.0/v7_*v8_*v24_-v25_*v52_); // 33 occurances
        double v59_ = (v53_/v6_+v54_/v6_-c22*1.0/v11_*v10_*v22_-c22*v9_*1.0/v11_*v21_); // 5 occurances
        double v60_ = (v25_-1.0/v7_*v10_*v34_+c22*1.0/v14_*v8_*v30_+1.0/v7_*v35_*v52_-c22*v9_*1.0/v14_*v52_*v36_); // 5 occurances
        double v61_ = (v55_/v6_+v54_/v6_+c11*1.0/v11_*v10_*v22_+c11*v9_*1.0/v11_*v21_); // 5 occurances
        double v62_ = (-v25_+1.0/v7_*v10_*v34_-c11*1.0/v14_*v8_*v30_+1.0/v7_*v35_*v52_+c11*v9_*1.0/v14_*v52_*v36_); // 5 occurances
        double v63_ = (v56_/v6_+c12*1.0/v11_*v10_*v30_+c12*v9_*1.0/v11_*v32_); // 5 occurances
        double v64_ = (v25_*4.0-1.0/v7_*v10_*v49_+c12*1.0/v14_*v8_*v34_-c12*v9_*1.0/v14_*v52_*3.0); // 5 occurances
        /* #endregion */  

        energy  =-K*(log(v6_/burgers)-v6_/burgers)+beta*(v17_*9.46969696969697E-4-v18_*(4.1E1/9.9E1)+v57_*1.0/v11_*v12_*v23_-v58_*v57_*1.0/v20_*v5_*v27_)+v17_*(1.7E1/5.28E2)+v18_*(4.0/1.1E1)-v58_*1.0/v7_*v13_*v28_+v58_*v57_*1.0/v20_*v5_*v29_;

        Singelton& s = Singelton::getInstance();
        energy -= s.zeroing_energy;
        r_s[0][0]  = -beta*(v19_*v59_*v33_+v58_*v60_*v37_-v57_*1.0/v11_*v13_*v38_-1.0/v11_*v12_*v59_*v23_+v58_*v57_*1.0/v20_*v27_+c22*v57_*1.0/v15_*v12_*v39_-v57_*1.0/v20_*v5_*v60_*v27_+v58_*1.0/v20_*v5_*v59_*v27_-c22*v58_*v57_*1.0/v7_*v5_*v40_)+K*(c22/burgers-c22/v6_)+v19_*v59_*v41_-v58_*v60_*v42_+1.0/v7_*v13_*v60_*v28_-v58_*1.0/v7_*v16_*v24_+v58_*v57_*1.0/v20_*v29_+c22*v58_*1.0/v14_*v13_*v43_-v57_*1.0/v20_*v5_*v60_*v29_+v58_*1.0/v20_*v5_*v59_*v29_-c22*v58_*v57_*1.0/v7_*v5_*v44_;

        r_s[1][1]   = beta*(v19_*v61_*v33_+v58_*v62_*v37_+v57_*1.0/v11_*v13_*v38_-1.0/v11_*v12_*v61_*v23_-v58_*v57_*1.0/v20_*v27_-c11*v57_*1.0/v15_*v12_*v39_-v57_*1.0/v20_*v5_*v62_*v27_+v58_*1.0/v20_*v5_*v61_*v27_+c11*v58_*v57_*1.0/v7_*v5_*v40_)+K*(c11/burgers-c11/v6_)-v19_*v61_*v41_+v58_*v62_*v42_-1.0/v7_*v13_*v62_*v28_-v58_*1.0/v7_*v16_*v24_+v58_*v57_*1.0/v20_*v29_+c11*v58_*1.0/v14_*v13_*v43_+v57_*1.0/v20_*v5_*v62_*v29_-v58_*1.0/v20_*v5_*v61_*v29_-c11*v58_*v57_*1.0/v7_*v5_*v44_;

        r_s[0][1]   = v19_*v63_*v41_-K*(v48_/burgers-v48_/v6_)+v58_*v64_*v42_-beta*(v19_*v63_*v33_-v58_*v64_*v37_+v57_*1.0/v11_*v13_*v38_-1.0/v11_*v63_*v12_*v23_-v58_*v57_*1.0/v20_*v27_-c12*v57_*1.0/v15_*v12_*v38_+v57_*1.0/v20_*v5_*v64_*v27_+v58_*1.0/v20_*v63_*v5_*v27_+c12*v58_*v57_*1.0/v7_*v5_*v27_)-1.0/v7_*v13_*v64_*v28_+v58_*1.0/v7_*v16_*v24_-v58_*v57_*1.0/v20_*v29_-c12*v58_*1.0/v14_*v13_*v24_+v57_*1.0/v20_*v5_*v64_*v29_+v58_*1.0/v20_*v63_*v5_*v29_+c12*v58_*v57_*1.0/v7_*v5_*v29_;
    }


    // TODO get better name for this
    // Does this perhaps convert the reduced stress to "real" stress?
    Matrix2x2<double> Cell::riduci_matrix_reductions_stress(){
        
        // c=r_s
        Matrix2x2<double> d;

        d[0][0] = r_s[0][0]*m[0][0]*m[0][0] + r_s[0][1]*m[0][0]*m[0][1] + r_s[1][1]*m[0][1]*m[0][1];
        d[1][1] = r_s[0][0]*m[1][0]*m[1][0] + r_s[0][1]*m[1][0]*m[1][1] + r_s[1][1]*m[1][1]*m[1][1];
        d[0][1] = m[1][0]*(2 * r_s[0][0]*m[0][0] + r_s[0][1]*m[0][1]) + m[1][1]*(r_s[0][1]*m[0][0] + 2 * r_s[1][1]*m[0][1]);

        return d;

    }

    // What is Piola stress really? ... Is it a way to calculate "normal" stress?
    triangle Cell::calculate_Piola_stress(){
        Matrix2x2<double> temp = riduci_matrix_reductions_stress();
        Matrix2x2<double> temp = r_s.sym_orth_conjugate(m);
        // v1=deformation

        //Temp is a bit strange

        // What kind of triangle is this?
        triangle  t;
        // You can directly change the original triangle
        
        //derivative of energy with respect to f
        //P11 = phi,11
        t.a2->x =  (temp[0][0] * 2 * D[0][0] + temp[0][1] * D[1][0]);
        //P21 = phi,21
        t.a2->y =  (temp[0][0] * 2 * D[0][1] + temp[0][1] * D[1][1]);
        //P12 = phi,12
        t.a3->x =  (temp[1][1] * 2 * D[1][0] + temp[0][1] * D[0][0]);
        //P22 = phi,22
        t.a3->y =  (temp[1][1] * 2 * D[1][1] + temp[0][1] * D[0][1]);

        //forces  on the node 0     
        t.a1->x = t.a2->x - t.a3->x;
        t.a1->y = t.a2->y - t.a3->y;

        return t;
    }
};
/*
void calc_energy_forces(struct Grid& g,  alglib::real_1d_array &grad){
	
    // This cell variable is reused for each triangle
	Cell cell;
    // This is the total energy from all the triangles
    double total_energy;

    //#pragma omp parallel
	//#pragma omp for reduction(+:energy_thread)
	for (const triangle t : g.triangles) {

            cell = Cell(t, g.a);	
			
			total_energy += cell.energy;

			d =  cell.fptr4(cell.r_s, metrics_reduced.m);

            //-------FORCES ON  NODES-------		


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
    Singelton& s = Singelton::getInstance();
	calc_energy_forces(s.g, grad);
}

void run_simulation(){

    int nx, ny = 5;
    int n=nx*ny;

    Singelton& s = Singelton::getInstance();
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

    }

    
    
// 			std::cout<<"---"<<std::endl;
    //this copies the output of ALGLIB into the singletion c.p = y = x+u = c.pfix + starting_point
// 			transfer_from_alglib_to_std_v2(singleton.c,starting_point);

    std::copy(starting_point.getcontent(),
    starting_point.getcontent()+starting_point.length(), 
    starting_point_keep.getcontent());

    save_results(singleton.c,starting_point);


    save_elastic_state_gl2z(singleton.C_.current_metrics,current_metrics_t0);



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

*/