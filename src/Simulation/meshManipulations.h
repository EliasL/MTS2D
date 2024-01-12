#ifndef MESHMANIPULATIONS_H
#define MESHMANIPULATIONS_H
#pragma once

#include "settings.h"
#include "Matrix/matrix2x2.h"
#include "Mesh/mesh.h"


Matrix2x2<double> getShear(double load, double theta=0)
{
    // perturb is currently unused. If it will be used, it should be implemeted
    // propperly.
    double perturb = 0;

    Matrix2x2<double> trans;
    trans[0][0] = (1. - load * cos(theta + perturb) * sin(theta + perturb));
    trans[1][1] = (1. + load * cos(theta - perturb) * sin(theta - perturb));
    trans[0][1] = load * pow(cos(theta), 2.);
    trans[1][0] = -load * pow(sin(theta - perturb), 2.);

    return trans;
}

#endif