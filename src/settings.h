#ifndef SETTINGS_H
#define SETTINGS_H
#pragma once

/*
SETTINGS
*/

// This determines how the lagrange_reduction function behaves
// (And some other things) TODO: Understand better
#define LINEARITY true

// This determines what function is used to calculate metrics in cells
#define METRICFUNCTION MetricFunction::faicella

// This determines what function is used to set the forces on the boundary nodes
#define BOUNDARYCONDITIONFUNCTION BoundaryConditionFunction::macro_shear;

/*
DEFINITIONS
*/
enum class MetricFunction
{
    faicella,        // TODO, describe better
    epsilon_lineaire // TODO, describe better
};
enum class BoundaryConditionFunction
{
    macro_shear // TODO, describe better 
};

#endif