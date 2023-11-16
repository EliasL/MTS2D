#ifndef SETTINGS_H
#define SETTINGS_H
#pragma once


/*
SETTINGS
*/

// Uncomment this define to disable logging
//#define NO_LOG
#define LOG_LEVEL DEBUG_LEVEL
#define LOG_DETAIL MINIMAL
// It is important that this include comes after these definitions
#include "Utility/macrologger.h"

// This determines how the lagrangeReduction function behaves
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


/**
 * @brief A Boundary Condition Function.
 * 
 * Determines how the transformation matrix to deform the boundary nodes is calculated
 */
enum class BCF
{
    macroShear // TODO, describe better 
};

#endif