#ifndef SETTINGS_H
#define SETTINGS_H
#pragma once

// This determines how the lagrange_reduction function behaves
// (And some other things) TODO: Understand better
#define LINEARITY true



enum class MetricFunction {faicella, epsilon_lineaire};
// This determines what function is used to calculate metrics in cells
#define METRICFUNCTION MetricFunction::faicella


enum class BoundaryConditionFunction {macro_shear};
// This determines what function is used to set the forces on the boundary nodes
#define BOUNDARYCONDITIONFUNCTION BoundaryConditionFunction::macro_shear;

#endif