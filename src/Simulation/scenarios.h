#ifndef SCENARIOS_H
#define SCENARIOS_H
#include <omp.h>
#pragma once

#include "../Data/data_export.h"
#include "simulation.h"

void handleInputArgs(int argc, char *argv[]);
void runSimulationScenario(
    Config config, std::string dataPath,
    std::shared_ptr<Simulation> loadedSimulation = nullptr);

// Helper functions

// This function is quite complicated.
// We want to be able to prepare a simulation, and then initialize it, but if we
// have a loadedSimulation from a dump, then we want to skip both the
// preparation and the initialization. When you use this function, you want to
// prepare the simulation using a lambda function inside the function call.
std::shared_ptr<Simulation>
initOrLoad(Config config, std::string dataPath,
           std::shared_ptr<Simulation> loadedSimulation,
           std::function<void(std::shared_ptr<Simulation> s)> prepFunction);

// This prepares a simulation OR loads an already started simulation
// It fixes the border nodes and applies the start transformation.
std::shared_ptr<Simulation>
getFixedBorderSimulation(Config config, std::string dataPath,
                         std::shared_ptr<Simulation> loadedSimulation);

// This prepares a simulation OR loads an already started simulation
// The prefFunctino only applies the start transformation.
std::shared_ptr<Simulation>
getPeriodicBorderSimulation(Config config, std::string dataPath,
                            std::shared_ptr<Simulation> loadedSimulation);

#endif