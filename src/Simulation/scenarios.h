#ifndef SCENARIOS_H
#define SCENARIOS_H
#pragma once

#include "Data/dataExport.h"
#include "simulation.h"

void handleInputArgs(int argc, char *argv[]);
void runSimulationScenario(
    Config config, std::string dataPath,
    std::shared_ptr<Simulation> loadedSimulation = nullptr);

#endif