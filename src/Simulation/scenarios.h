#ifndef SCENARIOS_H
#define SCENARIOS_H
#pragma once

#include "simulation.h"
#include "Data/dataExport.h"

void simpleShearFixedBoundary(Config config, std::string dataPath)
{
    // False here means that we do not have periodic boundary conditions
    Simulation s = Simulation(config, dataPath, false);

    // Boundary conditon transformation
    Matrix2x2<double> loadStepTransform = getShear(s.loadIncrement);
    Matrix2x2<double> startLoadTransform = getShear(s.startLoad);

    // We need to fix the border nodes
    s.mesh.fixBorderNodes();
    // Prepare initial load condition
    // Note that this transformation is applied to the ENTIRE mesh, not just
    // the fixed nodes
    s.mesh.applyTransformation(startLoadTransform);

    s.initialize();

    for (double load = s.startLoad; load <= s.maxLoad; load += s.loadIncrement)
    {
        s.mesh.applyTransformationToFixedNodes(loadStepTransform);

        s.setInitialGuess(loadStepTransform);

        // If it is the first step of the simulation
        if (load == s.startLoad)
        {
            s.addNoiseToGuess();
        }
        // minimizes the energy by moving the positions of the nodes in the mesh
        s.minimize_with_alglib();

        // Updates progress and writes to file
        s.finishStep(load);
    }
    s.finishSimulation();
}

void simpleShearPeriodicBoundary(Config config, std::string dataPath)
{
    // False here means that we do not have periodic boundary conditions
    Simulation s = Simulation(config, dataPath, true);

    // Boundary conditon transformation
    Matrix2x2<double> loadStepTransform = getShear(s.loadIncrement);
    Matrix2x2<double> startLoadTransform = getShear(s.startLoad);

    s.mesh.applyTransformation(startLoadTransform);

    s.initialize();

    for (double load = s.startLoad; load <= s.maxLoad; load += s.loadIncrement)
    {
        s.mesh.applyTransformationToSystemDeformation(loadStepTransform);

        s.setInitialGuess(loadStepTransform);

        // If it is the first step of the simulation
        if (load == s.startLoad)
        {
            s.addNoiseToGuess();
        }
        // minimizes the energy by moving the positions of the nodes in the mesh
        s.minimize_with_alglib();

        // Updates progress and writes to file
        s.finishStep(load);
    }
    s.finishSimulation();
}

void runSimunationScenario(int argc, char *argv[])
{
    // Check for minimum number of arguments upfront
    if (argc < 2)
    {
        std::cerr << "Error! Usage: " << argv[0] << " <Config File> [Output Path]" << '\n';
        if (argc < 2)
            std::cerr << " - Configuration file NOT provided!" << '\n';
        return; // Exit the function if required arguments are not provided
    }

    std::string configFile = argv[1];
    std::string dataPath = (argc >= 3) ? argv[2] : findOutputPath(); // Use findOutputPath as default if not specified

    // Proceed with the simulation setup using the parsed arguments
    std::cout << "Running simulation with:" << '\n';
    std::cout << " - Config File: " << configFile << '\n';
    std::cout << " - Data Path: " << dataPath << '\n';

    Config config = getConf(configFile);

    // Choose scenario
    if (config.scenario == "simpleShearFixedBoundary")
    {
        simpleShearFixedBoundary(config, dataPath);
    }
    else if (config.scenario == "simpleShearPeriodicBoundary")
    {
        simpleShearPeriodicBoundary(config, dataPath);
    }
}

#endif