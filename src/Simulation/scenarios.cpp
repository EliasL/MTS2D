#include "scenarios.h"

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

        // Modifies the nodeDisplacements
        s.setInitialGuess(loadStepTransform);

        // If it is the first step of the simulation
        if (load == s.startLoad)
        {
            s.addNoiseToGuess();
        }
        // Minimizes the energy by moving the positions of the free nodes in the mesh
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

        // Modifies the nodeDisplacements
        s.setInitialGuess(loadStepTransform);

        // If it is the first step of the simulation
        if (load == s.startLoad)
        {
            s.addNoiseToGuess();
        }
        // Minimizes the energy by moving the positions of the free nodes in the mesh
        s.minimize_with_alglib();

        // Updates progress and writes to file
        s.finishStep(load);
    }
    s.finishSimulation();
}

void periodicBoundaryTest(Config config, std::string dataPath)
{
    /*
    Note that the periodic boundary is not sheared! Only the fixed
    nodes are. The difference is that node (0,0) will be duplicated (as a ghost)
    to (0,n-1), not (alpha*(n-1), n-1). (You may think of two different alphas,
    one describing the loading of the fixed particles, and one describing the
    translation over the periodic boundary.)
    */

    // False here means that we do not have periodic boundary conditions
    Simulation s = Simulation(config, dataPath, true);

    // We fix two of the rows
    s.mesh.fixNodesInRow(0);
    int fixedMiddleRow = std::floor(s.rows / 2);
    s.mesh.fixNodesInRow(fixedMiddleRow);

    // We also fix the first column so that we can compare with fixed boundaries
    // later
    s.mesh.fixNodesInColumn(0);

    s.initialize();
    // Boundary conditon transformation
    Matrix2x2<double> loadStepTransform = getShear(s.loadIncrement);

    for (double load = s.startLoad; load <= s.maxLoad; load += s.loadIncrement)
    {
        // Moves the fixed nodes
        s.mesh.applyTransformationToFixedNodes(loadStepTransform);
        s.mesh.applyTransformationToSystemDeformation(loadStepTransform);

        // Modifies the nodeDisplacements
        s.setInitialGuess(loadStepTransform);

        // Minimizes the energy by moving the positions of the free nodes in the mesh
        s.minimize_with_alglib();

        // Updates progress and writes to file
        s.finishStep(load);
    }
    s.finishSimulation();
}

void periodicBoundaryFixedComparisonTest(Config config, std::string dataPath)
{
    // Now we try to do the same as in the periodic boundary test, but with fixed
    // boundaries
    // False here means that we do not have periodic boundary conditions
    Simulation s = Simulation(config, dataPath, false);

    // Boundary conditon transformation
    Matrix2x2<double> loadStepTransform = getShear(s.loadIncrement);

    s.mesh.fixBorderNodes();
    int fixedMiddleRow = std::floor(s.rows / 2);
    s.mesh.fixNodesInRow(fixedMiddleRow);

    s.initialize();

    for (double load = s.startLoad; load <= s.maxLoad; load += s.loadIncrement)
    {
        // Moves the fixed nodes
        s.mesh.applyTransformationToFixedNodes(loadStepTransform);
        s.mesh.applyTransformationToSystemDeformation(loadStepTransform);

        s.setInitialGuess(loadStepTransform);

        // Minimizes the energy by moving the positions of the free nodes in the mesh
        s.minimize_with_alglib();

        // Updates progress and writes to file
        s.finishStep(load);
    }
    s.finishSimulation();
}

void failedSingleDislocation(Config config, std::string dataPath)
{
    // Now we try to do the same as in the periodic boundary test, but with fixed
    // boundaries
    // False here means that we do not have periodic boundary conditions
    Simulation s = Simulation(config, dataPath, false);

    // Boundary conditon transformation
    Matrix2x2<double> loadStepTransform = getShear(s.loadIncrement);

    int middleRow = std::floor(s.rows / 2);
    int middleCol = std::floor(s.cols / 2);
    int distance = 10;
    int radius = 5; // half width

    // We want to add a single dislocation
    // We do that by gradually shifting a section one unit length
    for (int row = middleRow; row < s.rows; row++)
    {
        for (int col = std::max(middleCol - radius, 0); col < s.cols; col++)
        {
            Vector2d disp = {std::clamp(static_cast<double>(col - (middleCol - radius)) / (radius * 2), 0.0, 1.0) * s.mesh.a, 0.0};
            s.mesh.nodes[row][col].setDisplacement(disp);
        }
    }
    s.initialize();
    s.finishStep(-1);
    for (double load = s.startLoad; load <= s.maxLoad; load += s.loadIncrement)
    {
        s.mesh.applyTransformationToSystemDeformation(loadStepTransform);

        s.setInitialGuess(loadStepTransform);

        // Minimizes the energy by moving the positions of the free nodes in the mesh
        s.minimize_with_alglib();

        // Updates progress and writes to file
        s.finishStep(load);
    }
    s.finishSimulation();
}

void runSimulationScenario(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Error! Usage: " << argv[0] << " <Config File> [Output Path]\n"
                  << " - Configuration file NOT provided!\n";
        return;
    }

    std::string configFile = argv[1];
    std::string dataPath = (argc >= 3) ? argv[2] : findOutputPath();

    std::cout << "Running simulation with:\n"
              << " - Config File: " << configFile << '\n'
              << " - Data Path: " << dataPath << '\n';

    Config config = getConf(configFile);
    // Here we have a map associating each scenario with a string you can use
    // in the config file
    std::unordered_map<std::string, std::function<void(const Config &, const std::string &)>> scenarioMap = {
        {"simpleShearFixedBoundary", simpleShearFixedBoundary},
        {"simpleShearPeriodicBoundary", simpleShearPeriodicBoundary},
        {"periodicBoundaryTest", periodicBoundaryTest},
        {"periodicBoundaryFixedComparisonTest", periodicBoundaryFixedComparisonTest},
        {"failedSingleDislocation", failedSingleDislocation},
    };

    // We now search though the map until we find a match.
    auto it = scenarioMap.find(config.scenario);
    // If we didn't find anything, we throw an error.
    if (it != scenarioMap.end())
    {
        it->second(config, dataPath);
    }
    else
    {
        std::cerr << "No matching scenario!\n";
    }
}
