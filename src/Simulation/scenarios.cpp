#include "scenarios.h"
#include "Data/dataExport.h"
#include "Data/paramParser.h"
#include "Simulation/simulation.h"
#include <iostream>
#include <memory>
#include <unistd.h>

using SimPtr = std::shared_ptr<Simulation>;

void simpleShear(Config config, std::string dataPath, SimPtr loadedSimulation) {
  Matrix2d loadStepTransform = getShear(config.loadIncrement);

  SimPtr s = initOrLoadPeriodic(config, dataPath, loadedSimulation);

  while (s->mesh.load < s->maxLoad) {
    s->mesh.addLoad(s->loadIncrement);
    s->mesh.applyTransformationToSystemDeformation(loadStepTransform);

    // Modifies the nodeDisplacements
    s->setInitialGuess(loadStepTransform);

    // If it is the first step of the simulation
    if (s->mesh.loadSteps == 1) {
      s->addNoiseToGuess();
    }
    // Minimizes the energy by moving the positions of the free nodes in the
    // mesh
    s->minimize();
    // Updates progress and writes to file
    s->finishStep();
  }
  s->finishSimulation();
}

void simpleShearFixedBoundary(Config config, std::string dataPath,
                              SimPtr loadedSimulation) {
  Matrix2d loadStepTransform = getShear(config.loadIncrement);

  SimPtr s = initOrLoadFixed(config, dataPath, loadedSimulation);

  while (s->mesh.load < s->maxLoad) {
    s->mesh.addLoad(s->loadIncrement);
    s->mesh.applyTransformationToFixedNodes(loadStepTransform);

    // Modifies the nodeDisplacements
    s->setInitialGuess(loadStepTransform);

    // If it is the first step of the simulation
    if (s->mesh.loadSteps == 1) {
      s->addNoiseToGuess();
    }
    // Minimizes the energy by moving the positions of the free nodes in the
    // mesh
    s->minimize();
    // Updates progress and writes to file
    s->finishStep();
  }
  s->finishSimulation();
}

void simpleShearWithNoise(Config config, std::string dataPath,
                          SimPtr loadedSimulation) {
  Matrix2d loadStepTransform = getShear(config.loadIncrement);

  SimPtr s = initOrLoadPeriodic(config, dataPath, loadedSimulation);

  while (s->mesh.load < s->maxLoad) {
    s->mesh.addLoad(s->loadIncrement);
    s->mesh.applyTransformationToSystemDeformation(loadStepTransform);

    // Modifies the nodeDisplacements
    s->setInitialGuess(loadStepTransform);

    // If it is the first step of the simulation
    if (s->mesh.loadSteps == 1) {
      s->addNoiseToGuess();
    } else {
      s->addNoiseToGuess(0.0000008);
    }
    // Minimizes the energy by moving the positions of the free nodes in the
    // mesh
    s->minimize();
    // Updates progress and writes to file
    s->finishStep();
  }
  s->finishSimulation();
}

void cyclicSimpleShear(Config config, std::string dataPath,
                       SimPtr loadedSimulation) {
  Matrix2d loadStepTransform = getShear(config.loadIncrement);

  SimPtr s = initOrLoadPeriodic(config, dataPath, loadedSimulation);

  while (s->mesh.load < s->maxLoad) {
    s->mesh.addLoad(s->loadIncrement);
    s->mesh.applyTransformationToSystemDeformation(loadStepTransform);

    // We keep loading until we reach extremes
    if (s->mesh.currentDeformation(0, 1) > 0.2 ||
        s->mesh.currentDeformation(0, 1) < 0.0) {
      loadStepTransform(0, 1) *= -1;
      // Now we take a step to go back to where we were, and then another one to
      // solve for
      s->mesh.applyTransformationToSystemDeformation(loadStepTransform);
      s->mesh.applyTransformationToSystemDeformation(loadStepTransform);
    }

    // Modifies the nodeDisplacements
    s->setInitialGuess(loadStepTransform);

    // If it is the first step of the simulation
    if (s->mesh.loadSteps == 1) {
      s->addNoiseToGuess();
    }
    // Minimizes the energy by moving the positions of the free nodes in the
    // mesh
    s->minimize();

    // Updates progress and writes to file
    s->finishStep();
  }
  s->finishSimulation();
}

void periodicBoundaryTest(Config config, std::string dataPath,
                          SimPtr loadedSimulation) {
  /*
  Note that the periodic boundary is not sheared! Only the fixed
  nodes are. The difference is that node (0,0) will be duplicated (as a ghost)
  to (0,n-1), not (alpha*(n-1), n-1). (You may think of two different alphas,
  one describing the loading of the fixed particles, and one describing the
  translation over the periodic boundary.)
  */

  Matrix2d loadStepTransform = getShear(config.loadIncrement);
  Matrix2d startLoadTransform = getShear(config.startLoad);

  SimPtr s = initOrLoad(config, dataPath, loadedSimulation,
                        [startLoadTransform](SimPtr s) {
                          // We fix two of the rows
                          s->mesh.fixNodesInRow(0);
                          int fixedMiddleRow = std::floor(s->rows / 2);
                          s->mesh.fixNodesInRow(fixedMiddleRow);

                          // We also fix the first column so that we can compare
                          // with fixed boundaries later
                          s->mesh.fixNodesInColumn(0);
                        });

  s->writeToFile(true);
  while (s->mesh.load < s->maxLoad) {
    s->mesh.addLoad(s->loadIncrement);
    // Moves the fixed nodes
    s->mesh.applyTransformationToFixedNodes(loadStepTransform);
    s->mesh.applyTransformationToSystemDeformation(loadStepTransform);

    // Modifies the nodeDisplacements
    s->setInitialGuess(loadStepTransform);

    // Minimizes the energy by moving the positions of the free nodes in the
    // mesh
    s->minimize();

    // Updates progress and writes to file
    s->finishStep();
  }
  s->finishSimulation();
}

void periodicBoundaryFixedComparisonTest(Config config, std::string dataPath,
                                         SimPtr loadedSimulation) {
  Matrix2d loadStepTransform = getShear(config.loadIncrement);
  Matrix2d startLoadTransform = getShear(config.startLoad);

  SimPtr s = initOrLoad(config, dataPath, loadedSimulation,
                        [startLoadTransform](SimPtr s) {
                          s->mesh.fixBorderNodes();
                          int fixedMiddleRow = std::floor(s->rows / 2);
                          s->mesh.fixNodesInRow(fixedMiddleRow);
                        });

  while (s->mesh.load < s->maxLoad) {
    s->mesh.addLoad(s->loadIncrement);
    s->mesh.applyTransformationToSystemDeformation(loadStepTransform);

    // Modifies the nodeDisplacements
    s->setInitialGuess(loadStepTransform);

    // If it is the first step of the simulation
    if (s->mesh.loadSteps == 1) {
      s->addNoiseToGuess();
    }
    // Minimizes the energy by moving the positions of the free nodes in the
    // mesh
    s->minimize();

    // Updates progress and writes to file
    s->finishStep();
  }
  s->finishSimulation();
}

void failedSingleDislocation(Config config, std::string dataPath,
                             SimPtr loadedSimulation) {
  Matrix2d loadStepTransform = getShear(config.loadIncrement);
  Matrix2d startLoadTransform = getShear(config.startLoad);

  SimPtr s = initOrLoad(
      config, dataPath, loadedSimulation, [startLoadTransform](SimPtr s) {
        int middleRow = std::floor(s->rows / 2);
        int middleCol = std::floor(s->cols / 2);
        int radius = 5; // half width

        // We want to add a single dislocation
        // We do that by gradually shifting a section one unit length
        for (int row = middleRow; row < s->rows; row++) {
          for (int col = std::max(middleCol - radius, 0); col < s->cols;
               col++) {
            Vector2d disp = {
                std::clamp(static_cast<double>(col - (middleCol - radius)) /
                               (radius * 2),
                           0.0, 1.0) *
                    s->mesh.a,
                0.0};
            s->mesh.nodes(row, col).setDisplacement(disp);
          }
        }
      });

  writeMeshToVtu(s->mesh, s->name, s->dataPath);
  while (s->mesh.load < s->maxLoad) {
    s->mesh.addLoad(s->loadIncrement);
    s->mesh.applyTransformationToSystemDeformation(loadStepTransform);

    s->setInitialGuess(loadStepTransform);

    // Minimizes the energy by moving the positions of the free nodes in the
    // mesh
    s->minimize();

    // Updates progress and writes to file
    s->finishStep();
  }
  s->finishSimulation();
}

void createDumpBeforeEnergyDrop(Config config, std::string dataPath,
                                SimPtr loadedSimulation) {
  /*
  Note that the periodic boundary is not sheared! Only the fixed
  nodes are. The difference is that node (0,0) will be duplicated (as a ghost)
  to (0,n-1), not (alpha*(n-1), n-1). (You may think of two different alphas,
  one describing the loading of the fixed particles, and one describing the
  translation over the periodic boundary.)
  */

  Matrix2d loadStepTransform = getShear(config.loadIncrement);

  SimPtr s = initOrLoadPeriodic(config, dataPath, loadedSimulation);

  std::string dumps[] = {"dump1", "dump2"};
  int dumpInUse = 0;
  int saveInterval = 10;
  double lastEnergy = 0;
  int step = 0;
  while (s->mesh.load < s->maxLoad) {
    s->mesh.addLoad(s->loadIncrement);
    s->mesh.applyTransformationToSystemDeformation(loadStepTransform);

    // Modifies the nodeDisplacements
    s->setInitialGuess(loadStepTransform);

    // If it is the first step of the simulation
    if (s->mesh.loadSteps == 1) {
      s->addNoiseToGuess();
    }

    // Minimizes the energy by moving the positions of the free nodes in the
    // mesh
    s->minimize();

    // Updates progress and writes to file
    s->finishStep();

    if (step == 0) {
      lastEnergy = s->mesh.averageEnergy;
    }

    // Keep always one updated dump
    step++;
    if (step % saveInterval == 0) {
      s->saveSimulation(dumps[dumpInUse]);
      if (step % (saveInterval * 2) == 0) {
        dumpInUse = (dumpInUse + 1) % 2;
      }
    }

    // Check if we have had a big drop
    if (lastEnergy - s->mesh.averageEnergy < -0.0004) {
      s->saveSimulation(dumps[dumpInUse] + "_EnergyFall");
      break;
    }
    lastEnergy = s->mesh.averageEnergy;
  }
  s->finishSimulation();
}

std::string trim(const std::string &str) {
  size_t first = str.find_first_not_of(' ');
  if (std::string::npos == first) {
    return str;
  }
  size_t last = str.find_last_not_of(' ');
  return str.substr(first, (last - first + 1));
}

void handleInputArgs(int argc, char *argv[]) {
  std::string configPath;
  std::string dumpPath;
  std::string outputPath;
  // The program will check if the folder already contains a completed
  // simulation If the simulation is complete, the program will terminate and
  // not rerun the simulation unless forceReRun is true
  bool forceReRun = false;

  int opt;
  while ((opt = getopt(argc, argv, "c:d:o:r")) != -1) {
    switch (opt) {
    case 'c':
      configPath = trim(optarg);
      break;
    case 'd':
      dumpPath = trim(optarg);
      break;
    case 'o':
      outputPath = trim(optarg);
      break;
    case 'r':
      forceReRun = true;
      break;
    default:
      std::cerr << "Usage: " << argv[0]
                << " -c <Config File> -d <Dump File> [-o <Output Path>] [-r]\n";
      return;
    }
  }

  // Check if dumpPath is provided
  if (!dumpPath.empty()) {
    auto sPtr = std::make_shared<Simulation>(); // Create a shared pointer to a
                                                // new Simulation object

    if (!configPath.empty()) { // Check if configPath is provided
      std::cout << "Overwriting simulation settings using\n - " << configPath
                << '\n';
    } else {
      configPath = searchForConfig(
          dumpPath); // Try to find configPath in the same folder as dumpPath
      if (configPath.empty()) { // If no configPath is found
        std::cerr << "Error! No config provided or found in the same folder as "
                     "the dump file.\n";
        return; // Exit the function
      }
    }

    // Load and resume the simulation using the provided or found configPath
    Simulation::loadSimulation(*sPtr, dumpPath, configPath, forceReRun);
    std::cout << "Resuming simulation using " << dumpPath << '\n'
              << " - Config File: " << sPtr->config.name << '\n'
              << " - Data Path: " << sPtr->dataPath << '\n'
              << sPtr->config << '\n'
              << " - Current Load: " << sPtr->mesh.load << '\n';
    runSimulationScenario(sPtr->config, sPtr->dataPath,
                          sPtr); // Run the simulation scenario

    // If dumpPath is not provided but configPath is
  } else if (!configPath.empty()) {
    Config config = parseConfigFile(configPath); // Parse the configuration file
    config.forceReRun = forceReRun; // Set the forceReRun flag in the config

    if (outputPath.empty()) {        // If outputPath is not provided
      outputPath = findOutputPath(); // Find the output path
    }

    // Run the simulation with the provided configuration
    std::cout << "Running simulation with:\n"
              << " - Config File: " << configPath << '\n'
              << " - Data Path: " << outputPath << '\n'
              << config << '\n';
    runSimulationScenario(config, outputPath); // Run the simulation scenario
  }
}

void runSimulationScenario(Config config, std::string dataPath,
                           SimPtr loadedSimulation) {
  // Here we have a map associating each scenario with a string you can use
  // in the config file
  std::unordered_map<
      std::string,
      std::function<void(const Config &, const std::string &, SimPtr)>>
      scenarioMap = {
          {"simpleShear", simpleShear},
          {"simpleShearFixedBoundary", simpleShearFixedBoundary},
          {"simpleShearWithNoise", simpleShearWithNoise},
          {"periodicBoundaryTest", periodicBoundaryTest},
          {"periodicBoundaryFixedComparisonTest",
           periodicBoundaryFixedComparisonTest},
          {"failedSingleDislocation", failedSingleDislocation},
          {"cyclicSimpleShear", cyclicSimpleShear},
          {"createDumpBeforeEnergyDrop", createDumpBeforeEnergyDrop},
      };

  // We now search though the map until we find a match.
  auto it = scenarioMap.find(config.scenario);
  // If we didn't find anything, we throw an error.
  if (it != scenarioMap.end()) {
    it->second(config, dataPath, loadedSimulation);
  } else {
    std::cerr << "No matching scenario!\n";
  }
}

// This function is quite complicated.
// We want to be able to prepare a simulation, and then initialize it, but if we
// have a loadedSimulation from a dump, then we want to skip both the
// preparation and the initialization. When you use this function, you want to
// prepare the simulation using a lambda function inside the function call.
SimPtr initOrLoad(Config config, std::string dataPath, SimPtr loadedSimulation,
                  std::function<void(SimPtr s)> prepFunction) {
  SimPtr s;
  if (loadedSimulation) {
    s = loadedSimulation; // Use the loaded simulation directly
  } else {
    s = std::make_shared<Simulation>(config, dataPath);
    prepFunction(s); // Call the initialization lambda
    s->initialize(); // Call initialize after the custom init function
  }
  return s;
}
SimPtr initOrLoadFixed(Config config, std::string dataPath,
                       SimPtr loadedSimulation) {
  auto startLoadTransform = getShear(config.startLoad);
  return initOrLoad(config, dataPath, loadedSimulation,
                    [startLoadTransform](SimPtr s) {
                      s->mesh.fixBorderNodes();
                      s->mesh.applyTransformation(startLoadTransform);
                    });
}

SimPtr initOrLoadPeriodic(Config config, std::string dataPath,
                          SimPtr loadedSimulation) {
  auto startLoadTransform = getShear(config.startLoad);
  return initOrLoad(config, dataPath, loadedSimulation,
                    [startLoadTransform](SimPtr s) {
                      // Assuming some periodic boundary-specific operations
                      s->mesh.applyTransformation(startLoadTransform);
                    });
}
// Utility function to setup transformations based on the config
std::pair<Matrix2d, Matrix2d> getTransformations(const Config &config) {
  return {getShear(config.loadIncrement), getShear(config.startLoad)};
}