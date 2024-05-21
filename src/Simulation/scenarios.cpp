#include "scenarios.h"
#include "Data/dataExport.h"
#include "Data/paramParser.h"

// This function is quite complicated.
// We want to be able to prepare a simulation, and then initialize it, but if we
// have a loadedSimulation from a dump, then we want to skip both the
// preparation and the initialization. When you use this function, you want to
// prepare the simulation using a lambda function inside the function call. That
// is at least how I would do it.
std::shared_ptr<Simulation>
initOrLoad(Config config, std::string dataPath,
           std::shared_ptr<Simulation> loadedSimulation,
           std::function<void(std::shared_ptr<Simulation> s)> prepFunction) {
  std::shared_ptr<Simulation> s;
  if (loadedSimulation) {
    s = loadedSimulation; // Use the loaded simulation directly
  } else {
    s = std::make_shared<Simulation>(config, dataPath);
    prepFunction(s); // Call the initialization lambda
    s->initialize(); // Call initialize after the custom init function
  }
  return s;
}

void simpleShearFixedBoundary(Config config, std::string dataPath,
                              std::shared_ptr<Simulation> loadedSimulation) {
  // False here means that we do not have periodic boundary conditions
  // Boundary conditon transformation
  Matrix2d loadStepTransform = getShear(config.loadIncrement);
  Matrix2d startLoadTransform = getShear(config.startLoad);

  // This prepares a simulation OR loads an already started simulation
  // The last argument in this function is itself a function which makes some
  // changes to the simulation before it is initialized.
  // The initialization must happen after we have told the simulation what
  // nodes should be fixed or not for example
  auto s = initOrLoad(config, dataPath, loadedSimulation,
                      [startLoadTransform](std::shared_ptr<Simulation> s) {
                        // We need to fix the border nodes
                        s->mesh.fixBorderNodes();
                        // Prepare initial load condition
                        // Note that this transformation is applied to the
                        // ENTIRE mesh, not just the fixed nodes
                        s->mesh.applyTransformation(startLoadTransform); //
                      });

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

void simpleShear(Config config, std::string dataPath,
                 std::shared_ptr<Simulation> loadedSimulation) {
  Matrix2d loadStepTransform = getShear(config.loadIncrement);
  Matrix2d startLoadTransform = getShear(config.startLoad);

  auto s = initOrLoad(config, dataPath, loadedSimulation,
                      [startLoadTransform](std::shared_ptr<Simulation> s) {
                        // Prepare initial load condition
                        // Note that this transformation is applied to the
                        // ENTIRE mesh, not just the fixed nodes
                        s->mesh.applyTransformation(startLoadTransform);
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

void simpleShearWithNoise(Config config, std::string dataPath,
                          std::shared_ptr<Simulation> loadedSimulation) {
  Matrix2d loadStepTransform = getShear(config.loadIncrement);
  Matrix2d startLoadTransform = getShear(config.startLoad);

  auto s = initOrLoad(config, dataPath, loadedSimulation,
                      [startLoadTransform](std::shared_ptr<Simulation> s) {
                        // Prepare initial load condition
                        // Note that this transformation is applied to the
                        // ENTIRE mesh, not just the fixed nodes
                        s->mesh.applyTransformation(startLoadTransform);
                      });

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
                       std::shared_ptr<Simulation> loadedSimulation) {
  Matrix2d loadStepTransform = getShear(config.loadIncrement);
  Matrix2d startLoadTransform = getShear(config.startLoad);

  auto s = initOrLoad(config, dataPath, loadedSimulation,
                      [startLoadTransform](std::shared_ptr<Simulation> s) {
                        // Prepare initial load condition
                        // Note that this transformation is applied to the
                        // ENTIRE mesh, not just the fixed nodes
                        s->mesh.applyTransformation(startLoadTransform);
                      });

  while (s->mesh.load < s->maxLoad) {
    s->mesh.addLoad(s->loadIncrement);
    s->mesh.applyTransformationToSystemDeformation(loadStepTransform);

    // We keep loading until we reach extremes
    if (s->mesh.currentDeformation(0, 1) > 0.5 ||
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
                          std::shared_ptr<Simulation> loadedSimulation) {
  /*
  Note that the periodic boundary is not sheared! Only the fixed
  nodes are. The difference is that node (0,0) will be duplicated (as a ghost)
  to (0,n-1), not (alpha*(n-1), n-1). (You may think of two different alphas,
  one describing the loading of the fixed particles, and one describing the
  translation over the periodic boundary.)
  */

  Matrix2d loadStepTransform = getShear(config.loadIncrement);
  Matrix2d startLoadTransform = getShear(config.startLoad);

  auto s = initOrLoad(config, dataPath, loadedSimulation,
                      [startLoadTransform](std::shared_ptr<Simulation> s) {
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

void periodicBoundaryFixedComparisonTest(
    Config config, std::string dataPath,
    std::shared_ptr<Simulation> loadedSimulation) {
  Matrix2d loadStepTransform = getShear(config.loadIncrement);
  Matrix2d startLoadTransform = getShear(config.startLoad);

  auto s = initOrLoad(config, dataPath, loadedSimulation,
                      [startLoadTransform](std::shared_ptr<Simulation> s) {
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
                             std::shared_ptr<Simulation> loadedSimulation) {
  Matrix2d loadStepTransform = getShear(config.loadIncrement);
  Matrix2d startLoadTransform = getShear(config.startLoad);

  auto s = initOrLoad(
      config, dataPath, loadedSimulation,
      [startLoadTransform](std::shared_ptr<Simulation> s) {
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
                                std::shared_ptr<Simulation> loadedSimulation) {
  /*
  Note that the periodic boundary is not sheared! Only the fixed
  nodes are. The difference is that node (0,0) will be duplicated (as a ghost)
  to (0,n-1), not (alpha*(n-1), n-1). (You may think of two different alphas,
  one describing the loading of the fixed particles, and one describing the
  translation over the periodic boundary.)
  */

  Matrix2d loadStepTransform = getShear(config.loadIncrement);
  Matrix2d startLoadTransform = getShear(config.startLoad);

  auto s = initOrLoad(config, dataPath, loadedSimulation,
                      [startLoadTransform](std::shared_ptr<Simulation> s) {
                        s->mesh.applyTransformation(startLoadTransform);
                      });

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

void handleInputArgs(int argc, char *argv[]) {
  // If we are given less than two arguments, it means that there is no config
  // file or simulation dump
  if (argc < 2) {
    std::cerr << "Error! Usage: " << argv[0]
              << " <Config File / Dump.mtsb> [Output Path]\n"
              << " - Configuration file NOT provided!\n"
              << " - Simulation dump NOT provided!\n";
    return;
  }

  std::string arg1 = argv[1];
  // Lets check if we are given a dump instead of a config file
  if (arg1.size() >= 5 && arg1.substr(arg1.size() - 5) == ".mtsb") {
    auto sPtr = std::make_shared<Simulation>();

    if (argc > 2) {
      std::string configPath = argv[2];
      std::cout << "Overwriting simulation settings using\n - " << configPath
                << '\n';
      // Load directly into the shared_ptr managed Simulation, using new
      // settings
      Simulation::loadSimulation(*sPtr, arg1, configPath);
    } else {
      // Load directly into the shared_ptr managed Simulation
      Simulation::loadSimulation(*sPtr, arg1);
    }
    std::cout << "Resuming simulation using " << arg1 << '\n'
              << " - Config File: " << sPtr->config.name << '\n'
              << " - Data Path: " << sPtr->dataPath << '\n'
              << " - Current Load: " << sPtr->mesh.load << '\n';
    runSimulationScenario(sPtr->config, sPtr->dataPath, sPtr);
  } else {
    // We extract the configFile and dataPath we will be storing the data in.
    std::string configPath = argv[1];
    Config config = parseConfigFile(configPath);
    // If there is no given dataPath, we find one.
    std::string dataPath = (argc >= 3) ? argv[2] : findOutputPath();

    std::cout << "Running simulation with:\n"
              << " - Config File: " << configPath << '\n'
              << " - Data Path: " << dataPath << '\n';
    runSimulationScenario(config, dataPath);
  }
}

void runSimulationScenario(Config config, std::string dataPath,
                           std::shared_ptr<Simulation> loadedSimulation) {
  // Here we have a map associating each scenario with a string you can use
  // in the config file
  std::unordered_map<std::string,
                     std::function<void(const Config &, const std::string &,
                                        std::shared_ptr<Simulation>)>>
      scenarioMap = {
          {"simpleShearFixedBoundary", simpleShearFixedBoundary},
          {"simpleShear", simpleShear},
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
