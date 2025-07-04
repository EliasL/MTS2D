#include "scenarios.h"
#include "Data/data_export.h"
#include "Data/param_parser.h"
#include "Eigen/Core"
#include "Simulation/simulation.h"
#include <cassert>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <unistd.h>

using SimPtr = std::shared_ptr<Simulation>;

void simpleShear(Config config, std::string dataPath, SimPtr loadedSimulation) {
  Matrix2d loadStepTransform = getShear(config.loadIncrement);

  SimPtr s = getPeriodicBorderSimulation(config, dataPath, loadedSimulation);

  while (s->keepLoading()) {
    s->mesh.addLoad(s->loadIncrement);
    s->mesh.applyTransformationToSystemDeformation(loadStepTransform);

    // Modifies the nodeDisplacements
    s->setInitialGuess(loadStepTransform);

    // Minimizes the energy by moving the free nodes in the mesh
    s->minimize();

    // Updates progress and writes to file
    s->finishStep(true);
  }
  s->finishSimulation();
}

void simpleShearFixedBoundary(Config config, std::string dataPath,
                              SimPtr loadedSimulation) {
  Matrix2d loadStepTransform = getShear(config.loadIncrement);

  SimPtr s = getFixedBorderSimulation(config, dataPath, loadedSimulation);

  while (s->keepLoading()) {
    s->mesh.addLoad(s->loadIncrement);
    s->mesh.applyTransformationToFixedNodes(loadStepTransform);

    // Modifies the nodeDisplacements
    s->setInitialGuess(loadStepTransform);

    // Minimizes the energy by moving the free nodes in the mesh
    s->minimize();
    // Updates progress and writes to file
    s->finishStep();
  }
  s->finishSimulation();
}

void simpleShearWithNoise(Config config, std::string dataPath,
                          SimPtr loadedSimulation) {
  Matrix2d loadStepTransform = getShear(config.loadIncrement);

  SimPtr s = getPeriodicBorderSimulation(config, dataPath, loadedSimulation);

  while (s->keepLoading()) {
    s->mesh.addLoad(s->loadIncrement);
    s->mesh.applyTransformationToSystemDeformation(loadStepTransform);

    // Modifies the nodeDisplacements
    s->setInitialGuess(loadStepTransform);

    s->addNoiseToGuess(0.0000008);
    // Minimizes the energy by moving the free nodes in the mesh
    s->minimize();
    // Updates progress and writes to file
    s->finishStep();
  }
  s->finishSimulation();
}

void cyclicSimpleShear(Config config, std::string dataPath,
                       SimPtr loadedSimulation) {
  Matrix2d loadStepTransform = getShear(config.loadIncrement);

  SimPtr s = getPeriodicBorderSimulation(config, dataPath, loadedSimulation);

  while (s->keepLoading()) {
    s->mesh.addLoad(s->loadIncrement);
    s->mesh.applyTransformationToSystemDeformation(loadStepTransform);

    // We keep loading until we reach extremes
    double shear = s->mesh.currentDeformation(0, 1);
    if ((shear > 0.3 && s->loadIncrement > 0) ||
        (shear < 0.16 && s->loadIncrement < 0)) {

      loadStepTransform(0, 1) *= -1;
      s->loadIncrement *= -1;
      // Now we take a step to go back to where we were, and then another one to
      // solve for
      s->mesh.load += 2 * s->loadIncrement;
      s->mesh.applyTransformationToSystemDeformation(loadStepTransform);
      s->mesh.applyTransformationToSystemDeformation(loadStepTransform);
    }

    // Modifies the nodeDisplacements
    s->setInitialGuess(loadStepTransform);

    // Minimizes the energy by moving the free nodes in the mesh
    s->minimize(false);

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

  while (s->keepLoading()) {
    s->mesh.addLoad(s->loadIncrement);
    // Moves the fixed nodes
    s->mesh.applyTransformationToFixedNodes(loadStepTransform);
    s->mesh.applyTransformationToSystemDeformation(loadStepTransform);

    // Modifies the nodeDisplacements
    s->setInitialGuess(loadStepTransform);

    // Minimizes the energy by moving the free nodes in the mesh
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
  // TODO looks like startLoadTransform is not used, so probably something here
  // doesn't work
  SimPtr s = initOrLoad(config, dataPath, loadedSimulation,
                        [startLoadTransform](SimPtr s) {
                          s->mesh.fixBorderNodes();
                          int fixedMiddleRow = std::floor(s->rows / 2);
                          s->mesh.fixNodesInRow(fixedMiddleRow);
                        });

  while (s->keepLoading()) {
    s->mesh.addLoad(s->loadIncrement);
    s->mesh.applyTransformationToSystemDeformation(loadStepTransform);

    // Modifies the nodeDisplacements
    s->setInitialGuess(loadStepTransform);

    // Minimizes the energy by moving the free nodes in the mesh
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

  SimPtr s = getPeriodicBorderSimulation(config, dataPath, loadedSimulation);

  std::string dumps[] = {"dump1", "dump2"};
  int dumpInUse = 0;
  int saveInterval = 10;
  double lastEnergy = 0;
  int step = 0;
  while (s->keepLoading()) {
    s->mesh.addLoad(s->loadIncrement);
    s->mesh.applyTransformationToSystemDeformation(loadStepTransform);

    // Modifies the nodeDisplacements
    s->setInitialGuess(loadStepTransform);

    // Minimizes the energy by moving the free nodes in the mesh
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

void doubleDislocationTest(Config config, std::string dataPath,
                           SimPtr loadedSimulation) {
  // SimPtr s = getFixedBorderSimulation(config, dataPath, loadedSimulation);

  SimPtr s = initOrLoad(config, dataPath, loadedSimulation, [](SimPtr s) {
    s->mesh.fixNodesInRow(0);
    s->mesh.fixNodesInColumn(0);
  });

  while (s->mesh.load < 1) {
    s->mesh.addLoad(s->loadIncrement);
    s->mesh.moveMeshSection(0.0, s->mesh.a * config.rows / 2.0 - 0.5,
                            Vector2d{config.loadIncrement, 0});

    // Minimizes the energy by moving the free nodes in the mesh
    s->minimize();

    // Updates progress and writes to file
    s->finishStep();
  }

  s->mesh.reconnect();

  while (s->mesh.load < config.maxLoad) {
    s->mesh.addLoad(s->loadIncrement);
    s->mesh.moveMeshSection(s->mesh.a * config.cols / 2.0 - 0.5, 0.0,
                            Vector2d{0, config.loadIncrement});

    // Minimizes the energy by moving the free nodes in the mesh
    s->minimize();

    // Updates progress and writes to file
    s->finishStep();
  }

  s->finishSimulation();
}

void singleDislocationFixedBoundaryTest(Config config, std::string dataPath,
                                        SimPtr loadedSimulation) {
  // SimPtr s = getFixedBorderSimulation(config, dataPath, loadedSimulation);

  SimPtr s = initOrLoad(config, dataPath, loadedSimulation, [](SimPtr s) {
    s->mesh.fixNodesInRow(0);
    s->mesh.fixNodesInColumn(0);
    s->mesh.fixNodesInColumn(-1);
  });

  while (s->mesh.load < 1) {
    s->mesh.addLoad(s->loadIncrement);
    s->mesh.moveMeshSection(0.0, s->mesh.a * config.rows / 2.0 - 0.5,
                            Vector2d{config.loadIncrement, 0}, true, false, 2);

    // Minimizes the energy by moving the free nodes in the mesh
    s->minimize(false);

    // Updates progress and writes to file
    s->finishStep();
  }

  // Reconnect in the middle row
  for (int i = 0; i < s->mesh.elements.size(); i += 2) {

    TElement &e = s->mesh.elements[i];
    // Check if the element com is in the middle row
    int middleRow = std::round(s->mesh.rows / 2);
    if (floor(e.getCom()[1]) == middleRow - 1) {

      int twinIndex = e.getElementTwin(s->mesh);
      // If we found a twin
      if (twinIndex != -1) {
        std::cout << "Element " << i << " has twin " << twinIndex << '\n';
        TElement &twin = s->mesh.elements[twinIndex];

        s->mesh.fixElementPair(e, twin);
        writeMeshToVtu(s->mesh, s->mesh.simName, dataPath,
                       std::to_string(twin.eIndex));
      } else {
        std::cout << "No twin found for element " << i << '\n';
      }
    }
    s->finishSimulation();
  }
}

// Alternative implementation: step along a set of points with fixed step size
// Returns {current displacement, finished flag}
std::pair<Vector2d, bool>
getTargetDisplacement(const std::vector<Vector2d> &points) {
  static const double stepSize = 0.002; // Fixed step size
  static Vector2d currentDisplacement(0, 0);
  static size_t currentTargetIdx = 0;

  // All points visited → finished
  if (currentTargetIdx >= points.size()) {
    return {currentDisplacement, true};
  }

  Vector2d target = points[currentTargetIdx];
  Vector2d toTarget = target - currentDisplacement;
  double distance = toTarget.norm();

  if (distance < stepSize) {
    // Snap exactly onto the target and advance
    currentDisplacement = target;
    ++currentTargetIdx;
    bool finished = (currentTargetIdx >= points.size());
    return {currentDisplacement, finished}; // Spend one frame at target
  } else {
    // Move a single fixed‑length step toward the target
    currentDisplacement += toTarget.normalized() * stepSize;
    return {currentDisplacement, false};
  }
}

void reconnectTest(Config config, std::string dataPath,
                   SimPtr loadedSimulation) {
  SimPtr s = getFixedBorderSimulation(config, dataPath, loadedSimulation);

  // We assume 3x3 mesh
  assert(s->mesh.cols == 3);
  assert(s->mesh.rows == 3);

  static const std::vector<Vector2d> points = {
      Vector2d(0, 0),      Vector2d(-0.2, 0.2),  Vector2d(0.2, -0.2),
      Vector2d(0, 0),      Vector2d(0.3, 0),     Vector2d(0.3, 0.3),
      Vector2d(-0.3, 0.3), Vector2d(-0.3, -0.3), Vector2d(0.3, -0.3),
      Vector2d(0.3, 0),    Vector2d(0, 0)};

  bool finished = false;
  while (!finished) {
    // Advance geometry along the preset path
    auto res = getTargetDisplacement(points);
    Vector2d targetDisp = res.first;
    finished = res.second;

    // Update loading (optional – keeps load progressing visually)
    s->mesh.addLoad(s->loadIncrement);

    // Apply the displacement to the middle node
    s->mesh.nodes(1, 1).setDisplacement(targetDisp);
    s->mesh.updateMesh();

    // Bookkeeping and output
    s->finishStep(false);
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
  // simulation If the simulation is complete, the program will terminate
  // and not rerun the simulation unless forceReRun is true
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
    auto sPtr = std::make_shared<Simulation>(); // Create a shared pointer to
                                                // a new Simulation object

    if (!configPath.empty()) { // Check if configPath is provided
      std::cout << "Overwriting simulation settings using\n - " << configPath
                << '\n';
    } else {
      // Try to find configPath in the same folder as dumpPath
      configPath = searchForConfig(dumpPath);
      if (configPath.empty()) { // If no configPath is found
        // Technically, we don't need a config file, since we can use the
        // settings from the dump, but to force the user to always keep a
        // config file in the folder, we throw an error here by design.
        std::cerr << "Error! No config provided or found in the same folder as "
                     "the dump file.\n";
        return; // Exit the function
      }
    }
    // Load and resume the simulation using the provided or found configPath
    Simulation::loadSimulation(*sPtr, dumpPath, configPath, outputPath,
                               forceReRun);

    std::cout << "Resuming simulation using " << dumpPath << '\n'
              << " - Config File: " << sPtr->config.name << '\n'
              << " - Data Path: " << sPtr->dataPath << '\n'
              << sPtr->config << '\n'
              << " - Current Load: " << sPtr->mesh.load << std::endl;
    std::cout << std::endl;
    runSimulationScenario(sPtr->config, sPtr->dataPath,
                          sPtr); // Run the simulation scenario

    // If dumpPath is not provided but configPath is
  } else {
    if (configPath.empty()) {

      // Try to find configPath in the same folder as dumpPath
      configPath = searchForConfig(dumpPath);
      if (configPath.empty()) { // If no configPath is found
        // Technically, we don't need a config file, since we can use the
        // settings from the dump, but to force the user to always keep a
        // config file in the folder, we throw an error here by design.
        std::cerr << "Error! No config provided or found in the same folder as "
                     "the dump file.\n";
        return; // Exit the function
      }
    }
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

    runSimulationScenario(config,
                          outputPath); // Run the simulation scenario
  }
}

void runSimulationScenario(Config config, std::string dataPath,
                           SimPtr loadedSimulation) {
  static const std::unordered_map<
      std::string,
      std::function<void(const Config &, const std::string &, SimPtr)>>
      scenarioMap = {
          {"simpleShear", simpleShear},
          {"simpleShearFixedBoundary", simpleShearFixedBoundary},
          {"simpleShearWithNoise", simpleShearWithNoise},
          {"periodicBoundaryTest", periodicBoundaryTest},
          {"periodicBoundaryFixedComparisonTest",
           periodicBoundaryFixedComparisonTest},
          {"cyclicSimpleShear", cyclicSimpleShear},
          {"createDumpBeforeEnergyDrop", createDumpBeforeEnergyDrop},
          {"doubleDislocationTest", doubleDislocationTest},
          {"singleDislocationFixedBoundaryTest",
           singleDislocationFixedBoundaryTest},
          {"reconnectTest", reconnectTest},
      };

  auto it = scenarioMap.find(config.scenario);
  if (it == scenarioMap.end()) {
    throw std::invalid_argument("No matching scenario: " + config.scenario);
  }

  it->second(config, dataPath, loadedSimulation);
}

// This function is quite complicated.
// We want to be able to prepare a simulation, and then initialize it, but
// if we have a loadedSimulation from a dump, then we want to skip both the
// preparation and the initialization. When you use this function, you want
// to prepare the simulation using a lambda function inside the function
// call.
SimPtr initSimulation(Config config, std::string dataPath,
                      std::function<void(SimPtr)> prepFunction) {
  // Construct shared simulation pointer
  SimPtr s = std::make_shared<Simulation>(config, dataPath, true);

  // This is where we would fix the border nodes in fixed boundary
  // conditions and/or apply the initial load transformation The Reason this
  // is a bit convoluted is because the prep function needs to occur before
  // the initialization function.
  prepFunction(s);

  // Now we initialize which involves creating the elements in the mesh and
  // initialzing the minimization solvers.
  s->initialize();

  // This first step function is also special, as it attemps to give most
  // simulations of the same system size and seed the exact same starting
  // conditions.
  s->firstStep();
  return s;
}

SimPtr initOrLoad(Config config, std::string dataPath, SimPtr loadedSimulation,
                  std::function<void(SimPtr)> prepFunction) {
  // If loadedSimulation is not a nullptr, then we already have a simulation
  // to use, otherwise, we need to run initSimulation.
  return loadedSimulation ? loadedSimulation
                          : initSimulation(config, dataPath, prepFunction);
}

SimPtr getFixedBorderSimulation(Config config, std::string dataPath,
                                SimPtr loadedSimulation) {
  if (config.usingPBC) {
    throw std::logic_error(
        "Should not fix boarder nodes if we use PBC. Check config file.");
  }
  return initOrLoad(config, dataPath, loadedSimulation, [](SimPtr s) {
    s->mesh.fixBorderNodes();
    auto startLoadTransform = getShear(s->startLoad);
    s->mesh.applyTransformation(startLoadTransform);
  });
}

SimPtr getPeriodicBorderSimulation(Config config, std::string dataPath,
                                   SimPtr loadedSimulation) {
  return initOrLoad(config, dataPath, loadedSimulation, [](SimPtr s) {
    // This fixed node avoids translation (and maybe rotation?)
    // s->mesh.fixBottomLeftCorner();
    auto startLoadTransform = getShear(s->startLoad);
    // Assuming some periodic boundary-specific operations
    s->mesh.applyTransformation(startLoadTransform);
  });
}