#include "experiments.h"
#include "Data/data_export.h"
#include "Eigen/Core"
#include "Eigen/LU"
#include "Simulation/simulation.h"
#include <cassert>
#include <iostream>
#include <memory>
#include <stdexcept>

using SimPtr = std::shared_ptr<Simulation>;

namespace {
thread_local double activeMakeDumpAt = -1.0;

void applyRuntimeOptions(SimPtr s) {
  if (s != nullptr) {
    s->makeDumpAt = activeMakeDumpAt;
  }
}
} // namespace

void simpleShear(Config config, std::string dataPath, SimPtr loadedSimulation) {
  Matrix2d loadStepTransform = getShear(config.loadIncrement);

  SimPtr s = getPeriodicBorderSimulation(config, dataPath, loadedSimulation);

  while (s->keepLoading()) {
    // Modifies the nodeDisplacements
    s->applyAffineStep(loadStepTransform);

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
    // Modifies the nodeDisplacements
    s->applyAffineStep(loadStepTransform);

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
    // Modifies the nodeDisplacements
    s->applyAffineStep(loadStepTransform);

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
    s->applyLoadStepToGuess(loadStepTransform);

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

  SimPtr s = initOrLoad(config, dataPath, loadedSimulation, [](SimPtr s) {
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
    s->applyLoadStepToGuess(loadStepTransform);

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
  SimPtr s = initOrLoad(config, dataPath, loadedSimulation, [](SimPtr s) {
    s->mesh.fixBorderNodes();
    int fixedMiddleRow = std::floor(s->rows / 2);
    s->mesh.fixNodesInRow(fixedMiddleRow);
  });

  while (s->keepLoading()) {
    // Modifies the nodeDisplacements
    s->applyAffineStep(loadStepTransform);

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
    // Modifies the nodeDisplacements
    s->applyAffineStep(loadStepTransform);

    // Minimizes the energy by moving the free nodes in the mesh
    s->minimize();

    // Updates progress and writes to file
    s->finishStep();

    if (step == 0) {
      lastEnergy = s->mesh.totalEnergy;
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
    if (lastEnergy - s->mesh.totalEnergy < -0.0004 * s->mesh.nrElements) {
      s->saveSimulation(dumps[dumpInUse] + "_EnergyFall");
      break;
    }
    lastEnergy = s->mesh.totalEnergy;
  }
  s->finishSimulation();
}

enum class DoubleDislocationLoadDirection { Horizontal, Vertical };

static bool useLastDoubleDislocationBoundary(const Config &config) {
  return config.GP1 > 0.5;
}

static bool useVerticalFirstDoubleDislocationLoading(const Config &config) {
  return config.GP2 > 0.5;
}

static void
applyDoubleDislocationLoadStep(Simulation &simulation,
                               DoubleDislocationLoadDirection direction,
                               bool useLastBoundary) {
  Mesh &mesh = simulation.mesh;
  const double step = simulation.loadIncrement;
  const double sign = useLastBoundary ? -1.0 : 1.0;
  const double midpointX = mesh.a * mesh.cols / 2.0 - 0.5;
  const double midpointY = mesh.a * mesh.rows / 2.0 - 0.5;
  const double maxX = mesh.a * mesh.cols;

  if (direction == DoubleDislocationLoadDirection::Vertical) {
    const Vector2d displacement{0.0, sign * step};
    if (useLastBoundary) {
      mesh.moveMeshSection(0.0, 0.0, displacement, true, false, midpointX);
    } else {
      mesh.moveMeshSection(midpointX, 0.0, displacement);
    }
    return;
  }

  const Vector2d displacement{sign * step, 0.0};
  if (useLastBoundary) {
    mesh.moveMeshSection(0.0, 0.0, displacement, true, false, maxX, midpointY);
  } else {
    mesh.moveMeshSection(0.0, midpointY, displacement);
  }
}

static void
runDoubleDislocationLoadingPhase(Simulation &simulation, double targetLoad,
                                 DoubleDislocationLoadDirection direction,
                                 bool useLastBoundary) {
  while (simulation.mesh.load < targetLoad) {
    simulation.mesh.addLoad(simulation.loadIncrement);
    applyDoubleDislocationLoadStep(simulation, direction, useLastBoundary);

    // Minimizes the energy by moving the free nodes in the mesh.
    simulation.minimize();

    // Updates progress and writes to file.
    simulation.finishStep();
  }
}

void doubleDislocationTest(Config config, std::string dataPath,
                           SimPtr loadedSimulation) {
  /*
  doubleDislocationTest general parameters:
    GP1 <= 0.5 fixes the first row and first column, and pushes inward with
      positive vertical/horizontal displacements.
    GP1 > 0.5 fixes the last row and last column, and pushes inward with
      negative vertical/horizontal displacements.
    GP2 <= 0.5 loads horizontally first, then vertically.
    GP2 > 0.5 loads vertically first, then horizontally.
  */
  const bool useLastBoundary = useLastDoubleDislocationBoundary(config);
  const bool verticalFirst = useVerticalFirstDoubleDislocationLoading(config);

  SimPtr s = initOrLoad(config, dataPath, loadedSimulation,
                        [useLastBoundary](SimPtr s) {
                          s->mesh.compareEdgeFlipOptions = true;
                          if (useLastBoundary) {
                            s->mesh.fixNodesInRow(-1);
                            s->mesh.fixNodesInColumn(-1);
                          } else {
                            s->mesh.fixNodesInRow(0);
                            s->mesh.fixNodesInColumn(0);
                          }
                        });
  s->mesh.compareEdgeFlipOptions = true;

  const DoubleDislocationLoadDirection firstDirection =
      verticalFirst ? DoubleDislocationLoadDirection::Vertical
                    : DoubleDislocationLoadDirection::Horizontal;
  const DoubleDislocationLoadDirection secondDirection =
      verticalFirst ? DoubleDislocationLoadDirection::Horizontal
                    : DoubleDislocationLoadDirection::Vertical;

  runDoubleDislocationLoadingPhase(*s, 1.0, firstDirection, useLastBoundary);
  runDoubleDislocationLoadingPhase(*s, config.maxLoad, secondDirection,
                                   useLastBoundary);

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

        s->mesh.flipEdge(e, twin);
        s->mesh.ensureFull();
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
    s->mesh.markDirty();
    s->mesh.updateMesh();

    // Bookkeeping and output
    s->finishStep(false);
  }
  s->finishSimulation();
}

void reconnectSSTest(Config config, std::string dataPath,
                     SimPtr loadedSimulation) {
  SimPtr s = getPeriodicBorderSimulation(config, dataPath, loadedSimulation);
  Matrix2d refTransform = getShear(config.GP1);
  // We assume 3x3 mesh
  assert(s->mesh.cols == 3);
  assert(s->mesh.rows == 3);

  // We fix all nodes
  s->mesh.fixNodesInColumn(0);
  s->mesh.fixNodesInColumn(1);
  s->mesh.fixNodesInColumn(2);
  // Set reference config
  s->mesh.applyTransformation(refTransform);
  s->mesh.setRefConfiguration();
  s->mesh.applyTransformation(refTransform.inverse());
  // We put the center node out of equilibrium
  s->mesh.nodes(1, 1).addDisplacement({0, config.GP2});

  Matrix2d loadStepTransform = getShear(config.loadIncrement);
  while (s->keepLoading()) {
    // Modifies the nodeDisplacements
    // s->applyAffineStep(loadStepTransform);
    s->mesh.addLoad(config.loadIncrement);
    s->mesh.applyTransformation(loadStepTransform);
    if (config.GP3 > 0.5) {
      s->mesh.applyTransformation(loadStepTransform.transpose());
    }
    if (s->config.reconnectionMethod == "edgeFlip") {
      s->mesh.reconnect();
    }
    // Updates progress and writes to file
    s->finishStep(true);
  }
  s->finishSimulation();
}

void reversibilityProtocolTest(Config config, std::string dataPath,
                               SimPtr loadedSimulation) {
  SimPtr s = initOrLoad(config, dataPath, loadedSimulation,
                        [](SimPtr s) { s->addReversibilityCsvColumns(); });

  s->addReversibilityCsvColumns();

  const double eps = 1e-4;

  while (s->keepLoading()) {
    Matrix2d loadStepTransform = getShear(s->loadIncrement);
    // This function steps the simulation and fills the reversibilityState
    // struct with the results.
    s->checkReversibility(loadStepTransform, eps);

    // Updates progress and writes to file
    s->finishStep();
  }
  s->finishSimulation();
}

void simpleShearReferenceTest(Config config, std::string dataPath,
                              SimPtr loadedSimulation) {
  // This is an experiment where we test how the influence of the reference
  // position affects the simulation. (We can then decouple the influence of the
  // reference configuration and the geometry of the elements) GP1 serves as the
  // shear used to set the reference configuration state.
  Matrix2d loadStepTransform = getShear(config.loadIncrement);
  Matrix2d refTransform = getShear(config.GP1);

  SimPtr s = std::make_shared<Simulation>(config, dataPath, true);
  applyRuntimeOptions(s);
  s->initialize();
  // Set reference config
  s->mesh.applyTransformation(refTransform);
  s->mesh.setRefConfiguration();
  s->mesh.applyTransformation(refTransform.inverse());

  while (s->keepLoading()) {
    // Modifies the nodeDisplacements
    s->applyAffineStep(loadStepTransform);

    // Minimizes the energy by moving the free nodes in the mesh
    s->minimize();

    // Updates progress and writes to file
    s->finishStep(true);
  }
  s->finishSimulation();
}

void runSimulationExperiment(Config config, std::string dataPath,
                             SimPtr loadedSimulation, double makeDumpAt) {
  static const std::unordered_map<
      std::string,
      std::function<void(const Config &, const std::string &, SimPtr)>>
      experimentMap = {
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
          {"reconnectSSTest", reconnectSSTest},
          {"reversibilityProtocolTest", reversibilityProtocolTest},
          {"simpleShearReferenceTest", simpleShearReferenceTest},
      };

  auto it = experimentMap.find(config.experiment);
  if (it == experimentMap.end()) {
    throw std::invalid_argument("No matching experiment: " + config.experiment);
  }

  const double previousMakeDumpAt = activeMakeDumpAt;
  activeMakeDumpAt = makeDumpAt;
  try {
    it->second(config, dataPath, loadedSimulation);
  } catch (...) {
    activeMakeDumpAt = previousMakeDumpAt;
    throw;
  }
  activeMakeDumpAt = previousMakeDumpAt;
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
  applyRuntimeOptions(s);

  // This is where we would fix the border nodes in fixed boundary
  // conditions and/or apply the initial load transformation The Reason this
  // is a bit convoluted is because the prep function needs to occur before
  // the initialization function.
  prepFunction(s);

  // This first step function is special, as it attempts to give most
  // simulations of the same system size and seed the exact same starting
  // conditions. It also handles initialization when needed.
  s->firstStep();
  return s;
}

SimPtr initOrLoad(Config config, std::string dataPath, SimPtr loadedSimulation,
                  std::function<void(SimPtr)> prepFunction) {
  // If loadedSimulation is not a nullptr, then we already have a simulation
  // to use, otherwise, we need to run initSimulation.
  if (loadedSimulation) {
    applyRuntimeOptions(loadedSimulation);
    return loadedSimulation;
  }
  return initSimulation(config, dataPath, prepFunction);
}

SimPtr getFixedBorderSimulation(Config config, std::string dataPath,
                                SimPtr loadedSimulation) {
  if (config.usingPBC) {
    throw std::logic_error(
        "Should not fix boarder nodes if we use PBC. Check config file.");
  }
  return initOrLoad(config, dataPath, loadedSimulation,
                    [](SimPtr s) { s->mesh.fixBorderNodes(); });
}

SimPtr getPeriodicBorderSimulation(Config config, std::string dataPath,
                                   SimPtr loadedSimulation) {
  return initOrLoad(config, dataPath, loadedSimulation, [](SimPtr s) {
    // This fixed node avoids translation (and maybe rotation?)
    // s->mesh.fixBottomLeftCorner();
  });
}
