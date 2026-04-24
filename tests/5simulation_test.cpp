#include "../src/Simulation/scenarios.h"
#include "../src/Simulation/simulation.h" // Include the header for your surface struct
#include "Data/data_export.h"
#include "Mesh/mesh.h"
#include "Mesh/tElement.h"
#include "Simulation/randomUtils.h"
#include "run/doctest.h"
#include "settings.h"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>

// To speed up testing, we use smaller simulations. This makes the tests less
// valid, so set it to true for a more thorough test.
#define FULLTESTS false

TEST_CASE("Simulation Save/Load mesh Test") {
  // Create a simple config
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 2;
  testConfig.cols = 2;
  testConfig.loadIncrement = 0.1;
  testConfig.maxLoad = 0.2;
  testConfig.usingPBC = true;
  testConfig.name = "2x2PBCSaveLoadTest";

  // Create and initialize simulation
  std::string dataPath = "test_data";
  Simulation sim(testConfig, dataPath, true);
  sim.initialize();

  Matrix2d loadStepTransform = getShear(testConfig.loadIncrement);
  sim.applyLoadStepToGuess(loadStepTransform);

  // Run a simulation step
  sim.mesh.addLoad(sim.loadIncrement);
  sim.mesh.applyTransformationToSystemDeformation(loadStepTransform);
  sim.minimize();

  sim.mesh.updateMesh();
  sim.mesh.updateAngles();
  // Save simulation to file
  std::string saveFileName = "test_sim_save";
  std::string pathToDump = sim.saveSimulation(saveFileName);

  // Load simulation into a new object
  Simulation loadedSim;
  Simulation::loadSimulation(loadedSim, pathToDump, "", dataPath, true);

  CHECK(loadedSim.mesh == sim.mesh);
  if (loadedSim.mesh != sim.mesh) {
    std::cout << debugCompare(loadedSim.mesh, sim.mesh) << std::endl;
  }
}

TEST_CASE("Mesh Snapshot Save/Restore Test") {
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 3;
  testConfig.cols = 4;
  testConfig.usingPBC = false;
  testConfig.name = "SnapshotSaveRestoreTest";

  std::string dataPath = "test_data";
  Simulation sim(testConfig, dataPath, true);
  sim.initialize();

  sim.mesh.fixBorderNodes();
  size_t fixedCount = 0;
  for (int row = 0; row < sim.mesh.rows; row++) {
    for (int col = 0; col < sim.mesh.cols; col++) {
      if (sim.mesh.nodes(row, col).fixedNode) {
        fixedCount++;
      }
    }
  }
  const int rows = sim.mesh.rows;
  const int cols = sim.mesh.cols;
  const size_t expectedFixed =
      (rows <= 1 || cols <= 1) ? static_cast<size_t>(rows * cols)
                               : static_cast<size_t>(2 * rows + 2 * cols - 4);
  CHECK(fixedCount == expectedFixed);

  const double loadValue = 1.23;
  const int loadSteps = 7;
  const Matrix2d deformation = getShear(0.15);
  sim.mesh.load = loadValue;
  sim.mesh.loadSteps = loadSteps;
  sim.mesh.currentDeformation = deformation;

  auto expectedDisp = [](int row, int col) {
    return Vector2d{0.1 * (row + 1), -0.05 * (col + 1)};
  };
  for (int row = 0; row < sim.mesh.rows; row++) {
    for (int col = 0; col < sim.mesh.cols; col++) {
      sim.mesh.nodes(row, col).setDisplacement(expectedDisp(row, col));
    }
  }

  Mesh::Snapshot snap;
  sim.mesh.saveSnapshot(snap);
  const size_t nodeCount = static_cast<size_t>(sim.mesh.rows * sim.mesh.cols);
  CHECK(snap.displacements.size() == 2 * nodeCount);

  sim.mesh.load = 0.0;
  sim.mesh.loadSteps = 0;
  sim.mesh.currentDeformation = Matrix2d::Identity();
  for (int row = 0; row < sim.mesh.rows; row++) {
    for (int col = 0; col < sim.mesh.cols; col++) {
      sim.mesh.nodes(row, col).setDisplacement(Vector2d::Zero());
    }
  }

  sim.mesh.restoreState(snap);

  CHECK(sim.mesh.load == doctest::Approx(loadValue));
  CHECK(sim.mesh.loadSteps == loadSteps);
  for (int r = 0; r < 2; r++) {
    for (int c = 0; c < 2; c++) {
      CHECK(sim.mesh.currentDeformation(r, c) ==
            doctest::Approx(deformation(r, c)));
    }
  }

  for (int row = 0; row < sim.mesh.rows; row++) {
    for (int col = 0; col < sim.mesh.cols; col++) {
      const Vector2d expected = expectedDisp(row, col);
      const Vector2d actual = sim.mesh.nodes(row, col).u();
      CHECK(actual[0] == doctest::Approx(expected[0]));
      CHECK(actual[1] == doctest::Approx(expected[1]));
    }
  }
}

TEST_CASE("Simulation Save/Load Triangular Mesh") {
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 3;
  testConfig.cols = 3;
  testConfig.loadIncrement = 0.1;
  testConfig.maxLoad = 0.2;
  testConfig.usingPBC = true;
  testConfig.energyFunction = "contiTriangular";
  testConfig.name = "TriangularSaveLoadTest";

  std::string dataPath = "test_data";
  Simulation sim(testConfig, dataPath, true);
  sim.initialize();
  sim.mesh.updateMesh();
  sim.mesh.updateAveragesAndPlasticEvents();
  sim.mesh.updateAngles();

  std::string saveFileName = "triangular_sim_save";
  std::string pathToDump = sim.saveSimulation(saveFileName);

  Simulation loadedSim;
  Simulation::loadSimulation(loadedSim, pathToDump, "", dataPath, true);

  CHECK(loadedSim.mesh == sim.mesh);
  if (loadedSim.mesh != sim.mesh) {
    std::cout << debugCompare(loadedSim.mesh, sim.mesh) << std::endl;
  }

  const double h = std::sqrt(3.0) * 0.5 * loadedSim.mesh.a;
  CHECK(loadedSim.mesh.nodes(1, 0).pos()[0] ==
        doctest::Approx(0.5 * loadedSim.mesh.a));
  CHECK(loadedSim.mesh.nodes(1, 0).pos()[1] == doctest::Approx(h));
}

TEST_CASE("Simulation Save/Load Exact Node Positions") {
  Config testConfig;
  testConfig.setDefaultValues();
  const int meshSize = FULLTESTS ? 30 : 5;
  testConfig.rows = meshSize;
  testConfig.cols = meshSize;
  testConfig.usingPBC = true;
  testConfig.name = "SaveLoadExactPositions";
  testConfig.forceReRun = true;

  std::string dataPath = "test_data";
  Simulation sim(testConfig, dataPath, true);
  sim.initialize();

  // Add deterministic noise directly to node displacements.
  setSeed(1234);
  for (int row = 0; row < sim.mesh.rows; ++row) {
    for (int col = 0; col < sim.mesh.cols; ++col) {
      Node &n = sim.mesh.nodes(row, col);
      const Vector2d disp(sampleNormal(0.0, 0.01), sampleNormal(0.0, 0.01));
      n.addDisplacement(disp);
    }
  }
  sim.mesh.markDirty();
  sim.mesh.updateMesh();
  // Apply a load to make the minimization less trivial.
  const Matrix2d loadTransform = getShear(0.4);
  sim.mesh.addLoad(0.4);
  sim.mesh.applyTransformationToSystemDeformation(loadTransform);
  sim.applyLoadStepToGuess(loadTransform);

  std::string dumpPath = sim.saveSimulation("ExactPositions");

  Simulation loadedSim;
  Simulation::loadSimulation(loadedSim, dumpPath, "", dataPath, true);

  REQUIRE(loadedSim.mesh.rows == sim.mesh.rows);
  REQUIRE(loadedSim.mesh.cols == sim.mesh.cols);

  bool nodesDiffer = false;
  for (int row = 0; row < sim.mesh.rows; ++row) {
    for (int col = 0; col < sim.mesh.cols; ++col) {
      const Node &lhs = sim.mesh.nodes(row, col);
      const Node &rhs = loadedSim.mesh.nodes(row, col);
      if (lhs.pos()[0] != rhs.pos()[0] || lhs.pos()[1] != rhs.pos()[1]) {
        nodesDiffer = true;
        break;
      }
    }
    if (nodesDiffer) {
      break;
    }
  }
  CHECK(nodesDiffer == false);
}

TEST_CASE("Simulation Save/Load Minimize Determinism") {
  Config testConfig;
  testConfig.setDefaultValues();
  const int meshSize = FULLTESTS ? 30 : 5;
  testConfig.rows = meshSize;
  testConfig.cols = meshSize;
  testConfig.usingPBC = true;
  testConfig.reconnectionMethod = "none";
  testConfig.name = "SaveLoadMinimize";
  testConfig.forceReRun = true;

  std::string dataPath = "test_data";
  clearOutputFolder(testConfig.name, dataPath);

  Simulation sim(testConfig, dataPath, true);
  sim.initialize();

  // Add deterministic noise directly to node displacements.
  setSeed(1234);
  for (int row = 0; row < sim.mesh.rows; ++row) {
    for (int col = 0; col < sim.mesh.cols; ++col) {
      Node &n = sim.mesh.nodes(row, col);
      const Vector2d disp(sampleNormal(0.0, 0.01), sampleNormal(0.0, 0.01));
      n.addDisplacement(disp);
    }
  }
  sim.mesh.markDirty();
  sim.mesh.updateMesh();
  // Apply a load to make the minimization less trivial.
  const Matrix2d loadTransform = getShear(0.4);
  sim.mesh.addLoad(0.4);
  sim.mesh.applyTransformationToSystemDeformation(loadTransform);
  sim.applyLoadStepToGuess(loadTransform);

  std::string dumpPath = sim.saveSimulation("MinimizeDeterminism");

  Simulation loadedSim;
  Simulation::loadSimulation(loadedSim, dumpPath, "", dataPath, true);
  loadedSim.applyLoadStepToGuess(Matrix2d::Identity());

  // Minimize both simulations from the same starting state.
  sim.minimize(false);
  loadedSim.minimize(false);

  REQUIRE(loadedSim.mesh.rows == sim.mesh.rows);
  REQUIRE(loadedSim.mesh.cols == sim.mesh.cols);

  bool nodesDiffer = false;
  bool forcesDiffer = false;
  for (int row = 0; row < sim.mesh.rows; ++row) {
    for (int col = 0; col < sim.mesh.cols; ++col) {
      const Node &lhs = sim.mesh.nodes(row, col);
      const Node &rhs = loadedSim.mesh.nodes(row, col);
      if (lhs.pos()[0] != rhs.pos()[0] || lhs.pos()[1] != rhs.pos()[1]) {
        nodesDiffer = true;
        break;
      }
      if (lhs.f[0] != rhs.f[0] || lhs.f[1] != rhs.f[1]) {
        forcesDiffer = true;
        break;
      }
    }
    if (nodesDiffer) {
      break;
    }
    if (forcesDiffer) {
      break;
    }
  }
  CHECK(nodesDiffer == false);
  CHECK(forcesDiffer == false);
}

TEST_CASE("Simulation Save/Load Minimize Determinism With Reconnect") {
  Config testConfig;
  testConfig.setDefaultValues();
  const int meshSize = FULLTESTS ? 30 : 5;
  testConfig.rows = meshSize;
  testConfig.cols = meshSize;
  testConfig.usingPBC = true;
  testConfig.reconnectionMethod = "edgeFlip";
  testConfig.name = "SaveLoadMinimizeReconnect";
  testConfig.forceReRun = true;

  std::string dataPath = "test_data";
  clearOutputFolder(testConfig.name, dataPath);

  Simulation sim(testConfig, dataPath, true);
  sim.initialize();

  // Add deterministic noise directly to node displacements.
  setSeed(1234);
  for (int row = 0; row < sim.mesh.rows; ++row) {
    for (int col = 0; col < sim.mesh.cols; ++col) {
      Node &n = sim.mesh.nodes(row, col);
      const Vector2d disp(sampleNormal(0.0, 0.01), sampleNormal(0.0, 0.01));
      n.addDisplacement(disp);
    }
  }
  sim.mesh.markDirty();
  sim.mesh.updateMesh();
  // Apply a load to make the minimization less trivial.
  const Matrix2d loadTransform = getShear(0.4);
  sim.mesh.addLoad(0.4);
  sim.mesh.applyTransformationToSystemDeformation(loadTransform);
  sim.applyLoadStepToGuess(loadTransform);

  std::string dumpPath = sim.saveSimulation("MinimizeReconnect");

  Simulation loadedSim;
  Simulation::loadSimulation(loadedSim, dumpPath, "", dataPath, true);
  loadedSim.applyLoadStepToGuess(Matrix2d::Identity());

  // Minimize both simulations from the same starting state (with reconnection).
  sim.minimize(true);
  loadedSim.minimize(true);

  REQUIRE(loadedSim.mesh.rows == sim.mesh.rows);
  REQUIRE(loadedSim.mesh.cols == sim.mesh.cols);

  bool nodesDiffer = false;
  bool forcesDiffer = false;
  for (int row = 0; row < sim.mesh.rows; ++row) {
    for (int col = 0; col < sim.mesh.cols; ++col) {
      const Node &lhs = sim.mesh.nodes(row, col);
      const Node &rhs = loadedSim.mesh.nodes(row, col);
      if (lhs.pos()[0] != rhs.pos()[0] || lhs.pos()[1] != rhs.pos()[1]) {
        nodesDiffer = true;
        break;
      }
      if (lhs.f[0] != rhs.f[0] || lhs.f[1] != rhs.f[1]) {
        forcesDiffer = true;
        break;
      }
    }
    if (nodesDiffer || forcesDiffer) {
      break;
    }
  }
  CHECK(nodesDiffer == false);
  CHECK(forcesDiffer == false);
}

TEST_CASE("Simulation Save/Revert Minimize Determinism") {
  Config testConfig;
  testConfig.setDefaultValues();
  const int meshSize = FULLTESTS ? 30 : 5;
  testConfig.rows = meshSize;
  testConfig.cols = meshSize;
  testConfig.usingPBC = true;
  testConfig.name = "SaveRevertTest";
  testConfig.reconnectionMethod = "none";
  testConfig.epsR = 0.1;
  testConfig.forceReRun = true;

  std::string dataPath = "test_data";
  // Recompute derived values from the current state so we can compare meshes.
  auto normalizeForCompare = [](Mesh &m) {
    // These are history-style maxima, so after a revert they may reflect the
    // (now discarded) path. Reset them so they only reflect the current state.
    m.maxM3Nr = 0;
    m.maxPlasticJump = 0;
    m.updateAveragesAndPlasticEvents();
    m.updateBoundingBox();
    m.updateAngles();
  };

  Simulation sim(testConfig, dataPath, true);
  sim.initialize();

  const double load = 0.4;
  const Matrix2d loadTransform = getShear(load);

  sim.mesh.addLoad(load);
  sim.mesh.applyTransformationToSystemDeformation(loadTransform);

  setSeed(static_cast<unsigned int>(testConfig.seed));
  sim.applyLoadStepToGuess(loadTransform);
  sim.addNoiseToGuess(0.1);
  sim.mesh.updateMesh();
  sim.mesh.updateAveragesAndPlasticEvents();
  Mesh::Snapshot beforeMinSnapshot = sim.mesh.snapshotState();
  normalizeForCompare(sim.mesh);
  Mesh meshBeforeMin = sim.mesh;
  meshBeforeMin.writeToVtu("SavedMeshBeforeMinimize");

  sim.minimize();

  normalizeForCompare(sim.mesh);
  Mesh meshAfterMin = sim.mesh;
  meshAfterMin.writeToVtu("SavedMeshAfterMinimize");

  sim.mesh.restoreState(beforeMinSnapshot);
  normalizeForCompare(sim.mesh);
  Mesh meshRevertBeforeMin = sim.mesh;
  meshRevertBeforeMin.writeToVtu("SavedMeshAfterRevertBeforeMinimize");
  // We expect the mesh to be exactly the same as before minimize
  CHECK(meshRevertBeforeMin == meshBeforeMin);
  if (sim.mesh != meshBeforeMin) {
    std::cout << debugCompare(sim.mesh, meshBeforeMin) << std::endl;
  }
  CHECK(sim.mesh == meshBeforeMin);
  if (sim.mesh != meshBeforeMin) {
    std::cout << debugCompare(sim.mesh, meshBeforeMin) << std::endl;
  }

  setSeed(static_cast<unsigned int>(testConfig.seed));
  // After revert, the mesh already includes the load transform. We only want
  // to sync the guess arrays to the current mesh state (no extra deformation).
  sim.applyLoadStepToGuess(Matrix2d::Identity());
  sim.addNoiseToGuess(0.1);
  sim.minimize();
  normalizeForCompare(sim.mesh);
  sim.mesh.writeToVtu("SavedMeshAfterSecondMinimize");
  CHECK(sim.mesh == meshAfterMin);
  if (sim.mesh != meshAfterMin) {
    std::cout << debugCompare(sim.mesh, meshAfterMin) << std::endl;
  }
}

TEST_CASE("Simulation Save/Load mesh Test") {
  // Create a simple config
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 3;
  testConfig.cols = 3;
  testConfig.loadIncrement = 0.5;
  // testConfig.scenario = "simpleShearFixedBoundary";
  testConfig.maxLoad = 1;
  testConfig.name = "3x3LoadingTestSaveLoad";
  testConfig.forceReRun = true;

  // Create a data path and file paths
  std::string dataPath = "test_data/";
  std::string dumpPath = dataPath + testConfig.name + "/dumps/dump_l1.0.xml.gz";

  // Remove old data
  clearOutputFolder(testConfig.name, dataPath);
  std::shared_ptr<Simulation> s =
      std::make_shared<Simulation>(testConfig, dataPath);
  // s->mesh.fixBorderNodes();
  s->firstStep();

  // Run the scenario and check CSV
  runSimulationScenario(testConfig, dataPath, s);
  Mesh oldMeshAtOne = s->mesh;
  s->maxLoad = 2;
  runSimulationScenario(testConfig, dataPath, s);
  Mesh oldMeshAtTwo = s->mesh;

  // Load simulation into a new object
  using SimPtr = std::shared_ptr<Simulation>;
  SimPtr loadedSim = std::make_shared<Simulation>(testConfig, dataPath);
  Simulation::loadSimulation(*loadedSim, dumpPath, "", dataPath);

  CHECK(loadedSim->mesh == oldMeshAtOne);
  if (loadedSim->mesh != oldMeshAtOne) {
    std::cout << debugCompare(loadedSim->mesh, oldMeshAtOne) << std::endl;
  }

  // Rerun
  loadedSim->maxLoad = 2;
  runSimulationScenario(testConfig, dataPath, loadedSim);
  // std::cout << oldMeshAtTwo << std::endl;
  // std::cout << loadedSim->mesh << std::endl;
  CHECK(loadedSim->mesh == oldMeshAtTwo);
  if (loadedSim->mesh != oldMeshAtTwo) {
    std::cout << debugCompare(loadedSim->mesh, oldMeshAtTwo) << std::endl;
  }
}

TEST_CASE("Simulation Save/Load Energy Test") {
  // Create a simple config
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 10;
  testConfig.cols = 10;
  testConfig.loadIncrement = 0.1;
  testConfig.maxLoad = 0.2;
  testConfig.name = "10x10PBCSaveLoadTest";
  testConfig.usingPBC = true;
  // Create and initialize simulation
  std::string dataPath = "test_data";
  Simulation sim(testConfig, dataPath, true);
  sim.initialize();

  Matrix2d loadStepTransform = getShear(testConfig.loadIncrement);
  sim.applyLoadStepToGuess(loadStepTransform);

  // Run a simulation step
  sim.mesh.addLoad(sim.loadIncrement);
  sim.mesh.applyTransformationToSystemDeformation(loadStepTransform);
  sim.minimize();

  bool allOnes = true;
  // sim.minState.scale is alglib_impl::ae_vector *
  // Check that all values are 1.0
  for (int i = 0; i < sim.minState.scale->cnt; i++) {
    if (sim.minState.scale->ptr.p_double[i] != 1.0) {
      allOnes = false;
      break;
    }
  }
  CHECK(allOnes == true);

  sim.mesh.updateMesh();
  double originalEnergy = sim.mesh.totalEnergy;
  // Save simulation to file
  std::string saveFileName = "test_sim_save";
  std::string pathToDump = sim.saveSimulation(saveFileName);

  // Load simulation into a new object
  Simulation loadedSim;
  Simulation::loadSimulation(loadedSim, pathToDump, "", dataPath, true);

  double loadedEnergy = sim.mesh.totalEnergy;
  CHECK(doctest::Approx(loadedEnergy).epsilon(1e-12) == originalEnergy);
  // Update properties
  loadedSim.mesh.updateMesh();

  CHECK(loadedSim.mesh == sim.mesh);
  if (loadedSim.mesh != sim.mesh) {
    std::cout << debugCompare(loadedSim.mesh, sim.mesh) << std::endl;
  }
}

TEST_CASE("Simulation Dump Serialization Side Effects") {
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 3;
  testConfig.cols = 3;
  testConfig.usingPBC = true;
  testConfig.name = "DumpSerializationSideEffects";
  testConfig.forceReRun = true;

  std::string dataPath = "test_data";
  clearOutputFolder(testConfig.name, dataPath);

  Simulation sim(testConfig, dataPath, true);
  sim.initialize();

  // Seed energy history with non-default values so we can detect mutations.
  sim.energyHistory.prevLoadStepTotalEnergy = 1.1;
  sim.energyHistory.prevLoadStepAverageEnergy = 2.2;
  sim.energyHistory.loadStepTotalEnergyChange = 3.3;
  sim.energyHistory.loadStepAverageEnergyChange = 4.4;
  sim.energyHistory.initialGuessTotalEnergy = 5.5;
  sim.energyHistory.initialGuessAverageEnergy = 6.6;
  sim.energyHistory.totalEnergyChangeFromInitialGuess = 7.7;
  sim.energyHistory.averageEnergyChangeFromInitialGuess = 8.8;
  sim.energyHistory.initialGuessAverageSigma12 = 9.9;
  sim.energyHistory.averageSigma12ChangeFromInitialGuess = 10.1;
  sim.energyHistory.prevMinIterTotalEnergy = 11.2;
  sim.energyHistory.prevMinIterAverageEnergy = 12.3;
  sim.energyHistory.minIterTotalEnergyChange = 13.4;
  sim.energyHistory.minIterAverageEnergyChange = 14.5;

  const SimulationEnergyHistory historySnapshot = sim.energyHistory;

  // Seed internal mesh buffers so we can detect serialization side effects.
  Mesh::Snapshot snapshotBeforeSave = sim.mesh.snapshotState();

  sim.saveSimulation("DumpSerializationSideEffects");

  auto checkHistoryEqual = [](const SimulationEnergyHistory &lhs,
                              const SimulationEnergyHistory &rhs) {
    CHECK(lhs.prevLoadStepTotalEnergy == rhs.prevLoadStepTotalEnergy);
    CHECK(lhs.prevLoadStepAverageEnergy == rhs.prevLoadStepAverageEnergy);
    CHECK(lhs.loadStepTotalEnergyChange == rhs.loadStepTotalEnergyChange);
    CHECK(lhs.loadStepAverageEnergyChange == rhs.loadStepAverageEnergyChange);
    CHECK(lhs.initialGuessTotalEnergy == rhs.initialGuessTotalEnergy);
    CHECK(lhs.initialGuessAverageEnergy == rhs.initialGuessAverageEnergy);
    CHECK(lhs.totalEnergyChangeFromInitialGuess ==
          rhs.totalEnergyChangeFromInitialGuess);
    CHECK(lhs.averageEnergyChangeFromInitialGuess ==
          rhs.averageEnergyChangeFromInitialGuess);
    CHECK(lhs.initialGuessAverageSigma12 == rhs.initialGuessAverageSigma12);
    CHECK(lhs.averageSigma12ChangeFromInitialGuess ==
          rhs.averageSigma12ChangeFromInitialGuess);
    CHECK(lhs.prevMinIterTotalEnergy == rhs.prevMinIterTotalEnergy);
    CHECK(lhs.prevMinIterAverageEnergy == rhs.prevMinIterAverageEnergy);
    CHECK(lhs.minIterTotalEnergyChange == rhs.minIterTotalEnergyChange);
    CHECK(lhs.minIterAverageEnergyChange == rhs.minIterAverageEnergyChange);
  };

  checkHistoryEqual(sim.energyHistory, historySnapshot);
  CHECK(sim.mesh.rmsDistanceToSnapshot(snapshotBeforeSave) ==
        doctest::Approx(0.0));
}

// We create a helper function to read the first column of the CSV file
// and compare it to expected values. This function:
//  1) Opens the CSV file and ensures it is open
//  2) Reads the header line and discards it
//  3) Reads each subsequent line, parses the first column as an integer
//  4) Checks if the read integers match the expected values
static void checkMacroDataCsv(const std::string &csvPath,
                              const std::vector<int> &expectedValues) {
  // Open file
  std::ifstream file(csvPath);
  REQUIRE(file.is_open());

  // Discard header
  std::string line;
  REQUIRE(std::getline(file, line));

  // Read lines and parse first column
  std::vector<int> actualValues;
  while (std::getline(file, line)) {
    std::stringstream ss(line);
    int val;
    ss >> val;
    actualValues.push_back(val);
  }

  // Check size
  REQUIRE(actualValues.size() == expectedValues.size());

  // Check each value
  for (size_t i = 0; i < expectedValues.size(); i++) {
    CHECK(actualValues[i] == expectedValues[i]);
  }
}

static std::vector<std::string> splitCsvLine(const std::string &line) {
  std::vector<std::string> cells;
  std::stringstream ss(line);
  std::string cell;
  while (std::getline(ss, cell, ',')) {
    cells.push_back(cell);
  }
  return cells;
}

static int findHeaderIndex(const std::vector<std::string> &headers,
                           const std::string &name) {
  for (size_t i = 0; i < headers.size(); ++i) {
    if (headers[i] == name) {
      return static_cast<int>(i);
    }
  }
  return -1;
}

// Here, in the main test, we run the simulation in steps and check CSV results
TEST_CASE("Simulation Save/Load Macro Data Test") {
  // Create a simple config
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 3;
  testConfig.cols = 3;
  testConfig.loadIncrement = 0.1;
  testConfig.maxLoad = 0.3;
  testConfig.usingPBC = true;
  testConfig.name = "3x3PBCLoadingTestWithSaveLoad";
  testConfig.forceReRun = true;

  // Create a data path and file paths
  std::string dataPath = "test_data/";
  std::string dumpPath =
      dataPath + testConfig.name + "/dumps/dump_l0.30.xml.gz";
  std::string csvPath = dataPath + testConfig.name + "/macroData.csv";

  // Remove old data
  clearOutputFolder(testConfig.name, dataPath);
  std::shared_ptr<Simulation> s =
      std::make_shared<Simulation>(testConfig, dataPath);
  s->firstStep();

  // Run the scenario and check CSV
  runSimulationScenario(testConfig, dataPath, s);

  // Check that the first column is 1, 2, 3, 4
  checkMacroDataCsv(csvPath, {1, 2, 3, 4});

  // Load simulation into a new object
  using SimPtr = std::shared_ptr<Simulation>;
  SimPtr loadedSim = std::make_shared<Simulation>(testConfig, dataPath);
  Simulation::loadSimulation(*loadedSim, dumpPath, "", dataPath, true);

  CHECK(loadedSim->mesh == s->mesh);
  if (loadedSim->mesh != s->mesh) {
    std::cout << debugCompare(loadedSim->mesh, s->mesh) << std::endl;
  }

  // After loading from l0.2, check that the first column is 1, 2, 3, 4
  // (Corresponding to load 0.0, 0.1, 0.2 and 0.3)
  checkMacroDataCsv(csvPath, {1, 2, 3, 4});

  // Increase max load, run again, and check appended results
  loadedSim->maxLoad = 0.4;
  runSimulationScenario(testConfig, dataPath, loadedSim);

  // Now, the first column should be 1, 2, 3, 4, 5
  checkMacroDataCsv(csvPath, {1, 2, 3, 4, 5});
}

TEST_CASE("Simulation Macro CSV Sanity") {
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 3;
  testConfig.cols = 3;
  testConfig.loadIncrement = 0.1;
  testConfig.maxLoad = 0.2;
  testConfig.usingPBC = true;
  testConfig.name = "3x3PBCLoadingCsvSanity";
  testConfig.forceReRun = true;
  testConfig.logDuringMinimization = false;

  std::string dataPath = "test_data/";
  std::string csvPath = dataPath + testConfig.name + "/macroData.csv";

  clearOutputFolder(testConfig.name, dataPath);
  std::shared_ptr<Simulation> s =
      std::make_shared<Simulation>(testConfig, dataPath);
  s->firstStep();
  runSimulationScenario(testConfig, dataPath, s);

  std::ifstream file(csvPath);
  REQUIRE(file.is_open());

  std::string line;
  REQUIRE(std::getline(file, line));
  const auto headers = splitCsvLine(line);

  const int idxTotalEnergy = findHeaderIndex(headers, "total_energy");
  const int idxTotalEnergyChange =
      findHeaderIndex(headers, "total_energy_change");
  const int idxTotalInitEnergy = findHeaderIndex(headers, "total_init_energy");
  const int idxTotalChangeFromInit =
      findHeaderIndex(headers, "total_e_change_from_init");
  const int idxAvgEnergy = findHeaderIndex(headers, "avg_energy");
  const int idxAvgEnergyChange = findHeaderIndex(headers, "avg_energy_change");
  const int idxAvgInitEnergy = findHeaderIndex(headers, "avg_init_energy");
  const int idxAvgChangeFromInit =
      findHeaderIndex(headers, "avg_e_change_from_init");
  const int idxMinIterTotalChange =
      findHeaderIndex(headers, "min_iter_total_energy_change");
  const int idxMinIterAvgChange =
      findHeaderIndex(headers, "min_iter_avg_energy_change");
  const int idxNrIterations = findHeaderIndex(headers, "nr_iterations");
  const int idxNrFuncEvals = findHeaderIndex(headers, "nr_func_evals");
  const int idxRedQ1 = findHeaderIndex(headers, "nr_red_q1");
  const int idxRedQ2 = findHeaderIndex(headers, "nr_red_q2");
  const int idxRedQ3 = findHeaderIndex(headers, "nr_red_q3");
  const int idxRedQ4 = findHeaderIndex(headers, "nr_red_q4");
  const int idxRedQ1Fixed = findHeaderIndex(headers, "nr_red_q1_fixed");
  const int idxRedQ2Fixed = findHeaderIndex(headers, "nr_red_q2_fixed");
  const int idxRedQ3Fixed = findHeaderIndex(headers, "nr_red_q3_fixed");
  const int idxRedQ4Fixed = findHeaderIndex(headers, "nr_red_q4_fixed");

  REQUIRE(idxTotalEnergy >= 0);
  REQUIRE(idxTotalEnergyChange >= 0);
  REQUIRE(idxTotalInitEnergy >= 0);
  REQUIRE(idxTotalChangeFromInit >= 0);
  REQUIRE(idxAvgEnergy >= 0);
  REQUIRE(idxAvgEnergyChange >= 0);
  REQUIRE(idxAvgInitEnergy >= 0);
  REQUIRE(idxAvgChangeFromInit >= 0);
  REQUIRE(idxMinIterTotalChange >= 0);
  REQUIRE(idxMinIterAvgChange >= 0);
  REQUIRE(idxNrIterations >= 0);
  REQUIRE(idxNrFuncEvals >= 0);
  REQUIRE(idxRedQ1 >= 0);
  REQUIRE(idxRedQ2 >= 0);
  REQUIRE(idxRedQ3 >= 0);
  REQUIRE(idxRedQ4 >= 0);
  REQUIRE(idxRedQ1Fixed >= 0);
  REQUIRE(idxRedQ2Fixed >= 0);
  REQUIRE(idxRedQ3Fixed >= 0);
  REQUIRE(idxRedQ4Fixed >= 0);

  const int nrElements =
      testConfig.usingPBC ? 2 * testConfig.rows * testConfig.cols
                          : 2 * (testConfig.rows - 1) * (testConfig.cols - 1);

  std::vector<double> totalEnergyValues;
  std::vector<double> avgEnergyValues;
  std::vector<double> totalEnergyChangeValues;
  std::vector<double> avgEnergyChangeValues;

  while (std::getline(file, line)) {
    const auto cells = splitCsvLine(line);
    REQUIRE(cells.size() == headers.size());

    const double totalEnergy = std::stod(cells[idxTotalEnergy]);
    const double totalEnergyChange = std::stod(cells[idxTotalEnergyChange]);
    const double totalInitEnergy = std::stod(cells[idxTotalInitEnergy]);
    const double totalChangeFromInit = std::stod(cells[idxTotalChangeFromInit]);
    const double avgEnergy = std::stod(cells[idxAvgEnergy]);
    const double avgEnergyChange = std::stod(cells[idxAvgEnergyChange]);
    const double avgInitEnergy = std::stod(cells[idxAvgInitEnergy]);
    const double avgChangeFromInit = std::stod(cells[idxAvgChangeFromInit]);
    const double minIterTotalChange = std::stod(cells[idxMinIterTotalChange]);
    const double minIterAvgChange = std::stod(cells[idxMinIterAvgChange]);
    const double nrIterations = std::stod(cells[idxNrIterations]);
    const double nrFuncEvals = std::stod(cells[idxNrFuncEvals]);
    const double redQ1 = std::stod(cells[idxRedQ1]);
    const double redQ2 = std::stod(cells[idxRedQ2]);
    const double redQ3 = std::stod(cells[idxRedQ3]);
    const double redQ4 = std::stod(cells[idxRedQ4]);
    const double redQ1Fixed = std::stod(cells[idxRedQ1Fixed]);
    const double redQ2Fixed = std::stod(cells[idxRedQ2Fixed]);
    const double redQ3Fixed = std::stod(cells[idxRedQ3Fixed]);
    const double redQ4Fixed = std::stod(cells[idxRedQ4Fixed]);

    totalEnergyValues.push_back(totalEnergy);
    avgEnergyValues.push_back(avgEnergy);
    totalEnergyChangeValues.push_back(totalEnergyChange);
    avgEnergyChangeValues.push_back(avgEnergyChange);

    CHECK(totalEnergy >= 0.0);
    CHECK(avgEnergy >= 0.0);
    CHECK(std::abs(totalEnergy - avgEnergy * nrElements) < 1e-8);
    CHECK(std::abs((totalEnergy - totalInitEnergy) - totalChangeFromInit) <
          1e-8);
    CHECK(std::abs((avgEnergy - avgInitEnergy) - avgChangeFromInit) < 1e-8);

    // When we are not logging during minimization, min-iter deltas should be 0
    CHECK(std::abs(minIterTotalChange) < 1e-12);
    CHECK(std::abs(minIterAvgChange) < 1e-12);

    // If there were at least two function evaluations, we expect iterations.
    if (nrFuncEvals >= 2) {
      CHECK(nrIterations > 0);
    }

    const double sumRed = redQ1 + redQ2 + redQ3 + redQ4;
    const double sumRedFixed =
        redQ1Fixed + redQ2Fixed + redQ3Fixed + redQ4Fixed;
    CHECK(std::abs(sumRed - nrElements) < 1e-8);
    CHECK(std::abs(sumRedFixed - nrElements) < 1e-8);
  }

  REQUIRE(totalEnergyValues.size() >= 2);
  for (size_t i = 0; i < totalEnergyValues.size(); ++i) {
    const double totalEnergy = totalEnergyValues[i];
    const double avgEnergy = avgEnergyValues[i];
    const double totalEnergyChange = totalEnergyChangeValues[i];
    const double avgEnergyChange = avgEnergyChangeValues[i];

    if (i == 0) {
      CHECK(std::abs(totalEnergyChange) < 1e-12);
      CHECK(std::abs(avgEnergyChange) < 1e-12);
      continue;
    }

    // Energy change should not be zero after the first line.
    const double totalDelta = totalEnergy - totalEnergyValues[i - 1];
    const double avgDelta = avgEnergy - avgEnergyValues[i - 1];
    CHECK(std::abs(totalDelta) > 1e-12);
    CHECK(std::abs(avgDelta) > 1e-12);

    // CSV deltas should match the step-to-step change.
    CHECK(std::abs(totalDelta - totalEnergyChange) < 1e-8);
    CHECK(std::abs(avgDelta - avgEnergyChange) < 1e-8);
  }
}

TEST_CASE("CSV Header Change Comment on Resume") {
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 3;
  testConfig.cols = 3;
  testConfig.usingPBC = true;
  testConfig.loadIncrement = 0.1;
  testConfig.maxLoad = 0.5;
  testConfig.scenario = "simpleShear";
  testConfig.name = "CsvHeaderChangeCommentTest";
  testConfig.forceReRun = true;

  std::string dataPath = "test_data/";
  clearOutputFolder(testConfig.name, dataPath);

  auto sim = std::make_shared<Simulation>(testConfig, dataPath, true);
  sim->firstStep();

  const std::string csvPath =
      getOutputPath(testConfig.name, dataPath) + MACRODATANAME + ".csv";

  const Matrix2d loadStepTransform = getShear(testConfig.loadIncrement);
  std::string dumpPath;
  while (sim->keepLoading()) {
    sim->applyAffineStep(loadStepTransform);
    sim->minimize();
    sim->finishStep(true);
    if (std::abs(sim->mesh.load - 0.4) < 1e-6) {
      dumpPath = sim->saveSimulation("dump_l0.40");
    }
  }
  sim->finishSimulation();
  REQUIRE(!dumpPath.empty());

  auto loadedSim = std::make_shared<Simulation>(testConfig, dataPath);
  loadedSim->addCsvColumn("test_new_column",
                          [](const Simulation &) { return 0.0; });
  Simulation::loadSimulation(*loadedSim, dumpPath, "", dataPath, true);

  {
    std::ifstream file(csvPath);
    REQUIRE(file.is_open());
    std::string line;
    bool foundComment = false;
    size_t commentCols = 0;
    while (std::getline(file, line)) {
      if (line.rfind("#HEADER:", 0) == 0) {
        foundComment = true;
        const std::string headerLine =
            line.substr(std::string("#HEADER:").size());
        commentCols = splitCsvLine(headerLine).size();
        break;
      }
    }
    CHECK(foundComment);
    CHECK(commentCols == loadedSim->getCsvColumns().size());
  }

  loadedSim->maxLoad = 0.8;
  const Matrix2d resumeTransform = getShear(loadedSim->loadIncrement);
  while (loadedSim->keepLoading()) {
    loadedSim->applyAffineStep(resumeTransform);
    loadedSim->minimize();
    loadedSim->finishStep(true);
  }
  loadedSim->finishSimulation();
}

TEST_CASE("Simulation Min CSV Sanity") {
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 3;
  testConfig.cols = 3;
  testConfig.loadIncrement = 0.1;
  testConfig.maxLoad = 0.1;
  testConfig.usingPBC = true;
  testConfig.name = "3x3PBCLoadingMinCsvSanity";
  testConfig.forceReRun = true;
  testConfig.logDuringMinimization = true;
  testConfig.epsR = 1e-6;
  testConfig.initialGuessNoise = 0.1;

  std::string dataPath = "test_data/";

  clearOutputFolder(testConfig.name, dataPath);
  std::shared_ptr<Simulation> s =
      std::make_shared<Simulation>(testConfig, dataPath);
  s->firstStep();
  runSimulationScenario(testConfig, dataPath, s);

  const std::filesystem::path minRoot = std::filesystem::path(dataPath) /
                                        testConfig.name / DATAFOLDERPATH /
                                        MINDATAFOLDERPATH;
  REQUIRE(std::filesystem::exists(minRoot));

  std::vector<std::filesystem::path> minCsvPaths;
  for (const auto &entry :
       std::filesystem::recursive_directory_iterator(minRoot)) {
    if (entry.is_regular_file() &&
        entry.path().filename() == std::string(MACRODATANAME) + ".csv") {
      minCsvPaths.push_back(entry.path());
    }
  }
  REQUIRE(!minCsvPaths.empty());

  for (const auto &csvPath : minCsvPaths) {
    std::ifstream file(csvPath);
    REQUIRE(file.is_open());

    std::string line;
    REQUIRE(std::getline(file, line));
    const auto headers = splitCsvLine(line);

    const int idxTotalEnergy = findHeaderIndex(headers, "total_energy");
    const int idxMinIterTotalChange =
        findHeaderIndex(headers, "min_iter_total_energy_change");
    const int idxMinIterAvgChange =
        findHeaderIndex(headers, "min_iter_avg_energy_change");
    const int idxNrIterations = findHeaderIndex(headers, "nr_iterations");
    const int idxNrFuncEvals = findHeaderIndex(headers, "nr_func_evals");

    REQUIRE(idxTotalEnergy >= 0);
    REQUIRE(idxMinIterTotalChange >= 0);
    REQUIRE(idxMinIterAvgChange >= 0);
    REQUIRE(idxNrIterations >= 0);
    REQUIRE(idxNrFuncEvals >= 0);

    std::vector<double> totalEnergyValues;
    std::vector<double> minIterTotalChangeValues;
    std::vector<double> minIterAvgChangeValues;
    std::vector<double> nrIterationsValues;
    std::vector<double> nrFuncEvalsValues;
    size_t dataLineIndex = 0;

    while (std::getline(file, line)) {
      const auto cells = splitCsvLine(line);
      REQUIRE(cells.size() == headers.size());

      const double totalEnergy = std::stod(cells[idxTotalEnergy]);
      const double minIterTotalChange = std::stod(cells[idxMinIterTotalChange]);
      const double minIterAvgChange = std::stod(cells[idxMinIterAvgChange]);
      const double nrIterations = std::stod(cells[idxNrIterations]);
      const double nrFuncEvals = std::stod(cells[idxNrFuncEvals]);

      totalEnergyValues.push_back(totalEnergy);
      minIterTotalChangeValues.push_back(minIterTotalChange);
      minIterAvgChangeValues.push_back(minIterAvgChange);
      nrIterationsValues.push_back(nrIterations);
      nrFuncEvalsValues.push_back(nrFuncEvals);

      // Function evaluations should increase by one per line (starting at 0).
      CHECK(nrFuncEvals == static_cast<double>(dataLineIndex));
      dataLineIndex++;
    }

    INFO("Min CSV path: " << csvPath.string());
    REQUIRE(totalEnergyValues.size() >= 1);

    bool anyMinIterChange = false;
    bool anyEnergyChange = false;
    for (size_t i = 1; i < totalEnergyValues.size(); ++i) {
      if (std::abs(totalEnergyValues[i] - totalEnergyValues[i - 1]) > 1e-12) {
        anyEnergyChange = true;
      }
      if (std::abs(minIterTotalChangeValues[i]) > 1e-12 ||
          std::abs(minIterAvgChangeValues[i]) > 1e-12) {
        anyMinIterChange = true;
        break;
      }
    }
    if (anyEnergyChange) {
      CHECK(anyMinIterChange);
    }

    // On the final line, if more than one evaluation happened,
    // we expect at least one iteration.
    const double lastFuncEvals = nrFuncEvalsValues.back();
    const double lastIterations = nrIterationsValues.back();
    if (lastFuncEvals > 1) {
      CHECK(lastIterations >= 1);
    }
  }
}

TEST_CASE("Simulation Final Dump Marks Completion") {
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 3;
  testConfig.cols = 3;
  testConfig.usingPBC = true;
  testConfig.name = "3x3FinalDumpCompleteTest";
  testConfig.startLoad = 0.9990;
  testConfig.loadIncrement = 1e-5;
  testConfig.maxLoad = 1.0;
  testConfig.scenario = "simpleShear";

  std::string dataPath = "test_data/";
  clearOutputFolder(testConfig.name, dataPath);

  runSimulationScenario(testConfig, dataPath);

  std::string dumpDir = getDumpPath(testConfig.name, dataPath);
  std::string latestDump;
  std::filesystem::file_time_type latestTime;

  for (const auto &entry : std::filesystem::directory_iterator(dumpDir)) {
    if (!entry.is_regular_file()) {
      continue;
    }
    const auto path = entry.path();
    const auto fileName = path.filename().string();
    const auto filePath = path.string();
    if (fileName.rfind("dump_l", 0) != 0) {
      continue;
    }
    if (filePath.size() < 7 ||
        filePath.substr(filePath.size() - 7) != ".xml.gz") {
      continue;
    }
    auto writeTime = entry.last_write_time();
    if (latestDump.empty() || writeTime > latestTime) {
      latestTime = writeTime;
      latestDump = filePath;
    }
  }

  REQUIRE(!latestDump.empty());

  Simulation loadedSim;
  REQUIRE_THROWS_AS(
      Simulation::loadSimulation(loadedSim, latestDump, "", dataPath, false),
      SimulationAlreadyComplete);
}

TEST_CASE("Simulation Load Handles Old Dumps") {
  // Esure backwards compatibility by loading old dumps from the test folder.
  // Load all dumps in the test folder and check that they can be loaded without
  // error
  std::string testDumpFolderPath = "tests/oldDumps/";

  size_t filesFound = 0;
  for (const auto &entry :
       std::filesystem::recursive_directory_iterator(testDumpFolderPath)) {
    if (!entry.is_regular_file()) {
      continue;
    }
    const auto filePath = entry.path().string();
    const bool isXml =
        filePath.size() >= 4 && filePath.substr(filePath.size() - 4) == ".xml";
    const bool isXmlGz = filePath.size() >= 7 &&
                         filePath.substr(filePath.size() - 7) == ".xml.gz";
    if (!isXml && !isXmlGz) {
      continue;
    }

    filesFound++;
    Simulation loadedSim;
    REQUIRE_NOTHROW(Simulation::loadSimulation(loadedSim, filePath, "",
                                               "test_data/", true));
  }

  REQUIRE(filesFound > 0);
}

// Here, in the main test, we run the simulation in steps and check CSV results
TEST_CASE("Small Simulation Test") {
  // Create a simple config
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.scenario = "simpleShearFixedBoundary";
  testConfig.rows = 2;
  testConfig.cols = 2;
  testConfig.loadIncrement = 0.1;
  testConfig.maxLoad = 0.3;
  testConfig.usingPBC = false;
  testConfig.name = "2x2FixedLoading";

  // Create a data path and file paths
  std::string dataPath = "test_data/";

  // Remove old data
  clearOutputFolder(testConfig.name, dataPath);
  std::shared_ptr<Simulation> s;

  // Run the scenario and check CSV
  runSimulationScenario(testConfig, dataPath, s);
}

TEST_CASE("3x3 Periodic Mesh Simple Shear Stress Test") {
  // Create a simple config
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 3;
  testConfig.cols = 3;
  testConfig.loadIncrement = 0.05;
  testConfig.startLoad = -testConfig.loadIncrement;
  testConfig.maxLoad = 3.0;
  testConfig.initialGuessNoise = 0;
  testConfig.usingPBC = true;
  testConfig.name = "3x3PeriodicShearStressTest";
  testConfig.forceReRun = true;

  // Create a data path and file paths
  std::string dataPath = "test_data/";

  // Remove old data
  clearOutputFolder(testConfig.name, dataPath);

  // Create simulation
  Simulation sim(testConfig, dataPath);
  sim.firstStep();

  struct StressSample {
    double load = 0.0;
    double p21 = 0.0;
    double p12 = 0.0;
    double sigma21 = 0.0;
    double sigma12 = 0.0;
    double avgPXy = 0.0;
    double avgSigmaXy = 0.0;
  };
  std::vector<StressSample> samples;

  const double parityTol = 1e-8;
  bool p21ParityOk = true;
  bool p12ParityOk = true;
  bool sigma21ParityOk = true;
  bool sigma12ParityOk = true;
  // Run simulation manually to track stress at each step
  while (sim.keepLoading()) {
    // Apply load step
    Matrix2d loadTransform = getShear(testConfig.loadIncrement);
    sim.mesh.addLoad(testConfig.loadIncrement);
    sim.mesh.applyTransformationToSystemDeformation(loadTransform);
    sim.applyLoadStepToGuess(loadTransform);

    // Update mesh properties
    sim.mesh.updateMesh();
    sim.mesh.updateAveragesAndPlasticEvents();

    // Check that all elements of the same parity have the same P and sigma
    const auto &refEven = sim.mesh.elements[0];
    const auto &refOdd = sim.mesh.elements[1];
    for (size_t i = 0; i < sim.mesh.elements.size(); ++i) {
      const auto &e = sim.mesh.elements[i];
      const auto &ref = (i % 2 == 0) ? refEven : refOdd;
      if (std::abs(e.P(1, 0) - ref.P(1, 0)) >= parityTol) {
        p21ParityOk = false;
      }
      if (std::abs(e.P(0, 1) - ref.P(0, 1)) >= parityTol) {
        p12ParityOk = false;
      }
      if (std::abs(e.sigma(1, 0) - ref.sigma(1, 0)) >= parityTol) {
        sigma21ParityOk = false;
      }
      if (std::abs(e.sigma(0, 1) - ref.sigma(0, 1)) >= parityTol) {
        sigma12ParityOk = false;
      }
    }

    const auto &e0 = sim.mesh.elements[0];
    samples.push_back({sim.mesh.load, e0.P(1, 0), e0.P(0, 1), e0.sigma(1, 0),
                       e0.sigma(0, 1), sim.mesh.averageP12,
                       sim.mesh.averageSigma12});
    // std::cout << "load " << sim.mesh.load << " avgPxy " <<
    // sim.mesh.averageP12
    //           << " avgSigmaXy " << sim.mesh.averageSigma12 << '\n';
  }
  CHECK(p21ParityOk);
  CHECK(p12ParityOk);
  CHECK(sigma21ParityOk);
  CHECK(sigma12ParityOk);

  auto stepIndex = [&](double load) -> int {
    return static_cast<int>(std::llround(load / testConfig.loadIncrement));
  };

  auto findSample = [&](double target) -> const StressSample * {
    const int targetStep = stepIndex(target);
    for (const auto &sample : samples) {
      if (stepIndex(sample.load) == targetStep) {
        return &sample;
      }
    }
    return nullptr;
  };

  const double zeroTol = 1e-6;
  for (double targetLoad : {0.0, 1.0, 2.0, 3.0}) {
    const StressSample *sample = findSample(targetLoad);
    REQUIRE(sample);
    CHECK(std::abs(sample->p21) < zeroTol);
    CHECK(std::abs(sample->p12) < zeroTol);
    CHECK(std::abs(sample->sigma21) < zeroTol);
    CHECK(std::abs(sample->sigma12) < zeroTol);
  }

  struct ExpectedStress {
    double load = 0.0;
    double p21 = 0.0;
    double p12 = 0.0;
    double sigma21 = 0.0;
    double sigma12 = 0.0;
  };

  // Generated using MTMath/miniMTM.py
  const std::vector<ExpectedStress> expected = {
      {0.3, 0.179775, 0.179628, 0.179628, 0.179628},
      {0.5, -0.0231268, 0.0, 0.0, 0.0},
      {1.3, 0.180266, 0.179628, 0.179628, 0.179628},
      {1.5, -0.0693803, 0.0, 0.0, 0.0},
      {2.3, 0.180757, 0.179628, 0.179628, 0.179628},
      {2.5, -0.115634, 0.0, 0.0, 0.0},
  };

  for (const auto &exp : expected) {
    const StressSample *sample = findSample(exp.load);
    REQUIRE(sample);
    CHECK(sample->p21 == doctest::Approx(exp.p21).epsilon(1e-6));
    CHECK(sample->p12 == doctest::Approx(exp.p12).epsilon(1e-6));
    CHECK(sample->sigma21 == doctest::Approx(exp.sigma21).epsilon(1e-6));
    CHECK(sample->sigma12 == doctest::Approx(exp.sigma12).epsilon(1e-6));
  }
}

struct SimulationTestAccess {
  static Mesh::Snapshot &before(Simulation &s) { return s.beforeMinimization; }
  static Mesh::Snapshot &after(Simulation &s) { return s.afterMinimization; }
};

TEST_CASE("Participation fraction Test") {
  // Create a simple config
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 3;
  testConfig.cols = 3;
  testConfig.loadIncrement = 0.05;
  testConfig.maxLoad = 1.0;
  testConfig.initialGuessNoise = 0;
  testConfig.usingPBC = true;
  testConfig.name = "participationFractionTest";
  testConfig.forceReRun = true;

  // Create a data path and file paths
  std::string dataPath = "test_data/";

  // Remove old data
  clearOutputFolder(testConfig.name, dataPath);

  // Create simulation
  Simulation sim(testConfig, dataPath);
  sim.initialize();

  const size_t nodeCount = static_cast<size_t>(sim.mesh.rows * sim.mesh.cols);
  REQUIRE(nodeCount > 0);

  auto &before = SimulationTestAccess::before(sim);
  auto &after = SimulationTestAccess::after(sim);
  sim.mesh.saveSnapshot(before);
  sim.mesh.saveSnapshot(after);

  // Uniform displacement across all free nodes -> participation fraction = 1
  const double ux = 0.1;
  const double uy = -0.2;
  for (size_t i = 0; i < nodeCount; ++i) {
    after.displacements[i] = before.displacements[i] + ux;
    after.displacements[i + nodeCount] =
        before.displacements[i + nodeCount] + uy;
  }
  sim.computeParticipationFraction();
  const double uniformPf = sim.participationFraction;
  CHECK(uniformPf == doctest::Approx(1.0));

  // Displacement only on one node -> participation fraction = 1/N
  after.displacements = before.displacements;
  after.displacements[0] = before.displacements[0] + ux;
  after.displacements[nodeCount] = before.displacements[nodeCount] + uy;
  sim.computeParticipationFraction();
  const double localizedPf = sim.participationFraction;
  CHECK(localizedPf == doctest::Approx(1.0 / static_cast<double>(nodeCount)));
}
