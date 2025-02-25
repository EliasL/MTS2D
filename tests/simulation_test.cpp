#include "../src/Simulation/scenarios.h"
#include "../src/Simulation/simulation.h" // Include the header for your surface struct
#include "Data/data_export.h"
#include "Mesh/mesh.h"
#include "Mesh/tElement.h"
#include "run/doctest.h"
#include <iostream>
#include <string>

TEST_CASE("Simulation Save/Load mesh Test") {
  // Create a simple config
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 2;
  testConfig.cols = 2;
  testConfig.loadIncrement = 0.1;
  testConfig.maxLoad = 0.2;
  testConfig.usingPBC = true;

  // Create and initialize simulation
  std::string dataPath = "test_data";
  Simulation sim(testConfig, dataPath, true);
  sim.initialize();

  Matrix2d loadStepTransform = getShear(testConfig.loadIncrement);
  sim.setInitialGuess(loadStepTransform);

  // Run a simulation step
  sim.mesh.addLoad(sim.loadIncrement);
  sim.minimize();

  sim.mesh.updateMesh();
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

  // Update properties (should not change the result)
  loadedSim.mesh.updateMesh();

  CHECK(loadedSim.mesh == sim.mesh);
  if (loadedSim.mesh != sim.mesh) {
    std::cout << debugCompare(loadedSim.mesh, sim.mesh) << std::endl;
  }
  // abort further testing
}

TEST_CASE("Simulation Save/Load Energy Test") {
  // Create a simple config
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 10;
  testConfig.cols = 10;
  testConfig.loadIncrement = 0.1;
  testConfig.maxLoad = 0.2;

  // Create and initialize simulation
  std::string dataPath = "test_data";
  Simulation sim(testConfig, dataPath, true);
  sim.initialize();

  Matrix2d loadStepTransform = getShear(testConfig.loadIncrement);
  sim.setInitialGuess(loadStepTransform);

  // Run a simulation step
  sim.mesh.addLoad(sim.loadIncrement);
  sim.minimize();

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
  loadedEnergy = loadedSim.mesh.totalEnergy;

  CHECK(doctest::Approx(loadedEnergy).epsilon(1e-12) == originalEnergy);
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

// Here, in the main test, we run the simulation in steps and check CSV results
TEST_CASE("Simulation Save/Load Macro Data Test") {
  // Create a simple config
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 3;
  testConfig.cols = 3;
  testConfig.loadIncrement = 0.1;
  testConfig.maxLoad = 0.3;

  // Create a data path and file paths
  std::string dataPath = "test_data/";
  std::string dumpPath = dataPath + "default_name/dumps/dump_l0.20.xml.gz";
  std::string csvPath = dataPath + "default_name/macroData.csv";

  // Create and initialize simulation
  Simulation sim(testConfig, dataPath);

  // Run the scenario and check CSV
  runSimulationScenario(testConfig, dataPath);

  // Check that the first column is 1, 2, 3
  checkMacroDataCsv(csvPath, {1, 2, 3, 4});

  // Load simulation into a new object
  using SimPtr = std::shared_ptr<Simulation>;
  SimPtr loadedSim = std::make_shared<Simulation>(testConfig, dataPath);
  Simulation::loadSimulation(*loadedSim, dumpPath, "", dataPath, true);

  // After loading from l0.2, check that the first column is 1, 2
  checkMacroDataCsv(csvPath, {1, 2, 3});

  // Increase max load, run again, and check appended results
  loadedSim->maxLoad = 0.4;
  runSimulationScenario(testConfig, dataPath, loadedSim);

  // Now, the first column should be 1, 2, 3, 4
  checkMacroDataCsv(csvPath, {1, 2, 3, 4, 5});
}
