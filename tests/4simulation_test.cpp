#include "../src/Simulation/scenarios.h"
#include "../src/Simulation/simulation.h" // Include the header for your surface struct
#include "Data/data_export.h"
#include "Mesh/mesh.h"
#include "Mesh/tElement.h"
#include "run/doctest.h"
#include <filesystem>
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
  testConfig.name = "2x2PBCSaveLoadTest";

  // Create and initialize simulation
  std::string dataPath = "test_data";
  Simulation sim(testConfig, dataPath, true);
  sim.initialize();

  Matrix2d loadStepTransform = getShear(testConfig.loadIncrement);
  sim.setInitialGuess(loadStepTransform);

  // Run a simulation step
  sim.mesh.addLoad(sim.loadIncrement);
  sim.mesh.applyTransformationToSystemDeformation(loadStepTransform);
  sim.minimize();

  sim.mesh.updateMesh(true);
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
  testConfig.showProgress = 0;
  testConfig.forceReRun = true;

  // Create a data path and file paths
  std::string dataPath = "test_data/";
  std::string dumpPath = dataPath + testConfig.name + "/dumps/dump_l1.0.xml.gz";

  // Remove old data
  clearOutputFolder(testConfig.name, dataPath);
  std::shared_ptr<Simulation> s =
      std::make_shared<Simulation>(testConfig, dataPath);
  // s->mesh.fixBorderNodes();
  s->initialize();
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
  sim.setInitialGuess(loadStepTransform);

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
  testConfig.usingPBC = true;
  testConfig.name = "3x3PBCLoadingTestWithSaveLoad";
  testConfig.showProgress = 0;

  // Create a data path and file paths
  std::string dataPath = "test_data/";
  std::string dumpPath =
      dataPath + testConfig.name + "/dumps/dump_l0.20.xml.gz";
  std::string csvPath = dataPath + testConfig.name + "/macroData.csv";

  // Remove old data
  clearOutputFolder(testConfig.name, dataPath);
  std::shared_ptr<Simulation> s =
      std::make_shared<Simulation>(testConfig, dataPath);
  s->initialize();
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

  // After loading from l0.2, check that the first column is 1, 2, 3
  // (Corresponding to load 0.0, 0.1, and 0.2)
  checkMacroDataCsv(csvPath, {1, 2, 3});

  // Increase max load, run again, and check appended results
  loadedSim->maxLoad = 0.4;
  runSimulationScenario(testConfig, dataPath, loadedSim);

  // Now, the first column should be 1, 2, 3, 4, 5
  checkMacroDataCsv(csvPath, {1, 2, 3, 4, 5});
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
  testConfig.showProgress = 0;

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
  std::string testDumpFolderPath = "../tests/oldDumps/";

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
  testConfig.showProgress = 0;

  // Create a data path and file paths
  std::string dataPath = "test_data/";

  // Remove old data
  clearOutputFolder(testConfig.name, dataPath);
  std::shared_ptr<Simulation> s;

  // Run the scenario and check CSV
  runSimulationScenario(testConfig, dataPath, s);
}
