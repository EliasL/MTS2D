#include "../src/Simulation/simulation.h" // Include the header for your surface struct
#include "run/doctest.h"

TEST_CASE("Simulation Save/Load Test") {
  // Create a simple config
  Config testConfig;
  testConfig.setDefaultValues();
  testConfig.rows = 10;
  testConfig.cols = 10;
  testConfig.loadIncrement = 0.1;
  testConfig.maxLoad = 0.2;

  // Create and initialize simulation
  std::string dataPath = "test_data/";
  Simulation sim(testConfig, dataPath);
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
  std::string pathToSave = sim.saveSimulation(saveFileName);

  // Load simulation into a new object
  Simulation loadedSim;
  Simulation::loadSimulation(loadedSim, pathToSave, "", dataPath, true);

  double loadedEnergy = sim.mesh.totalEnergy;
  CHECK(doctest::Approx(loadedEnergy).epsilon(1e-12) == originalEnergy);

  // Update properties
  sim.mesh.updateMesh();
  loadedEnergy = sim.mesh.totalEnergy;

  CHECK(doctest::Approx(loadedEnergy).epsilon(1e-12) == originalEnergy);
}