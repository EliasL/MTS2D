#include "Data/cereal_help.h"
#include "Data/data_export.h"
#include "Simulation/protocols.h"
#include "element_history_export.h"
#include "run/doctest.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cereal/archives/json.hpp>

namespace {
Config protocolConfig(const std::string &name) {
  Config config;
  config.setDefaultValues();
  config.name = name;
  config.rows = config.cols = 4;
  config.loading = "relaxation";
  config.deformation = "simpleShear";
  config.boundary = "periodic";
  config.startLoad = 0.01;
  config.loadIncrement = 0.01;
  config.totalLoad = 0.02;
  config.initialGuessNoise = 0;
  config.nrThreads = 1;
  config.showProgress = -1;
  config.writeMeshVTUs = false;
  config.writeDumps = false;
  config.forceReRun = true;
  return config;
}

double firstVtuScalar(const std::string &xml, const std::string &name) {
  const auto field = xml.find("Name=\"" + name + "\"");
  if (field == std::string::npos) {
    throw std::runtime_error("Missing VTU field " + name);
  }
  return std::stod(xml.substr(xml.find('>', field) + 1));
}
} // namespace

TEST_CASE("Protocol options have shared defaults and reject removed selectors") {
  Config defaults;
  defaults.setDefaultValues();
  const Config parsed = initializeConfig({});
  CHECK(parsed.loadIncrement == defaults.loadIncrement);
  CHECK(parsed.rows == defaults.rows);
  CHECK(parsed.totalLoad == defaults.totalLoad);
  CHECK(initializeConfig({{"repeat", "false"}}).totalLoad == defaults.totalLoad);
  CHECK(initializeConfig({{"startLoad", "0.15"}, {"totalLoad", "0.01"}}).totalLoad == 0.01);
  CHECK_THROWS_AS(initializeConfig({{"maxLoad", "1"}}), std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"loading", "affineOnly"}}), std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"deformation", "directDisplacement"},
                                   {"boundary", "periodic"}}), std::invalid_argument);
  const Config config = initializeConfig({
      {"loading", "reversibility"}, {"deformation", "simpleShear"},
      {"boundary", "fixed"}, {"directions", "vertical"},
      {"minimum_transport_distance", "false"},
  });
  CHECK(config.loading == "reversibility");
  CHECK(config.deformation == "simpleShear");
  CHECK(config.directions == "vertical");
  CHECK_FALSE(config.usesPeriodicBoundary());
  CHECK_FALSE(config.minimum_transport_distance);
  CHECK_THROWS_AS(initializeConfig({{"protocol", "standard"}}), std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"experiment", "doubleDislocationTest"}}), std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"loading", "reversibilityProtocolTest"}}), std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"directions", "horizontal,vertical"}}), std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"directions", "horizontal,"}}), std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"repeat", "true"}}), std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"loadIncrement", "0"}}), std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"boundary", "unknown"}}), std::invalid_argument);
  CHECK_THROWS_AS(initializeConfig({{"directions", "diagonal"}}), std::invalid_argument);
  CHECK_FALSE(initializeConfig({{"loading", "noRelaxation"},
                                   {"affineInitialGuess", "false"}}).affineInitialGuess);
  CHECK(initializeConfig({{"loading", "reversibility"},
                                   {"stepNoise", "0.01"}}).stepNoise == 0.01);
}

TEST_CASE("Protocol selections and maximum flip-distance survive a dump") {
  Config config = protocolConfig("protocolSerialization");
  config.loading = "reversibility";
  config.directions = "vertical,-horizontal";
  config.segmentLoad = 0.02;
  config.repeat = true;
  config.checkpointInterval = 7;
  config.stopAfterEnergyDrops = 1;
  config.reversibilityTolerance = 2e-5;
  config.minimum_transport_distance = false;
  const std::string xml = serializeToXml(config);
  Config loaded;
  deserializeFromXml(xml, loaded);
  CHECK(loaded.totalLoad == config.totalLoad);
  CHECK(loaded.deformation == config.deformation);
  CHECK(loaded.loading == config.loading);
  CHECK(loaded.boundary == config.boundary);
  CHECK(loaded.directions == config.directions);
  CHECK(loaded.segmentLoad == config.segmentLoad);
  CHECK(loaded.repeat == config.repeat);
  CHECK(loaded.checkpointInterval == config.checkpointInterval);
  CHECK(loaded.stopAfterEnergyDrops == config.stopAfterEnergyDrops);
  CHECK(loaded.reversibilityTolerance == config.reversibilityTolerance);
  CHECK(loaded.usesPeriodicBoundary());
  CHECK_FALSE(loaded.minimum_transport_distance);
}

TEST_CASE("Shared protocols compose boundaries and reconnection choices") {
  for (const std::string protocol : {"relaxation", "reversibility", "noRelaxation"}) {
    for (const std::string boundary : {"periodic", "fixed"}) {
      for (const std::string reconnection : {"none", "edgeFlip"}) {
        for (const bool minimum : {false, true}) {
          CAPTURE(protocol);
          CAPTURE(boundary);
          CAPTURE(reconnection);
          CAPTURE(minimum);
          Config config = protocolConfig("protocols_" + protocol + "_" +
                                         boundary + "_" + reconnection +
                                         (minimum ? "_min" : "_max"));
          config.loading = protocol;
          config.boundary = boundary;
          config.reconnectionMethod = reconnection;
          config.minimum_transport_distance = minimum;
          auto simulation = runSimulationProtocol(config, "test_data");
          CHECK(simulation->accumulatedLoad() == doctest::Approx(config.totalLoad));
          CHECK(simulation->mesh.usingPBC == (boundary == "periodic"));
          CHECK(simulation->mesh.fixedNodeIds.size() ==
                (boundary == "periodic" ? 0 : 12));
          CHECK(simulation->mesh.minimum_transport_distance == minimum);
          simulation->mesh.ensureFull();
          const Matrix2d expected = getShear(config.startLoad + config.totalLoad);
          bool affine = true;
          for (const auto &element : simulation->mesh.elements) {
            affine &= element.F.isApprox(expected, 1e-8);
          }
          CHECK(affine);
        }
      }
    }
  }
}

TEST_CASE("Direction schedules split steps and compose initial and final deformations") {
  struct Example {
    std::string directions;
    double segmentLoad;
    bool repeat;
    double startLoad;
    double endLoad;
    double increment;
    Matrix2d expected;
  };
  const auto vertical = [](double load) -> Matrix2d {
    return getShear(load).transpose();
  };
  const std::vector<Example> examples = {
      {"vertical", 0, false, 0.01, 0.035, -0.017, vertical(-0.035)},
      {"horizontal,vertical", 0.02, false, 0.025, 0.065, 0.017,
       vertical(0.045) * getShear(0.02)},
      {"horizontal,-horizontal", 0.02, true, 0, 0.07, 0.017, getShear(0.01)},
      {"horizontal,vertical", 0.02, true, 0.085, 0.095, 0.017,
       getShear(0.015) * vertical(0.02) * getShear(0.02) * vertical(0.02) * getShear(0.02)},
  };
  for (size_t i = 0; i < examples.size(); ++i) {
    CAPTURE(i);
    const Example &example = examples[i];
    Config config = protocolConfig("directionSchedule" + std::to_string(i));
    config.loading = "noRelaxation";
    config.directions = example.directions;
    config.segmentLoad = example.segmentLoad;
    config.repeat = example.repeat;
    config.startLoad = example.startLoad;
    config.totalLoad = example.endLoad - example.startLoad;
    config.loadIncrement = example.increment;
    auto simulation = runSimulationProtocol(config, "test_data");
    CHECK(simulation->mesh.load == example.endLoad);
    CHECK(simulation->mesh.currentDeformation.isApprox(example.expected, 1e-12));
    simulation->mesh.ensureFull();
    for (const auto &element : simulation->mesh.elements) {
      CHECK(element.F.isApprox(example.expected, 1e-10));
    }
    CHECK(protocolDeformation(config, config.startLoad + config.totalLoad).isApprox(example.expected, 1e-12));
    if (i == 1) CHECK(simulation->mesh.loadSteps == 4);
  }
}

TEST_CASE("An energy-drop stop keeps the preceding checkpoint") {
  Config config = protocolConfig("energyDropCheckpoint");
  config.loading = "noRelaxation";
  config.directions = "horizontal,-horizontal";
  config.segmentLoad = 0.01;
  config.loadIncrement = 0.005;
  config.checkpointInterval = 1;
  config.stopAfterEnergyDrops = 1;
  config.stepEnergyTolerance = 0;
  auto simulation = runSimulationProtocol(config, "test_data");
  CHECK(simulation->mesh.load == doctest::Approx(0.015));
  const auto folder = std::filesystem::path("test_data") / config.name / "dumps";
  CHECK(std::filesystem::exists(folder / "checkpoint_1.xml.gz"));
  CHECK(simulation->termination.energyDrops == 1);
  CHECK(simulation->termination.stopReason == "sample_targets");
  CHECK(std::filesystem::exists(folder.parent_path() / "termination.json"));
  Simulation before;
  Simulation::loadSimulation(before, (folder / "checkpoint_1.xml.gz").string(),
                             "", "test_data", true);
  CHECK(before.mesh.load == doctest::Approx(0.01));
}

TEST_CASE("Unified all and any targets stop after completed loading steps") {
  for (const std::string mode : {"all", "any"}) {
    Config config = protocolConfig("stepTargets" + mode);
    config.loading = "noRelaxation";
    config.totalLoad = 0.05;
    config.stopAfterSteps = 2;
    config.stopAfterEnergyIncreases = 3;
    config.stopWhen = mode;
    auto simulation = runSimulationProtocol(config, "test_data");
    const int expected = mode == "all" ? 3 : 2;
    CHECK(simulation->termination.steps == expected);
    CHECK(simulation->termination.energyIncreases == expected);
    CHECK(simulation->mesh.loadSteps == expected + 1); // Exclude initialization.
    CHECK(simulation->termination.targetsMet);
    CHECK(simulation->termination.stopReason == "sample_targets");
    std::ifstream file(std::filesystem::path("test_data") / config.name / "termination.json");
    cereal::JSONInputArchive archive(file);
    StepTermination saved;
    archive(cereal::make_nvp("termination", saved));
    CHECK(saved.steps == expected);
    CHECK(saved.targetsMet);
  }
}

TEST_CASE("A step target caps an unmet plastic target in any mode") {
  Config config = protocolConfig("plasticStepCap");
  config.startLoad = 0.15;
  config.loadIncrement = 1e-5;
  config.totalLoad = 1e-4;
  config.stopAfterSteps = 3;
  config.stopAfterPlasticSteps = 20;
  config.stopWhen = "any";
  auto simulation = runSimulationProtocol(config, "test_data");
  CHECK(simulation->termination.steps == 3);
  CHECK(simulation->termination.plasticSteps == 0);
  CHECK(simulation->termination.elasticSteps == 3);
  CHECK(simulation->termination.targetsMet);
  CHECK(simulation->mesh.load == doctest::Approx(0.15003));
}

TEST_CASE("Strain exhaustion and resume preserve cumulative step targets") {
  Config config = protocolConfig("stepTargetResume");
  config.loading = "noRelaxation";
  config.totalLoad = 0.045;
  config.stopAfterSteps = 2;
  auto simulation = runSimulationProtocol(config, "test_data");
  const std::string dump = simulation->saveSimulation("resumeStepTargets");
  auto resumed = std::make_shared<Simulation>();
  Simulation::loadSimulation(*resumed, dump, "", "test_data", true);
  resumed->config.stopAfterSteps = 100;
  resumed = runSimulationProtocol(resumed->config, "test_data", resumed);
  CHECK(resumed->termination.steps == 5); // Includes the shortened final step.
  CHECK(resumed->mesh.loadSteps == 6);
  CHECK_FALSE(resumed->termination.targetsMet);
  CHECK(resumed->termination.stopReason == "strain_limit");
  CHECK(resumed->mesh.load == doctest::Approx(0.055));
}

TEST_CASE("Resuming the new reversibility protocol preserves its trajectory") {
  Config config = protocolConfig("newProtocolResume");
  config.loading = "reversibility";
  config.startLoad = 0;
  config.directions = "horizontal,vertical,-horizontal,-vertical";
  config.segmentLoad = 0.02;
  config.repeat = true;
  config.loadIncrement = 0.007;
  config.totalLoad = 0.034;
  auto first = runSimulationProtocol(config, "test_data");
  const std::string dump = first->saveSimulation("protocol_checkpoint");
  first->csvFile.close();

  auto resumed = std::make_shared<Simulation>();
  Simulation::loadSimulation(*resumed, dump, "", "test_data", true);
  resumed->config.totalLoad = resumed->totalLoad = 0.091;
  runSimulationProtocol(resumed->config, "test_data", resumed);

  config.name = "newProtocolUninterrupted";
  config.totalLoad = 0.091;
  auto uninterrupted = runSimulationProtocol(config, "test_data");
  CHECK(resumed->mesh.load == doctest::Approx(uninterrupted->mesh.load));
  CHECK(resumed->mesh.rmsDistanceToMesh(uninterrupted->mesh) < 1e-10);
  CHECK(resumed->mesh.totalEnergy == doctest::Approx(uninterrupted->mesh.totalEnergy));
  bool hasReversibilityColumn = false;
  for (const auto &column : resumed->getCsvColumns()) {
    if (column.name == "is_reversible") {
      hasReversibilityColumn = true;
    }
  }
  CHECK(hasReversibilityColumn);
}

TEST_CASE("Resume completes a final step shorter than half the configured increment") {
  Config config = protocolConfig("shortFinalStepResume");
  config.loading = "noRelaxation";
  config.startLoad = 0;
  config.loadIncrement = 0.02;
  config.totalLoad = 0.03;
  auto first = runSimulationProtocol(config, "test_data");
  first->config.totalLoad = first->totalLoad = 0.031;
  const std::string dump = first->saveSimulation("short_final_checkpoint");
  first->csvFile.close();

  auto resumed = std::make_shared<Simulation>();
  Simulation::loadSimulation(*resumed, dump, "", "test_data", false);
  runSimulationProtocol(resumed->config, "test_data", resumed);
  CHECK(resumed->mesh.load == 0.031);
  CHECK(resumed->mesh.currentDeformation.isApprox(getShear(0.031), 1e-12));
}

TEST_CASE("VTU and history exports distinguish total and plastic deformation") {
  Mesh mesh(2, 2, false, "major");
  mesh.applyTransformation(getShear(0.2));
  mesh.ensureFull();
  for (auto &element : mesh.elements) {
    element.T_p = getShear(2.0);
  }
  const auto row = captureElementMatrixHistoryRow(mesh, 0);
  CHECK(row.T_p.isApprox(getShear(2.0), 1e-12));
  CHECK(row.T.isApprox(getShear(2.2), 1e-12));
  CHECK(row.T.isApprox(row.F_E * row.T_p, 1e-12));
  const std::string file = writeMeshToVtu(mesh, "tensorNames", "test_data", "fields");
  std::ifstream input(file);
  std::ostringstream contents;
  contents << input.rdbuf();
  CHECK(firstVtuScalar(contents.str(), "T_p12") == doctest::Approx(2.0));
  CHECK(firstVtuScalar(contents.str(), "T12") == doctest::Approx(2.2));
}

TEST_CASE("Total loading excludes initialization and permits an initialization-only run") {
  Config config = protocolConfig("loadBudget");
  config.loading = "noRelaxation";
  config.startLoad = 0.15;
  config.totalLoad = 0.025;
  config.loadIncrement = -0.017;
  auto simulation = runSimulationProtocol(config, "test_data");
  CHECK(simulation->accumulatedLoad() == doctest::Approx(0.025));
  CHECK(simulation->mesh.currentDeformation.isApprox(getShear(-0.175), 1e-12));
  config.name = "initializationOnly";
  config.totalLoad = 0;
  config.writeDumps = config.writeMeshVTUs = true;
  simulation = runSimulationProtocol(config, "test_data");
  CHECK(simulation->accumulatedLoad() == 0);
  CHECK(simulation->mesh.currentDeformation.isApprox(getShear(-0.15), 1e-12));
  CHECK(std::filesystem::exists(simulation->saveSimulation()));
}

TEST_CASE("Direct displacement uses stable driven groups and the shared direction schedule") {
  Config config = protocolConfig("directDisplacementSchedule");
  config.loading = "noRelaxation";
  config.deformation = "directDisplacement";
  config.boundary = "periodicAnchored";
  config.startLoad = 0.085;
  config.totalLoad = 0.01;
  config.loadIncrement = 0.007;
  config.segmentLoad = 0.02;
  config.directions = "horizontal,vertical";
  config.repeat = true;
  auto simulation = runSimulationProtocol(config, "test_data");
  CHECK(simulation->mesh.fixedNodeIds.size() == 7);
  CHECK(simulation->drivenNodeGroups[0].size() == 2);
  CHECK(simulation->drivenNodeGroups[1].size() == 2);
  CHECK(simulation->accumulatedLoad() == doctest::Approx(config.totalLoad));
  for (int row = 0; row < config.rows; ++row) {
    for (int col = 0; col < config.cols; ++col) {
      const Node &node = simulation->mesh.nodes(row, col);
      const Vector2d expected(col == 0 && row >= 2 ? 0.055 : 0,
                              row == 0 && col >= 2 ? 0.04 : 0);
      CHECK(node.u().isApprox(expected, 1e-12));
    }
  }
  const Mesh forward = simulation->mesh;
  simulation->loadIncrement = 2;
  simulation->applyDeformation({Direction::Horizontal, 2});
  simulation->loadIncrement = -2;
  simulation->applyDeformation({Direction::Horizontal, -2});
  CHECK(simulation->mesh.rmsDistanceToMesh(forward) < 1e-12);
}

TEST_CASE("Direct displacement composes loading modes and resumes with saved driven groups") {
  for (const std::string loading : {"relaxation", "reversibility", "noRelaxation"}) {
    for (const std::string boundary : {"fixed", "anchored", "periodicAnchored"}) {
      CAPTURE(loading);
      CAPTURE(boundary);
      Config config = protocolConfig("directResume_" + loading + "_" + boundary);
      // Resume rebuilds the factor. Resolve the physical state tightly enough
      // to compare paths that use different cached preconditioners.
      config.epsR = 1e-9;
      config.loading = loading;
      config.deformation = "directDisplacement";
      config.boundary = boundary;
      config.startLoad = 0.01;
      config.totalLoad = 0.024;
      config.loadIncrement = 0.007;
      config.segmentLoad = 0.02;
      config.directions = "horizontal,vertical,-horizontal,-vertical";
      config.repeat = true;
      // Force the backward check even for these small elastic steps.
      config.saveElasticReversibilityStates = loading == "reversibility";
      config.maximumSavedElasticReversibilityStates = 1;
      auto partial = runSimulationProtocol(config, "test_data");
      const auto groups = partial->drivenNodeGroups;
      const std::string dump = partial->saveSimulation("direct_checkpoint");
      partial->csvFile.close();
      auto resumed = std::make_shared<Simulation>();
      Simulation::loadSimulation(*resumed, dump, "", "test_data", true);
      for (size_t axis = 0; axis < groups.size(); ++axis) {
        REQUIRE(resumed->drivenNodeGroups[axis].size() == groups[axis].size());
        for (size_t i = 0; i < groups[axis].size(); ++i) {
          CHECK(resumed->drivenNodeGroups[axis][i].i == groups[axis][i].i);
        }
      }
      resumed->config.totalLoad = resumed->totalLoad = 0.081;
      runSimulationProtocol(resumed->config, "test_data", resumed);
      config.name += "_uninterrupted";
      config.totalLoad = 0.081;
      auto uninterrupted = runSimulationProtocol(config, "test_data");
      CHECK(resumed->accumulatedLoad() == doctest::Approx(0.081));
      CHECK(resumed->mesh.rmsDistanceToMesh(uninterrupted->mesh) < 1e-8);
      CHECK(resumed->mesh.totalEnergy == doctest::Approx(uninterrupted->mesh.totalEnergy));
      if (loading == "reversibility") {
        CHECK(resumed->mesh.nodes(3, 0).u().isApprox(Vector2d(0.011, 0), 1e-10));
        CHECK(resumed->mesh.nodes(0, 3).u().norm() < 1e-10);
        const auto output = std::filesystem::path("test_data") / config.name;
        size_t snapshots = 0;
        for (const auto &entry : std::filesystem::recursive_directory_iterator(output)) {
          if (entry.path().filename().string().find("state3_deformed") != std::string::npos) ++snapshots;
        }
        CHECK(snapshots > 0);
      }
    }
  }
}
