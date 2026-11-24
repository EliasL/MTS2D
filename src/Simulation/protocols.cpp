#include "protocol_deformation.h"
#include "Data/data_export.h"
#include "finite.h"
#include "randomUtils.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <stdexcept>

namespace {
enum class Loading { Relaxation, Reversibility, NoRelaxation };

Loading loadingChoice(const std::string &name) {
  if (name == "relaxation") {
    return Loading::Relaxation;
  }
  if (name == "reversibility") {
    return Loading::Reversibility;
  }
  if (name == "noRelaxation") {
    return Loading::NoRelaxation;
  }
  throw std::invalid_argument("Unknown loading: " + name +
      ". Choose relaxation, reversibility, or noRelaxation.");
}

struct ScheduledStep {
  SignedDirection direction;
  double amount;
};

ScheduledStep nextStep(const Config &config,
                       const std::vector<SignedDirection> &directions, double load,
                       double endLoad) {
  double amount = std::min(std::abs(config.loadIncrement), endLoad - load);
  if (config.segmentLoad == 0) {
    return {directions.front(), amount};
  }

  // Find the segment from the saved path coordinate. The small tolerance keeps
  // roundoff at a segment boundary from selecting the segment just completed.
  const double position = load / config.segmentLoad;
  const double segment = std::floor(position + 16 * DBL_EPSILON *
                                               std::max(1.0, std::abs(position)));
  const double last = static_cast<double>(directions.size() - 1);
  const size_t index = config.repeat
      ? static_cast<size_t>(std::fmod(std::fmod(segment, directions.size()) + directions.size(), directions.size()))
      : static_cast<size_t>(std::max(0.0, std::min(segment, last)));
  if (config.repeat || segment < last) {
    amount = std::min(amount, (segment + 1) * config.segmentLoad - load);
  }
  return {directions[index], amount};
}

} // namespace

void validateProtocolConfig(const Config &config) {
  if (!config.initialMeshSha256.empty()) {
    if (config.initialMeshSha256.size() != 64 ||
        config.initialMeshSha256.find_first_not_of("0123456789abcdef") != std::string::npos) {
      throw std::invalid_argument("initialMeshSha256 must be a lowercase SHA-256 checksum.");
    }
  }
  loadingChoice(config.loading);
  if (config.deformation != "simpleShear" && config.deformation != "directDisplacement" &&
      config.deformation != "sequentialSlip" && !isAxialDeformation(config) &&
      !config.usesIndentation()) {
    throw std::invalid_argument("Unknown deformation: " + config.deformation);
  }
  if (config.boundary != "periodic" && config.boundary != "fixed" &&
      config.boundary != "periodicAnchored" && config.boundary != "anchored" &&
      !config.usesRectangleBoundary() && config.boundary != "fixedBottom") {
    throw std::invalid_argument("Unknown boundary: " + config.boundary);
  }
  if (config.constitutiveReference != "areaNormalizedTriangular") {
    throw std::invalid_argument("Unsupported constitutiveReference.");
  }
  validateRectangleBoundary(config);
  if (config.usesIndentation()) {
    if (config.boundary != "fixedBottom" || config.rows < 2 || config.cols < 2 ||
        config.loading != "relaxation" ||
        config.segmentLoad != 0 || config.repeat) {
      throw std::invalid_argument("Nanoindentation requires bottom fixation, relaxation, at least "
                                  "two rows and columns, and an unsegmented loading path.");
    }
    const double height = (config.rows - 1) *
        (config.energyFunction == "contiTriangular" ? std::sqrt(3.0) / 2 : 1);
    if (!finiteValue(config.indenterRadius) || config.indenterRadius <= 0 ||
        config.startLoad < 0 || config.loadIncrement <= 0 || config.startLoad + config.totalLoad >= height) {
      throw std::invalid_argument("Nanoindentation requires a positive radius and increment, a nonnegative "
                                  "initial depth, and a final depth below the sample height.");
    }
    if (config.indentationContact != "penalty" && config.indentationContact != "projection") {
      throw std::invalid_argument("indentationContact must be penalty or projection.");
    }
  } else if (config.boundary == "fixedBottom" || config.indenterRadius != 25 ||
             config.indentationContact != "penalty") {
    throw std::invalid_argument("Indenter settings and fixedBottom require nanoindentation.");
  }
  if (config.deformation == "directDisplacement" && config.boundary == "periodic") {
    throw std::invalid_argument("directDisplacement needs driven fixed nodes: use anchored, periodicAnchored, or fixed.");
  }
  const auto directions = directionSequence(config);
  if (config.deformation == "sequentialSlip") {
    if (config.boundary != "anchored" || config.rows < 2 || config.cols < 2) {
      throw std::invalid_argument("sequentialSlip requires an anchored mesh of at least 2 by 2 nodes.");
    }
    if (config.startLoad != 0 || config.segmentLoad != 1 || config.repeat ||
        config.totalLoad != directions.size() || config.loadIncrement <= 0 ||
        config.stopAfterSteps || config.stopAfterPlasticSteps || config.stopAfterElasticSteps ||
        config.stopAfterEnergyDrops || config.stopAfterEnergyIncreases) {
      throw std::invalid_argument("sequentialSlip requires zero startLoad, positive loadIncrement, "
          "one unit per entry, and a full sequence without repetition or early stops.");
    }
  }
  if (!(config.loadIncrement != 0 && std::abs(config.loadIncrement) <= DBL_MAX)) {
    throw std::invalid_argument("loadIncrement must be nonzero and bounded.");
  }
  if (!(finiteValue(config.startLoad) && config.totalLoad >= 0 &&
        config.startLoad <= DBL_MAX - config.totalLoad)) {
    throw std::invalid_argument("startLoad must be finite and totalLoad nonnegative with a bounded sum.");
  }
  if (config.totalLoad > 0 && config.startLoad + config.totalLoad <= config.startLoad) {
    throw std::invalid_argument("totalLoad is too small to advance from startLoad.");
  }
  if (isAxialDeformation(config) && !config.usesRectangleBoundary() &&
      ((config.boundary != "periodic" && config.boundary != "fixed") ||
       directions.size() != 1 || directions.front().sign < 0 || config.repeat ||
       config.segmentLoad != 0 || config.loadIncrement <= 0 || config.startLoad < 0 ||
       axialStretchAtLoad(config, config.startLoad + config.totalLoad) <= 0)) {
    throw std::invalid_argument("Affine compression and tension require periodic or fixed boundaries, "
        "one unsigned direction, and nonnegative engineering strain with positive stretch and increment.");
  }
  if (!(config.segmentLoad >= 0 && config.segmentLoad <= DBL_MAX)) {
    throw std::invalid_argument("segmentLoad must be nonnegative and bounded.");
  }
  if (config.repeat && !(config.segmentLoad * directions.size() <= DBL_MAX)) {
    throw std::invalid_argument("The loading sequence has an invalid period.");
  }
  if ((directions.size() > 1 || config.repeat) && config.segmentLoad == 0) {
    throw std::invalid_argument("Multiple/repeated directions need segmentLoad > 0.");
  }
  if (!(config.stepNoise >= 0 && config.stepNoise <= DBL_MAX))
    throw std::invalid_argument("stepNoise must be nonnegative and bounded.");
  if (!(config.reversibilityTolerance > 0 && config.reversibilityTolerance <= DBL_MAX)) {
    throw std::invalid_argument("reversibilityTolerance must be positive and bounded.");
  }
  if (config.stopWhen != "all" && config.stopWhen != "any") {
    throw std::invalid_argument("stopWhen must be all or any.");
  }
  for (int target : {config.stopAfterSteps, config.stopAfterPlasticSteps,
                     config.stopAfterElasticSteps, config.stopAfterEnergyDrops,
                     config.stopAfterEnergyIncreases}) {
    if (target < 0) throw std::invalid_argument("Step-count targets must be nonnegative.");
  }
  if (!finiteValue(config.stepEnergyTolerance) || config.stepEnergyTolerance < 0) {
    throw std::invalid_argument("stepEnergyTolerance must be finite and nonnegative.");
  }
  if ((config.stopAfterPlasticSteps || config.stopAfterElasticSteps) &&
      (config.loading != "relaxation" || config.deformation != "simpleShear" ||
       config.energyFunction != "contiSquare")) {
    throw std::invalid_argument("Plastic/elastic targets require relaxation, simpleShear and contiSquare.");
  }
  if (config.checkpointInterval < 0 ||
      !(config.energyDropThreshold >= 0 && config.energyDropThreshold <= DBL_MAX)) {
    throw std::invalid_argument("Checkpoint interval and energy-drop threshold must be nonnegative.");
  }
  if (!finiteValue(config.vtuFilesPerUnitLoad) || config.vtuFilesPerUnitLoad < 0 ||
      (config.writeMeshVTUs && config.vtuFilesPerUnitLoad == 0) ||
      !finiteValue(config.vtuFilesPerUnitLoad * config.totalLoad)) {
    throw std::invalid_argument("VTU files per unit load must be finite and positive when mesh output is enabled.");
  }
  if (config.reconnectionMethod != "none" && config.reconnectionMethod != "edgeFlip" &&
      config.reconnectionMethod != "delaunay") {
    throw std::invalid_argument("Unknown reconnection method: " + config.reconnectionMethod);
  }
}

void previewInitialMesh(const Config &config, const std::string &output) {
  validateProtocolConfig(config);
  validateDrivenEdgeConstraint(config);
  omp_set_num_threads(1);
  rng::setSeed(config.seed);
  Simulation simulation;
  simulation.config = config;
  simulation.mesh = makeInitialMesh(config);
  Mesh &mesh = simulation.mesh;
  const auto directions = directionSequence(config);
  prepareSimulation(simulation, directions);
  if (config.usesIndentation()) {
    mesh.initializeIndenter(config.indenterRadius, config.indentationContact == "projection");
  } else {
    initializeProtocolDeformation(simulation);
  }
  std::vector<Eigen::Vector2d> current;
  current.reserve(mesh.nrNodes);
  for (int i = 0; i < mesh.nrNodes; ++i) current.push_back(mesh.nodes(i).pos());
  const Matrix2d currentDeformation = mesh.currentDeformation;

  constexpr double previewIncrement = 0.1;
  const double proposedLoad = config.startLoad + previewIncrement;
  if (config.deformation == "simpleShear" && config.segmentLoad > 0) {
    applyPreviewShear(mesh, config, directions, proposedLoad);
  } else if (!config.usesIndentation() &&
      !(config.deformation == "compression" && !config.usesRectangleBoundary() && proposedLoad >= 1)) {
    const auto scheduled = nextStep(config, directions, config.startLoad, proposedLoad);
    const double sign = config.loadIncrement < 0 ? -1 : 1;
    const double stretch = isAxialDeformation(config) && !config.usesRectangleBoundary()
        ? axialStretchAtLoad(config, proposedLoad) / axialStretchAtLoad(config, config.startLoad) : 1;
    applyDeformationToMesh(mesh, config,
                           {scheduled.direction.axis, sign * scheduled.direction.sign * previewIncrement,
                            scheduled.direction.slip, stretch}, simulation.drivenNodeGroups);
  }
  writeInitialMeshPreview(mesh, current, currentDeformation, output, config.startLoad + config.totalLoad);
}

std::shared_ptr<Simulation> runSimulationProtocol(
    Config config, std::string dataPath, std::shared_ptr<Simulation> simulation,
    double makeDumpAt) {
  validateProtocolConfig(config);
  if (!simulation) validateDrivenEdgeConstraint(config);
  if (config.energyFunction == "contiTriangular" && config.reconnectionMethod != "none") {
    std::cerr << "Warning: T is unreliable for deformation history after triangular reconnection.\n";
  }
  if (config.energyFunction == "conti2Phase") {
    std::cerr << "Scientific issues: Conti2Phase is experimental. Interpret the elastoplastic "
                 "decomposition and deformation history T with caution.\n";
  }
  const Loading loading = loadingChoice(config.loading);
  const auto directions = directionSequence(config);

  if (simulation) {
    // Constraints and reference geometry are part of the saved physical state.
    // A resume restores them; preparation must not be applied a second time.
    validateResumedBoundary(simulation->config, config);
    if (simulation->mesh.usingPBC != config.usesPeriodicBoundary()) {
      throw std::invalid_argument("Cannot change boundary topology on resume.");
    }
    if (loading == Loading::Reversibility) {
      simulation->addReversibilityCsvColumns();
    }
    simulation->termination.initialize(simulation->mesh, config);
    simulation->initialize();
  } else {
    simulation = std::make_shared<Simulation>(config, dataPath, true);
    simulation->makeDumpAt = makeDumpAt;
    prepareSimulation(*simulation, directions);
    if (loading == Loading::Reversibility) {
      simulation->addReversibilityCsvColumns();
    }
    simulation->termination.initialize(simulation->mesh, config);
    if (loading == Loading::NoRelaxation) {
      // Record the starting geometry without adding noise or relaxing it.
      simulation->initialize();
      initializeProtocolDeformation(*simulation);
      simulation->mesh.addLoad(0);
      simulation->finishStep();
    } else {
      // firstStep applies the starting deformation, adds initial noise, and
      // relaxes it using the shared initialization procedure.
      simulation->firstStep();
    }
  }
  Simulation &s = *simulation;
  s.makeDumpAt = makeDumpAt;
  // Callers can extend a run by changing its total budget before resuming.
  s.config.totalLoad = s.totalLoad;

  // This is the main loading loop where the actual simulation takes place.
  // endLoad() is startLoad + totalLoad, so initialization uses none of the budget.
  while (s.mesh.load < s.endLoad() && !s.termination.targetsMet) {
    // Choose a direction and shorten the step at a segment boundary or the end.
    const ScheduledStep scheduled = nextStep(config, directions, s.mesh.load, s.endLoad());
    const double nextLoad = s.mesh.load + scheduled.amount;
    if (!(scheduled.amount > 0 && nextLoad > s.mesh.load)) {
      throw std::runtime_error("Loading step is too small to advance the saved load.");
    }
    if (config.checkpointInterval > 0 && s.mesh.loadSteps % config.checkpointInterval == 0) {
      // Keep alternating checkpoints from before the next deformation step.
      const int slot = (s.mesh.loadSteps / config.checkpointInterval) % 2;
      s.saveSimulation("checkpoint_" + std::to_string(slot));
    }
    s.mesh.ensureForces();
    const double beforeEnergy = s.mesh.totalEnergy;
    const double sign = config.loadIncrement < 0 ? -1 : 1;
    const double stretch = isAxialDeformation(config) && !config.usesRectangleBoundary() ?
        axialStretchAtLoad(config, nextLoad) / axialStretchAtLoad(config, s.mesh.load) : 1;
    const DeformationStep step{scheduled.direction.axis,
                               sign * scheduled.direction.sign * scheduled.amount,
                               scheduled.direction.slip, stretch};
    // Progress is positive; the deformation step contains its physical sign.
    s.loadIncrement = scheduled.amount;
    const auto addStepNoise = [&]() {
      if (config.stepNoise > 0) s.addStepNoise();
    };
    switch (loading) {
    case Loading::Relaxation:
      s.applyDeformation(step);
      if (config.stepNoise > 0) s.addNoiseToGuess(config.stepNoise);
      s.minimize();
      break;
    case Loading::Reversibility:
      // This helper reverses this exact step and restores the forward state.
      // It never asks the schedule to choose a second direction for unloading.
      addStepNoise();
      s.checkReversibility(step, config.reversibilityTolerance);
      break;
    case Loading::NoRelaxation:
      s.applyDeformation(step);
      addStepNoise();
      s.reconnectWithoutMinimization();
      break;
    }
    // Record the planned endpoint, including any shortened final step, and
    // update the measured quantities and output for this accepted state.
    s.mesh.load = nextLoad;
    s.termination.record(config, beforeEnergy, s.mesh.totalEnergy,
                         s.energyHistory.initialGuessTotalEnergy, scheduled.amount);
    if (!s.termination.targetsMet && nextLoad >= s.endLoad()) {
      s.termination.stopReason = "strain_limit";
    }
    s.finishStep();
  }
  s.loadIncrement = std::abs(config.loadIncrement);
  if (!s.termination.targetsMet) s.termination.stopReason = "strain_limit";
  s.finishSimulation();
  return simulation;
}
