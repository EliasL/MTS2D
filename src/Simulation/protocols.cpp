/*
A protocol is the whole experiment described by the configuration.

  loading      says how we examine each step: relax, test reversibility,
               or apply the deformation without relaxation.
  deformation  says what displacement pattern we impose: simpleShear or
               directDisplacement.
  boundary     says which nodes or periodic images are driven.
  directions   says which way we shear, and in which order.
  reconnectionMethod is handled by Simulation during relaxation/reconnection.

The runner below follows this pseudocode:

  validate the choices
  prepare a new simulation, or resume its saved state
  while there is more of the loading path to follow:
      choose the next direction and step size
      shorten the step if it would cross a segment boundary or exhaust totalLoad
      apply the selected loading behavior
      record the state and check whether to stop
  finish the output

To add a deformation, extend applyDeformationToMesh() and initialization.
To add a loading behavior, extend the switch in runSimulationProtocol().
Boundary preparation belongs in prepareSimulation(). For direct displacement,
that is also where we choose which fixed nodes move in each direction.
All these choices share one simulation loop.

For example, loading=relaxation, deformation=simpleShear, boundary=periodic,
directions=horizontal,vertical, segmentLoad=0.2, repeat=true means:

  relax after each small shear step
  switch between horizontal and vertical shear every 0.2 of path progress
  repeat until totalLoad has been accumulated
  use the chosen reconnection method throughout

totalLoad is the loading budget after initialization. The initial deformation
at startLoad does not consume it. The saved mesh.load is the path coordinate:
accumulated loading is mesh.load - startLoad. It always increases, even when a
deformation step is negative, so backward/forward sequences still finish.
Shear amounts are strains; direct-displacement amounts are distances.
A minus sign on a direction reverses that segment; a negative loadIncrement
reverses every segment. Neither changes how progress is counted.

segmentLoad is the amount of path traversed per direction. Zero means one
unswitched direction. A finite sequence continues in its final direction;
repeat=true returns to its first direction instead. The schedule is determined
from the saved load, so it needs no separate counter to resume correctly.
*/
#include "protocols.h"
#include "finite.h"
#include <algorithm>
#include <cfloat>
#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {
enum class Loading { Relaxation, Reversibility, NoRelaxation };
struct SignedDirection {
  Direction axis;
  double sign;
};

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

std::vector<SignedDirection> directionSequence(const std::string &text) {
  if (text.empty() || text.back() == ',') {
    throw std::invalid_argument("directions must contain horizontal or vertical.");
  }
  std::vector<SignedDirection> directions;
  std::istringstream input(text);
  std::string name;
  while (std::getline(input, name, ',')) {
    name.erase(std::remove_if(name.begin(), name.end(),
        [](unsigned char c) { return std::isspace(c); }), name.end());
    double sign = 1;
    if (!name.empty() && (name.front() == '-' || name.front() == '+')) {
      sign = name.front() == '-' ? -1 : 1;
      name.erase(name.begin());
    }
    if (name == "horizontal") {
      directions.push_back({Direction::Horizontal, sign});
    } else if (name == "vertical") {
      directions.push_back({Direction::Vertical, sign});
    } else {
      throw std::invalid_argument("Unknown direction: " + name);
    }
  }
  return directions;
}

Matrix2d deformationIncrement(const Config &config, SignedDirection direction,
                               double amount) {
  // Horizontal shear moves in x with a gradient in y. Vertical shear is its
  // transpose. The schedule and the sign of loadIncrement choose the sense.
  const double sign = config.loadIncrement < 0 ? -1 : 1;
  Matrix2d step = getShear(sign * direction.sign * amount);
  if (direction.axis == Direction::Vertical) {
    step.transposeInPlace();
  }
  return step;
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

Matrix2d repeatedTransform(Matrix2d transform, double count) {
  // Apply whole repeated cycles without walking every small loading increment.
  Matrix2d result = Matrix2d::Identity();
  while (count >= 1) {
    if (std::fmod(count, 2) >= 1) {
      result = (transform * result).eval();
    }
    count = std::floor(count / 2);
    if (count >= 1) {
      transform = (transform * transform).eval();
    }
  }
  return result;
}

void prepareSimulation(Simulation &simulation) {
  // Prepare constraints before creating the solver's free-node arrays.
  // A resumed simulation already carries these constraints in its saved mesh.
  Mesh &mesh = simulation.mesh;
  if (simulation.config.boundary == "fixed") {
    mesh.fixBorderNodes();
  }
  if (simulation.config.boundary == "anchored" ||
      simulation.config.boundary == "periodicAnchored") {
    mesh.fixNodesInRow(0);
    mesh.fixNodesInColumn(0);
  }
  if (simulation.config.deformation == "directDisplacement") {
    // Select by original grid indices, once. These same IDs are saved in dumps
    // and reused for forward steps, reverse checks, and crash replay.
    for (const NodeId &id : mesh.fixedNodeIds) {
      if (id.row() >= mesh.rows / 2) {
        simulation.drivenNodeGroups[0].push_back(id);
      }
      if (id.col() >= mesh.cols / 2) {
        simulation.drivenNodeGroups[1].push_back(id);
      }
    }
  }
  if (simulation.config.loading == "reversibility") {
    simulation.addReversibilityCsvColumns();
  }
}
} // namespace

void validateProtocolConfig(const Config &config) {
  const Loading loading = loadingChoice(config.loading);
  if (config.deformation != "simpleShear" && config.deformation != "directDisplacement") {
    throw std::invalid_argument("Unknown deformation: " + config.deformation);
  }
  if (config.boundary != "periodic" && config.boundary != "fixed" &&
      config.boundary != "periodicAnchored" && config.boundary != "anchored") {
    throw std::invalid_argument("boundary must be periodic, periodicAnchored, anchored, or fixed.");
  }
  if (config.deformation == "directDisplacement" && config.boundary == "periodic") {
    throw std::invalid_argument("directDisplacement needs driven fixed nodes: use anchored, periodicAnchored, or fixed.");
  }
  const auto directions = directionSequence(config.directions);
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
  if (config.reconnectionMethod != "none" && config.reconnectionMethod != "edgeFlip" &&
      config.reconnectionMethod != "delaunay") {
    throw std::invalid_argument("Unknown reconnection method: " + config.reconnectionMethod);
  }
}

Matrix2d protocolDeformation(const Config &config, double load) {
  // Combine the shear segments up to this path coordinate. Horizontal and
  // vertical shears must be multiplied in their original loading order.
  auto directions = directionSequence(config.directions);
  if (load < 0) {
    for (auto &direction : directions) direction.sign *= -1;
    load = -load;
  }
  if (config.segmentLoad == 0) {
    return deformationIncrement(config, directions.front(), load);
  }
  Matrix2d result = Matrix2d::Identity();
  double remaining = load;
  if (config.repeat) {
    const double period = config.segmentLoad * directions.size();
    if (!(period > 0 && period <= DBL_MAX)) {
      throw std::invalid_argument("The loading sequence has an invalid period.");
    }
    Matrix2d cycle = Matrix2d::Identity();
    for (const SignedDirection direction : directions) {
      cycle = (deformationIncrement(config, direction, config.segmentLoad) * cycle).eval();
    }
    result = repeatedTransform(cycle, std::floor(load / period));
    remaining = std::fmod(load, period);
  }
  for (size_t i = 0; i < directions.size() && remaining > 0; ++i) {
    const bool holdLast = !config.repeat && i + 1 == directions.size();
    const double amount = holdLast ? remaining : std::min(remaining, config.segmentLoad);
    result = (deformationIncrement(config, directions[i], amount) * result).eval();
    remaining -= amount;
  }
  return result;
}

void applyDeformationToMesh(Mesh &mesh, const Config &config,
                            const DeformationStep &step,
                            const DrivenNodeGroups &drivenNodes) {
  // Apply one imposed deformation. The step already contains its direction
  // and signed amount; this function does not choose the next segment.
  // Loading counters, relaxation, and output are handled by the caller. This
  // lets reversibility snapshots use the same deformation as the actual run.
  if (config.deformation == "simpleShear") {
    // Horizontal shear moves in x with a gradient in y. Transposing the matrix
    // gives vertical shear: movement in y with a gradient in x.
    Matrix2d transform = getShear(step.amount);
    if (step.direction == Direction::Vertical) {
      transform.transposeInPlace();
    }
    // For periodic boundaries we deform the cell. For open boundaries we move
    // the fixed nodes that impose the deformation on the free interior.
    if (mesh.usingPBC) {
      mesh.applyTransformationToSystemDeformation(transform);
    } else {
      mesh.applyTransformationToFixedNodes(transform);
    }
    // Free nodes can start the next relaxation at affinely predicted positions.
    // With this option disabled, their previous positions remain the guess.
    if (config.affineInitialGuess) {
      const Matrix2d change = transform - Matrix2d::Identity();
      for (const NodeId &id : mesh.freeNodeIds) {
        Node &node = *mesh[id];
        // Write the new displacement as (transform - I) X + transform * u.
        node.setDisplacement(change * node.ref_pos() + transform * node.u());
      }
    }
  } else if (config.deformation == "directDisplacement") {
    // Every node in the selected group receives the same translation. The
    // groups were chosen during preparation: upper-half fixed nodes move
    // horizontally, and right-half fixed nodes move vertically.
    const bool horizontal = step.direction == Direction::Horizontal;
    const auto &nodes = drivenNodes[horizontal ? 0 : 1];
    if (nodes.empty()) {
      throw std::invalid_argument("No driven nodes prepared for directDisplacement.");
    }
    const Vector2d displacement = horizontal ? Vector2d(step.amount, 0)
                                             : Vector2d(0, step.amount);
    for (const NodeId &id : nodes) {
      mesh[id]->addDisplacement(displacement);
    }
    // Free nodes stay at their previous positions until relaxation moves them.
  } else {
    throw std::invalid_argument("Unknown deformation: " + config.deformation);
  }
  // Element geometry and forces must be refreshed before they are used again.
  mesh.markDirty();
}

void Simulation::applyDeformation(const DeformationStep &step,
                                  bool saveReplayCheckpoint) {
  // Save the state before changing the load or positions, so a failed
  // minimization can retry the exact same step.
  if (saveReplayCheckpoint && !debugReplayActive) {
    saveLoadingStepReplayCheckpoint(step);
  }
  mesh.addLoad(loadIncrement);
  applyDeformationToMesh(mesh, config, step, drivenNodeGroups);
  // Copy the resulting free-node displacements into both solver arrays.
  applyLoadStepToGuess();
}

void initializeProtocolDeformation(Simulation &simulation) {
  // Build the starting geometry for a new simulation at startLoad. Earlier
  // segments are combined before the first relaxation; we do not run their
  // intermediate relaxation or reconnection steps here.
  // This does not advance the load counters or consume any of totalLoad.
  // Resumed simulations already have their geometry and skip this function.
  const Config &config = simulation.config;
  if (config.deformation == "simpleShear") {
    // The combined shear acts on the whole initial mesh, including free nodes.
    // The caller then relaxes this starting state if the loading mode asks for it.
    simulation.mesh.applyTransformation(protocolDeformation(config, config.startLoad));
    return;
  }
  // Translations of the two fixed groups add, so whole initial cycles can be
  // collected without simulating every increment or performing relaxation.
  const auto directions = directionSequence(config.directions);
  // Collect horizontal and vertical displacement separately. A minus sign in
  // a direction, or a negative loadIncrement, reverses the corresponding amount.
  std::array<double, 2> amounts{0, 0};
  const double sign = (config.loadIncrement < 0 ? -1 : 1) * (config.startLoad < 0 ? -1 : 1);
  auto collect = [&](SignedDirection direction, double amount) {
    amounts[direction.axis == Direction::Horizontal ? 0 : 1] +=
        sign * direction.sign * amount;
  };
  double remaining = std::abs(config.startLoad);
  if (config.segmentLoad == 0) {
    collect(directions.front(), remaining);
    remaining = 0;
  } else if (config.repeat) {
    // Count all complete cycles first, then handle the remaining partial cycle.
    const double period = config.segmentLoad * directions.size();
    const double cycles = std::floor(remaining / period);
    for (const auto &direction : directions) {
      collect(direction, cycles * config.segmentLoad);
    }
    remaining = std::fmod(remaining, period);
  }
  for (size_t i = 0; i < directions.size() && remaining > 0; ++i) {
    // A non-repeating sequence stays in its final direction after earlier
    // segments have finished, just as it does in the main loading loop.
    const bool holdLast = !config.repeat && i + 1 == directions.size();
    const double amount = holdLast ? remaining : std::min(remaining, config.segmentLoad);
    collect(directions[i], amount);
    remaining -= amount;
  }
  // Apply the two collected translations to the same saved node groups that
  // will be driven during loading. Their translations add even if groups overlap.
  for (size_t i = 0; i < amounts.size(); ++i) {
    applyDeformationToMesh(simulation.mesh, config,
        {i == 0 ? Direction::Horizontal : Direction::Vertical, amounts[i]},
        simulation.drivenNodeGroups);
  }
}

std::shared_ptr<Simulation> runSimulationProtocol(
    Config config, std::string dataPath, std::shared_ptr<Simulation> simulation,
    double makeDumpAt) {
  validateProtocolConfig(config);
  const Loading loading = loadingChoice(config.loading);
  const auto directions = directionSequence(config.directions);

  if (simulation) {
    // Constraints and reference geometry are part of the saved physical state.
    // A resume restores them; preparation must not be applied a second time.
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
    prepareSimulation(*simulation);
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
    const DeformationStep step{scheduled.direction.axis,
                               sign * scheduled.direction.sign * scheduled.amount};
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
