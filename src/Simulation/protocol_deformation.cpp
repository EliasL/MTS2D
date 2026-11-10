#include "protocol_deformation.h"
#include <algorithm>
#include <cfloat>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <stdexcept>

namespace {
SignedDirection slipDirection(const Config &config, const std::string &token,
                              double sign) {
  if (token.empty() || (token[0] != '1' && token[0] != '2') ||
      (token.size() > 1 &&
       (token[1] != '.' || token.size() < 3 || token.size() > 11 ||
        token.find_first_not_of("0123456789", 2) != std::string::npos))) {
    throw std::invalid_argument("Invalid slip entry: " + token);
  }
  const bool horizontal = token[0] == '1';
  const int size = horizontal ? config.rows : config.cols;
  int index = size / 2;
  if (token.size() > 1) {
    // Convert the fraction exactly so a value such as .29 selects index 29
    // on a 100-node edge, independently of floating-point rounding.
    std::int64_t denominator = 1;
    for (size_t digit = 2; digit < token.size(); ++digit) denominator *= 10;
    index = size * std::stoll(token.substr(2)) / denominator;
  }
  if (index < 1 || index >= size) {
    throw std::invalid_argument("Slip entry must select an interior split on the edge: " + token);
  }
  return {horizontal ? Direction::Horizontal : Direction::Vertical, sign, {index}};
}

Matrix2d shearTransform(Direction axis, double amount) {
  Matrix2d transform = getShear(amount);
  if (axis == Direction::Vertical) transform.transposeInPlace();
  return transform;
}

Matrix2d axialTransform(Direction axis, double stretch) {
  Matrix2d transform = Matrix2d::Identity();
  const int index = axis == Direction::Horizontal ? 0 : 1;
  transform(index, index) = stretch;
  return transform;
}

Matrix2d deformationIncrement(const Config &config, SignedDirection direction, double amount) {
  const double sign = config.loadIncrement < 0 ? -1 : 1;
  return shearTransform(direction.axis, sign * direction.sign * amount);
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

} // namespace

bool isAxialDeformation(const Config &config) {
  return config.deformation == "compression" || config.deformation == "tension";
}

double axialStretchAtLoad(const Config &config, double load) {
  return 1 + (config.deformation == "compression" ? -load : load);
}

std::vector<SignedDirection> directionSequence(const Config &config) {
  std::string text = config.directions;
  text.erase(std::remove_if(text.begin(), text.end(),
      [](unsigned char c) { return std::isspace(c); }), text.end());
  if (text.empty() || text.back() == ',') {
    throw std::invalid_argument("directions requires a nonempty comma-separated sequence.");
  }
  std::vector<SignedDirection> directions;
  std::array<double, 2> slipSigns{};
  std::istringstream input(text);
  std::string name;
  while (std::getline(input, name, ',')) {
    double sign = 1;
    if (!name.empty() && (name.front() == '-' || name.front() == '+')) {
      sign = name.front() == '-' ? -1 : 1;
      name.erase(name.begin());
    }
    if (config.deformation == "sequentialSlip") {
      const auto direction = slipDirection(config, name, sign);
      double &systemSign = slipSigns[direction.axis == Direction::Horizontal ? 0 : 1];
      if (systemSign != 0 && systemSign != sign) {
        throw std::invalid_argument("Each slip system must use a single sign.");
      }
      systemSign = sign;
      directions.push_back(direction);
    } else if (name == "horizontal") {
      directions.push_back({Direction::Horizontal, sign});
    } else if (name == "vertical") {
      directions.push_back({Direction::Vertical, sign});
    } else {
      throw std::invalid_argument("Unknown direction: " + name);
    }
  }
  for (auto &direction : directions) {
    direction.slip.lastRow = slipSigns[1] < 0;
    direction.slip.lastColumn = slipSigns[0] < 0;
  }
  return directions;
}

void prepareSimulation(Simulation &simulation, const std::vector<SignedDirection> &directions) {
  // Prepare constraints before creating the solver's free-node arrays.
  // A resumed simulation already carries these constraints in its saved mesh.
  Mesh &mesh = simulation.mesh;
  if (simulation.config.boundary == "fixed") {
    mesh.fixBorderNodes();
  }
  if (simulation.config.boundary == "anchored" ||
      simulation.config.boundary == "periodicAnchored") {
    const SlipRegion &slip = directions.front().slip;
    mesh.fixNodesInRow(slip.lastRow ? mesh.rows - 1 : 0);
    mesh.fixNodesInColumn(slip.lastColumn ? mesh.cols - 1 : 0);
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
}

Matrix2d protocolDeformation(const Config &config, double load) {
  auto directions = directionSequence(config);
  if (isAxialDeformation(config)) {
    return axialTransform(directions.front().axis, axialStretchAtLoad(config, load));
  }
  // Combine the shear segments up to this path coordinate. Horizontal and
  // vertical shears must be multiplied in their original loading order.
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
  if (config.deformation == "simpleShear" || isAxialDeformation(config)) {
    const Matrix2d transform = isAxialDeformation(config) ? axialTransform(step.direction, step.axialStretch)
                                                        : shearTransform(step.direction, step.amount);
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
  } else if (config.deformation == "sequentialSlip") {
    const bool horizontal = step.direction == Direction::Horizontal;
    if (step.slip.index < 1 || step.slip.index >= (horizontal ? mesh.rows : mesh.cols)) {
      throw std::invalid_argument("Missing or invalid sequential slip index.");
    }
    const int edge = horizontal ? (step.slip.lastColumn ? mesh.cols - 1 : 0)
                                : (step.slip.lastRow ? mesh.rows - 1 : 0);
    // Translate the side of the split away from the other anchored edge.
    const bool beforeSplit = horizontal ? step.slip.lastRow : step.slip.lastColumn;
    const Vector2d displacement = horizontal ? Vector2d(step.amount, 0) : Vector2d(0, step.amount);
    bool moved = false;
    for (const NodeId &id : mesh.fixedNodeIds) {
      const int alongEdge = horizontal ? id.row() : id.col();
      if ((horizontal ? id.col() : id.row()) == edge &&
          (alongEdge < step.slip.index) == beforeSplit) {
        mesh[id]->addDisplacement(displacement);
        moved = true;
      }
    }
    if (!moved) throw std::invalid_argument("The slip region contains no constrained nodes.");
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
  if (config.deformation == "simpleShear" || isAxialDeformation(config)) {
    // The starting transformation acts on the whole mesh, including free nodes.
    // The caller then relaxes this starting state if the loading mode asks for it.
    simulation.mesh.applyTransformation(protocolDeformation(config, config.startLoad));
    return;
  }
  if (config.deformation == "sequentialSlip") {
    // A fresh slip sequence starts from the unloaded geometry.
    return;
  }
  // Prescribed translations add, including on overlapping slip regions.
  const auto directions = directionSequence(config);
  std::array<double, 2> amounts{};
  const double sign = (config.loadIncrement < 0 ? -1 : 1) * (config.startLoad < 0 ? -1 : 1);
  auto collect = [&](SignedDirection direction, double amount) {
    amounts[direction.axis == Direction::Horizontal ? 0 : 1] += sign * direction.sign * amount;
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
  for (size_t i = 0; i < amounts.size(); ++i) {
    applyDeformationToMesh(simulation.mesh, config,
        {i == 0 ? Direction::Horizontal : Direction::Vertical, amounts[i]},
        simulation.drivenNodeGroups);
  }
}
