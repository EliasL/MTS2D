#include "logging.h"
#include "Mesh/mesh.h"
#include "Mesh/node.h"
#include "Mesh/tElement.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>

Timer::Timer() {}

void Timer::addKey(const std::string &key) {
  if (runtimes.find(key) == runtimes.end()) {
    runtimes[key] = std::chrono::milliseconds(0);
    runningStatus[key] = false;
    checkpoints[key] = std::chrono::high_resolution_clock::now();
  }
}

void Timer::Start(const std::string &key) {
  // Ensure that each key is initialized properly if it doesn't exist
  if (runtimes.find(key) == runtimes.end()) {
    runtimes[key] = std::chrono::milliseconds(0);
    runningStatus[key] =
        false; // Ensures that the timer is initially not running
  }

  // Start the timer only if it's not already running
  if (!runningStatus[key]) {
    checkpoints[key] = std::chrono::high_resolution_clock::now();
    runningStatus[key] = true;
  }
}

// Adds the time from the checkpoint to the runtime
// This prevents the difference between now and the checkpoint from
// getting very large. This means that if the program crashes, we don't loose
// much progress. (Since when we restart and calculate the difference between
// the checkpoint and now, it will be huge, so we don't want that)
void Timer::Save(const std::string &key) {
  auto now = std::chrono::high_resolution_clock::now();
  auto timeInterval = std::chrono::duration_cast<std::chrono::milliseconds>(
      now - checkpoints[key]);
  runtimes[key] += timeInterval;
  checkpoints[key] = now;
}

void Timer::SaveAll() {
  for (const auto &pair : runtimes) {
    Save(pair.first);
  }
}

size_t Timer::Stop(const std::string &key) {
  if (runtimes.find(key) == runtimes.end()) {
    std::cerr << "Timer key " << key << " not found!\n";
  }

  if (runningStatus[key]) {
    auto now = std::chrono::high_resolution_clock::now();
    auto timeInterval = std::chrono::duration_cast<std::chrono::milliseconds>(
        now - checkpoints[key]);
    runtimes[key] += timeInterval;
    runningStatus[key] = false;
    return timeInterval.count();
  }
  return 0;
}

std::chrono::milliseconds Timer::RunTime(const std::string &key) const {
  // Check that the key exists in all required maps.
  if (runningStatus.find(key) == runningStatus.end() ||
      checkpoints.find(key) == checkpoints.end() ||
      runtimes.find(key) == runtimes.end()) {
    // throw std::invalid_argument("Key '" + key + "' not found in Timer
    // maps.");
    std::cerr << "Warning: Key '" + key + "' not found in Timer maps.\n";
    return std::chrono::milliseconds::zero();
  }

  if (runningStatus.at(key)) {
    auto now = std::chrono::high_resolution_clock::now();
    auto currentDuration =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            now - checkpoints.at(key));
    return (runtimes.at(key) + currentDuration);
  } else {
    return runtimes.at(key);
  }
}

std::chrono::milliseconds Timer::ETR(double completion) {
  // Completion is expected to be in [0, 1]
  assert(completion >= 0.0 && completion <= 1.0);

  // Add the current runtime and completion to the deques
  average_time.push_back(RunTime());
  average_completion.push_back(completion);

  // Keep the deque size fixed at 1024 by removing the oldest entry if necessary
  // At 20 seconds per update, that makes for an average period of 5.7 hours
  if (average_time.size() > 1024) {
    average_time.pop_front();
    average_completion.pop_front();
  }

  // Need at least three data points so we can ignore the first sample
  // (the first time step is unnaturally long)
  if (average_time.size() < 3) {
    return std::chrono::milliseconds(0);
  }

  // Use the 2nd sample as the baseline, and the last sample as the latest
  const size_t base_idx = 1;
  const auto t0 = average_time[base_idx];
  const auto t1 = average_time.back();
  const double c0 = average_completion[base_idx];
  const double c1 = average_completion.back();

  const double total_time_change_ms =
      std::chrono::duration<double, std::milli>(t1 - t0).count();
  const double total_completion_change = c1 - c0;

  // Guard against non-positive changes (no progress or non-monotone signals)
  if (!(total_time_change_ms > 0.0) || !(total_completion_change > 0.0)) {
    return std::chrono::milliseconds(0);
  }

  // Calculate the average rate of completion per millisecond
  const double average_rate = total_completion_change / total_time_change_ms;
  if (!(average_rate > 0.0)) {
    return std::chrono::milliseconds(0);
  }

  const double remaining_completion = 1.0 - completion;
  if (!(remaining_completion > 0.0)) {
    return std::chrono::milliseconds(0);
  }

  // Estimate the time remaining in milliseconds
  const double estimated_time_remaining_ms =
      remaining_completion / average_rate;

  // Cap / sanity checks
  if (!(estimated_time_remaining_ms >= 0.0)) {
    return std::chrono::milliseconds(0);
  }

  const double max_ms =
      static_cast<double>(std::chrono::milliseconds::max().count());
  if (!(estimated_time_remaining_ms < max_ms)) {
    return std::chrono::milliseconds::max();
  }

  return std::chrono::milliseconds(
      static_cast<int64_t>(estimated_time_remaining_ms));
}

std::string Timer::RTString(const std::string &key, int precision) const {
  // Check if key is present
  return FormatDuration(RunTime(key), precision);
}

std::string Timer::ETRString(double progress) {
  oldETRString = FormatDuration(ETR(progress));
  return oldETRString;
}

void Timer::Reset(const std::string &key) {
  runtimes[key] = std::chrono::milliseconds(0);
  runningStatus[key] = false;
}

// Implementation of PrintAllRuntimes function
void Timer::PrintAllRuntimes() const {
  // Collect all keys and their durations in a vector of pairs
  std::vector<std::pair<std::string, std::chrono::milliseconds>>
      keyDurationPairs;
  for (const auto &pair : runtimes) {
    // We use RunTime instead of pair.second so that we can get updated values
    // even if the timers are still runing.
    keyDurationPairs.emplace_back(pair.first, RunTime(pair.first));
  }

  // Sort the vector by duration
  std::sort(keyDurationPairs.begin(), keyDurationPairs.end(),
            [](const std::pair<std::string, std::chrono::milliseconds> &a,
               const std::pair<std::string, std::chrono::milliseconds> &b) {
              return a.second > b.second;
            });

  // Find the default timer's runtime for percentage calculation
  std::chrono::milliseconds defaultRuntime = RunTime(DEFAULT_KEY);
  double defaultRuntimeInMs = defaultRuntime.count();

  // Print the formatted duration for each key and their percentage compared to
  // the default
  for (const auto &pair : keyDurationPairs) {
    double durationInMs = pair.second.count();
    double percentage = 100.0 * (durationInMs / defaultRuntimeInMs);
    std::cout << pair.first << ": " << FormatDuration(pair.second) << " ("
              << std::fixed << std::setprecision(2) << percentage << "%)"
              << std::endl;
  }
}

std::string FormatDuration(std::chrono::milliseconds duration, int precision) {

  bool useMilliseconds = duration.count() < 1e4 || precision > 3;
  std::ostringstream stream;
  auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
  duration -= hours;
  auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
  duration -= minutes;
  auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
  auto milliseconds = duration.count() % 1000;

  if (hours.count() >= 24) {
    auto days = hours.count() / 24;
    hours = std::chrono::hours(hours.count() % 24);
    stream << days << "d ";
  }

  if (hours.count() > 0) {
    stream << hours.count() << "h ";
  }

  if (minutes.count() > 0) {
    stream << minutes.count() << "m ";
  }
  if (useMilliseconds) {
    stream << std::fixed << std::setprecision(3)
           << seconds.count() + milliseconds / 1000.0 << "s";
  } else {
    stream << seconds.count() << "s";
  }
  return stream.str();
}

// Function to calculate the Estimated Time Remaining (ETR) using progress
// fraction
std::chrono::milliseconds calculateETR(std::chrono::milliseconds elapsed,
                                       float progressFraction) {
  if (progressFraction <= 0) {
    return std::chrono::milliseconds(
        0); // Avoid division by zero if no progress
  }
  double elapsedSeconds = elapsed.count() / 1000.0;
  if (!(elapsedSeconds > 0.0)) {
    return std::chrono::milliseconds(0);
  }
  double rate = progressFraction / elapsedSeconds;
  if (!(rate > 0.0)) {
    return std::chrono::milliseconds(0);
  }
  long long etrInMilliseconds =
      static_cast<long long>(((1 - progressFraction) / rate) * 1000);

  // Ignore negative values
  if (etrInMilliseconds < 0) {
    etrInMilliseconds = 0;
  }
  return std::chrono::milliseconds(etrInMilliseconds);
}

namespace {

std::string formatVector(const Eigen::Vector2d &v) {
  std::ostringstream oss;
  oss << std::setprecision(std::numeric_limits<double>::max_digits10) << v;
  return oss.str();
}

std::string formatVector(const Eigen::Vector2i &v) {
  std::ostringstream oss;
  oss << "(" << v[0] << ", " << v[1] << ")";
  return oss.str();
}

std::string updateStateName(Mesh::UpdateState state) {
  switch (state) {
  case Mesh::UpdateState::Dirty:
    return "Dirty";
  case Mesh::UpdateState::Forces:
    return "Forces";
  case Mesh::UpdateState::Geometry:
    return "Geometry";
  case Mesh::UpdateState::Full:
    return "Full";
  }
  return "Unknown";
}

} // namespace

namespace DebugLog {

GhostNodeSnapshot snapshot(const GhostNode &node) {
  GhostNodeSnapshot s;
  s.referenceId = node.referenceId.i;
  s.id = node.id;
  s.periodicShift = node.periodicShift;
  s.pos = node.pos;
  s.refPos = node.ref_pos;
  return s;
}

ElementSnapshot snapshot(const TElement &element) {
  ElementSnapshot s;
  s.eIndex = element.eIndex;
  s.m3Nr = element.m3Nr;
  s.pastM3Nr = element.pastM3Nr;
  s.pastStepM3Nr = element.pastStepM3Nr;
  s.redQuadrant = element.red_quadrant;
  s.angleNode = element.angleNode;
  if (element.angleNode >= 0 &&
      element.angleNode < static_cast<int>(element.ghostNodes.size())) {
    const GhostNode &first =
        element.ghostNodes[(element.angleNode + 1) % 3];
    const GhostNode &second =
        element.ghostNodes[(element.angleNode + 2) % 3];
    s.angleEdgeNodeIdA = std::min(first.referenceId.i, second.referenceId.i);
    s.angleEdgeNodeIdB = std::max(first.referenceId.i, second.referenceId.i);
  }
  s.noise = element.noise;
  for (size_t i = 0; i < element.ghostNodes.size(); ++i) {
    s.ghostNodes[i] = snapshot(element.ghostNodes[i]);
  }
  return s;
}

std::string formatGhostNodeSnapshot(const GhostNodeSnapshot &snapshot,
                                    std::string_view name) {
  std::ostringstream oss;
  oss << name << ": referenceId=" << snapshot.referenceId
      << ", id=" << formatVector(snapshot.id)
      << ", periodicShift=" << formatVector(snapshot.periodicShift)
      << ", pos=" << formatVector(snapshot.pos)
      << ", ref_pos=" << formatVector(snapshot.refPos);
  return oss.str();
}

std::string formatElementSnapshot(const ElementSnapshot &snapshot,
                                  std::string_view name) {
  std::ostringstream oss;
  oss << name << ":\n"
      << "  eIndex=" << snapshot.eIndex << ", m3Nr=" << snapshot.m3Nr
      << ", pastM3Nr=" << snapshot.pastM3Nr
      << ", pastStepM3Nr=" << snapshot.pastStepM3Nr
      << ", redQuadrant=" << snapshot.redQuadrant
      << ", angleNode=" << snapshot.angleNode << ", angleEdge=";
  if (snapshot.angleEdgeNodeIdA >= 0 && snapshot.angleEdgeNodeIdB >= 0) {
    oss << formatEdgeKey(snapshot.angleEdgeNodeIdA, snapshot.angleEdgeNodeIdB);
  } else {
    oss << "invalid";
  }
  oss << ", noise=" << std::setprecision(std::numeric_limits<double>::max_digits10)
      << snapshot.noise << "\n";
  for (size_t i = 0; i < snapshot.ghostNodes.size(); ++i) {
    oss << "  " << formatGhostNodeSnapshot(snapshot.ghostNodes[i],
                                           "ghost[" + std::to_string(i) + "]")
        << "\n";
  }
  return oss.str();
}

std::string formatElementState(const TElement &element, std::string_view name) {
  std::ostringstream oss;
  oss << formatElementSnapshot(snapshot(element), name);
  oss << "  currentArea="
      << std::setprecision(std::numeric_limits<double>::max_digits10)
      << element.area()
      << ", referenceArea=" << tElementInitialArea(element.ghostNodes) << "\n"
      << "  F:\n"
      << element.F << "\n"
      << "  C:\n"
      << element.C << "\n"
      << "  C_R:\n"
      << element.C_R << "\n"
      << "  G:\n"
      << element.G << "\n"
      << "  M_l:\n"
      << element.M_l << "\n";
  return oss.str();
}

std::string formatInvalidAngleNode(std::string_view context,
                                   const TElement &element) {
  std::ostringstream oss;
  oss << context << ": invalid angleNode=" << element.angleNode
      << ". Expected an index in [0, 2].\n"
      << "This usually means the element geometry was marked stale or modified "
         "after construction, but updateGeometry()/updateAngleNode() was not "
         "run before angle-node dependent access.\n"
      << formatElementState(element, "Element");
  return oss.str();
}

std::string formatEdgeKey(int nodeIdA, int nodeIdB) {
  std::ostringstream oss;
  oss << "(node " << nodeIdA << ", node " << nodeIdB << ")";
  return oss.str();
}

std::string formatEdgeTwinLookupOverflow(int nodeIdA, int nodeIdB,
                                         const TElement &firstElement,
                                         const TElement &secondElement,
                                         const TElement &newElement) {
  std::ostringstream oss;
  oss << "addToEdgeTwinLookup: more than two elements share angle edge "
      << formatEdgeKey(nodeIdA, nodeIdB) << ".\n"
      << "Existing element indices: " << firstElement.eIndex << ", "
      << secondElement.eIndex << "; new element index: " << newElement.eIndex
      << ".\n"
      << "Situation: this lookup key is built from the real node ids stored in "
         "ghostNode.referenceId.i. If these elements are in different periodic "
         "images, the same real-node edge may appear in multiple ghost "
         "placements. If they are from the same placement, the lookup may have "
         "gone stale after an edge flip.\n"
      << formatElementState(firstElement, "Existing first element")
      << formatElementState(secondElement, "Existing second element")
      << formatElementState(newElement, "New element");
  return oss.str();
}

std::string formatEdgeTwinLookupMismatch(
    int nodeIdA, int nodeIdB, const TElement &queryElement,
    const TElement *firstElement, int firstElementIndex,
    const TElement *secondElement, int secondElementIndex) {
  std::ostringstream oss;
  oss << "findTwinFromLookup: edge lookup mismatch.\n"
      << "  Query element index: " << queryElement.eIndex << "\n"
      << "  Query angle edge node ids: " << formatEdgeKey(nodeIdA, nodeIdB)
      << "\n"
      << "  Lookup entry element indices: " << firstElementIndex << ", "
      << secondElementIndex << "\n"
      << "  Situation: the queried element has an angle edge that exists in "
         "the lookup, but the lookup entry points to two different elements. "
         "This edge key is not a simple twin pair for the queried element. In "
         "PBC cases this may indicate that the same real-node edge appears in "
         "multiple ghost placements, or that the lookup was not updated "
         "consistently after an edge flip.\n"
      << formatElementState(queryElement, "Query element");

  if (firstElement != nullptr) {
    oss << formatElementState(*firstElement, "Lookup first element");
  }
  if (secondElement != nullptr) {
    oss << formatElementState(*secondElement, "Lookup second element");
  }
  return oss.str();
}

std::string formatMeshContext(const Mesh &mesh) {
  std::ostringstream oss;
  oss << "Mesh context:\n"
      << "  load=" << std::setprecision(std::numeric_limits<double>::max_digits10)
      << mesh.load << ", loadSteps=" << mesh.loadSteps
      << ", nrMinFunctionCalls=" << mesh.nrMinFunctionCalls
      << ", nrMinItterations=" << mesh.nrMinItterations << "\n"
      << "  totalEnergy=" << mesh.totalEnergy << ", maxForce=" << mesh.maxForce
      << ", rows=" << mesh.rows << ", cols=" << mesh.cols
      << ", nrElements=" << mesh.nrElements << ", nrNodes=" << mesh.nrNodes
      << ", usingPBC=" << mesh.usingPBC
      << ", updateState=" << updateStateName(mesh.getUpdateState()) << "\n";
  return oss.str();
}

std::string formatReductionExplosion(const TElement &element, const Mesh &mesh,
                                     std::string_view context) {
  std::ostringstream oss;
  oss << "Reduction exploded in " << context << ": eIndex=" << element.eIndex
      << ", m3Nr=" << element.m3Nr << ", load="
      << std::setprecision(std::numeric_limits<double>::max_digits10)
      << mesh.load << ", loadSteps=" << mesh.loadSteps << ".\n"
      << formatMeshContext(mesh) << formatElementState(element, "Element");
  return oss.str();
}

std::string formatLagrangeReductionFailure(const TElement &element) {
  std::ostringstream oss;
  oss << "Lagrange reduction failed for fixed-reference state.\n"
      << formatElementState(element, "Element");
  return oss.str();
}

std::string formatMinimizationFailure(std::string_view caughtType,
                                      std::string_view message,
                                      const Mesh &mesh,
                                      std::string_view minimizer, bool rough) {
  std::ostringstream oss;
  oss << "Minimization failed: " << caughtType;
  if (!message.empty()) {
    oss << ": " << message;
  }
  oss << "\n"
      << "Minimization context:\n"
      << "  minimizer=" << minimizer << ", rough=" << rough << "\n"
      << formatMeshContext(mesh);
  return oss.str();
}

std::string formatDisplacementSnapshotSizeError(std::string_view context,
                                                size_t actualSize,
                                                size_t expectedSize, int rows,
                                                int cols) {
  std::ostringstream oss;
  oss << context << ": displacement snapshot size does not match node count."
      << " actualSize=" << actualSize << ", expectedSize=" << expectedSize
      << ", rows=" << rows << ", cols=" << cols << ".";
  return oss.str();
}

std::string formatDisplacementSnapshotPairSizeError(
    std::string_view context, size_t beforeSize, size_t afterSize,
    size_t expectedSize, int rows, int cols) {
  std::ostringstream oss;
  oss << context << ": displacement snapshot sizes do not match node count."
      << " beforeSize=" << beforeSize << ", afterSize=" << afterSize
      << ", expectedSize=" << expectedSize << ", rows=" << rows
      << ", cols=" << cols << ".";
  return oss.str();
}

std::string exceptionMessage(std::exception_ptr exception) {
  if (!exception) {
    return "No exception.";
  }
  try {
    std::rethrow_exception(exception);
  } catch (const alglib::ap_error &e) {
    return "ALGLIB error: " + std::string(e.msg);
  } catch (const std::exception &e) {
    return "Standard exception: " + std::string(e.what());
  } catch (...) {
    return "Unknown exception.";
  }
}

std::string formatDebugReplayFailure(std::string_view originalError,
                                     std::string_view replayError) {
  std::ostringstream oss;
  oss << "Debug replay reproduced an error while logDuringMinimization was "
         "enabled.\n"
      << "Original error:\n"
      << originalError << "\n"
      << "Replay error:\n"
      << replayError;
  return oss.str();
}

std::string formatDebugReplayDidNotReproduce(std::string_view originalError) {
  std::ostringstream oss;
  oss << "Debug replay did not reproduce the original error after enabling "
         "logDuringMinimization.\n"
      << "Original error:\n"
      << originalError;
  return oss.str();
}

} // namespace DebugLog

void printReport(const SimReport &report) {
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgsresults
  std::cout << "Optimization Report:\n";
  std::cout << "\tIterations Count: " << report.nrIter << '\n';
  std::cout << "\tNumber of Function Evaluations: " << report.nfev << '\n';
  std::cout << "\tTermination Reason: ";
  switch (report.termType) {
  case -8:
    std::cout << "Infinite or NAN values in function/gradient";
    break;
  case -2:
    std::cout << "Rounding errors prevent further improvement";
    break;
  case -1:
    std::cout << "Incorrect parameters were specified";
    break;
  case 1:
    std::cout << "Relative function improvement is no more than EpsF";
    break;
  case 2:
    std::cout << "Relative step is no more than EpsX";
    break;
  case 4:
    std::cout << "Gradient norm is no more than EpsG";
    break;
  case 5:
    std::cout << "MaxIts steps was taken";
    break;
  case 7:
    std::cout << "Stopping conditions are too stringent, further improvement "
                 "is impossible";
    break;
  case 8:
    std::cout << "Terminated by user request";
    break;
  default:
    std::cout << "Unknown termination reason";
  }
  std::cout << std::endl;
}

void printReport(const alglib::minlbfgsreport &report) {
  printReport(SimReport(report));
}
void printReport(const alglib::mincgreport &report) {
  printReport(SimReport(report));
}

// New method to print nodeDisplacements in (x, y) pairs
void printNodeDisplacementsGrid(alglib::real_1d_array nodeDisplacements) {
  int nr_x_values = nodeDisplacements.length() / 2;

  std::cout << "Node Displacements (x, y):" << std::endl;

  // Calculate the grid size for printing, assuming a rectangular (not
  // necessarily square) layout
  int gridSizeX = std::ceil(std::sqrt(nr_x_values)); // Width of the grid
  int gridSizeY =
      std::ceil(double(nr_x_values) /
                gridSizeX); // Height of the grid, ensuring all nodes fit

  for (int y = gridSizeY - 1; y >= 0;
       y--) { // Start from the bottom row to have (0,0) in the bottom left
    for (int x = 0; x < gridSizeX; x++) {
      int index = y * gridSizeX + x;
      if (index < nr_x_values) { // Ensure index is within the range of node
                                 // displacements
        std::cout << std::setw(10) << "(" << nodeDisplacements[index] << ", "
                  << nodeDisplacements[index + nr_x_values] << ") ";
      } else {
        // Print placeholders for grid positions without a corresponding node
        std::cout << std::setw(10) << "(--, --) ";
      }
    }
    std::cout << std::endl; // New line for each row of the grid
  }
}
