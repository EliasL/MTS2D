#include "logging.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

Timer::Timer() {}

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
  if (runningStatus.at(key)) {
    auto now = std::chrono::high_resolution_clock::now();
    auto currentDuration =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            now - checkpoints.at(key));
    return (runtimes.at(key) + currentDuration);
  }
  return runtimes.at(key);
}

std::string Timer::RTString(const std::string &key) const {
  return FormatDuration(RunTime(key));
}

std::string Timer::ETRString(double progress) const {
  return FormatDuration(calculateETR(RunTime(), progress));
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

std::string FormatDuration(std::chrono::milliseconds duration) {
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

  stream << std::fixed << std::setprecision(3)
         << seconds.count() + milliseconds / 1000.0 << "s";
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
  double rate = progressFraction / elapsedSeconds;
  if (rate == 0) {
    return std::chrono::milliseconds::min(); // Avoid infinity if rate
                                             // calculates to zero
  }
  long long etrInMilliseconds =
      static_cast<long long>(((1 - progressFraction) / rate) * 1000);

  // Ignore negative values
  if (etrInMilliseconds < 0) {
    etrInMilliseconds = 0;
  }
  return std::chrono::milliseconds(etrInMilliseconds);
}

void printReport(const alglib::minlbfgsreport &report) {
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgsresults
  std::cout << "Optimization Report:\n";
  std::cout << "\tIterations Count: " << report.iterationscount << '\n';
  std::cout << "\tNumber of Function Evaluations: " << report.nfev << '\n';
  std::cout << "\tTermination Reason: ";
  switch (report.terminationtype) {
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