#include "logging.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

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
  } else {
    return runtimes.at(key);
  }
}

std::chrono::milliseconds Timer::ETR(double completion) {
  // Add the current runtime and completion to the deques
  average_time.push_back(RunTime());
  average_completion.push_back(completion);

  // Keep the deque size fixed at 1024 by removing the oldest entry if necessary
  // At 20 seconds per update, that makes for an average period of 5.7 hours
  if (average_time.size() > 1024) {
    average_time.pop_front();
    average_completion.pop_front();
  }

  // Need at least two data points to calculate the rate
  if (average_time.size() > 1) {
    // Initialize variables to calculate total time and completion changes
    double total_time_change = 0.0;       // in milliseconds
    double total_completion_change = 0.0; // in percentage (0.0 to 1.0)

    // Iterate over the deque to calculate total changes
    for (size_t i = 1; i < average_time.size(); ++i) {
      // Get the time and completion values for the current and previous indices
      auto time_current = average_time[i];
      auto time_previous = average_time[i - 1];
      double completion_current = average_completion[i];
      double completion_previous = average_completion[i - 1];

      // Calculate the differences
      auto delta_time = std::chrono::duration<double, std::milli>(time_current -
                                                                  time_previous)
                            .count(); // in milliseconds
      double delta_completion = completion_current - completion_previous;

      // Accumulate the total changes
      total_time_change += delta_time;
      total_completion_change += delta_completion;
    }

    // Avoid division by zero
    if (total_completion_change > 0.0) {
      // Calculate the average rate of completion per millisecond
      double average_rate = total_completion_change / total_time_change;

      // Calculate the remaining completion needed
      double remaining_completion = 1.0 - completion;

      // Estimate the time remaining in milliseconds
      double estimated_time_remaining_ms = remaining_completion / average_rate;

      // Return the estimated time remaining as std::chrono::milliseconds
      return std::chrono::milliseconds(
          static_cast<int64_t>(estimated_time_remaining_ms));
    } else {
      // If no progress was made, return zero or handle accordingly
      return std::chrono::milliseconds(0);
    }
  } else {
    // Not enough data points to estimate
    return std::chrono::milliseconds(0);
  }
}

std::string Timer::RTString(const std::string &key) const {
  return FormatDuration(RunTime(key));
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

std::string FormatDuration(std::chrono::milliseconds duration) {

  bool useMilliseconds = duration.count() < 1e4;
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