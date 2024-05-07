#include "logging.h"
#include <chrono>
#include <iomanip>
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

void Timer::Stop(const std::string &key) {
  if (runningStatus[key]) {
    auto now = std::chrono::high_resolution_clock::now();
    runtimes[key] += std::chrono::duration_cast<std::chrono::milliseconds>(
        now - checkpoints[key]);
    runningStatus[key] = false;
  }
}

std::string Timer::RunTimeString(const std::string &key) const {
  return FormatDuration(RunTime(key));
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
