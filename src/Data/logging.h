#ifndef LOGGING_H
#define LOGGING_H
#pragma once

#include <cereal/access.hpp>
#include <cereal/types/chrono.hpp> // Include Cereal support for std::chrono types
#include <cereal/types/map.hpp>
#include <chrono>
#include <cstddef>
#include <map>
#include <optimization.h>

#define DEFAULT_KEY "main"

// A timer using a key to keep track of multiple timers at the same time
class Timer {
public:
  Timer();
  void Start(const std::string &key = DEFAULT_KEY);
  // Adds the time from the checkpoint to the runtime
  void Save(const std::string &key = DEFAULT_KEY);
  // Returns number of milliseconds taken for minimization
  size_t Stop(const std::string &key = DEFAULT_KEY);
  std::chrono::milliseconds RunTime(const std::string &key = DEFAULT_KEY) const;
  std::string RunTimeString(const std::string &key = DEFAULT_KEY) const;
  void Reset(const std::string &key = DEFAULT_KEY);
  void PrintAllRuntimes() const;

private:
  std::map<std::string, std::chrono::milliseconds> runtimes;
  std::map<std::string, std::chrono::high_resolution_clock::time_point>
      checkpoints;
  std::map<std::string, bool> runningStatus;
  friend class cereal::access;
  template <class Archive> void serialize(Archive &ar) {
    ar(runningStatus, checkpoints, runtimes);
  }
};

std::string FormatDuration(std::chrono::milliseconds duration);

struct SimReport {
  int terminationType;
  size_t nrIter; // Number of itterations
  size_t nfev;   // Number of function evaluations
  size_t nms;    // Number of milliseconds taken
  SimReport(alglib::minlbfgsreport r) {
    terminationType = r.terminationtype;
    nrIter = r.iterationscount;
    nfev = r.nfev;
  }
  SimReport()
      : terminationType(0), nrIter(0), nfev(0) {
  } // Explicit default constructor
};

void printReport(const alglib::minlbfgsreport &report);

// Function to calculate the Estimated Time Remaining (ETR) using progress
// fraction
std::chrono::milliseconds calculateETR(std::chrono::milliseconds elapsed,
                                       float progressFraction);

// Debug function to see nodeDisplacements
void printNodeDisplacementsGrid(alglib::real_1d_array nodeDisplacements);

#endif