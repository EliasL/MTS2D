#ifndef LOGGING_H
#define LOGGING_H
#include "cereal/cereal.hpp"
#include <ap.h>
#pragma once
#include "cereal_help.h"
#include <cereal/access.hpp>
#include <cereal/types/chrono.hpp> // Include Cereal support for std::chrono types
#include <cereal/types/map.hpp>
#include <chrono>
#include <cstddef>
#include <deque>
#include <map>
#include <optimization.h>
#include <string>

#define DEFAULT_KEY "main"

// A timer using a key to keep track of multiple timers at the same time
class Timer {
public:
  Timer();
  void Start(const std::string &key = DEFAULT_KEY);
  // Adds the time from the checkpoint to the runtime
  void Save(const std::string &key = DEFAULT_KEY);
  void SaveAll();
  // Returns number of milliseconds taken for minimization
  size_t Stop(const std::string &key = DEFAULT_KEY);
  std::chrono::milliseconds RunTime(const std::string &key = DEFAULT_KEY) const;
  // Calculates an estimated time remaining based on an average sored in the
  // timer
  std::chrono::milliseconds ETR(double);
  // Runtime string
  std::string RTString(const std::string &key = DEFAULT_KEY,
                       int precision = 3) const;
  // Estimated runtime string
  std::string ETRString(double progress);

  // used to access the ETR without updating it
  std::string oldETRString;

  void Reset(const std::string &key = DEFAULT_KEY);
  void PrintAllRuntimes() const;

private:
  std::map<std::string, std::chrono::milliseconds> runtimes;
  std::map<std::string, std::chrono::high_resolution_clock::time_point>
      checkpoints;
  std::map<std::string, bool> runningStatus;

  // We keep one array that we can average over to calculate ETR
  // We use a pair of (time, completion) where completion goes from 0 to 1
  std::deque<std::chrono::milliseconds> average_time;
  std::deque<double> average_completion;

  friend class cereal::access;
  template <class Archive> void save(Archive &ar) const {
    // Convert runtimes to a serializable format (milliseconds as integers)
    std::map<std::string, int64_t> tempRuntimes;
    for (const auto &pair : runtimes) {
      tempRuntimes[pair.first] =
          pair.second.count(); // Convert chrono::milliseconds to int64_t
    }

    // Adjust runtimes for timers still running
    for (const auto &pair : runningStatus) {
      if (pair.second) { // Timer is running
        auto now = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
            now - checkpoints.at(pair.first));
        tempRuntimes[pair.first] += elapsed.count();
      }
    }

    // Serialize using MAKE_NVP for proper XML structure
    ar(MAKE_NVP(runningStatus), cereal::make_nvp("runtimes", tempRuntimes));
  }

  template <class Archive> void load(Archive &ar) {
    // Temporary map to load runtime data
    std::map<std::string, int64_t> tempRuntimes;

    // Deserialize
    ar(MAKE_NVP(runningStatus), cereal::make_nvp("runtimes", tempRuntimes));

    // Convert runtimes back to std::chrono::milliseconds
    runtimes.clear();
    for (const auto &pair : tempRuntimes) {
      runtimes[pair.first] = std::chrono::milliseconds(pair.second);
    }

    // Restart timers that were running
    for (const auto &pair : runningStatus) {
      if (pair.second) {
        checkpoints[pair.first] = std::chrono::high_resolution_clock::now();
      }
    }
  }
};

std::string FormatDuration(std::chrono::milliseconds duration,
                           int precision = 3);

struct MinState {
  double *energy;
  double *epsf;
  double *epsg;
  double *epsx;
  bool *userStop;
  long *preConType;
  long *n;
  alglib_impl::ae_int_t *nrIt;
  alglib_impl::ae_vector *gradient;
  alglib_impl::ae_vector *work;
  alglib_impl::ae_vector *scale;
  long *termType;

  MinState() {};
  MinState(const alglib::minlbfgsstate &s) {
    auto state = s.c_ptr();
    energy = &s.f;
    epsf = &state->epsf;
    epsg = &state->epsg;
    epsx = &state->epsx;
    userStop = &state->userterminationneeded;
    preConType = &state->prectype;
    n = &state->n;
    nrIt = &state->repiterationscount;
    gradient = &state->g;
    work = &state->work;
    scale = &state->s;
    termType = &state->repterminationtype;
  };

  double gNorm() {

    double v = 0;
    long i = 0;
    v = (double)(0);
    for (i = 0; i <= *n - 1; i++) {
      v = v + pow(gradient->ptr.p_double[i] * scale->ptr.p_double[i], 2);
    }
    return sqrt(v);
  }
};

struct SimReport {
  int termType;  // Termination type
  size_t nrIter; // Number of itterations
  size_t nfev;   // Number of function evaluations
  size_t nms;    // Number of milliseconds taken

  SimReport(const alglib::minlbfgsreport &r) {
    termType = r.terminationtype;
    nrIter = r.iterationscount;
    nfev = r.nfev;
  }
  SimReport(const alglib::mincgreport &r) {
    termType = r.terminationtype;
    nrIter = r.iterationscount;
    nfev = r.nfev;
  }
  SimReport()
      : termType(0), nrIter(0), nfev(0) {} // Explicit default constructor
};

void printReport(const SimReport &report);
void printReport(const alglib::minlbfgsreport &report);
void printReport(const alglib::mincgreport &report);

// Function to calculate the Estimated Time Remaining (ETR) using progress
// fraction
std::chrono::milliseconds calculateETR(std::chrono::milliseconds elapsed,
                                       float progressFraction);

// Debug function to see nodeDisplacements
void printNodeDisplacementsGrid(alglib::real_1d_array nodeDisplacements);

#endif