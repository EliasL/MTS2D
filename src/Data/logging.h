#ifndef LOGGING_H
#define LOGGING_H
#pragma once

#include <cereal/access.hpp>
#include <cereal/types/chrono.hpp> // Include Cereal support for std::chrono types
#include <cereal/types/map.hpp>
#include <chrono>
#include <map>

#define DEFAULT_KEY "main"

// A timer using a key to keep track of multiple timers at the same time
class Timer {
public:
  Timer();
  void Start(const std::string &key = DEFAULT_KEY);
  void Stop(const std::string &key = DEFAULT_KEY);
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

#endif