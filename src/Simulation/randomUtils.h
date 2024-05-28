#ifndef RANDOM_UTILS_H
#define RANDOM_UTILS_H

#include <random>

// Declare the global random number generator
inline std::mt19937 &globalGenerator() {
  static std::mt19937 gen;
  return gen;
}

// Function to set the global seed
inline void setSeed(unsigned int seed) { globalGenerator().seed(seed); }

// Function to sample from a normal distribution
inline double sampleNormal(double mean, double stddev) {
  std::normal_distribution<> d(mean, stddev);
  return d(globalGenerator());
}

#endif // RANDOM_UTILS_H