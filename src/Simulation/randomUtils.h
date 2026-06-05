#ifndef RANDOM_UTILS_H
#define RANDOM_UTILS_H

#include <cmath>
#include <cstdint>
#include <random>
#include <sstream>
#include <string>

// Rationale:
// - Avoid std::normal_distribution because its caching/state and implementation
//   differ across libstdc++/libc++/MSVC and between versions. That can yield
//   different sequences in Debug vs Release or across toolchains.
// - Implement our own deterministic transforms on top of std::mt19937 so that
//   given the same seed and the same call order, results are identical across
//   build types.
// - Single-thread assumption (caller guarantees no concurrent use).
//
// NOTE1: If a simulation is loaded from a dump, the RNG state is not restored.
// NOTE2: Reproducibility requires that the *number and order* of calls stay
//       the same between Debug/Release. If optimization changes code paths that
//       call the RNG, results will still diverge. This header only guarantees
//       that each primitive call is deterministic and toolchain-stable.

namespace rng {

// Global engine (single-threaded use only)
inline std::mt19937 &engine() {
  static std::mt19937 gen; // default-seeded; call setSeed() before use
  return gen;
}

// Set deterministic seed (32-bit to match std::mt19937 expectation)
inline void setSeed(uint32_t seed) { engine().seed(seed); }

// --- Engine state I/O (optional helpers for checkpointing) ---
inline std::string dumpState() {
  std::ostringstream oss;
  oss << engine();
  return oss.str();
}

inline void loadState(const std::string &state) {
  std::istringstream iss(state);
  iss >> engine();
}

// --- Stable uniform primitives ---
// Returns U in [0,1) using a fixed mapping from 32-bit output to double.
// This avoids library-specific uniform_real_distribution differences.
inline double uniform01() {
  // Pull 32 bits and scale by 2^-32 to get [0,1). This never returns 1.0.
  // Use static_cast to avoid implicit promotions that could vary with flags.
  uint32_t x = engine()();
  // 1.0 / 2^32
  constexpr double inv_2_32 = 1.0 / 4294967296.0; // 2^32
  return static_cast<double>(x) * inv_2_32;
}

// Uniform in (0,1]; avoids exact 0 which would cause log(0) in Box–Muller.
inline double uniform_open01_inclusive() {
  // Map x in [0, 2^32-1] to (0,1] by (x+1)/2^32
  uint32_t x = engine()();
  constexpr double inv_2_32 = 1.0 / 4294967296.0;
  return (static_cast<double>(x) + 1.0) * inv_2_32; // in (0,1]
}

// --- Stable normal via Marsaglia polar method (Box–Muller variant) ---
// Deterministic across platforms given our uniform01() above.
inline double standardNormal() {
  // Cache one spare sample to match typical distribution behavior while being
  // fully deterministic under fixed call order.
  static bool hasSpare = false;
  static double spare = 0.0;

  if (hasSpare) {
    hasSpare = false;
    return spare;
  }

  double u, v, s;
  do {
    // u, v in (-1, 1)
    u = 2.0 * uniform01() - 1.0;
    v = 2.0 * uniform01() - 1.0;
    s = u * u + v * v;
  } while (s >= 1.0 || s == 0.0);

  const double m = std::sqrt(-2.0 * std::log(s) / s);
  spare = v * m;
  hasSpare = true;
  return u * m;
}

inline double sampleNormal(double mean, double stddev) {
  return mean + stddev * standardNormal();
}

inline double sampleLogNormal(double mean, double stddev) {
  return std::exp(sampleNormal(mean, stddev));
}

} // namespace rng

// Backward-compatible API (optional): keep previous function names
inline std::mt19937 &globalGenerator() { return rng::engine(); }
inline void setSeed(unsigned int seed) { rng::setSeed(seed); }
inline double sampleNormal(double mean, double stddev) {
  return rng::sampleNormal(mean, stddev);
}
inline double sampleLogNormal(double mean, double stddev) {
  return rng::sampleLogNormal(mean, stddev);
}

#endif // RANDOM_UTILS_H