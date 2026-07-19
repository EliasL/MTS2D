#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif

#include "Simulation/simulation.h"

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#include <unistd.h>
#elif defined(__APPLE__)
#include <unistd.h>
#endif

namespace {

struct Options {
  int rows = 32;
  int cols = 32;
  int threads = 1;
  int seed = 1729;
  long long calls = 1000;
  long long warmupCalls = 10;
  double shear = 0.15;
  double perturbation = 1e-3;
  double disorder = 0.0;
  std::string mode = "force-evaluation";
  std::string scenario = "controlled";
  std::string reconnection = "none";
};

struct WorkloadMetrics {
  long long functionEvaluations = 0;
  long long minimizerIterations = 0;
  long long edgeFlips = 0;
  int reconnectCycles = 0;
  int plasticElements = 0;
  int terminationType = 0;
  double finalEnergy = 0.0;
  double finalMaxForce = 0.0;
};

struct ThreadPlacement {
  int thread = -1;
  int cpu = -1;
  int place = -1;
};

std::string jsonEscape(const std::string &value) {
  std::ostringstream out;
  for (const unsigned char c : value) {
    switch (c) {
    case '\\':
      out << "\\\\";
      break;
    case '"':
      out << "\\\"";
      break;
    case '\n':
      out << "\\n";
      break;
    case '\r':
      out << "\\r";
      break;
    case '\t':
      out << "\\t";
      break;
    default:
      if (c < 0x20) {
        out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
            << static_cast<int>(c) << std::dec << std::setfill(' ');
      } else {
        out << static_cast<char>(c);
      }
    }
  }
  return out.str();
}

const char *environmentValue(const char *name) {
  const char *value = std::getenv(name);
  return value == nullptr ? "" : value;
}

long long parseInteger(const char *text, const std::string &name) {
  std::size_t consumed = 0;
  long long value = 0;
  try {
    value = std::stoll(text, &consumed);
  } catch (const std::exception &) {
    throw std::invalid_argument("Invalid integer for " + name + ": " + text);
  }
  if (consumed != std::strlen(text)) {
    throw std::invalid_argument("Invalid integer for " + name + ": " + text);
  }
  return value;
}

double parseDouble(const char *text, const std::string &name) {
  std::size_t consumed = 0;
  double value = 0.0;
  try {
    value = std::stod(text, &consumed);
  } catch (const std::exception &) {
    throw std::invalid_argument("Invalid number for " + name + ": " + text);
  }
  if (consumed != std::strlen(text)) {
    throw std::invalid_argument("Invalid number for " + name + ": " + text);
  }
  return value;
}

void printHelp() {
  std::cout
      << "Usage: benchmark_MTS2D [options]\n\n"
      << "Times output-free force or minimization workloads and writes one JSON "
         "record to stdout.\n\n"
      << "Options:\n"
      << "  --rows N                 Mesh rows (default: 32)\n"
      << "  --cols N                 Mesh columns (default: 32)\n"
      << "  --threads N              OpenMP threads (default: 1)\n"
      << "  --seed N                 Deterministic fixture seed (default: 1729)\n"
      << "  --calls N                Timed evaluations (default: 1000)\n"
      << "  --warmup-calls N         Untimed evaluations (default: 10)\n"
      << "  --mode NAME              force-kernel, force-evaluation, or minimization\n"
      << "  --scenario NAME          controlled, quiet, or event\n"
      << "  --reconnection NAME      none or edgeFlip (minimization only)\n"
      << "  --shear X                Affine shear in the fixture (default: 0.15)\n"
      << "  --perturbation X         Non-affine amplitude (default: 0.001)\n"
      << "  --disorder X             Quenched-disorder SD (default: 0)\n"
      << "  --help                    Show this help\n";
}

Options parseOptions(int argc, char **argv) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string argument = argv[i];
    if (argument == "--help" || argument == "-h") {
      printHelp();
      std::exit(0);
    }
    if (i + 1 >= argc) {
      throw std::invalid_argument("Missing value after " + argument);
    }
    const char *value = argv[++i];
    if (argument == "--rows") {
      options.rows = static_cast<int>(parseInteger(value, argument));
    } else if (argument == "--cols") {
      options.cols = static_cast<int>(parseInteger(value, argument));
    } else if (argument == "--threads") {
      options.threads = static_cast<int>(parseInteger(value, argument));
    } else if (argument == "--seed") {
      options.seed = static_cast<int>(parseInteger(value, argument));
    } else if (argument == "--calls") {
      options.calls = parseInteger(value, argument);
    } else if (argument == "--warmup-calls") {
      options.warmupCalls = parseInteger(value, argument);
    } else if (argument == "--mode") {
      options.mode = value;
    } else if (argument == "--scenario") {
      options.scenario = value;
    } else if (argument == "--reconnection") {
      options.reconnection = value;
    } else if (argument == "--shear") {
      options.shear = parseDouble(value, argument);
    } else if (argument == "--perturbation") {
      options.perturbation = parseDouble(value, argument);
    } else if (argument == "--disorder") {
      options.disorder = parseDouble(value, argument);
    } else {
      throw std::invalid_argument("Unknown option: " + argument);
    }
  }

  if (options.rows < 2 || options.cols < 2) {
    throw std::invalid_argument("Rows and columns must both be at least 2.");
  }
  if (options.threads < 1) {
    throw std::invalid_argument("Thread count must be positive.");
  }
  if (options.calls < 1 || options.warmupCalls < 0) {
    throw std::invalid_argument(
        "Calls must be positive and warm-up calls must be non-negative.");
  }
  if (options.disorder < 0.0) {
    throw std::invalid_argument("Disorder must be non-negative.");
  }
  if (options.mode != "force-kernel" && options.mode != "force-evaluation" &&
      options.mode != "minimization") {
    throw std::invalid_argument(
        "Mode must be force-kernel, force-evaluation, or minimization.");
  }
  if (options.scenario != "controlled" && options.scenario != "quiet" &&
      options.scenario != "event") {
    throw std::invalid_argument("Scenario must be controlled, quiet, or event.");
  }
  if (options.reconnection != "none" && options.reconnection != "edgeFlip") {
    throw std::invalid_argument("Reconnection must be none or edgeFlip.");
  }
  if (options.mode == "minimization" && options.calls != 1) {
    throw std::invalid_argument(
        "Minimization samples must use exactly one independent minimization.");
  }
  if (options.mode == "minimization" && options.warmupCalls != 0) {
    throw std::invalid_argument(
        "Minimization samples do not use warm-up minimizations; pass "
        "--warmup-calls 0.");
  }
  if (options.mode == "minimization" && options.scenario == "controlled") {
    throw std::invalid_argument(
        "Minimization mode requires the quiet or event scenario.");
  }
  if (options.mode != "minimization" && options.reconnection != "none") {
    throw std::invalid_argument(
        "Reconnection is only supported by minimization mode.");
  }
  return options;
}

std::string compilerDescription() {
#if defined(__clang__)
  return std::string("Clang ") + __clang_version__;
#elif defined(__GNUC__)
  return std::string("GCC ") + __VERSION__;
#else
  return "unknown";
#endif
}

std::string hostName() {
  char buffer[256] = {};
#if defined(__linux__) || defined(__APPLE__)
  if (gethostname(buffer, sizeof(buffer) - 1) == 0) {
    return buffer;
  }
#endif
  return "unknown";
}

std::vector<ThreadPlacement> captureThreadPlacement(int requestedThreads,
                                                    int &observedThreads) {
  std::vector<ThreadPlacement> placement(
      static_cast<std::size_t>(requestedThreads));
  observedThreads = 0;
#pragma omp parallel
  {
    const int thread = omp_get_thread_num();
#pragma omp single
    observedThreads = omp_get_num_threads();

    if (thread < requestedThreads) {
      ThreadPlacement item;
      item.thread = thread;
      item.place = omp_get_place_num();
#if defined(__linux__)
      item.cpu = sched_getcpu();
#endif
      placement[static_cast<std::size_t>(thread)] = item;
    }
  }
  placement.resize(static_cast<std::size_t>(observedThreads));
  return placement;
}

std::vector<double> makeDisplacements(const Mesh &mesh, double amplitude) {
  constexpr double twoPi = 6.283185307179586476925286766559;
  const std::size_t count = mesh.freeNodeIds.size();
  std::vector<double> displacement(2 * count, 0.0);
  for (std::size_t i = 0; i < count; ++i) {
    const NodeId &id = mesh.freeNodeIds[i];
    const double x = static_cast<double>(id.col()) /
                     static_cast<double>(std::max(1, mesh.cols));
    const double y = static_cast<double>(id.row()) /
                     static_cast<double>(std::max(1, mesh.rows));
    displacement[i] = amplitude * std::sin(twoPi * x) * std::cos(twoPi * y);
    displacement[count + i] =
        amplitude * std::cos(twoPi * x) * std::sin(twoPi * y);
  }
  return displacement;
}

double copyGradient(const Mesh &mesh, std::vector<double> &gradient) {
  const std::size_t count = mesh.freeNodeIds.size();
  double probe = 0.0;
  for (std::size_t i = 0; i < count; ++i) {
    const Node *node = mesh[mesh.freeNodeIds[i]];
    gradient[i] = -node->f.x();
    gradient[count + i] = -node->f.y();
  }
  if (!gradient.empty()) {
    probe = gradient.front() + gradient[gradient.size() / 2] + gradient.back();
  }
  return probe;
}

double evaluate(Mesh &mesh, const Options &options,
                const std::vector<double> &displacement,
                std::vector<double> &gradient) {
  if (options.mode == "force-evaluation") {
    mesh.updateNodePositions(displacement.data(), displacement.size());
  } else {
    // ensureForces() intentionally caches an unchanged force state. Marking the
    // mesh dirty is required to benchmark the arithmetic rather than that
    // software-level cache hit.
    mesh.markDirty();
  }
  mesh.updateMesh();

  double checksum = mesh.totalEnergy + mesh.maxForce;
  if (options.mode == "force-evaluation") {
    checksum += copyGradient(mesh, gradient);
  }
  return checksum;
}

template <typename Duration>
double seconds(Duration duration) {
  return std::chrono::duration<double>(duration).count();
}

void printJson(const Options &options, const Mesh &mesh, int observedThreads,
               const std::vector<ThreadPlacement> &placement,
               double initializationSeconds, double warmupSeconds,
               double elapsedSeconds, double checksum,
               const WorkloadMetrics &workload) {
  const double secondsPerCall = elapsedSeconds / options.calls;
  const double callsPerSecond = options.calls / elapsedSeconds;
  const double elementsPerSecond =
      static_cast<double>(mesh.nrElements) * callsPerSecond;
  const unsigned long long scratchBytes =
      static_cast<unsigned long long>(observedThreads) *
      static_cast<unsigned long long>(mesh.nrNodes) * sizeof(Vector2d);

  std::ostringstream out;
  out << std::setprecision(17);
  out << "{";
  out << "\"schema_version\":2";
  out << ",\"status\":\"ok\"";
  out << ",\"mode\":\"" << jsonEscape(options.mode) << "\"";
  out << ",\"scenario\":\"" << jsonEscape(options.scenario) << "\"";
  out << ",\"reconnection\":\"" << jsonEscape(options.reconnection)
      << "\"";
  out << ",\"rows\":" << options.rows;
  out << ",\"cols\":" << options.cols;
  out << ",\"nodes\":" << mesh.nrNodes;
  out << ",\"elements\":" << mesh.nrElements;
  out << ",\"threads_requested\":" << options.threads;
  out << ",\"threads_observed\":" << observedThreads;
  out << ",\"fixture_seed\":" << options.seed;
  out << ",\"calls\":" << options.calls;
  out << ",\"warmup_calls\":" << options.warmupCalls;
  out << ",\"shear\":" << options.shear;
  out << ",\"perturbation\":" << options.perturbation;
  out << ",\"disorder\":" << options.disorder;
  out << ",\"initialization_seconds\":" << initializationSeconds;
  out << ",\"warmup_seconds\":" << warmupSeconds;
  out << ",\"elapsed_seconds\":" << elapsedSeconds;
  out << ",\"seconds_per_call\":" << secondsPerCall;
  out << ",\"calls_per_second\":" << callsPerSecond;
  out << ",\"elements_per_second\":" << elementsPerSecond;
  out << ",\"checksum\":" << checksum;
  out << ",\"function_evaluations\":" << workload.functionEvaluations;
  out << ",\"minimizer_iterations\":" << workload.minimizerIterations;
  out << ",\"seconds_per_function_evaluation\":"
      << (workload.functionEvaluations > 0
              ? elapsedSeconds / workload.functionEvaluations
              : 0.0);
  out << ",\"edge_flips\":" << workload.edgeFlips;
  out << ",\"reconnect_cycles\":" << workload.reconnectCycles;
  out << ",\"plastic_elements\":" << workload.plasticElements;
  out << ",\"termination_type\":" << workload.terminationType;
  out << ",\"final_energy\":" << workload.finalEnergy;
  out << ",\"final_max_force\":" << workload.finalMaxForce;
  out << ",\"force_scratch_bytes_estimate\":" << scratchBytes;
  out << ",\"host\":\"" << jsonEscape(hostName()) << "\"";
  out << ",\"compiler\":\"" << jsonEscape(compilerDescription()) << "\"";
#ifdef _OPENMP
  out << ",\"openmp_version\":" << _OPENMP;
#else
  out << ",\"openmp_version\":0";
#endif
  out << ",\"openmp_num_places\":" << omp_get_num_places();
  out << ",\"openmp_proc_bind\":"
      << static_cast<int>(omp_get_proc_bind());
  out << ",\"omp_places_env\":\""
      << jsonEscape(environmentValue("OMP_PLACES")) << "\"";
  out << ",\"omp_proc_bind_env\":\""
      << jsonEscape(environmentValue("OMP_PROC_BIND")) << "\"";
  out << ",\"thread_placement\":[";
  for (std::size_t i = 0; i < placement.size(); ++i) {
    if (i != 0) {
      out << ',';
    }
    out << "{\"thread\":" << placement[i].thread
        << ",\"cpu\":" << placement[i].cpu
        << ",\"place\":" << placement[i].place << '}';
  }
  out << "]}";
  std::cout << out.str() << '\n';
}

} // namespace

int main(int argc, char **argv) {
  try {
    const Options options = parseOptions(argc, argv);

    omp_set_dynamic(0);
    omp_set_num_threads(options.threads);
    omp_set_max_active_levels(1);

    int observedThreads = 0;
    const std::vector<ThreadPlacement> placement =
        captureThreadPlacement(options.threads, observedThreads);
    if (observedThreads != options.threads) {
      throw std::runtime_error(
          "OpenMP created " + std::to_string(observedThreads) +
          " threads, but " + std::to_string(options.threads) +
          " were requested. Check the scheduler allocation, OMP_THREAD_LIMIT, "
          "and affinity settings.");
    }

    if (options.mode == "minimization") {
      const auto initializationStart = std::chrono::steady_clock::now();
      Config config;
      config.setDefaultValues();
      config.name = "benchmark";
      config.rows = options.rows;
      config.cols = options.cols;
      config.usingPBC = true;
      config.reconnectionMethod = options.reconnection;
      config.nrThreads = options.threads;
      config.seed = options.seed;
      config.QDSD = options.disorder;
      config.showProgress = -1;
      config.logDuringMinimization = false;
      config.fullMinimizationLogging = false;
      config.writeDumps = false;
      config.writeDebugVTUs = false;
      config.forceReRun = true;

      Simulation simulation(config, "", false);
      simulation.initSolver();

      // Relax the disordered zero-load fixture before constructing the timed
      // load step. This preparation is deliberately outside the timed region:
      // a production simulation arrives at both quiet and eventful steps from
      // an already-relaxed state rather than the freshly constructed mesh.
      simulation.applyLoadStepToGuess();
      simulation.minimize(false);
      simulation.mesh.resetCounters();

      simulation.mesh.addLoad(options.shear);
      const Matrix2d loadTransform = getShear(options.shear);
      simulation.mesh.applyTransformationToSystemDeformation(loadTransform);
      simulation.applyLoadStepToGuess(loadTransform);
      simulation.addNoiseToGuess(options.perturbation);
      const auto initializationEnd = std::chrono::steady_clock::now();

      const auto timedStart = std::chrono::steady_clock::now();
      simulation.minimize(options.reconnection == "edgeFlip");
      const auto timedEnd = std::chrono::steady_clock::now();

      WorkloadMetrics workload;
      workload.functionEvaluations = simulation.mesh.nrMinFunctionCalls;
      workload.minimizerIterations = simulation.mesh.nrMinItterations;
      workload.edgeFlips =
          static_cast<long long>(simulation.mesh.totalEdgeFlipsInStep);
      workload.reconnectCycles = simulation.reconnectingCycles();
      workload.plasticElements =
          simulation.mesh.nr_elements_with_m3_changeInStep;
      workload.terminationType = simulation.LBFGSRep.termType;
      workload.finalEnergy = simulation.mesh.totalEnergy;
      workload.finalMaxForce = simulation.mesh.maxForce;
      const double checksum =
          workload.finalEnergy + workload.finalMaxForce +
          static_cast<double>(workload.functionEvaluations) +
          static_cast<double>(workload.edgeFlips);

      printJson(options, simulation.mesh, observedThreads, placement,
                seconds(initializationEnd - initializationStart), 0.0,
                seconds(timedEnd - timedStart), checksum, workload);
    } else {
      const auto initializationStart = std::chrono::steady_clock::now();
      Mesh mesh(options.rows, options.cols, 1.0, 0.0, true, "major",
                "contiSquare", 4.0);
      Matrix2d shear = Matrix2d::Identity();
      shear(0, 1) = options.shear;
      mesh.applyTransformationToSystemDeformation(shear);
      const std::vector<double> displacement =
          makeDisplacements(mesh, options.perturbation);
      std::vector<double> gradient(displacement.size(), 0.0);
      mesh.updateNodePositions(displacement.data(), displacement.size());
      mesh.updateMesh();
      const auto initializationEnd = std::chrono::steady_clock::now();

      double checksum = 0.0;
      const auto warmupStart = std::chrono::steady_clock::now();
      for (long long i = 0; i < options.warmupCalls; ++i) {
        checksum += evaluate(mesh, options, displacement, gradient);
      }
      const auto warmupEnd = std::chrono::steady_clock::now();

      checksum = 0.0;
      const auto timedStart = std::chrono::steady_clock::now();
      for (long long i = 0; i < options.calls; ++i) {
        checksum += evaluate(mesh, options, displacement, gradient);
      }
      const auto timedEnd = std::chrono::steady_clock::now();

      WorkloadMetrics workload;
      workload.finalEnergy = mesh.totalEnergy;
      workload.finalMaxForce = mesh.maxForce;
      printJson(options, mesh, observedThreads, placement,
                seconds(initializationEnd - initializationStart),
                seconds(warmupEnd - warmupStart),
                seconds(timedEnd - timedStart), checksum, workload);
    }
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "benchmark_MTS2D: " << error.what() << '\n';
    return 2;
  }
}
