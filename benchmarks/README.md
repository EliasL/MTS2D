# MTS2D benchmark suite

The benchmark suite measures controlled force evaluations and output-free full
minimizations. It is portable across macOS development machines and Linux
compute nodes and has no required Python packages beyond the standard library.

The suite has two layers:

- `benchmark_MTS2D` is the native timed executable. It initializes one mesh,
  performs untimed warm-up calls for force modes, or prepares one relaxed
  minimization fixture, times with `steady_clock`, verifies the requested
  OpenMP team size, consumes the numerical result, and emits one JSON record.
- `tools/run_benchmarks.py` builds the native target, generates system/strong/
  weak-scaling cases, calibrates safe call counts, runs one child process at a
  time, randomizes repetitions, enforces budgets, and writes statistics,
  plots, and hardware-specific configuration candidates.

## Benchmark modes

- `force-kernel` marks an unchanged mesh dirty and recomputes element forces,
  energy, per-thread scratch, and nodal force accumulation. This is the most
  controlled kernel measurement.
- `force-evaluation` additionally copies a deterministic displacement field
  into the nodes and copies nodal forces into an optimizer-style gradient.
  This is the preferred approximation to one minimizer function evaluation.

Both modes repeatedly use one deterministic mesh state after warm-up. They do
not include mesh construction, random-state generation, logging, VTU output,
relaxation, or reconnection in the timed region.

The `minimization` mode measures the entire relaxation call, including solver
overhead and optional edge-flip reconnection. Mesh construction and an initial
zero-load relaxation prepare the state outside the timed region. Each timed
sample runs one independent minimization in a fresh process.

## Quiet, event, and reconnection workloads

The workload experiment uses a compact three-cell design:

- `quiet + none`: a tiny load step from a relaxed state. It normally needs only
  a few function evaluations and exposes OpenMP barrier/team overhead.
- `event + none`: a heterogeneous, high-activity synthetic step with many
  function evaluations.
- `event + edgeFlip`: the same event fixture with only the production
  `edgeFlip` reconnection method enabled.

`quiet + edgeFlip` is deliberately omitted. The simulation itself skips
reconnection when minimization produces no plastic-state change, so this cell
would duplicate `quiet + none`. Delaunay reconnection is never selected by the
suite.

Workload repetitions use consecutive deterministic seeds, matched across every
thread count and binding policy. Their standard deviation therefore includes
meaningful fixture-to-fixture variation as well as timing noise. Controlled
force repetitions reuse the same data, so their standard deviation measures
timing variation only.

This is a fractional factorial screening design, not a formal Taguchi array.
It keeps the interactions most likely to matter—size × threads, workload ×
threads, and event × reconnection—while avoiding a combinatorial sweep of all
simulation parameters. A generic orthogonal array would be less useful here
because it can alias precisely those interactions. Promising settings are
confirmed with repeated matched fixtures.

The built-in event is intentionally a synthetic stress fixture (default load
step `0.72`, quenched-disorder SD `0.3`, perturbation `0.02`). It is useful for
portable screening and reliably exercises nontrivial minimizations, but it is
not a claim that every scientific configuration behaves that way. For final
paper results, preserve dumps immediately before representative quiet and
late-load events (for example with `--makeDumpAt`) so a future replay fixture
can be compared against this screening result.

## Local use

Run the tiny validation suite while developing:

```sh
python3 tools/run_benchmarks.py --preset smoke
```

Run the default suite, designed to finish within one hour:

```sh
python3 tools/run_benchmarks.py --preset quick
```

Run the higher-accuracy suite, with a six-hour hard total budget and a
20-minute budget for each benchmark case:

```sh
python3 tools/run_benchmarks.py --preset full
```

The full preset includes square sizes through `1024`, but it calibrates the
number of calls independently for each case and never exceeds 10,000 calls per
sample. Slow or memory-heavy cases are reduced, timed out, or skipped rather
than breaking the global budget.

Useful customizations include:

```sh
# Exact thread sweep and selected mesh sizes
python3 tools/run_benchmarks.py --preset custom \
  --threads 1-16 --sizes 64,128,256,512 \
  --strong-size 512 --repetitions 8

# Request 10,000 calls when the time budget permits
python3 tools/run_benchmarks.py --preset custom --calls 10000

# Inspect the generated matrix without building or running it
python3 tools/run_benchmarks.py --preset full --dry-run

# Workload-only screening with matched seeds and selected sizes
python3 tools/run_benchmarks.py --preset custom --experiments workload \
  --workload-sizes 32,64 --threads 1,2,4,8 --repetitions 8
```

By default the benchmark gets a clean optimized build in `build-benchmark`
with `ENABLE_NOINLINE=OFF`. Use `--keep-noinline` only when matching a
profiler-oriented build is intentional. Use `--no-build --exe PATH` to test a
specific existing executable.

## Linux compute nodes and affinity

Resource allocation and thread placement solve different parts of the problem:

1. The batch scheduler grants one task an exclusive set of CPUs.
2. The runner places OpenMP threads inside that granted CPU set.

Force-scaling cases use the first selected policy. With the default
`--affinity auto`, minimization workloads compare both `close` and `spread` so
the recommendation is measured on the current host. The portable constraints
are:

```text
OMP_PLACES=cores
OMP_PROC_BIND=close
OMP_DYNAMIC=FALSE
```

`close` packs threads into nearby OpenMP core places, which usually keeps a
small team within a socket or cache domain. It is not universally optimal:
memory-bandwidth-limited or multi-socket cases may prefer `--affinity spread`.
The suite also supports `--affinity none` for an explicit unbound comparison.
Passing `--affinity close` or `--affinity spread` disables the comparison and
uses only that policy.

For Slurm, request a whole node (or an exclusive CPU allocation), start one
task, and give that task all physical cores that the suite may use. For a
128-core node, a typical job script is:

```sh
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --time=06:15:00

srun --cpu-bind=cores --mem-bind=local \
  python3 tools/run_benchmarks.py --preset full --threads auto
```

Cluster policies differ, so confirm the site's meaning of `--exclusive`, CPU
counts, simultaneous multithreading, and memory binding. Do not make the runner
start multiple benchmark processes concurrently.

On Linux, each native result records the logical CPU ID and OpenMP place of
every thread. The run manifest also records the allowed CPU set, package, die,
core, and NUMA-node topology read from sysfs. These records make it possible to
verify that `close` actually produced the intended placement. On runtimes that
do not expose OpenMP places, the report warns that affinity is unverified.

The runner uses physical cores by default. Pass an explicit `--threads` list if
simultaneous multithreading is itself part of the experiment.

## Time and memory budgets

`quick` defaults to a one-hour global budget. `full` defaults to six hours.
Every case in either preset has a 20-minute total budget covering calibration,
process startup, warm-up, and all repetitions. These are hard subprocess
timeouts, not estimates.

Before running a case, the runner estimates the current dense force scratch as

```text
threads × nodes × sizeof(Vector2d)
```

and skips it if that estimate alone exceeds half of currently available
memory. Change the guard with `--max-memory-fraction` when an allocation makes
a different limit appropriate.

## Results

Results are written below `tmp/benchmark_suite/<preset>_<timestamp>/`:

- `samples.jsonl` is appended after every successful independent process, so
  partial results survive interruption.
- `samples.csv` contains all raw observations and thread-placement records.
- `summary.csv` and `summary.json` contain the arithmetic mean, sample standard
  deviation, median, range, coefficient of variation, standard error,
  throughput, strong/weak efficiency, and checksum validation.
- `manifest.json` records the executable hash, git state, effective build
  flags, system and scheduler metadata, topology, budgets, skips, and failures.
- `report.md` provides a compact human-readable table.
- `plots/` shows mean time per force evaluation or full minimization with
  standard-deviation error bars. Matplotlib produces PNGs when available;
  otherwise the runner writes dependency-free SVG plots.
- `recommendations.json` records conservative candidates separately for quiet,
  event/no-reconnection, and event/edge-flip workloads at every tested size.
- `recommended_configs/` contains one OpenMP environment file and one
  `nrThreads` config fragment per workload profile, selected from the largest
  measured size. Config fragments must be merged into an existing simulation
  config; they are not standalone configurations.

The directory also contains `mixed_none` and `mixed_edgeFlip` profiles for
complete simulations that contain both quiet and event-heavy steps. These use
a weight-free minimax rule: choose the setting whose worst relative slowdown
from either workload's individual optimum is smallest. This avoids assuming an
unknown quiet/event frequency. The per-workload profiles remain useful when the
actual workload mix is known.

Recommendations prefer fewer threads when their mean is within 10% of the
fastest and the difference is no larger than 1.96 combined standard errors.
An edge-flip recommendation is not emitted if the fixture never flips an edge.
The JSON records flip trigger rates, means, standard deviations, fixture seeds,
selection rule, and confidence. Affinity confidence is capped at
`exploratory` when the runtime cannot expose enough placement information to
verify binding.

The existing `tools/benchmark_reconnect.py` remains the end-to-end benchmark
for relaxation/reconnection, output, wall time, memory, and profiling. It is
complementary to this controlled benchmark suite.
