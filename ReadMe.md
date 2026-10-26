# Mesoscopic Tensorial Simulation (MTS)

## About The Project

This project presents a novel simulation using the Mesoscopic Tensorial Model (MTM) to study dislocation nucleation in defect-free crystals. This simulation aims to explore the formation and behaviour of dislocation nucleation patterns.

## Getting Started

### Prerequisites

Required:
- Git
- CMake (>= 3.14)
- A C++17 compiler (Clang or GCC)
- Make or Ninja
- OpenMP support (on macOS you may need `libomp`)
- Zlib development package
- LLVM/LLDB debugger support

Optional:
- Python 3 (for config generation and plotting scripts)
- CGAL dependencies (CMake will tell you if anything is missing)

Notes:
- ALGLIB, Eigen, and Cereal are included in `libs/` and built by CMake.
- CGAL is not fetched automatically. Place it in `libs/cgal` or update `CGAL_DIR` in `CMakeLists.txt`.
- CGAL is currently optional (used only for Delaunay reconnection). You can build without it using `-DIDE_LIGHTWEIGHT=ON`.
- On macOS, `workingZProfile.txt` provides a minimal Homebrew setup and selects Apple Clang through `CC` and `CXX`. CMake configures OpenMP using Homebrew `libomp`; global compiler/linker flags and a system-directory `PATH` override are unnecessary. This example is not a complete replacement for your personal `~/.zprofile`.

### Installation

1. **Clone the repository**

   ```sh
   git clone https://github.com/EliasL/MTS2D.git
   cd MTS2D
   ```

2. **CGAL setup**

   Place CGAL in `libs/cgal`, or update the `CGAL_DIR` path in `CMakeLists.txt`.

### Build

Release builds live in `build/release`; optional Debug builds and scientific
tests live in `build/debug`. CMake presets require CMake 3.21 or newer.

In the complete MTS workspace, open the parent MTS folder in VS Code. Select
`MTS2D/CMakeLists.txt` and the `release` preset. Use **Build managed simulator**
from **Terminal → Run Build Task** to build the Release executable and record
its identity for managed simulations. From the parent MTS folder, the same
operation is `bash setup.sh --local-sim`, or use the GUI Setup page to rebuild.
Personal VS Code settings are ignored by Git.

For a standalone simulator build, run from this folder:

```sh
cmake --preset release
cmake --build --preset release
```

This standalone build does not create the manager's build manifest. To prepare
it for managed runs, use the managed build task above.

To build and run scientific tests:

```sh
cmake --preset debug
cmake --build --preset debug
./build/debug/test_MTS2D
```

With CMake 3.14–3.20, configure explicitly instead of using presets:
`cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release`, then
`cmake --build build/release --target MTS2D --parallel 2`.

### Run

Basic usage:
```sh
./build/release/MTS2D -c path/to/config.conf -o path/to/output
```

Resume from a dump:
```sh
./build/release/MTS2D -d path/to/dump.xml.gz -o path/to/output
```

Common flags:
- `-c` config file
- `-d` dump file (resume)
- `-o` output path
- `-r` force re-run even if output folder looks complete
- `--makeDumpAt <load>` write a debug dump near the given load

### Tips

LBFGS and CG use Sylvain's sparse tangent-stiffness preconditioning by default.
It changes the coordinates seen by the optimizer; physical energy, forces, and
the `epsR` force tolerance keep their usual meaning. FIRE does not use it.

For comparisons, change just this config entry:

```ini
preconditioner=none
```

Use `preconditioner=exactTangent` to enable it explicitly, or omit the entry.
Use `preconditioner=diagonal` for ALGLIB's native diagonal preconditioning. It
assembles only the tangent diagonal and keeps physical coordinates, without a
sparse matrix or Cholesky factorization. Very small or negative diagonal entries
are floored at `1e-8` times the larger of one and the largest absolute diagonal
entry, since ALGLIB requires positive values.

Both modes use the same rebuild policy. By default, the preconditioner is reused
until nodal displacement since its construction
reaches `preconditionerDisplacementThreshold=0.05`. Set
`preconditionerUpdate=every` to rebuild before each minimization, or `frozen` to
reuse it until the mesh connectivity or free-node layout changes. It is always
held fixed during a minimization. If exact-tangent factorization fails even after
diagonal shifts, that solve uses ordinary minimization.

When debugging a problem that happens late in a simulation, pass
`--makeDumpAt <load>` on the command line with a load value just before the issue.
The simulation can then be resumed with `-d` to restart closer to the problem.

For performance benchmarks or short load windows where VTU output dominates the
runtime, reduce `nrVTUFrames` in the config. The default is `200`, matching the
old hard-coded VTU spacing.

Benchmark runners, the native benchmark driver, plotting scripts, and fixture
catalog live in the sibling SimulationScripts repository. From its root, run
`.venv/bin/python Benchmarks/run_benchmarks.py --preset smoke`; see
[the benchmark documentation](../SimulationScripts/Benchmarks/README.md).
Animation drivers and generated figures also belong in SimulationScripts.

### Example Config

Minimal LBFGS example:
```ini
name=smallSimulation
rows=10
cols=10
boundary=periodic
reconnectionMethod=none
loading=relaxation
deformation=simpleShear
nrThreads=1
seed=0
QDSD=0.0
initialGuessNoise=0.05
meshDiagonal=major
energyFunction=contiSquare
bulkModulus=4.0

startLoad=0.15
loadIncrement=1e-5
totalLoad=0.001

minimizer=LBFGS
epsR=1e-5

LBFGSNrCorrections=10
LBFGSScale=1.0
LBFGSEpsg=1e-8
LBFGSEpsf=0
LBFGSEpsx=0
LBFGSMaxIterations=0

logDuringMinimization=false
fullMinimizationLogging=false
writeDumps=false
nrVTUFrames=200
plasticityEventThreshold=0.01
energyDropThreshold=0.001
showProgress=1
```

### Designing a protocol

A protocol is the whole recipe in the config. One runner combines the choices
below; there is no separate `protocol=standard` selector.

```ini
loading=relaxation
deformation=simpleShear
boundary=periodic
reconnectionMethod=edgeFlip
minimum_transport_distance=true

# One direction; no switching or repetition needed.
directions=horizontal
startLoad=0
loadIncrement=0.001
totalLoad=1
```

- `loading=relaxation` relaxes the mesh after each step.
  `loading=reversibility` relaxes forward, reverses the exact step to check
  reversibility, then continues from the saved forward state.
  `loading=noRelaxation` applies the deformation and reconnection without
  minimization or noise. For simple shear it requires
  `affineInitialGuess=true`, the default.
- `deformation=simpleShear` applies shear. `deformation=directDisplacement`
  translates selected fixed nodes by the same displacement, while free nodes
  are left for relaxation. Its load amounts are distances; shear load amounts
  are strains.
- `boundary=periodic` uses a periodic cell with all nodes free.
  `boundary=fixed`, the default, fixes the outer border.
  `boundary=anchored` fixes the first row and column of an open mesh,
  as in the double-dislocation construction. `boundary=periodicAnchored`
  applies those same constraints to a periodic cell.
- `directions=horizontal` shears in x with a gradient in y; `vertical` shears
  in y with a gradient in x. A minus sign reverses an individual direction,
  for example `horizontal,-vertical`. A negative `loadIncrement` reverses
  every direction in the sequence.
- `reconnectionMethod=none`, `edgeFlip`, or `delaunay` chooses reconnection.
  Delaunay requires a CGAL-enabled build.
- `minimum_transport_distance=true` chooses minimum flip-distance inheritance;
  `false` chooses maximum. This affects donor inheritance during edge flips.

For a direction sequence, set `segmentLoad` to the amount of loading per
segment. For example, these settings alternate horizontal and vertical shear
in segments of 0.2 until the `totalLoad` budget is used:

```ini
directions=horizontal,vertical
segmentLoad=0.2
repeat=true
```

With `repeat=false` (the default), the last direction uses the remaining budget.
For a single direction, omit `segmentLoad`: its default of zero means no
switching. Steps are shortened to reach segment boundaries and use the exact
remaining budget. `directions=horizontal,-horizontal` with `repeat=true` gives cyclic
shear.

`totalLoad` is the cumulative loading **after initialization**, regardless of
sign or direction. Its default is 1. `startLoad` initializes the deformation
at a position along the loading path without consuming this budget. For
example, `startLoad=0.15` and `totalLoad=0.85` finish at shear 1 for positive
horizontal simple shear. Both values are nonnegative; zero `totalLoad` records
only the initialized state. There is no `maxLoad` setting.

The saved `load` is the path coordinate, starting at `startLoad`. The CSV also
records `accumulated_load = load - startLoad`, starting at zero. A restart
continues counting toward the original total budget. The direction sequence
is anchored at path coordinate zero, so a nonzero `startLoad` can initialize
partway through a segment or cycle.

The initial state is relaxed with seeded `initialGuessNoise` unless loading
is `noRelaxation`. For simple shear, `mesh.currentDeformation` records the
imposed affine deformation; a local direct displacement cannot be described
by one global matrix.

A double-dislocation-style recipe uses the same runner:

```ini
loading=relaxation
deformation=directDisplacement
boundary=anchored
directions=horizontal,vertical
segmentLoad=1
repeat=false
startLoad=0
totalLoad=2
loadIncrement=0.002
reconnectionMethod=edgeFlip
```

Horizontal direct displacement moves fixed nodes in the upper half of the
original grid in x. Vertical displacement moves fixed nodes in the right half
in y. The driven node groups are selected once, saved in dumps, and reused
when reversing steps. `directDisplacement` requires `fixed`, `anchored`, or
`periodicAnchored` boundaries. To design another displacement region, edit the
selection in `prepareSimulation()` in the protocol file.

Boundary constraints, reference geometry, and driven node groups are restored
from the saved state. The new dump schema requires `totalLoad`, the saved
groups, and each element's accumulated plastic history; dumps from before this change are not supported.

Optional settings replace the remaining specialized loops:

| Setting | Default | Purpose |
| --- | --- | --- |
| `stepNoise` | `0` | Add noise before each relaxation, independently of initial noise. |
| `reversibilityTolerance` | `0.0001` | Displacement tolerance for the reversibility check. |
| `checkpointInterval` | `0` | Save a checkpoint before every Nth load step; zero disables it. Two files alternate. |
| `stopOnEnergyDrop` | `false` | Save an `energyDrop` dump and stop when average energy falls by more than `energyDropThreshold` between steps. |

An energy decrease during deliberate unloading can also trigger the stop.
Checkpoints can be enabled independently to keep states before the decrease.

The former simple-shear, fixed-boundary, noisy, cyclic, and energy-drop
experiments now use this shared runner. All experiment entries containing
`Test` have been removed. General parameters `GP1`, `GP2`, and `GP3` and old
config selector keys (`experiment`, `scenario`, `protocol`, and `usingPBC`)
are no longer accepted. Missing settings use the same defaults as
`Config::setDefaultValues()`; unknown keys cause an error.

The implementation is in [`src/Simulation/protocols.cpp`](src/Simulation/protocols.cpp).
Its opening comment gives pseudocode and explains where new deformation
patterns, loading behaviors, and boundaries belong. The numerical solvers and
reconnection methods remain in `Simulation` and `Mesh`.

### Deformation fields

Each element stores its accumulated plastic deformation directly as `T_p`.
The current element's elastic and plastic factors remain `F_E` and `F_P`:

```text
F = F_E F_P
T = F_E T_p
```

Between reconnections, a change of plastic well updates the stored history as
`T_p_new = F_P_new * inverse(F_P_old) * T_p_old`. An edge flip copies the donor's
`T_p` to its inheriting element (Rule I in the thesis); the new geometry supplies
`F_E`. Setting a new reference configuration resets the history.

VTU exports use `T_p11`, `T_p12`, `T_p21`, `T_p22` for accumulated plastic
deformation and `T11`, `T12`, `T21`, `T22` for history-corrected total deformation.
Matrix-history JSON and reports use these same meanings. No separate `H` matrix
is stored or exported, and old report-field aliases have been removed.

Dumps store complete state and require the current field names. Old reference
position names, plastic-counter aliases, and missing-field fallbacks have been
removed. Text configuration files still use `Config::setDefaultValues()` for
omitted options; removing old dump support does not change those defaults.

Non-improving reconnection candidates are a normal stopping condition. They
are restored without writing VTU files during ordinary runs. When
`logDuringMinimization=true`, the rejected candidate is retained as a
`rejectedReconnect` VTU pair for diagnosis. The regular CSV records
`nr_reconnect_cycles`, `reconnect_stop_reason`, and
`rejected_reconnect_energy_delta` without requiring VTU logging.

### Python Config Generator (Optional)

There is a Python config generator in `SimulationScripts/Management/configGenerator.py`. VSCode task `generateDefaultSettings` runs it.

There is also a complementary set of Python scripts here:
[https://github.com/EliasL/MTMSimulationScripts](https://github.com/EliasL/MTMSimulationScripts)
