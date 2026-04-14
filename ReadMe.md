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

Optional:
- Python 3 (for config generation and plotting scripts)
- CGAL dependencies (CMake will tell you if anything is missing)

Notes:
- ALGLIB, Eigen, and Cereal are included in `libs/` and built by CMake.
- CGAL is not fetched automatically. Place it in `libs/cgal` or update `CGAL_DIR` in `CMakeLists.txt`.
- CGAL is currently optional (used only for Delaunay reconnection). You can build without it using `-DIDE_LIGHTWEIGHT=ON`.

### Installation

1. **Clone the repository**

   ```sh
   git clone https://github.com/EliasL/MTS2D.git
   cd MTS2D
   ```

2. **CGAL setup**

   Place CGAL in `libs/cgal`, or update the `CGAL_DIR` path in `CMakeLists.txt`.

### Build

Using Visual Studio Code (recommended):
- Open the folder in VSCode.
- Use the built-in tasks in `.vscode/tasks.json` such as:
  - `build`
  - `build-release`
  - `test`
  - `generateDefaultSettings`

From a terminal:
- Debug:
  ```sh
  mkdir -p build && cd build && cmake .. && make
  ```
- Release:
  ```sh
  mkdir -p build-release && cd build-release && cmake -DCMAKE_BUILD_TYPE=Release .. && make
  ```

### Run

Basic usage:
```sh
./build-release/MTS2D -c path/to/config.conf -o path/to/output
```

Resume from a dump:
```sh
./build-release/MTS2D -d path/to/dump.xml.gz -o path/to/output
```

Common flags:
- `-c` config file
- `-d` dump file (resume)
- `-o` output path
- `-r` force re-run even if output folder looks complete

### Example Config

Minimal LBFGS example:
```ini
name=smallSimulation
rows=10
cols=10
usingPBC=true
reconnectionMethod=none
scenario=simpleShear
nrThreads=1
seed=0
QDSD=0.0
initialGuessNoise=0.05
meshDiagonal=major
energyFunction=contiSquare
bulkModulus=4.0

startLoad=0.15
loadIncrement=1e-5
maxLoad=0.151

minimizer=LBFGS
epsR=1e-5

LBFGSNrCorrections=10
LBFGSScale=1.0
LBFGSEpsg=1e-8
LBFGSEpsf=0
LBFGSEpsx=0
LBFGSMaxIterations=0

logDuringMinimization=false
writeDumps=false
plasticityEventThreshold=0.01
energyDropThreshold=0.001
showProgress=1
```

### Python Config Generator (Optional)

There is a Python config generator in `SimulationScripts/Management/configGenerator.py`. VSCode task `generateDefaultSettings` runs it.

There is also a complementary set of Python scripts here:
```
https://github.com/EliasL/MTMSimulationScripts
```
