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
- CGAL is included in `libs/`, but some systems still require CGAL dependencies.

### Installation

1. **Clone the repository**

   ```sh
   git clone https://github.com/EliasL/MTS2D.git
   cd MTS2D
   ```

2. **Initialize submodules (safe and pinned)**

   This checks out the exact submodule commits referenced by the repo.
   ```sh
   git submodule update --init --recursive
   ```

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

### Python Config Generator (Optional)

There is a Python config generator in `SimulationScripts/Management/configGenerator.py`. VSCode task `generateDefaultSettings` runs it.

There is also a complementary set of Python scripts here:
```
https://github.com/EliasL/MTMSimulationScripts
```
