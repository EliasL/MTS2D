# Mesoscopic Tensorial Model (MTM) Simulation

## About The Project

This project presents a novel simulation using the Mesoscopic Tensorial Model (MTM) to study dislocation nucleation in defect-free crystals. This simulation aims to explore the formation and behavioure of dislocation nucleation patterns.

## Getting Started

Instructions to set up your project locally. To get a local copy up and running, follow these simple steps.

### Prerequisites

This project requires the following tools and libraries:

- **CMake**: Essential for building the project. Download from the [official CMake website](https://cmake.org/download/).
- **ALGLIB**: A library for numerical analysis and data processing. More information can be found on the [ALGLIB website](https://www.alglib.net/download.php).
- **OpenMP**: For parallel computing. OpenMP is usually included with compilers like GCC or Clang. Check your compiler's documentation for OpenMP support.
- **VSCode**: Recommended for leveraging pre-configured build settings. Download VSCode from the [official VSCode website](https://code.visualstudio.com/Download).


### Installation

Follow these steps to set up the project on your local machine:

1. **Clone the Repository**

   Open a terminal and run the following command to clone the repo:

   ```sh
   git clone https://github.com/EliasL/1D-version1.git
   ```

2. **Using Visual Studio Code (Recommended)**

   It's recommended to use Visual Studio Code (VSCode) for this project to leverage the pre-configured build settings in the '.vscode' folder. 

   - If you haven't already, [download and install VSCode](https://code.visualstudio.com/Download).
   - Open the cloned project folder in VSCode.

3. **Alternative Setup without VSCode**

   If you prefer not to use VSCode, you can manually build the project using the following commands in your terminal:

   - For a Debug Build:
     ```sh
     mkdir -p build-debug && cd build-debug && cmake .. && make
     ```
   - For a Release Build:
     ```sh
     mkdir -p build-release && cd build-release && cmake -DCMAKE_BUILD_TYPE=Release .. && make
     ```

   These commands create a new directory for the build (either debug or release), configure the project with CMake, and then compile it with Make.