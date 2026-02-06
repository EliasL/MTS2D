#include "Simulation/scenarios.h"
#include "Simulation/simulation.h"

int main(int argc, char *argv[]) {
  try {
    // This starts the simulation.
    handleInputArgs(argc, argv);
    // If the simulation is already complete,
    // the program will terminate:
  } catch (const SimulationAlreadyComplete &e) {
    std::cout << e.what() << std::endl;
    return EXIT_SUCCESS;
  }
}
