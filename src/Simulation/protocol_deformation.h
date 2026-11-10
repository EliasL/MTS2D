#ifndef PROTOCOL_DEFORMATION_H
#define PROTOCOL_DEFORMATION_H

#include "protocols.h"
#include <vector>

struct SignedDirection {
  Direction axis;
  double sign;
  SlipRegion slip{};
};

std::vector<SignedDirection> directionSequence(const Config &config);
void prepareSimulation(Simulation &simulation, const std::vector<SignedDirection> &directions);
bool isAxialDeformation(const Config &config);
double axialStretchAtLoad(const Config &config, double load);

#endif
