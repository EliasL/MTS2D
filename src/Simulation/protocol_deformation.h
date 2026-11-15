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
void applyPreviewShear(Mesh &mesh, const Config &config,
                       const std::vector<SignedDirection> &directions, double endLoad);
bool isAxialDeformation(const Config &config);
void validateRectangleBoundary(const Config &config);
void validateDrivenEdgeConstraint(const Config &config);
void validateResumedBoundary(const Config &saved, const Config &requested);
double axialStretchAtLoad(const Config &config, double load);

#endif
