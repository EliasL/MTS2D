#ifndef PROTOCOLS_H
#define PROTOCOLS_H

#include "simulation.h"
#include <memory>
#include <string>

void validateProtocolConfig(const Config &config);
Matrix2d protocolDeformation(const Config &config, double load);
void initializeProtocolDeformation(Simulation &simulation);
void applyDeformationToMesh(Mesh &mesh, const Config &config,
                            const DeformationStep &step,
                            const DrivenNodeGroups &drivenNodes);
std::shared_ptr<Simulation> runSimulationProtocol(
    Config config, std::string dataPath,
    std::shared_ptr<Simulation> loadedSimulation = nullptr,
    double makeDumpAt = -1.0);

#endif
