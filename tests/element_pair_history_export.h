#ifndef ELEMENT_PAIR_HISTORY_EXPORT_H
#define ELEMENT_PAIR_HISTORY_EXPORT_H

#include "Eigen/Core"
#include "Simulation/simulation.h"
#include <array>
#include <filesystem>
#include <string>
#include <vector>

struct ElementPairMatrixHistoryRow {
  double gamma = 0.0;
  int loadStep = 0;
  std::array<Matrix2d, 2> F = {Matrix2d::Identity(), Matrix2d::Identity()};
  std::array<Matrix2d, 2> F_E = {Matrix2d::Identity(), Matrix2d::Identity()};
  std::array<Matrix2d, 2> F_P = {Matrix2d::Identity(), Matrix2d::Identity()};
  std::array<Matrix2d, 2> H = {Matrix2d::Identity(), Matrix2d::Identity()};
  std::array<Matrix2d, 2> T = {Matrix2d::Identity(), Matrix2d::Identity()};
  std::array<std::array<Vector2d, 3>, 2> currentNodes = {};
  std::array<std::array<Vector2d, 3>, 2> referenceNodes = {};
};

struct ElementPairStepReconnectSnapshot {
  double gamma = 0.0;
  int loadStep = 0;
  std::array<bool, 2> elementChanged = {false, false};
  ElementPairMatrixHistoryRow before;
  ElementPairMatrixHistoryRow after;
};

struct ElementPairStepReconnectLoggerContext {
  struct Target {
    std::array<int, 2> elementIndices = {0, 1};
    std::vector<ElementPairStepReconnectSnapshot> *rows = nullptr;
    bool hasPending = false;
    ElementPairStepReconnectSnapshot pending;
  };
  std::vector<Target> targets;
};

struct ElementPairStepReconnectJsonSpec {
  std::array<int, 2> elementIndices = {0, 1};
  const std::vector<ElementPairStepReconnectSnapshot> *rows = nullptr;
  std::string label;
};

ElementPairMatrixHistoryRow
captureElementPairMatrixHistoryRow(const Mesh &mesh,
                                   const std::array<int, 2> &elementIndices);

ElementPairMatrixHistoryRow
captureElementPairMatrixHistoryRow(Simulation &simulation,
                                   const std::array<int, 2> &elementIndices);

void recordElementPairStepReconnectSnapshot(Simulation &simulation,
                                            ReconnectStepStage stage,
                                            void *context);

void writeElementPairStepReconnectJson(
    const std::filesystem::path &path,
    const std::vector<ElementPairStepReconnectJsonSpec> &tables,
    const Config &config);

#endif
