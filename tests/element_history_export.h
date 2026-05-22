#ifndef ELEMENT_HISTORY_EXPORT_H
#define ELEMENT_HISTORY_EXPORT_H

#include "Eigen/Core"
#include "Simulation/simulation.h"
#include <array>
#include <filesystem>
#include <string>
#include <vector>

struct ElementMatrixHistoryRow {
  double gamma = 0.0;
  int loadStep = 0;
  Matrix2d F = Matrix2d::Identity();
  Matrix2d F_E = Matrix2d::Identity();
  Matrix2d F_P = Matrix2d::Identity();
  Matrix2d H = Matrix2d::Identity();
  Matrix2d T = Matrix2d::Identity();
  std::array<Vector2d, 3> currentNodes = {};
  std::array<Vector2d, 3> referenceNodes = {};
};

struct ElementStepReconnectSnapshot {
  double gamma = 0.0;
  int loadStep = 0;
  bool elementChanged = false;
  ElementMatrixHistoryRow before;
  ElementMatrixHistoryRow after;
};

struct ElementStepReconnectLoggerContext {
  int elementIndex = -1;
  std::vector<ElementStepReconnectSnapshot> *rows = nullptr;
  bool hasPending = false;
  ElementStepReconnectSnapshot pending;
};

ElementMatrixHistoryRow captureElementMatrixHistoryRow(const Mesh &mesh,
                                                       int elementIndex);

ElementMatrixHistoryRow captureElementMatrixHistoryRow(Simulation &simulation,
                                                       int elementIndex);

void recordElementStepReconnectSnapshot(Simulation &simulation,
                                        ReconnectStepStage stage,
                                        void *context);

void writeElementStepReconnectJson(const std::filesystem::path &path,
                                   int elementIndex, const std::string &label,
                                   const std::vector<ElementStepReconnectSnapshot> &rows,
                                   const Config &config);

#endif
