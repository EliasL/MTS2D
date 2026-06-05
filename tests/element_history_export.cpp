#include "element_history_export.h"
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace {

constexpr double kReportElementChangeTolerance = 0.1;

bool vectorsApproxEqual(const Vector2d &lhs, const Vector2d &rhs,
                        double tol = kReportElementChangeTolerance) {
  return (lhs - rhs).norm() <= tol;
}

bool matricesApproxEqual(const Matrix2d &lhs, const Matrix2d &rhs,
                         double tol = kReportElementChangeTolerance) {
  return (lhs - rhs).cwiseAbs().maxCoeff() <= tol;
}

bool elementEntryChanged(const ElementMatrixHistoryRow &before,
                         const ElementMatrixHistoryRow &after,
                         double tol = kReportElementChangeTolerance) {
  if (!matricesApproxEqual(before.F, after.F, tol) ||
      !matricesApproxEqual(before.F_E, after.F_E, tol) ||
      !matricesApproxEqual(before.F_P, after.F_P, tol) ||
      !matricesApproxEqual(before.H, after.H, tol) ||
      !matricesApproxEqual(before.T, after.T, tol)) {
    return true;
  }
  for (size_t nodeIndex = 0; nodeIndex < before.currentNodes.size();
       ++nodeIndex) {
    if (!vectorsApproxEqual(before.currentNodes[nodeIndex],
                            after.currentNodes[nodeIndex], tol) ||
        !vectorsApproxEqual(before.referenceNodes[nodeIndex],
                            after.referenceNodes[nodeIndex], tol)) {
      return true;
    }
  }
  return false;
}

void writeJsonEscaped(std::ostream &out, const std::string &text) {
  out << '"';
  for (char c : text) {
    switch (c) {
    case '\\':
      out << "\\\\";
      break;
    case '"':
      out << "\\\"";
      break;
    case '\n':
      out << "\\n";
      break;
    case '\r':
      out << "\\r";
      break;
    case '\t':
      out << "\\t";
      break;
    default:
      out << c;
      break;
    }
  }
  out << '"';
}

void writeJsonDouble(std::ostream &out, double value) {
  out << std::setprecision(17) << value;
}

void writeJsonVector2d(std::ostream &out, const Vector2d &value) {
  out << "[";
  writeJsonDouble(out, value.x());
  out << ",";
  writeJsonDouble(out, value.y());
  out << "]";
}

void writeJsonMatrix2d(std::ostream &out, const Matrix2d &matrix) {
  out << "[[";
  writeJsonDouble(out, matrix(0, 0));
  out << ",";
  writeJsonDouble(out, matrix(0, 1));
  out << "],[";
  writeJsonDouble(out, matrix(1, 0));
  out << ",";
  writeJsonDouble(out, matrix(1, 1));
  out << "]]";
}

void writeJsonNodes(std::ostream &out, const std::array<Vector2d, 3> &nodes) {
  out << "[";
  for (size_t nodeIndex = 0; nodeIndex < nodes.size(); ++nodeIndex) {
    if (nodeIndex > 0) {
      out << ",";
    }
    writeJsonVector2d(out, nodes[nodeIndex]);
  }
  out << "]";
}

void writeJsonHistoryRow(std::ostream &out, const ElementMatrixHistoryRow &row) {
  out << "{";
  out << "\"gamma\":";
  writeJsonDouble(out, row.gamma);
  out << ",\"load_step\":" << row.loadStep;
  out << ",\"F\":";
  writeJsonMatrix2d(out, row.F);
  out << ",\"F_E\":";
  writeJsonMatrix2d(out, row.F_E);
  out << ",\"F_P\":";
  writeJsonMatrix2d(out, row.F_P);
  out << ",\"H\":";
  writeJsonMatrix2d(out, row.H);
  out << ",\"T\":";
  writeJsonMatrix2d(out, row.T);
  out << ",\"current_nodes\":";
  writeJsonNodes(out, row.currentNodes);
  out << ",\"reference_nodes\":";
  writeJsonNodes(out, row.referenceNodes);
  out << "}";
}

void writeJsonSnapshot(std::ostream &out,
                       const ElementStepReconnectSnapshot &snapshot) {
  out << "{";
  out << "\"gamma\":";
  writeJsonDouble(out, snapshot.gamma);
  out << ",\"load_step\":" << snapshot.loadStep;
  out << ",\"element_changed\":"
      << (snapshot.elementChanged ? "true" : "false");
  out << ",\"before\":";
  writeJsonHistoryRow(out, snapshot.before);
  out << ",\"after\":";
  writeJsonHistoryRow(out, snapshot.after);
  out << "}";
}

} // namespace

ElementMatrixHistoryRow captureElementMatrixHistoryRow(const Mesh &mesh,
                                                       int elementIndex) {
  if (mesh.F_P_H.size() != mesh.elements.size()) {
    throw std::runtime_error(
        "captureElementMatrixHistoryRow: branch-history size mismatch.");
  }
  if (elementIndex < 0 ||
      elementIndex >= static_cast<int>(mesh.elements.size())) {
    throw std::runtime_error(
        "captureElementMatrixHistoryRow: requested element index out of "
        "range.");
  }

  const TElement &element = mesh.elements[static_cast<size_t>(elementIndex)];
  if (element.eIndex != elementIndex) {
    throw std::runtime_error(
        "captureElementMatrixHistoryRow: element index/eIndex mismatch.");
  }

  ElementMatrixHistoryRow row;
  row.gamma = mesh.load;
  row.loadStep = mesh.loadSteps;
  row.F = element.F;
  row.F_E = element.F_E;
  row.F_P = element.F_P;
  row.H = mesh.F_P_H.at(static_cast<size_t>(elementIndex));
  row.T = element.totalBranch(row.H);
  for (int nodeIndex = 0; nodeIndex < 3; ++nodeIndex) {
    row.currentNodes[static_cast<size_t>(nodeIndex)] =
        element.ghostNodes[static_cast<size_t>(nodeIndex)].pos;
    row.referenceNodes[static_cast<size_t>(nodeIndex)] =
        element.ghostNodes[static_cast<size_t>(nodeIndex)].ref_pos;
  }
  return row;
}

ElementMatrixHistoryRow captureElementMatrixHistoryRow(Simulation &simulation,
                                                       int elementIndex) {
  simulation.mesh.ensureFull();
  return captureElementMatrixHistoryRow(simulation.mesh, elementIndex);
}

void recordElementStepReconnectSnapshot(Simulation &simulation,
                                        ReconnectStepStage stage,
                                        void *context) {
  auto *loggerContext = static_cast<ElementStepReconnectLoggerContext *>(context);
  if (loggerContext == nullptr) {
    throw std::runtime_error(
        "recordElementStepReconnectSnapshot: missing logger context.");
  }
  if (loggerContext->rows == nullptr) {
    throw std::runtime_error(
        "recordElementStepReconnectSnapshot: missing row storage.");
  }

  if (stage == ReconnectStepStage::BeforeReconnect) {
    if (loggerContext->hasPending) {
      throw std::runtime_error(
          "recordElementStepReconnectSnapshot: before-reconnect snapshot "
          "arrived while another step was pending.");
    }
    loggerContext->pending = ElementStepReconnectSnapshot();
    loggerContext->pending.gamma = simulation.mesh.load;
    loggerContext->pending.loadStep = simulation.mesh.loadSteps;
    loggerContext->pending.before = captureElementMatrixHistoryRow(
        simulation, loggerContext->elementIndex);
    loggerContext->hasPending = true;
    return;
  }

  if (!loggerContext->hasPending) {
    throw std::runtime_error(
        "recordElementStepReconnectSnapshot: after-reconnect snapshot arrived "
        "without a matching before-reconnect snapshot.");
  }

  loggerContext->pending.after =
      captureElementMatrixHistoryRow(simulation, loggerContext->elementIndex);
  if (std::abs(loggerContext->pending.gamma - simulation.mesh.load) > 1e-10 ||
      loggerContext->pending.loadStep != simulation.mesh.loadSteps) {
    throw std::runtime_error(
        "recordElementStepReconnectSnapshot: load/loadStep mismatch between "
        "before and after reconnect snapshots.");
  }
  loggerContext->pending.elementChanged = elementEntryChanged(
      loggerContext->pending.before, loggerContext->pending.after);
  loggerContext->rows->push_back(loggerContext->pending);
  loggerContext->hasPending = false;
}

void writeElementStepReconnectJson(
    const std::filesystem::path &path, int elementIndex,
    const std::string &label, const std::vector<ElementStepReconnectSnapshot> &rows,
    const Config &config) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not open step reconnect JSON file: " +
                             path.string());
  }

  out << std::setprecision(17);
  out << "{";
  out << "\"config\":{";
  out << "\"rows\":" << config.rows;
  out << ",\"cols\":" << config.cols;
  out << ",\"load_increment\":";
  writeJsonDouble(out, config.loadIncrement);
  out << ",\"experiment\":";
  writeJsonEscaped(out, config.experiment);
  out << ",\"name\":";
  writeJsonEscaped(out, config.name);
  out << "},\"element_index\":" << elementIndex;
  out << ",\"label\":";
  writeJsonEscaped(out, label);
  out << ",\"rows\":[";
  for (size_t rowIndex = 0; rowIndex < rows.size(); ++rowIndex) {
    if (rowIndex > 0) {
      out << ",";
    }
    writeJsonSnapshot(out, rows[rowIndex]);
  }
  out << "]}";
}
