#include "element_pair_history_export.h"
#include "Mesh/tElement.h"
#include <cmath>
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

bool elementPairEntryChanged(const ElementPairMatrixHistoryRow &before,
                             const ElementPairMatrixHistoryRow &after,
                             size_t elementIndex,
                             double tol = kReportElementChangeTolerance) {
  if (!matricesApproxEqual(before.F[elementIndex], after.F[elementIndex], tol) ||
      !matricesApproxEqual(before.F_E[elementIndex], after.F_E[elementIndex],
                           tol) ||
      !matricesApproxEqual(before.F_P[elementIndex], after.F_P[elementIndex],
                           tol) ||
      !matricesApproxEqual(before.H[elementIndex], after.H[elementIndex], tol) ||
      !matricesApproxEqual(before.T[elementIndex], after.T[elementIndex], tol)) {
    return true;
  }
  for (size_t nodeIndex = 0; nodeIndex < 3; ++nodeIndex) {
    if (!vectorsApproxEqual(before.currentNodes[elementIndex][nodeIndex],
                            after.currentNodes[elementIndex][nodeIndex], tol) ||
        !vectorsApproxEqual(before.referenceNodes[elementIndex][nodeIndex],
                            after.referenceNodes[elementIndex][nodeIndex],
                            tol)) {
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
  if (!std::isfinite(value)) {
    throw std::runtime_error(
        "writeElementPairStepReconnectJson: encountered non-finite value.");
  }
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

void writeJsonNodes(std::ostream &out,
                    const std::array<std::array<Vector2d, 3>, 2> &nodes) {
  out << "[";
  for (size_t elementIndex = 0; elementIndex < nodes.size(); ++elementIndex) {
    if (elementIndex > 0) {
      out << ",";
    }
    out << "[";
    for (size_t nodeIndex = 0; nodeIndex < nodes[elementIndex].size();
         ++nodeIndex) {
      if (nodeIndex > 0) {
        out << ",";
      }
      writeJsonVector2d(out, nodes[elementIndex][nodeIndex]);
    }
    out << "]";
  }
  out << "]";
}

void writeJsonMatrixPair(std::ostream &out,
                         const std::array<Matrix2d, 2> &matrices) {
  out << "[";
  for (size_t i = 0; i < matrices.size(); ++i) {
    if (i > 0) {
      out << ",";
    }
    writeJsonMatrix2d(out, matrices[i]);
  }
  out << "]";
}

void writeJsonHistoryRow(std::ostream &out,
                         const ElementPairMatrixHistoryRow &row) {
  out << "{";
  out << "\"gamma\":";
  writeJsonDouble(out, row.gamma);
  out << ",\"load_step\":" << row.loadStep;
  out << ",\"F\":";
  writeJsonMatrixPair(out, row.F);
  out << ",\"F_E\":";
  writeJsonMatrixPair(out, row.F_E);
  out << ",\"F_P\":";
  writeJsonMatrixPair(out, row.F_P);
  out << ",\"H\":";
  writeJsonMatrixPair(out, row.H);
  out << ",\"T\":";
  writeJsonMatrixPair(out, row.T);
  out << ",\"current_nodes\":";
  writeJsonNodes(out, row.currentNodes);
  out << ",\"reference_nodes\":";
  writeJsonNodes(out, row.referenceNodes);
  out << "}";
}

void writeJsonBoolPair(std::ostream &out,
                       const std::array<bool, 2> &values) {
  out << "[";
  out << (values[0] ? "true" : "false");
  out << ",";
  out << (values[1] ? "true" : "false");
  out << "]";
}

void writeJsonSnapshot(std::ostream &out,
                       const ElementPairStepReconnectSnapshot &snapshot) {
  out << "{";
  out << "\"gamma\":";
  writeJsonDouble(out, snapshot.gamma);
  out << ",\"load_step\":" << snapshot.loadStep;
  out << ",\"element_changed\":";
  writeJsonBoolPair(out, snapshot.elementChanged);
  out << ",\"before\":";
  writeJsonHistoryRow(out, snapshot.before);
  out << ",\"after\":";
  writeJsonHistoryRow(out, snapshot.after);
  out << "}";
}

} // namespace

ElementPairMatrixHistoryRow
captureElementPairMatrixHistoryRow(const Mesh &mesh,
                                   const std::array<int, 2> &elementIndices) {
  if (mesh.F_P_H.size() != mesh.elements.size()) {
    throw std::runtime_error(
        "captureElementPairMatrixHistoryRow: branch-history size mismatch.");
  }

  ElementPairMatrixHistoryRow row;
  row.gamma = mesh.load;
  row.loadStep = mesh.loadSteps;

  for (size_t i = 0; i < elementIndices.size(); ++i) {
    const int elementIndex = elementIndices[i];
    if (elementIndex < 0 ||
        elementIndex >= static_cast<int>(mesh.elements.size())) {
      throw std::runtime_error(
          "captureElementPairMatrixHistoryRow: requested element index out of "
          "range.");
    }
    const TElement &element = mesh.elements[static_cast<size_t>(elementIndex)];
    if (element.eIndex != elementIndex) {
      throw std::runtime_error(
          "captureElementPairMatrixHistoryRow: element index/eIndex mismatch.");
    }
    const Matrix2d &H = mesh.F_P_H.at(static_cast<size_t>(elementIndex));
    row.F[i] = element.F;
    row.F_E[i] = element.F_E;
    row.F_P[i] = element.F_P;
    row.H[i] = H;
    row.T[i] = element.totalBranch(H);
    for (int nodeIndex = 0; nodeIndex < 3; ++nodeIndex) {
      row.currentNodes[i][static_cast<size_t>(nodeIndex)] =
          element.ghostNodes[static_cast<size_t>(nodeIndex)].pos;
      row.referenceNodes[i][static_cast<size_t>(nodeIndex)] =
          element.ghostNodes[static_cast<size_t>(nodeIndex)].ref_pos;
    }
  }

  return row;
}

ElementPairMatrixHistoryRow
captureElementPairMatrixHistoryRow(Simulation &simulation,
                                   const std::array<int, 2> &elementIndices) {
  simulation.mesh.ensureFull();
  return captureElementPairMatrixHistoryRow(simulation.mesh, elementIndices);
}

void recordElementPairStepReconnectSnapshot(Simulation &simulation,
                                            ReconnectStepStage stage,
                                            void *context) {
  auto *loggerContext =
      static_cast<ElementPairStepReconnectLoggerContext *>(context);
  if (loggerContext == nullptr) {
    throw std::runtime_error(
        "recordElementPairStepReconnectSnapshot: missing logger context.");
  }

  for (ElementPairStepReconnectLoggerContext::Target &target :
       loggerContext->targets) {
    if (target.rows == nullptr) {
      throw std::runtime_error(
          "recordElementPairStepReconnectSnapshot: missing target row "
          "storage.");
    }

    if (stage == ReconnectStepStage::BeforeReconnect) {
      if (target.hasPending) {
        throw std::runtime_error(
            "recordElementPairStepReconnectSnapshot: before-reconnect "
            "snapshot arrived while another step was pending.");
      }
      target.pending = ElementPairStepReconnectSnapshot();
      target.pending.gamma = simulation.mesh.load;
      target.pending.loadStep = simulation.mesh.loadSteps;
      target.pending.before =
          captureElementPairMatrixHistoryRow(simulation, target.elementIndices);
      target.hasPending = true;
      continue;
    }

    if (!target.hasPending) {
      throw std::runtime_error(
          "recordElementPairStepReconnectSnapshot: after-reconnect snapshot "
          "arrived without a matching before-reconnect snapshot.");
    }

    target.pending.after =
        captureElementPairMatrixHistoryRow(simulation, target.elementIndices);
    if (std::abs(target.pending.gamma - simulation.mesh.load) > 1e-10 ||
        target.pending.loadStep != simulation.mesh.loadSteps) {
      throw std::runtime_error(
          "recordElementPairStepReconnectSnapshot: load/loadStep mismatch "
          "between before and after reconnect snapshots.");
    }
    for (size_t i = 0; i < target.pending.elementChanged.size(); ++i) {
      target.pending.elementChanged[i] =
          elementPairEntryChanged(target.pending.before, target.pending.after,
                                  i);
    }
    target.rows->push_back(target.pending);
    target.hasPending = false;
  }
}

void writeElementPairStepReconnectJson(
    const std::filesystem::path &path,
    const std::vector<ElementPairStepReconnectJsonSpec> &tables,
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
  out << ",\"scenario\":";
  writeJsonEscaped(out, config.scenario);
  out << ",\"name\":";
  writeJsonEscaped(out, config.name);
  out << "},\"tables\":[";

  for (size_t tableIndex = 0; tableIndex < tables.size(); ++tableIndex) {
    if (tableIndex > 0) {
      out << ",";
    }
    const ElementPairStepReconnectJsonSpec &table = tables[tableIndex];
    if (table.rows == nullptr) {
      throw std::runtime_error(
          "writeElementPairStepReconnectJson: missing table row data.");
    }
    out << "{";
    out << "\"label\":";
    writeJsonEscaped(out, table.label);
    out << ",\"element_indices\":[";
    out << table.elementIndices[0] << "," << table.elementIndices[1] << "]";
    out << ",\"rows\":[";
    for (size_t rowIndex = 0; rowIndex < table.rows->size(); ++rowIndex) {
      if (rowIndex > 0) {
        out << ",";
      }
      writeJsonSnapshot(out, table.rows->at(rowIndex));
    }
    out << "]";
    out << "}";
  }

  out << "]}";
}
