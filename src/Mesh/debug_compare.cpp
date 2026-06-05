#include "debug_compare.h"
#include "compare_macros.h"
#include <string>

bool compareNodeIdsInternal(const NodeId &lhs, const NodeId &rhs,
                            std::string *debugMsg = nullptr,
                            int tabNumber = 0) {
  bool equal = true;
  COMPARE_FIELD(i);
  COMPARE_FIELD(idPos);
  COMPARE_FIELD(cols);
  return equal;
}

bool operator==(const NodeId &lhs, const NodeId &rhs) {
  return compareNodeIdsInternal(lhs, rhs, nullptr);
}

bool operator!=(const NodeId &lhs, const NodeId &rhs) { return !(lhs == rhs); }

bool compareNodesInternal(const Node &lhs, const Node &rhs,
                          std::string *debugMsg, int tabNumber) {
  bool equal = true;

  COMPARE_FIELD(id);
  COMPARE_FIELD(f);
  COMPARE_FIELD(fixedNode);
  COMPARE_FIELD(connectedElements);
  COMPARE_FIELD(nodeIndexInElement);
  COMPARE_FIELD(elementCount);
  COMPARE_FIELD(m_pos);
  COMPARE_FIELD(m_ref_pos);
  COMPARE_FIELD(m_u);

  return equal;
}

bool compareNodesInternal(const GhostNode &lhs, const GhostNode &rhs,
                          std::string *debugMsg = nullptr,
                          int tabNumber = 0) {
  bool equal = true;

  COMPARE_FIELD(referenceId);
  COMPARE_FIELD(id);
  COMPARE_FIELD(f);
  COMPARE_FIELD(periodicShift);
  COMPARE_FIELD(pos);
  COMPARE_FIELD(ref_pos);
  COMPARE_FIELD(u);

  return equal;
}

bool operator==(const Node &lhs, const Node &rhs) {
  return compareNodesInternal(lhs, rhs, nullptr, 0);
}

bool operator!=(const Node &lhs, const Node &rhs) { return !(lhs == rhs); }

bool operator==(const GhostNode &lhs, const GhostNode &rhs) {
  return compareNodesInternal(lhs, rhs, nullptr, 0);
}

bool operator!=(const GhostNode &lhs, const GhostNode &rhs) {
  return !(lhs == rhs);
}

std::string debugCompare(const Node &lhs, const Node &rhs, int tabNumber) {
  std::string diff;
  compareNodesInternal(lhs, rhs, &diff, tabNumber);
  return diff;
}

std::string debugCompare(const GhostNode &lhs, const GhostNode &rhs,
                         int tabNumber) {
  std::string diff;
  compareNodesInternal(lhs, rhs, &diff, tabNumber);
  return diff;
}

bool compareTElementsInternal(const TElement &lhs, const TElement &rhs,
                              std::string *debugMsg, int tabNumber) {
  bool equal = true;

  COMPARE_FIELD(ghostNodes);
  COMPARE_FIELD(C);
  COMPARE_FIELD(C_R);
  COMPARE_FIELD(M_l);
  COMPARE_FIELD(S);
  COMPARE_FIELD(P);
  COMPARE_FIELD(energy);
  COMPARE_FIELD(dN_dX);
  COMPARE_FIELD(m3Nr);
  COMPARE_FIELD(pastM3Nr);
  COMPARE_FIELD(pastStepM3Nr);
  COMPARE_FIELD(red_quadrant);
  COMPARE_FIELD(eIndex);
  COMPARE_FIELD(noise);
  COMPARE_FIELD(largestAngle);
  COMPARE_FIELD(angleNode);
  COMPARE_FIELD(initArea);
  COMPARE_FIELD(beta);
  COMPARE_FIELD(K);

  return equal;
}

bool operator==(const TElement &lhs, const TElement &rhs) {
  return compareTElementsInternal(lhs, rhs, nullptr, 0);
}

bool operator!=(const TElement &lhs, const TElement &rhs) {
  return !(lhs == rhs);
}

std::string debugCompare(const TElement &lhs, const TElement &rhs,
                         int tabNumber) {
  std::string diff;

  for (size_t i = 0; i < lhs.ghostNodes.size(); i++) {
    if (!(lhs.ghostNodes[i] == rhs.ghostNodes[i])) {
      diff += std::string(tabNumber, '\t') + "tElementNodes[" +
              std::to_string(i) + "] differs -> \n";
      diff += debugCompare(lhs.ghostNodes[i], rhs.ghostNodes[i], tabNumber + 1);
    }
  }
  compareTElementsInternal(lhs, rhs, &diff, tabNumber);
  return diff;
}

bool compareconnectesInternal(const Mesh &lhs, const Mesh &rhs,
                              std::string *debugMsg = nullptr,
                              int tabNumber = 0) {
  bool equal = true;

  COMPARE_FIELD(nodes);
  COMPARE_FIELD(elements);
  COMPARE_FIELD(F_P_H);
  COMPARE_FIELD(fixedNodeIds);
  COMPARE_FIELD(freeNodeIds);
  COMPARE_FIELD(a);
  COMPARE_FIELD(rows);
  COMPARE_FIELD(cols);
  COMPARE_FIELD(load);
  COMPARE_FIELD(loadSteps);
  COMPARE_FIELD(currentDeformation);
  COMPARE_FIELD(referenceDeformation);
  COMPARE_FIELD(nrElements);
  COMPARE_FIELD(nrNodes);
  COMPARE_FIELD(totalEnergy);
  COMPARE_FIELD(averageEnergy);
  COMPARE_FIELD(maxEnergy);
  COMPARE_FIELD(maxForce);
  COMPARE_FIELD(averageP11);
  COMPARE_FIELD(averageP12);
  COMPARE_FIELD(averageP21);
  COMPARE_FIELD(averageP22);
  COMPARE_FIELD(averageSigma11);
  COMPARE_FIELD(averageSigma12);
  COMPARE_FIELD(averageSigma22);
  COMPARE_FIELD(averageSigmaTrace);
  COMPARE_FIELD(maxM3Nr);
  COMPARE_FIELD(sumM3Nr);
  COMPARE_FIELD(maxPlasticJump);
  COMPARE_FIELD(minPlasticJump);
  COMPARE_FIELD(redQuadrantCounts);
  COMPARE_FIELD(redQuadrantFixedCounts);
  COMPARE_FIELD(nr_elements_with_m3_change);
  COMPARE_FIELD(nr_elements_with_m3_changeInStep);
  COMPARE_FIELD(QDSD);
  COMPARE_FIELD(energyFunction);
  COMPARE_FIELD(bulkModulus);
  COMPARE_FIELD(usingPBC);
  COMPARE_FIELD(simName);
  COMPARE_FIELD(dataPath);
  COMPARE_FIELD(bounds);

  return equal;
}

bool operator==(const Mesh &lhs, const Mesh &rhs) {
  return compareconnectesInternal(lhs, rhs, nullptr, 0);
}

bool operator!=(const Mesh &lhs, const Mesh &rhs) { return !(lhs == rhs); }

static std::string limitDebugLines(const std::string &text, size_t maxLines) {
  if (maxLines == 0 || text.empty()) {
    return text;
  }

  size_t lines = 0;
  size_t pos = 0;
  while (pos < text.size()) {
    size_t lineEnd = text.find('\n', pos);
    lines++;
    if (lines == maxLines) {
      if (lineEnd == std::string::npos) {
        return text;
      }
      if (lineEnd + 1 >= text.size()) {
        return text;
      }
      std::string out = text.substr(0, lineEnd);
      out += " ... (truncated after " + std::to_string(maxLines) + " lines)";
      return out;
    }
    if (lineEnd == std::string::npos) {
      return text;
    }
    pos = lineEnd + 1;
  }

  return text;
}

std::string debugCompare(const Mesh &lhs, const Mesh &rhs, size_t maxLines) {
  std::string diff;

  if (lhs.elements.size() != rhs.elements.size()) {
    diff += "elements size differs; ";
  } else {
    for (size_t i = 0; i < lhs.elements.size(); i++) {
      if (!(lhs.elements[i] == rhs.elements[i])) {
        diff += "elements[" + std::to_string(i) + "] differs -> \n";
        diff += debugCompare(lhs.elements[i], rhs.elements[i], 1);
      }
    }
  }
  if (lhs.nodes.size() != rhs.nodes.size()) {
    diff += "nodes size differs; ";
  } else {
    for (size_t i = 0; i < lhs.nodes.size(); i++) {
      if (!(lhs.nodes(i) == rhs.nodes(i))) {
        diff += "nodes[" + std::to_string(i) + "] differs -> \n";
        diff += debugCompare(lhs.nodes(i), rhs.nodes(i), 1);
      }
    }
  }

  compareconnectesInternal(lhs, rhs, &diff, 0);
  return limitDebugLines(diff, maxLines);
}
