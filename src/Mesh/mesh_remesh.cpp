#include "mesh.h"
#include "Data/logging.h"
#include "Mesh/node.h"
#include "Mesh/tElement.h"
#include "Simulation/randomUtils.h"
#include <Eigen/Core>
#include <Eigen/LU>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <vector>

using Eigen::Vector2i;
// CGAL is a heavy library and slows down clangd and other tools a lot.
// We use a seperate build without CGAL.
#if !defined(IDE_LIGHTWEIGHT)
#include <CGAL/Delaunay_triangulation_2.h>
// I had some issues with inexact constructions giving different triangulations
// where they should be the same.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/centroid.h>

using K = CGAL::Exact_predicates_exact_constructions_kernel; // exact predicates
                                                             // & constructions
using Point = K::Point_2;
using Delaunay = CGAL::Delaunay_triangulation_2<K>;
// Add this just above the CGAL aliases
struct VInfo {
  int refNodeIndex; // index of the reference node in the MTS mesh
                    // (0..nrNodes-1)
};
using Vb = CGAL::Triangulation_vertex_base_with_info_2<VInfo, K>;
using Fb = CGAL::Triangulation_face_base_2<K>;
using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using DelaunayInfo = CGAL::Delaunay_triangulation_2<K, Tds>;
#endif

inline bool inRegion(const Eigen::Matrix2d &G) {
  const double a = G(0, 0);
  const double b = G(0, 1);
  assert(G(0, 1) == G(1, 0));
  const double c = G(1, 1);
  if (b < 0) {
    return false;
  }
  if (b > std::min(a, c)) {
    return false;
  }
  return true;
}

inline bool inExtendedRegion(const Eigen::Matrix2d &G) {
  const double a = G(0, 0);
  const double b = 0.5 * (G(0, 1) + G(1, 0)); // robust off-diagonal
  const double c = G(1, 1);

  if (a <= c) {
    const double b_lo = -0.5 * a;
    const double b_hi = std::min(1.5 * a, (3.0 * a + c) / 4.0);
    return (b >= b_lo) && (b <= b_hi);
  } else {
    const double b_lo = -0.5 * c;
    const double b_hi = std::min(1.5 * c, (a + 3.0 * c) / 4.0);
    return (b >= b_lo) && (b <= b_hi);
  }
}

inline double maxCosineForMinAngle(const TElement &e) {
  // We want to compare the smalest angle, but
  // avoid expensive trigonometric functions.
  double maxCos = -1.0;
  for (int i = 0; i < 3; ++i) {
    const int next = (i + 1) % 3;
    const int prev = (i + 2) % 3;
    const Vector2d v1 = e.ghostNodes[next].pos - e.ghostNodes[i].pos;
    const Vector2d v2 = e.ghostNodes[prev].pos - e.ghostNodes[i].pos;
    const double denom = v1.norm() * v2.norm();
    if (denom <= 1e-12) {
      return 1.0;
    }
    double cosA = v1.dot(v2) / denom;
    cosA = std::clamp(cosA, -1.0, 1.0);
    if (cosA > maxCos) {
      maxCos = cosA;
    }
  }
  return maxCos;
}

static Mesh::EdgeKey getSharedEdgeKey(const TElement &e1, const TElement &e2) {
  std::array<int, 2> sharedNodeIds = {-1, -1};
  int sharedNodeCount = 0;
  for (const GhostNode &g1 : e1.ghostNodes) {
    for (const GhostNode &g2 : e2.ghostNodes) {
      if (g1.referenceId.i != g2.referenceId.i) {
        continue;
      }
      if (sharedNodeCount >= 2) {
        throw std::runtime_error(
            "getSharedEdgeKey: element pair shares more than two nodes.");
      }
      sharedNodeIds[sharedNodeCount++] = g1.referenceId.i;
      break;
    }
  }

  if (sharedNodeCount != 2) {
    throw std::runtime_error(
        "getSharedEdgeKey: element pair does not share a common edge.");
  }
  return Mesh::EdgeKey(sharedNodeIds[0], sharedNodeIds[1]);
}

static Mesh::EdgeKey getCurrentAngleEdgeKey(const TElement &e) {
  const auto coNodes = e.getCoNodesByIndex();
  return Mesh::EdgeKey(coNodes[0]->referenceId.i, coNodes[1]->referenceId.i);
}

// Stores the TElement::eIndex values of the elements using one edge.
struct EdgeTwinEntry {
  int firstElementIndex = -1;
  int secondElementIndex = -1;
};

using EdgeTwinLookup =
    std::unordered_map<Mesh::EdgeKey, EdgeTwinEntry, Mesh::EdgeKeyHash>;

static bool isValidElementIndex(const std::vector<TElement> &elements,
                                int elementIndex) {
  return elementIndex >= 0 && elementIndex < static_cast<int>(elements.size());
}

static const TElement *elementPtrOrNull(const std::vector<TElement> &elements,
                                        int elementIndex) {
  if (!isValidElementIndex(elements, elementIndex)) {
    return nullptr;
  }
  return &elements[static_cast<size_t>(elementIndex)];
}

static bool angleEdgesAreCompatible(const TElement &e1, const TElement &e2) {
  const auto e1Co = e1.getCoNodesByIndex();
  const auto e2Co = e2.getCoNodesByIndex();

  if (e1Co[0]->referenceId.i != e2Co[0]->referenceId.i ||
      e1Co[1]->referenceId.i != e2Co[1]->referenceId.i) {
    return false;
  }

  const Vector2i deltaShift0 = e1Co[0]->periodicShift - e2Co[0]->periodicShift;
  const Vector2i deltaShift1 = e1Co[1]->periodicShift - e2Co[1]->periodicShift;
  return (deltaShift0.array() == deltaShift1.array()).all();
}

static void addToEdgeTwinLookup(EdgeTwinLookup &lookup,
                                const std::vector<TElement> &elements,
                                const TElement &e) {
  const Mesh::EdgeKey angleEdge = getCurrentAngleEdgeKey(e);
  EdgeTwinEntry &entry = lookup[angleEdge];

  if (entry.firstElementIndex == e.eIndex ||
      entry.secondElementIndex == e.eIndex) {
    throw std::runtime_error(
        "addToEdgeTwinLookup: duplicate element in edge lookup.");
  }
  if (entry.firstElementIndex == -1) {
    entry.firstElementIndex = e.eIndex;
    return;
  }
  if (entry.secondElementIndex == -1) {
    entry.secondElementIndex = e.eIndex;
    return;
  }

  const TElement *firstElement =
      elementPtrOrNull(elements, entry.firstElementIndex);
  const TElement *secondElement =
      elementPtrOrNull(elements, entry.secondElementIndex);
  if (firstElement == nullptr || secondElement == nullptr) {
    throw std::runtime_error(
        "addToEdgeTwinLookup: edge lookup contains invalid element indices.");
  }

  // This is a rare case where the same edge appears more than twice in the
  // mesh. This can easily happen in a 2x2 mesh, but otherwise, it should not
  // and is probably a bug
  if (elements.size() == 8) {
    const bool keepExistingPair =
        angleEdgesAreCompatible(*firstElement, *secondElement);
    const bool keepFirstAndNew = angleEdgesAreCompatible(*firstElement, e);
    const bool keepSecondAndNew = angleEdgesAreCompatible(*secondElement, e);
    const int compatiblePairCount = static_cast<int>(keepExistingPair) +
                                    static_cast<int>(keepFirstAndNew) +
                                    static_cast<int>(keepSecondAndNew);

    if (compatiblePairCount == 1) {
      if (keepFirstAndNew) {
        entry.secondElementIndex = e.eIndex;
      } else if (keepSecondAndNew) {
        entry.firstElementIndex = entry.secondElementIndex;
        entry.secondElementIndex = e.eIndex;
      }
      return;
    }
  }

  throw std::runtime_error(DebugLog::formatEdgeTwinLookupOverflow(
      angleEdge.nodeIdA, angleEdge.nodeIdB, *firstElement, *secondElement, e));
}

static void removeFromEdgeTwinLookup(EdgeTwinLookup &lookup,
                                     const TElement &e) {
  const Mesh::EdgeKey angleEdge = getCurrentAngleEdgeKey(e);
  auto it = lookup.find(angleEdge);
  if (it == lookup.end()) {
    throw std::runtime_error(
        "removeFromEdgeTwinLookup: missing angle edge in lookup.");
  }

  EdgeTwinEntry &entry = it->second;
  if (entry.firstElementIndex == e.eIndex) {
    entry.firstElementIndex = entry.secondElementIndex;
    entry.secondElementIndex = -1;
  } else if (entry.secondElementIndex == e.eIndex) {
    entry.secondElementIndex = -1;
  } else {
    throw std::runtime_error(
        "removeFromEdgeTwinLookup: element not found on angle edge.");
  }

  if (entry.firstElementIndex == -1) {
    lookup.erase(it);
  }
}

static EdgeTwinLookup
buildEdgeTwinLookup(const std::vector<TElement> &elements) {
  EdgeTwinLookup lookup;
  lookup.reserve(elements.size());
  for (const TElement &e : elements) {
    addToEdgeTwinLookup(lookup, elements, e);
  }
  return lookup;
}

static int findTwinFromLookup(const EdgeTwinLookup &lookup,
                              const std::vector<TElement> &elements,
                              const TElement &e) {
  const Mesh::EdgeKey angleEdge = getCurrentAngleEdgeKey(e);
  auto it = lookup.find(angleEdge);
  if (it == lookup.end()) {
    return -1;
  }

  const EdgeTwinEntry &entry = it->second;
  if (entry.firstElementIndex == e.eIndex) {
    return entry.secondElementIndex;
  }
  if (entry.secondElementIndex == e.eIndex) {
    return entry.firstElementIndex;
  }

  const TElement *firstElement =
      elementPtrOrNull(elements, entry.firstElementIndex);
  const TElement *secondElement =
      elementPtrOrNull(elements, entry.secondElementIndex);
  if (firstElement != nullptr && secondElement != nullptr &&
      !angleEdgesAreCompatible(e, *firstElement) &&
      !angleEdgesAreCompatible(e, *secondElement)) {
    return -1;
  }

  throw std::runtime_error(DebugLog::formatEdgeTwinLookupMismatch(
      angleEdge.nodeIdA, angleEdge.nodeIdB, e, firstElement,
      entry.firstElementIndex, secondElement, entry.secondElementIndex));
}

static void toggleEdgeKey(Mesh::EdgeSet &edgeSet, const Mesh::EdgeKey &edge) {
  auto [it, inserted] = edgeSet.insert(edge);
  if (!inserted) {
    edgeSet.erase(it);
  }
}

bool Mesh::reconnect(bool onlyCheck, EdgeSet *lockedEdges) {
  ensureGeometry();
  meshReconnected = false;
  bool changedInAnySweep = false;
  EdgeTwinLookup edgeTwinLookup = buildEdgeTwinLookup(elements);

  while (true) {
    bool changedThisSweep = false;
    for (int i = 0; i < elements.size(); i++) {
      TElement &e = elements[i];

      // Is the element deformed enough to be reconnected?
      if (inRegion(e.G)) {
        continue;
      }
      // Does the element have a matching twin to reconnect with?
      const int twinIndex = findTwinFromLookup(edgeTwinLookup, elements, e);
      if (twinIndex == -1) {
        continue;
      }
      // Is the twin also deformed enough to be reconnected?
      TElement &twin = elements[twinIndex];
      if (inRegion(twin.G)) {
        continue;
      }

      // if (e.F_P != twin.F_P) {
      //   std::cout << "skipping different F_P\n";
      //   continue;
      // }

      // Is the edge being flipped locked?
      const EdgeKey sharedEdge = getSharedEdgeKey(e, twin);
      if (lockedEdges != nullptr &&
          lockedEdges->find(sharedEdge) != lockedEdges->end()) {
        continue;
      }

      if (onlyCheck) {
        return true;
      }
      removeFromEdgeTwinLookup(edgeTwinLookup, e);
      removeFromEdgeTwinLookup(edgeTwinLookup, twin);
      try {
        checkAndFixPeriodicElementPair(e, twin);
      } catch (const std::exception &ex) {
        throw std::runtime_error("Mesh::reconnect: could not align periodic "
                                 "twin before edge flip: " +
                                 std::string(ex.what()));
      } catch (...) {
        throw std::runtime_error(
            "Mesh::reconnect: could not align periodic twin before edge flip.");
      }

      flipEdge(e, twin);
      addToEdgeTwinLookup(edgeTwinLookup, elements, elements[i]);
      addToEdgeTwinLookup(edgeTwinLookup, elements, elements[twinIndex]);

      const EdgeKey newSharedEdge =
          getSharedEdgeKey(elements[i], elements[twinIndex]);
      toggleEdgeKey(edgeFlipDeltaSinceLastStep, sharedEdge);
      toggleEdgeKey(edgeFlipDeltaSinceLastStep, newSharedEdge);
      if (lockedEdges != nullptr) {
        lockedEdges->insert(newSharedEdge);
      }

      changedThisSweep = true;
      changedInAnySweep = true;
      meshReconnected = true;
    }

    if (onlyCheck || !changedThisSweep) {
      break;
    }
  }

  if (!onlyCheck && changedInAnySweep) {
    markDirty();
    nrMinItterationsSinceLastReconnect = 0;
    nrMinFunctionCallsSinceLastReconnect = 0;
  }
  return changedInAnySweep;
}

// Delauney Helpers

#if !defined(IDE_LIGHTWEIGHT)
inline Eigen::Vector2d toEigen(const Point &p) {
  return {CGAL::to_double(p.x()), CGAL::to_double(p.y())};
}
// Unique key for a triangle based on its three real node ids.
struct TriKey {
  uint64_t nodeIdA;
  uint64_t nodeIdB;
  uint64_t nodeIdC;
  bool operator==(const TriKey &o) const noexcept {
    return nodeIdA == o.nodeIdA && nodeIdB == o.nodeIdB && nodeIdC == o.nodeIdC;
  }
  // Define less-than operator for use in std::set
  bool operator<(const TriKey &o) const noexcept {
    if (nodeIdA != o.nodeIdA)
      return nodeIdA < o.nodeIdA;
    if (nodeIdB != o.nodeIdB)
      return nodeIdB < o.nodeIdB;
    return nodeIdC < o.nodeIdC;
  }
};

static inline TriKey makeTriKey(const DelaunayInfo::Face_handle &f) {
  std::array<uint64_t, 3> nodeIds{0, 0, 0};
  for (int k = 0; k < 3; ++k) {
    nodeIds[k] = f->vertex(k)->info().refNodeIndex;
  }
  std::sort(nodeIds.begin(), nodeIds.end());
  return TriKey{nodeIds[0], nodeIds[1], nodeIds[2]};
}

void Mesh::reconnectDelaunay() {
  // We will refer to the native MTS mesh as the "MTSMesh", and the Delaunay
  // triangulation as the "Dmesh".
  // In order to deal with periodic boundaries, we will extend the mesh
  // by three node layers in each direction.
  int extension = 3;
  int extendedRows = rows + extension * 2;
  int extendedCols = cols + extension * 2;
  int nrNodesInDMesh = extendedRows * extendedCols;
  // Note that the same NodeId will appear multiple times
  std::vector<NodeId> DNodeToMTSNode; // index in Dmesh -> NodeId in MTSMesh
  DNodeToMTSNode.reserve(nrNodesInDMesh);

  // 2) Build Delaunay with vertex info = VInfo (baseId, dr, dc)
  DelaunayInfo dt;
  for (int dr = -extension; dr < rows + extension; ++dr) {
    for (int dc = -extension; dc < cols + extension; ++dc) {
      const GhostNode gn = getGhostNode(Vector2i{dc, dr});
      const Point p(gn.pos.x(), gn.pos.y());
      auto vh = dt.insert(p);
      VInfo vi;
      vi.refNodeIndex = gn.referenceId.i;
      vh->info() = vi;
    }
  }

  // 3) Remove excess elements
  // We only keep faces with two or more vertices inside the window
  // And we also make sure to only keep unique triangles (no duplicates)
  std::set<TriKey> seenTriangles;

  auto keep_element = [&](const DelaunayInfo::Face_handle &f) {
    // Deduplicate by canonical triangle key based on reference vertex indices
    TriKey key = makeTriKey(f);
    if (seenTriangles.find(key) != seenTriangles.end()) {
      return false; // duplicate triangle
    }

    const Point &p0 = f->vertex(0)->point();
    const Point &p1 = f->vertex(1)->point();
    const Point &p2 = f->vertex(2)->point();

    // For some reason, CGAL makes flat triagnles sometimes.
    const K::FT area = CGAL::area(p0, p1, p2); // non-negative area in K::FT
    if (area < K::FT(1e-5)) {
      return false; // near-flat triangle
    }

    // Half-open base window test via exact comparisons on K::FT
    const Point c = CGAL::centroid(p0, p1, p2);
    const double cx = CGAL::to_double(c.x());
    const double cy = CGAL::to_double(c.y());

    // We keep an extra row and column if using PBC. These are the elements
    // that wrap around.
    int extra = usingPBC ? 1 : 0;
    // Keep iff 0 <= cx < Lx and 0 <= cy < Ly (half-open [0,L))
    if (!(0 <= cx && cx < (cols - 1 + extra) * a && //
          0 <= cy && cy < (rows - 1 + extra) * a)) {
      return false;
    }

    seenTriangles.insert(key);
    return true;
  };

  // Collect faces to keep
  std::vector<DelaunayInfo::Face_handle> kept_faces;
  kept_faces.reserve(
      std::distance(dt.finite_faces_begin(), dt.finite_faces_end()));
  for (auto f = dt.finite_faces_begin(); f != dt.finite_faces_end(); ++f) {
    if (keep_element(f)) {
      kept_faces.push_back(f);
    }
  }

  // 4) Sanity: number of faces should match our preallocated elements
  const int nFaces = static_cast<int>(kept_faces.size());
  if (nFaces != nrElements) {
    elements.resize(nFaces);
    F_P_H.resize(nFaces, Matrix2d::Identity());
    std::cerr
        << "Warning: reconnectDelaunay(): triangle count mismatch (expected "
        << nrElements << ", got " << nFaces << ")." << std::endl;
    nrElements = nFaces;
  }

  // 5) Reset per-node connectivity
  for (long i = 0; i < nodes.size(); ++i) {
    nodes(i).elementCount = 0;
  }

  // 6) Refill elements from CGAL faces (each face -> one TElement)
  int eIdx = 0;
  for (const auto &f : kept_faces) {
    std::array<GhostNode, 3> g;
    for (int k = 0; k < 3; ++k) {
      auto v = f->vertex(k);
      const VInfo &vi = v->info();
      const Node *n = &nodes(vi.refNodeIndex); // reference node in MTS mesh
      // Build ghost from DT vertex coordinate. Use the targetPos to set the
      // periodic shift.
      g[k] = m_gn(n, toEigen(v->point()));
    }

    double noise = (eIdx < (int)elements.size()) ? elements[eIdx].noise
                                                 : sampleNormal(1, QDSD);
    elements[eIdx] = TElement((*this), g[0], g[1], g[2], eIdx, noise,
                              energyFunction, bulkModulus);
    ++eIdx;
  }

  markDirty();
  nrMinItterationsSinceLastReconnect = 0;
  nrMinFunctionCallsSinceLastReconnect = 0;
}
#else
void Mesh::reconnectDelaunay() {
  throw std::runtime_error(
      "reconnectDelaunay requires CGAL; build without IDE_LIGHTWEIGHT.");
}
#endif

std::array<GhostNode, 4> Mesh::getElementPairNodes(const TElement &e1,
                                                   const TElement &e2) {
  auto e1Co = e1.getCoNodesByIndex();

  const GhostNode *el1AngleNode = e1.getAngleNode();
  const GhostNode *el2AngleNode = e2.getAngleNode();

  return {*el1AngleNode, *e1Co[0], *e1Co[1], *el2AngleNode};
}

void Mesh::removeElementsFromNodes(
    const std::array<const GhostNode *, 4> &gNodes, int e1i, int e2i) {
  const std::array<Node *, 4> nodes = {
      (*this)[gNodes[0]->referenceId], (*this)[gNodes[1]->referenceId],
      (*this)[gNodes[2]->referenceId], (*this)[gNodes[3]->referenceId]};
  removeElementsFromNodes(nodes, e1i, e2i);
}

void Mesh::removeElementsFromNodes(const std::array<Node *, 4> &nodes, int e1i,
                                   int e2i) {
  for (Node *n : nodes) {
    int tempElementIndices[MAX_ELEMENTS_PER_NODE];
    int tempNodeIndexInElement[MAX_ELEMENTS_PER_NODE];
    int newCount = 0;

    for (int i = 0; i < n->elementCount; i++) {
      const int el = n->connectedElements[i];
      if (el == e1i || el == e2i) {
        continue;
      }
      tempElementIndices[newCount] = n->connectedElements[i];
      tempNodeIndexInElement[newCount] = n->nodeIndexInElement[i];
      newCount++;
    }

    for (int i = 0; i < newCount; i++) {
      n->connectedElements[i] = tempElementIndices[i];
      n->nodeIndexInElement[i] = tempNodeIndexInElement[i];
    }
    for (int i = newCount; i < n->elementCount; i++) {
      n->connectedElements[i] = -1;
      n->nodeIndexInElement[i] = -1;
    }
    n->elementCount = newCount;
  }
}

void Mesh::flipEdge(TElement &e1, TElement &e2) {
  // This function takes two elements that should both have large angles, and
  // reconfigures the 4 nodes into two new elements that have smaller angles.

  const int e1i = e1.eIndex;
  const int e2i = e2.eIndex;
  if (e1i < 0 || e2i < 0 || e1i >= nrElements || e2i >= nrElements) {
    throw std::runtime_error("Mesh::flipEdge: invalid element indices.");
  }
  if (F_P_H.size() != elements.size()) {
    throw std::runtime_error(
        "Mesh::flipEdge: branch history size does not match elements.");
  }
  const auto oe1 = e1.getAngleCo1Co2Nodes();
  const auto oe2 = e2.getAngleCo1Co2Nodes();

  // The two reconnect options keep both old angle nodes and choose which
  // counterclockwise co-node stays with each new child element.
  // We either need to replace the first or the second node with the angle node
  // of the twin. We can compute the distances from the old position to the new
  // position as follows:

  double distA =
      (oe1[1].pos - oe2[0].pos).norm() + (oe2[1].pos - oe1[0].pos).norm();
  double distB =
      (oe1[2].pos - oe2[0].pos).norm() + (oe2[2].pos - oe1[0].pos).norm();
  const std::array<GhostNode, 3> n1a = {oe1[0], oe1[1], oe2[0]};
  const std::array<GhostNode, 3> n2a = {oe2[0], oe2[1], oe1[0]};
  const std::array<GhostNode, 3> n1b = {oe1[0], oe1[2], oe2[0]};
  const std::array<GhostNode, 3> n2b = {oe2[0], oe2[2], oe1[0]};

  const bool chooseA = distA <= distB;
  const auto &n1 = chooseA ? n1a : n1b;
  const auto &n2 = chooseA ? n2a : n2b;

  if (compareEdgeFlipOptions) {
    const TElement::EdgeFlipRemeshState state1 =
        e1.evaluateEdgeFlipRemeshState(n1, F_P_H[static_cast<size_t>(e1i)]);
    const TElement::EdgeFlipRemeshState state2 =
        e2.evaluateEdgeFlipRemeshState(n2, F_P_H[static_cast<size_t>(e2i)]);
    if (!state1.valid || !state2.valid) {
      throw std::runtime_error("Mesh::flipEdge: chosen remesh evaluation "
                               "returned an invalid state.");
    }

    const auto &otherN1 = chooseA ? n1b : n1a;
    const auto &otherN2 = chooseA ? n2b : n2a;
    const TElement::EdgeFlipRemeshState otherState1 =
        e1.evaluateEdgeFlipRemeshState(otherN1,
                                       F_P_H[static_cast<size_t>(e1i)]);
    const TElement::EdgeFlipRemeshState otherState2 =
        e2.evaluateEdgeFlipRemeshState(otherN2,
                                       F_P_H[static_cast<size_t>(e2i)]);
    if (!otherState1.valid || !otherState2.valid) {
      throw std::runtime_error("Mesh::flipEdge: other remesh option returned "
                               "an invalid state.");
    }
    const double chosenEnergy = state1.energy_new + state2.energy_new;
    const double otherEnergy = otherState1.energy_new + otherState2.energy_new;
    edgeFlipChosenMinusOtherEnergyInStep += chosenEnergy - otherEnergy;
    if (chosenEnergy > otherEnergy) {
      edgeFlipAlwaysChoseLowerEnergyInStep = false;
    }
  }

  const Matrix2d oldP1 = elements[static_cast<size_t>(e1i)].F_P;
  const Matrix2d oldP2 = elements[static_cast<size_t>(e2i)].F_P;
  const Matrix2d oldH1 = F_P_H[static_cast<size_t>(e1i)];
  const Matrix2d oldH2 = F_P_H[static_cast<size_t>(e2i)];
  const std::array<GhostNode, 3> newGhostNodes1 =
      prepareEdgeFlipCandidate(n1, "Mesh::flipEdge element 1");
  const std::array<GhostNode, 3> newGhostNodes2 =
      prepareEdgeFlipCandidate(n2, "Mesh::flipEdge element 2");

  removeElementFromNodes(elements[e1i]);
  removeElementFromNodes(elements[e2i]);

  createElementPair(newGhostNodes1, newGhostNodes2, e1i, e2i, true);

  const Matrix2d newP1 = elements[static_cast<size_t>(e1i)].F_P;
  const Matrix2d newP2 = elements[static_cast<size_t>(e2i)].F_P;
  F_P_H[static_cast<size_t>(e1i)] = newP1.inverse() * oldP1 * oldH1;
  F_P_H[static_cast<size_t>(e2i)] = newP2.inverse() * oldP2 * oldH2;

  totalEdgeFlipsInStep++;
  throwIfReductionExploded(elements[e1i], "Mesh::fixElementPair");
  throwIfReductionExploded(elements[e2i], "Mesh::fixElementPair");
}

bool Mesh::checkAndFixPeriodicElementPair(TElement &e1, TElement &e2) {
  // Check if the coAngleNodes are actually the same. If they are in different
  // periodic images, we need to move them together.
  auto e1Co = e1.getCoNodesByIndex();
  auto e2Co = e2.getCoNodesByIndex();
  if (e1Co[0]->id == e2Co[0]->id && e1Co[1]->id == e2Co[1]->id) {
    // The element are already together, no need to move elements
    return false;
  }
  // The situation here is like this: We would usually fix an element pair
  // that looks like this:
  // c2____a3
  //  |  \  |
  // a0____c1
  // But now, this element pair is accros the periodic boundary, like this:
  //  c2____a3           c2
  //      \  |    ...    |  \
  //        c1          a0____c1
  // What we need to do now, is to move one of the elements to the other
  // element so they can be reconnected properly. We will do this by considering
  // which element is closest to the center of the system. We will then move
  // the other one.

  // We move the one that is furthest away from the center of the mesh

  double distE1 = (e1.getCom() - com).norm();
  double distE2 = (e2.getCom() - com).norm();
  if (distE1 > distE2) {
    moveElementToTwin(e1, e2);
  } else {
    moveElementToTwin(e2, e1);
  }
  return true;
}

void Mesh::moveElementToTwin(TElement &elementToMove,
                             const TElement &fixedElement) {
  auto fixedCoAngleNodes = fixedElement.getCoNodesByIndex();
  auto movingCoAngleNodes = elementToMove.getCoNodesByIndex();

  const Vector2i deltaShift = fixedCoAngleNodes[0]->periodicShift -
                              movingCoAngleNodes[0]->periodicShift;
  const Vector2i secondDeltaShift = fixedCoAngleNodes[1]->periodicShift -
                                    movingCoAngleNodes[1]->periodicShift;
  if ((secondDeltaShift.array() != deltaShift.array()).any()) {
    throw std::runtime_error(
        "moveElementToTwin: shared edge nodes disagree on periodic shift.");
  }

  std::array<GhostNode, 3> shiftedGhostNodes = elementToMove.ghostNodes;
  for (GhostNode &gn : shiftedGhostNodes) {
    gn.applyPeriodicShift(deltaShift, latticeBasis, currentDeformation);
  }

  // We remove the element from all the nodes it is connected to.
  const std::vector<Node *> oldNodes = {
      (*this)[elementToMove.ghostNodes[0].referenceId],
      (*this)[elementToMove.ghostNodes[1].referenceId],
      (*this)[elementToMove.ghostNodes[2].referenceId]};
  removeElementFromNodes(elementToMove);

  // overwrite the element
  elementToMove = TElement{(*this),
                           shiftedGhostNodes[0],
                           shiftedGhostNodes[1],
                           shiftedGhostNodes[2],
                           elementToMove.eIndex,
                           elementToMove.noise,
                           energyFunction,
                           bulkModulus};
}

int Mesh::countConnectionsInGhostNode(const GhostNode &gn) {
  int connections = 0;

  // First we find the reference node
  Node *n = (*this)[gn.referenceId];

  // Then we need to go through each element it is connected to, and compare
  // with the ghost node in those elements to see if they are the same ghost
  // node as the one we are considering. We check if they are the same by
  // comparing their row, col AND periodicShift.

  for (int i = 0; i < n->elementCount; i++) {
    // Get one of the elements connected to the reference node.
    TElement &e = elements[n->connectedElements[i]];
    // Find the ghost node in that element that represents the reference node.
    GhostNode &gnInElement = e.ghostNodes[n->nodeIndexInElement[i]];
    // Check if this node is the same node as our input ghost node.
    // If it is, it means that our input ghost node is connected to this
    // element. (That means that the minimum number of connections will almost
    // always be 1, since we count the element that the input ghost node comes
    // from.)
    if (gnInElement.id == gn.id) {
      connections += 1;
    }
  }

  return connections;
}

void Mesh::setDiagonal(int row, int col, bool majorDiagonalOrder) {
  // get the 4 ghost nodes of the selected section of the mesh
  const std::array<GhostNode, 4> ghosts = getSquareGhostNodes(row, col);
  // Get the indexes of the elements
  auto [e1i, e2i] = getElementIndices(row, col);

  // Now we need to remove e1i and e2i from the 4 nodes, since they will be
  // added back in the tElement constructor
  removeElementsFromNodes(row, col, {e1i, e2i});

  // Update elements, preserving existing noise values
  createElementPair(ghosts, e1i, e2i, majorDiagonalOrder, true);
}

void Mesh::removeElementFromNodes(const TElement &element) {

  std::vector<const GhostNode *> ghostNodesVector;
  ghostNodesVector.reserve(3);

  for (const auto &gn : element.ghostNodes) {
    ghostNodesVector.push_back(&gn);
  }

  removeElementsFromNodes(ghostNodesVector, {element.eIndex});
}

void Mesh::removeElementsFromNodes(int row, int col,
                                   const std::vector<int> elIndexToRemove) {
  std::array<Node *, 4> nodes = getSquareNodes(row, col);
  removeElementsFromNodes({nodes[0], nodes[1], nodes[2], nodes[3]},
                          elIndexToRemove);
}

void Mesh::removeElementsFromNodes(const std::vector<const GhostNode *> gNodes,
                                   const std::vector<int> elIndexToRemove) {
  std::vector<Node *> nodes;
  nodes.reserve(gNodes.size()); // Reserve space for efficiency

  std::transform(
      gNodes.begin(), gNodes.end(), std::back_inserter(nodes),
      [this](const GhostNode *g) { return (*this)[g->referenceId]; });

  removeElementsFromNodes(nodes, elIndexToRemove);
}

void Mesh::removeElementsFromNodes(std::vector<Node *> nodes,
                                   const std::vector<int> elIndexToRemove) {

  // We check each node
  for (Node *n : nodes) {
    std::vector<int> indexesToRemove;

    // Find indices of elements to remove
    for (int i = 0; i < n->elementCount; i++) {
      for (int j = 0; j < elIndexToRemove.size(); j++) {
        if (n->connectedElements[i] == elIndexToRemove[j]) {
          indexesToRemove.push_back(i);
          break; // Found a match, no need to check other elements
        }
      }
    }

    // If there's nothing to remove, continue to next node
    if (indexesToRemove.empty()) {
      continue;
    }

    // Create temporary arrays to store elements we want to keep
    int tempElementIndices[MAX_ELEMENTS_PER_NODE];
    int tempNodeIndexInElement[MAX_ELEMENTS_PER_NODE];
    int newCount = 0;

    // Copy only the elements we want to keep
    for (int i = 0; i < n->elementCount; i++) {
      // Check if this index should be removed
      bool shouldRemove = false;
      for (int removeIdx : indexesToRemove) {
        if (i == removeIdx) {
          shouldRemove = true;
          break;
        }
      }

      // If not marked for removal, keep it
      if (!shouldRemove) {
        tempElementIndices[newCount] = n->connectedElements[i];
        tempNodeIndexInElement[newCount] = n->nodeIndexInElement[i];
        newCount++;
      }
    }

    // Update the node with the new arrays
    for (int i = 0; i < newCount; i++) {
      n->connectedElements[i] = tempElementIndices[i];
      n->nodeIndexInElement[i] = tempNodeIndexInElement[i];
    }
    for (int i = newCount; i < n->elementCount; i++) {
      n->connectedElements[i] = -1;
      n->nodeIndexInElement[i] = -1;
    }

    // Update the count
    n->elementCount = newCount;
  }
}

// Helper function to update positions using a generic buffer and its size
