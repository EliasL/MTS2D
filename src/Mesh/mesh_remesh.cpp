#include "mesh.h"
#include "Data/logging.h"
#include "Mesh/node.h"
#include "Mesh/tElement.h"
#include "Simulation/randomUtils.h"
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
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

using K = CGAL::Exact_predicates_exact_constructions_kernel; // exact predicates
                                                             // & constructions
using Point = K::Point_2;
using Delaunay = CGAL::Delaunay_triangulation_2<K>;
// Keep the original node id and the exact integer periodic image used to
// create the Delaunay vertex.  The latter must be carried through the
// triangulation instead of reconstructed from a floating-point position.
struct VInfo {
  int refNodeIndex = -1; // index of the reference node in the MTS mesh
                         // (0..nrNodes-1)
  Vector2i periodicShift = Vector2i::Zero();
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

  // Lev's double triangle
  bool useLevsDoubleTriangle = false; // Should always be false
  if (useLevsDoubleTriangle) {
    if (abs(b) > std::min(a, c)) {
      return false;
    }
  } else {
    if (b < 0) {
      return false;
    }
    if (b > std::min(a, c)) {
      return false;
    }
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
// Unique key for a triangle based on its three real node ids.
struct TriKey {
  uint64_t nodeIdA;
  uint64_t nodeIdB;
  uint64_t nodeIdC;
  bool operator==(const TriKey &o) const noexcept {
    return nodeIdA == o.nodeIdA && nodeIdB == o.nodeIdB && nodeIdC == o.nodeIdC;
  }
  // Define less-than operator for use in ordered representative selection.
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

static Vector2d pointVector(const Point &p) {
  return {CGAL::to_double(p.x()), CGAL::to_double(p.y())};
}

static double pointSegmentDistanceSquared(const Vector2d &p, const Vector2d &a,
                                          const Vector2d &b) {
  const Vector2d edge = b - a;
  const double edgeLengthSquared = edge.squaredNorm();
  if (edgeLengthSquared <= 0.0) {
    return (p - a).squaredNorm();
  }
  const double t = std::clamp((p - a).dot(edge) / edgeLengthSquared, 0.0, 1.0);
  return (p - (a + t * edge)).squaredNorm();
}

static double pointCellDistanceSquared(const Vector2d &p,
                                       const Matrix2d &cellBasis,
                                       const Matrix2d &cellBasisInverse) {
  const Vector2d q = cellBasisInverse * p;
  if (q.x() >= 0.0 && q.x() <= 1.0 && q.y() >= 0.0 && q.y() <= 1.0) {
    return 0.0;
  }

  const Vector2d origin = Vector2d::Zero();
  const Vector2d a = cellBasis.col(0);
  const Vector2d b = cellBasis.col(1);
  const Vector2d ab = a + b;
  return std::min({pointSegmentDistanceSquared(p, origin, a),
                   pointSegmentDistanceSquared(p, origin, b),
                   pointSegmentDistanceSquared(p, a, ab),
                   pointSegmentDistanceSquared(p, b, ab)});
}

static int checkedPeriodicShift(double latticeCoordinate, int period) {
  if (!std::isfinite(latticeCoordinate) || period <= 0) {
    throw std::runtime_error(
        "checkedPeriodicShift: invalid lattice coordinate or period.");
  }
  const double cell = std::floor(latticeCoordinate / period);
  const double shift = -cell * period;
  if (shift < std::numeric_limits<int>::min() ||
      shift > std::numeric_limits<int>::max()) {
    throw std::runtime_error(
        "checkedPeriodicShift: periodic shift exceeds integer range.");
  }
  return static_cast<int>(shift);
}

static Vector2i canonicalPeriodicShift(const Node &node,
                                       const Matrix2d &physicalLatticeInverse,
                                       int rows, int cols) {
  const Vector2d q = physicalLatticeInverse * node.pos();
  return {checkedPeriodicShift(q.x(), cols), checkedPeriodicShift(q.y(), rows)};
}

static Vector2i checkedCellShift(int cellX, int cellY, int rows, int cols) {
  const long long x = static_cast<long long>(cellX) * cols;
  const long long y = static_cast<long long>(cellY) * rows;
  if (x < std::numeric_limits<int>::min() ||
      x > std::numeric_limits<int>::max() ||
      y < std::numeric_limits<int>::min() ||
      y > std::numeric_limits<int>::max()) {
    throw std::runtime_error("checkedCellShift: shift exceeds integer range.");
  }
  return {static_cast<int>(x), static_cast<int>(y)};
}

static bool ownsCircumcenter(const Point &center,
                             const Matrix2d &cellBasisInverse) {
  Vector2d q = cellBasisInverse * pointVector(center);
  constexpr double boundaryTolerance = 1e-12;
  for (int i = 0; i < 2; ++i) {
    if (std::abs(q[i]) <= boundaryTolerance) {
      q[i] = 0.0;
    } else if (std::abs(q[i] - 1.0) <= boundaryTolerance) {
      q[i] = 1.0;
    }
  }
  return q.x() >= 0.0 && q.x() < 1.0 && q.y() >= 0.0 && q.y() < 1.0;
}

static void insertDelaunayPoint(DelaunayInfo &dt, const Vector2d &position,
                                int refNodeIndex,
                                const Vector2i &periodicShift) {
  const size_t oldVertexCount = dt.number_of_vertices();
  auto vertex = dt.insert(Point(position.x(), position.y()));
  if (dt.number_of_vertices() == oldVertexCount) {
    const VInfo &existing = vertex->info();
    if (existing.refNodeIndex != refNodeIndex ||
        (existing.periodicShift.array() != periodicShift.array()).any()) {
      throw std::runtime_error(
          "Mesh::reconnectDelaunay: distinct periodic nodes occupy exactly "
          "the same current position.");
    }
    return;
  }
  vertex->info().refNodeIndex = refNodeIndex;
  vertex->info().periodicShift = periodicShift;
}

static void
validatePeriodicTopology(const std::vector<DelaunayInfo::Face_handle> &faces,
                         int expectedFaces) {
  if (static_cast<int>(faces.size()) != expectedFaces) {
    throw std::runtime_error(
        "Mesh::reconnectDelaunay: geometric periodic ownership produced an "
        "unexpected triangle count (expected " +
        std::to_string(expectedFaces) + ", got " +
        std::to_string(faces.size()) + ").");
  }

  std::set<TriKey> triangleKeys;
  std::unordered_map<Mesh::EdgeKey, int, Mesh::EdgeKeyHash> edgeIncidence;
  for (const auto &face : faces) {
    if (!triangleKeys.insert(makeTriKey(face)).second) {
      throw std::runtime_error(
          "Mesh::reconnectDelaunay: duplicate canonical triangle after "
          "geometric periodic ownership.");
    }
    std::array<int, 3> ids{};
    for (int k = 0; k < 3; ++k) {
      ids[k] = face->vertex(k)->info().refNodeIndex;
    }
    if (ids[0] == ids[1] || ids[1] == ids[2] || ids[2] == ids[0]) {
      std::string detail;
      for (int k = 0; k < 3; ++k) {
        const VInfo &info = face->vertex(k)->info();
        detail += " [node=" + std::to_string(info.refNodeIndex) + ", shift=(" +
                  std::to_string(info.periodicShift.x()) + "," +
                  std::to_string(info.periodicShift.y()) + ")]";
      }
      throw std::runtime_error("Mesh::reconnectDelaunay: periodic triangle "
                               "repeats a real node:" +
                               detail);
    }
    for (int k = 0; k < 3; ++k) {
      edgeIncidence[Mesh::EdgeKey(ids[k], ids[(k + 1) % 3])]++;
    }
  }
  for (const auto &[edge, count] : edgeIncidence) {
    (void)edge;
    if (count != 2) {
      throw std::runtime_error(
          "Mesh::reconnectDelaunay: canonical periodic edge incidence is "
          "not two.");
    }
  }
}

void Mesh::reconnectDelaunay() {
  const Matrix2d physicalLattice = currentDeformation * latticeBasis;
  if (!physicalLattice.allFinite() ||
      std::abs(physicalLattice.determinant()) <= 1e-14) {
    throw std::runtime_error(
        "Mesh::reconnectDelaunay: singular current periodic lattice.");
  }
  const Matrix2d physicalLatticeInverse = physicalLattice.inverse();
  Matrix2d cellBasis = physicalLattice;
  cellBasis.col(0) *= cols;
  cellBasis.col(1) *= rows;
  const Matrix2d cellBasisInverse = cellBasis.inverse();

  Eigen::JacobiSVD<Matrix2d> svd(cellBasis);
  const double minimumCellStretch = svd.singularValues().minCoeff();
  if (!std::isfinite(minimumCellStretch) || minimumCellStretch <= 0.0) {
    throw std::runtime_error(
        "Mesh::reconnectDelaunay: invalid periodic-cell singular value.");
  }

  double haloWidth = 2.0 * std::max(physicalLattice.col(0).norm(),
                                    physicalLattice.col(1).norm());
  constexpr int maxHaloAttempts = 8;
  DelaunayInfo dt;
  std::vector<DelaunayInfo::Face_handle> kept_faces;
  bool haloCertified = false;
  for (int attempt = 0; attempt < maxHaloAttempts; ++attempt) {
    dt.clear();
    kept_faces.clear();

    const int translationRange =
        static_cast<int>(
            std::ceil(std::sqrt(2.0) + haloWidth / minimumCellStretch)) +
        1;
    const double haloWidthSquared = haloWidth * haloWidth;
    for (int nodeIndex = 0; nodeIndex < nrNodes; ++nodeIndex) {
      const Node &node = nodes(nodeIndex);
      const Vector2i canonicalShift =
          canonicalPeriodicShift(node, physicalLatticeInverse, rows, cols);
      const Vector2d canonicalPosition =
          node.pos() + physicalLattice * canonicalShift.cast<double>();
      if (pointCellDistanceSquared(canonicalPosition, cellBasis,
                                   cellBasisInverse) >
          1e-18 * cellBasis.squaredNorm()) {
        throw std::runtime_error(
            "Mesh::reconnectDelaunay: failed to wrap a node into the current "
            "periodic cell.");
      }

      for (int cellY = -translationRange; cellY <= translationRange; ++cellY) {
        for (int cellX = -translationRange; cellX <= translationRange;
             ++cellX) {
          const Vector2i cellShift = checkedCellShift(cellX, cellY, rows, cols);
          const Vector2d position =
              canonicalPosition + physicalLattice * cellShift.cast<double>();
          if (cellX != 0 || cellY != 0) {
            const double distanceSquared =
                pointCellDistanceSquared(position, cellBasis, cellBasisInverse);
            if (distanceSquared > haloWidthSquared) {
              continue;
            }
          }
          insertDelaunayPoint(dt, position, nodeIndex,
                              canonicalShift + cellShift);
        }
      }
    }

    double maxCircumradius = 0.0;
    for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end();
         ++face) {
      const Point &p0 = face->vertex(0)->point();
      const Point &p1 = face->vertex(1)->point();
      const Point &p2 = face->vertex(2)->point();
      const K::FT area = CGAL::area(p0, p1, p2);
      if (area < K::FT(1e-5)) {
        continue;
      }
      const Point center = CGAL::circumcenter(p0, p1, p2);
      if (!ownsCircumcenter(center, cellBasisInverse)) {
        continue;
      }
      kept_faces.push_back(face);
      maxCircumradius = std::max(
          maxCircumradius, (pointVector(center) - pointVector(p0)).norm());
    }

    if (maxCircumradius > haloWidth * (1.0 - 1e-10)) {
      haloWidth = std::max(2.0 * haloWidth, 1.05 * maxCircumradius);
      continue;
    }
    validatePeriodicTopology(kept_faces, nrElements);
    haloCertified = true;
    break;
  }
  if (!haloCertified) {
    throw std::runtime_error(
        "Mesh::reconnectDelaunay: adaptive geometric halo did not converge.");
  }

  // Reset per-node connectivity.
  for (long i = 0; i < nodes.size(); ++i) {
    nodes(i).elementCount = 0;
  }

  // Refill elements from CGAL faces (each face -> one TElement).
  int eIdx = 0;
  for (const auto &f : kept_faces) {
    std::array<GhostNode, 3> g;
    for (int k = 0; k < 3; ++k) {
      auto v = f->vertex(k);
      const VInfo &vi = v->info();
      const Node *n = &nodes(vi.refNodeIndex); // reference node in MTS mesh
      // Build the ghost from the exact integer image stored on the Delaunay
      // vertex. Do not infer the image from the floating-point point position.
      g[k] = GhostNode(n, vi.periodicShift, cols, latticeBasis,
                       currentDeformation, referenceDeformation);
    }

    // Delaunay is defined in current space, so its triangle can connect
    // collinear reference-grid nodes. Reuse the edge-flip remeshing rule to
    // assign a valid square-lattice reference triangle before construction.
    g = prepareEdgeFlipCandidate(g, "Mesh::reconnectDelaunay eIndex=" +
                                        std::to_string(eIdx));

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
