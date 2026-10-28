#include "mesh.h"
#include "Simulation/randomUtils.h"
#include <CDT.h>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using Eigen::Vector2i;

namespace {
// Keep the real node and exact periodic image associated with each input point.
struct DelaunayVertex {
  int refNodeIndex;
  Vector2i periodicShift;
};
using DelaunayFace = std::array<CDT::VertInd, 3>;

// Unique key for a triangle based on its three real node ids.
using TriKey = std::array<int, 3>;

TriKey makeTriKey(const DelaunayFace &face,
                  const std::vector<DelaunayVertex> &vertices) {
  TriKey ids{};
  for (int k = 0; k < 3; ++k) {
    ids[k] = vertices[face[k]].refNodeIndex;
  }
  std::sort(ids.begin(), ids.end());
  return ids;
}

// Use relative coordinates and extended precision for circumcenters; the
// periodic ownership and halo certification below still use physical doubles.
Vector2d circumcenter(const Vector2d &a, const Vector2d &b,
                            const Vector2d &c) {
  const long double bx = static_cast<long double>(b.x()) - a.x();
  const long double by = static_cast<long double>(b.y()) - a.y();
  const long double cx = static_cast<long double>(c.x()) - a.x();
  const long double cy = static_cast<long double>(c.y()) - a.y();
  const long double divisor = 2 * (bx * cy - by * cx);
  if (divisor == 0) {
    throw std::runtime_error("Mesh::reconnectDelaunay: collinear triangle.");
  }
  const long double bb = bx * bx + by * by;
  const long double cc = cx * cx + cy * cy;
  const Vector2d center{static_cast<double>(a.x() + (cy * bb - by * cc) / divisor),
                        static_cast<double>(a.y() + (bx * cc - cx * bb) / divisor)};
  if (!center.allFinite()) {
    throw std::runtime_error("Mesh::reconnectDelaunay: invalid circumcenter.");
  }
  return center;
}

double pointSegmentDistanceSquared(const Vector2d &p, const Vector2d &a,
                                          const Vector2d &b) {
  const Vector2d edge = b - a;
  const double edgeLengthSquared = edge.squaredNorm();
  if (edgeLengthSquared <= 0.0) {
    return (p - a).squaredNorm();
  }
  const double t = std::clamp((p - a).dot(edge) / edgeLengthSquared, 0.0, 1.0);
  return (p - (a + t * edge)).squaredNorm();
}

double pointCellDistanceSquared(const Vector2d &p,
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

int checkedPeriodicShift(double latticeCoordinate, int period) {
  if (!std::isfinite(latticeCoordinate) || period <= 0) {
    throw std::runtime_error(
        "checkedPeriodicShift: invalid lattice coordinate or period.");
  }
  const double cell = std::floor(latticeCoordinate / period);
  const double shift = -cell * period;
  if (shift < INT_MIN ||
      shift > INT_MAX) {
    throw std::runtime_error(
        "checkedPeriodicShift: periodic shift exceeds integer range.");
  }
  return static_cast<int>(shift);
}

Vector2i canonicalPeriodicShift(const Node &node,
                                       const Matrix2d &physicalLatticeInverse,
                                       int rows, int cols) {
  const Vector2d q = physicalLatticeInverse * node.pos();
  return {checkedPeriodicShift(q.x(), cols), checkedPeriodicShift(q.y(), rows)};
}

Vector2i checkedCellShift(int cellX, int cellY, int rows, int cols) {
  const long long x = static_cast<long long>(cellX) * cols;
  const long long y = static_cast<long long>(cellY) * rows;
  if (x < INT_MIN ||
      x > INT_MAX ||
      y < INT_MIN ||
      y > INT_MAX) {
    throw std::runtime_error("checkedCellShift: shift exceeds integer range.");
  }
  return {static_cast<int>(x), static_cast<int>(y)};
}

bool ownsCircumcenter(const Vector2d &center,
                             const Matrix2d &cellBasisInverse) {
  Vector2d q = cellBasisInverse * center;
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

void validatePeriodicTopology(const std::vector<DelaunayFace> &faces,
                         const std::vector<DelaunayVertex> &vertices,
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
    if (!triangleKeys.insert(makeTriKey(face, vertices)).second) {
      throw std::runtime_error(
          "Mesh::reconnectDelaunay: duplicate canonical triangle after "
          "geometric periodic ownership.");
    }
    std::array<int, 3> ids{};
    for (int k = 0; k < 3; ++k) {
      ids[k] = vertices[face[k]].refNodeIndex;
    }
    if (ids[0] == ids[1] || ids[1] == ids[2] || ids[2] == ids[0]) {
      std::string detail;
      for (int k = 0; k < 3; ++k) {
        const auto &info = vertices[face[k]];
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

} // namespace

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
  std::vector<Vector2d> points;
  std::vector<DelaunayVertex> vertices;
  std::vector<DelaunayFace> kept_faces;
  bool haloCertified = false;
  for (int attempt = 0; attempt < maxHaloAttempts; ++attempt) {
    points.clear();
    vertices.clear();
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
          points.push_back(position);
          vertices.push_back({nodeIndex, canonicalShift + cellShift});
        }
      }
    }

    // Ordinary Delaunay only: no constraints, refinement or added mesh nodes.
    // CDT rejects duplicate input positions instead of merging physical nodes.
    CDT::Triangulation<double> dt(CDT::VertexInsertionOrder::Auto);
    dt.insertVertices(points.begin(), points.end(),
                      [](const Vector2d &p) { return p.x(); },
                      [](const Vector2d &p) { return p.y(); });
    dt.eraseSuperTriangle();
    if (dt.vertices.size() != points.size()) {
      throw std::runtime_error(
          "Mesh::reconnectDelaunay: triangulation changed the vertex count.");
    }
    double maxCircumradius = 0.0;
    for (const auto &face : dt.triangles) {
      const Vector2d &p0 = points[face.vertices[0]];
      const Vector2d &p1 = points[face.vertices[1]];
      const Vector2d &p2 = points[face.vertices[2]];
      const Vector2d b = p1 - p0, c = p2 - p0;
      const double area = 0.5 * (b.x() * c.y() - b.y() * c.x());
      if (area < 1e-5) {
        continue;
      }
      const Vector2d center = circumcenter(p0, p1, p2);
      if (!ownsCircumcenter(center, cellBasisInverse)) {
        continue;
      }
      kept_faces.push_back(face.vertices);
      maxCircumradius = std::max(maxCircumradius, (center - p0).norm());
    }

    if (maxCircumradius > haloWidth * (1.0 - 1e-10)) {
      haloWidth = std::max(2.0 * haloWidth, 1.05 * maxCircumradius);
      continue;
    }
    validatePeriodicTopology(kept_faces, vertices, nrElements);
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

  // Refill elements from triangle indices (each face -> one TElement).
  int eIdx = 0;
  for (const auto &f : kept_faces) {
    std::array<GhostNode, 3> g;
    for (int k = 0; k < 3; ++k) {
      const auto &vi = vertices[f[k]];
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
                                                 : rng::sampleNormal(1, QDSD);
    const Matrix2d previousPlastic = elements[eIdx].F_P;
    const Matrix2d previousHistory = elements[eIdx].T_p;
    elements[eIdx] = TElement((*this), g[0], g[1], g[2], eIdx, noise,
                              energyFunction, bulkModulus);
    elements[eIdx].T_p = elements[eIdx].F_P * previousPlastic.inverse() * previousHistory;
    ++eIdx;
  }

  markDirty();
  nrMinItterationsSinceLastReconnect = 0;
  nrMinFunctionCallsSinceLastReconnect = 0;
}
