#include "node.h"
#include "Eigen/Core"
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

NodeId::NodeId() : i(0), idPos(0, 0), cols(0) {}

NodeId::NodeId(int row_, int col_, int cols_)
    : i(row_ * cols_ + col_), idPos(col_, row_), cols(cols_) {}

NodeId::NodeId(int i_, int cols_)
    : i(i_), idPos(i_ % cols_, i_ / cols_), cols(cols_) {}

std::ostream &operator<<(std::ostream &os, const NodeId &nodeId) {
  os << "Node " << nodeId.i << ", row: " << nodeId.row()
     << ", col: " << nodeId.col();
  return os;
}

Node::Node(double x, double y) {
  m_pos = {x, y};
  m_ref_pos = {x, y};
  m_u = {0, 0};
  f = {0, 0};
  fixedNode = false;

  connectedElements.fill(-1);
  nodeIndexInElement.fill(-1);
}

Node::Node(int row, int col, int cols) : Node(1, row, col, cols) {}
Node::Node(double a, int row, int col, int cols) : Node(a * col, a * row) {
  id = NodeId(row, col, cols);
}

void Node::setPos(const Vector2d &pos) {
  m_pos = pos;
  updateDisplacement();
}

void Node::addPos(const Vector2d &pos) {
  m_pos += pos;
  updateDisplacement();
}

void Node::setRefPos(const Vector2d &ref_pos) {
  m_ref_pos = ref_pos;
  updateDisplacement();
}

// Function to update displacement based on the current and initial positions.
void Node::updateDisplacement() { m_u = m_pos - m_ref_pos; }

void Node::setDisplacement(const Vector2d &disp) {
  m_pos = m_ref_pos + disp;
  m_u = disp;
}

void Node::addDisplacement(const Vector2d &dispChange) {
  m_u += dispChange;
  m_pos = m_ref_pos + m_u;
}

void Node::addForce(const Vector2d &_f) { f += _f; }

void Node::resetForce() {
  // This sets all the values in f to 0
  f.setZero();
}
void Node::applyDeformation(const Matrix2d &deformation) {
  setPos(deformation * m_pos);
}

Node::Node() : Node(0, 0) {}
GhostNode::GhostNode(const Node *referenceNode, Vector2i periodicShift,
                     int cols, const Matrix2d &latticeBasis,
                     const Matrix2d &deformation,
                     const Matrix2d &referenceDeformation)
    : referenceId(referenceNode->id), id(periodicShift + referenceId.idPos),
      periodicShift(periodicShift) {
  Vector2d shift = latticeBasis * periodicShift.cast<double>();
  pos = referenceNode->pos() + deformation * shift;
  ref_pos = referenceNode->ref_pos() + referenceDeformation * shift;
  u = pos - ref_pos;
}

GhostNode::GhostNode(const Node *referenceNode, int row, int col, int cols,
                     const Matrix2d &latticeBasis, const Matrix2d &deformation,
                     const Matrix2d &referenceDeformation)
    : GhostNode(referenceNode, Vector2i{col, row} - referenceNode->id.idPos,
                cols, latticeBasis, deformation, referenceDeformation) {}

GhostNode::GhostNode(const Node *referenceNode, const Matrix2d &latticeBasis,
                     const Matrix2d &deformation,
                     const Matrix2d &referenceDeformation)
    : GhostNode(referenceNode, Vector2i::Zero(), referenceNode->id.cols,
                latticeBasis, deformation, referenceDeformation) {}

GhostNode::GhostNode(const Node *referenceNode, int row, int col, int cols,
                     const Matrix2d &deformation)
    : GhostNode(referenceNode, row, col, cols, Matrix2d::Identity(),
                deformation) {}

GhostNode::GhostNode(const Node *referenceNode, const Matrix2d &deformation)
    : GhostNode(referenceNode, Matrix2d::Identity(), deformation) {}

// Return integer shift (in grid units) that brings ref near target.
inline Vector2i nearestShift(const Eigen::Vector2d &refNodePos,
                             const Eigen::Vector2d &targetPos,
                             const Eigen::Matrix2d &Ainv, // inverse deformation
                             int cols, int rows) {

  // We undo the deformation to work in the undeformed lattice
  const Eigen::Vector2d rL = Ainv * refNodePos;
  const Eigen::Vector2d tL = Ainv * targetPos;
  Eigen::Vector2d dL = tL - rL;

  // Scale by cell counts and round to nearest integer cell translation
  // (componentwise minimal image)
  const double sx = std::round(dL.x() / cols);
  const double sy = std::round(dL.y() / rows);
  if (abs(sx) > 100) {
    std::cerr << "Warning: large periodic shift sx=" << sx
              << " dL.x()=" << dL.x() << " cols=" << cols << "\n";
  }
  if (abs(sy) > 100) {
    std::cerr << "Warning: large periodic shift sx=" << sx
              << " dL.x()=" << dL.x() << " cols=" << cols << "\n";
  }

  return Vector2i{int(sx) * cols, int(sy) * rows};
}

GhostNode::GhostNode(const Node *referenceNode, Vector2d targetPos, int rows,
                     int cols, const Matrix2d &latticeBasis,
                     const Matrix2d &deformation,
                     const Matrix2d &referenceDeformation)
    : referenceId(referenceNode->id) {
  const Matrix2d Ainv = (deformation * latticeBasis).inverse();
  periodicShift =
      nearestShift(referenceNode->pos(), targetPos, Ainv, cols, rows);
  id = periodicShift + referenceId.idPos;
  Vector2d shift = latticeBasis * periodicShift.cast<double>();
  pos = referenceNode->pos() + deformation * shift;
  ref_pos = referenceNode->ref_pos() + referenceDeformation * shift;
  u = pos - ref_pos;
}

GhostNode::GhostNode(Vector2d currentPos, Vector2d referencePos)
    : id(Vector2i::Zero()), f(Vector2d::Zero()), periodicShift(Vector2i::Zero()),
      pos(currentPos), ref_pos(referencePos), u(currentPos - referencePos) {
  referenceId = NodeId();
  referenceId.i = -1;
}

void GhostNode::updateCurrentPosition(const Node *referenceNode,
                                      const Matrix2d &deformation,
                                      const Matrix2d &latticeBasis,
                                      const Matrix2d &referenceDeformation) {
  Vector2d shift = latticeBasis * periodicShift.cast<double>();
  pos = referenceNode->pos() + deformation * shift;
  u = pos - ref_pos;
}

void GhostNode::updateReferencePosition(Vector2d new_ref_pos) {

  // Note that this function does not account for periodic shifts!
  // That is because updatePosition extracts information from its
  // reference node, which might be in a different periodic image.
  // Here, the correct possition should be given directly.
  ref_pos = new_ref_pos;
  u = pos - ref_pos;
}

void GhostNode::transformReferencePosition(Matrix2d trans, Vector2d oldAnchor,
                                           Vector2d newAnchor) {
  ref_pos = newAnchor + trans * (ref_pos - oldAnchor);
  u = pos - ref_pos;
}

void GhostNode::applyPeriodicShift(const Vector2i &deltaShift,
                                   const Matrix2d &latticeBasis,
                                   const Matrix2d &currentDeformation) {
  periodicShift += deltaShift;
  id = referenceId.idPos + periodicShift;

  const Vector2d shift = latticeBasis * deltaShift.cast<double>();
  pos += currentDeformation * shift;
  u = pos - ref_pos;
}

std::ostream &operator<<(std::ostream &os, const Node &node) {
  // This implementation is confusing because (3,4) resembles vector notation
  // where x=3 and y=4.
  // os << "Node " << nodeId.i << "(" << nodeId.col << ", " << nodeId.row <<
  // ")";

  // This implementation, while less compact, is clearer.
  os << "Node " << node.id.i << ", pos: " << node.pos()
     << " disp: " << node.u();

  return os;
}
std::ostream &operator<<(std::ostream &os, const GhostNode &node) {
  // This implementation is confusing because (3,4) resembles vector notation
  // where x=3 and y=4.
  // os << "Node " << nodeId.i << "(" << nodeId.col << ", " << nodeId.row <<
  // ")";

  // This implementation, while less compact, is clearer.
  os << "GNode " << node.referenceId.i << ", pos: " << node.pos
     << " disp: " << node.u;
  // NOTE This only holds when the system deformation is identity
  os << " pShift: " << node.periodicShift.transpose();

  return os;
}

void transformInPlace(const Matrix2d &matrix, Node &n) {
  auto test = n.pos();
  n.setPos(matrix * test);
}

Node transform(const Matrix2d &matrix, const Node &n) {
  Node result = n;
  transformInPlace(matrix, result);
  return result;
}

void translateInPlace(Node &n, const Vector2d &disp, double multiplier) {
  // Update the node's position
  n.setPos(n.pos() + disp * multiplier);
}

void translateInPlace(Node &n, double dx, double dy, double multiplier) {
  translateInPlace(n, Vector2d{dx, dy}, multiplier);
}

void translateInPlace(Node &n, const Node &delta, double multiplier) {
  translateInPlace(n, delta.pos(), multiplier);
}

Node translate(const Node &n, const Node &delta, double multiplier) {
  Node result = n;
  translateInPlace(result, delta, multiplier);
  return result;
}

// Overload the << operator for Vector2d
std::ostream &operator<<(std::ostream &os, const Vector2d &arr) {
  os << "(";
  for (long i = 0; i < arr.size(); ++i) {
    os << arr[i];
    if (i < arr.size() - 1) {
      os << ", ";
    }
  }
  os << ")";
  return os;
}
