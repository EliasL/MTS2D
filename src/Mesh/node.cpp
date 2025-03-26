#include "node.h"
#include "Eigen/Core"
#include "Eigen/src/Core/Matrix.h"
#include <cassert>

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
  m_init_pos = {x, y};
  m_u = {0, 0};
  f = {0, 0};
  fixedNode = false;

  elementIndices.fill(-1);
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

void Node::setInitPos(const Vector2d &init_pos) {
  m_init_pos = init_pos;
  updateDisplacement();
}

// Function to update displacement based on the current and initial positions.
void Node::updateDisplacement() { m_u = m_pos - m_init_pos; }

void Node::setDisplacement(const Vector2d &disp) {
  m_pos = m_init_pos + disp;
  m_u = disp;
}

void Node::addDisplacement(const Vector2d &dispChange) {
  m_u += dispChange;
  m_pos = m_init_pos + m_u;
}

void Node::addForce(const Vector2d &_f) { f += _f; }

void Node::resetForce() {
  // This sets all the values in f to 0
  f = {0, 0};
}
void Node::applyDeformation(const Matrix2d &deformation) {
  setPos(deformation * m_pos);
}

Node::Node() : Node(0, 0) {}

GhostNode::GhostNode(const Node *referenceNode, int row, int col, int cols,
                     double a, const Matrix2d &deformation)
    : GhostNode(referenceNode, Vector2i{col, row} - referenceNode->id.idPos,
                cols, a, deformation) {}

GhostNode::GhostNode(const Node *referenceNode, Vector2i periodicShift,
                     int cols, double a, const Matrix2d &deformation)
    : referenceId(referenceNode->id), id(periodicShift + referenceId.idPos),
      periodShift(periodicShift) {
  updatePosition(referenceNode, deformation, a);
}

GhostNode::GhostNode(const Node *referenceNode, double a,
                     const Matrix2d &deformation)
    : GhostNode(referenceNode, referenceNode->id.idPos, referenceNode->id.cols,
                a, deformation) {}

GhostNode::GhostNode(const Node *referenceNode, int row, int col, int cols,
                     const Matrix2d &deformation)
    : GhostNode(referenceNode, row, col, cols, 1, deformation) {}

GhostNode::GhostNode(const Node *referenceNode, const Matrix2d &deformation)
    : GhostNode(referenceNode, 1, deformation) {}

void GhostNode::updatePosition(const Node *referenceNode,
                               const Matrix2d &deformation, double a) {
  Vector2d shift = periodShift.cast<double>() * a;
  pos = referenceNode->pos() + deformation * shift;
  init_pos = referenceNode->init_pos() + shift;
  u = pos - init_pos;
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
  os << " pShift: " << node.periodShift.transpose();

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
