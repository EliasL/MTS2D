#include "node.h"
#include "mesh.h"

NodeId::NodeId() : i(0), row(0), col(0) {}
NodeId::NodeId(int row_, int col_, int cols)
    : i(row_ * cols + col_), row(row_), col(col_), cols(cols) {}

NodeId::NodeId(int i_, int cols)
    : i(i_), row(i_ / cols), col(i_ % cols), cols(cols) {}

std::ostream &operator<<(std::ostream &os, const NodeId &nodeId) {
  // This implementation is confusing because (3,4) resembles vector notation
  // where x=3 and y=4.
  // os << "Node " << nodeId.i << "(" << nodeId.col << ", " << nodeId.row <<
  // ")";

  // This implementation, while less compact, is clearer.
  os << "Node " << nodeId.i << ", row: " << nodeId.row
     << ", col: " << nodeId.col;
  return os;
}

Node::Node(double x, double y) {
  m_pos = {x, y};
  m_init_pos = {x, y};
  m_u = {0, 0};
  f = {0, 0};
  fixedNode = false;
  isGhostNode = false;
  ghostShift = {0, 0};
}

Node::Node(double a, int row, int col, int cols) : Node(a * col, a * row) {
  id = NodeId(row, col, cols);
  ghostId = NodeId(row, col, cols + 1);
}

Vector2d Node::pos() const { return m_pos; }
Vector2d Node::init_pos() const { return m_init_pos; }
Vector2d Node::u() const { return m_u; }

void Node::setPos(Vector2d pos) {
  m_pos = pos;
  updateDisplacement();
}

void Node::addPos(Vector2d pos) {
  m_pos += pos;
  updateDisplacement();
}

void Node::setInitPos(Vector2d init_pos) {
  m_init_pos = init_pos;
  updateDisplacement();
}

// Function to update displacement based on the current and initial positions.
void Node::updateDisplacement() { m_u = m_pos - m_init_pos; }

void Node::setDisplacement(Vector2d disp) {
  m_pos = m_init_pos + disp;
  m_u = disp;
}

void Node::addForce(Vector2d _f) { f += _f; }

void Node::resetForce() {
  // This sets all the values in f to 0
  f = {0, 0};
}

Node::Node() : Node(0, 0) {}

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

double tElementArea(const Node &A, const Node &B, const Node &C) {
  // Accessing private member m_pos for each Node, which holds the position.
  Vector2d posA = A.m_pos;
  Vector2d posB = B.m_pos;
  Vector2d posC = C.m_pos;

  double area = 0.5 * std::abs(posA[0] * (posB[1] - posC[1]) +
                               posB[0] * (posC[1] - posA[1]) +
                               posC[0] * (posA[1] - posB[1]));
  return area;
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

void translateInPlace(Node &n, Vector2d disp, double multiplier) {
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

// Operator overloading for Vector2d scalar multiplication
Vector2d operator*(const Vector2d &v, double scalar) {
  return {v[0] * scalar, v[1] * scalar};
}

Vector2d operator*(const Vector2d &lhs, const Vector2d &rhs) {
  return {lhs[0] * rhs[0], lhs[1] * rhs[1]};
}

// Operator overloading for Vector2d addition
Vector2d operator+(const Vector2d &v1, const Vector2d &v2) {
  return {v1[0] + v2[0], v1[1] + v2[1]};
}
// Element-wise vector addition
Vector2d &operator+=(Vector2d &lhs, const Vector2d &rhs) {
  for (long i = 0; i < lhs.size(); ++i) {
    lhs[i] += rhs[i];
  }
  return lhs;
}

// Element-wise vector subtraction
Vector2d operator-(const Vector2d &lhs, const Vector2d &rhs) {
  Vector2d result;
  for (long i = 0; i < lhs.size(); ++i) {
    result[i] = lhs[i] - rhs[i];
  }
  return result;
}
