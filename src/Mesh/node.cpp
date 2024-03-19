#include "node.h"
#include "mesh.h"

NodeId::NodeId() : row(0), col(0), i(0) {}
NodeId::NodeId(int row_, int col_, int cols)
    : row(row_),
      col(col_),
      cols(cols),
      i(row_ * cols + col_) {}
NodeId::NodeId(int i_, int cols)
    : row(i_ / cols),
      col(i_ % cols),
      cols(cols),
      i(i_) {}

std::ostream &operator<<(std::ostream &os, const NodeId &nodeId)
{
    // This implementation is confusing because (3,4) resembles vector notation
    // where x=3 and y=4.
    // os << "Node " << nodeId.i << "(" << nodeId.col << ", " << nodeId.row << ")";

    // This implementation, while less compact, is clearer.
    os << "Node " << nodeId.i << ", row: " << nodeId.row << ", col: " << nodeId.col;
    return os;
}

Node::Node(double x, double y)
{
    m_pos = {x, y};
    m_init_pos = {0, 0};
    m_u = m_pos;
    f = {0, 0};
    fixedNode = false;
}

void Node::setPos(VArray pos)
{
    m_pos = pos;
    updateDisplacement();
}

void Node::addPos(VArray pos)
{
    m_pos += pos;
    updateDisplacement();
}

void Node::setInitPos(VArray init_pos)
{
    m_init_pos = init_pos;
    updateDisplacement();
}

// Function to update displacement based on the current and initial positions.
void Node::updateDisplacement()
{
    m_u = m_pos - m_init_pos;
}

void Node::setDisplacement(VArray disp)
{
    m_pos = m_init_pos + disp;
    m_u = disp;
}

void Node::addForce(VArray _f)
{
    f += _f;
}

void Node::resetForce()
{
    // This sets all the values in f to 0
    f = 0;
}

void Node::copyForceAndDisplacement(const Node &node)
{
    // Note that we do NOT copy the position or initial position
    // This is so that we can copy the values from a ghost node
    setDisplacement(node.u());
    f = node.f;
}

Node::Node() : Node(0, 0) {}

void transformInPlace(const Matrix2x2<double> &matrix, Node &n)
{
    auto test = n.pos();
    n.setPos(matrix * test);
}

Node transform(const Matrix2x2<double> &matrix, const Node &n)
{
    Node result = n;
    transformInPlace(matrix, result);
    return result;
}

void translateInPlace(Node &n, VArray disp, double multiplier)
{
    // Update the node's position
    n.setPos(n.pos() + disp * multiplier);
}

void translateInPlace(Node &n, double dx, double dy, double multiplier)
{
    translateInPlace(n, VArray{dx, dy}, multiplier);
}

void translateInPlace(Node &n, const Node &delta, double multiplier)
{
    translateInPlace(n, delta.pos(), multiplier);
}

Node translate(const Node &n, const Node &delta, double multiplier)
{
    Node result = n;
    translateInPlace(result, delta, multiplier);
    return result;
}
