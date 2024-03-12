#include "node.h"

NodeId::NodeId() : row(0), col(0), i(0) {}
NodeId::NodeId(int row_, int col_, int cols) : row(row_), col(col_), i(row_ * cols + col_) {}
NodeId::NodeId(int i_, int cols) : row(i_ / cols), col(i_ % cols), i(i_) {}

std::ostream &operator<<(std::ostream &os, const NodeId &nodeId)
{
    // This implementation is confusing because (3,4) resembles vector notation
    // where x=3 and y=4.
    // os << "Node " << nodeId.i << "(" << nodeId.col << ", " << nodeId.row << ")";

    // This implementation, while less compact, is clearer.
    os << "Node " << nodeId.i << ", row: " << nodeId.row << ", col: " << nodeId.col;
    return os;
}

Node::Node(double x_, double y_)
{
    m_x = x_;
    m_y = y_;
    f_x = f_y = 0;
    fixedNode = false;
}

void Node::setPos(double x_, double y_)
{
    m_x = x_;
    m_y = y_;
    updateDisplacement();
}

void Node::setInitPos(double x, double y)
{
    m_init_x = x;
    m_init_y = y;
    updateDisplacement();
}

// Function to update displacement based on the current and initial positions.
void Node::updateDisplacement()
{
    m_u_x = m_x - m_init_x;
    m_u_y = m_y - m_init_y;
}

void Node::addForce(std::array<double, 2> f)
{
    f_x += f[0];
    f_y += f[1];
}

void Node::resetForce()
{
    f_x = 0;
    f_y = 0;
}

void Node::update(PeriodicNode node)
{
    setInitPos(node.init_x(), node.init_y());
    setPos(node.x(), node.y());
    f_x = node.f_x();
    f_y = node.f_y();
}

Node::Node() : Node(0, 0) {}

PeriodicNode::PeriodicNode() {}

PeriodicNode::PeriodicNode(Node *realNode, double dx, double dy, int rows, int cols) : realNode(realNode), dx(dx), dy(dy)
{
    // If we are looping around, we know we will end up at the last row/col
    int newRow = (dy > 0 ? rows : realNode->id.row);
    int newCol = (dx > 0 ? cols : realNode->id.col);
    periodicId = NodeId(
        newRow,
        newCol,
        cols + 1);
}

void PeriodicNode::addForce(std::array<double, 2> f)
{
    realNode->addForce(f);
}

void transformInPlace(const Matrix2x2<double> &matrix, Node &n)
{
    double newX = matrix[0][0] * n.x() + matrix[0][1] * n.y();
    double newY = matrix[1][0] * n.x() + matrix[1][1] * n.y();
    n.setPos(newX, newY);
}

Node transform(const Matrix2x2<double> &matrix, const Node &n)
{
    Node result = n;
    transformInPlace(matrix, result);
    return result;
}

void translateInPlace(Node &n, double dx, double dy, double multiplier)
{
    // Calculate the new position
    double newX = n.x() + dx * multiplier;
    double newY = n.y() + dy * multiplier;

    // Update the node's position
    n.setPos(newX, newY);
}

void translateInPlace(Node &n, const Node &delta, double multiplier)
{
    translateInPlace(n, delta.x(), delta.y(), multiplier);
}

Node translate(const Node &n, const Node &delta, double multiplier)
{
    Node result = n;
    translateInPlace(result, delta, multiplier);
    return result;
}
