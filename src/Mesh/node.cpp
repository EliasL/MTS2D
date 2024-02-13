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
    x = x_;
    y = y_;
    f_x = f_y = 0;
    fixedNode = false;
}

void Node::setPos(double x_, double y_)
{
    x = x_;
    y = y_;
    updateDisplacement();
}

void Node::setInitPos(double x, double y)
{
    init_x = x;
    init_y = y;
    updateDisplacement();
}

// Function to update displacement based on the current and initial positions.
void Node::updateDisplacement()
{
    u_x = x - init_x;
    u_y = y - init_y;
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

Node::Node() : Node(0, 0) {}

void transformInPlace(const Matrix2x2<double> &matrix, Node &n)
{
    double newX = matrix[0][0] * n.X() + matrix[0][1] * n.Y();
    double newY = matrix[1][0] * n.X() + matrix[1][1] * n.Y();
    n.setPos(newX, newY);
}

Node transform(const Matrix2x2<double> &matrix, const Node &n)
{
    Node result = n;
    transformInPlace(matrix, result);
    return result;
}

void translateInPlace(Node& n, double dx, double dy, double multiplier)
{
    // Calculate the new position
    double newX = n.X() + dx * multiplier;
    double newY = n.Y() + dy * multiplier;
    
    // Update the node's position
    n.setPos(newX, newY);
}

void translateInPlace(Node &n, const Node &delta, double multiplier)
{
    translateInPlace(n, delta.X(), delta.Y(), multiplier);
}

Node translate(const Node &n, const Node &delta, double multiplier)
{
    Node result = n;
    translateInPlace(result, delta, multiplier);
    return result;
}