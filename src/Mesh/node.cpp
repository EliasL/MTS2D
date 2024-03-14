#include "node.h"
#include "mesh.h"

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

void Node::copyValues(PeriodicNode node)
{
    setInitPos(node.init_pos());
    setPos(node.pos());
    f = node.f();
}

Node::Node() : Node(0, 0) {}

PeriodicNode::PeriodicNode(NodeId nodeId, std::weak_ptr<Mesh> mesh, bool shiftX, bool shiftY)
    : realId(nodeId), mesh(mesh)
{
    updatePeriodicity(shiftX, shiftY);
}

void PeriodicNode::addForce(VArray f)
{
    mesh[realId]->addForce(f);
}

VArray PeriodicNode::pos() const
{
    return mesh.periodicTransformation * displacement + mesh[realId]->pos();
}

VArray PeriodicNode::init_pos() const
{
    return mesh.periodicTransformation * displacement + mesh[realId]->init_pos();
}

VArray PeriodicNode::u() const
{
    return mesh[realId]->u();
}
VArray PeriodicNode::f() const
{
    return mesh[realId]->f;
}

void PeriodicNode::updatePeriodicity(bool shiftX, bool shiftY)
{
    // If we are supposed to shift in the x direction, we add a*cols to the
    // displacement. Same for y.
    displacement = {shiftX ? mesh.a * mesh.nodes.cols : 0,
                    shiftY ? mesh.a * mesh.nodes.rows : 0};
    // The id will be different since the number of columns in the non-periodic
    // mesh is one greater. If we shift, we always end up on the last row or column
    periodicId = {shiftY ? mesh.nodes.rows : mesh[realId]->id.row,
                  shiftX ? mesh.nodes.cols : mesh[realId]->id.col,
                  mesh.nodes.cols + 1};
}

void transformInPlace(const Matrix2x2<double> &matrix, Node &n)
{
    n.setPos(matrix * n.pos());
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
