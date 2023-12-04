#include "triangle.h"

// e1 and e2 form basis vectors for the triangle
std::array<double, 2> Triangle::e12() const
{ // a is the lattice spacing of the gird
    return {
        a2->x - a1->x,
        a2->y - a1->y,
    };
}

std::array<double, 2> Triangle::e13() const
{ // a is the lattice spacing of the gird
    return {
        a3->x - a1->x,
        a3->y - a1->y,
    };
}

std::array<double, 2> Triangle::e23() const
{ // a is the lattice spacing of the gird
    return {
        a3->x - a2->x,
        a3->y - a2->y,
    };
}

// Provices a metric tensor for the triangle
Matrix2x2<double> Triangle::metric(MetricFunction f) const
{
    // Symetric matricies would be faster, but only slightly for 2x2 matrix
    Matrix2x2<double> m;
    auto e1_ = e12();
    auto e2_ = e13();
    // There are many ways to calculate a metric. The user can specify which
    // to use.
    switch (f)
    {
    case MetricFunction::faicella:
        m[0][0] = e1_[0] * e1_[0] + e1_[1] * e1_[1];
        m[1][1] = e2_[0] * e2_[0] + e2_[1] * e2_[1];
        m[1][0] = m[0][1] = e1_[0] * e2_[0] + e1_[1] * e2_[1];
        return m;
    case MetricFunction::epsilon_lineaire: // This function looks super strange
        m[0][0] = e1_[0] - 1;
        m[1][1] = e2_[1] - 1;
        m[1][0] = m[0][1] = e2_[0];
    default:
        throw std::invalid_argument("Invalid metric function");
        break;
    }
}

std::ostream &operator<<(std::ostream &os, const Triangle &triangle)
{
    os << "node 1: (" << triangle.a1->x << ", " << triangle.a1->y << "), "
       << "node 2: (" << triangle.a2->x << ", " << triangle.a2->y << "), "
       << "node 3: (" << triangle.a3->x << ", " << triangle.a3->y << ")";
    return os;
}

std::ostream &operator<<(std::ostream &os, const Triangle *trianglePtr)
{
    if (trianglePtr)
    {
        // Use the existing Triangle operator<< overload
        os << *trianglePtr;
    }
    else
    {
        // Handle the nullptr case
        os << "nullptr";
    }
    return os;
}