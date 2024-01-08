#include "matrix2x2.h"

template <typename T>
Matrix2x2<T>::Matrix2x2(T a)
{
    // Assuming data is a member variable that is some kind of container
    data.fill(a);
}

template <typename T>
Matrix2x2<T>::Matrix2x2()
{
    (*this)[0][0] = (*this)[1][1] = static_cast<T>(1);
    (*this)[0][1] = (*this)[1][0] = static_cast<T>(0);
}

// Depricated because of ambiguity and confusion of whether c12 or c21 should come first
// template <typename T>
// Matrix2x2<T>::Matrix2x2(T c11, T c12, T c21, T c22)
// {
//     (*this)[0][0] = c11;
//     (*this)[0][1] = c12;
//     (*this)[1][0] = c21;
//     (*this)[1][1] = c22;
// }

// Takes a list of ROWS, not columns
template <typename T>
Matrix2x2<T>::Matrix2x2(std::initializer_list<std::initializer_list<T>> list)
{
    auto it = list.begin();
    for (int i = 0; i < 2; ++i)
    {
        auto itRow = it->begin();
        for (int j = 0; j < 2; ++j)
        {
            // This is technically slightly slower than [i][j], but I prefer
            // this, and for a 2x2 matrix, it really doesn't matter
            (*this)[i][j] = *itRow;
            ++itRow;
        }
        ++it;
    }
}

template <typename T>
Matrix2x2<T>::Matrix2x2(std::array<T, 4> d)
{
    data = d;
}

// Change the sign of the value in [x][y]
template <typename T>
void Matrix2x2<T>::flip(int row, int col)
{
    (*this)[row][col] *= -1;
}

template <typename T>
void Matrix2x2<T>::swap(int row1, int col1, int row2, int col2)
{
    T temp;
    temp = (*this)[row1][col1];
    (*this)[row1][col1] = (*this)[row2][col2];
    (*this)[row2][col2] = temp;
}

template <typename T>
void Matrix2x2<T>::swapCols()
{
    swap(0, 0, 0, 1);
    swap(1, 0, 1, 1);
}

template <typename T>
void Matrix2x2<T>::setCol(std::array<T, 2> column, int col)
{
    (*this)[0][col] = column[0];
    (*this)[1][col] = column[1];
}

template <typename T>
void Matrix2x2<T>::setCols(std::array<T, 2> column1, std::array<T, 2> column2)
{
    setCol(column1, 0);
    setCol(column2, 1);
}

// Matrix multiplication optimized for 2x2 matrices
template <typename T>
Matrix2x2<T> Matrix2x2<T>::operator*(const Matrix2x2<T> &other) const
{
    Matrix2x2 result;
    result[0][0] = (*this)[0][0] * other[0][0] + (*this)[0][1] * other[1][0];
    result[0][1] = (*this)[0][0] * other[0][1] + (*this)[0][1] * other[1][1];
    result[1][0] = (*this)[1][0] * other[0][0] + (*this)[1][1] * other[1][0];
    result[1][1] = (*this)[1][0] * other[0][1] + (*this)[1][1] * other[1][1];
    return result;
}

template <typename T>
Matrix2x2<T> Matrix2x2<T>::operator*(T scalar) const
{
    Matrix2x2<T> result;
    for (size_t i = 0; i < data.size(); i++)
    {
        result.data[i] = this->data[i] * scalar;
    }
    return result;
} 

template <typename T>
Matrix2x2<T> operator*(T scalar, const Matrix2x2<T>& matrix)
{
    return matrix*scalar;
}

// Check if two matrixies are equal
template <typename T>
bool Matrix2x2<T>::operator==(const Matrix2x2 &other) const
{
    for (size_t i = 0; i < 4; ++i)
    {
        if (data[i] != other.data[i])
        {
            return false; // If any elements are not equal, return false
        }
    }
    return true; // All elements are equal
}

// Calculate the determinant of a 2x2 matrix
template <typename T>
T Matrix2x2<T>::det() const
{
    return (*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0];
}

// Matrix addition for 2x2 matrices
template <typename T>
Matrix2x2<T> Matrix2x2<T>::operator+(const Matrix2x2 &other) const
{
    Matrix2x2 result;
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            result[i][j] = (*this)[i][j] + other[i][j];
        }
    }
    return result;
}

// Transpose the 2x2 matrix
template <typename T>
Matrix2x2<T> Matrix2x2<T>::transpose() const
{
    Matrix2x2 result;
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            result[i][j] = (*this)[j][i];
        }
    }
    return result;
}

// Calculate the inverse of the 2x2 matrix
template <typename T>
Matrix2x2<T> Matrix2x2<T>::inverse() const
{
    T detVal = det();
    if (detVal == 0)
        throw std::runtime_error("Matrix is not invertible.");

    T invDet = 1 / detVal;
    Matrix2x2 result;
    result[0][0] = (*this)[1][1] * invDet;
    result[0][1] = -(*this)[0][1] * invDet;
    result[1][0] = -(*this)[1][0] * invDet;
    result[1][1] = (*this)[0][0] * invDet;
    return result;
}

// Calculates the similarityTransform of (this) matrix
template <typename T>
Matrix2x2<T> Matrix2x2<T>::similarityTransform(const Matrix2x2 &g) const
{
    return g.inverse() * (*this) * g;
}

template <typename T>
bool Matrix2x2<T>::isOrthogonal() const
{
    // First check determinant
    if (abs(abs(det()) - 1) > __DBL_EPSILON__)
    {
        return false;
    }
    // Then do a more thurough test
    Matrix2x2 inv = Matrix2x2(data);
    Matrix2x2 tran = Matrix2x2(data);
    return tran.transpose() == inv.inverse();
}

// Calculate the similarityTransform provided that g is orthogonal
template <typename T>
Matrix2x2<T> Matrix2x2<T>::congruenceTransform(const Matrix2x2 &g) const
{
#if defined(__OPTIMIZE__) && __OPTIMIZE__ < 3
    if (!g.isOrthogonal())
    {
        throw std::runtime_error("Provided matrix is not orthogonal");
    }
#endif

    return g.transpose() * (*this) * g;
}

// Calculate the similarityTransform provided (this) matrix is symetric, and g is orthogonal
template <typename T>
Matrix2x2<T> Matrix2x2<T>::sym_orth_conjugate(const Matrix2x2 &g) const
{
#if __OPTIMIZE__ < 3
    if ((*this)[0][1] != (*this)[1][0])
    {
        throw std::runtime_error("Matrix is not symmetric.");
    }
    if (!g.isOrthogonal())
    {
        throw std::runtime_error("Provided matrix is not orthogonal");
    }
#endif

    return congruenceTransform(g);
    // This does not work TODO!
    Matrix2x2<double> d;
    d[0][0] = (*this)[0][0] * g[0][0] * g[0][0] + (*this)[0][1] * g[0][0] * g[0][1] + (*this)[1][1] * g[0][1] * g[0][1];
    d[1][1] = (*this)[0][0] * g[1][0] * g[1][0] + (*this)[0][1] * g[1][0] * g[1][1] + (*this)[1][1] * g[1][1] * g[1][1];
    d[0][1] = g[1][0] * (2 * (*this)[0][0] * g[0][0] + (*this)[0][1] * g[0][1]) + g[1][1] * ((*this)[0][1] * g[0][0] + 2 * (*this)[1][1] * g[0][1]);
    d[1][0] = d[0][1];
    return d;
}

// The two next functions used to be static. I hope I don't break anything by removing that
template <typename T>
Matrix2x2<T> Matrix2x2<T>::rotation_matrix(double theta)
{
    return Matrix2x2({{cos(theta), -sin(theta)}, {sin(theta), cos(theta)}});
}
template <typename T>
Matrix2x2<T> Matrix2x2<T>::reflection_matrix(double theta)
{
    return Matrix2x2({{cos(theta), sin(theta)}, {sin(theta), -cos(theta)}});
}


template <typename T>
std::ostream &operator<<(std::ostream &os, const Matrix2x2<T> &m)
{
    // Helper lambda to format output based on type
    auto formatOutput = [&os](const T &value) {
        if constexpr (std::is_floating_point<T>::value) {
            os << std::fixed << std::setprecision(1) << value;
        } else {
            os << value;
        }
    };

    os << "{{";
    formatOutput(m.data[0]);
    os << ", ";
    formatOutput(m.data[1]);
    os << "}, {";
    formatOutput(m.data[2]);
    os << ", ";
    formatOutput(m.data[3]);
    os << "}}";

    return os;
}

// Lagrange mutipliers used in the lagrange reduction algoritm.
// Homogeneous nucleation of dislocations as a pattern formation phenomenon - R. Baggio

// Lagrange multiplier 1
template <typename T>
void Matrix2x2<T>::lag_m1()
{
    // Multiply by 1  0
    //             0 -1
    flip(0, 1);
    flip(1, 1);
}

// Lagrange multiplier 2
template <typename T>
void Matrix2x2<T>::lag_m2()
{
    // Multiply by 0 1
    //             1 0
    swapCols();
}

// Lagrange multiplier 3
template <typename T>
void Matrix2x2<T>::lag_m3(int n)
{
    // Multiply by 1 -1
    //             0  1

    // (*this) -> data[]
    (*this) = (*this) * Matrix2x2<T>{{{static_cast<T>(1), static_cast<T>(-n)},
                                      {static_cast<T>(0), static_cast<T>(1)}}};
}
