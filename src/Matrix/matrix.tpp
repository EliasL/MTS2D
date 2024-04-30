#include "matrix.h"

template <typename T>
Matrix<T>::Matrix() {}

template <typename T>
Matrix<T>::Matrix(int rows, int cols) : data(rows * cols), rows(rows), cols(cols)
{
    if (rows <= 0 || cols <= 0)
    {
        throw std::invalid_argument("Rows and cols must be greater than 0.");
    }
}

template <typename T>
// Check if two matrices are equal
bool Matrix<T>::operator==(const Matrix<T> &other) const
{
    if (rows != other.rows || cols != other.cols)
    {
        return false;
    }
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            if ((*this)[i][j] != other[i][j])
            {
                return false;
            }
        }
    }
    return true;
}

// Fill the matrix with a specific value
template <typename T>
void Matrix<T>::fill(T value)
{
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            (*this)[i][j] = value;
        }
    }
}

// Operator+ to add two matrices
template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &other) const
{
    Matrix<T> result(rows, cols); // Create a new Matrix for the result

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result[i][j] = (*this)[i][j] + other[i][j];
        }
    }

    return result; // Return the result Matrix
}

// Operator* to ELEMENTWISE MULTIPLICATION, not matrix multiplication
template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &other) const
{
    if (cols != other.cols || rows != other.rows)
    {
        throw std::invalid_argument("Size of matrices do not match.");
    }
    Matrix<T> result(rows, cols); // Create a new Matrix for the result

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result[i][j] = (*this)[i][j] * other[i][j];
        }
    }

    return result; // Return the result Matrix
}

// Define matrix multiplication
template <typename T>
Matrix<T> Matrix<T>::mat_mult(const Matrix<T> &other) const
{
    if (cols != other.rows)
    {
        throw std::invalid_argument("Number of columns in the first matrix must match the number of rows in the second matrix for multiplication.");
    }

    int resultRows = rows;
    int resultCols = other.cols;
    Matrix<T> result(resultRows, resultCols);

    for (int i = 0; i < resultRows; ++i)
    {
        for (int j = 0; j < resultCols; ++j)
        {
            result[i][j] = 0;
            for (int k = 0; k < cols; ++k)
            {
                result[i][j] += (*this)[i][k] * other[k][j];
            }
        }
    }
    return result;
}

// Recursive function to find the determinant of a square matrix
template <typename T>
T Matrix<T>::det() const
{
    if (rows != cols)
    {
        throw std::invalid_argument("Determinant can only be calculated for square matrices.");
    }

    if (rows == 1)
    {
        return (*this)[0][0];
    }
    else if (rows == 2)
    {
        return (*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0];
    }
    else
    {
        T det = 0;
        for (int i = 0; i < cols; ++i)
        {
            Matrix<T> submatrix(rows - 1, cols - 1);
            for (int j = 1; j < rows; ++j)
            {
                for (int k = 0; k < cols; ++k)
                {
                    if (k < i)
                    {
                        submatrix[j - 1][k] = (*this)[j][k];
                    }
                    else if (k > i)
                    {
                        submatrix[j - 1][k - 1] = (*this)[j][k];
                    }
                }
            }
            det += (i % 2 == 0 ? 1 : -1) * (*this)[0][i] * submatrix.det();
        }
        return det;
    }
}

template <typename T>
void Matrix<T>::transpose()
{
    int oldRows = rows;
    int oldCols = cols;

    // Swap rows and columns
    rows = oldCols;
    cols = oldRows;

    std::vector<T> newData(rows * cols);

    for (int i = 0; i < oldRows; ++i)
    {
        for (int j = 0; j < oldCols; ++j)
        {
            newData[j * oldRows + i] = data[i * oldCols + j];
        }
    }
    data = newData;
}

template <typename T>
Matrix<T> transpose(const Matrix<T> &matrix)
{
    int rows = matrix.cols;
    int cols = matrix.rows;

    Matrix<T> result(cols, rows);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result[j][i] = matrix[i][j];
        }
    }

    return result;
}

// Overload << operator for Matrix
template <typename T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &matrix)
{
    Matrix<T> m = transpose(matrix);
    os << std::endl;
    for (int i = m.rows - 1; i >= 0; --i)
    {
        for (int j = 0; j < m.cols; ++j)
        {
            os << m[i][j] << "\t";
        }
        os << std::endl;
    }
    return os;
}