#include <vector>
#include <iostream>

template <typename T>
class Matrix {
public:
    std::vector<T> data;
    int rows;
    int cols;
    Matrix(int rows, int cols) : data(rows*cols), cols(cols), rows(rows) {}

    // This is complicated: https://stackoverflow.com/questions/36123452/statically-declared-2-d-array-c-as-this-member-of-a-class/36123944#36123944
    // It essentially returns a pointer to, not the beginning of the array,
    // but the array row*col elements in. 
    T* operator[](int row) { return &data[row*cols]; }

      // Const version of operator[]
    const T* operator[](int row) const {
        return &data[row * cols];
    }

    // Check if two matrices are equal
    bool operator==(const Matrix<T>& other) const {
        if (rows != other.rows || cols != other.cols) {
            return false;
        }
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if ((*this)[i][j] != other[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    // Fill the matrix with a specific value
    void fill(T value) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                (*this)[i][j] = value;
            }
        }
    }

    // Operator+ to add two matrices
    Matrix<T> operator+(const Matrix<T>& other) const {
        Matrix<T> result(rows, cols); // Create a new Matrix for the result

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result[i][j] = (*this)[i][j] + other[i][j];
            }
        }

        return result; // Return the result Matrix
    }

    // Overload * operator for matrix multiplication
    Matrix<T> operator*(const Matrix<T>& other) const {
        if (cols != other.rows) {
            throw std::invalid_argument("Number of columns in the first matrix must match the number of rows in the second matrix for multiplication.");
        }

        int resultRows = rows;
        int resultCols = other.cols;
        Matrix<T> result(resultRows, resultCols);

        for (int i = 0; i < resultRows; ++i) {
            for (int j = 0; j < resultCols; ++j) {
                result[i][j] = 0;
                for (int k = 0; k < cols; ++k) {
                    result[i][j] += (*this)[i][k] * other[k][j];
                }
            }
        }
        return result;
    }

    // Recursive function to find the determinant of a square matrix
    T det() const {
        if (rows != cols) {
            throw std::invalid_argument("Determinant can only be calculated for square matrices.");
        }

        if (rows == 1) {
            return (*this)[0][0];
        } else if (rows == 2) {
            return (*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0];
        } else {
            T det = 0;
            for (int i = 0; i < cols; ++i) {
                Matrix<T> submatrix(rows - 1, cols - 1);
                for (int j = 1; j < rows; ++j) {
                    for (int k = 0; k < cols; ++k) {
                        if (k < i) {
                            submatrix[j - 1][k] = (*this)[j][k];
                        } else if (k > i) {
                            submatrix[j - 1][k - 1] = (*this)[j][k];
                        }
                    }
                }
                det += (i % 2 == 0 ? 1 : -1) * (*this)[0][i] * submatrix.det();
            }
            return det;
        }
    }
};

// Overload << operator for Matrix
template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix) {
    for (int i = 0; i < matrix.rows; ++i) {
        for (int j = 0; j < matrix.cols; ++j) {
            os << matrix.data[i + j * matrix.cols] << "\t";
        }
        os << std::endl;
    }
    return os;
}
