#ifndef MATRIX2X2_H
#define MATRIX2X2_H


#include <stdexcept>

template <typename T>
class Matrix2x2 {
public:
    std::array<T, 4> data;
    Matrix2x2(T a) {
        data.fill(a);
    }
    // This is complicated: https://stackoverflow.com/questions/36123452/statically-declared-2-d-array-c-as-(*this)-member-of-a-class/36123944#36123944
    // It essentially returns a pointer to, not the beginning of the array,
    // but the array row*col elements in. 
    T* operator[](int row) { return &data[row*2]; }

      // Const version of operator[]
    const T* operator[](int row) const { return &data[row * 2]; }

    Matrix2x2() {
        (*this)[0][0] = (*this)[1][1] = static_cast<T>(1);
        (*this)[0][1] = (*this)[1][0] = static_cast<T>(0);
    }

    Matrix2x2(T a, T b, T c, T d) {
        (*this)[0][0] = a;
        (*this)[0][1] = b;
        (*this)[1][0] = c;
        (*this)[1][1] = d;
    }



    // Matrix multiplication optimized for 2x2 matrices
    Matrix2x2 operator*(const Matrix2x2& other) const {
        Matrix2x2 result;
        result[0][0] = (*this)[0][0] * other[0][0] + (*this)[0][1] * other[1][0];
        result[0][1] = (*this)[0][0] * other[0][1] + (*this)[0][1] * other[1][1];
        result[1][0] = (*this)[1][0] * other[0][0] + (*this)[1][1] * other[1][0];
        result[1][1] = (*this)[1][0] * other[0][1] + (*this)[1][1] * other[1][1];
        return result;
    }

    // Calculate the determinant of a 2x2 matrix
    T det() const {
        return (*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0];
    }

    // Matrix addition for 2x2 matrices
    Matrix2x2 operator+(const Matrix2x2& other) const {
        Matrix2x2 result;
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                result[i][j] = (*this)[i][j] + other[i][j];
            }
        }
        return result;
    }

    // Transpose the 2x2 matrix
    Matrix2x2 transpose() const {
        Matrix2x2 result;
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                result[i][j] = (*this)[j][i];
            }
        }
        return result;
    }

    // Calculate the inverse of the 2x2 matrix
    Matrix2x2 inverse() const {
        T detVal = det();
        if (detVal == 0)
            throw std::runtime_error("Matrix is not invertible");

        T invDet = 1 / detVal;
        Matrix2x2 result;
        result[0][0] = (*this)[1][1] * invDet;
        result[0][1] = -(*this)[0][1] * invDet;
        result[1][0] = -(*this)[1][0] * invDet;
        result[1][1] = (*this)[0][0] * invDet;
        return result;
    }
};

#endif