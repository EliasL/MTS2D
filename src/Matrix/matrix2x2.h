#ifndef MATRIX2X2_H
#define MATRIX2X2_H
#pragma once

#include <stdexcept>
#include <cmath>
#include <array>

template <typename T>
class Matrix2x2
{
    /**
     * @brief A row, column accessed 2x2 matrix.
     *
     * A = {{1,2},
     *      {3,4}}
     *  -> A[0][1] returns 2, not 3
     *
     * To avoid confusion, be very careful about using x or y as indexing variables
     */
public:
    std::array<T, 4> data;

    Matrix2x2(T a);
    Matrix2x2();
    // Matrix2x2(T c11, T c12, T c21, T c22);
    Matrix2x2(std::array<T, 4> d);

    Matrix2x2(std::initializer_list<std::initializer_list<T>> list);

    // This is complicated: https://stackoverflow.com/questions/36123452/statically-declared-2-d-array-c-as-(*this)-member-of-a-class/36123944#36123944
    // It essentially returns a pointer to, not the beginning of the array,
    // but the array row*col elements in.
    T *operator[](int row) { return &data[row * 2]; }

    // Const version of operator[]
    const T *operator[](int row) const { return &data[row * 2]; }

    // flipp the sign of specified value
    void flip(int row, int col);
    // swap the possition of two values in the matrix
    void swap(int row1, int col1, int row2, int col2);
    // swap the two columns
    void swapCols();
    // set a column of the matrix
    void setCol(std::array<T, 2> column, int col);
    // set the matrix using two column vectors
    void setCols(std::array<T, 2> column1, std::array<T, 2> column2);

    Matrix2x2 operator*(const Matrix2x2 &other) const;
    Matrix2x2 operator*(T scalar) const;

    bool operator==(const Matrix2x2 &other) const;

    T det() const;
    Matrix2x2 operator+(const Matrix2x2 &other) const;
    Matrix2x2 transpose() const;
    Matrix2x2 inverse() const;
    Matrix2x2 similarityTransform(const Matrix2x2 &g) const;
    bool isOrthogonal() const;
    Matrix2x2 congruenceTransform(const Matrix2x2 &g) const;
    Matrix2x2 sym_orth_conjugate(const Matrix2x2 &g) const;

    static Matrix2x2 rotation_matrix(double theta);
    static Matrix2x2 reflection_matrix(double theta);

    void lag_m1();
    void lag_m2();
    void lag_m3();
};

template <typename T>
Matrix2x2<T> operator*(T scalar, const Matrix2x2<T>& matrix);

#include "matrix2x2.tpp"

// extern template class Matrix2x2<int>;
// extern template class Matrix2x2<double>;

#endif // MATRIX2X2_H
