#ifndef MATRIX2X2_H
#define MATRIX2X2_H
#pragma once

#include <stdexcept>
#include <cmath>
#include <array>

template <typename T>
class Matrix2x2
{
public:
    std::array<T, 4> data;

    Matrix2x2(T a);
    Matrix2x2();
    Matrix2x2(T c11, T c12, T c21, T c22);
    Matrix2x2(std::array<T, 4> d);

    // This is complicated: https://stackoverflow.com/questions/36123452/statically-declared-2-d-array-c-as-(*this)-member-of-a-class/36123944#36123944
    // It essentially returns a pointer to, not the beginning of the array,
    // but the array row*col elements in.
    T *operator[](int row) { return &data[row * 2]; }

    // Const version of operator[]
    const T *operator[](int row) const { return &data[row * 2]; }

    void flip(int x, int y);
    void swap(int x1, int y1, int x2, int y2);
    void swapCols();

    Matrix2x2 operator*(const Matrix2x2 &other) const;
    bool operator==(const Matrix2x2 &other) const;

    T det() const;
    Matrix2x2 operator+(const Matrix2x2 &other) const;
    Matrix2x2 transpose() const;
    Matrix2x2 inverse() const;
    Matrix2x2 conjugate(const Matrix2x2 &g) const;
    bool isOrthogonal() const;
    Matrix2x2 orth_conjugate(const Matrix2x2 &g) const;
    Matrix2x2 sym_orth_conjugate(const Matrix2x2 &g) const;

    static Matrix2x2 rotation_matrix(double theta);
    static Matrix2x2 reflection_matrix(double theta);

    void lag_m1();
    void lag_m2();
    void lag_m3();
};

#include "matrix2x2.tpp"

//extern template class Matrix2x2<int>;
//extern template class Matrix2x2<double>;

#endif // MATRIX2X2_H
