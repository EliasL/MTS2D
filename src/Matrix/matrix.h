#ifndef MATRIX_H
#define MATRIX_H
#pragma once

#include <vector>
#include <iostream>
#include <type_traits> // For static_assert

#include <cereal/types/vector.hpp> // Cereal serialization for std::vector

template <typename T>
class Matrix
{
public:
    // Constructors
    Matrix();
    Matrix(int rows, int cols);

    // This is complicated: https://stackoverflow.com/questions/36123452/statically-declared-2-d-array-c-as-this-member-of-a-class/36123944#36123944
    // It essentially returns a pointer to, not the beginning of the array,
    // but the array row*col elements in.
    T *operator[](int row) { return &data[row * cols]; }
    // Const version of operator[]
    const T *operator[](int row) const { return &data[row * cols]; }

    void fill(T value);
    bool operator==(const Matrix<T> &other) const;
    Matrix<T> operator+(const Matrix<T> &other) const;
    Matrix<T> operator*(const Matrix<T> &other) const;
    Matrix<T> mat_mult(const Matrix<T> &other) const;
    T det() const;
    void transpose();

    // Members
    std::vector<T> data;
    int rows;
    int cols;

    template <class Archive>
    void serialize(Archive &ar)
    {
        ar(data, rows, cols);
    }
};

// Non-member functions
template <typename T>
Matrix<T> transpose(const Matrix<T> &matrix);

template <typename T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &matrix);

#include "matrix.tpp"

// Explicit instantiation declarations
// extern template class Matrix<int>;
// extern template class Matrix<double>;

#endif // MATRIX_H
