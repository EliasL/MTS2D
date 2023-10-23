#ifndef MATRIX2X2_H
#define MATRIX2X2_H


#include <stdexcept>
#include <cmath>

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

    Matrix2x2(T c11, T c12, T c21, T c22) {
        (*this)[0][0] = c11;
        (*this)[0][1] = c12;
        (*this)[1][0] = c21;
        (*this)[1][1] = c22;
    }

    Matrix2x2(std::array<T, 4> d){
        data = d;
    }

    // Change the sign of the value in [x][y]
    void flip(int x, int y){
        (*this)[x][y] *= -1;
    }

    void swap(int x1, int y1, int x2, int y2){
        T temp;
        temp = (*this)[x1][y1];
        (*this)[x1][y1] = (*this)[x2][y2];
        (*this)[x2][y2] = temp;
    }

    void swapCols() {
        swap(0,0,0,1);
        swap(1,0,1,1);
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

    // Check if two matrixies are equal
    bool operator==(const Matrix2x2& other) const {
        for (size_t i = 0; i < 4; ++i) {
            if (data[i] != other.data[i]) {
                return false; // If any elements are not equal, return false
            }
        }
        return true; // All elements are equal
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
            throw std::runtime_error("Matrix is not invertible.");

        T invDet = 1 / detVal;
        Matrix2x2 result;
        result[0][0] = (*this)[1][1] * invDet;
        result[0][1] = -(*this)[0][1] * invDet;
        result[1][0] = -(*this)[1][0] * invDet;
        result[1][1] = (*this)[0][0] * invDet;
        return result;
    }

    // Calculates the conjugate of (this) matrix 
    Matrix2x2 conjugate(const Matrix2x2& g) const {
        return g * (*this) * g.inverse();
    }

    bool isOrthogonal() const {
        // First check determinant
        if(abs(abs(det())-1) > __DBL_EPSILON__){
            return false;
        }
        // Then do a more thurough test
        Matrix2x2 inv = Matrix2x2(data);
        Matrix2x2 tran = Matrix2x2(data);
        return tran.transpose() == inv.inverse();
    }
    
    // Calculate the conjugate provided that g is orthogonal
    Matrix2x2 orth_conjugate(const Matrix2x2& g) const {
        #if defined(__OPTIMIZE__) && __OPTIMIZE__ < 3
            if(!g.isOrthogonal()){
                throw std::runtime_error("Provided matrix is not orthogonal");
            }
        #endif

        return g * (*this) * g.transpose();
    }
    
    // Calculate the conjugate provided (this) matrix is symetric, and g is orthogonal
    Matrix2x2 sym_orth_conjugate(const Matrix2x2& g) const {
        #if __OPTIMIZE__ < 3
            if ((*this)[0][1] != (*this)[1][0]) {
                throw std::runtime_error("Matrix is not symmetric.");
            }
            if(!g.isOrthogonal()){
                throw std::runtime_error("Provided matrix is not orthogonal");
            }
        #endif

        return orth_conjugate(g);
        // This does not work TODO!
        Matrix2x2<double> d;
        d[0][0] = (*this)[0][0]*g[0][0]*g[0][0] + (*this)[0][1]*g[0][0]*g[0][1] + (*this)[1][1]*g[0][1]*g[0][1];
        d[1][1] = (*this)[0][0]*g[1][0]*g[1][0] + (*this)[0][1]*g[1][0]*g[1][1] + (*this)[1][1]*g[1][1]*g[1][1];
        d[0][1] = g[1][0]*(2 * (*this)[0][0]*g[0][0] + (*this)[0][1]*g[0][1]) + g[1][1]*((*this)[0][1]*g[0][0] + 2 * (*this)[1][1]*g[0][1]);
        d[1][0] = d[0][1];
        return d;
    }

    static Matrix2x2 rotation_matrix(double theta) {
        return Matrix2x2(cos(theta), -sin(theta), sin(theta), cos(theta));
    }
    static Matrix2x2 reflection_matrix(double theta) {
        return Matrix2x2(cos(theta), sin(theta), sin(theta), -cos(theta));
    }



// Lagrange mutipliers used in the lagrange reduction algoritm.
    // Lagrange multiplier 1
    void lag_m1(){
        // Multiply by 1  0
        //             0 -1
        flip(1,0);
        flip(1,1);
    }

    // Lagrange multiplier 2
    void lag_m2(){
        // Multiply by 0 1
        //             1 0
        swapCols();
    }

    // Lagrange multiplier 3
    void lag_m3(){
        // Multiply by 1 -1
        //             1  1

        // (*this) = {
        //     data[]
        // }
        (*this) = (*this) * Matrix2x2<T>{static_cast<T>(1), static_cast<T>(-1),
                                static_cast<T>(1), static_cast<T>(1)};
    }

};

#endif