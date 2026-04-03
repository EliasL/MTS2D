#include "reduction.h"
#include <algorithm>
#include <cassert>
#include <cmath>

using Eigen::Matrix2d;

// Right-multiply by U_m = [[1, m],[0,1]]: col1 += m*col0
inline void mul_U(Matrix2d &m, int n) {
  m(0, 1) += n * m(0, 0);
  m(1, 1) += n * m(1, 0);
}

// Right-multiply by V_m = [[1,0],[m,1]]: col0 += m*col1
inline void mul_V(Matrix2d &m, int n) {
  m(0, 0) += n * m(0, 1);
  m(1, 0) += n * m(1, 1);
}

// Right-multiply by m1 = diag(1,-1): col1 -> -col1
inline void mul_m1(Matrix2d &m) {
  m(0, 1) = -m(0, 1);
  m(1, 1) = -m(1, 1);
}

// Right-multiply by m2 = [[0,1],[1,0]]: swap columns
inline void mul_m2(Matrix2d &m) {
  std::swap(m(0, 0), m(0, 1));
  std::swap(m(1, 0), m(1, 1));
}
inline void mul_m3(Matrix2d &m, int n = -1) { mul_U(m, n); }

inline bool in_elastic_domain(double a, double b, double c) {
  const double cmin = (a < c) ? a : c;
  const double abs_b = (b < 0.0) ? -b : b;
  return (cmin > 0.0) && (abs_b <= 0.5 * cmin);
}

int elastic_quadrant(double a, double b, double c) {
  if (!in_elastic_domain(a, b, c)) {
    return 0;
  }
  if (a < c) {
    return (b >= 0.0) ? 1 : 2;
  }
  return (b >= 0.0) ? 3 : 4;
}

inline void elastic_to_fundamental(double &a, double &b, double &c,
                                   Matrix2d &m) {
  if (b < 0.0) {
    b = -b;
    mul_m1(m);
  }
  if (c < a) {
    std::swap(a, c);
    mul_m2(m);
  }
}

/**
 * Elastic reduction for a 2x2 symmetric metric in-place.
 *
 * @param C_R           Output reduced metric (also used as the working matrix).
 * @param C_in          Input metric to start from (read-only).
 * @param m             Accumulated right-multiply matrix (updated in-place).
 * @param m3Nr           Counter for shear ops.
 * @param q             Elastic-domain quadrant label (1..4), or 0 if outside.
 * @param maxLoops      Safety cap on iterations.
 * @param fullReduction If true, apply final m1/m2 to map into fundamental
 * domain.
 * @return              true if converged before maxLoops; false otherwise.
 */
bool elasticReduction(Matrix2d &C_R, const Matrix2d &C_in, Matrix2d &m,
                      int &m3Nr, int &q, int maxLoops, bool fullReduction) {
  C_R = C_in;
  m.setIdentity();
  m3Nr = 0;
  q = 0;

  double &a = C_R(0, 0);
  double &b = C_R(0, 1);
  double &c = C_R(1, 1);

  // Fast-path: already in elastic domain before any iterations.
  if (in_elastic_domain(a, b, c)) {
    q = elastic_quadrant(a, b, c);
    if (fullReduction) {
      elastic_to_fundamental(a, b, c, m);
    }
    C_R(1, 0) = C_R(0, 1);
    return true;
  }

  bool converged = false;
  for (int iter = 0; iter < maxLoops; ++iter) {
    const bool useU = (a < c);
    const double denom = useU ? a : c; // == min(a, c)
    const double abs_b = (b < 0.0) ? -b : b;
    if ((denom > 0.0) && (abs_b <= 0.5 * denom)) {
      converged = true;
      break;
    }
    if (b == 0.0) {
      break; // no progress possible
    }
    // const int n = sgn(-b / denom);
    // const int n = -static_cast<int>(b / denom + 0.5);
    const double x = -b / denom;
    const int n =
        (x >= 0.0) ? static_cast<int>(x + 0.5) : static_cast<int>(x - 0.5);
    if (n == 0) {
      break; // no progress possible
    }

    // It might seem like we could use abs(n) larger than 1 to do multiple steps
    // at once, but none of round, floor or ceil work.
    const double old_b = b;
    if (useU) {
      // U_n: a' = a, b' = b + n a, c' = c + 2 n b + n^2 a
      const double n_a = n * a;
      b = old_b + n_a;
      c = c + 2.0 * n * old_b + n * n_a;
      mul_U(m, n);
    } else {
      // V_n: c' = c, b' = b + n c, a' = a + 2 n b + n^2 c
      const double n_c = n * c;
      b = old_b + n_c;
      a = a + 2.0 * n * old_b + n * n_c;
      mul_V(m, n);
    }
    m3Nr += (n < 0) ? -n : n;
  }

  q = elastic_quadrant(a, b, c);

  if (fullReduction) {
    elastic_to_fundamental(a, b, c, m);
  }

  C_R(1, 0) = C_R(0, 1); // enforce symmetry
  return converged;
}

/**
 * Lagrange reduction for a 2x2 symmetric metric in-place.
 * All arguments are passed by reference; no dynamic allocations.
 *
 * @param C_R       Output reduced metric (also used as the working matrix).
 * @param C_in      Input metric to start from (read-only).
 * @param m         Accumulated integer-transform matrix (updated in-place).
 * @param m3Nr      Counter for op3 (+= by ref).
 * @param q         Elastic-domain quadrant label (1..4), or 0 if outside.
 * @param maxLoops  Safety cap on iterations.
 * @param eps       Tolerance to avoid numerical chattering.
 * @return          true if converged before maxLoops; false otherwise.
 */
bool lagrangeReduction(Matrix2d &C_R,        // work/output: [[a,b],[b,c]]
                       const Matrix2d &C_in, // input metric
                       Matrix2d &m,          // accumulated transform
                       int &m3Nr, int &q, int maxLoops) {
  return elasticReduction(C_R, C_in, m, m3Nr, q, maxLoops, true);
}

/**
 * Lagrange reduction for a 2x2 symmetric metric in-place.
 * All arguments are passed by reference; no dynamic allocations.
 *
 * @param C_R       Output reduced metric (also used as the working matrix).
 * @param C_in      Input metric to start from (read-only).
 * @param m         Accumulated integer-transform matrix (updated in-place).
 * @param m3Nr      Counter for op3 (+= by ref).
 * @param maxLoops  Safety cap on iterations.
 * @param eps       Tolerance to avoid numerical chattering.
 * @return          true if converged before maxLoops; false otherwise.
 */
bool legacy_lagrangeReduction(Matrix2d &C_R, // work/output: [[a,b],[b,c]]
                              const Matrix2d &C_in, // input metric
                              Matrix2d &m,          // accumulated transform
                              int &m3Nr, int maxLoops) {
  C_R = C_in;
  m.setIdentity();
  m3Nr = 0;
  double &a = C_R(0, 0);
  double &b = C_R(0, 1);
  double &c = C_R(1, 1);
  assert(a > 0 && c > 0);

  for (int iter = 0; iter < maxLoops; ++iter) {
    bool changed = false;

    // 1) make off-diagonal non-negative
    if (std::signbit(b)) {
      b = -b;
      mul_m1(m);
      changed = true;
    }

    // 2) enforce a <= c
    if (c < a) {
      std::swap(a, c);
      mul_m2(m);
      changed = true;
    }

    // 3) single-shot shear: choose n so that -a/2 <= b + n a < a/2
    // Traditionally, n is always -1, but here, we avoid repeated m3 steps
    // Instead of looping several times and doing b += -a each time,
    // we calculate how many times we would need to do that, and do it all at
    // once.
    const double two_b = b + b;
    if (two_b > a) {
      // b >= 0 here, so this is a fast nearest-integer round (ties away from 0)
      const int n = -static_cast<int>(b / a + 0.5);
      if (n != 0) {
        const double n_a = n * a;
        c += n * (two_b + n_a); // uses old b
        b += n_a;
        mul_m3(m, n);
        m3Nr += (n < 0) ? -n : n;
        changed = true;
      }
    }

    if (!changed) {
      C_R(1, 0) = C_R(0, 1); // enforce symmetry once at the end
      return true;
    }
  }

  C_R(1, 0) = C_R(0, 1);
  return false; // hit loop cap (shouldn't happen in practice)
}
