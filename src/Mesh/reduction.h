#include <Eigen/Dense>
using Eigen::Matrix2d;

/**
 * Return elastic-domain quadrant label (1..4), or 0 if outside.
 */
int elastic_quadrant(double a, double b, double c);
inline int elastic_quadrant(const Matrix2d &C) {
  return elastic_quadrant(C(0, 0), C(0, 1), C(1, 1));
}

/**
 * Elastic reduction for a 2x2 symmetric metric in-place.
 *
 * @param C_R           Output reduced metric (also used as the working matrix).
 * @param C_in          Input metric to start from (read-only).
 * @param m             Accumulated right-multiply matrix (updated in-place).
 * @param m1Nr,m2Nr,m3Nr Counters for final m1/m2 ops and shear ops.
 * @param q             Elastic-domain quadrant label (1..4), or 0 if outside.
 * @param maxLoops      Safety cap on iterations.
 * @param fullReduction If true, apply final m1/m2 to map into fundamental
 * domain.
 * @return              true if converged before maxLoops; false otherwise.
 */
bool elasticReduction(Matrix2d &C_R, const Matrix2d &C_in, Matrix2d &m,
                      int &m1Nr, int &m2Nr, int &m3Nr, int &q,
                      int maxLoops = 100'000, bool fullReduction = false);

/**
 * Lagrange reduction for a 2x2 symmetric metric in-place.
 * All arguments are passed by reference; no dynamic allocations.
 *
 * @param C_R       Output reduced metric (also used as the working matrix).
 * @param C_in      Input metric to start from (read-only).
 * @param m         Accumulated integer-transform matrix (updated in-place).
 * @param m1Nr      Counter for op1 (+= by ref).
 * @param m2Nr      Counter for op2 (+= by ref).
 * @param m3Nr      Counter for op3 (+= by ref).
 * @param q         Elastic-domain quadrant label (1..4), or 0 if outside.
 * @param maxLoops  Safety cap on iterations.
 * @param eps       Tolerance to avoid numerical chattering.
 * @return          true if converged before maxLoops; false otherwise.
 */
bool lagrangeReduction(Matrix2d &C_R,        // work/output: [[a,b],[b,c]]
                       const Matrix2d &C_in, // input metric
                       Matrix2d &m,          // accumulated transform
                       int &m1Nr, int &m2Nr, int &m3Nr, int &q,
                       int maxLoops = 100'000);

/**
 * Lagrange reduction for a 2x2 symmetric metric in-place.
 * All arguments are passed by reference; no dynamic allocations.
 *
 * @param C_R       Output reduced metric (also used as the working matrix).
 * @param C_in      Input metric to start from (read-only).
 * @param m         Accumulated integer-transform matrix (updated in-place).
 * @param m1Nr      Counter for op1 (+= by ref).
 * @param m2Nr      Counter for op2 (+= by ref).
 * @param m3Nr      Counter for op3 (+= by ref).
 * @param maxLoops  Safety cap on iterations.
 * @param eps       Tolerance to avoid numerical chattering.
 * @return          true if converged before maxLoops; false otherwise.
 */
bool legacy_lagrangeReduction(Matrix2d &C_R, // work/output: [[a,b],[b,c]]
                              const Matrix2d &C_in, // input metric
                              Matrix2d &m,          // accumulated transform
                              int &m1Nr, int &m2Nr, int &m3Nr,
                              int maxLoops = 100'000);
