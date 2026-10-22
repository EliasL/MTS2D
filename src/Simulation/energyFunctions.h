#ifndef ENERGYFUNCTIONS_H
#define ENERGYFUNCTIONS_H
#pragma once

// #include "itensor/itensor.h"
#include <Eigen/Core>
#include <cmath>

using namespace Eigen;

namespace ContiPotential {
/*
 * =====================================================================================
 *        Discontinuous yielding of pristine micro-crystals, page 10
 *        Energy density function Φ(C) is defined as:
 *
 *                  Φ(C) = Φ₀_d(Ĉ/sqrt(detC)) + Φᵥ(detĈ)
 *
 *        where C is the metric tensor and
 *        Ĉ = MᵀCM is the lagrange reduced metric tensor.
 *
 *        Φᵥ(det(Ĉ)) = K(detĈ - log(detĈ))
 *        Φ₀_d(Ĉ) = βψ₁(Ĉ) + ψ₂(Ĉ)
 *
 *        ψ₁ = I₁⁴I₂ - 41I₂³/99 + 7I₁I₂I₃/66 + I₃²/1056
 *        ψ₂ = 4I₂³/11 + I₁³I₃ - 8I₁I₂I₃/11 + 17I₃²/528
 *
 *        I₁ = 1/3 (Ĉ₁₁ + Ĉ₂₂ - Ĉ₁₂)
 *        I₂ = 1/4 (Ĉ₁₁ - Ĉ₂₂)² + 1/12 (Ĉ₁₁ + Ĉ₂₂ - 4Ĉ₁₂)²
 *        I₃ = (Ĉ₁₁ - Ĉ₂₂)²(Ĉ₁₁ + Ĉ₂₂ - 4Ĉ₁₂) - 1/9 (Ĉ₁₁ + Ĉ₂₂ - Ĉ₁₂)³
 *
 * =====================================================================================
 */
double energyDensity(double C11, double C22, double C12, double beta, double K,
                     double noise = 1);

// Shared energy definition for ordinary numbers and the tangent's derivative
// scalar. Unqualified math calls also find that scalar's pow/log overloads.
template <typename Scalar>
Scalar energyDensityFormula(const Scalar &C11, const Scalar &C22,
                            const Scalar &C12, double beta, double K,
                            double noise) {
  using std::log;
  using std::pow;

  const Scalar x0 = C11 * C22 - pow(C12, 2);
  const Scalar x1 = noise * x0;
  const Scalar x2 = pow(x0, -1.0 / 2.0);
  const Scalar x3 = C11 * x2;
  const Scalar x4 = C22 * x2;
  const Scalar x5 = pow(x3 - x4, 2);
  const Scalar x6 = C12 * x2;
  const Scalar x7 = x3 + x4 - 4 * x6;
  const Scalar x8 = 0.25 * x5 + 0.083333333333333329 * pow(x7, 2);
  const Scalar x9 = pow(x8, 3);
  const Scalar x10 = x5 * x7 - 0.1111111111111111 * pow(x7, 3);
  const Scalar x11 = pow(x10, 2);
  const Scalar x12 = 0.33333333333333331 * x3 + 0.33333333333333331 * x4 -
                     0.33333333333333331 * x6;
  const Scalar x13 = x10 * x8;

  const Scalar phi =
      K * (x1 - log(x1)) +
      beta * (0.000946969696969697 * x11 + pow(x12, 4) * x8 +
              0.035353535353535352 * x13 * (x3 + x4 - x6) -
              0.4141414141414142 * x9) +
      x10 * pow(x12, 3) + 0.032196969696969696 * x11 -
      0.24242424242424243 * x13 * (x3 + x4 - x6) +
      0.36363636363636365 * x9;

  return phi;
}

/*
 * =====================================================================================
 *        Discontinuous yielding of pristine micro-crystals, page 15
 *        Stress tensor function S(C) derivation from the energy density
 *        function Φ(C):
 *
 *                  Σ = ∂Φ/∂C = [  ∂Φ/∂C₁₁   1/2 ∂Φ/∂C₁₂]
 *                              [1/2 ∂Φ/∂C₁₂   ∂Φ/∂C₂₂  ]
 * =====================================================================================
 */
Matrix2d stress(double C11, double C22, double C12, double beta, double K,
                double noise = 1);

/*
 * =====================================================================================
 *        Combined energy density and stress calculation.
 *        Returns Φ(C) and optionally writes Σ = ∂Φ/∂C into stress_out.
 *        You might think this is faster, but my testing indicates that it's
 *        somehow significantly slower than calling energyDensity and stress
 *        separately.
 * =====================================================================================
 */
double energyDensityAndStress(double C11, double C22, double C12, double beta,
                              double K, double noise, Matrix2d *stress_out);
/*
 * =====================================================================================
 *        Discontinuous yielding of pristine micro-crystals, page 15
 *        Hessian tensor function H(C) derivation from the energy density
 *        function Φ(C):
 *
  [ [ ∂²Φ⁰/∂C₁₁²        1/2 ∂²Φ⁰/∂C₁₂∂C₁₁]  1[ ∂²Φ⁰/∂C₁₁∂C₁₂   1/2 ∂²Φ⁰/∂C²₁₂] ]
  | [ 1/2 ∂²Φ⁰/∂C₁₂∂C₁₁ ∂²Φ⁰/∂C₂₂∂C₁₁    ]  2[ 1/2 ∂²Φ⁰/∂C²₁₂  ∂²Φ⁰/∂C₂₂∂C₁₂ ] |
H=|                                                                            |
  | 1[ 1/2 ∂²Φ⁰/∂C₁₁∂C₁₂ 1/2 ∂²Φ⁰/∂C²₁₂] [ ∂²Φ⁰/∂C₁₁∂C₂₂    1/2 ∂²Φ⁰/∂C₁₂∂C₁₁] |
  [ 2[ 1/2 ∂²Φ⁰/∂C²₁₂    ∂²Φ⁰/∂C₂₂∂C₁₂ ] [ 1/2 ∂²Φ⁰/∂C₁₂∂C₁₁∂²Φ⁰/∂C₂²        ] ]
 *
 *        Note that the diagonal of H is symetric. We can therefore use the
 *        bottom right, element to store the first derivative and only return
 *        a single Matrix3d object.
 *
 * =====================================================================================
 */
// itensor::ITensor hessian(double c11, double c22, double c12, double beta,
//                          double K, double noise);

} // namespace ContiPotential
#endif
