#ifndef ENERGYFUNCTIONS_H
#define ENERGYFUNCTIONS_H
#pragma once

#include "itensor/itensor.h"
#include <Eigen/Core>

using namespace Eigen;

namespace ContiPotential {
using namespace itensor;
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
double energyDensity(double c11, double c22, double c12, double beta, double K,
                     double noise = 0);

/*
 * =====================================================================================
 *        Discontinuous yielding of pristine micro-crystals, page 15
 *        Stress tensor function S(C) derivation from the energy density
 * function Φ(C):
 *
 *                  Σ = ∂Φ/∂C = [  ∂Φ/∂C₁₁   1/2 ∂Φ/∂C₁₂]
 *                              [1/2 ∂Φ/∂C₁₂   ∂Φ/∂C₂₂  ]
 * =====================================================================================
 */
Matrix2d stress(double c11, double c22, double c12, double beta, double K);

/*
 * =====================================================================================
 *        Discontinuous yielding of pristine micro-crystals, page 15
 *        Hessian tensor function H(C) derivation from the energy density
 * function Φ(C):
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
ITensor hessian(double c11, double c22, double c12, double beta, double K);

} // namespace ContiPotential
#endif