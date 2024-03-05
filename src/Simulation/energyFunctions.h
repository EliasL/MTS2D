#ifndef ENERGYFUNCTIONS_H
#define ENERGYFUNCTIONS_H
#pragma once

#include <cmath>
#include "Matrix/matrix2x2.h"
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
 *        Φᵥ(det(Ĉ)) = μ(detĈ - log(detĈ))
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
double polynomialEnergy(double c11, double c22, double c12, double beta, double mu);

/*
 * =====================================================================================
 *        Discontinuous yielding of pristine micro-crystals, page 15
 *        Stress tensor function S(C) derivation from the energy density function Φ(C):
 *
 *                  Σ = ∂Φ/∂C = [  ∂Φ/∂C₁₁   1/2 ∂Φ/∂C₁₂]
 *                              [1/2 ∂Φ/∂C₁₂   ∂Φ/∂C₂₂  ]
 * =====================================================================================
 */
Matrix2x2<double> polynomialStress(double c11, double c22, double c12, double beta, double mu);
#endif