/***********************************************************************************************
 * Copyright (C) 2018 Dillon Cislo
 *
 * This file is part of FIRE++.
 *
 * FIRE++ is free software: you can redistribute it and/or modify it under the
 *terms of the GNU General Public License as published by the Free Software
 *Foundation, either version 3 of the License, or (at your option) any later
 *version.
 *
 * This program is distributed in the hope that it will by useful, but WITHOUT
 *ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or FITNESS
 *FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 *this program. If not, see <http://www.gnu.org/licenses/>
 *
 ***********************************************************************************************/

/**
 * 	\file FIRE.h
 * 	\brief A FIRE solver for numerical optimization
 *
 * 	I should write a detailed explanation at some point.
 *
 * 	\author Dillon Cislo
 * 	\date 11/16/2018
 * 	\copyright GNU Public License
 *
 **/

#ifndef _FIRE_H_
#define _FIRE_H_

#include <cfloat>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <stdexcept>

#include <Eigen/Core>

#include "MolecularDynamics.h"
#include "Param.h"

namespace FIREpp {

///
/// FIRE solver for unconstrained numerical optimization
///
template <typename Scalar> class FIRESolver {

private:
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

  const FIREParam<Scalar> &m_param; // Parameters to control the FIRE algorithm

  Vector m_fx;     // History of objective function values
  Vector m_xp;     // Old x
  Vector m_xStart; // Initial x
  Vector m_grad;   // New gradient
  Vector m_gradp;  // Old gradient
  Vector m_v;      // New velocity
  Vector m_vp;     // Old velocity
  Vector m_mass;   // Particle masses
  bool m_userMass; // Check whether the user supplied masses

  inline void reset(int n) {

    m_xp.resize(n);
    m_xStart.resize(n);
    m_grad.resize(n);
    m_gradp.resize(n);

    m_v.resize(n);
    m_v = Vector::Zero(n);
    m_vp.resize(n);
    m_vp = Vector::Zero(n);

    if (m_param.past > 0)
      m_fx.resize(m_param.past);
  };

public:
  ///
  /// Constructor for FIRE solver
  ///
  /// \param param An object of \ref FIREParam to store parameters
  /// 	for the algorithm
  ///
  FIRESolver(const FIREParam<Scalar> &param) : m_param(param) {
    m_param.check_param();
    m_userMass = false;
  };

  ///
  /// Constructor for FIRE solver with mass
  ///
  /// \param param An object of \ref FIREParam to store parameters
  /// 	for the algorithm
  /// \param mass A vector of particle masses
  ///
  FIRESolver(const FIREParam<Scalar> &param, Vector mass)
      : m_param(param), m_mass(mass) {

    m_param.check_param();
    m_userMass = true;
  };

  ///
  /// Minimizing a multivariate function using the FIRE algorithm
  /// Exceptions will be thrown if error occurs.
  ///
  /// \param f	A function object such that `f(x, grad)` returns the
  ///		objective function value at `x`, and overwrites `grad`
  /// 		with the gradient.
  /// \param x 	IN: An initial guress of the optimal point.
  ///		OUT: The best point found
  /// \param fx 	OUT: The objective function value at `x`.
  ///
  /// \param optPtr IN: A pointer to an object passed to the f function
  /// 		Use reinterpret_cast to use as any object.
  ///
  /// \param termT OUT: An integer that is set to the termination type
  ///     Numbers are chosen to match with LBFGS from alglib
  ///     1: Gradient norm less than relative epsilon
  ///     2: Objective function value has small change
  ///     4: Gradient norm less than epsilon
  ///     3: Guess is good enough
  ///    -3: Energy is too high. Solution is not likely to be found
  /// \return 	Number of iterations used
  ///
  template <typename Foo>
  inline int minimize(Foo &f, Vector &x, Scalar &fx, void *optPtr, int &termT) {

    const int n = x.size();
    const int fpast = m_param.past;
    // Save the initial point so we can revert if minimization fails.
    m_xStart.noalias() = x;

    reset(n);
    if (!m_userMass) {
      m_mass.resize(n);
      m_mass = Vector::Constant(n, 1, Scalar(1.0));
    }

    // Evaluate function and compute gradient
    fx = f(x, m_grad, optPtr);
    Scalar xnorm = x.norm();
    Scalar gnorm = m_grad.norm();
    if (fpast > 0)
      m_fx[0] = fx;

    int k = 0;

    // Display iterative updates
    if (m_param.iter_display) {

      std::cout << "( 0)"
                << " ||dx|| = " << std::scientific << std::setw(10)
                << std::setprecision(5) << gnorm
                << " ||x|| = " << std::setprecision(5) << std::setw(10) << xnorm
                << " f(x) = " << std::setw(10) << std::setprecision(5) << fx
                << std::endl;
    }

    // Handle NaNs produced by the initial guess
    if ((fx != fx) || (xnorm != xnorm) || (gnorm != gnorm))
      throw std::invalid_argument("Initial guess generates NaNs");

    // Handle Infs produced by the initial guess
    if ((std::isinf(fx)) || std::isinf(xnorm) || std::isinf(gnorm))
      throw std::invalid_argument("Initial guess generates Infs");

    // Early exit if the initial x is already a minimizer
    if (gnorm <= m_param.epsilon || gnorm <= m_param.epsilon_rel * xnorm) {

      if (m_param.iter_display)
        std::cout << "EARLY EXIT CONDITION: Gradient Norm" << std::endl;
      termT = 3;
      return 1;
    }

    k++;
    int downhillStepsInARow = 0;
    int uphillStepsInARow = 0;
    Scalar dt = m_param.dt_start;
    Scalar alpha = m_param.alpha_start;

    while (true) {

      // Save the current x, v, and gradient
      m_xp.noalias() = x;
      m_gradp.noalias() = m_grad;
      m_vp.noalias() = m_v;

      // Molecular dynamics update to current system configuration
      MolecularDynamics<Scalar>::VelocityVerlet(f, fx, x, m_v, m_grad, dt,
                                                alpha, m_mass, m_param, optPtr);

      // New x norm and gradient norm
      xnorm = x.norm();
      gnorm = m_grad.norm();

      // Handle NaNs produced by the current iterate
      if ((fx != fx) || (xnorm != xnorm) || (gnorm != gnorm))
        throw std::invalid_argument("Current iterate generates NaNs");

      // Handle Infs produced by the current interate
      if ((std::isinf(fx)) || std::isinf(xnorm) || std::isinf(gnorm))
        throw std::invalid_argument("Current iterate generates Infs");

      // Display iterative updates
      if (m_param.iter_display) {

        std::cout << "(" << std::setw(2) << k << ")"
                  << " ||dx|| = " << std::scientific << std::setw(10)
                  << std::setprecision(5) << gnorm
                  << " ||x|| = " << std::setprecision(5) << std::setw(10)
                  << xnorm << " f(x) = " << std::setw(10)
                  << std::setprecision(5) << fx << std::endl;
      }

      // Convergence test -- gradient
      if (gnorm <= m_param.epsilon || gnorm <= m_param.epsilon_rel * xnorm) {

        if (m_param.iter_display)
          std::cout << "CONVERGENCE CRITERION: Gradient Norm" << std::endl;
        termT = 1;
        return k;
      }

      // Convergence test -- objective function value
      if (fpast > 0) {

        if ((k >= fpast) &&
            (std::abs((m_fx[k % fpast] - fx) / fx) < m_param.delta)) {

          if (m_param.iter_display)
            std::cout << "CONVERGENCE CRITERION: Objective Function Value"
                      << std::endl;
          termT = 2;
          return k;
        }

        m_fx[k % fpast] = fx;
      }

      // Maximum number of iterations
      if (m_param.max_iterations != 0 && k >= m_param.max_iterations) {

        if (m_param.iter_display)
          std::cout << "CONVERGENCE CRITERION: Maximum Iteration Count"
                    << std::endl;
        termT = 4;
        return k;
      }

      // Abort simulation
      if (k >= 10 && fx > 3e3) { // Expected energy would be around 0.1
        // std::runtime_error("Energy diverging.");
        std::cout << "CONVERGENCE FAILED: Energy is too high: " << fx
                  << std::endl;
        termT = -3;

        // We want to go back to how things were before
        f(m_xStart, m_grad, optPtr);

        return k;
      }

      if (uphillStepsInARow > m_param.max_uphillSteps) {
        std::cout << "CONVERGENCE FAILED: Always uphill" << std::endl;
        termT = 7;
        return k;
      }

      // This indicates whether or not the velocity is pointing in the same
      // direction as the gradient
      Scalar P = m_v.dot(-m_grad);

      m_v = (Scalar(1) - alpha) * m_v +
            alpha * m_v.norm() * (-m_grad) / m_grad.norm();

      // If the current velocity is 'downhill' (good)
      if (P > Scalar(0)) {

        uphillStepsInARow = 0;
        downhillStepsInARow++;

        if (downhillStepsInARow > m_param.nmin) {
          alpha = alpha * m_param.falpha;
          dt = std::min(dt * m_param.finc, m_param.dt_max);
        }

        // If the current velocity is 'uphill' (bad)
      } else {

        downhillStepsInARow = 0;
        uphillStepsInARow++;
        // We don't decrease the timestep if we have just started minimizing, we
        // will make a few more attempts first
        if (k > m_param.nmin) {
          dt = std::max(dt * m_param.fdec, m_param.dt_min);
          alpha = m_param.alpha_start;
        }
        // If we are moving uphill, we should take half a step back before
        // trying again
        x -= Scalar(0.5) * dt * m_v;
        m_v = Vector::Zero(n);
      }

      k++;
    }

    return k;
  };
};

} // namespace FIREpp

#endif // _FIRE_H_
