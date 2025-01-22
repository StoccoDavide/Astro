/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Astro project is distributed under the BSD 2-Clause License.                              *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef ASTRO_UTILITIES_HXX
#define ASTRO_UTILITIES_HXX

namespace Astro {

  /*\
   |   _   _ _   _ _ _ _   _
   |  | | | | |_(_) (_) |_(_) ___  ___
   |  | | | | __| | | | __| |/ _ \/ __|
   |  | |_| | |_| | | | |_| |  __/\__ \
   |   \___/ \__|_|_|_|\__|_|\___||___/
   |
  \*/

  //! Compute the square of a number.
  //! \param[in] x Number to be squared.
  //! \return Square of the number.
  Real power2(Real x) {return x*x;}

  //! Compute the cube of a number.
  //! \param[in] x Number to be cubed.
  //! \return Cube of the number.
  Real power3(Real x) {return x*x*x;}

  //! Convert an angle in degrees to radiants using the formula \f$ \text{rad} =
  //! \pi/180 \text{deg} \f$.
  //! \param[in] x Angle in degrees.
  //! \return Angle in radiants.
  Real degrees_to_radiants(Real x) {return DEG2RAD*x;}

  //! Convert an angle in radiants to degrees using the formula \f$ \text{deg} =
  //! 180/\pi\text{rad} \f$.
  //! \param[in] x Angle in radiants.
  //! \return Angle in degrees.
  Real radiants_to_degrees(Real x) {return RAD2DEG*x;}

  //! Add or remove multiple of \f$ 2\pi \f$ to an angle in order to clamp it in
  //! the range \f$ [0, 2\pi] \f$.
  //! \param[in] x Angle to be normalized.
  void angle_in_range(Real x) {
    x = std::fmod(x, PIMUL2);
    while (x < Real(0.0)) {x += PIMUL2;}
    while (x > PIMUL2) {x -= PIMUL2;}
  }

  //! Add or remove multiple of \f$ 2\pi \f$ to an angle in order to clamp it in
  //! the range \f$ [-\pi, \pi] \f$.
  //! \param[in] x Angle to be normalized.
  void angle_in_range_sym(Real x) {
    x = std::fmod(x, PIMUL2);
    while (x < -PI) {x += PIMUL2;}
    while (x > PI) {x -= PIMUL2;}
  }

} // namespace Astro

#endif // ASTRO_UTILITIES_HXX
