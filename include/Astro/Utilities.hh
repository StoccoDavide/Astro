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

#ifndef ASTRO_UTILITIES_HH
#define ASTRO_UTILITIES_HH

#include "Astro.hh"

namespace Astro {

  /*\
   |   ____
   |  |  _ \ _____      _____ _ __
   |  | |_) / _ \ \ /\ / / _ \ '__|
   |  |  __/ (_) \ V  V /  __/ |
   |  |_|   \___/ \_/\_/ \___|_|
   |
  \*/

  /**
  * Compute the square of a number.
  * \param[in] x Number to be squared.
  * \return Square of the number.
  */
  Real Power2(Real x) {return x*x;}

  /**
  * Compute the cube of a number.
  * \param[in] x Number to be cubed.
  * \return Cube of the number.
  */
  template<typename Real>
  Real Power3(Real x) {return x*x*x;}

  /**
  * Compute the fourth power of a number.
  * \param[in] x Number to be raised to the fourth power.
  * \return Number raised to the fourth power.
  */
  Real Power4(Real x) {return x*x*x*x;}

  /*\
   |      _                _
   |     / \   _ __   __ _| | ___
   |    / _ \ | '_ \ / _` | |/ _ \
   |   / ___ \| | | | (_| | |  __/
   |  /_/   \_\_| |_|\__, |_|\___|
   |                 |___/
  \*/

  static Real const PI         = Real(3.141592653589793238462643383279502884197); /**< Pi static constant value. */
  static Real const PIMUL2     = Real(6.283185307179586476925286766559005768394); /**< The value of \f$ 2\pi \f$. */
  static Real const PIDIV2     = Real(1.570796326794896619231321691639751442098); /**< The value of \f$ \pi/2 \f$. */
  static Real const DEG_TO_RAD = Real(0.017453292519943295769236907684886127134); /**< The value of \f$ \pi/180 \f$. */
  static Real const RAD_TO_DEG = Real(57.29577951308232087679815481410517033240); /**< The value of \f$ 180/\pi \f$. */

  /**
  * Convert an angle in degrees to radiants using the formula \f$ \text{rad} =
  * \pi/180 \text{deg} \f$.
  * \param[in] x Angle in degrees.
  * \return Angle in radiants.
  */
  Real Deg_To_Rad(Real x) {return DEG_TO_RAD*x;}

  /**
  * Convert an angle in radiants to degrees using the formula \f$ \text{deg} =
  * 180/\pi\text{rad} \f$.
  * \param[in] x Angle in radiants.
  * \return Angle in degrees.
  */
  Real Rad_To_Deg(Real x) {return RAD_TO_DEG*x;}

  /**
  * Add or remove multiple of \f$ 2\pi \f$ to an angle in order to clamp it in
  * the range \f$ [0, 2\pi] \f$.
  * \param[in] x Angle to be normalized.
  * \return Normalized angle.
  */
  Real AngleInRange(Real x)
  {
    x = std::fmod(x, PIMUL2);
    while (x < Real(0.0)) {x += PIMUL2;}
    while (x > PIMUL2) {x -= PIMUL2;}
    return x;
  }

  /**
  * Add or remove multiple of \f$ 2\pi \f$ to an angle in order to clamp it in
  * the range \f$ [-\pi, \pi] \f$.
  * \param[in] x Angle to be normalized.
  * \return Normalized angle.
  */
  Real AngleInRangeSym(Real x)
  {
    x = std::fmod(x, PIMUL2);
    while (x < -PI) {x += PIMUL2;}
    while (x > PI) {x -= PIMUL2;}
    return x;
  }

  /*\
   |   ____  _     _
   |  |  _ \(_)___| |_ __ _ _ __   ___ ___
   |  | | | | / __| __/ _` | '_ \ / __/ _ \
   |  | |_| | \__ \ || (_| | | | | (_|  __/
   |  |____/|_|___/\__\__,_|_| |_|\___\___|
   |
  \*/

  static Real const AU_TO_KM{1.49597870707e+08}; /**< One astronomical unit in kilometers /f$ 1 \text{AU} = 1.49597870707 \times 10^8 \text{km} \f$. */
  static Real const KM_TO_AU{1.0/AU_TO_KM};      /**< One kilometer in astronomical units /f$ 1 \text{km} = 1.49597870707 \times 10^{-8} \text{AU} \f$. */
  static Real const KM_TO_M{1000.0};             /**< One kilometer in meters /f$ 1 \text{km} = 1000 \text{m} \f$. */
  static Real const M_TO_KM{1.0/KM_TO_M};        /**< One meter in kilometers /f$ 1 \text{m} = 10^{-3} \text{km} \f$. */
  static Real const AU_TO_M{AU_TO_KM*KM_TO_M};   /**< One astronomical unit in meters /f$ 1 \text{AU} = 1.49597870707 \times 10^8 \times 1000 \text{m} \f$. */
  static Real const M_TO_AU{AU_TO_KM*KM_TO_M};   /**< One meter in astronomical units /f$ 1 \text{m} = 1.49597870707 \times 10^{-8} \text{AU} \f$. */


  // Time units
  static Real const DAY_TO_SEC{86400.0};        /**< One day in seconds /f$ 1 \text{day} = 86400 \text{s} \f$. */
  static Real const SEC_TO_DAY{1.0/DAY_TO_SEC}; /**< One second in days /f$ 1 \text{s} = 86400 \times 10^{-1} \text{day} \f$. */

  static Real const KG_M_SEC2_TO_KG_AU_DAY2{M_TO_AU/(SEC_TO_DAY*SEC_TO_DAY)}; // Kg * m/s^2 => Kg * UA / day^2

  static Real const gravity_kg_m_s2{9.80665};
  static Real const gravity_kg_AU_DAY2{gravity_kg_m_s2*KG_M_SEC2_TO_KG_AU_DAY2};


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Real KM_To_M(Real x) {return x * KM_TO_M;}
  Real M_To_KM(Real x) {return x * M_TO_KM;}

  Real KM_To_AU(Real x) {return x * KM_TO_AU;}
  Real AU_To_KM(Real x) {return x * AU_TO_KM;}

  Real M_To_AU(Real x) {return x * M_TO_AU;}
  Real AU_To_M(Real x) {return x * AU_TO_M;}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Real KM_S_To_M_S(Real x) {return x * KM_TO_M;}
  Real M_S_To_KM_S(Real x) {return x * M_TO_KM;}

  Real KM_S_To_AU_S(Real x) {return x * (KM_TO_AU/SEC_TO_DAY);}
  Real AU_S_To_KM_S(Real x) {return x * (AU_TO_KM*DAY_TO_SEC);}

  Real M_S_To_AU_S(Real x) {return x * (M_TO_AU/SEC_TO_DAY);}
  Real AU_S_To_M_S(Real x) {return x * (AU_TO_M*DAY_TO_SEC);}

  Real M_S_To_AU_DAY(Real x) {return x * (M_TO_AU/DAY_TO_SEC);}
  Real AU_DAY_To_M_S(Real x) {return x * (AU_TO_M*DAY_TO_SEC);}

  Real KM_S_To_AU_DAY(Real x) {return x * (KM_TO_AU/DAY_TO_SEC);}
  Real AU_DAY_To_KM_S(Real x) {return x * (AU_TO_KM*DAY_TO_SEC);}

  Real M_S_To_KM_DAY(Real x) {return x * (M_TO_KM/DAY_TO_SEC);}
  Real KM_DAY_To_M_S(Real x) {return x * (KM_TO_M*DAY_TO_SEC);}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Real AU_DAY2_To_KM_S2(Real x) {return x * (AU_TO_KM/(DAY_TO_SEC*DAY_TO_SEC));}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Real AU_DAY3_To_KM_S3(Real x) {return x * (AU_TO_KM/(DAY_TO_SEC*DAY_TO_SEC*DAY_TO_SEC));}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Real KM3_S2_To_KM3_DAY2(Real x) {return x * (1.0/(SEC_TO_DAY*SEC_TO_DAY));}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Real KM3_S2_To_AU3_DAY2(Real x) {return x * (1.0/(AU_TO_KM*AU_TO_KM*AU_TO_KM)/(SEC_TO_DAY*SEC_TO_DAY));}

} // namespace Astro

#endif // ASTRO_UTILITIES_HH
