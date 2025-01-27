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
  Real power2(Real x) {return x*x;}

  /**
  * Compute the cube of a number.
  * \param[in] x Number to be cubed.
  * \return Cube of the number.
  */
  Real power3(Real x) {return x*x*x;}

  /**
  * Compute the fourth power of a number.
  * \param[in] x Number to be raised to the fourth power.
  * \return Number raised to the fourth power.
  */
  Real power4(Real x) {return x*x*x*x;}

  /*\
   |      _                _
   |     / \   _ __   __ _| | ___
   |    / _ \ | '_ \ / _` | |/ _ \
   |   / ___ \| | | | (_| | |  __/
   |  /_/   \_\_| |_|\__, |_|\___|
   |                 |___/
  \*/

  /**
  * Convert an angle in degrees to radiants using the formula \f$ \text{rad} =
  * \pi/180 \text{deg} \f$.
  * \param[in] x Angle in degrees.
  * \return Angle in radiants.
  */
  Real deg_to_rad(Real x) {return DEG2RAD*x;}

  /**
  * Convert an angle in radiants to degrees using the formula \f$ \text{deg} =
  * 180/\pi\text{rad} \f$.
  * \param[in] x Angle in radiants.
  * \return Angle in degrees.
  */
  Real rad_to_deg(Real x) {return RAD2DEG*x;}

  /**
  * Add or remove multiple of \f$ 2\pi \f$ to an angle in order to clamp it in
  * the range \f$ [0, 2\pi] \f$.
  * \param[in] x Angle to be normalized.
  * \return Normalized angle.
  */
  Real angle_in_range(Real x)
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
  Real angle_in_range_sym(Real x)
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

  static Real const AU_TO_KM{1.49597870707e+08};   /**< One astronomical unit in kilometers /f$ 1 \text{AU} = 1.49597870707 \times 10^8 \text{km} \f$. */
  static Real const LY_TO_KM{9.4607304725808e+12}; /**< One light year in meters /f$ 1 \text{ly} = 9.4607304725808 \times 10^{12} \text{km} \f$. */
  static Real const PC_TO_AU{6.48e+05/PI};         /**< One parsec in astronomical units /f$ 1 \text{pc} = 648000/\pi \text{AU} \f$. */
  static Real const PC_TO_KM{PC_TO_AU*AU_TO_KM};   /**< One parsec in kilometers /f$ 1 \text{pc} = 648000/\pi \times 1.49597870707 \times 10^{8} \text{km} \f$. */
  static Real const AU_TO_PC{1.0/PC_TO_AU};        /**< One astronomical unit in parsecs /f$ 1 \text{AU} = \pi/648000 \text{pc} \f$. */
  static Real const LY_TO_PC{1.0/PC_TO_AU};        /**< One light year in parsecs /f$ 1 \text{ly} = \pi/648000 \text{pc} \f$. */
  static Real const KM_TO_PC{1.0/PC_TO_KM};        /**< One kilometer in parsecs /f$ 1 \text{km} = 1.0/648000/\pi \text{pc} \f$. */
  static Real const PC_TO_LY{1.0/LY_TO_PC};        /**< One parsec in light years /f$ 1 \text{pc} = 648000/\pi \times 1.49597870707 \times 10^{-8} \text{ly} \f$. */
  static Real const KM_TO_LY{1.0/LY_TO_KM};        /**< One kilometer in light years /f$ 1 \text{km} = 1.0570008340246 \times 10^{-13} \text{ly} \f$. */
  static Real const LY_TO_AU{LY_TO_KM/AU_TO_KM};   /**< One light year in astronomical units /f$ 1 \text{ly} = 9.4607304725808 \times 10^{12} \text{km} \f$. */
  static Real const AU_TO_LY{1.0/LY_TO_AU};        /**< One astronomical unit in light years /f$ 1 \text{AU} = 1.49597870707 \times 10^{-8} \text{ly} \f$. */
  static Real const KM_TO_AU{1.0/AU_TO_KM};        /**< One kilometer in astronomical units /f$ 1 \text{km} = 1.49597870707 \times 10^{-8} \text{AU} \f$. */

  static Real const KM_TO_M{1000.0};           /**< One kilometer in meters /f$ 1 \text{km} = 1000 \text{m} \f$. */
  static Real const M_TO_KM{1.0/KM_TO_M};      /**< One meter in kilometers /f$ 1 \text{m} = 10^{-3} \text{km} \f$. */
  static Real const AU_TO_M{AU_TO_KM*KM_TO_M}; /**< One astronomical unit in meters /f$ 1 \text{AU} = 1.49597870707 \times 10^8 \times 1000 \text{m} \f$. */
  static Real const LY_TO_M{LY_TO_KM*KM_TO_M}; /**< One light year in meters /f$ 1 \text{ly} = 9.4607304725808 \times 10^{12} \times 1000 \text{m} \f$. */
  static Real const PC_TO_M{LY_TO_KM*KM_TO_M}; /**< One parsec in meters /f$ 1 \text{pc} = 648000/\pi \times 1.49597870707 \times 10^{8} \times 1000 \text{m} \f$. */
  static Real const M_TO_PC{LY_TO_KM*KM_TO_M}; /**< One meter in parsecs /f$ 1 \text{m} = 1.0/648000/\pi \times 1.49597870707 \times 10^{8} \text{pc} \f$. */
  static Real const M_TO_LY{LY_TO_KM*KM_TO_M}; /**< One meter in light years /f$ 1 \text{m} = 1.0570008340246 \times 10^{-13} \text{ly} \f$. */
  static Real const M_TO_AU{LY_TO_KM*KM_TO_M}; /**< One meter in astronomical units /f$ 1 \text{m} = 1.49597870707 \times 10^{-8} \text{AU} \f$. */


  // Time units
  static Real const DAY_TO_SEC{86400.0};        /**< One day in seconds /f$ 1 \text{day} = 86400 \text{s} \f$. */
  static Real const SEC_TO_DAY{1.0/DAY_TO_SEC}; /**< One second in days /f$ 1 \text{s} = 86400 \times 10^{-1} \text{day} \f$. */

  static Real const KG_M_SEC2_TO_KG_AU_DAY2{M_TO_AU/(SEC_TO_DAY*SEC_TO_DAY)}; // Kg * m/s^2 => Kg * UA / day^2

  static Real const gravity_kg_m_s2{9.80665};
  static Real const gravity_kg_AU_DAY2{gravity_kg_m_s2*KG_M_SEC2_TO_KG_AU_DAY2};


  static Real const muSun_km3s2{1.32712440018E11}; // Km^3/s^2
  static Real const muSun_AU3DAY2{muSun_km3s2*(KM_TO_AU*KM_TO_AU*KM_TO_AU)/(SEC_TO_DAY*SEC_TO_DAY)}; // Km^3/s^2

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Real AU_to_KM(Real x) {return x * AU_TO_KM;}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Real AU_by_DAY_to_km_by_s(Real x) {return x * (AU_TO_KM/DAY_TO_SEC);}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Real AU_by_DAY2_to_km_by_s2(Real x) {return x * (AU_TO_KM/(DAY_TO_SEC*DAY_TO_SEC));}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Real AU_by_DAY3_to_km_by_s3(Real x) {return x * (AU_TO_KM/(DAY_TO_SEC*DAY_TO_SEC*DAY_TO_SEC));}

} // namespace Astro

#endif // ASTRO_UTILITIES_HXX
