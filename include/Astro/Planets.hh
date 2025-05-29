/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Astro project is distributed under the GNU GPLv3.                     *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ASTRO_PLANETS_HH
#define ASTRO_PLANETS_HH

#include "Astro/Body.hh"
#include "Astro/Utilities.hh"

namespace Astro
{
  namespace Planets
  {

    static const Real KM3S2_TO_AU3DAY2{(KM_TO_AU*KM_TO_AU*KM_TO_AU)/(SEC_TO_DAY*SEC_TO_DAY)}; /**< Conversion factor from Km^3/s^2 to AU^3/day^2. */
    static const Real AU3DAY2_TO_KM3S2{1.0/KM3S2_TO_AU3DAY2}; /**< Conversion factor from AU^3/day^2 to Km^3/s^2. */

    /*\
     |   ____
     |  / ___| _   _ _ __
     |  \___ \| | | | '_ \
     |   ___) | |_| | | | |
     |  |____/ \__,_|_| |_|
     |
    \*/

    static Real const Sun_mass_KG{1.9885E30}; /** Mass of the Sun in Kg. */
    static Real const Sun_radius_M{695700000.0}; /**< Radius of the Sun in m. */
    static Real const Sun_radius_KM{Sun_radius_M/1.0E3}; /**< Radius of the Sun in Km. */
    static Real const Sun_radius_AU{Sun_radius_KM*KM_TO_AU}; /**< Radius of the Sun in AU. */
    static Real const Sun_mu_M3S2{1.32712440018E20}; /**< Gravitational constant of the Sun in m^3/s^2. */
    static Real const Sun_mu_KM3S2{Sun_mu_M3S2/1.0E9}; /**< Gravitational constant of the Sun in Km^3/s^2. */
    static Real const Sun_mu_AU3DAY2{Sun_mu_KM3S2*KM3S2_TO_AU3DAY2}; /**< Gravitational constant of the Sun in AU^3/day^2. */

    /*\
     |   _____           _   _
     |  | ____|__ _ _ __| |_| |__
     |  |  _| / _` | '__| __| '_ \
     |  | |__| (_| | |  | |_| | | |
     |  |_____\__,_|_|   \__|_| |_|
     |
    \*/

    static Real const Earth_mass_KG{5.97219E24}; /**< Mass of the Earth in Kg. */
    static Real const Earth_radius_M{6371000.0}; /**< Radius of the Earth in m. */
    static Real const Earth_radius_KM{Earth_radius_M/1.0E3}; /**< Radius of the Earth in Km. */
    static Real const Earth_radius_AU{Earth_radius_KM*KM_TO_AU}; /**< Radius of the Earth in AU. */
    static Real const Earth_mu_M3S2{3.986004418E14}; /**< Gravitational constant of the Earth in m^3/s^2. */
    static Real const Earth_mu_KM3S2{Earth_mu_M3S2/1.0E9}; /**< Gravitational constant of the Earth in Km^3/s^2. */
    static Real const Earth_mu_AU3DAY2{Earth_mu_KM3S2*KM3S2_TO_AU3DAY2}; /**< Gravitational constant of the Earth in AU^3/day^2. */

    /**
    * \brief Structure container for the Earth J2000 Keplerian orbital elements for orbit about the Sun.
    *
    * Structure container for the Earth J2000 Keplerian orbital elements for orbit about the Sun, which are:
    *   - the semi-major axis \f$ a = 1.00000011 \f$ (AU),
    *   - the eccentricity \f$ e = 0.01671022 \f$ (-),
    *   - the inclination \f$ i = 0.00005 \f$ (deg),
    *   - the longitude of the ascending node \f$ \Omega = 0.0 \f$ (deg),
    *   - the argument of periapsis \f$ \omega = 102.93768193 \f$ (deg).
    */
    struct KeplerianEarth : public OrbitalElements::Keplerian
    {
      KeplerianEarth()
      {
        this->a     = 1.00000011; // Semi-major axis (AU)
        this->e     = 0.01671022; // Eccentricity (-)
        this->i     = 0.00005*DEG_TO_RAD; // Inclination (rad)
        this->Omega = 0.0*DEG_TO_RAD; // Longitude of the ascending node (rad)
        this->omega = 102.93768193*DEG_TO_RAD; // Argument of periapsis (rad)
      }
    };

    /**
    * \brief Create a Earth object with J2000 Keplerian orbital elements.
    *
    * \return A new Earth object with J2000 Keplerian orbital elements.
    */
    inline Body Earth()
    {
      Body earth("Earth", Earth_mass_KG, Earth_radius_AU);
      earth.factor(Factor::POSIGRADE);
      earth.mu(Earth_mu_M3S2);
      earth.keplerian(KeplerianEarth());
      return earth;
    }

    /*\
     |   __  __
     |  |  \/  | ___   ___  _ __
     |  | |\/| |/ _ \ / _ \| '_ \
     |  | |  | | (_) | (_) | | | |
     |  |_|  |_|\___/ \___/|_| |_|
     |
    \*/

    static Real const Moon_mass_KG{7.34767309E22}; /**< Mass of the Moon in Kg. */
    static Real const Moon_radius_M{1737400.0}; /**< Radius of the Moon in m. */
    static Real const Moon_radius_KM{Moon_radius_M/1.0E3}; /**< Radius of the Moon in Km. */
    static Real const Moon_radius_AU{Moon_radius_KM*KM_TO_AU}; /**< Radius of the Moon in AU. */
    static Real const Moon_mu_M3S2{4.9048695E12}; /**< Gravitational constant of the Moon in m^3/s^2. */
    static Real const Moon_mu_KM3S2{Moon_mu_M3S2/1.0E9}; /**< Gravitational constant of the Moon in Km^3/s^2. */
    static Real const Moon_mu_AU3DAY2{Moon_mu_KM3S2*KM3S2_TO_AU3DAY2}; /**< Gravitational constant of the Moon in AU^3/day^2. */

    /**
    * \brief Structure container for the Moon Keplerian orbital elements for orbit about the Earth.
    *
    * Structure container for the Moon Keplerian orbital elements for orbit about the Earth, which are:
    *   - the semi-major axis \f$ a = 0.002569555 \f$ (AU),
    *   - the eccentricity \f$ e = 0.0549006 \f$ (-),
    *   - the inclination \f$ i = 5.145396 \f$ (deg),
    *   - the longitude of the ascending node \f$ \Omega = 125.1228 \f$ (deg),
    *   - the argument of periapsis \f$ \omega = 318.0634 \f$ (deg).
    */
    struct KeplerianMoon : public OrbitalElements::Keplerian
    {
      KeplerianMoon()
      {
        this->a     = 0.002569555; // Semi-major axis (AU)
        this->e     = 0.0549006; // Eccentricity (-)
        this->i     = 5.145396*DEG_TO_RAD; // Inclination (rad)
        this->Omega = 125.1228*DEG_TO_RAD; // Longitude of the ascending node (rad)
        this->omega = 318.0634*DEG_TO_RAD; // Argument of periapsis (rad)
      }
    };

  } // namespace Planets

} // namespace Astro

#endif // ASTRO_PLANETS_HH
