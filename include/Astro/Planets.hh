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

    static const Real KM3_S2_TO_AU3_DAY2{(KM_TO_AU*KM_TO_AU*KM_TO_AU)/(SEC_TO_DAY*SEC_TO_DAY)}; /**< Conversion factor from Km^3/s^2 to AU^3/day^2. */
    static const Real AU3_DAY2_TO_KM3_S2{1.0/KM3_S2_TO_AU3_DAY2}; /**< Conversion factor from AU^3/day^2 to Km^3/s^2. */

    /*\
     |   ____
     |  / ___| _   _ _ __
     |  \___ \| | | | '_ \
     |   ___) | |_| | | | |
     |  |____/ \__,_|_| |_|
     |
    \*/

    static Real const Sun_mass_KG{1.9885E30}; /** Mass of the Sun in Kg. */
    static Real const Sun_radius_KM{695700.0}; /**< Radius of the Sun in Km. */
    static Real const Sun_radius_AU{KM_To_AU(Sun_radius_KM)}; /**< Radius of the Sun in AU. */
    static Real const Sun_mu_KM3_S2{1.32712440018E11}; /**< Gravitational constant of the Sun in Km^3/s^2. */
    static Real const Sun_mu_AU3_DAY2{KM3_S2_To_AU3_DAY2(Sun_mu_KM3_S2)}; /**< Gravitational constant of the Sun in AU^3/day^2. */

    /*\
     |   _____           _   _
     |  | ____|__ _ _ __| |_| |__
     |  |  _| / _` | '__| __| '_ \
     |  | |__| (_| | |  | |_| | | |
     |  |_____\__,_|_|   \__|_| |_|
     |
    \*/

    static Real const Earth_mass_KG{5.97219E24}; /**< Mass of the Earth in Kg. */
    static Real const Earth_radius_KM{6378.1370}; /**< Radius of the Earth in Km. */
    static Real const Earth_radius_AU{KM_To_AU(Earth_radius_KM)}; /**< Radius of the Earth in AU. */
    static Real const Earth_mu_KM3_S2{398600.4418}; /**< Gravitational constant of the Earth in Km^3/s^2. */
    static Real const Earth_mu_KM3_DAY2{KM3_S2_To_KM3_DAY2(Earth_mu_KM3_S2)}; /**< Gravitational constant of the Earth in Km^3/day^2. */
    static Real const Earth_mu_AU3_DAY2{KM3_S2_To_AU3_DAY2(Earth_mu_KM3_S2)}; /**< Gravitational constant of the Earth in AU^3/day^2. */

    /**
    * \brief Structure container for the Earth J2000 Keplerian orbital elements for orbit about the Sun.
    *
    * Structure container for the Earth J2000 Keplerian orbital elements for orbit about the Sun, which are:
    *   - the semi-major axis \f$ a = 1.00000011 \f$ (AU),
    *   - the eccentricity \f$ e = 0.01671022 \f$ (-),
    *   - the inclination \f$ i = 0.00005 \f$ (deg),
    *   - the longitude of the ascending node \f$ \Omega = -11.26064 \f$ (deg),
    *   - the argument of periapsis \f$ \omega = 102.93768193 \f$ (deg).
    */
    struct KeplerianEarth : public OrbitalElements::Keplerian
    {
      KeplerianEarth()
      {
        this->a     = 1.00000011; // Semi-major axis (AU)
        this->e     = 0.01671022; // Eccentricity (-)
        this->i     = Deg_To_Rad(0.00005); // Inclination (rad)
        this->Omega = Deg_To_Rad(-11.26064); // Longitude of the ascending node (rad)
        this->omega = Deg_To_Rad(102.93768193); // Argument of periapsis (rad)
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
      earth.set_orbit().set_factor(Factor::POSIGRADE);
      earth.set_orbit().set_mu(Sun_mu_AU3_DAY2);
      earth.set_orbit().set_keplerian(KeplerianEarth(), 0.0);
      earth.set_epoch(51544.0); // J2000 epoch (MJD 2451545.0 = JD 51544.5)
      earth.set_epoch_anomaly().set_M(Deg_To_Rad(100.46435), earth.orbit().keplerian(), Factor::POSIGRADE);
      return earth;
    }

    static constexpr Real EARTH_MAG_MOMENT{7.96e15}; /**< Magnetic dipole moment of the Earth in T*m^3. */

    /**
    * \brief Compute the Eath magnetic field at a given positionusing the IGRF model.
    *
    * \param[in] position The position vector in Km.
    * \return The magnetic field vector at the given position in nT.
    */
    inline Vector3 EarthMagneticFieldDipole(Vector3 const & position)
    {
      // Compute radial distance
      Real r{position.norm()};
      if (r < EPSILON_LOW) {return ZEROS_VEC3;} // avoid division by zero

      // Unit vector
      Vector3 r_hat{position/r};

      // Magnetic dipole moment aligned with z-axis (simplified)
      Vector3 m(0.0, 0.0, EARTH_MAG_MOMENT);

      Vector3 B{(3.0 * r_hat * m.dot(r_hat) - m) / Power3(r)};

      return B;
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
    static Real const Moon_radius_KM{1737.400}; /**< Radius of the Moon in Km. */
    static Real const Moon_radius_AU{Moon_radius_KM*KM_TO_AU}; /**< Radius of the Moon in AU. */
    static Real const Moon_mu_M3S2{4.9048695E12}; /**< Gravitational constant of the Moon in m^3/s^2. */
    static Real const Moon_mu_KM3_S2{Moon_mu_M3S2/1.0E9}; /**< Gravitational constant of the Moon in Km^3/s^2. */
    static Real const Moon_mu_AU3_DAY2{Moon_mu_KM3_S2*KM3_S2_TO_AU3_DAY2}; /**< Gravitational constant of the Moon in AU^3/day^2. */

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
        this->i     = Deg_To_Rad(5.145396); // Inclination (rad)
        this->Omega = Deg_To_Rad(125.1228); // Longitude of the ascending node (rad)
        this->omega = Deg_To_Rad(318.0634); // Argument of periapsis (rad)
      }
    };

    /**
    * \brief Create a Moon object with Keplerian orbital elements.
    * \return A new Moon object with Keplerian orbital elements.
    */
    inline Body Moon()
    {
      Body moon("Moon", Moon_mass_KG, Moon_radius_AU);
      moon.set_orbit().set_factor(Factor::POSIGRADE);
      moon.set_orbit().set_mu(Earth_mu_AU3_DAY2);
      moon.set_orbit().set_keplerian(KeplerianMoon(), 0.0);
      moon.set_epoch(51544.0); // J2000 epoch (JD 2451545.0 = MJD 51544.5)
      moon.set_epoch_anomaly().set_M(Deg_To_Rad(134.96340251), moon.orbit().keplerian(), Factor::POSIGRADE);
      return moon;
    }

  } // namespace Planets

} // namespace Astro

#endif // ASTRO_PLANETS_HH
