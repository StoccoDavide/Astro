/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Astro project is distributed under the GNU GPLv3.                     *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef ASTRO_ORBITAL_ELEMENTS_HXX
#define ASTRO_ORBITAL_ELEMENTS_HXX

namespace Astro
{

  /**
  * \namespace OrbitalElements
  * \brief The namespace for the orbital elements definition and conversion.
  *
  * \includedoc OrbitalElements.md
  */
  namespace OrbitalElements
  {
    /** Orbit posigrade/retrograte factor. */
    using Factor = enum class Factor : Integer {POSIGRADE = 1, RETROGRADE = -1};

    /**< Orbit type. */
    using Type = enum class Type : Integer {HYPERBOLIC = 0, ELLIPTIC = 1, PARABOLIC = 2};

    /*\
     |     ____           _            _
     |    / ___|__ _ _ __| |_ ___  ___(_) __ _ _ __
     |   | |   / _` | '__| __/ _ \/ __| |/ _` | '_ \
     |   | |__| (_| | |  | ||  __/\__ \ | (_| | | | |
     |    \____\__,_|_|   \__\___||___/_|\__,_|_| |_|
     |
    \*/

    /**
    * \brief Class container for the cartesian orbital parameters.
    *
    * Class container for the cartesian orbit parameters, which are:
    *   - the position vector \f$ \mathbf{r} \f$ (UA),
    *   - the velocity vector \f$ \mathbf{v} \f$ (UA/day).
    *
    * \note For more information on orbital elements, refer to *"Survey of Orbital Elements"*, by
    * G. R. Hintz, Journal of Guidance, Control, and Dynamics, Vol. 31, No. 3, May-June 2008.
    */
    class Cartesian
    {
      Vector3 m_r{NAN_VEC3}; /**< Position vector \f$ \mathbf{r} \f$ (UA). */
      Vector3 m_v{NAN_VEC3}; /**< Velocity vector \f$ \mathbf{v} \f$ (UA/day). */

    public:

      /**
      * Class constructor for Cartesian orbit parameters.
      */
      Cartesian() {}

      /**
      * Class constructor for Cartesian orbit parameters.
      * \param r Position vector \f$ \mathbf{r} \f$.
      * \param v Velocity vector \f$ \mathbf{v} \f$.
      */
      Cartesian(Vector3 const &r, Vector3 const &v) : m_r(r), m_v(v) {}

      /**
      * Class constructor for Cartesian orbit parameters.
      * \param r_x Position vector \f$ x \f$-axis component.
      * \param r_y Position vector \f$ y \f$-axis component.
      * \param r_z Position vector \f$ z \f$-axis component.
      * \param v_x Velocity vector \f$ x \f$-axis component.
      * \param v_y Velocity vector \f$ y \f$-axis component.
      * \param v_z Velocity vector \f$ z \f$-axis component.
      */
      Cartesian(Real r_x, Real r_y, Real r_z, Real v_x, Real v_y, Real v_z) :
        m_r(r_x, r_y, r_z), m_v(v_x, v_y, v_z) {}

      /**
      * Enable the default Cartesian orbit parameters copy constructor.
      */
      Cartesian(Cartesian const &) = default;

      /**
      * Enable the default Cartesian orbit parameters move constructor.
      */
      Cartesian(Cartesian &&) = default;

      /**
      * Enable the default Cartesian orbit parameters assignment operator.
      */
      Cartesian & operator=(const Cartesian &) = default;

      /**
      * Enable the default Cartesian orbit parameters move assignment operator.
      */
      Cartesian & operator=(Cartesian &&) = default;

      /**
      * Get the position vector \f$ \mathbf{r} \f$.
      * \return The position vector \f$ \mathbf{r} \f$.
      */
      Vector3 const & r() const {return this->m_r;}

      /**
      * Set the position vector \f$ \mathbf{r} \f$.
      * \param[in] t_r The position vector \f$ \mathbf{r} \f$.
      */
      void r(Vector3 const & t_r) {this->m_r = t_r;}

      /**
      * Get the velocity vector \f$ \mathbf{v} \f$.
      * \return The velocity vector \f$ \mathbf{v} \f$.
      */
      Vector3 const & v() const {return this->m_v;}

      /**
      * Get the velocity vector \f$ \mathbf{v} \f$.
      * \param[in] t_v The velocity vector \f$ \mathbf{v} \f$.
      */
      void v(Vector3 const & t_v) {this->m_v = t_v;}

      /**
      * Print the cartesian orbit parameters on a string.
      * \return The cartesian orbit parameters string.
      */
      std::string info() const {
        std::ostringstream os;
        os <<
          "r = " << this->m_r.transpose() << " (UA)" << std::endl <<
          "v = " << this->m_v.transpose() << " (UA/day)" << std::endl;
          return os.str();
      }

      /**
      * Print the cartesian orbit parameters on a stream.
      * \param[in,out] os Output stream.
      */
      void info(std::ostream & os) {os << this->info();}

      /**
      * Reset the cartesian orbit parameters to NaN.
      */
      void reset() {
        this->m_r = NAN_VEC3;
        this->m_v = NAN_VEC3;
      }

      /**
      * Check if the cartesian orbit parameters are valid, *i.e.*, finite.
      */
      bool sanity_check() const {
        #define CMD "Astro::OrbitalElements::Cartesian::sanity_check(...) "

        if (!(this->m_r.allFinite())) {
          ASTRO_ERROR(CMD "invalid position vector detected.");
          return false;
        }
        if (!(this->m_v.allFinite())) {
          ASTRO_ERROR(CMD "invalid velocity vector detected.");
          return false;
        }
        return true;

        #undef CMD
      }

      /**
      * Compute the orbital momentum vector \f$ \mathbf{h} = \mathbf{r} \times \mathbf{v} \f$ (UA^2/day).
      * \return The orbital momentum vector \f$ \mathbf{h} \f$.
      */
      Vector3 h() const {return this->m_r.cross(this->m_v);}

      /**
      * Compute the eccentricity \f$ e = \left\| \mathbf{e} \right\| \f$ (-), with
      * \f[ \mathbf{e} = \displaystyle\frac{\mathbf{v} \times \mathbf{h}}{\mu} - \displaystyle
      * \frac{\mathbf{r}}{r} \f$
      * \param[in] mu The gravitational parameter \f$ \mu \f$.
      * \return The eccentricity \f$ e \f$.
      */
      Real e(Real mu) const
      {
        return ((this->m_v.cross(this->h()) / mu) - (this->m_r / this->m_r.norm())).norm();
      }

    }; // class Cartesian

    /*\
     |    _  __          _           _
     |   | |/ /___ _ __ | | ___ _ __(_) __ _ _ __
     |   | ' // _ \ '_ \| |/ _ \ '__| |/ _` | '_ \
     |   | . \  __/ |_) | |  __/ |  | | (_| | | | |
     |   |_|\_\___| .__/|_|\___|_|  |_|\__,_|_| |_|
     |            |_|
    \*/

    /**
    * \brief Class container for the (modified) Keplerian (or classical) orbital parameters.
    *
    * Class container for the (modified) Keplerian (or classical) orbit parameters, which are:
    *   - the semi-major axis \f$ a \f$ (UA),
    *   - the eccentricity \f$ e \in [0, 1] \f$ (-),
    *   - the inclination \f$ i \f$ (rad),
    *   - the right ascension of the ascending node \f$ \Omega \f$ (rad),
    *   - the argument of periapsis \f$ \omega \f$ (rad).
    *   - the anomaly, which can be:
    * Singular at \f$ e = 0 \f$, or \f$ 1 \f$, \f$ i = 0 \f$, or \f$ \pi \f$.
    *
    * \note For more information on orbital elements, refer to *"Survey of Orbital Elements"*, by
    * G. R. Hintz, Journal of Guidance, Control, and Dynamics, Vol. 31, No. 3, May-June 2008.
    */
    class Keplerian
    {
      Real m_a{QUIET_NAN};     /**< Semi-major axis \f$ a \f$ (UA). */
      Real m_e{QUIET_NAN};     /**< Eccentricity \f$ e \in [0, 1] \f$ (-). */
      Real m_i{QUIET_NAN};     /**< Inclination \f$ i \f$ (rad). */
      Real m_Omega{QUIET_NAN}; /**< Right ascension of the ascending node \f$ \Omega \f$ (rad). */
      Real m_omega{QUIET_NAN}; /**< Argument of periapsis \f$ \omega \f$ (rad). */

    public:

      /**
      * Class constructor for the (modified) Keplerian orbit parameters.
      */
      Keplerian() {}

      /**
      * Class constructor for the (modified) Keplerian orbit parameters.
      * \param[in] t_a The semi-major axis \f$ a \f$.
      * \param[in] t_e The eccentricity \f$ e \f$.
      * \param[in] t_i The inclination \f$ i \f$.
      * \param[in] t_Omega The longitude of the ascending node \f$ \Omega \f$.
      * \param[in] t_omega The argument of periapsis \f$ \omega \f$.
      */
      Keplerian(Real t_a, Real t_e, Real t_i, Real t_Omega, Real t_omega)
        : m_a(t_a), m_e(t_e), m_i(t_i), m_Omega(t_Omega), m_omega(t_omega)
      {}

      /**
      * Enable the default Keplerian orbit parameters copy constructor.
      */
      Keplerian(Keplerian const &) = default;

      /**
      * Enable the default Keplerian orbit parameters move constructor.
      */
      Keplerian(Keplerian &&) = default;

      /**
      * Enable the default Keplerian orbit parameters assignment operator.
      */
      Keplerian & operator=(const Keplerian &) = default;

      /**
      * Enable the default Keplerian orbit parameters move assignment operator.
      */
      Keplerian & operator=(Keplerian &&) = default;

      /**
      * Get the semi-major axis \f$ a \f$.
      * \return The semi-major axis \f$ a \f$.
      */
      Real a() const {return this->m_a;}

      /**
      * Set the semi-major axis \f$ a \f$.
      * \param[in] t_a The semi-major axis \f$ a \f$.
      */
      void a(Real t_a) {this->m_a = t_a;}

      /**
      * Get the eccentricity \f$ e \f$.
      * \return The eccentricity \f$ e \f$.
      */
      Real e() const {return this->m_e;}

      /**
      * Set the eccentricity \f$ e \f$.
      * \param[in] t_e The eccentricity \f$ e \f$.
      */
      void e(Real t_e) {this->m_e = t_e;}

      /**
      * Get the inclination \f$ i \f$.
      * \return The inclination \f$ i \f$.
      */
      Real i() const {return this->m_i;}

      /**
      * Set the inclination \f$ i \f$.
      * \param[in] t_i The inclination \f$ i \f$.
      */
      void i(Real t_i) {this->m_i = t_i;}

      /**
      * Get the longitude of the ascending node \f$ \Omega \f$.
      * \return The longitude of the ascending node \f$ \Omega \f$.
      */
      Real Omega() const {return this->m_Omega;}

      /**
      * Set the longitude of the ascending node \f$ \Omega \f$.
      * \param[in] t_Omega The longitude of the ascending node \f$ \Omega \f$.
      */
      void Omega(Real t_Omega) {this->m_Omega = t_Omega;}

      /**
      * Get the argument of periapsis \f$ \Omega \f$
      * \return The argument of periapsis \f$ \Omega \f$
      */
      Real omega() const {return this->m_omega;}

      /**
      * Set the argument of periapsis \f$ \omega \f$.
      * \param[in] t_omega The argument of periapsis \f$ \omega \f$.
      */
      void omega(Real t_omega) {this->m_omega = t_omega;}

      /**
      * Print the Keplerian orbit parameters on a string.
      * \return The Keplerian orbit parameters string.
      */
      std::string info() const {
        std::ostringstream os;
        os <<
          "a : semi-major axis   = " << this->m_a << std::endl <<
          "e : eccentricity      = " << this->m_e << std::endl <<
          "i : inclination       = " << this->m_i << " (rad) = " << rad_to_deg(this->m_i) << " (deg)" << std::endl <<
          "Ω : right ascension … = " << this->m_Omega << " (rad) = " << rad_to_deg(this->m_Omega) << " (deg)" << std::endl <<
          "ω : arg. of periapsis = " << this->m_omega << " (rad) = " << rad_to_deg(this->m_omega) << " (deg)" << std::endl;
          return os.str();
      }

      /**
      * Print the keplerian orbit parameters on a stream.
      * \param[in,out] os Output stream.
      */
      void info(std::ostream & os) {os << this->info();}

      /**
      * Reset the Keplerian orbit parameters to NaN.
      */
      void reset() {
        this->m_a = QUIET_NAN;
        this->m_e = QUIET_NAN;
        this->m_i = QUIET_NAN;
        this->m_Omega = QUIET_NAN;
        this->m_omega = QUIET_NAN;
      }

      /**
      * Check if the keplerian orbit parameters are valid, *i.e.*, finite, and with \f$ e > 0 \f$,
      * \f$ a > 0 \f$.
      * \param[in] tol_i Tolerance \f$ \varepsilon_i \f$ for the singularity check on the inclination.
      * \param[in] tol_e Tolerance \f$ \varepsilon_e \f$ for the singularity check on the eccentricity.
      * \return True if the keplerian orbit parameters are valid, false otherwise.
      */
      bool sanity_check() const
      {
        #define CMD "Astro::OrbitalElements::Keplerian::sanity_check(...) "

        if (!std::isfinite(this->m_a) && this->m_a > 0.0) {
          ASTRO_ERROR(CMD "invalid semi-major axis, a = " << this->m_a << ".");
          return false;
        }
        if (!std::isfinite(this->m_e) && this->m_e > 0.0 && this->m_e < 1.0) {
          ASTRO_ERROR(CMD "invalid eccentricity, e = " << this->m_e << ".");
          return false;
        }
        if (!std::isfinite(this->m_i)) {
          ASTRO_ERROR(CMD "invalid inclination, i = " << this->m_i << ".");
          return false;
        }
        if (!std::isfinite(this->m_Omega)) {
          ASTRO_ERROR(CMD "invalid right ascension of the ascending node, Ω = " << this->m_Omega << ".");
          return false;
        }
        if (!std::isfinite(this->m_omega)) {
          ASTRO_ERROR(CMD "invalid argument of periapsis, ω = " << this->m_omega << ".");
          return false;
        }
        return true;

        #undef CMD
      }

      /**
      * Check if the keplerian orbit is singular, *i.e.*, *\f$ i < \pi-\varepsilon_i \f$, \f$ i >
      * \pi+\varepsilon_i\f$, or \f$ e \geq 1 - \varepsilon_e\f$.
      * \param[in] tol_i Tolerance \f$ \varepsilon_i \f$ for the singularity check on the inclination.
      * \param[in] tol_e Tolerance \f$ \varepsilon_e \f$ for the singularity check on the eccentricity.
      * \return True if the keplerian orbit is singular, false otherwise.
      */
      bool is_singular(Real tol_i = EPSILON_LOW, Real tol_e = EPSILON_LOW) const
      {
        #define CMD "Astro::OrbitalElements::Keplerian::is_singular(...) "

        Real i{angle_in_range(this->m_i)};
        if (i < PI + tol_i && i > PI - tol_i) {
          ASTRO_WARNING(CMD "singular inclination detected.");
          return true;
        }
        if (this->m_e > 1.0 - tol_e) {
          ASTRO_WARNING(CMD "singular eccentricity detected.");
          return true;
        }
        return false;

        #undef CMD
      }

      /**
      * Check if the keplerian orbit is singular, *i.e.*, *\f$ \pi-\varepsilon_i < i < \pi+\varepsilon_i\f$,
      * and \f$ e < 1 - \varepsilon_e\f$.
      * \param[in] tol_i Tolerance \f$ \varepsilon_i \f$ for the singularity check on the inclination.
      * \param[in] tol_e Tolerance \f$ \varepsilon_e \f$ for the singularity check on the eccentricity.
      * \return True if the keplerian orbit is nonsingular, false otherwise.
      */
      bool is_nonsingular(Real tol_i = EPSILON_LOW, Real tol_e = EPSILON_LOW) const
      {
        return !this->is_singular(tol_i, tol_e);
      }

      /**
      * Compute the argument of latitude \f$ u = \omega + \Omega \f$.
      * \return The argument of latitude \f$ u \f$.
      */
      Real u() const {return this->m_omega + this->m_Omega;}

      /**
      * Compute the semi-latus rectum \f$ p = a(1 - e^2) \f$ (UA).
      * \return The semi-latus rectum \f$ p \f$.
      */
      Real p() const {return this->m_a * (1.0 - this->m_e * this->m_e);}

    }; // class Keplerian

    /*\
     |   _____            _                  _   _       _
     |  | ____|__ _ _   _(_)_ __   ___   ___| |_(_) __ _| |
     |  |  _| / _` | | | | | '_ \ / _ \ / __| __| |/ _` | |
     |  | |__| (_| | |_| | | | | | (_) | (__| |_| | (_| | |
     |  |_____\__, |\__,_|_|_| |_|\___/ \___|\__|_|\__,_|_|
     |           |_|
    \*/

    /**
    * \brief Class container for the (modified) Equinoctical orbital parameters.
    *
    * The (modified) Equinoctical orbital parameters are defined as follows:
    *   - the semi-latus rectum \f$ p \f$ (UA).
    *   - the \f$ x \f$-axis component of the eccentricity vector in the orbital frame \f$ f \f$ (-).
    *   - the \f$ y \f$-axis component of the eccentricity vector in the orbital frame \f$ g \f$ (-).
    *   - the \f$ x \f$-axis component of the node vector in the orbital frame \f$ h \f$ (-).
    *   - the \f$ y \f$-axis component of the node vector in the orbital frame \f$ k \f$ (-).
    *   - the posigrade (+1)/retrograde (-1) factor \f$ I = +1 \f$ (posigrade) or \f$ I = -1 \f$ (retrograde) (-).
    * Singular at \f$ i = \pi \f$.
    *
    * \note For more information on orbital elements, refer to *"Survey of Orbital Elements"*, by
    * G. R. Hintz, Journal of Guidance, Control, and Dynamics, Vol. 31, No. 3, May-June 2008.
    */
    class Equinoctical
    {
      Real m_p{QUIET_NAN}; /**< Semi-latus rectum \f$ p \f$ (UA). */
      Real m_f{QUIET_NAN}; /**< \f$ X \f$-axis component of the eccentricity vector in the orbital frame (-). */
      Real m_g{QUIET_NAN}; /**< \f$ Y \f$-axis component of the eccentricity vector in the orbital frame (-). */
      Real m_h{QUIET_NAN}; /**< \f$ X \f$-axis component of the node vector in the orbital frame (-). */
      Real m_k{QUIET_NAN}; /**< \f$ Y \f$-axis component of the node vector in the orbital frame (-). */

    public:

      /**
      * Class constructor for the (modified) Equinoctical orbit parameters.
      */
      Equinoctical() {}

      /**
      * Class constructor for the (modified) Equinoctical orbit parameters.
      * \param[in] t_p The semi-latus rectum \f$ p \f$.
      * \param[in] t_f The \f$ x \f$-axis component of the eccentricity vector in the orbital frame.
      * \param[in] t_g The \f$ y \f$-axis component of the eccentricity vector in the orbital frame.
      * \param[in] t_h The \f$ x \f$-axis component of the node vector in the orbital frame.
      * \param[in] t_k The \f$ y \f$-axis component of the node vector in the orbital frame.
      */
      Equinoctical(Real t_p, Real t_f, Real t_g, Real t_h, Real t_k)
        : m_p(t_p), m_f(t_f), m_g(t_g), m_h(t_h), m_k(t_k) {}

      /**
      * Enable the default equinoctial orbit parameters copy constructor.
      */
      Equinoctical(Equinoctical const &) = default;

      /**
      * Enable the default equinoctial orbit parameters move constructor.
      */
      Equinoctical(Equinoctical &&) = default;

      /**
      * Enable the default equinoctial orbit parameters assignment operator.
      */
      Equinoctical & operator=(const Equinoctical &) = default;

      /**
      * Enable the default equinoctial orbit parameters move assignment operator.
      */
      Equinoctical & operator=(Equinoctical &&) = default;

      /**
      * Get the semi-latus rectum \f$ p \f$.
      * \return The semi-latus rectum \f$ p \f$.
      */
      Real p() const {return this->m_p;}

      /**
      * Set the semi-latus rectum \f$ p \f$.
      * \param[in] t_p The semi-latus rectum \f$ p \f$.
      */
      void p(Real t_p) {this->m_p = t_p;}

      /**
      * Get the \f$ x \f$-axis component of the eccentricity vector in the orbital frame \f$ f \f$.
      * \return The \f$ x \f$-axis component of the eccentricity vector in the orbital frame \f$ f \f$.
      */
      Real f() const {return this->m_f;}

      /**
      * Set the \f$ x \f$-axis component of the eccentricity vector in the orbital frame \f$ f \f$.
      * \param[in] t_f The \f$ x \f$-axis component of the eccentricity vector in the orbital frame \f$ f \f$.
      */
      void f(Real t_f) {this->m_f = t_f;}

      /**
      * Get the \f$ y \f$-axis component of the eccentricity vector in the orbital frame \f$ g \f$.
      * \return The \f$ y \f$-axis component of the eccentricity vector in the orbital frame \f$ g \f$.
      */
      Real g() const {return this->m_g;}

      /**
      * Set the \f$ y \f$-axis component of the eccentricity vector in the orbital frame \f$ g \f$.
      * \param[in] t_g The \f$ y \f$-axis component of the eccentricity vector in the orbital frame \f$ g \f$.
      */
      void g(Real t_g) {this->m_g = t_g;}

      /**
      * Get the \f$ x \f$-axis component of the node vector in the orbital frame \f$ h \f$.
      * \return The \f$ x \f$-axis component of the node vector in the orbital frame \f$ h \f$.
      */
      Real h() const {return this->m_h;}

      /**
      * Set the \f$ x \f$-axis component of the node vector in the orbital frame \f$ h \f$.
      * \param[in] t_h The \f$ x \f$-axis component of the node vector in the orbital frame \f$ h \f$.
      */
      void h(Real t_h) {this->m_h = t_h;}

      /**
      * Get the \f$ y \f$-axis component of the node vector in the orbital frame \f$ k \f$.
      * \return The \f$ y \f$-axis component of the node vector in the orbital frame \f$ k \f$.
      */
      Real k() const {return this->m_k;}

      /**
      * Set the \f$ y \f$-axis component of the node vector in the orbital frame \f$ k \f$.
      * \param[in] t_k The \f$ y \f$-axis component of the node vector in the orbital frame \f$ k \f$.
      */
      void k(Real t_k) {this->m_k = t_k;}

      /**
      * Print the equinoctial orbit parameters on a string.
      * \return Equinoctical orbit parameters string.
      */
      std::string info() const {
        std::ostringstream os;
        os <<
          "p : semi-latus rectum  = " << this->m_p << " (UA)" << std::endl <<
          "f : x-axis ecc. vector = " << this->m_f << " (-)" << std::endl <<
          "g : y-axis ecc. vector = " << this->m_g << " (-)" << std::endl <<
          "h : x-axis node vector = " << this->m_h << " (-)" << std::endl <<
          "k : y-axis node vector = " << this->m_k << " (-)" << std::endl;
          return os.str();
      }

      /**
      * Print the equinoctial orbit parameters on a stream.
      * \param[in,out] os Output stream.
      */
      void info(std::ostream & os) {os << this->info();}

      /**
      * Reset the equinoctial orbit parameters to NaN.
      */
      void reset() {
        this->m_p = QUIET_NAN;
        this->m_f = QUIET_NAN;
        this->m_g = QUIET_NAN;
        this->m_h = QUIET_NAN;
        this->m_k = QUIET_NAN;
      }

      /**
      * Check if the equinoctial orbit parameters are valid, *i.e.*, finite, and with \f$ p > 0 \f$.
      * \return True if the equinoctial orbit parameters are valid, false otherwise.
      */
      bool sanity_check() const
      {
        #define CMD "Astro::OrbitalElements::Keplerian::sanity_check(...) "

        if (!std::isfinite(this->m_p) && this->m_p > 0.0) {
          ASTRO_ERROR(CMD "invalid semi-latus rectum, p = " << this->m_p << ".");
          return false;
        }
        if (!std::isfinite(this->m_f)) {
          ASTRO_ERROR(CMD "invalid x-axis component of the eccentricity vector, f = " << this->m_f << ".");
          return false;
        }
        if (!std::isfinite(this->m_g)) {
          ASTRO_ERROR(CMD "invalid y-axis component of the eccentricity vector, g = " << this->m_g << ".");
          return false;
        }
        if (!std::isfinite(this->m_h)) {
          ASTRO_ERROR(CMD "invalid x-axis component of the node vector, h = " << this->m_h << ".");
          return false;
        }
        if (!std::isfinite(this->m_k)) {
          ASTRO_ERROR(CMD "invalid y-axis component of the node vector, k = " << this->m_k << ".");
          return false;
        }

        // Check if it is singular
        ASTRO_ASSERT_WARNING(this->is_nonsingular(), CMD "singular orbit detected.");

        return true;

        #undef CMD
      }

      /**
      * Check if the equinoctial orbit is singular, *i.e.*, *\f$ i < \pi-\varepsilon_i \f$, or
      * \f$ i > \pi+\varepsilon_i\f$.
      * \param[in] tol_i Tolerance \f$ \varepsilon_i \f$ for the singularity check on the inclination.
      * \return True if the equinoctial orbit is singular, false otherwise.
      */
      bool is_singular(Real tol_i = EPSILON_LOW) const
      {
        #define CMD "Astro::OrbitalElements::Keplerian::is_singular(...) "

        Real i{angle_in_range(this->i())};
        if (i < PI + tol_i && i > PI - tol_i) {
          ASTRO_WARNING(CMD "singular inclination detected.");
          return true;
        }
        return false;

        #undef CMD
      }

      /**
      * Check if the equinoctial orbit is singular, *i.e.*, *\f$ \pi-\varepsilon_i < i < \pi+\varepsilon_i\f$.
      * \param[in] tol_i Tolerance \f$ \varepsilon_i \f$ for the singularity check on the inclination.
      * \return True if the equinoctial orbit is nonsingular, false otherwise.
      */
      bool is_nonsingular(Real tol_i = EPSILON_LOW) const
      {
        return !this->is_singular(tol_i);
      }

      /**
      * Compute the argument of latitude \f$ u = \omega + \Omega = \arctan\left(h\sin(L) - k\cos(L),
      * h\cos(L) + k\sin(L)\right) \f$.
      * \return The argument of latitude \f$ u \f$.
      */
      Real u() const
      {
        Real s_L{std::sin(this->L())}, c_L{std::cos(this->L())};
        return std::atan2(this->m_h*s_L - this->m_k*c_L, this->m_h*c_L + this->m_k*s_L);
      }

      /**
      * Compute the semi-major axis \f$ a = \displaystyle\frac{p}{1 - f^2 - g^2} \f$ (UA).
      * \return The semi-major axis \f$ a \f$.
      */
      Real a() const {return this->m_p / (1.0 - this->m_f*this->m_f - this->m_g*this->m_g);}

      /**
      * Compute the eccentricity \f$ \mathbf{e} = [f, g, 0]^\top \f$ (-).
      * \return The eccentricity \f$ \mathbf{e} \f$.
      */
      Real e() const {return std::sqrt(this->m_f*this->m_f + this->m_g*this->m_g);}

      /**
      * Compute the inclination \f$ i = \arctan\left(2\sqrt{h^2 + k^2}, 1 - h^2 - k^2\right) \f$ (rad).
      * \return The inclination \f$ i \f$.
      */
      Real i() const
      {
        Real h2{this->m_h * this->m_h}, k2{this->m_k * this->m_k};
        return std::atan2(2.0*std::sqrt(h2 + k2), 1.0 - h2 - k2);
      }

      /**
      * Compute the argument of periapsis \f$ \omega = \arctan\left(gh - fk, fh + gk\right) \f$ (rad).
      * \return The argument of periapsis \f$ \omega \f$.
      */
      Real omega() const
      {
        return std::atan2(
          this->m_g*this->m_h - this->m_f*this->m_k,
          this->m_f*this->m_h + this->m_g*this->m_k
        );
      }

      /**
      * Compute the right ascension of the ascending node \f$ \Omega = \arctan\left(k, h\right) \f$ (rad).
      * \return The right ascension of the ascending node \f$ \Omega \f$.
      */
      Real Omega() const {return std::atan2(this->m_k, this->m_h);}

    }; // class Equinoctical


    /*\
     |    ___              _                  _             _
     |   / _ \ _   _  __ _| |_ ___ _ __ _ __ (_) ___  _ __ (_) ___
     |  | | | | | | |/ _` | __/ _ \ '__| '_ \| |/ _ \| '_ \| |/ __|
     |  | |_| | |_| | (_| | ||  __/ |  | | | | | (_) | | | | | (__
     |   \__\_\\__,_|\__,_|\__\___|_|  |_| |_|_|\___/|_| |_|_|\___|
     |
    \*/

    /**
    * \brief Class container for the quaternionic orbital parameters.
    *
    * Class container for the quaternionic orbit parameters, which is made of four quaternionic
    * parameters: \f$ [q¹, q², q³, q⁴] \f$.
    *
    * \note For more information on the quaternionic orbital elements, refer to *"Alternative Set of
    * Nonsingular Quaternionic Orbital Elements"*, by J. Roa and J. Kasdin, Journal of Guidance,
    * Control, and Dynamics, Vol. 40, No. 11, November 2017.
    */
    class Quaternionic
    {
      Vector4 m_q{NAN_VEC4}; /**< Quaternionic orbit parameters. */

    public:

      /**
      * Class constructor for the quaternionic orbit parameters.
      */
      Quaternionic() {}

      /**
      * Class constructor for the quaternionic orbit parameters.
      * \param[in] t_q The quaternionic orbit parameters.
      */
      Quaternionic(Vector4 const &t_q) : m_q(t_q) {}

      /**
      * Class constructor for the quaternionic orbit parameters.
      * \param[in] t_q_1 The first quaternionic orbit parameter.
      * \param[in] t_q_2 The second quaternionic orbit parameter.
      * \param[in] t_q_3 The third quaternionic orbit parameter.
      * \param[in] t_q_4 The fourth quaternionic orbit parameter.
      */
      Quaternionic(Real t_q_1, Real t_q_2, Real t_q_3, Real t_q_4)
        : m_q(t_q_1, t_q_2, t_q_3, t_q_4) {}

      /**
      * Enable the default quaternionic orbit parameters copy constructor.
      */
      Quaternionic(Quaternionic const &) = default;

      /**
      * Enable the default quaternionic orbit parameters move constructor.
      */
      Quaternionic(Quaternionic &&) = default;

      /**
      * Enable the default quaternionic orbit parameters assignment operator.
      */
      Quaternionic & operator=(const Quaternionic &) = default;

      /**
      * Enable the default quaternionic orbit parameters move assignment operator.
      */
      Quaternionic & operator=(Quaternionic &&) = default;

      /**
      * Get the quaternionic orbit parameters.
      * \return The quaternionic orbit parameters.
      */
      Vector4 const & q() const {return this->m_q;}

      /**
      * Set the quaternionic orbit parameters.
      * \param[in] t_q The quaternionic orbit parameters.
      */
      void q(Vector4 const & t_q) {this->m_q = t_q;}

      /**
      * Print the quaternionic orbit parameters on a string.
      * \return The quaternionic orbit parameters string.
      */
      std::string info() const {
        std::ostringstream os;
        os <<
          "q¹ : 1st parameter = " << this->m_q(0) << " (-)" << std::endl <<
          "q² : 2nd parameter = " << this->m_q(1) << " (-)" << std::endl <<
          "q³ : 3rd parameter = " << this->m_q(2) << " (-)" << std::endl <<
          "q⁴ : 4th parameter = " << this->m_q(3) << " (-)" << std::endl;
          return os.str();
      }

      /**
      * Print the quaternionic orbit parameters on a stream.
      * \param[in,out] os Output stream.
      */
      void info(std::ostream & os) {os << this->info();}

      /**
      * Reset the quaternionic orbit parameters to NaN.
      */
      void reset() {this->m_q = NAN_VEC4;}

      /**
      * Check if the quaternionic orbit parameters are valid, *i.e.*, finite.
      * \return True if the quaternionic orbit parameters are valid, false otherwise.
      */
      bool sanity_check() const
      {
        #define CMD "Astro::OrbitalElements::Quaternionic::sanity_check(...) "

        if (!(this->m_q.allFinite())) {
          ASTRO_ERROR(CMD "invalid quaternionic orbit parameters.");
          return false;
        }
        return true;

        #undef CMD
      }

    }; // class Quaternionic

    /*\
     |    ____                          _
     |   / ___|___  _ ____   _____ _ __| |_
     |  | |   / _ \| '_ \ \ / / _ \ '__| __|
     |  | |__| (_) | | | \ V /  __/ |  | |_
     |   \____\___/|_| |_|\_/ \___|_|   \__|
     |
    \*/

    /**
    * Namespace to collect the orbital elements and anomalies conversions\.
    */
    namespace Convert
    {

      /**
      * Compute the mean anomaly \f$ M \f$ from the eccentric anomaly \f$ E \f$ as
      * \f[ M = \nu - 2\arctan\left(\frac{\sin(\nu)}{\beta + \cos(\nu)}\right) -
      * \displaystyle\frac{e\sqrt{1 - e^2}\sin(\nu)}{1 + e\cos(\nu)} \text{,} \f]
      * with \f$ \beta = \displaystyle\frac{1 + \sqrt{1 - e^2}}{e} \f$. That is equivalent to
      *
      * \note Refer to *"A note on the relations between true and eccentric anomalies in the two-body
      * problem"* by R. Broucke and P. Cefola, Celestial Mechanics, Vol. 7, pp. 300-389, 1973.
      * \param[in] nu The true anomaly \f$ \nu \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \return The mean anomaly \f$ M \f$.
      */
      Real nu_to_M(Real nu, Keplerian const & kepl)
      {
        Real e{kepl.e()};
        Real beta{(1.0 + std::sqrt(1.0 - e*e)) / e};
        return nu - 2.0 * std::atan(std::sin(nu) / (beta + std::cos(nu))) -
          (e * std::sqrt(1.0 - e*e) * std::sin(nu)) / (1.0 + e * std::cos(nu));
      }

      /**
      * Compute the true anomaly \f$ \nu \f$ from the eccentric anomaly \f$ E \f$ as
      * \f[ \tan\left(\frac{\nu - E}{2}\right) = \displaystyle\frac{\sin(\nu)}{\beta + \cos(\nu)} \text{,} \f]
      * with \f$ \beta = \displaystyle\frac{1 + \sqrt{1 - e^2}}{e} \f$. That is equivalent to
      * \f[ E = -\nu - 2\arctan\left(\frac{\sin(\nu)}{\beta + \cos(\nu)}\right) \text{.} \f]
      *
      \note Refer to *"A note on the relations between true and eccentric anomalies in the two-body
      * problem"* by R. Broucke and P. Cefola, Celestial Mechanics, Vol. 7, pp. 300-389, 1973.
      * \param[in] nu The true anomaly \f$ \nu \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \return The eccentric anomaly \f$ E \f$.
      */
      Real nu_to_E(Real nu, Keplerian const & kepl)
      {
        Real e{kepl.e()};
        Real beta{(1.0 + std::sqrt(1.0 - e*e)) / e};
        return -nu - 2.0 * std::atan(std::sin(nu) / (beta + std::cos(nu)));
      }

      /**
      * Compute the true longitude \f$ L \f$ from the true anomaly \f$ \nu \f$ as
      * \f[ L = \nu + \omega + I\Omega \text{.} \f]
      * \param[in] nu The true anomaly \f$ \nu \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \param[in] I The posigrade (+1)/retrograde (-1) factor.
      * \return The true longitude \f$ L \f$.
      */
      Real nu_to_L(Real nu, Keplerian const & kepl, Factor I)
      {
        return nu + kepl.omega() + kepl.Omega()*static_cast<Real>(I);
      }

      /**
      * Compute the eccentric anomaly \f$ E \f$ from the mean anomaly \f$ M \f$ through a the solution
      * of the nonlinear equation \f$ E = M + e\sin(E) \f$ for \f$ E \f$. The solution is found through
      * a basic Newton method.
      * \param[in] M The mean anomaly \f$ M \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \return The eccentric anomaly \f$ E \f$.
      */
      Real M_to_E(Real M, Keplerian const & kepl)
      {
        #define CMD "Astro::OrbitalElements::Anomaly::M_to_E(...): "

        // Clamp the mean anomaly in the range [0, 2\pi]
        M = angle_in_range(M);

        // Solve Kepler equation through a basic Newton method
        Real dE{0.0}, E{M}, e{kepl.e()};
        for (Integer k{0}; k < 100; ++k)
        {
          dE = (E - e * std::sin(E) - M) / (1.0 - e * std::cos(E));
          // Saturate if steps are too largex
          E -= std::max(std::min(dE, 0.1), -0.1);;
          // Break if the error is small enough
          if (std::abs(dE) < EPSILON_HIGH || std::abs(E) < EPSILON_HIGH) {break;}
        }

        ASTRO_ASSERT(std::abs(dE) < EPSILON_HIGH, CMD "convergence not reached: E = " << E <<
          ", dE = " << dE << ", M = " << M << ", e = " << e << ".");

        return E;

        #undef CMD
      }

      /**
      * Compute the mean anomaly \f$ M \f$ from the mean anomaly \f$ M \f$ as
      * \f$ M = e\sinh(H) - H \text{.} \f$. The solution is found through a basic Newton
      * method.
      *
      * \note The mean anomaly \f$ M \f$ must not be clamped in the range \f$ [0, 2\pi] \f$!
      * \param[in] M The mean anomaly \f$ M \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \return The mean anomaly \f$ M \f$.
      */
      Real M_to_H(Real M, Keplerian const & kepl)
      {

        #define CMD "Astro::OrbitalElements::Anomaly::M_to_H(...): "

        // Solve Kepler equation through a basic Newton method
        Real e{kepl.e()}, abs_M{M > 0 ? M : -M};
        Real dH{0.0}, H{5.0*e-2.5 > abs_M ? std::pow(6.0*abs_M/e, 1.0/3.0) : std::log(2.0*abs_M/e)};
        for (Integer k{0}; k < 100; ++k)
        {
          dH = (e * std::sinh(H) - H - abs_M) / (e * std::cosh(H) - 1.0);
          H -= dH;
          // Break if the error is small enough
          if (std::abs(dH) < EPSILON_HIGH || std::abs(H) < EPSILON_HIGH) {break;}
        }

        ASTRO_ASSERT(std::abs(dH) < EPSILON_HIGH, CMD "convergence not reached: H = " << H <<
          ", dH = " << dH << ", M = " << M << ", e = " << e << ".");

        return H > 0 ? H : -H; // Return the positive value

        #undef CMD
      }

      /**
      * Compute the mean longitude \f$ \lambda \f$ from the mean anomaly \f$ M \f$ as
      * \f[ \lambda = M + \omega + I\Omega \text{.} \f]
      * \param[in] M The mean anomaly \f$ M \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \param[in] I The posigrade (+1)/retrograde (-1) factor.
      * \return The mean longitude \f$ \lambda \f$.
      */
      Real M_to_lambda(Real M, Keplerian const & kepl, Factor I)
      {
        return M + kepl.omega() + kepl.Omega()*static_cast<Real>(I);
      }

      /**
      * Compute the true anomaly \f$ \nu \f$ from the eccentric anomaly \f$ E \f$ as
      * \f[ \tan\left(\frac{\nu - E}{2}\right) = \displaystyle\frac{\sin(E)}{\beta - \cos(E)} \text{,} \f]
      * with \f$ \beta = \displaystyle\frac{1 + \sqrt{1 - e^2}}{e} \f$. That is equivalent to
      * \f[ \nu = E + 2\arctan\left(\frac{\sin(E)}{\beta - \cos(E)}\right) \text{.} \f]
      *
      \note Refer to *"A note on the relations between true and eccentric anomalies in the two-body
      * problem"* by R. Broucke and P. Cefola, Celestial Mechanics, Vol. 7, pp. 300-389, 1973.
      * \param[in] E The eccentric anomaly \f$ E \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \return The true anomaly \f$ \nu \f$.
      */
      Real E_to_nu(Real E, Keplerian const & kepl)
      {
        Real beta{(1.0 + std::sqrt(1.0 - kepl.e()*kepl.e())) / kepl.e()};
        return E + 2.0 * std::atan(std::sin(E) / (beta - std::cos(E)));
      }

      /**
      * Compute the mean anomaly \f$ M \f$ from the eccentric anomaly \f$ E \f$ as
      * \f[ M = E - e\sin(E) \text{.} \f]
      * \param[in] E The eccentric anomaly \f$ E \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \return The true anomaly \f$ \nu \f$.
      */
      Real E_to_M(Real E, Keplerian const & kepl) {return E - kepl.e() * std::sin(E);}

      /**
      * Compute the true anomaly \f$ \nu \f$ from the hyperbolic anomaly \f$ H \f$ as
      * \f[ \nu = 2\arctan\left(\sqrt{\frac{e + 1}{e - 1}}\tanh\left(\frac{H}{2}\right)\right) \text{.} \f]
      *
      \note Refer to (p. 167) *"An Introduction to the Mathematics and Methods of Astrodynamics"* by
      * R. H. Battin, AIAA Education Series, 1999.
      * \param[in] H The hyperbolic anomaly \f$ H \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \return The true anomaly \f$ \nu \f$.
      */
      Real H_to_nu(Real H, Keplerian const & kepl)
      {
        return 2.0*std::atan(std::sqrt((kepl.e()+1.0)/(kepl.e()-1.0)) * std::tanh(H/2.0));
      }

      /**
      * Compute the mean anomaly \f$ M \f$ from the hyperbolic anomaly \f$ H \f$ as
      * \f[ M = e\sinh(H) - H \text{.} \f]
      * \param[in] H The hyperbolic anomaly \f$ H \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \return The mean anomaly \f$ M \f$.
      */
      Real H_to_M(Real H, Keplerian const & kepl) {return kepl.e() * std::sinh(H) - H;}

      /**
      * Compute the true anomaly \f$ \nu \f$ from the true longitude \f$ L \f$ as
      * \f[ \nu = L - \omega - I\Omega \text{.} \f]
      * \param[in] L The true longitude \f$ L \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \param[in] I The posigrade (+1)/retrograde (-1) factor.
      * \return The true anomaly \f$ \nu \f$.
      */
      Real L_to_nu(Real L, Keplerian const & kepl, Factor I)
      {
        return L - kepl.omega() - kepl.Omega()*static_cast<Real>(I);
      }

      /**
      * Compute the mean longitude \f$ \lambda \f$ from the true longitude \f$ L \f$ as
      * \f[ \lambda = L - \nu + M \text{.} \f]
      * \param[in] L The true longitude \f$ L \f$.
      * \param[in] nu The true anomaly \f$ \nu \f$.
      * \param[in] M The mean anomaly \f$ M \f$.
      * \return The mean longitude \f$ \lambda \f$.
      */
      Real L_to_lambda(Real L, Real nu, Real M) {return L - nu + M;}

      /**
      * Compute the mean anomaly \f$ M \f$ from the mean longitude \f$ \lambda \f$ as
      * \f[ M = \lambda - \omega - I\Omega \text{.} \f]
      * \param[in] lambda The mean longitude \f$ \lambda \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \param[in] I The posigrade (+1)/retrograde (-1) factor.
      * \return The mean anomaly \f$ M \f$.
      */
      Real lambda_to_M(Real lambda, Keplerian const & kepl, Factor I)
      {
        return lambda - kepl.omega() - static_cast<Real>(I)*kepl.Omega();
      }

      /**
      * Compute the true longitude \f$ L \f$ from the mean longitude \f$ \lambda \f$ as
      * \f[ L = \lambda + \nu - M \text{.} \f]
      * \param[in] lambda The mean longitude \f$ \lambda \f$.
      * \param[in] nu The true anomaly \f$ \nu \f$.
      * \param[in] M The mean anomaly \f$ M \f$.
      * \return The true longitude \f$ L \f$.
      */
      Real lambda_to_L(Real lambda, Real nu, Real M) {return lambda + nu - M;}

    } // namespace Convert

    /*\
     |      _                                _
     |     / \   _ __   ___  _ __ ___   __ _| |_   _
     |    / _ \ | '_ \ / _ \| '_ ` _ \ / _` | | | | |
     |   / ___ \| | | | (_) | | | | | | (_| | | |_| |
     |  /_/   \_\_| |_|\___/|_| |_| |_|\__,_|_|\__, |
     |                                         |___/
    \*/

    /**
    * \brief Class container for the orbital anomalies.
    *
    * Class container for the orbital anomalies, which is made of the following anomalies:
    *   - the true anomaly \f$ \nu \f$ (rad),
    *   - the mean anomaly \f$ M \f$ (rad),
    *   - the eccentric anomaly \f$ E \f$ (rad).
    *   - the hyperbolic anomaly \f$ H \f$ (rad).
    *   - the true longitude \f$ L \f$ (rad).
    *   - the mean longitude \f$ \lambda \f$ (rad).
    */
    class Anomaly
    {
      Real m_nu{QUIET_NAN};     /**< True anomaly \f$ \nu \f$ (rad). */
      Real m_M{QUIET_NAN};      /**< Mean anomaly \f$ M \f$ (rad). */
      Real m_E{QUIET_NAN};      /**< Eccentric anomaly \f$ E \f$ (rad). */
      Real m_H{QUIET_NAN};      /**< Hyperbolic anomaly \f$ H \f$ (rad). */
      Real m_L{QUIET_NAN};      /**< True longitude \f$ L \f$ (rad). */
      Real m_lambda{QUIET_NAN}; /**< Mean longitude \f$ \lambda \f$ (rad). */

    public:

      /**
      * Class constructor for the orbit anomalies.
      */
      Anomaly() {}

      /**
      * Enable the default orbit anomalies copy constructor.
      */
      Anomaly(Anomaly const &) = default;

      /**
      * Enable the default orbit anomalies move constructor.
      */
      Anomaly(Anomaly &&) = default;

      /**
      * Enable the default orbit anomalies assignment operator.
      */
      Anomaly & operator=(const Anomaly &) = default;

      /**
      * Enable the default orbit anomalies move assignment operator.
      */
      Anomaly & operator=(Anomaly &&) = default;

      /**
      * Get the true anomaly \f$ \nu \f$.
      * \return The true anomaly \f$ \nu \f$.
      */
      Real nu() const {return this->m_nu;}

      /**
      * Set the true anomaly \f$ \nu \f$ and compute the other anomalies accordingly.
      * \param[in] nu The true anomaly \f$ \nu \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \param[in] I The posigrade (+1)/retrograde (-1) factor.
      */
      void nu(Real nu, Keplerian const & kepl, Factor I)
      {
        this->m_nu     = nu;
        this->m_M      = Convert::nu_to_M(nu, kepl);
        this->m_E      = Convert::nu_to_E(nu, kepl);
        this->m_H      = Convert::M_to_H(this->m_M, kepl);
        this->m_L      = Convert::nu_to_L(nu, kepl, I);
        this->m_lambda = Convert::M_to_lambda(this->m_M, kepl, I);
      }

      /**
      * Get the mean anomaly \f$ M \f$.
      * \return The mean anomaly \f$ M \f$.
      */
      Real M() const {return this->m_M;}

      /**
      * Set the mean anomaly \f$ M \f$.
      * \param[in] t_M The mean anomaly \f$ M \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \param[in] I The posigrade (+1)/retrograde (-1) factor.
      */
      void M(Real t_M, Keplerian const & kepl, Factor I)
      {
        this->m_E      = Convert::M_to_E(t_M, kepl);
        this->m_nu     = Convert::E_to_nu(this->m_E, kepl);
        this->m_M      = t_M;
        this->m_H      = Convert::M_to_H(t_M, kepl);
        this->m_L      = Convert::nu_to_L(this->m_nu, kepl, I);
        this->m_lambda = Convert::M_to_lambda(this->m_M, kepl, I);
      }

      /**
      * Get the eccentric anomaly \f$ E \f$.
      * \return The eccentric anomaly \f$ E \f$.
      */
      Real E() const {return this->m_E;}

      /**
      \brief Set the eccentric anomaly \f$ E \f$.
      *
      * Set the eccentric anomaly \f$ E \f$, and compute the mean anomaly \f$ M \f$ and the true anomaly
      * \f$ \nu \f$ accordingly.
      *
      \note Refer to *"A note on the relations between true and eccentric anomalies in the two-body
      * problem"* by R. Broucke and P. Cefola, Celestial Mechanics, Vol. 7, pp. 300-389, 1973.
      * \param[in] t_E The eccentric anomaly \f$ E \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \param[in] I The posigrade (+1)/retrograde (-1) factor.
      */
      void E(Real t_E, Keplerian const & kepl, Factor I)
      {
        this->m_nu     = Convert::E_to_nu(t_E, kepl);
        this->m_M      = Convert::E_to_M(t_E, kepl);
        this->m_E      = t_E;
        this->m_H      = Convert::M_to_H(this->m_M, kepl);
        this->m_L      = Convert::nu_to_L(this->m_nu, kepl, I);
        this->m_lambda = Convert::M_to_lambda(this->m_M, kepl, I);
      }

      /**
      * Get the hyperbolic anomaly \f$ H \f$.
      * \return The hyperbolic anomaly \f$ H \f$.
      */
      Real H() const {return this->m_H;}

      /**
      * \brief Set the hyperbolic anomaly \f$ H \f$.
      *
      * Set the hyperbolic anomaly \f$ H \f$, and compute the mean anomaly \f$ M \f$ and the true anomaly
      * \f$ \nu \f$ accordingly.
      * \param[in] t_H The hyperbolic anomaly \f$ H \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \param[in] I The posigrade (+1)/retrograde (-1) factor.
      */
      void H(Real t_H, Keplerian const & kepl, Factor I)
      {
        #define CMD "Astro::OrbitalElements::Anomaly::H(...): "

        this->nu(Convert::H_to_nu(t_H, kepl), kepl, I);

        ASTRO_ASSERT(std::abs(this->m_H - t_H) < EPSILON_HIGH,
          CMD "conversion error, H = " << t_H << " ≠ " << this->m_H << ".");

        #undef CMD
      }

      /**
      * Get the true longitude \f$ L \f$.
      * \return The true longitude \f$ L \f$.
      */
      Real L() const {return this->m_L;}

      /**
      \brief Set the true longitude \f$ L \f$.
      *
      * Set the true longitude \f$ L \f$, and compute the true anomaly \f$ \nu \f$, the mean anomaly
      * \f$ M \f$, and the eccentric anomaly \f$ E \f$ accordingly.
      * \param[in] t_L The true longitude \f$ L \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \param[in] I The posigrade (+1)/retrograde (-1) factor.
      */
      void L(Real t_L, Keplerian const & kepl, Factor I)
      {
        #define CMD "Astro::OrbitalElements::Anomaly::L(...): "

        this->nu(Convert::L_to_nu(t_L, kepl, I), kepl, I);

        ASTRO_ASSERT(std::abs(this->m_L - t_L) < EPSILON_HIGH,
          CMD "conversion error, L = " << t_L << " ≠ " << this->m_L << ".");

        #undef CMD
      }

      /**
      * Get the mean longitude \f$ \lambda \f$.
      * \return The mean longitude \f$ \lambda \f$.
      */
      Real lambda() const {return this->m_lambda;}

      /**
      \brief Set the mean longitude \f$ \lambda \f$.
      *
      * Set the mean longitude \f$ \lambda \f$, and compute the true anomaly \f$ \nu \f$, the mean anomaly
      * \f$ M \f$, and the eccentric anomaly \f$ E \f$ accordingly.
      * \param[in] t_lambda The mean longitude \f$ \lambda \f$.
      * \param[in] kepl The keplerian orbital elements.
      * \param[in] I The posigrade (+1)/retrograde (-1) factor.
      */
      void lambda(Real t_lambda, Keplerian const & kepl, Factor I)
      {
        #define CMD "Astro::OrbitalElements::Anomaly::lambda(...): "

        this->M(Convert::lambda_to_M(t_lambda, kepl, I), kepl, I);

        ASTRO_ASSERT(std::abs(this->m_lambda - t_lambda) < EPSILON_HIGH,
          CMD "conversion error, lambda = " << t_lambda << " ≠ " << this->m_lambda << ".");

        #undef CMD
      }

      /**
      * Print the orbit anomialies on a string.
      * \return The orbit anomialies parameters string.
      */
      std::string info() const {
        std::ostringstream os;
        os <<
          "v : true anomaly       = " << this->m_nu << " (rad) = " << rad_to_deg(this->m_nu) << " (deg)" << std::endl <<
          "M : mean anomaly       = " << this->m_M << " (rad) = " << rad_to_deg(this->m_M) << " (deg)" << std::endl <<
          "E : eccentric anomaly  = " << this->m_E << " (rad) = " << rad_to_deg(this->m_E) << " (deg)" << std::endl <<
          "H : hyperbolic anomaly = " << this->m_H << " (rad) = " << rad_to_deg(this->m_H) << " (deg)" << std::endl <<
          "L : true longitude     = " << this->m_L << " (rad) = " << rad_to_deg(this->m_L) << " (deg)" << std::endl <<
          "λ : mean longitude     = " << this->m_lambda << " (rad) = " << rad_to_deg(this->m_lambda) << " (deg)" << std::endl;
          return os.str();
      }

      /**
      * Print the orbit anomialies on a stream.
      * \param[in,out] os Output stream.
      */
      void info(std::ostream & os) {os << this->info();}

    }; // class Anomaly


    /**
    * Convert the cartesian state (position and velocity) vectors to the keplerian orbital elements.
    * \param[in] cart The cartesian state (position and velocity) vectors.
    * \param[in] mu The gravitational parameter.
    * \param[out] kepl The keplerian orbital elements.
    */
    void cartesian_to_keplerian(Cartesian const & cart, Real mu, Keplerian & kepl)
    {
      #define CMD "Astro::OrbitalElements::Convert::cartesian_to_keplerian(...): "

      // Reset the keplerian orbital elements
      kepl.reset();

      // Compute the orbital momentum vector
      Vector3 h{cart.h()};

      // Compute the eccentricity vector
      Real r_norm{cart.r().norm()};
      Vector3 e_vec{cart.v().cross(h)/mu - cart.r()/r_norm};

      // Determine the vector pointing towards the ascending node
      Vector3 n{Vector3::UnitZ().cross(h)};

      ASTRO_ASSERT(std::abs(h.z()) > EPSILON_MEDIUM,
        CMD "conversion error h = " << h << " ≠ 0.");

      // Compute the eccentricity
      Real e{e_vec.norm()};

      // Compute the true anomaly
      Real e_dot_r{e_vec.transpose()*cart.r()};
      Real nu{std::acos(e_dot_r / (e * r_norm))};
      if (cart.r().transpose()*cart.v() < 0.0) {nu = 2.0*PI - nu;}

      // Compute the inclination
      Real i{std::acos(h.z() / h.norm())};

      // Compute the longitude of the ascending node
      Real n_norm{n.norm()};
      Real Omega{std::acos(n.x() / n_norm)};
      if (n.y() < 0.0) {Omega = 2.0*PI - Omega;}

      // Compute the argument of the periapsis
      Real n_dot_e_vec{n.transpose()*e_vec};
      Real omega{std::acos(n_dot_e_vec / (n_norm * e))};
      if (e_vec.z() < 0.0) {omega = 2.0*PI - omega;}

      // Compute the semi-major axis
      Real a{1.0 / (2.0 / r_norm - cart.v().squaredNorm() / mu)};

      // Assign the keplerian orbital elements
      kepl.a(a);
      kepl.e(e);
      kepl.i(i);
      kepl.Omega(Omega);
      kepl.omega(omega);

      ASTRO_ASSERT(kepl.sanity_check(), CMD "conversion error.");

      #undef CMD
    }

    /**
    * Convert the keplerian orbital elements to the cartesian state (position and velocity) vectors.
    * \param[in] kepl The keplerian orbital elements.
    * \param[in] anom The orbital anomalies.
    * \param[in] mu The gravitational parameter.
    * \param[out] cart The cartesian state (position and velocity) vectors.
    */
    void keplerian_to_cartesian(Keplerian const & kepl, Anomaly const & anom, Real mu, Cartesian & cart)
    {
      #define CMD "Astro::OrbitalElements::Convert::keplerian_to_cartesian(...): "

      // Reset the cartesian state vectors
      cart.reset();

      // Compute the distance to the central body
      Real r_c{kepl.a() * (1.0 - kepl.e()*std::cos(anom.E()))};

      // Obtain the position and velocity vector r_of and v_of in the orbital frame
      Vector2 r_of(
        r_c * std::cos(anom.nu()),
        r_c * std::sin(anom.nu())
        // 0.0 third component
      );
      Vector2 v_of(
        -std::sin(anom.E()),
        std::sqrt(1.0 - kepl.e()*kepl.e()) * std::cos(anom.E())
        // 0.0 third componen
      );
      v_of *= std::sqrt(mu / kepl.a());

      // Compute the position and velocity vectors in the inertial frame
      Real c_Omega{std::cos(kepl.Omega())}, s_Omega{std::sin(kepl.Omega())};
      Real c_omega{std::cos(kepl.omega())}, s_omega{std::sin(kepl.omega())};
      Real c_i{std::cos(kepl.i())}, s_i{std::sin(kepl.i())};
      Vector3 r(
        r_of.x()*(c_Omega*c_omega - s_Omega*s_omega*c_i) - r_of.y()*(c_Omega*s_omega + s_Omega*c_omega*c_i),
        r_of.x()*(s_Omega*c_omega + c_Omega*s_omega*c_i) - r_of.y()*(s_Omega*s_omega - c_Omega*c_omega*c_i),
        r_of.x()*(s_omega*s_i)                           + r_of.y()*(c_omega*s_i)
      );
      cart.r(r);
      Vector3 v(
        v_of.x()*(c_Omega*c_omega - s_Omega*s_omega*c_i) - v_of.y()*(c_Omega*s_omega + s_Omega*c_omega*c_i),
        v_of.x()*(s_Omega*c_omega + c_Omega*s_omega*c_i) - v_of.y()*(s_Omega*s_omega - c_Omega*c_omega*c_i),
        v_of.x()*(s_omega*s_i)                           + v_of.y()*(c_omega*s_i)
      );
      cart.v(v);

      ASTRO_ASSERT(cart.sanity_check(), CMD "conversion error.");

      #undef CMD
    }

    /**
    * Convert the equinoctial orbital elements to the cartesian state (position and velocity) vectors.
    * \param[in] equi The equinoctial orbital elements.
    * \param[in] anom The orbital anomalies.
    * \param[in] mu The gravitational parameter.
    * \param[out] cart The cartesian state (position and velocity) vectors.
    */
    void equinoctial_to_cartesian(Equinoctical const & equi, Anomaly const & anom, Real mu, Cartesian & cart)
    {
      #define CMD "Astro::OrbitalElements::Convert::equinoctial_to_cartesian(...): "

      // Reset the cartesian state vectors
      cart.reset();

      // Compute the distance to the central body
      Real c_L{std::cos(anom.L())}, s_L{std::sin(anom.L())};
      Real alpha_2{equi.h()*equi.h() - equi.k()*equi.k()};
      Real s_2{1.0 + equi.h()*equi.h() + equi.k()*equi.k()};
      Real w{1.0 + equi.f()*c_L + equi.g()*s_L};
      Real r_c{equi.p() / w};

      // Compute the position and velocity vectors in the inertial frame
      Real tmp{r_c/s_2};
      Vector3 r(
        tmp * (c_L + alpha_2*c_L + 2.0*equi.h()*equi.k()*s_L),
        tmp * (s_L - alpha_2*s_L + 2.0*equi.h()*equi.k()*c_L),
        2.0*tmp * (equi.h()*s_L - equi.k()*c_L)
      );
      cart.r(r);
      tmp = -1.0/s_2 * std::sqrt(mu/equi.p());
      Vector3 v(
        tmp * ( s_L + alpha_2*s_L - 2.0*equi.h()*equi.k()*c_L + equi.g() - 2.0*equi.f()*equi.h()*equi.k() + alpha_2*equi.g()),
        tmp * (-c_L + alpha_2*c_L + 2.0*equi.h()*equi.k()*s_L - equi.f() + 2.0*equi.g()*equi.h()*equi.k() - alpha_2*equi.f()),
        2.0*tmp * (-equi.h()*c_L - equi.k()*s_L + equi.f()*equi.h() - equi.g()*equi.k())
      );
      cart.v(v);

      ASTRO_ASSERT(cart.sanity_check(), CMD "conversion error.");

      #undef CMD
    }z

    /**
    * Convert the keplerian orbital elements to the equinoctial orbital elements.
    * \param[in] kepl The keplerian orbital elements.
    * \param[in] I The posigrade (+1)/retrograde (-1) factor.
    * \param[out] equi The equinoctial orbital elements.
    */
    void keplerian_to_equinoctial(Keplerian const & kepl, Equinoctical & equi)
    {
      #define CMD "Astro::OrbitalElements::Convert::keplerian_to_equinoctial(...): "

      // Reset the equinoctial orbital elements
      equi.reset();

      // Compute the semi-latus rectum
      Real p{kepl.p()};

      // Compute the x-axis component of the eccentricity vector in the orbital frame f
      Real omega_plus_Omega{kepl.omega() + static_cast<Real>(I)*kepl.Omega()};
      Real f{kepl.e() * std::cos(omega_plus_Omega)};

      // Compute the y-axis component of the eccentricity vector in the orbital frame g
      Real g{kepl.e() * std::sin(omega_plus_Omega)};

      // Compute the x- and y-axis components of the node vector in the orbital frame h
      Real tmp;
      tmp = std::tan(kepl.i() / 2.0);
      if (I == Factor::RETROGRADE) {tmp = 1.0/tmp;} // 1/tan(x) = cot(x)
      Real h{tmp * std::cos(kepl.Omega())};
      Real k{tmp * std::sin(kepl.Omega())};

      // Assign the equinoctial orbital elements
      equi.p(p);
      equi.f(f);
      equi.g(g);
      equi.h(h);
      equi.k(k);

      ASTRO_ASSERT(equi.sanity_check(), CMD "conversion error.");

      #undef CMD
    }

    /**
    * Convert the equinoctial orbital elements to the keplerian orbital elements.
    * \param[in] equi The equinoctial orbital elements.
    * \param[in] I The posigrade (+1)/retrograde (-1) factor.
    * \param[out] kepl The keplerian orbital elements.
    */
    void equinoctial_to_keplerian(Equinoctical const & equi, Keplerian & kepl)
    {
      #define CMD "Astro::OrbitalElements::Convert::equinoctial_to_keplerian(...): "

      // Reset the keplerian orbital elements
      kepl.reset();

      // Assign the keplerian orbital elements
      kepl.a(equi.a());
      kepl.e(equi.e());
      kepl.i(equi.i());
      kepl.Omega(Omegequi.a());
      kepl.omega(omegequi.a());

      ASTRO_ASSERT(kepl.sanity_check(), CMD "conversion error.");

      #undef CMD
    }

  } // namespace OrbitalElements

  /*\
   |    ___       _     _ _
   |   / _ \ _ __| |__ (_) |_
   |  | | | | '__| '_ \| | __|
   |  | |_| | |  | |_) | | |_
   |   \___/|_|  |_.__/|_|\__|
   |
  \*/
//
//  /**
//  * \brief Orbit class container.
//  *
//  * \includedoc Orbit.md
//  */
//  class Orbit
//  {
//    Keplerian    m_keplerian;    /**< Keplerian orbit parameters. */
//    Equinoctical m_equinoctial;  /**< Equinoctical orbit parameters. */
//    Quaternionic m_quaternionic; /**< Quaternionic orbit parameters. */
//    Anomaly      m_anomaly;      /**< Orbit anomalies. */
//    Type         m_type;         /**< Orbit type. */
//    Factor       m_factor;       /**< Orbit posigrade (+1)/retrograde (-1) factor. */
//
//  public:
//
//    /**
//    * Get the posigrade (+1)/retrograde (-1) factor \f$ I \f$.
//    * \return The posigrade (+1)/retrograde (-1) factor \f$ I \f$.
//    */
//    Factor retrograde() const {return this->m_factor;}
//
//    /**
//    * Set the posigrade (+1)/retrograde (-1) factor \f$ I \f$.
//    * \param[in] t_factor The posigrade (+1)/retrograde (-1) factor \f$ I \f$.
//    */
//    void retrograde(Factor t_factor) {this->m_factor = t_factor;}
//
//    /**
//    * Check if the orbit is posigrade.
//    * \return True if the orbit is posigrade, false otherwise.
//    */
//    bool is_posigrade() const {return this->m_factor == Factor::POSIGRADE;}
//
//    /**
//    * Check if the orbit is retrograde.
//    * \return True if the orbit is retrograde, false otherwise.
//    */
//    bool is_retrograde() const {return this->m_factor == Factor::RETROGRADE;}
//
//    /**
//    * Set the orbit as posigrade.
//    */
//    void set_posigrade() {this->m_factor = Factor::POSIGRADE;}
//
//    /**
//    * Set the orbit as retrograde.
//    */
//    void set_retrograde() {this->m_factor = Factor::RETROGRADE;}
//        // "I : factor             = " << (this->m_factor == Factor::POSIGRADE ?
//        //   "+1 (posigrade)" : "-1 (retrograde)") << std::endl;
//
//  }; // class Orbit

} // namespace Astro

#endif // ASTRO_ORBITAL_ELEMENTS_HXX
