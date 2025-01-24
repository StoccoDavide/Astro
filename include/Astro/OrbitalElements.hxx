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
    using Factor = enum class Factor : Integer {POSIGRADE = 1, RETROGRADE = -1}; /**< Orbit type. */

    /*\
     |     ____           _            _
     |    / ___|__ _ _ __| |_ ___  ___(_) __ _ _ __
     |   | |   / _` | '__| __/ _ \/ __| |/ _` | '_ \
     |   | |__| (_| | |  | ||  __/\__ \ | (_| | | | |
     |    \____\__,_|_|   \__\___||___/_|\__,_|_| |_|
     |
    \*/

    /**
    * \brief Class container for the Cartesian orbital parameters.
    *
    * Class container for the Cartesian orbit parameters, which are:
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
      * Print the Cartesian orbit parameters on a string.
      * \return The Cartesian orbit parameters string.
      */
      std::string info() const {
        std::ostringstream os;
        os <<
          "r = " << this->m_r.transpose() << " (UA)" << std::endl <<
          "v = " << this->m_v.transpose() << " (UA/day)" << std::endl;
          return os.str();
      }

      /**
      * Print the Cartesian orbit parameters on a stream.
      * \param[in,out] os Output stream.
      */
      void info(std::ostream &os) {os << this->info();}

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
    *     - the true anomaly \f$ \nu \f$ (rad),
    *     - the mean anomaly \f$ M \f$ (rad),
    *     - the eccentric anomaly \f$ E \f$ (rad).
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
      Real m_nu{QUIET_NAN};    /**< True anomaly \f$ \nu \f$ (rad). */
      Real m_M{QUIET_NAN};     /**< Mean anomaly \f$ M \f$ (rad). */
      Real m_E{QUIET_NAN};     /**< Eccentric anomaly \f$ E \f$ (rad). */

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
      * \brief Compute the mean anomaly \f$ M \f$ from the eccentric anomaly \f$ E \f$.
      *
      * Compute the mean anomaly \f$ M \f$ from the eccentric anomaly \f$ E \f$ as
      * \f[ M = \nu - 2\arctan\left(\frac{\sin(\nu)}{\beta + \cos(\nu)}\right) -
      * \frac{e\sqrt{1 - e^2}\sin(\nu)}{1 + e\cos(\nu)} \text{,} \f]
      * with \f$ \beta = \frac{1 + \sqrt{1 - e^2}}{e} \f$. That is equivalent to
      *
      * \note Refer to *"A note on the relations between true and eccentric anomalies in the two-body
      * problem"* by R. Broucke and P. Cefola, Celestial Mechanics, Vol. 7, pp. 300-389, 1973.
      * \param[in] t_nu The true anomaly \f$ \nu \f$.
      */
      Real nu_to_M(Real t_nu) const
      {
        Real beta{(1.0 + std::sqrt(1.0 - this->m_e*this->m_e)) / this->m_e};
        return t_nu - 2.0 * std::atan(std::sin(t_nu) / (beta + std::cos(t_nu))) -
          (this->m_e * std::sqrt(1.0 - this->m_e*this->m_e) * std::sin(t_nu)) /
          (1.0 + this->m_e * std::cos(t_nu));
      }

      /**
      * \brief Compute the true anomaly \f$ \nu \f$ from the eccentric anomaly \f$ E \f$.
      *
      * Compute the true anomaly \f$ \nu \f$ from the eccentric anomaly \f$ E \f$ as
      * \f[ \tan\left(\frac{\nu - E}{2}\right) = \frac{\sin(\nu)}{\beta + \cos(\nu)} \text{,} \f]
      * with \f$ \beta = \frac{1 + \sqrt{1 - e^2}}{e} \f$. That is equivalent to
      * \f[ E = -\nu - 2\arctan\left(\frac{\sin(\nu)}{\beta + \cos(\nu)}\right) \text{.} \f]
      *
      \note Refer to *"A note on the relations between true and eccentric anomalies in the two-body
      * problem"* by R. Broucke and P. Cefola, Celestial Mechanics, Vol. 7, pp. 300-389, 1973.
      * \param[in] t_nu The true anomaly \f$ \nu \f$.
      */
      Real nu_to_E(Real t_nu) const
      {
        Real beta{(1.0 + std::sqrt(1.0 - this->m_e*this->m_e)) / this->m_e};
        return -t_nu - 2.0 * std::atan(std::sin(t_nu) / (beta + std::cos(t_nu)));
      }

      /**
      * Get the true anomaly \f$ \nu \f$.
      * \return The true anomaly \f$ \nu \f$.
      */
      Real nu() const {return this->m_nu;}

      /**
      * Set the true anomaly \f$ \nu \f$.
      * \param[in] t_nu The true anomaly \f$ \nu \f$.
      */
      void nu(Real t_nu)
      {
        this->m_nu = t_nu;
        this->m_M  = this->nu_to_M(t_nu);
        this->m_E  = this->nu_to_E(t_nu);
      }

      /**
      * Get the mean anomaly \f$ M \f$.
      * \return The mean anomaly \f$ M \f$.
      */
      Real M() const {return this->m_M;}

      /**
      * \brief Compute the mean anomaly \f$ M \f$ from the true anomaly \f$ \nu \f$.
      *
      * Compute the mean anomaly \f$ M \f$ from the true anomaly \f$ \nu \f$ by solving the Kepler
      * equation 1f$ M = E - e\sin(E) \f$ for \f$ E \f$. The solution is found through a basic Newton
      * method. Then, the mean anomaly is computed as \f$ M = E - e\sin(E) \f$.
      * \param[in] t_M The mean anomaly \f$ M \f$.
      * \param[in] t_E The eccentric anomaly \f$ E \f$ (default is computed from \f$ M \f$).
      */
      Real M_to_nu(Real t_M, Real t_E = FIXME) const
      {
        return this->E_to_nu(t_E);
      }

      /**
      * \brief Compute the eccentric anomaly \f$ E \f$ from the mean anomaly \f$ M \f$.
      *
      * Compute the eccentric anomaly \f$ E \f$ from the mean anomaly \f$ M \f$ through a the solution
      * of the nonlinear equation \f$ E = M + e\sin(E) \f$ for \f$ E \f$. The solution is found through
      * a basic Newton method.
      * \param[in] t_M The mean anomaly \f$ M \f$.
      */
      Real M_to_E(Real t_M) const
      {

        #define CMD "Astro::OrbitalElements::Keplerian::M_to_E(...): "

        angle_in_range(t_M);

        // Solve Kepler equation through a basic Newton method
        Real dE{0.0}, E{t_M};
        for (Integer k{0}; k < 100; ++k)
        {
          dE = (E - this->m_e * std::sin(E) - t_M) / (1.0 - this->m_e * std::cos(E));
          // Saturate if steps are too large
          dE = std::min(dE,  0.1);
          dE = std::max(dE, -0.1);
          E -= dE;
          // Break if the error is small enough
          if (std::abs(dE) < EPSILON_LOW) {break;}
          // Covergence should be reached at E ~ 1
        }

        ASTRO_ASSERT(std::abs(dE) < EPSILON_MEDIUM,
          CMD "convergence not reached: E = " << E << ", dE = " << dE << ", M = " << t_M << ", e = "
          << this->m_e << ".");

        return E;

        #undef CMD
      }

      /**
      * Set the mean anomaly \f$ M \f$.
      * \param[in] t_M The mean anomaly \f$ M \f$.
      */
      void M(Real t_M)
      {
        this->m_M  = t_M;
        this->m_E  = this->M_to_E(t_M);
        this->m_nu = this->M_to_nu(t_M, this->m_E);
      }

      /**
      * Get the eccentric anomaly \f$ E \f$.
      * \return The eccentric anomaly \f$ E \f$.
      */
      Real E() const {return this->m_E;}

      /**
      * \brief Compute the true anomaly \f$ \nu \f$ from the eccentric anomaly \f$ E \f$.
      *
      * Compute the true anomaly \f$ \nu \f$ from the eccentric anomaly \f$ E \f$ as
      * \f[ \tan\left(\frac{\nu - E}{2}\right) = \frac{\sin(E)}{\beta - \cos(E)} \text{,} \f]
      * with \f$ \beta = \frac{1 + \sqrt{1 - e^2}}{e} \f$. That is equivalent to
      * \f[ \nu = E + 2\arctan\left(\frac{\sin(E)}{\beta - \cos(E)}\right) \text{.} \f]
      *
      \note Refer to *"A note on the relations between true and eccentric anomalies in the two-body
      * problem"* by R. Broucke and P. Cefola, Celestial Mechanics, Vol. 7, pp. 300-389, 1973.
      * \param[in] t_E The eccentric anomaly \f$ E \f$.
      */
      Real E_to_nu(Real t_E) const
      {
        Real beta{(1.0 + std::sqrt(1.0 - this->m_e*this->m_e)) / this->m_e};
        return t_E + 2.0 * std::atan(std::sin(t_E) / (beta - std::cos(t_E)));
      }

      /**
      * \brief Compute the mean anomaly \f$ M \f$ from the eccentric anomaly \f$ E \f$.
      *
      * Compute the mean anomaly \f$ M \f$ from the eccentric anomaly \f$ E \f$ as
      * \f[ M = E - e\sin(E) \text{.} \f]
      * \param[in] t_E The eccentric anomaly \f$ E \f$.
      */
      Real E_to_M(Real t_E) const {return t_E - this->m_e * std::sin(t_E);}

      /**
      \brief Set the eccentric anomaly \f$ E \f$.
      *
      * Set the eccentric anomaly \f$ E \f$, and compute the mean anomaly \f$ M \f$ and the true anomaly
      * \f$ \nu \f$ accordingly.
      *
      \note Refer to *"A note on the relations between true and eccentric anomalies in the two-body
      * problem"* by R. Broucke and P. Cefola, Celestial Mechanics, Vol. 7, pp. 300-389, 1973.
      * \param[in] t_E The eccentric anomaly \f$ E \f$.
      */
      void E(Real t_E)
      {
        this->m_E  = t_E;
        this->m_M  = this->E_to_M(t_E);
        this->m_nu = this->E_to_nu(t_E);
      }

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
          "ω : arg. of periapsis = " << this->m_omega << " (rad) = " << rad_to_deg(this->m_omega) << " (deg)" << std::endl <<
          "θ : true anomaly      = " << this->m_nu << " (rad) = " << rad_to_deg(this->m_nu) << " (deg)" << std::endl <<
          "M : mean anomaly      = " << this->m_M << " (rad) = " << rad_to_deg(this->m_M) << " (deg)" << std::endl <<
          "E : eccentric anomaly = " << this->m_E << " (rad) = " << rad_to_deg(this->m_E) << " (deg)" << std::endl;
          return os.str();
      }

      /**
      * Print the Keplerian orbit parameters on a stream.
      * \param[in,out] os Output stream.
      */
      void info(std::ostream &os) {os << this->info();}

      /**
      * Check if the Keplerian orbit parameters are nonsingular (*i.e.*, far from  singularities).
      * \param[in] tol Tolerance for the singulary check.
      * \return True if the Keplerian orbit parameters are nonsingular, false otherwise.
      */
      bool is_nonsingular(Real tol = EPSILON_LOW) const {
        return this->m_e > tol && this->m_e < 1.0 - tol &&
               this->m_i > tol && this->m_i < PI - tol;
      }

      /**
      * Check if the Keplerian orbit parameters are singular (*i.e.*, close to singularities).
      * \param[in] tol Tolerance for the singulary check.
      * \return True if the Keplerian orbit parameters are singular, false otherwise.
      */
      bool is_singular(Real tol = EPSILON_LOW) const {return !this->is_nonsingular(tol);}

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
    *   - the true longitude \f$ L \f$ (rad).
    *   - the posigrade/retrograde factor \f$ I = +1 \f$ (posigrade) or \f$ I = -1 \f$ (retrograde) (-).
    * Singular at \f$ i = \pi \f$.
    *
    * \note For more information on orbital elements, refer to *"Survey of Orbital Elements"*, by
    * G. R. Hintz, Journal of Guidance, Control, and Dynamics, Vol. 31, No. 3, May-June 2008.
    */
    class Equinoctical
    {
      Real   m_p{QUIET_NAN};              /**< Semi-latus rectum \f$ p \f$ (UA). */
      Real   m_f{QUIET_NAN};              /**< \f$ X \f$-axis component of the eccentricity vector in the orbital frame (-). */
      Real   m_g{QUIET_NAN};              /**< \f$ Y \f$-axis component of the eccentricity vector in the orbital frame (-). */
      Real   m_h{QUIET_NAN};              /**< \f$ X \f$-axis component of the node vector in the orbital frame (-). */
      Real   m_k{QUIET_NAN};              /**< \f$ Y \f$-axis component of the node vector in the orbital frame (-). */
      Real   m_true_longitude{QUIET_NAN}; /**< True longitude \f$ L \f$ (rad). */
      Factor m_factor{Factor::POSIGRADE}; /**< Posigrade/retrograde factor (-). */

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
      * \param[in] t_factor The posigrade/retrograde factor.
      */
      Equinoctical(Real t_p, Real t_f, Real t_g, Real t_h, Real t_k, Factor t_factor = Factor::POSIGRADE)
        : m_p(t_p), m_f(t_f), m_g(t_g), m_h(t_h), m_k(t_k), m_factor(t_factor) {}

      /**
      * Enable the default Equinoctical orbit parameters copy constructor.
      */
      Equinoctical(Equinoctical const &) = default;

      /**
      * Enable the default Equinoctical orbit parameters move constructor.
      */
      Equinoctical(Equinoctical &&) = default;

      /**
      * Enable the default Equinoctical orbit parameters assignment operator.
      */
      Equinoctical & operator=(const Equinoctical &) = default;

      /**
      * Enable the default Equinoctical orbit parameters move assignment operator.
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
      * Get the true longitude \f$ L \f$.
      * \return The true longitude \f$ L \f$.
      */
      Real true_longitude() const {return this->m_true_longitude;}

      /**
      * Set the true longitude \f$ L \f$.
      * \param[in] t_true_long The true longitude \f$ L \f$.
      */
      void true_longitude(Real t_true_long) {this->m_true_longitude = t_true_long;}

      /**
      * Get the posigrade/retrograde factor \f$ I \f$.
      * \return The posigrade/retrograde factor \f$ I \f$.
      */
      Factor retrograde() const {return this->m_factor;}

      /**
      * Set the posigrade/retrograde factor \f$ I \f$.
      * \param[in] t_factor The posigrade/retrograde factor \f$ I \f$.
      */
      void retrograde(Factor t_factor) {this->m_factor = t_factor;}

      /**
      * Check if the orbit is posigrade.
      * \return True if the orbit is posigrade, false otherwise.
      */
      bool is_posigrade() const {return this->m_factor == Factor::POSIGRADE;}

      /**
      * Check if the orbit is retrograde.
      * \return True if the orbit is retrograde, false otherwise.
      */
      bool is_retrograde() const {return this->m_factor == Factor::RETROGRADE;}

      /**
      * Set the orbit as posigrade.
      */
      void set_posigrade() {this->m_factor = Factor::POSIGRADE;}

      /**
      * Set the orbit as retrograde.
      */
      void set_retrograde() {this->m_factor = Factor::RETROGRADE;}

      /**
      * Print the Equinoctical orbit parameters on a string.
      * \return Equinoctical orbit parameters string.
      */
      std::string info() const {
        std::ostringstream os;
        os <<
          "p : semi-latus rectum  = " << this->m_p << " (UA)" << std::endl <<
          "f : x-axis ecc. vector = " << this->m_f << " (-)" << std::endl <<
          "g : y-axis ecc. vector = " << this->m_g << " (-)" << std::endl <<
          "h : x-axis node vector = " << this->m_h << " (-)" << std::endl <<
          "k : y-axis node vector = " << this->m_k << " (-)" << std::endl <<
          "L : true longitude     = " << this->m_true_longitude << " (rad)" << std::endl <<
          "I = " << (this->m_factor == Factor::POSIGRADE ? "POSIGRADE" : "RETROGRADE") << std::endl;
          return os.str();
      }

      /**
      * Print the Equinoctical orbit parameters on a stream.
      * \param[in,out] os Output stream.
      */
      void info(std::ostream &os) {os << this->info();}

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
      * Class constructor for quaternionic orbit parameters.
      */
      Quaternionic() {}

      /**
      * Class constructor for quaternionic orbit parameters.
      * \param[in] t_q The quaternionic orbit parameters.
      */
      Quaternionic(Vector4 const &t_q) : m_q(t_q) {}

      /**
      * Class constructor for quaternionic orbit parameters.
      * \param[in] t_q_1 The first quaternionic orbit parameter.
      * \param[in] t_q_2 The second quaternionic orbit parameter.
      * \param[in] t_q_3 The third quaternionic orbit parameter.
      * \param[in] t_q_4 The fourth quaternionic orbit parameter.
      */
      Quaternionic(Real t_q_1, Real t_q_2, Real t_q_3, Real t_q_4) : m_q(t_q_1, t_q_2, t_q_3, t_q_4) {}

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
          "q⁴ : 4yh parameter = " << this->m_q(3) << " (-)" << std::endl;
          return os.str();
      }

      /**
      * Print the quaternionic orbit parameters on a stream.
      * \param[in,out] os Output stream.
      */
      void info(std::ostream &os) {os << this->info();}

    }; // class Quaternionic


    //! Compute the Keplerian semi-major axis \f$ a \f$ from Equinoctical elements.
    //! \return Semi-major axis \f$ a \f$ (UA).
    Real keplerian_a() const;

    //! Compute the Keplerian eccentricity \f$ e \f$ from Equinoctical elements.
    //! \return Eccentricity \f$ e \f$ (-).
    Real keplerian_e() const;

    //! Compute the Keplerian inclination \f$ i \f$ from Equinoctical elements.
    //! \return Inclination \f$ i \f$ (rad).
    Real keplerian_i() const;

    //! Compute the Keplerian longitude of the ascending node \f$ \Omega \f$ from
    //! Equinoctical elements.
    Real keplerian_Omega() const;

    //! Compute the Keplerian longitude of the ascending node \f$ \omega \f$ from
    //! Equinoctical elements.
    //! \return Longitude of the ascending node \f$ \omega \f$ (rad).
    Real keplerian_omega() const;

    //! Compute the Keplerian angular momentum vector \f$ \mathbf{L} \f$ from
    //! Equinoctical elements.
    //! \param mu_s Standard gravitational parameter of the Sun.
    //! \param L Angular momentum vector \f$ \mathbf{L} \f$.
    void keplerian_L_vec(Real mu_s, Vector3 L) const;

    //! Compute the Keplerian vector \f$ \mathbf{A} \f$ from Equinoctical elements.
    //! \param mu_s Standard gravitational parameter of the Sun.
    //! \param A Vector \f$ \mathbf{A} \f$.
    void keplerian_A_vec(Real mu_s, Vector3 A) const;

    //! Compute the orbit energy.
    //! \param[in] mu_s Standard gravitational parameter of the Sun.
    //! \return Orbit energy.
    Real orbit_energy(Real mu_s) const{
      return mu_s * (this->m_g * this->m_g + this->m_f * this->m_f - 1.0) / this->m_p/2.0;
    }

  #if 0
    // da rivedere
    void
    invariantLAEtoEquinoctical(
      Real const L[3],
      Real const A[3],
      Real       E,
      Real       muS,
      Equinoctical &   EQ
    );
  #endif

    Real
    equinoctial_to_apoapsis(Equinoctical const &EQ)
    {
      Real e = hypot(EQ.f, EQ.g);
      return EQ.p / (1 - e);
    }

    Real
    equinoctial_to_periapsis(Equinoctical const &EQ)
    {
      Real e = hypot(EQ.f, EQ.g);
      return EQ.p / (1 + e);
    }

    Real to_cartesian_x(Real L);
    Real to_cartesian_y(Real L);
    Real to_cartesian_z(Real L);
    Vector3 to_cartesian_xyz(Real L);
    Real to_cartesian_vx(Real L, Real mu_s);
    Real to_cartesian_vy(Real L, Real mu_s);
    Real to_cartesian_vz(Real L, Real mu_s);
    Vector3 to_cartesian_vxyz(Real L, Real mu_s);

    Real
    equinoctial_to_radius(Equinoctical const &EQ, Real L);
    Real
    equinoctial_to_velocity(Equinoctical const &EQ, Real L, Real mu_s);

  } // namespace OrbitalElements

  /*\
   |    ___       _     _ _
   |   / _ \ _ __| |__ (_) |_
   |  | | | | '__| '_ \| | __|
   |  | |_| | |  | |_) | | |_
   |   \___/|_|  |_.__/|_|\__|
   |
  \*/

  /**
  * \brief Orbit class container.
  *
  * \includedoc Orbit.md
  */
  class Orbit
  {
    using Type = enum class Type : integer {HYPERBOLIC = 0, ELLIPTIC = 1, PARABOLIC = 2};

    Keplerian    m_keplerian;
    Equinoctical m_equinoctial;
    Quaternionic m_quaternionic;
    Anomaly      m_anomaly;
    Type         m_type;
    bool         m_retrograde;

  public:

  }; // class Orbit

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    Real
    mean_anomaly_to_E(Real M, Real e);
    Real
    mean_anomaly_to_H(Real M, Real e);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    Real
    eccentric_anomaly_to_true_anomaly(Real E, Real e);
    Real
    E_to_true_anomaly(Real E, Real e);
    Real
    H_to_true_anomaly(Real H, Real e);
    Real
    true_anomaly_to_M(Real theta, Real e);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    from_Keplerian_to_Equinoctical(
      Keplerian const &K,
      Equinoctical     &EQ);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    from_Keplerian_to_Equinoctical(
      Keplerian const &K,
      Real             theta,
      Equinoctical     &EQ,
      Real            &L);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    from_equinoctial_to_Keplerian(
      Equinoctical const &EQ,
      Keplerian         &K);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    from_equinoctial_to_Keplerian(
      Equinoctical const &EQ,
      Real               L,
      Keplerian         &K,
      Real              &theta);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    equinoctial_to_reference(
      Equinoctical const &EQ,
      Real               f[3],
      Real               g[3],
      Real               w[3]);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    Keplerian_to_reference(
      Keplerian const &K,
      Real             X[3],
      Real             Y[3],
      Real             Z[3]);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    point_and_velocity_to_Equinoctical_and_Keplerian(
      Vector3 const  &P,
      Vector3 const  &V,
      Real         muS,
      Equinoctical &EQ,
      Real        &L,
      Keplerian   &K,
      Real        &theta,
      Real        &M0);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    equinoctial_to_point(
      Equinoctical const &EQ,
      Real               L,
      Real               P[3]);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    equinoctial_to_velocity(
      Equinoctical const &EQ,
      Real               L,
      Real               muS,
      Real               V[3]);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    equinoctial_to_point_and_velocity(
      Equinoctical const &EQ,
      Real               L,
      Real               muS,
      Real               P[3],
      Real               V[3]);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    point_and_velocity_to_Frenet_RTN(
      Vector3 const &P,
      Vector3 const &V,
      Vector3       &Dr,
      Vector3       &Dt,
      Vector3       &Dn);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    equinoctial_to_Frenet_RTN(
      Equinoctical const &EQ,
      Real               L,
      Vector3              &Dr,
      Vector3              &Dt,
      Vector3              &Dn);

    /*
    //                   _ _____
    //    _____   ____ _| |_   _|
    //   / _ \ \ / / _` | | | |
    //  |  __/\ V / (_| | | | |
    //   \___| \_/ \__,_|_| |_|
    */ //Evaluate Trust

    void
    equinoctial_Trtn_to_Txyz(
      Equinoctical const &EQ,
      Real               L,
      Vector3 const        &Trtn,
      Vector3              &Txyz);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    point_and_velocity_Trtn_to_Txyz(
      Vector3 const &P,
      Vector3 const &V,
      Vector3 const &Trtn,
      Vector3       &Txyz);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    equinoctial_Txyz_to_Trtn(
      Equinoctical const &EQ,
      Real               L,
      Vector3 const        &Txyz,
      Vector3              &Trtn);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    point_and_velocity_Txyz_to_Trtn(
      Vector3 const &P,
      Vector3 const &V,
      Vector3 const &Txyz,
      Vector3       &Trtn);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // matrice
    Real // value of b^T
    equinoctial_matrix(
      Equinoctical const &EQ,
      Real               L,
      Real               muS,
      Real               A[6][3]);
    }

} // namespace Astro

#endif // ASTRO_ORBITAL_ELEMENTS_HXX
