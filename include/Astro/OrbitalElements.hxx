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
    * Class container for the Cartesian orbit parameters, which are the position vector \f$
    * \mathbf{r} \f$ and the velocity vector \f$ \mathbf{v} \f$.
    */
    class Cartesian
    {
      Vector3 m_r{NAN_VEC3}; /**< Position vector \f$ \mathbf{r} \f$ (UA). */
      Vector3 m_v{NAN_VEC3}; /**< Velocity vector \f$ \mathbf{v} \f$ (UA/day). */

    public:

      /**
      * Class constructor for Cartesian orbit parameters.
      */
      Cartesian(){}

      /**
      * Class constructor from Cartesian orbit parameters.
      * \param r Position vector \f$ \mathbf{r} \f$.
      * \param v Velocity vector \f$ \mathbf{v} \f$.
      */
      Cartesian(Vector3 const &r, Vector3 const &v) : m_r(r), m_v(v) {}

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
      * \return Position vector \f$ \mathbf{r} \f$.
      */
      Vector3 const & r() const {return this->m_r;}

      /**
      * Set the position vector \f$ \mathbf{r} \f$.
      * \param[in] t_r Position vector \f$ \mathbf{r} \f$.
      */
      void r(Vector3 const & t_r) {this->m_r = t_r;}

      /**
      * Get the velocity vector \f$ \mathbf{v} \f$.
      * \return Velocity vector \f$ \mathbf{v} \f$.
      */
      Vector3 const & v() const {return this->m_v;}

      /**
      * Get the velocity vector \f$ \mathbf{v} \f$.
      * \param[in] t_v Velocity vector \f$ \mathbf{v} \f$.
      */
      void v(Vector3 const & t_v) {this->m_v = t_v;}

      /**
      * Print the orbit parameters.
      * \param[in,out] os Output stream.
      */
      void info(std::ostream &os) const {
        os <<
          "r = " << this->m_r.transpose() << " (UA)" << std::endl <<
          "v = " << this->m_v.transpose() << " (UA/day)" << std::endl;
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
    * \brief Class container for the Keplerian orbital parameters.
    */
    class Keplerian
    {
      Real m_a{QUIET_NAN};        /**< Semi-major axis \f$ a \f$ (UA). */
      Real m_e{QUIET_NAN};        /**< Eccentricity \f$ e \in [0, 1] \f$ (-). */
      Real m_i{QUIET_NAN};        /**< Inclination \f$ i \f$ (rad). */
      Real m_uc_Omega{QUIET_NAN}; /**< Longitude of the ascending node \f$ \Omega \f$ (rad). */
      Real m_lc_omega{QUIET_NAN}; /**< Argument of periapsis \f$ \omega \f$ (rad). */

    public:

      /**
      * Class constructor for Keplerian orbit parameters.
      */
      Keplerian(){}

      /**
      * Class constructor from Keplerian orbit parameters.
      * \param[in] a Semi-major axis \f$ a \f$.
      * \param[in] e Eccentricity \f$ e \f$.
      * \param[in] i Inclination \f$ i \f$.
      * \param[in] uc_Omega Longitude of the ascending node \f$ \Omega \f$.
      * \param[in] lc_omega Argument of periapsis \f$ \omega \f$.
      */
      Keplerian(Real a, Real e, Real i, Real uc_Omega, Real lc_omega)
        : m_a(a), m_e(e), m_i(i), m_uc_Omega(uc_Omega), m_lc_omega(lc_omega)
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
      * \return Semi-major axis \f$ a \f$.
      */
      Real a() const {return this->m_a;}

      /**
      * Set the semi-major axis \f$ a \f$.
      * \param[in] t_a Semi-major axis \f$ a \f$.
      */
      void a(Real t_a) {this->m_a = a;}

      /**
      * Get the eccentricity \f$ e \f$.
      * \return Eccentricity \f$ e \f$.
      */
      Real e() const {return this->m_e;}

      /**
      * Set the eccentricity \f$ e \f$.
      * \param[in] t_e Eccentricity \f$ e \f$.
      */
      void e(Real t_e) {this->m_e = e;}

      /**
      * Get the inclination \f$ i \f$.
      * \return Inclination \f$ i \f$.
      */
      Real i() const {return this->m_e;}

      /**
      * Set the inclination \f$ i \f$.
      * \param[in] t_i Inclination \f$ i \f$.
      */
      void i(Real t_i) {this->m_e = i;}

      /**
      * Get the longitude of the ascending node \f$ \Omega \f$.
      * \return Longitude of the ascending node \f$ \Omega \f$.
      */
      Real uc_Omega() const {return this->m_uc_Omega;}

      /**
      * Set the longitude of the ascending node \f$ \Omega \f$.
      * \param[in] t_uc_Omega Longitude of the ascending node \f$ \Omega \f$.
      */
      void uc_Omega(Real t_uc_Omega) {this->m_uc_Omega = uc_Omega;}

      /**
      * Get the argument of periapsis \f$ \Omega \f$
      * \return Argument of periapsis \f$ \Omega \f$
      */
      Real lc_omega() const {return this->m_lc_omega;}

      /**
      * Set the argument of periapsis \f$ \omega \f$.
      * \param[in] t_lc_omega Argument of periapsis \f$ \omega \f$.
      */
      void lc_omega(Real t_lc_omega) {this->m_lc_omega = lc_omega;}

      /**
      * Print the orbit parameters.
      * \param[in,out] os Output stream.
      */
      void info(std::ostream & os) const {
        os <<
          "a = " << this->m_a << std::endl <<
          "e = " << this->m_e << std::endl <<
          "i = " << this->m_i        << " (rad) = " << radiants_to_degrees(this->m_i)        << " (deg)" << std::endl <<
          "Ω = " << this->m_uc_Omega << " (rad) = " << radiants_to_degrees(this->m_uc_Omega) << " (deg)" << std::endl <<
          "ω = " << this->m_lc_omega << " (rad) = " << radiants_to_degrees(this->m_lc_omega) << " (deg)" << std::endl;
      }

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
    * \brief Class container for the Equinoctical orbital parameters.
    *
    * The Equinoctical orbital parameters are defined as follows:
    * - \f$ p \f$: Semi-latus rectum (UA).
    * - \f$ f \f$: Equinoctical element.
    * - \f$ g \f$: Equinoctical element.
    * - \f$ h \f$: Equinoctical element.
    * - \f$ k \f$: Equinoctical element.
    * - Retrograde flag.
    */
    class Equinoctical
    {
      Real m_p{QUIET_NAN}; /**< Semi-latus rectum \f$ p \f$ (UA). */
      Real m_f{QUIET_NAN}; /**< Equinoctical element \f$ f \f$ (-). */
      Real m_g{QUIET_NAN}; /**< Equinoctical element \f$ g \f$ (-). */
      Real m_h{QUIET_NAN}; /**< Equinoctical element \f$ h \f$ (-). */
      Real m_k{QUIET_NAN}; /**< Equinoctical element \f$ k \f$ (-). */
      bool m_r{false};     /**< Retrograde flag. */

    public:
      /**
      * Class constructor for Equinoctical orbit parameters.
      */
      Equinoctical(){}

      /**
      * Class constructor for Equinoctical orbit parameters.
      * \param p Semi-latus rectum.
      * \param f Equinoctical element.
      * \param g Equinoctical element.
      * \param h Equinoctical element.
      * \param k Equinoctical element.
      * \param r Retrograde flag.
      */
      Equinoctical(Real p, Real f, Real g, Real h, Real k, bool r = false)
        : m_p(p), m_f(f), m_g(g), m_h(h), m_k(k), m_r(r) {}

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
      * Get the semi-latus rectum.
      * \return Semi-latus rectum.
      */
      Real p() const {return this->m_p;}

      /**
      * Set the semi-latus rectum.
      * \param[in] t_p Semi-latus rectum.
      */
      void p(Real t_p) {this->m_p = t_p;}

      /**
      * Get the Equinoctical element \f$ f \f$.
      * \return Equinoctical element \f$ f \f$.
      */
      Real f() const {return this->m_f;}

      /**
      * Set the Equinoctical element \f$ f \f$.
      * \param[in] t_f Equinoctical element \f$ f \f$.
      */
      void f(Real t_f) {this->m_f = t_f;}

      /**
      * Get the Equinoctical element \f$ g \f$.
      * \return Equinoctical element \f$ g \f$.
      */
      Real g() const {return this->m_g;}

      /**
      * Set the Equinoctical element \f$ g \f$.
      * \param[in] t_g Equinoctical element \f$ g \f$.
      */
      void g(Real t_g) {this->m_g = t_g;}

      /**
      * Get the Equinoctical element \f$ h \f$.
      * \return Equinoctical element \f$ h \f$.
      */
      Real h() const {return this->m_h;}

      /**
      * Set the Equinoctical element \f$ h \f$.
      * \param[in] t_h Equinoctical element \f$ h \f$.
      */
      void h(Real t_h) {this->m_h = t_h;}

      /**
      * Get the Equinoctical element \f$ k \f$.
      * \return Equinoctical element \f$ k \f$.
      */
      Real k() const {return this->m_k;}

      /**
      * Set the Equinoctical element \f$ k \f$.
      * \param[in] t_k Equinoctical element \f$ k \f$.
      */
      void k(Real t_k) {this->m_k = t_k;}

      /**
      * Get the Equinoctical retrograde flag.
      * \return Retrograde flag.
      */
      bool retrograde() const {return this->m_r;}

      /**
      * Set the Equinoctical retrograde flag.
      * \param[in] t_r Retrograde flag.
      */
      void retrograde(bool t_r) {this->m_r = t_r;}

      /**
      * Check if the Equinoctical orbit is retrograde.
      * \return True if the Equinoctical orbit is retrograde, false otherwise.
      */
      bool is_retrograde() const {return this->m_r;}

      /**
      * Check if the Equinoctical orbit is posigrade.
      * \return True if the Equinoctical orbit is posigrade, false otherwise.
      */
      bool is_posigrade() const {return !this->m_r;}

      /**
      * Print the Equinoctical orbit parameters.
      * \param[in,out] os Output stream.
      */
      void info(std::ostream & os) const {
        os <<
          "p = " << this->m_p << " (UA)" << std::endl <<
          "f = " << this->m_f << " (-)" << std::endl <<
          "g = " << this->m_g << " (-)" << std::endl <<
          "h = " << this->m_h << " (-)" << std::endl <<
          "k = " << this->m_k << " (-)" << std::endl <<
          "r = " << std::boolalpha << this->m_r << std::endl;
      }

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
    * \brief Class container for the Quaternionic orbital parameters.
    */
    class Quaternionic
    {

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
    Real keplerian_uc_Omega() const;

    //! Compute the Keplerian longitude of the ascending node \f$ \omega \f$ from
    //! Equinoctical elements.
    //! \return Longitude of the ascending node \f$ \omega \f$ (rad).
    Real keplerian_lc_omega() const;

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
   |      _                                _
   |     / \   _ __   ___  _ __ ___   __ _| |_   _
   |    / _ \ | '_ \ / _ \| '_ ` _ \ / _` | | | | |
   |   / ___ \| | | | (_) | | | | | | (_| | | |_| |
   |  /_/   \_\_| |_|\___/|_| |_| |_|\__,_|_|\__, |
   |                                         |___/
  \*/

  /**
  * \brief Anomaly class container.
  *
  * \includedoc Anomaly.md
  */
  struct Anomaly
  {
    Real m_theta{QUIET_NAN}; /**< True anomaly \f$ \theta \f$ (rad). */
    Real M{QUIET_NAN};       /**< Mean anomaly \f$ M \f$ (rad). */
    Real E{QUIET_NAN};       /**< Eccentric anomaly \f$ E \f$ (rad). */
    Real H{QUIET_NAN};       /**< Hyperbolic anomaly \f$ H \f$ (rad). */
  }

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
    true_anomaly_to_mean_anomaly(Real theta, Real e);

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
