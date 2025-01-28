/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Orbit project is distributed under the GNU GPLv3.                     *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef ASTRO_ORBIT_HXX
#define ASTRO_ORBIT_HXX

namespace Orbit
{

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
  * The Orbit class container is used to store the orbital elements of an object in space. The class
  * provides methods for setting and getting the orbital elements in various representations, such as
  * cartesian, keplerian, equinoctial, and quaternionic. The class also provides methods for converting
  * between these representations.
  */
  class Orbit
  {
    using Cartesian    = OrbitalElements::Cartesian;
    using Keplerian    = OrbitalElements::Keplerian;
    using Equinoctial  = OrbitalElements::Equinoctial;
    using Quaternionic = OrbitalElements::Quaternionic;
    using Anomaly      = OrbitalElements::Anomaly;

    Cartesian    m_cart;                      /**< Cartesian orbit parameters. */
    Keplerian    m_kepl;                      /**< Keplerian orbit parameters. */
    Equinoctial  m_equi;                      /**< Equinoctial orbit parameters. */
    Quaternionic m_quat;                      /**< Quaternionic orbit parameters. */
    Anomaly      m_anom;                      /**< Orbit anomalies. */
    Type         m_type{Type::UNDEFINED};     /**< Orbit type. */
    Factor       m_factor{Factor::UNDEFINED}; /**< Orbit posigrade (+1)/retrograde (-1) factor. */
    Real         m_mu{QUIET_NAN};             /**< Gravitational constant of the central body. */

  public:
    /**
    * Class constructor for the astro object.
    */
    Orbit() {}

    /**
    * Enable the default astro copy constructor.
    */
    Orbit(Orbit const &) = default;

    /**
    * Enable the default astro move constructor.
    */
    Orbit(Orbit &&) = default;

    /**
    * Enable the default astro assignment operator.
    */
    Orbit & operator=(const Orbit &) = default;

    /**
    * Enable the default astro move assignment operator.
    */
    Orbit & operator=(Orbit &&) = default;

    /**
    * Get the cartesian orbital elements.
    * \return The cartesian orbital elements.
    */
    Cartesian const & cartesian() const {return this->m_cart;}

    /**
    * Set the cartesian orbital elements.
    * \param[in] r_x Position vector \f$ x \f$-axis component.
    * \param[in] r_y Position vector \f$ y \f$-axis component.
    * \param[in] r_z Position vector \f$ z \f$-axis component.
    * \param[in] v_x Velocity vector \f$ x \f$-axis component.
    * \param[in] v_y Velocity vector \f$ y \f$-axis component.
    * \param[in] v_z Velocity vector \f$ z \f$-axis component.
    */
    void set_cartesian(Real r_x, Real r_y, Real r_z, Real v_x, Real v_y, Real v_z)
    {
      this->m_cart.set(r_x, r_y, r_z, v_x, v_y, v_z);
      OrbitalElements::cartesian_to_keplerian(this->m_cart, this->m_mu, this->m_kepl);
      OrbitalElements::keplerian_to_equinoctial(this->m_kepl, this->m_factor, this->m_equi);
      // TODO: OrbitalElements::cartesian_to_quaternionic(this->m_cart, this->m_quat);
      this->type(this->m_kepl.e());
    }

    /**
    * Set the cartesian orbital elements.
    * \param[in] t_r The position vector \f$ \mathbf{r} \f$.
    * \param[in] t_v The velocity vector \f$ \mathbf{v} \f$.
    */
    void set_cartesian(Vector3 const & t_r, Vector3 const & t_v)
    {
      this->m_cart.set(t_r, t_v);
      OrbitalElements::cartesian_to_keplerian(this->m_cart, this->m_mu, this->m_kepl);
      OrbitalElements::keplerian_to_equinoctial(this->m_kepl, this->m_factor, this->m_equi);
      // TODO: OrbitalElements::cartesian_to_quaternionic(this->m_cart, this->m_quat);
      this->type(this->m_kepl.e());
    }

    /**
    * Set the cartesian orbital elements.
    * \param[in] t_cart The cartesian orbital elements.
    */
    void set_cartesian(Cartesian const & t_cart)
    {
      this->m_cart = t_cart;
      OrbitalElements::cartesian_to_keplerian(this->m_cart, this->m_mu, this->m_kepl);
      OrbitalElements::keplerian_to_equinoctial(this->m_kepl, this->m_factor, this->m_equi);
      // TODO: OrbitalElements::cartesian_to_quaternionic(this->m_cart, this->m_quat);
      this->type(this->m_kepl.e());
    }

    /**
    * Get the keplerian orbital elements.
    * \return The keplerian orbital elements.
    */
    Keplerian const & keplerian() const {return this->m_kepl;}

    /**
    * Set the keplerian orbital elements.
    * \param[in] t_a The semi-major axis \f$ a \f$.
    * \param[in] t_e The eccentricity \f$ e \f$.
    * \param[in] t_i The inclination \f$ i \f$.
    * \param[in] t_Omega The longitude of the ascending node \f$ \Omega \f$.
    * \param[in] t_omega The argument of periapsis \f$ \omega \f$.
    */
    void set_keplerian(Real t_a, Real t_e, Real t_i, Real t_Omega, Real t_omega)
    {
      this->m_kepl.set(t_a, t_e, t_i, t_Omega, t_omega);
      OrbitalElements::keplerian_to_cartesian(this->m_kepl, this->m_anom, this->m_mu, this->m_cart);
      OrbitalElements::keplerian_to_equinoctial(this->m_kepl, this->m_factor, this->m_equi);
      // TODO: OrbitalElements::keplerian_to_quaternionic(this->m_kepl, this->m_quat);
      this->type(this->m_kepl.e());
    }

    /**
    * Set the keplerian orbital elements.
    * \param[in] t_kepl The keplerian orbital elements.
    */
    void set_keplerian(Keplerian const & t_kepl)
    {
      this->m_kepl = t_kepl;
      OrbitalElements::keplerian_to_cartesian(this->m_kepl, this->m_anom, this->m_mu, this->m_cart);
      OrbitalElements::keplerian_to_equinoctial(this->m_kepl, this->m_factor, this->m_equi);
      // TODO: OrbitalElements::keplerian_to_quaternionic(this->m_kepl, this->m_quat);
      this->type(this->m_kepl.e());
    }

    /**
    * Get the equinoctial orbital elements.
    * \return The equinoctial orbital elements.
    */
    Equinoctial const & equinoctial() const {return this->m_equi;}

    /**
      * Set the (modified) equinoctial orbit parameters.
      * \param[in] t_p The semi-latus rectum \f$ p \f$.
      * \param[in] t_f The \f$ x \f$-axis component of the eccentricity vector in the orbital frame.
      * \param[in] t_g The \f$ y \f$-axis component of the eccentricity vector in the orbital frame.
      * \param[in] t_h The \f$ x \f$-axis component of the node vector in the orbital frame.
      * \param[in] t_k The \f$ y \f$-axis component of the node vector in the orbital frame.
      */
    void set_equinoctial(Real t_p, Real t_f, Real t_g, Real t_h, Real t_k)
    {
      this->m_equi.set(t_p, t_f, t_g, t_h, t_k);
      OrbitalElements::equinoctial_to_keplerian(this->m_equi, this->m_kepl);
      OrbitalElements::equinoctial_to_cartesian(this->m_equi, this->m_anom, this->m_mu, this->m_cart);
      // TODO: OrbitalElements::equinoctial_to_quaternionic(this->m_equi, this->m_quat);
      this->type(this->m_kepl.e());
    }

    /**
    * Set the equinoctial orbital elements.
    * \param[in] t_equi The equinoctial orbital elements.
    */
    void set_equinoctial(Equinoctial const & t_equi)
    {
      this->m_equi = t_equi;
      OrbitalElements::equinoctial_to_keplerian(this->m_equi, this->m_kepl);
      OrbitalElements::equinoctial_to_cartesian(this->m_equi, this->m_anom, this->m_mu, this->m_cart);
      // TODO: OrbitalElements::equinoctial_to_quaternionic(this->m_equi, this->m_quat);
      this->type(this->m_kepl.e());
    }

    /**
    * Get the quaternionic orbital elements.
    * \return The quaternionic orbital elements.
    */
    Quaternionic const & quaternionic() const {return this->m_quat;}

    /**
    * Set the quaternionic orbital elements.
    * \param[in] t_q_1 The first quaternionic orbit parameter.
    * \param[in] t_q_2 The second quaternionic orbit parameter.
    * \param[in] t_q_3 The third quaternionic orbit parameter.
    * \param[in] t_q_4 The fourth quaternionic orbit parameter.
    */
    void set_quaternionic(Real t_q_1, Real t_q_2, Real t_q_3, Real t_q_4)
    {
      this->m_quat.set(t_q_1, t_q_2, t_q_3, t_q_4);
      // TODO: OrbitalElements::quaternionic_to_keplerian(this->m_quat, this->m_kepl);
      // TODO: OrbitalElements::quaternionic_to_equinoctial(this->m_quat, this->m_factor, this->m_equi);
      // TODO: OrbitalElements::quaternionic_to_cartesian(this->m_quat, this->m_anom, this->m_mu, this->m_cart);
      this->type(this->m_kepl.e());
    }

    /**
    * Set the quaternionic orbital elements.
    * \param[in] t_quat The quaternion vector.
    */
    void set_quaternionic(Quaternion const & t_quat)
    {
      this->m_quat.set(t_quat);
      // TODO: OrbitalElements::quaternionic_to_keplerian(this->m_quat, this->m_kepl);
      // TODO: OrbitalElements::quaternionic_to_equinoctial(this->m_quat, this->m_factor, this->m_equi);
      // TODO: OrbitalElements::quaternionic_to_cartesian(this->m_quat, this->m_anom, this->m_mu, this->m_cart);
      this->type(this->m_kepl.e());
    }

    /**
    * Set the quaternionic orbital elements.
    * \param[in] t_quat The quaternionic orbital elements.
    */
    void set_quaternionic(Quaternionic const & t_quat)
    {
      this->m_quat = t_quat;
      // TODO: OrbitalElements::quaternionic_to_keplerian(this->m_quat, this->m_kepl);
      // TODO: OrbitalElements::quaternionic_to_equinoctial(this->m_quat, this->m_factor, this->m_equi);
      // TODO: OrbitalElements::quaternionic_to_cartesian(this->m_quat, this->m_anom, this->m_mu, this->m_cart);
      this->type(this->m_kepl.e());
    }

    /**
    * Get the orbital anomalies.
    * \return The orbital anomalies.
    */
    Anomaly const & anomaly() const {return this->m_anom;}

    /**
    * Set the orbital anomalies.
    * \param[in] t_anom The orbital anomalies.
    */
    void anomaly(Anomaly const & t_anom)
    {
      ASTRO_ASSERT(t_anom.sanity_check(),
        "Orbit::Orbit::anomaly(...): invalid orbital anomalies detected.");
      this->m_anom = t_anom;
    }

    /**
    * Get the type of the orbit.
    * \return The type of the orbit.
    */
    Type type() const {return this->m_type;}

    /**
    * Set the type of the orbit.
    * \param[in] e The eccentricity of the orbit.
    */
    void type(Real e)
    {
      if (e > 0.0 && e < 1.0) {
        this->m_type = Type::ELLIPTIC;
      } else if (e == 1.0) {
        this->m_type = Type::PARABOLIC;
      } else if (e > 1.0) {
        this->m_type = Type::HYPERBOLIC;
      } else {
        ASTRO_ERROR("Orbit::Orbit::type(...): invalid eccentricity detected.");
      }
    }

    /**
    * Get the posigrade (+1)/retrograde (-1) factor \f$ I \f$.
    * \return The posigrade (+1)/retrograde (-1) factor \f$ I \f$.
    */
    Factor factor() const {return this->m_factor;}

    /**
    * Set the posigrade (+1)/retrograde (-1) factor \f$ I \f$.
    * \param[in] t_factor The posigrade (+1)/retrograde (-1) factor \f$ I \f$.
    */
    void factor(Factor t_factor) {this->m_factor = t_factor;}

    /**
    * Get the gravitational constant of the central body.
    * \return The gravitational constant of the central body.
    */
    Real mu() const {return this->m_mu;}

    /**
    * Set the gravitational constant of the central body.
    * \param[in] t_mu The gravitational constant of the central body.
    */
    void mu(Real t_mu) {this->m_mu = t_mu;}

    /**
    * Print the orbit information on a string.
    * \return The orbit information string.
    */
    std::string info() const {
      std::ostringstream os;
      os <<
        "Cartesian orbit parameters:" << std::endl <<
        this->m_cart.info() << std::endl <<
        "Keplerian orbit parameters:" << std::endl <<
        this->m_kepl.info() << std::endl <<
        "Equinoctial orbit parameters:" << std::endl <<
        this->m_equi.info() << std::endl <<
        "Quaternionic orbit parameters:" << std::endl <<
        this->m_quat.info() << std::endl <<
        "Orbital anomalies:" << std::endl <<
        this->m_anom.info() << std::endl <<
        "Orbit:" << std::endl <<
        "µ: grav. const. = " << this->m_mu << " (UA³/day²)" << std::endl <<
        "type            = " << (this->m_type == Type::UNDEFINED ? "undefined" :
          (this->m_type == Type::ELLIPTIC ? "elliptic" :
          (this->m_type == Type::PARABOLIC ? "parabolic" :
          (this->m_type == Type::HYPERBOLIC ? "hyperbolic" : "undefined")))) << std::endl <<
        "factor          = " << (this->m_factor == Factor::UNDEFINED ? "undefined" :
          (this->m_factor == Factor::POSIGRADE ? "posigrade" :
          (this->m_factor == Factor::RETROGRADE ? "retrograde" : "undefined"))) << std::endl;
        return os.str();
      }

      /**
      * Print the orbit information on a stream.
      * \param[in,out] os Output stream.
      */
      void info(std::ostream & os) {os << this->info();}

      /**
      * Reset the orbit parameters to undefined values.
      */
      void reset() {
        this->m_cart.reset();
        this->m_kepl.reset();
        this->m_equi.reset();
        this->m_quat.reset();
        this->m_anom.reset();
        this->m_type   = Type::UNDEFINED;
        this->m_factor = Factor::UNDEFINED;
      }

      /**
      * Compute the rotation matrix of the orbital plane through the keplerian orbital elements.
      * \return The rotation matrix of the orbital plane.
      */
      Rotation keplerian_to_reference() const
      {
        Real s_Omega{std::sin(this->m_kepl.Omega())};
        Real c_Omega{std::cos(this->m_kepl.Omega())};
        Real s_omega{std::sin(this->m_kepl.omega())};
        Real c_omega{std::cos(this->m_kepl.omega())};
        Real s_i{std::sin(this->m_kepl.i())};
        Real c_i{std::cos(this->m_kepl.i())};

        Rotation R;
        R <<
          c_Omega*c_omega - s_Omega*s_omega*c_i, -c_Omega*s_omega - s_Omega*c_omega*c_i, s_Omega*s_i,
          s_Omega*c_omega + c_Omega*s_omega*c_i, -s_Omega*s_omega + c_Omega*c_omega*c_i, -c_Omega*s_i,
          s_i*s_omega,                           s_i*c_omega,                            c_i;
        return R;
      }

      /**
      * Compute the rotation matrix of the orbital plane through the equinoctial orbital elements.
      * \return The rotation matrix of the orbital plane.
      */
      Rotation equinoctial_to_reference() const
      {
        Real h{this->m_equi.h()};
        Real k{this->m_equi.k()};
        Real h2{h*h};
        Real k2{k*k};
        Real hk{h*k};
        Real I{static_cast<Real>(this->m_factor)};

        Rotation R;
        R <<
          1.0-h2+k2, 2.0*I*hk,      2.0*h,
          2.0*hk,    I*(1.0+h2-k2), -2.0*k,
          -2.0*I*h,  2.0*k,         I*(1.0-h2-k2);
        return R / (1.0+h2+k2);
      }

      /**
      * Compute the Frenet-Serret frame of the orbit (radial, tangential, normal)
      * through the cartesian orbital elements.
      * \return The Frenet-Serret frame of the orbit.
      */
      Rotation cartesian_to_frenet_rtn() const
      {
        // Initialize the Frenet-Serret frame
        Rotation rtn;

        // Compute the radial vector
        rtn.col(0) = this->m_cart.r();
        rtn.col(0).normalize();

        // Compute the tangential vector
        rtn.col(1) = rtn.col(0).cross(this->m_cart.v());
        rtn.col(1).normalize();

        // Compute the binormal vector
        rtn.col(2) = rtn.col(0).cross(rtn.col(1));

        // Return the Frenet-Serret frame
        return rtn;
      }

      /**
      * Compute the Frenet-Serret frame of the orbit (radial, tangential, normal)
      * through the equinoctial orbital elements.
      * \return The Frenet-Serret frame of the orbit.
      */
      Rotation equinoctial_to_frenet_rtn() const
      {
        Real h{this->m_equi.h()};
        Real k{this->m_equi.k()};
        Real c_L{std::cos(this->m_anom.L())};
        Real s_L{std::sin(this->m_anom.L())};
        Real h2{h*h};
        Real k2{k*k};
        Real hk{h*k};
        Real bf{1.0+h2+k2};
        Real I{static_cast<Real>(this->m_factor)};

        Rotation R;
        R <<
          ((h2-k2+1)*c_L+2.0*I*hk*s_L)/bf,   ((k2-h2-1.0)*s_L+2.0*hk*I*c_L)/bf, 2.0*k/bf,
          (2.0*hk*c_L-I*(h2-k2-1.0)*s_L)/bf, (I*(1.0+k2-h2)*c_L-2.0*hk*s_L)/bf, -2.0*h/bf,
          (2.0*(h*s_L-c_L*k*I))/bf,          (2.0*k*I*s_L+2.0*h*c_L)/bf,        I*(1-h2-k2)/bf;
        return R;
      }

      /**
      * Transform a vector in the Frent-Serret frame of the orbit (radial, tangential, normal) to
      * a vector in the cartesian frame through the cartesian orbital elements.
      * \param[in] vec The vector in Frenet-Serret frame.
      * \return The vector in the cartesian frame.
      */
      Vector3 cartesian_rtn_to_xyz(Vector3 const & vec) const
      {
        return this->cartesian_to_frenet_rtn() * vec;
      }

      /**
      * Transform a vector in the Frent-Serret frame of the orbit (radial, tangential, normal) to
      * a vector in the cartesian frame through the equinoctial orbital elements.
      * \param[in] vec The vector in Frenet-Serret frame.
      * \return The vector in the cartesian frame.
      */
      Vector3 equinoctial_rtn_to_xyz(Vector3 const & vec) const
      {
        return this->equinoctial_to_frenet_rtn() * vec;
      }

      /**
      * Transform a vector in the cartesian frame to a vector in the Frent-Serret frame of the orbit
      * (radial, tangential, normal) through the cartesian orbital elements.
      * \param[in] vec The vector in the cartesian frame.
      * \return The vector in the Frenet-Serret frame.
      */
      Vector3 cartesian_xyz_to_rtn(Vector3 const & vec) const
      {
        return this->cartesian_to_frenet_rtn().transpose() * vec;
      }

      /**
      * Transform a vector in the cartesian frame to a vector in the Frent-Serret frame of the orbit
      * (radial, tangential, normal) through the equinoctial orbital elements.
      * \param[in] vec The vector in the cartesian frame.
      * \return The vector in the Frenet-Serret frame.
      */
      Vector3 equinoctial_xyz_to_rtn(Vector3 const & vec) const
      {
        return this->equinoctial_to_frenet_rtn().transpose() * vec;
      }

  }; // class Orbit

} // namespace Astro

#endif // ASTRO_ORBIT_HXX
