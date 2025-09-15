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

#ifndef ASTRO_BODY_HH
#define ASTRO_BODY_HH

#include "Astro/Orbit.hh"

// Sandals library
#include <Sandals/RungeKutta.hh>
#include <Sandals/Solution.hh>

namespace Astro
{

  /** State format for the orbit representation. */
  using Coordinates = enum class Coordinates : Integer {CARTESIAN = 0, KEPLERIAN = 1, EQUINOCTIAL = 2, QUATERNIONIC = 3};

  /*\
   |   ____            _
   |  | __ )  ___   __| |_   _
   |  |  _ \ / _ \ / _` | | | |
   |  | |_) | (_) | (_| | |_| |
   |  |____/ \___/ \__,_|\__, |
   |                     |___/
  \*/

  /**
  * \brief Astronomical body class container.
  *
  * The astronomical body class container is used to model the properties of an astronomical body
  * in space. The class provides methods for setting and getting the mass of the body, as well as
  * its orbit parameters.
  */
  class Body
  {
  public:
    using Anomaly = OrbitalElements::Anomaly;

  private:
    std::string m_name{"(undefined)"}; /**< Name of the astronomical body. */
    Real        m_epoch{0.0};          /**< Epoch of the orbital elements (in days). */
    Anomaly     m_epoch_anom;         /**< Orbital anomalies at the epoch. */
    Orbit       m_orbit;               /**< Orbit parameters of the astronomical body. */
    Real        m_radius{QUIET_NAN};   /**< Occupancy radius of the astronomical body (in AU). */
    Real        m_mass{QUIET_NAN};     /**< Mass of the astronomical body. */

  public:
    /**
    * Class constructor for the astronomical body object.
    */
    Body() {}

    /**
    * Class constructor for the astronomical body object with a given mass.
    * \param[in] t_mass The mass of the astronomical body.
    */
    Body(std::string const & t_name, Real const t_mass, Real const t_radius = 0.0)
    {
      this->name(t_name);
      this->mass(t_mass);
      this->radius(t_radius);
    }

    /**
    * Enable the default astronomical body copy constructor.
    */
    Body(Body const &) = default;

    /**
    * Enable the default astronomical body move constructor.
    */
    Body(Body &&) = default;

    /**
    * Enable the default astronomical body assignment operator.
    */
    Body & operator=(const Body &) = default;

    /**
    * Enable the default astronomical body move assignment operator.
    */
    Body & operator=(Body &&) = default;

    /**
    * Get the name of the astronomical body.
    * \return The name of the astronomical body.
    */
    std::string name() const {return this->m_name;}

    /**
    * Set the name of the astronomical body.
    * \param[in] t_name The name of the astronomical body.
    */
    void name(std::string const & t_name)
    {
      #define CMD "Astro::Body::name(...): "

      ASTRO_ASSERT(!t_name.empty(), CMD "name cannot be empty.");
      this->m_name = t_name;

      #undef CMD
    }

    /**
    * Get the orbit parameters of the astronomical body.
    * \return The orbit parameters of the astronomical body.
    */
    Orbit const & orbit() const {return this->m_orbit;}

    /**
    * Set the orbit parameters of the astronomical body.
    * \return The orbit parameters of the astronomical body.
    */
    Orbit & set_orbit() {return this->m_orbit;}

    /**
    * Get the orbital anomalies at the epoch.
    * \return The orbital anomalies at the epoch.
    */
    Anomaly const & epoch_anomaly() const {return this->m_epoch_anom;}

    /**
    * Set the orbital anomalies at the epoch.
    * \return The orbital anomalies at the epoch.
    */
    Anomaly & set_epoch_anomaly() {return this->m_epoch_anom;}

    /**
    * Set the orbital anomalies at the epoch.
    * \param[in] t_epoch_anomaly The orbital anomalies at the epoch.
    */
    void set_epoch_anomaly(Anomaly const & t_epoch_anomaly)
    {
      ASTRO_ASSERT(t_epoch_anomaly.sanity_check(),
        "Astro::Orbit::epoch_anomaly(...): invalid orbital anomalies detected.");
      this->m_epoch_anom = t_epoch_anomaly;
    }

    /**
    * Compute the anomalies at a given time.
    * \param[in] t The time \f$ t \f$.
    * \return The anomalies at time \f$ t \f$.
    */
    Anomaly anomaly(Real const t) const
    {
      Anomaly anom;
      anom.set_lambda(
        this->m_orbit.compute_lambda(this->m_epoch+t, this->m_epoch_anom), this->m_orbit.keplerian(), this->m_orbit.factor()
      );
      return anom;
    }

    /**
    * Get the epoch of the orbital elements (in days).
    * \return The epoch of the orbital elements (in days).
    */
    Real epoch() const {return this->m_epoch;}

    /**
    * Set the epoch of the orbital elements (in days).
    * \param[in] t_epoch The epoch of the orbital elements (in days).
    */
    void set_epoch(Real const t_epoch) {this->m_epoch = t_epoch;}

    /**
    * Get the mass of the astronomical body.
    * \return The mass of the astronomical body.
    */
    Real mass() const {return this->m_mass;}

    /**
    * Set the mass of the astronomical body.
    * \param[in] t_mass The mass of
    */
    void mass(Real const t_mass) {
      #define CMD "Astro::Body::mass(...): "

      ASTRO_ASSERT(t_mass > 0.0, CMD "mass must be positive.");
      this->m_mass = t_mass;

      #undef CMD
    }

    /**
    * Get the radius of the astronomical body.
    * \return The radius of the astronomical body.
    */
    Real radius() const {return this->m_radius;}

    /**
    * Set the radius of the astronomical body.
    * \param[in] t_radius The radius of the astronomical body.
    */
    void radius(Real const t_radius) {
      #define CMD "Astro::Body::radius(...): "

      ASTRO_ASSERT(t_radius >= 0.0, CMD "radius must be non-negative.");
      this->m_radius = t_radius;

      #undef CMD
    }

    /**
    * Get the cartesian state vector of the astronomical body.
    * \return The cartesian state vector of the astronomical body.
    */
    Vector6 cartesian_state() const
    {
      Vector6 cart(this->m_orbit.cartesian().vector());
      ASTRO_ASSERT(cart.allFinite(),
        "Astro::Body::cartesian(...): invalid cartesian state vector detected.");
      return cart;
    }

    /**
    * Get the cartesian state vector of the astronomical body.
    * \param[in] t The time \f$ t \f$.
    * \return The cartesian state vector of the astronomical body.
    */
    Vector6 cartesian_state(Real const t) const
    {
      // Compute the anomaly at time t
      Anomaly anom;
      anom.set_lambda(
        this->m_orbit.compute_lambda(this->m_epoch+t, this->m_epoch_anom), this->m_orbit.keplerian(), this->m_orbit.factor()
      );

      // Compute the cartesian state vector
      OrbitalElements::Cartesian cart;
      OrbitalElements::equinoctial_to_cartesian(this->m_orbit.equinoctial(), anom.L, this->m_orbit.mu(), cart);
      return cart.vector();
    }

    /**
    * Get the cartesian position vector of the astronomical body.
    * \param[in] t The time \f$ t \f$.
    * \return The cartesian position vector of the astronomical body.
    */
    Vector3 cartesian_position(Real const t) const {return this->cartesian_state(t).head<3>();}

    /**
    * Get the cartesian velocity vector of the astronomical body.
    * \param[in] t The time \f$ t \f$.
    * \return The cartesian velocity vector of the astronomical body.
    */
    Vector3 cartesian_velocity(Real const t) const {return this->cartesian_state(t).tail<3>();}

    /**
    * Set the cartesian state vector of the astronomical body.
    * \param[in] cart The cartesian state vector of the astronomical body.
    */
    void set_cartesian_state(Vector6 const & cart)
    {
      this->set_orbit().set_cartesian(cart);
    }

    /**
    * Get the keplerian state vector of the astronomical body at time \f$ t \f$.
    * \param[in] t The time \f$ t \f$.
    * \return The keplerian state vector of the astronomical body.
    */
    Vector6 keplerian_state(Real const t) const
    {
      Anomaly anom{this->anomaly(t)};
      Vector6 kepl;
      kepl <<
        this->orbit().keplerian().a,
        this->orbit().keplerian().e,
        this->orbit().keplerian().i,
        this->orbit().keplerian().Omega,
        this->orbit().keplerian().omega,
        anom.M;
      ASTRO_ASSERT(kepl.allFinite(),
        "Astro::Body::keplerian(...): invalid keplerian state vector detected.");
      return kepl;
    }

    /**
    * Set the keplerian state vector of the astronomical body.
    * \param[in] state The keplerian state vector of the astronomical body.
    */
    void set_keplerian_state(Vector6 const & state)
    {
      Real E{OrbitalElements::M_to_E(state(5), state(1))};
      Real nu{OrbitalElements::E_to_nu(E, state(1))};
      this->set_orbit().set_keplerian(state.head<5>(), nu);
    }

    /**
    * Get the equinoctial state vector of the astronomical body.
    * \param[in] t The time \f$ t \f$.
    * \return The equinoctial state vector of the astronomical body.
    */
    Vector6 equinoctial_state(Real const t) const
    {
      Anomaly anom{this->anomaly(t)};
      Vector6 equi;
      equi <<
        this->orbit().equinoctial().p,
        this->orbit().equinoctial().f,
        this->orbit().equinoctial().g,
        this->orbit().equinoctial().h,
        this->orbit().equinoctial().k,
        anom.L;
      ASTRO_ASSERT(equi.allFinite(),
        "Astro::Body::equinoctial(...): invalid equinoctial state vector detected.");
      return equi;
    }

    /**
    * Set the equinoctial state vector of the astronomical body.
    * \param[in] state The equinoctial state vector of the astronomical body.
    */
    void set_equinoctial_state(Vector6 const & state)
    {
      this->set_orbit().set_equinoctial(state.head<5>(), state(5));
    }

    /**
    * Compute the vector of the first-order cartesian equations of orbital motion given the
    * thrust vector components \f$ \mathbf{t} = [t_{\text{rad}}, t_{\text{tan}}, t_{\text{nor}}]^\top \f$.
    * \param[in] x The cartesian state vector \f$ \mathbf{x} = [r_x, r_y, r_z, v_x, v_y, v_z]^\top \f$.
    * \param[in] thrust_rtn The thrust components vector.
    * \return The vector of the first-order cartesian equations of orbital motion.
    */
    template<typename T>
    Eigen::Matrix<T, 6, 1> cartesian_eom(Eigen::Matrix<T, 6, 1> const & x, Eigen::Matrix<T, 3, 1> const & thrust_rtn) const
    {
      // Compute the cartesian perturbation
      Eigen::Matrix<T, 3, 1> thrust_xyz(this->orbit().cartesian_rtn_to_xyz(thrust_rtn));
      thrust_xyz /= this->m_mass;

      // Compute the equations of motion
      T mur3{this->orbit().mu()/Power3(x.template head<3>().norm())};
      return Eigen::Matrix<T, 6, 1>(
        /* dx/dt   */ x[3],
        /* dy/dt   */ x[4],
        /* dz/dt   */ x[5],
        /* dv_x/dt */ thrust_xyz[0] - mur3*x[0],
        /* dv_y/dt */ thrust_xyz[1] - mur3*x[1],
        /* dv_z/dt */ thrust_xyz[2] - mur3*x[2]
      );
    }

    /**
    * Compute the vector of the first-order keplerian equations of orbital motion given the
    * thrust vector components \f$ \mathbf{t} = [t_{\text{rad}}, t_{\text{tan}}, t_{\text{nor}}]^\top \f$.
    * \param[in] x The keplerian state vector \f$ \mathbf{x} = [a, e, i, \Omega, \omega, M]^\top \f$.
    * \param[in] thrust The components (radial, tangential, normal) of the thrust.
    * \return The vector of the first-order keplerian equations of orbital motion.
    */
    Vector6 keplerian_eom(Vector6 const & x, Vector3 const & thrust) const
    {
      // https://farside.ph.utexas.edu/teaching/celestial/Celestial/node164.html

      Real const & a{x[0]};
      Real const & e{x[1]};
      Real const & i{x[2]};
      // Real const & Omega{x[3]};
      Real const & omega{x[4]};
      Real const & M{x[5]};
      Real const & mu{this->orbit().mu()};

      Real const C_rad{thrust[0]/this->m_mass};
      Real const C_tan{thrust[1]/this->m_mass};
      Real const C_nor{thrust[2]/this->m_mass};

      // Mean motion
      Real const n{std::sqrt(mu / std::pow(a, 3))};
      Real const e2{e*e};
      Real const ecc{std::sqrt(1.0 - e2)};

      // Anomalies calculations
      Real const E{OrbitalElements::M_to_E(M, e)};
      Real const cos_E{std::cos(E)};
      Real const nu{OrbitalElements::E_to_nu(E, e)};
      Real const cos_nu{std::cos(nu)};
      Real const sin_nu{std::sin(nu)};

      // Other parameters
      Real const p{a*(1.0 - e2)};
      Real const r{p / (1.0 + e*cos_nu)};
      Real const u{omega + nu};
      Real const h{ecc * std::sqrt(mu*a)}; // Specific angular momentum

      Real const sin_i{std::sin(i)};
      Real const cos_i{std::cos(i)};
      Real const sin_u{std::sin(u)};
      Real const cos_u{std::cos(u)};

      Vector6 res;
      res <<
        /* da/dt */ 2.0*a*h/(mu*(1.0 - e2))*(e*sin_nu*C_rad + (1.0 + e*cos_nu)*C_tan),
        /* de/dt */ h/mu*(sin_nu*C_rad + (cos_nu + cos_E)*C_tan),
        /* di/dt */ cos_u*r*C_nor/h,
        /* dO/dt */ sin_u*r*C_nor/(h*sin_i),
        /* do/dt */ -h/(mu*e)*(cos_nu*C_rad - (2.0 + e*cos_nu)/(1.0 + e*cos_nu)*sin_nu*C_tan) - cos_i*sin_u*r*C_nor/(h*sin_i),
        /* dM/dt */ n + h*ecc/(mu*e) * ((cos_nu - 2.0*e*r/p)*C_rad - (1.0 + r/p)*sin_nu*C_tan);

      return res;
    }

    /**
    * Compute the system of first-order modified equinoctial equations of orbital motion given the
    * thrust vector components \f$ \mathbf{t} = [t_{\text{rad}}, t_{\text{tan}}, t_{\text{nor}}]^\top \f$.
    * \param[in] x The modified equinoctial state vector \f$ \mathbf{x} = [p, f, g, h, k, L]^\top \f$.
    * \param[in] thrust The components (radial, tangential, normal) of the thrust.
    * \return The system of first-order modified equinoctial equations of orbital motion.
    */
    Vector6 equinoctial_eom(Vector6 const & x, Vector3 const & thrust) const
    {
      // Source: https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf

      Real const & p{x[0]};
      Real const & f{x[1]};
      Real const & g{x[2]};
      Real const & h{x[3]};
      Real const & k{x[4]};
      Real const & L{x[5]};
      Real const & mu{this->orbit().mu()};

      Real sin_L{std::sin(L)};
      Real cos_L{std::cos(L)};

      Real w{1.0 + (f*cos_L + g*sin_L)};
      Real s2{1.0 + (h*h + k*k)};
      Real sqrt_p_mu{std::sqrt(p/mu)};

      Real C_rad{thrust[0]/this->m_mass};
      Real C_tan{thrust[1]/this->m_mass};
      Real C_nor{thrust[2]/this->m_mass};

      Vector6 res;
      res <<
        /* dp/dt */ 2.0*p*sqrt_p_mu/w*C_tan,
        /* df/dt */ sqrt_p_mu*(+C_rad*sin_L + ((w + 1.0)*cos_L + f)*C_tan/w - (h*sin_L - k*cos_L)*g/w*C_nor),
        /* dg/dt */ sqrt_p_mu*(-C_rad*cos_L + ((w + 1.0)*sin_L + g)*C_tan/w + (h*sin_L - k*cos_L)*g/w*C_nor),
        /* dh/dt */ sqrt_p_mu*s2*C_nor*cos_L/(2.0*w),
        /* dk/dt */ sqrt_p_mu*s2*C_nor*sin_L/(2.0*w),
        /* dL/dt */ std::sqrt(mu*p) * Power2(w/p) + sqrt_p_mu/w*(h*sin_L - k*cos_L)*C_nor;

      return res;
    }

    /**
    * Compute the system of first-order modified equinoctial equations of orbital motion given the
    * thrust vector components \f$ \mathbf{t} = [t_{\text{rad}}, t_{\text{tan}}, t_{\text{nor}}]^\top \f$.
    * \param[in] x The modified equinoctial state vector \f$ \mathbf{x} = [p, f, g, h, k, L]^\top \f$.
    * \param[in] thrust_rad The radial component of the thrust.
    * \param[in] thrust_tan The tangential component of the thrust.
    * \param[in] thrust_nor The normal component of the thrust.
    * \return The system of first-order modified equinoctial equations of orbital motion.
    */
    Vector6 equinoctial_eom(Vector6 const & x, Real thrust_rad, Real thrust_tan, Real thrust_nor) const
    {
      return this->equinoctial_eom(x, Vector3(thrust_rad, thrust_tan, thrust_nor));
    }

    /**
    * Integrate the equations of motion using the Runge-Kutta method.
    * \param[in] rk The Runge-Kutta integrator.
    * \param[in] t_mesh The time mesh for the integration.
    * \param[in] ics The initial conditions for the integration (initial modified equinoctial state).
    * \param[out] sol The solution object to store the results of the integration.
    * \tparam S The stages of the Runge-Kutta method.
    * \tparam IntCoords The coordinate system for the integration.
    * \tparam OutCoords The coordinate system for the output of the integration.
    * \return The integrated modified equinoctial state vectors.
    * \note The `IntCoords` and `OutCoords` parameters are used to specify the coordinate systems for
    * the integration and output, respectively. The `IntCoords` parameter is used to specify the
    * coordinate system in which the equations of motion are integrated, while the `OutCoords` parameter
    * is used to specify the coordinate system in which the results of the integration are output.
    */
    template <Coordinates IntCoords, Coordinates OutCoords = Coordinates::CARTESIAN, Integer S>
    bool integrate(Sandals::RungeKutta<Real, S, 6, 0> & rk, VectorX const &t_mesh, Vector6 const &ics,
      Sandals::Solution<Real, 6, 0> & sol, bool const adaptive = false)
    {
      #define CMD "Astro::Body::integrate(...): "

      // Create the solution object according to the integration coordinates
      if constexpr (IntCoords == Coordinates::CARTESIAN) {
        rk.explicit_system(
          [this](Vector6 const & x, Real) -> Vector6 {return this->cartesian_eom<Real>(x, 50.0*Vector3::UnitZ());}, // f(x, t)
          [](Vector6 const &, Real) -> Matrix6 {return ZEROS_MAT6;} // Jf_x(x, t)
        );
      } else if constexpr (IntCoords == Coordinates::KEPLERIAN) {
        rk.explicit_system(
          [this](Vector6 const & x, Real) -> Vector6 {return this->keplerian_eom(x, 50.0*Vector3::UnitZ());}, // f(x, t)
          [](Vector6 const &, Real) -> Matrix6 {return ZEROS_MAT6;} // Jf_x(x, t)
        );
      } else if constexpr (IntCoords == Coordinates::EQUINOCTIAL) {
        rk.explicit_system(
          [this](Vector6 const & x, Real) -> Vector6 {return this->equinoctial_eom(x, 50.0*Vector3::UnitZ());}, // f(x, t)
          [](Vector6 const &, Real) -> Matrix6 {return ZEROS_MAT6;} // Jf_x(x, t)
        );
      } else if constexpr (IntCoords == Coordinates::QUATERNIONIC) {
        ASTRO_ERROR(CMD "quaternionic integration not implemented yet.");
      } else {
        ASTRO_ERROR(CMD "unknown coordinate system for the integration.");
      }

      // Transform the internal state to the output coordinate system
      rk.step_callback([this, &sol](Integer const i, Vector6 const & x, Real const t) {

        (void)sol; // Avoid unused variable warning

        // Set the internal state according to the output coordinate system
        if constexpr (IntCoords == Coordinates::CARTESIAN) {
          this->set_cartesian_state(x);
        } else if constexpr (IntCoords == Coordinates::KEPLERIAN) {
          this->set_keplerian_state(x);
        } else if constexpr (IntCoords == Coordinates::EQUINOCTIAL) {
          this->set_equinoctial_state(x);
        } else if constexpr (IntCoords == Coordinates::QUATERNIONIC) {
          ASTRO_ERROR(CMD "quaternionic output not implemented yet.");
        } else {
          ASTRO_ERROR(CMD "unknown coordinate system for the output.");
        }

        // Set the output state in the solution object
        if constexpr (IntCoords != OutCoords) {
          if constexpr (OutCoords == Coordinates::CARTESIAN) {
            sol.x.col(i) << this->cartesian_state(t);
          } else if constexpr (OutCoords == Coordinates::KEPLERIAN) {
            sol.x.col(i) << this->keplerian_state(t);
          } else if constexpr (OutCoords == Coordinates::EQUINOCTIAL) {
            sol.x.col(i) << this->equinoctial_state(t);
          } else if constexpr (OutCoords == Coordinates::QUATERNIONIC) {
            ASTRO_ERROR(CMD "quaternionic output not implemented yet.");
          } else {
            ASTRO_ERROR(CMD "unknown coordinate system for the output.");
          }
        }

      });

      // Solve the equations of motion
      if (adaptive) {
        return rk.adaptive_solve(t_mesh, ics, sol);
      } else {
        return rk.solve(t_mesh, ics, sol);
      }

      #undef CMD
    }

  }; // class Body

} // namespace Astro

#endif // ASTRO_BODY_HH
