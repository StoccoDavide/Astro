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
  class Body : public Orbit
  {
    using MatrixA = Eigen::Matrix<Real, 6, 3>; /**< Matrix \f$ \mathbf{A} \in \mathbb{R}^{6 \times 3}\f$. */
    using VectorB = Eigen::Vector<Real, 6>;    /**< Vector \f$ \mathbf{b} \in \mathbb{R}^{6 \times 1}\f$. */

    std::string m_name{"(undefined)"}; /**< Name of the astronomical body. */
    Real        m_radius{QUIET_NAN}; /**< Occupancy radius of the astronomical body (in AU). */
    Real        m_mass{QUIET_NAN}; /**< Mass of the astronomical body. */

  public:
    /**
    * Class constructor for the astronomical body object.
    */
    Body() {}

    /**
    * Class constructor for the astronomical body object with a given mass.
    * \param[in] t_mass The mass of the astronomical body.
    */
    Body(std::string const & t_name, Real t_mass, Real t_radius = 0.0)
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
    * Get the mass of the astronomical body.
    * \return The mass of the astronomical body.
    */
    Real mass() const {return this->m_mass;}

    /**
    * Set the mass of the astronomical body.
    * \param[in] t_mass The mass of
    */
    void mass(Real t_mass) {
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
    void radius(Real t_radius) {
      #define CMD "Astro::Body::radius(...): "

      ASTRO_ASSERT(t_radius >= 0.0, CMD "radius must be non-negative.");
      this->m_radius = t_radius;

      #undef CMD
    }

    /**
    * Compute the vector of the first-order cartesian equations of orbital motion given the
    * thrust vector components \f$ \mathbf{t} = [t_{\text{rad}}, t_{\text{tan}}, t_{\text{nor}}]^\top \f$.
    * \param[in] x The cartesian state vector \f$ \mathbf{x} = [r_x, r_y, r_z, v_x, v_y, v_z]^\top \f$.
    * \param[in] thrust_rtn The thrust components vector.
    * \return The vector of the first-order cartesian equations of orbital motion.
    */
    Vector6 cartesian_eom(Vector6 const & x, Vector3 const & thrust_rtn) const
    {
      // Set the new cartesian state
      Vector3 r(x.head<3>());
      Vector3 v(x.tail<3>());

      // Compute the cartesian perturbation
      Vector3 thrust_xyz(this->cartesian_rtn_to_xyz(thrust_rtn));
      thrust_xyz /= this->m_mass;

      // Compute the equations of motion
      Real r_norm{r.norm()};
      Vector3 r_mur3(this->m_mu/Power3(r_norm) * r);
      return Vector6(
        v.x(),
        v.y(),
        v.z(),
        thrust_xyz.x() - r_mur3.x(),
        thrust_xyz.y() - r_mur3.y(),
        thrust_xyz.z() - r_mur3.z()
      );
    }

    /**
    * Compute the vector of the first-order cartesian equations of orbital motion given the
    * thrust vector components \f$ \mathbf{t} = [t_{\text{rad}}, t_{\text{tan}}, t_{\text{nor}}]^\top \f$.
    * \param[in] x The cartesian state vector \f$ \mathbf{x} = [r_x, r_y, r_z, v_x, v_y, v_z]^\top \f$.
    * \param[in] thrust_rad The radial thrust component.
    * \param[in] thrust_tan The tangential thrust component.
    * \param[in] thrust_nor The normal thrust component.
    * \return The vector of the first-order cartesian equations of orbital motion.
    */
    Vector6 cartesian_eom(Vector6 const & x, Real thrust_rad, Real thrust_tan, Real thrust_nor) const
    {
      return this->cartesian_eom(x, Vector3(thrust_rad, thrust_tan, thrust_nor));
    }

    /**
    * Compute the matrix of the first-order modified equinoctial equations of orbital motion as
    * \f[ A = \begin{bmatrix}
    *   0 & 2\sqrt{\displaystyle\frac{p}{\mu}} & 0 \\
    *   \sqrt{p}\sin(L) & \sqrt{\displaystyle\frac{p}{\mu}}\left[(w+1)\cos(L)+f\right] & -\sqrt{\displaystyle\frac{p}{\mu}}(h\sin(L)-k\cos(L)) \\
    *   -\sqrt{p}\cos(L) & \sqrt{\displaystyle\frac{p}{\mu}}\left[(w+1)\sin(L)+g\right] & \sqrt{\displaystyle\frac{p}{\mu}}(h\cos(L)+k\sin(L)) \\
    *   0 & 0 & \displaystyle\frac{\sqrt{p}}{2}\displaystyle\frac{s^2\cos(L)}{2} \\
    *   0 & 0 & \displaystyle\frac{\sqrt{p}}{2}\displaystyle\frac{s^2\sin(L)}{2} \\
    *   0 & 0 & \displaystyle\frac{\sqrt{p}}{2}(h\sin(L)-k\cos(L))
    * \end{bmatrix} \text{.} \f]
    * \return The matrix \f$ \mathbf{A} \f$.
    */
    MatrixA equinoctial_eom_A() const
    {
      Real const & p{this->m_equi.p};
      Real const & f{this->m_equi.f};
      Real const & g{this->m_equi.g};
      Real const & h{this->m_equi.h};
      Real const & k{this->m_equi.k};
      Real s_L{std::sin(this->m_anom.L)};
      Real c_L{std::cos(this->m_anom.L)};

      Real w{1.0 + (f*c_L+g*s_L)};
      Real s2{1.0 + (h*h+k*k)};
      Real bf{std::sqrt(p/this->m_mu)};
      Real bf1{bf/w};
      Real bf2{(h*s_L-k*c_L)*bf1};

      MatrixA A;
      A <<
        0.0,     bf1*2.0*p,           0.0,
        bf*s_L,  bf1*((w+1.0)*c_L+f), -bf2*g,
        -bf*c_L, bf1*((w+1.0)*s_L+g), bf2*f,
        0.0,     0.0,                 bf1*s2*c_L/2.0,
        0.0,     0.0,                 bf1*s2*s_L/2.0,
        0.0,     0.0,                 bf2;
      return A;
    }

    /**
    * Compute the vector of the first-order modified equinoctial equations of orbital motion as
    * \f[ \mathbf{b} = \left[0, 0, 0, 0, 0, \sqrt{\mu p} \left(\displaystyle\frac{w}{p}\right)^2
    * \right]^\top \text{.} \f]
    * \return The vector \f$ \mathbf{b} \f$.
    */
    VectorB equinoctial_eom_b() const
    {
      Real const & p{this->m_equi.p};
      Real const & f{this->m_equi.f};
      Real const & g{this->m_equi.g};
      Real w{1.0 + (f*std::cos(this->m_anom.L) + g*std::sin(this->m_anom.L))};
      return VectorB(0.0, 0.0, 0.0, 0.0, 0.0, std::sqrt(p*this->m_mu)*Power2(w/p));
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
      Real const & p{x(0)};
      Real const & f{x(1)};
      Real const & g{x(2)};
      Real const & h{x(3)};
      Real const & k{x(4)};
      Real const & L{x(5)};

      Real s_L{std::sin(L)};
      Real c_L{std::cos(L)};

      Real w{1.0 + (f*c_L + g*s_L)};
      Real s2{1.0 + (h*h + k*k)};
      Real r{std::sqrt(p/this->m_mu)};

      Real C_rad{thrust(0)*r/this->m_mass};
      Real C_tan{thrust(1)*r/this->m_mass/w};
      Real C_nor{thrust(2)*r/this->m_mass/w};

      Vector6 res;
      res <<
        /* dp/dt */ 0.0*C_rad*s2,
        /* df/dt */ 0.0*C_tan,
        /* dg/dt */ 0.0*C_nor,
        /* dh/dt */ 0.0,
        /* dk/dt */ 0.0,
        /* dL/dt */ std::sqrt(p*this->m_mu) * Power2(w/p);
      return res;

      //res(1) += C_rad*s_L;
      //res(2) -= C_rad*c_L;
//
      //res(0) += C_tan*2.0*p;
      //res(1) += C_tan*((w+1)*c_L+f);
      //res(2) += C_tan*((w+1)*s_L+g);
//
      //res(1) -= C_nor*(h*s_L-k*c_L)*g;
      //res(2) += C_nor*(h*s_L-k*c_L)*f;
      //res(3) += C_nor*s2*c_L/2.0;
      //res(4) += C_nor*s2*s_L/2.0;
      //res(5) += C_nor*(h*s_L-k*c_L);
//
      //res(5) += std::sqrt(this->m_mu / Power3(p)) * Power2(w);

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
    template < Coordinates IntCoords, Coordinates OutCoords = Coordinates::CARTESIAN, Integer S>
    bool integrate(Sandals::RungeKutta<Real, S, 6, 0> & rk, VectorX const &t_mesh, Vector6 const &ics,
      Sandals::Solution<Real, 6, 0> & sol)
    {
      #define CMD "Astro::Body::integrate(...): "

      // Create the solution object according to the integration coordinates
      if constexpr (IntCoords == Coordinates::CARTESIAN) {
        rk.explicit_system(
          [this](Vector6 const & x, Real) -> Vector6 {return this->cartesian_eom(x, ZEROS_VEC3);}, // f(x, t)
          [](Vector6 const &, Real) -> Matrix6 {return ZEROS_MAT6;} // Jf_x(x, t)
        );
      } else if constexpr (IntCoords == Coordinates::KEPLERIAN) {
        ASTRO_ERROR(CMD "keplerian integration not implemented yet.");
      } else if constexpr (IntCoords == Coordinates::EQUINOCTIAL) {
        rk.explicit_system(
          [this](Vector6 const & x, Real) -> Vector6 {return this->equinoctial_eom(x, ZEROS_VEC3);}, // f(x, t)
          [](Vector6 const &, Real) -> Matrix6 {return ZEROS_MAT6;} // Jf_x(x, t)
        );
      } else if constexpr (IntCoords == Coordinates::QUATERNIONIC) {
        ASTRO_ERROR(CMD "quaternionic integration not implemented yet.");
      } else {
        ASTRO_ERROR(CMD "unknown coordinate system for the integration.");
      }

      // Transform the internal state to the output coordinate system
      rk.step_callback([this, &sol](Integer const i, Vector6 const & x, Real const) {

        (void)sol; // Avoid unused variable warning

        // Set the internal state according to the output coordinate system
        if constexpr (IntCoords == Coordinates::CARTESIAN) {
          this->set_cartesian(x);
        } else if constexpr (IntCoords == Coordinates::KEPLERIAN) {
          ASTRO_ERROR(CMD "keplerian output not implemented yet.");
        } else if constexpr (IntCoords == Coordinates::EQUINOCTIAL) {
          this->set_equinoctial(x.head<5>());
          this->set_anomaly().set_L(x(5), this->keplerian(), this->factor());
        } else if constexpr (IntCoords == Coordinates::QUATERNIONIC) {
          ASTRO_ERROR(CMD "quaternionic output not implemented yet.");
        } else {
          ASTRO_ERROR(CMD "unknown coordinate system for the output.");
        }

        // Set the output state in the solution object
        if constexpr (IntCoords != OutCoords) {
          if constexpr (OutCoords == Coordinates::CARTESIAN) {
            sol.x.col(i) << this->cartesian().vector();
          } else if constexpr (OutCoords == Coordinates::KEPLERIAN) {
            sol.x.col(i).head<5>() = this->keplerian().vector();
            sol.x(5,i) = this->anomaly().nu;
          } else if constexpr (OutCoords == Coordinates::EQUINOCTIAL) {
            sol.x.col(i).head<5>() = this->equinoctial().vector();
            sol.x(5,i) = this->anomaly().L;
          } else if constexpr (OutCoords == Coordinates::QUATERNIONIC) {
            ASTRO_ERROR(CMD "quaternionic output not implemented yet.");
          } else {
            ASTRO_ERROR(CMD "unknown coordinate system for the output.");
          }
        }

      });

      // Solve the equations of motion
      return rk.solve(t_mesh, ics, sol);

      #undef CMD
    }

  }; // class Body

} // namespace Astro

#endif // ASTRO_BODY_HH
