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
    * Compute the vector of the first-order cartesian equations of orbital motion given the
    * thrust vector components \f$ \mathbf{t} = [t_{\text{rad}}, t_{\text{tan}}, t_{\text{nor}}]^\top \f$.
    * \param[in] x The cartesian state vector \f$ \mathbf{x} = [r_x, r_y, r_z, v_x, v_y, v_z]^\top \f$.
    * \param[in] thrust_rtn The thrust components vector.
    * \return The vector of the first-order cartesian equations of orbital motion.
    */
    Vector6 cartesian_eom(Vector6 const & x, Vector3 const & thrust_rtn) const
    {
      // Set the new cartesian state
      Vector3 r = x.head<3>();
      Vector3 v = x.tail<3>();

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

      Real w{1.0 + (f*c_L+g*s_L)};
      Real s2{1.0 + (h*h+k*k)};
      Real bf{std::sqrt(p/this->m_mu)};

      Real C_rad{thrust(0)*bf/this->m_mass};
      Real C_tan{thrust(1)*bf/this->m_mass/w};
      Real C_nor{thrust(2)*bf/this->m_mass/w};

      Vector6 res;
      res.setZero();

      res(1) += C_rad*s_L;
      res(2) -= C_rad*c_L;

      res(0) += C_tan*2.0*p;
      res(1) += C_tan*((w+1)*c_L+f);
      res(2) += C_tan*((w+1)*s_L+g);

      res(1) -= C_nor*(h*s_L-k*c_L)*g;
      res(2) += C_nor*(h*s_L-k*c_L)*f;
      res(3) += C_nor*s2*c_L/2.0;
      res(4) += C_nor*s2*s_L/2.0;
      res(5) += C_nor*(h*s_L-k*c_L);

      res(5) += std::sqrt(p*this->m_mu)*Power2(w/p);

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
    * Integrate the cartesian equations of motion using the Runge-Kutta method.
    * \param[in] rk The Runge-Kutta integrator.
    * \param[in] t_mesh The time mesh for the integration.
    * \param[in] ics The initial conditions for the integration (initial cartesian state).
    * \param[out] sol The solution object to store the results of the integration.
    * \tparam S The stages of the Runge-Kutta method.
    * \return The integrated cartesian state (position and velocity) vectors.
    */
    template <Integer S>
    bool integrate_cartesian(
      Sandals::RungeKutta<Real, S, 6, 0> & rk,
      VectorX const &t_mesh,
      Vector6 const &ics,
      Sandals::Solution<Real, 6, 0> & sol
    ) {
      sol.resize(t_mesh.size());
      rk.explicit_system(
        [this](Vector6 const & x, Real) -> Vector6 { // f(x, t)
          return this->cartesian_eom(x, ZEROS_VEC3);
        },
        [](Vector6 const &, Real) -> Matrix6 { // Jf_x(x, t)
          return ZEROS_MAT6;
        }
      );
      Vector6 x_step;
      Real h_new_step;
      sol.t(0)     = t_mesh(0);
      sol.x.col(0) = ics;
      for (Integer i = 0; i < t_mesh.size() - 1; ++i) {
        Real t_step{t_mesh(i)};
        Real h_step{t_mesh(i + 1) - t_step};
        if (!rk.advance(sol.x.col(i), t_step, h_step, x_step, h_new_step)) {
          return false; // Error in the integration
        } else {
          sol.x.col(i + 1) = x_step;
          sol.t(i + 1) = t_step + h_step;
          this->cartesian(x_step.head<3>(), x_step.tail<3>());
        }
      }
      return true;
    }

    /**
    * Integrate the modified equinoctial equations of motion using the Runge-Kutta method.
    * \param[in] rk The Runge-Kutta integrator.
    * \param[in] t_mesh The time mesh for the integration.
    * \param[in] ics The initial conditions for the integration (initial modified equinoctial state).
    * \param[out] sol The solution object to store the results of the integration.
    * \tparam S The stages of the Runge-Kutta method.
    * \return The integrated modified equinoctial state vectors.
    */
    template <Integer S>
    bool integrate_equinoctial(
      Sandals::RungeKutta<Real, S, 6, 0> & rk,
      VectorX const &t_mesh,
      Vector6 const &ics,
      Sandals::Solution<Real, 6, 0> & sol
    ) {
      Sandals::Solution<Real, 6, 0> sol_xyz; sol_xyz.resize(t_mesh.size());
      sol.resize(t_mesh.size());
      rk.explicit_system(
        [this](Vector6 const & x, Real) -> Vector6 { // f(x, t)
          return this->equinoctial_eom(x, ZEROS_VEC3);
        },
        [](Vector6 const &, Real) -> Matrix6 { // Jf_x(x, t)
          return ZEROS_MAT6;
        }
      );
      Vector6 x_step;
      Real h_new_step;
      sol.t(0)     = t_mesh(0);
      sol.x.col(0) = ics;
      for (Integer i = 0; i < t_mesh.size() - 1; ++i) {
        Real t_step{t_mesh(i)};
        Real h_step{t_mesh(i + 1) - t_step};
        if (!rk.advance(sol.x.col(i), t_step, h_step, x_step, h_new_step)) {
          return false; // Error in the integration
        } else {
          sol.x.col(i + 1) = x_step;
          sol.t(i + 1) = t_step + h_step;

          this->equinoctial(x_step(0), x_step(1), x_step(2), x_step(3), x_step(4));
          this->anomaly().set_L(x_step(5));

          sol_xyz.x.col(i + 1).head<3>() = this->cartesian().r;
          sol_xyz.x.col(i + 1).tail<3>() = this->cartesian().v;
          sol_xyz.t(i + 1) = t_step + h_step;
        }
      }
      sol = sol_xyz;
      return true;
    }

  }; // class Body

} // namespace Astro

#endif // ASTRO_BODY_HH
