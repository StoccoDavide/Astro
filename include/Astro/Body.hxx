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

#ifndef ASTRO_BODY_HXX
#define ASTRO_BODY_HXX

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

    Real m_mass{QUIET_NAN}; /**< Mass of the astronomical body. */

  public:
    /**
    * Class constructor for the astronomical body object.
    */
    Body() {}

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
    * Get the mass of the astronomical body.
    * \return The mass of the astronomical body.
    */
    Real mass() const {return this->m_mass;}

    /**
    * Set the mass of the astronomical body.
    * \param[in] t_mass The mass of
    */
    void mass(Real t_mass) {this->m_mass = t_mass;}

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
      Real p{this->m_equi.p()};
      Real f{this->m_equi.f()};
      Real g{this->m_equi.g()};
      Real h{this->m_equi.h()};
      Real k{this->m_equi.k()};
      Real s_L{std::sin(this->m_anom.L())};
      Real c_L{std::cos(this->m_anom.L())};

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
      Real p{this->m_equi.p()};
      Real f{this->m_equi.f()};
      Real g{this->m_equi.g()};
      Real w{1.0 + (f*std::cos(this->m_anom.L()) + g*std::sin(this->m_anom.L()))};
      return VectorB(0.0, 0.0, 0.0, 0.0, 0.0, std::sqrt(p*this->m_mu)*power2(w/p));
    }

    /**
    * Compute the vector of the first-order cartesian equations of orbital motion given the
    * thrust vector components \f$ \mathbf{t} = [t_{\text{rad}}, t_{\text{tan}}, t_{\text{nor}}]^\top \f$.
    * \param[in] thrust_rtn The thrust components vector.
    * \return The vector of the first-order cartesian equations of orbital motion.
    */
    Vector6 cartesian_eom(Vector3 const & thrust_rtn) const
    {
      // Compute the cartesian perturbation
      Vector3 thrust_xyz(this->cartesian_rtn_to_xyz(thrust_rtn));

      // Compute the equations of motion
      Vector3 r(this->m_cart.r());
      Vector3 v(this->m_cart.v());
      Real r_norm{r.norm()};
      Real mur3{this->m_mu/power3(r_norm)};
      return Vector6(
        v.x(),
        v.y(),
        v.z(),
        thrust_xyz.x()/this->m_mass - mur3*r.x(),
        thrust_xyz.y()/this->m_mass - mur3*r.y(),
        thrust_xyz.z()/this->m_mass - mur3*r.z()
      );
    }

    /**
    * Compute the vector of the first-order cartesian equations of orbital motion given the
    * thrust vector components \f$ \mathbf{t} = [t_{\text{rad}}, t_{\text{tan}}, t_{\text{nor}}]^\top \f$.
    * \param[in] thrust_rad The radial thrust component.
    * \param[in] thrust_tan The tangential thrust component.
    * \param[in] thrust_nor The normal thrust component.
    * \return The vector of the first-order cartesian equations of orbital motion.
    */
    Vector6 cartesian_eom(Real thrust_rad, Real thrust_tan, Real thrust_nor) const
    {
      return this->cartesian_eom(Vector3(thrust_rad, thrust_tan, thrust_nor));
    }

    /**
    * Compute the system of first-order modified equinoctial equations of orbital motion given the
    * thrust vector components \f$ \mathbf{t} = [t_{\text{rad}}, t_{\text{tan}}, t_{\text{nor}}]^\top \f$.
    * \param[in] thrust The components (radial, tangential, normal) of the thrust.
    * \return The system of first-order modified equinoctial equations of orbital motion.
    */
    Vector6 equinoctial_eom(Vector3 const & thrust) const
    {
      Real p{this->m_equi.p()};
      Real f{this->m_equi.f()};
      Real g{this->m_equi.g()};
      Real h{this->m_equi.h()};
      Real k{this->m_equi.k()};
      Real L{this->m_anom.L()};

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

      res(5) += std::sqrt(p*this->m_mu)*power2(w/p);

      return res;
    }

    /**
    * Compute the system of first-order modified equinoctial equations of orbital motion given the
    * thrust vector components \f$ \mathbf{t} = [t_{\text{rad}}, t_{\text{tan}}, t_{\text{nor}}]^\top \f$.
    * \param[in] thrust_rad The radial component of the thrust.
    * \param[in] thrust_tan The tangential component of the thrust.
    * \param[in] thrust_nor The normal component of the thrust.
    * \return The system of first-order modified equinoctial equations of orbital motion.
    */
    Vector6 equinoctial_eom(Real thrust_rad, Real thrust_tan, Real thrust_nor) const
    {
      return this->equinoctial_eom(Vector3(thrust_rad, thrust_tan, thrust_nor));
    }

  }; // class Orbit

} // namespace Astro

#endif // ASTRO_ORBIT_HXX
