/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Astro project is distributed under the GNU GPLv3.                                         *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "Astro.hh"

using namespace Astro;
using namespace Astro::OrbitalElements;

int main()
{


  auto a{1.0};
  auto e{0.1};
  auto i{0.1};
  auto Omega{0.1};
  auto omega{0.1};
  auto I{Factor::POSIGRADE};
  auto mu{1.0};
  Keplerian kepl(a, e, i, Omega, omega);
  std::cout << kepl.info() << std::endl;

  Anomaly anom;
  anom.set_nu(0.1, kepl, I);

  Equinoctial equi;
  Cartesian cart;
  Quaternionic quat;

  keplerian_to_equinoctial(kepl, I, equi);
  std::cout << equi.info() << std::endl;

  equinoctial_to_keplerian(equi, kepl);
  std::cout << kepl.info() << std::endl;

  equinoctial_to_cartesian(equi, anom, I, mu, cart);
  std::cout << cart.info() << std::endl;

  cartesian_to_keplerian(cart, mu, kepl);
  std::cout << kepl.info() << std::endl;


  return 0;
}
