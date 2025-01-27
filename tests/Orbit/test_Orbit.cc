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

int main()
{
  Astro::OrbitalElements::Cartesian cart;
  Astro::OrbitalElements::Keplerian kepl;
  Astro::OrbitalElements::Equinoctical equi;
  Astro::OrbitalElements::Quaternionic quat;

  std::cout << "Hello World!" << std::endl;

  return 0;
}
