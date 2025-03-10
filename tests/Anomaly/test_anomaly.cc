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
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators_range.hpp>

using namespace Catch::Matchers;
using namespace Astro;
using namespace Astro::OrbitalElements;

TEST_CASE("Anomaly conversions") {

  // Create a Keplerian orbital elements object
  Real nu_tmp, M_tmp, E_tmp, L_tmp, lambda_tmp;

  auto nu = GENERATE(range(0.1, 1.2, 0.1));
  auto a = GENERATE(range(0.1, 1.2, 0.1));
  auto e = GENERATE(range(0.2, 0.3, 0.1));
  auto i = GENERATE(range(0.3, 0.4, 0.1));
  auto Omega = GENERATE(range(0.4, 0.5, 0.1));
  auto omega = GENERATE(range(0.5, 0.6, 0.1));
  auto I = GENERATE(Factor::POSIGRADE, Factor::RETROGRADE);

  // Set the Keplerian orbital elements
  Anomaly anom;
  Keplerian kepl(a, e, i, Omega, omega);
  kepl.sanity_check();

  anom.set_nu(nu, kepl, I);

  REQUIRE_THAT(anom.nu, WithinRel(nu, EPSILON_HIGH));
  nu_tmp     = anom.nu;
  M_tmp      = anom.M;
  E_tmp      = anom.E;
  L_tmp      = anom.L;
  lambda_tmp = anom.lambda;

  anom.set_M(anom.M, kepl, I);

  REQUIRE_THAT(anom.nu, WithinRel(nu_tmp, EPSILON_HIGH));
  REQUIRE_THAT(anom.M, WithinRel(M_tmp, EPSILON_HIGH));
  REQUIRE_THAT(anom.E, WithinRel(E_tmp, EPSILON_HIGH));
  REQUIRE_THAT(anom.L, WithinRel(L_tmp, EPSILON_HIGH));
  REQUIRE_THAT(anom.lambda, WithinRel(lambda_tmp, EPSILON_HIGH));
  nu_tmp     = anom.nu;
  M_tmp      = anom.M;
  E_tmp      = anom.E;
  L_tmp      = anom.L;
  lambda_tmp = anom.lambda;

  anom.set_E(anom.E, kepl, I);

  REQUIRE_THAT(anom.nu, WithinRel(nu_tmp, EPSILON_HIGH));
  REQUIRE_THAT(anom.M, WithinRel(M_tmp, EPSILON_HIGH));
  REQUIRE_THAT(anom.E, WithinRel(E_tmp, EPSILON_HIGH));
  REQUIRE_THAT(anom.L, WithinRel(L_tmp, EPSILON_HIGH));
  REQUIRE_THAT(anom.lambda, WithinRel(lambda_tmp, EPSILON_HIGH));
  nu_tmp     = anom.nu;
  M_tmp      = anom.M;
  E_tmp      = anom.E;
  L_tmp      = anom.L;
  lambda_tmp = anom.lambda;

  anom.set_L(anom.L, kepl, I);

  REQUIRE_THAT(anom.nu, WithinRel(nu_tmp, EPSILON_HIGH));
  REQUIRE_THAT(anom.M, WithinRel(M_tmp, EPSILON_HIGH));
  REQUIRE_THAT(anom.E, WithinRel(E_tmp, EPSILON_HIGH));
  REQUIRE_THAT(anom.L, WithinRel(L_tmp, EPSILON_HIGH));
  REQUIRE_THAT(anom.lambda, WithinRel(lambda_tmp, EPSILON_HIGH));
  nu_tmp     = anom.nu;
  M_tmp      = anom.M;
  E_tmp      = anom.E;
  L_tmp      = anom.L;
  lambda_tmp = anom.lambda;

  anom.set_lambda(anom.lambda, kepl, I);

  REQUIRE_THAT(anom.nu, WithinRel(nu_tmp, EPSILON_HIGH));
  REQUIRE_THAT(anom.M, WithinRel(M_tmp, EPSILON_HIGH));
  REQUIRE_THAT(anom.E, WithinRel(E_tmp, EPSILON_HIGH));
  REQUIRE_THAT(anom.L, WithinRel(L_tmp, EPSILON_HIGH));
  REQUIRE_THAT(anom.lambda, WithinRel(lambda_tmp, EPSILON_HIGH));

}
