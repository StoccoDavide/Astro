/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Astro project is distributed under the GNU GPLv3.                                         *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "Astro/Body.hh"
#include "Astro/Planets.hh"
#ifdef ASTRO_ENABLE_PLOTTING
#include "Astro/Plotting.hh"
#endif // ASTRO_ENABLE_PLOTTING

#include <Sandals/RungeKutta/RK4.hh>

using namespace Astro;

#ifdef ASTRO_ENABLE_PLOTTING
int main(int argc, char** argv) {
  TApplication app("app", &argc, argv);

  // Create canvas
  TCanvas* canvas = new TCanvas("c1", "Plotting test", 800, 600);
  canvas->SetGrid();

  // Set up 3D axis
  TAxis3D rulers; rulers.Draw();
  TAxis3D::ToggleRulers();
  TAxis3D::ToggleZoom();

  // Draw the absolute reference frame
  TObjArray* abs_axes = Astro::Plotting::DrawAbsoluteAxes(0.25);
  abs_axes->SetOwner(kTRUE); // Ensure axes are deleted with the canvas
#else
  int main() {
#endif // ASTRO_ENABLE_PLOTTING

  // Create a Body object (e.g., Earth)
  Body earth(Planets::Earth());

  std::cout << "Cartesian elements of the Earth:\n";
  std::cout << earth.orbit().cartesian().vector().transpose() << std::endl;
  std::cout << "Keplerian elements of the Earth:\n";
  std::cout << earth.orbit().keplerian().vector().transpose() << std::endl;
  std::cout << "Equinoctial elements of the Earth:\n";
  std::cout << earth.orbit().equinoctial().vector().transpose() << std::endl;

  // Set the time mesh for the integration
  Real t_start = 0.0;
  Real t_end = 365.25/2; // One year in days
  Real dt = 1.0; // Time step in days
  VectorX t_mesh = VectorX::LinSpaced((t_end-t_start)/dt + 1, t_start, t_end);

  // Set the initial state of the orbit (e.g., cartesian, keplerian, equinoctial)
  Vector6 ics;
  ics << earth.cartesian_state(t_start);
  std::cout << "Initial state vector (equinoctial coordinates):\n";
  std::cout << ics.transpose() << std::endl;

  // Integrate the orbit using the Runge-Kutta solver
  Sandals::RK4<Real, 6, 0> rk4;
  Sandals::Solution<Real, 6, 0> sol;
  earth.integrate<
    Astro::Coordinates::CARTESIAN, // Integration in equinoctial coordinates
    Astro::Coordinates::CARTESIAN    // Output in cartesian coordinates
    >(rk4, t_mesh, ics, sol);

#ifdef ASTRO_ENABLE_PLOTTING
  // Plot the orbit trace selecting last 3 rows (position) of the solution
  TPolyLine3D* orbit_trace = Plotting::DrawTrace(sol.x.topRows<3>());
  orbit_trace->Draw("same L");

  canvas->Update();
  app.Run();
#endif // ASTRO_ENABLE_PLOTTING

  return 0;
}
