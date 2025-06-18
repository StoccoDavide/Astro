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
#include "Astro/Plotting.hh"

#include <Sandals/RungeKutta/RK4.hh>

using namespace Astro;


int main(int argc, char** argv) {
  TApplication app("app", &argc, argv);

  // Create canvas
  TCanvas* canvas = new TCanvas("c1", "Plotting test", 800, 600);
  canvas->SetGrid();

  // Set up 3D axis
  TAxis3D rulers; rulers.Draw();
  TAxis3D::ToggleRulers();
  TAxis3D::ToggleZoom();

  // Add axes labels
  rulers.GetXaxis()->SetTitle("X (Km)");
  rulers.GetYaxis()->SetTitle("Y (Km)");
  rulers.GetZaxis()->SetTitle("Z (Km)");

  // Create the satellite
  Real payload_mass{100.0*1.339E-6}; // Mass of the satellite in Kg
  Real payload_radius{1E-3}; // Radius of the satellite in Km
  Body sat("sat", payload_mass, payload_radius);
  sat.set_factor(Factor::POSIGRADE);
  sat.set_mu(Planets::Earth_mu_KM3_DAY2);
  sat.set_keplerian(OrbitalElements::Keplerian(
    Planets::Earth_radius_KM + 4000.0, // Semi-major axis in AU
    1.36E-2, // Eccentricity
    Deg_To_Rad(9.015), // Inclination in radians
    0.0, // Longitude of ascending node in radians
    0.0  // Argument of periapsis in radians
  ));

  // Set the time mesh for the integration
  Real const t_start = 0.0; // days
  Real const t_end = 5.0; // days
  Real const dt = 5.0/(24*60); // days
  VectorX const t_mesh = VectorX::LinSpaced((t_end-t_start)/dt + 1, t_start, t_end);

  // Integrate the orbit using the Runge-Kutta solver
  Sandals::RK4<Real, 6, 0> rk;
  Sandals::Solution<Real, 6, 0> sol_cart, sol_equi, sol_kepl;
  Vector6 ics_cart, ics_equi, ics_kepl;
  ics_cart << sat.cartesian_state();
  ics_equi << sat.equinoctial_state();
  ics_kepl << sat.keplerian_state();

  sat.integrate<
    Astro::Coordinates::CARTESIAN, // Integration coordinates
    Astro::Coordinates::CARTESIAN // Output coordinates
    >(rk, t_mesh, ics_cart, sol_cart);

  // Plot the orbit trace selecting last 3 rows (position) of the solution
  TPolyLine3D* orbit_trace_cart = Plotting::DrawTrace(sol_cart.x.topRows<3>());
  orbit_trace_cart->SetLineColor(kRed);
  orbit_trace_cart->SetLineWidth(1);
  orbit_trace_cart->Draw("same L");

  sat.integrate<
    Astro::Coordinates::KEPLERIAN, // Integration coordinates
    Astro::Coordinates::CARTESIAN // Output coordinates
    >(rk, t_mesh, ics_kepl, sol_kepl);

  // Plot the orbit trace selecting last 3 rows (position) of the solution
  TPolyLine3D* orbit_trace_kepl = Plotting::DrawTrace(sol_kepl.x.topRows<3>());
  orbit_trace_kepl->SetLineColor(kGreen);
  orbit_trace_kepl->SetLineWidth(1);
  orbit_trace_kepl->Draw("same L");

  sat.integrate<
    Astro::Coordinates::EQUINOCTIAL, // Integration coordinates
    Astro::Coordinates::CARTESIAN // Output coordinates
    >(rk, t_mesh, ics_equi, sol_equi);

  // Plot the orbit trace selecting last 3 rows (position) of the solution
  TPolyLine3D* orbit_trace_equi = Plotting::DrawTrace(sol_equi.x.topRows<3>());
  orbit_trace_equi->SetLineColor(kBlue);
  orbit_trace_equi->SetLineWidth(1);
  orbit_trace_equi->Draw("same L");

  // Plot the wireframe sphere representing the earth
  TGeoVolume* earth_sphere = Plotting::DrawSphere(ZEROS_VEC3, Planets::Earth_radius_KM, kBlack, 1.0);
  earth_sphere->Draw("same");

  canvas->Update();

  // Make another canvas to plot the magnetic field with respect to time
  TCanvas* canvas2 = new TCanvas("c2", "Magnetic Field vs Time", 800, 600);
  canvas2->SetGrid();
  canvas2->SetTitle("Magnetic Field vs Time");

  // Plot the last row of the solution to check the orbit
  TGraph* graph_cart = new TGraph(t_mesh.size(), sol_cart.t.data(), sol_cart.x.row(3).data());
  graph_cart->SetLineColor(kRed);
  graph_cart->SetLineWidth(1);
  graph_cart->Draw("AL");

  TGraph* graph_kepl = new TGraph(t_mesh.size(), sol_kepl.t.data(), sol_kepl.x.row(3).data());
  graph_kepl->SetLineColor(kGreen);
  graph_kepl->SetLineWidth(1);
  graph_kepl->Draw("L SAME");

  TGraph* graph_equi = new TGraph(t_mesh.size(), sol_equi.t.data(), sol_equi.x.row(3).data());
  graph_equi->SetLineColor(kBlue);
  graph_equi->SetLineWidth(1);
  graph_equi->Draw("L SAME");

  //// Prepare data for magnetic field plot
  //std::vector<Real> time_data;
  //std::vector<Real> magnetic_field_data;
  //for (Integer i = 0; i < t_mesh.size(); ++i) {
  //  Vector3 position = sol.x.col(i).head<3>(); // Get position at time t_mesh[i]
  //  Vector3 magnetic_field = Planets::EarthMagneticFieldDipole(position); // Compute magnetic field
  //  time_data.push_back(t_mesh[i]);
  //  magnetic_field_data.push_back(magnetic_field.norm());
  //}
  //// Create TGraph for magnetic field
  //TGraph* graph = new TGraph(time_data.size(), time_data.data(), magnetic_field_data.data());
  //graph->SetLineColor(kBlue);
  //graph->SetLineWidth(2);
  //graph->Draw("AL");
//
  //graph->GetXaxis()->SetTitle("Time (days)");
  //graph->GetYaxis()->SetTitle("Magnetic Field (nT)");
  //graph->SetTitle("Magnetic Field Strength vs Time");

  canvas2->Update();
  app.Run();

  return 0;
}
