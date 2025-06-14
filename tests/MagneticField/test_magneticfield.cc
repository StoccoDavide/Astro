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
  Real payload_mass{100.0}; // Mass of the satellite in Kg
  Real payload_radius{1E-3}; // Radius of the satellite in Km
  Body sat("sat", payload_mass, payload_radius);
  sat.set_factor(Factor::POSIGRADE);
  sat.set_mu(Planets::Earth_mu_KM3_DAY2);
  std::cout << "Satellite mass: " << sat.mu() << " Kg" << std::endl;
  sat.set_keplerian(OrbitalElements::Keplerian(
    Planets::Earth_radius_KM + 400.0, // Semi-major axis in AU
    1.36E-4, // Eccentricity
    Deg_To_Rad(89.015), // Inclination in radians
    0.0, // Longitude of ascending node in radians
    0.0  // Argument of periapsis in radians
  ));

  sat.info(); // Print the satellite's orbital elements

  // Set the time mesh for the integration
  Real t_start = 0.0;
  Real t_end = 1.0/30; // days
  Real dt = 1.0/(24*60); // days
  VectorX t_mesh = VectorX::LinSpaced((t_end-t_start)/dt + 1, t_start, t_end);

  // Integrate the orbit using the Runge-Kutta solver
  Sandals::RK4<Real, 6, 0> rk4;
  Sandals::Solution<Real, 6, 0> sol;
  Vector6 ics;
  //ics << sat.cartesian().vector();
  ics << sat.equinoctial().vector(), sat.anomaly().L;
  sat.integrate<
    Astro::Coordinates::EQUINOCTIAL, // Integration coordinates
    Astro::Coordinates::CARTESIAN // Output coordinates
    >(rk4, t_mesh, ics, sol);

  // Plot the orbit trace selecting last 3 rows (position) of the solution
  TPolyLine3D* orbit_trace = Plotting::DrawTrace(sol.x.topRows<3>());
  orbit_trace->SetLineColor(kRed);
  orbit_trace->SetLineWidth(1);
  orbit_trace->Draw("same L");

  // Plot the wireframe sphere representing the earth
  TGeoVolume* earth_sphere = Plotting::DrawSphere(ZEROS_VEC3, Planets::Earth_radius_KM, kBlack, 1.0);
  earth_sphere->Draw("same");

  canvas->Update();

  // Make another canvas to plot the magnetic field with respect to time
  TCanvas* canvas2 = new TCanvas("c2", "Magnetic Field vs Time", 800, 600);
  canvas2->SetGrid();
  canvas2->SetTitle("Magnetic Field vs Time");

  // plot the last row of the solution (anomaly L) to check the orbit
  //TGraph* graph = new TGraph(t_mesh.size(), sol.t.data(), sol.x.row(5).data());
  //graph->SetLineColor(kBlue);
  //graph->SetLineWidth(1);
  //graph->Draw("AL");

  // Prepare data for magnetic field plot
  std::vector<Real> time_data;
  std::vector<Real> magnetic_field_data;
  for (Integer i = 0; i < t_mesh.size(); ++i) {
    Vector3 position = sol.x.col(i).head<3>(); // Get position at time t_mesh[i]
    Vector3 magnetic_field = Planets::EarthMagneticFieldDipole(position); // Compute magnetic field
    time_data.push_back(t_mesh[i]);
    magnetic_field_data.push_back(magnetic_field.norm());
  }
  // Create TGraph for magnetic field
  TGraph* graph = new TGraph(time_data.size(), time_data.data(), magnetic_field_data.data());
  graph->SetLineColor(kBlue);
  graph->SetLineWidth(2);
  graph->Draw("AL");

  graph->GetXaxis()->SetTitle("Time (days)");
  graph->GetYaxis()->SetTitle("Magnetic Field (nT)");
  graph->SetTitle("Magnetic Field Strength vs Time");

  canvas2->Update();
  app.Run();

  return 0;
}
