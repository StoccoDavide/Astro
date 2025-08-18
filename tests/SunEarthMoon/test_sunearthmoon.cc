/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Astro project is distributed under the GNU GPLv3.                                         *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <fstream>

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

  // Draw the absolute reference frame
  TObjArray* abs_axes = Astro::Plotting::DrawAbsoluteAxes(0.25);
  abs_axes->SetOwner(kTRUE); // Ensure axes are deleted with the canvas

  // Create the Sun-Earth-Moon system
  Body earth(Planets::Earth());
  Body moon(Planets::Moon());

  // Satellite orbiting the Moon
  Real payload_mass{100.0*1.339E-6}; // Mass of the satellite in Kg
  Real payload_radius{1E-3}; // Radius of the satellite in Km
  Body sat("sat", payload_mass, payload_radius);
  sat.set_factor(Factor::POSIGRADE);
  sat.set_mu(Planets::Moon_mu_AU3_DAY2);
  sat.set_keplerian(OrbitalElements::Keplerian(
    Planets::Moon_radius_AU + KM_To_AU(400.0), // Semi-major axis in AU
    1.36E-2, // Eccentricity
    Deg_To_Rad(9.015), // Inclination in radians
    0.0, // Longitude of ascending node in radians
    0.0  // Argument of periapsis in radians
  ));

  // Set the anomaly for the Earth and Moon
  earth.set_anomaly().set_nu(0.0, earth.keplerian(), Factor::POSIGRADE);
  moon.set_anomaly().set_nu(0.0, earth.keplerian(), Factor::POSIGRADE);
  sat.set_anomaly().set_nu(0.0, sat.keplerian(), Factor::POSIGRADE);

  // Set the time mesh for the integration
  Real t_start{0.0};
  Real t_end{365.0}; // One year in days
  Real dt{0.1}; // Time step in days
  VectorX t_mesh = VectorX::LinSpaced((t_end-t_start)/dt + 1, t_start, t_end);

  // Prepare the solution object
  std::vector<Vector3> pos_earth, pos_moon, vel_earth, vel_moon, vel_sat, pos_sat;
  std::vector<Matrix3> rtn_earth, rtn_moon, rtn_sat;
  pos_earth.reserve(t_mesh.size());
  pos_moon.reserve(t_mesh.size());
  pos_sat.reserve(t_mesh.size());
  vel_earth.reserve(t_mesh.size());
  vel_moon.reserve(t_mesh.size());
  vel_sat.reserve(t_mesh.size());
  rtn_earth.reserve(t_mesh.size());
  rtn_moon.reserve(t_mesh.size());
  rtn_sat.reserve(t_mesh.size());

  // Propagate the orbits analytically
  pos_earth.emplace_back(earth.cartesian().r);
  pos_moon.emplace_back(moon.cartesian().r);
  pos_sat.emplace_back(sat.cartesian().r);
  vel_earth.emplace_back(earth.cartesian().v);
  vel_moon.emplace_back(moon.cartesian().v);
  vel_sat.emplace_back(sat.cartesian().v);
  rtn_earth.emplace_back(earth.cartesian_to_frenet_rtn(earth.cartesian().r, earth.cartesian().v));
  rtn_moon.emplace_back(moon.cartesian_to_frenet_rtn(moon.cartesian().r, moon.cartesian().v));
  rtn_sat.emplace_back(sat.cartesian_to_frenet_rtn(sat.cartesian().r, sat.cartesian().v));
  for (int i{1}; i < t_mesh.size(); ++i) {
    Real dti{t_mesh[i] - t_mesh[i-1]};
    earth.propagate(dti);
    moon.propagate(dti);
    sat.propagate(dti);

    pos_earth.emplace_back(earth.cartesian().r);
    pos_moon.emplace_back(moon.cartesian().r);
    pos_sat.emplace_back(sat.cartesian().r);
    vel_earth.emplace_back(earth.cartesian().v);
    vel_moon.emplace_back(moon.cartesian().v);
    vel_sat.emplace_back(sat.cartesian().v);
    rtn_earth.emplace_back(earth.cartesian_to_frenet_rtn(earth.cartesian().r, earth.cartesian().v));
    rtn_moon.emplace_back(moon.cartesian_to_frenet_rtn(moon.cartesian().r, moon.cartesian().v));
    rtn_sat.emplace_back(sat.cartesian_to_frenet_rtn(sat.cartesian().r, sat.cartesian().v));
  }

  // Plot the orbit trace selecting last 3 rows (position) of the solution
  std::ofstream ofs("satellite_positions.csv");
  ofs << "Time,X,Y,Z\n";
  TPolyLine3D* poly_line = new TPolyLine3D(pos_earth.size());
  for (int k{0}; k < pos_earth.size(); ++k) {
    Vector3 const & pos = pos_earth[k] + rtn_earth[k]*(pos_moon[k] + 0*rtn_moon[k]*pos_sat[k]);
    poly_line->SetPoint(k, pos.x(), pos.y(), pos.z());
    ofs << t_mesh[k] << "," << pos.x() << "," << pos.y() << "," << pos.z() << "\n";
  }
  poly_line->SetLineColor(kRed);
  poly_line->SetLineWidth(2);
  poly_line->Draw("same L");
  ofs.close();

  canvas->Update();
  app.Run();
  return 0;
}
