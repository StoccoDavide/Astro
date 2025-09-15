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
#ifdef ASTRO_ENABLE_PLOTTING
#include "Astro/Plotting.hh"
#endif // ASTRO_ENABLE_PLOTTING

using namespace Astro;


#ifdef ASTRO_ENABLE_PLOTTING
  int main(int argc, char** argv)

    // Initialize the ROOT application
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
  int main()
#endif // ASTRO_ENABLE_PLOTTING
  {

  // Create the Sun-Earth-Moon system
  Body earth(Planets::Earth());
  Body moon(Planets::Moon());

  // Satellite formation orbiting the Moon
  Real payload_mass{100.0*1.339E-6}; // Mass of the satellite in Kg
  Real payload_radius{1E-3}; // Radius of the satellite in Km
  Body sat("sat", payload_mass, payload_radius);
  sat.set_orbit().set_factor(Factor::POSIGRADE);
  sat.set_orbit().set_mu(Planets::Moon_mu_AU3_DAY2);
  sat.set_orbit().set_keplerian(OrbitalElements::Keplerian(
    Planets::Moon_radius_AU + KM_To_AU(400.0), // Semi-major axis in AU
    1.36E-2, // Eccentricity
    Deg_To_Rad(9.015), // Inclination in radians
    0.0, // Longitude of ascending node in radians
    0.0  // Argument of periapsis in radians
    ),
    0.0 // True anomaly (rad)
  );

  // Set the anomaly for the Earth and Moon
  sat.set_epoch_anomaly().set_nu(0.0, sat.orbit().keplerian(), Factor::POSIGRADE);

  // Set the time mesh for the integration
  Real t_start{0.0};
  Real t_end{365.25}; // One year in days
  Real dt{1.0}; // Time step in days
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
  for (int i{0}; i < t_mesh.size(); ++i) {
    pos_earth.emplace_back(earth.cartesian_position(t_mesh[i]));
    pos_moon.emplace_back(moon.cartesian_position(t_mesh[i]));
    pos_sat.emplace_back(sat.cartesian_position(t_mesh[i]));
    vel_earth.emplace_back(earth.cartesian_velocity(t_mesh[i]));
    vel_moon.emplace_back(moon.cartesian_velocity(t_mesh[i]));
    vel_sat.emplace_back(sat.cartesian_velocity(t_mesh[i]));
    rtn_earth.emplace_back(earth.orbit().cartesian_to_frenet_rtn(pos_earth[i], vel_earth[i]));
    rtn_moon.emplace_back(moon.orbit().cartesian_to_frenet_rtn(pos_moon[i], vel_moon[i]));
    rtn_sat.emplace_back(sat.orbit().cartesian_to_frenet_rtn(pos_sat[i], vel_sat[i]));
  }

  // Compute periods
  std::cout << "Earth period (days): " << earth.orbit().period() << "\n";
  std::cout << "Moon period (days): " << moon.orbit().period() << "\n";
  std::cout << "Sat. formation period (days): " << sat.orbit().period() << "\n";

  // Save the cartesian states to CSV files
  std::ofstream of_moon("satellite_positions_moon_coords.csv");
  std::ofstream of_earth("satellite_positions_earth_coords.csv");
  std::ofstream of_sun("satellite_positions_sun_coords.csv");
  of_moon << "T,X,Y,Z,RX,RY,RZ,TX,TY,TZ,NX,NY,NZ\n";
  of_earth << "T,X,Y,Z,RX,RY,RZ,TX,TY,TZ,NX,NY,NZ\n";
  of_sun << "T,X,Y,Z,RX,RY,RZ,TX,TY,TZ,NX,NY,NZ\n";
  Vector3 pos, vel;
  Matrix3 rtn;
  for (int i{0}; i < t_mesh.size(); ++i) {
    pos = pos_sat[i];
    rtn = rtn_sat[i];
    of_moon <<
      t_mesh[i] << "," << pos.x() << "," << pos.y() << "," << pos.z() << "," <<
      rtn.row(0).x() << "," << rtn.row(0).y() << "," << rtn.row(0).z() << "," <<
      rtn.row(1).x() << "," << rtn.row(1).y() << "," << rtn.row(1).z() << "," <<
      rtn.row(2).x() << "," << rtn.row(2).y() << "," << rtn.row(2).z() << "\n";
    pos = pos_earth[i] + rtn_earth[i] * pos_sat[i];
    rtn = rtn_earth[i] * rtn_sat[i];
    of_earth <<
      t_mesh[i] << "," << pos.x() << "," << pos.y() << "," << pos.z() << "," <<
      rtn.row(0).x() << "," << rtn.row(0).y() << "," << rtn.row(0).z() << "," <<
      rtn.row(1).x() << "," << rtn.row(1).y() << "," << rtn.row(1).z() << "," <<
      rtn.row(2).x() << "," << rtn.row(2).y() << "," << rtn.row(2).z() << "\n";
    pos = pos_earth[i] + rtn_earth[i] * (pos_moon[i] + rtn_moon[i] * pos_sat[i]);
    rtn = rtn_earth[i] * rtn_moon[i] * rtn_sat[i];
    of_sun <<
      t_mesh[i] << "," << pos.x() << "," << pos.y() << "," << pos.z() << "," <<
      rtn.row(0).x() << "," << rtn.row(0).y() << "," << rtn.row(0).z() << "," <<
      rtn.row(1).x() << "," << rtn.row(1).y() << "," << rtn.row(1).z() << "," <<
      rtn.row(2).x() << "," << rtn.row(2).y() << "," << rtn.row(2).z() << "\n";
  }
  of_moon.close();
  of_earth.close();
  of_sun.close();

#ifdef ASTRO_ENABLE_PLOTTING
  // Plot the orbit trace selecting last 3 rows (position) of the solution
  TPolyLine3D* poly_line = new TPolyLine3D(pos_earth.size());
  for (int i{0}; i < t_mesh.size(); ++i) {
    Vector3 const pos(pos_earth[i] + rtn_earth[i] * (pos_moon[i] + rtn_moon[i] * pos_sat[i]));
    poly_line->SetPoint(i, pos.x(), pos.y(), pos.z());
  }
  poly_line->SetLineColor(kRed);
  poly_line->SetLineWidth(2);
  poly_line->Draw("same L");

  canvas->Update();
  app.Run();
#endif // ASTRO_ENABLE_PLOTTING

  return 0;
}
