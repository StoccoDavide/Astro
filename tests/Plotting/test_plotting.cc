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

int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);

    // Create canvas
    TCanvas* canvas = new TCanvas("c1", "Plotting test", 800, 600);
    canvas->SetGrid();

    // Set up 3D axis
    TAxis3D rulers; rulers.Draw();
    TAxis3D::ToggleRulers();
    TAxis3D::ToggleZoom();

    // Create a Body object (e.g., Earth)
    Astro::Body earth(Astro::Planets::Earth());

    // Retrieve the orbit trace
    TPolyLine3D* orbitTrace = Astro::Plotting::Trace(earth, 1000);
    orbitTrace->SetLineColor(kBlack);
    orbitTrace->SetLineWidth(1);
    orbitTrace->Draw("L");

    // Draw the absolute reference frame
    TObjArray* abs_axes = Astro::Plotting::DrawAbsoluteAxes();
    abs_axes->SetOwner(kTRUE); // Ensure axes are deleted with the canvas

    // Draw the orbit plane axes
    TObjArray* orbit_plane_axes = Astro::Plotting::DrawOrbitalPlaneAxes(earth);
    orbit_plane_axes->SetOwner(kTRUE); // Ensure axes are deleted with the canvas

    // Draw the Frenet-Serret frame at the first point of the orbit
    TObjArray* fs_axes = Astro::Plotting::DrawFrenetSerretAxes(earth);
    fs_axes->SetOwner(kTRUE); // Ensure axes are deleted with the canvas

    // Draw the sphere representing the Earth
    Astro::Plotting::DrawMarker(earth, kBlue, 2.0);

    canvas->Update();
    app.Run();
    return 0;
}
