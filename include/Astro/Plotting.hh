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

#ifdef ASTRO_ENABLE_PLOTTING

#ifndef ASTRO_PLOTTING_HH
#define ASTRO_PLOTTING_HH

#include "Astro.hh"
#include "Astro/Orbit.hh"

// ROOT library
#include <TApplication.h>
#include <TCanvas.h>
#include <TAxis3D.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoSphere.h>
#include <TGeoMatrix.h>
#include <TMath.h>

namespace Astro
{
  namespace Plotting
  {
    /**
    * \brief Draw a referce frame in the ROOT canvas.
    *
    * This function draws a reference frame in the ROOT canvas, with the X-axis in red, Y-axis in green,
    * and Z-axis in blue. The axes are drawn as arrows originating from the specified origin point, with a
    * specified rotation and length.
    * \param[in] origin The origin point of the axes.
    * \param[in] rotation The rotation matrix for the axes.
    * \param[in] length The length of the axes.
    * \param[in] line_width The width of the lines.
    * \return A TObjArray containing the drawn axes.
    * \note The function creates a new TObjArray containing the drawn axes, which can be stored or deleted later.
    */
    TObjArray* DrawAxes(const Vector3 & origin, const Rotation & rotation, Real length = 1.0, Real line_width = 1.0) {
      TObjArray* arrows = new TObjArray();

      // X-axis (Red)
      TPolyLine3D* xAxis = new TPolyLine3D(2);
      xAxis->SetPoint(0, origin.x(), origin.y(), origin.z());
      Vector3 x_end(origin + rotation * Vector3(length, 0, 0));
      xAxis->SetPoint(1, x_end.x(), x_end.y(), x_end.z());
      xAxis->SetLineColor(kRed);
      xAxis->SetLineWidth(line_width);
      xAxis->Draw("same L");
      arrows->Add(xAxis);

      // Y-axis (Green)
      TPolyLine3D* yAxis = new TPolyLine3D(2);
      yAxis->SetPoint(0, origin.x(), origin.y(), origin.z());
      Vector3 y_end(origin + rotation * Vector3(0, length, 0));
      yAxis->SetPoint(1, y_end.x(), y_end.y(), y_end.z());
      yAxis->SetLineColor(kGreen + 2);
      yAxis->SetLineWidth(line_width);
      yAxis->Draw("same L");
      arrows->Add(yAxis);

      // Z-axis (Blue)
      TPolyLine3D* zAxis = new TPolyLine3D(2);
      zAxis->SetPoint(0, origin.x(), origin.y(), origin.z());
      Vector3 z_end(origin + rotation * Vector3(0, 0, length));
      zAxis->SetPoint(1, z_end.x(), z_end.y(), z_end.z());
      zAxis->SetLineColor(kBlue);
      zAxis->SetLineWidth(line_width);
      zAxis->Draw("same L");
      arrows->Add(zAxis);

      return arrows; // You can store/delete later
    }

    /**
    * \brief Draw the absolute reference frame in the ROOT canvas.
    *
    * This function draws the absolute reference frame in the ROOT canvas.
    * \param[in] length The length of the axes.
    * \param[in] line_width The width of the lines.
    * \return A TObjArray containing the drawn axes.
    */
    TObjArray* DrawAbsoluteAxes(Real length = 1.0, Real line_width = 1) {
      return DrawAxes(ZEROS_VEC3, IDENTITY_MAT3, length, line_width);
    }

    /**
    * \brief Draw the orbital plane axes in the ROOT canvas.
    *
    * This function draws the orbital plane axes in the ROOT canvas.
    * \param[in] orbit The orbit body whose orbital plane axes are to be drawn.
    * \param[in] length The length of the axes.
    * \param[in] line_width The width of the lines.
    * \return A TObjArray containing the drawn axes.
    */
    TObjArray* DrawOrbitalPlaneAxes(const Orbit & orbit, Real length = 1.0, Real line_width = 1.0) {
      // Get the orbital plane rotation matrix
      Rotation rotation{orbit.keplerian_to_reference()};

      // Draw the axes in the orbital plane
      return DrawAxes(ZEROS_VEC3, rotation, length, line_width);
    }

    /**
    * \brief Plot the orbit trace of an astronomical body.
    *
    * This function plots the orbit trace of an astronomical body using the ROOT library.
    *
    * \param[in] orbit The astronomical body whose orbit is to be plotted.
    */
    TPolyLine3D* Trace(const Orbit & orbit, Integer num_points = 360) {

      // Extract orbital elements
      Real a{orbit.keplerian().a};
      Real e{orbit.keplerian().e};
      Real i{orbit.keplerian().i};
      Real Omega{orbit.keplerian().Omega};
      Real omega{orbit.keplerian().omega};

      // Precompute rotation values
      Real cos_Omega{std::cos(Omega)};
      Real sin_Omega{std::sin(Omega)};
      Real cos_i{std::cos(i)};
      Real sin_i{std::sin(i)};

      TPolyLine3D* poly_line = new TPolyLine3D(num_points + 1);
      for (int k = 0; k <= num_points; ++k) {
        // Fictitious true anomaly in [0, 2pi]
        Real nu{2.0 * TMath::Pi() * k / num_points};

        // Radius from orbit equation
        Real r{a*(1.0 - e*e)/(1.0 + e * std::cos(nu))};

        // Rotate to inertial frame
        Real cos_w_nu{std::cos(omega + nu)};
        Real sin_w_nu{std::sin(omega + nu)};
        Real x{r*(cos_Omega*cos_w_nu - sin_Omega*sin_w_nu*cos_i)};
        Real y{r*(sin_Omega*cos_w_nu + cos_Omega*sin_w_nu*cos_i)};
        Real z{r*(sin_w_nu*sin_i)};

        poly_line->SetPoint(k, x, y, z);
      }

      poly_line->SetLineColor(kRed);
      poly_line->SetLineWidth(2);
      return poly_line;
    }

    /**
    * \brief Plot the Frenet-Serret frame of an orbit body.
    *
    * This function plots the Frenet-Serret frame of an orbit body using the ROOT library.
    * \param[in] orbit The orbit body whose Frenet-Serret frame is to be plotted.
    * \param[in] length The length of the axes.
    * \param[in] line_width The width of the lines.
    * \return A TObjArray containing the drawn axes of the Frenet-Serret frame.
    */
    TObjArray* DrawFrenetSerretAxes(const Orbit & orbit, Real length = 0.1, Real line_width = 1.0) {
      return DrawAxes(orbit.cartesian().r, orbit.cartesian_to_frenet_rtn(), length, line_width);
    }

    /**
    * \brief Plot a marker in 3D space.
    *
    * Plots a marker in 3D space at the specified position using the ROOT library.
    * \param[in] body The celestial body whose position is to be plotted.
    * \param[in] color The color of the marker.
    * \param[in] size The size of the marker.
    * \return A TPolyMarker3D pointer to the drawn marker.
    */
    TPolyMarker3D* DrawMarker(const Vector3 & postion, Color_t color, Real size = 1.0) {
      TPolyMarker3D* marker = new TPolyMarker3D(1);
      marker->SetPoint(0, postion.x(), postion.y(), postion.z());
      marker->SetMarkerColor(color);
      marker->SetMarkerSize(size);
      marker->SetMarkerStyle(20); // Default marker style
      marker->Draw("same P");
      return marker;
    }

    /**
    * \brief Plot a marker representing a celestial body at the origin.
    *
    * Plots a marker representing a celestial body at the origin using the ROOT library.
    * \param[in] body The celestial body whose position is to be plotted.
    * \param[in] color The color of the marker.
    * \param[in] size The size of the marker.
    * \return A TPolyMarker3D pointer to the drawn marker.
    */
    TPolyMarker3D* DrawMarker(const Body & body, Color_t color, Real size = 1.0) {
      return DrawMarker(body.cartesian().r, color, size);
    }

    /**
    * \brief Plot a sphere representing a celestial body at a specific position.
    *
    * This function plots a sphere representing a celestial body at a specific position using the ROOT library.
    * \param[in] position The position of the sphere in 3D space.
    * \param[in] radius The radius of the sphere.
    * \param[in] color The color of the sphere.
    * \param[in] line_width The width of the lines.
    * \return A TGeoVolume pointer to the drawn sphere.
    */
    TGeoVolume* DrawSphere(Vector3 position, Real radius, Color_t color, Real ) {
      static TGeoManager* geom = nullptr;

      // Create TGeoManager only once
      if (!geom) {
          geom = new TGeoManager("SphereGeom", "Sphere Geometry");

          // Material and medium
          TGeoMaterial *mat = new TGeoMaterial("Vacuum", 0, 0, 0);
          TGeoMedium *med = new TGeoMedium("VacuumMedium", 1, mat);

          // Top volume large enough to contain everything
          TGeoVolume *top = geom->MakeBox("TOP", med, 1000, 1000, 1000);
          geom->SetTopVolume(top);
      }

      // Ensure the geometry is initialized
      TGeoMedium* med = geom->GetMedium("VacuumMedium");
      if (!med) {med = new TGeoMedium("VacuumMedium", 1, geom->GetMaterial("Vacuum"));}

      // Create a unique sphere
      static int sphere_id = 0;
      std::string sphere_name = "Sphere_" + std::to_string(sphere_id++);
      TGeoSphere* sphere = new TGeoSphere(0.0, radius);
      TGeoVolume* volume = new TGeoVolume(sphere_name.c_str(), sphere, med);
      volume->SetLineColor(color);
      //volume->SetLineWidth(line_width);

      // Add to top volume
      geom->GetTopVolume()->AddNode(volume, 1, new TGeoTranslation(position.x(), position.y(), position.z()));

      // Draw top volume in wireframe
      geom->CloseGeometry();
      geom->GetTopVolume()->Draw("same");

      return volume;
    }

    /**
    * \brief Plot a sphere representing a celestial body at the origin.
    *
    * This function plots a sphere representing a celestial body at the origin using the ROOT library.
    * \param[in] body The celestial body whose position is to be plotted.
    * \param[in] color The color of the sphere.
    * \param[in] line_width The width of the lines.
    * \return A TGeoVolume pointer to the drawn sphere.
    */
    TGeoVolume* DrawSphere(const Body & body, Color_t color, Real line_width = 1.0) {
      return DrawSphere(body.cartesian().r, body.radius(), color, line_width);
    }


  } // namespace Plotting

} // namespace Astro

#endif // ASTRO_PLOTTING_HH

#endif // ASTRO_ENABLE_PLOTTING
