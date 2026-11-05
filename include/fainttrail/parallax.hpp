// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Curved tracks of near-Earth objects.
//
// Seen from a point on the rotating Earth, an object at geocentric distance D is
// displaced from its geocentric direction by the observer's offset from the
// Earth's centre, projected on the sky and divided by D: up to 8.8" at 1 AU, 15'
// at 0.01 AU. Over a sequence of hours the displacement changes with the Earth's
// rotation, so the track of a close object is a straight line plus 1/D times a
// curve that depends only on the site, the time and the field: the parallax
// model. One extra search dimension (q = 1/D) therefore covers all curvatures
// that diurnal parallax can produce, instead of two free acceleration terms.
// The curvature of the geocentric orbit itself is not modelled; over a few hours
// it is much smaller than the parallax for the objects this is meant for.
#pragma once

#include <string>
#include <vector>

#include "fainttrail/fits.hpp"
#include "fainttrail/search.hpp"
#include "fainttrail/types.hpp"

namespace fainttrail {

struct ObserverSite {
    bool valid = false;
    double lat_deg = 0.0;   // geodetic latitude, north positive
    double lon_deg = 0.0;   // longitude, EAST positive
    double height_m = 0.0;  // above the WGS84 ellipsoid
};

// From FITS keywords: OBSGEO-B/OBSGEO-L/OBSGEO-H (FITS standard), SITELAT/SITELONG/
// SITEELEV (NINA, Ekos, MaxIm), LAT-OBS/LONG-OBS/ALT-OBS. Values in degrees or
// sexagesimal strings; longitudes east positive.
bool site_from_header(const FitsHeader& header, ObserverSite& site);
// "51.0,5.3" or "51.0,5.3,40" (lat, east lon, height in m); throws on bad input.
ObserverSite parse_site(const std::string& text);

// Geocentric position of the observer in AU, equatorial (of date) axes, at a UTC
// time (MJD); UT1 = UTC is assumed (the difference is below a second).
void observer_position_au(const ObserverSite& site, double mjd_utc, double xyz[3]);

// Apparent displacement in pixels of an object at 1 AU in the direction of
// pixel (x, y), for each time (MJD UTC): the raw parallax, straight part included.
void parallax_pixels(const Wcs& wcs, double x, double y, const ObserverSite& site, const std::vector<double>& mjd,
                     std::vector<double>& px, std::vector<double>& py);

// Fill p.par_x / p.par_y for the field centre, with the straight-line part over
// the sequence removed (it is absorbed by the velocity grid). Needs p.wcs.
void set_parallax_model(Prepared& p, const ObserverSite& site);

// Largest departure from a straight line, in pixels, of an object at 1 AU.
double parallax_amplitude_px(const Prepared& p);

// Trial inverse distances (1/AU) of the near-Earth search: 0 (straight tracks)
// first, then steps that move a track by at most `grid_tol_px` up to
// 1 / params.neo_min_dist_au. {0} when the search is off. Throws if the search is
// on but there is no parallax model.
std::vector<double> inverse_distance_grid(const Prepared& p, const SearchParams& params, double* dq = nullptr);

}  // namespace fainttrail
