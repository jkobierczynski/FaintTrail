// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// fainttrail — synthetic tracking engine for faint asteroids and comets.
// Umbrella header for the C++ API.
#pragma once

#include "fainttrail/astrometry.hpp"
#include "fainttrail/backend.hpp"
#include "fainttrail/fits.hpp"
#include "fainttrail/io.hpp"
#include "fainttrail/search.hpp"
#include "fainttrail/simulate.hpp"
#include "fainttrail/types.hpp"
#include "fainttrail/validate.hpp"

namespace fainttrail {
// Library version, e.g. "0.1.0".
const char* version();
}  // namespace fainttrail
