// tcrad_regions.cpp - Properties of horizontal regions in each layer
//
// Copyright (C) 2025 ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
//
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.
//
// Author:  Robin Hogan
// Email:   r.j.hogan@ecmwf.int

#ifndef TCRAD_REGIONS_H
#define TCRAD_REGIONS_H 1

#include "tcrad_adept.hpp"

namespace tcrad {

// ---------------------------------------------------------------------
// Compute the optical depth scalings for the optically "thick" and
// "thin" regions of a Tripleclouds representation of a sub-grid PDF
// of cloud optical depth. Following Shonk and Hogan (2008), the 16th
// percentile is used for the thin region, and the formulas estimate
// this for both lognormal and gamma distributions. However, an
// adjustment is needed for the gamma distribution at large fractional
// standard deviations.
template <bool IsActive> 
void
calc_region_properties(int nlev,
		       int nreg,
		       const Array<1,IsActive>& cloud_fraction,
		       const Array<1>& fractional_std,
		       Array<2,IsActive> region_fracs,
		       Array<2,IsActive> od_scaling,
		       Real cloud_fraction_threshold);

};

#endif
