// tcrad__regions.cpp -- Properties of horizontal regions in Tripleclouds & SPARTACUS
//
// (C) Copyright 2016- ECMWF.
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
//
// This is adapted from ecRad's radiation_regions.F90

#include "tcrad_regions.hpp"

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
		       Real cloud_fraction_threshold)
{
  // See appendix of Hogan et al. (JAS, 2019)  
  static const Real MIN_GAMMA_OD_SCALING = 0.025;
  static const Real MIN_LOWER_FRAC = 0.5;
  static const Real MAX_LOWER_FRAC = 0.9;
  static const Real FSD_AT_MIN_LOWER_FRAC = 1.5;
  static const Real FSD_AT_MAX_LOWER_FRAC = 3.725;
  static const Real LOWER_FRAC_FSD_GRADIENT 
    = (MAX_LOWER_FRAC - MIN_LOWER_FRAC) / (FSD_AT_MAX_LOWER_FRAC - FSD_AT_MIN_LOWER_FRAC);
  static const Real LOWER_FRAC_FSD_INTERCEPT
    = MIN_LOWER_FRAC - FSD_AT_MIN_LOWER_FRAC * LOWER_FRAC_FSD_GRADIENT;

  if (nreg == 2) {
    // Only one clear-sky and one cloudy region: cloudy region is
    // homogeneous
    region_fracs(__,1) = cloud_fraction;
    region_fracs(__,0) = 1.0 - cloud_fraction;
    od_scaling(__,0) = 1.0; // First cloudy region
  }
  else {
    // We treat the distribution as a gamma then the 16th
    // percentile is close to the following.  Note that because it
    // becomes vanishingly small for FSD >~ 2, we have a lower
    // limit of 1/40, and for higher FSDs reduce the fractional
    // cover of the denser region - see region_fractions routine
    // below
    for (int jlev = 0; jlev < nlev; ++jlev) {
      if (cloud_fraction(jlev) < cloud_fraction_threshold) {
	region_fracs(jlev,0) = 1.0; // Clear region
	region_fracs(jlev,1) = 0.0; // First cloudy region
	region_fracs(jlev,2) = 0.0; // Second cloudy region
	od_scaling(jlev,0)   = 1.0; // First cloudy region
	od_scaling(jlev,1)   = 1.0; // Second cloudy region
      }
      else {
	region_fracs(jlev,0) = 1.0 - cloud_fraction(jlev);
	// Fraction and optical-depth scaling of the first of the
	// two cloudy regions (least optical depth)
	region_fracs(jlev,1) = cloud_fraction(jlev)
	  * max(MIN_LOWER_FRAC,
		min(MAX_LOWER_FRAC,
		    LOWER_FRAC_FSD_INTERCEPT
		    + fractional_std(jlev)*LOWER_FRAC_FSD_GRADIENT));
	od_scaling(jlev,0) = MIN_GAMMA_OD_SCALING
	  + (1.0 - MIN_GAMMA_OD_SCALING)
	  * exp(-fractional_std(jlev)
		*(1.0 + 0.5*fractional_std(jlev)
		  *(1.0 + 0.5*fractional_std(jlev))));
	// Fraction of the upper of the two cloudy regions
	region_fracs(jlev,2)
	  = 1.0 - region_fracs(jlev,0) - region_fracs(jlev,1);
	// Ensure conservation of the mean optical depth
	od_scaling(jlev,1)
	  = (cloud_fraction(jlev)
	     -region_fracs(jlev,1)*od_scaling(jlev,0))
	  / region_fracs(jlev,2);
      }
    } // levels
  }
}
  
}; // namespace tcrad

// Instantiate the direct function
template
void
tcrad::calc_region_properties<false>(int nlev,
				     int nreg,
				     const Array<1,false>& cloud_fraction,
				     const Array<1>& fractional_std,
				     Array<2,false> region_fracs,
				     Array<2,false> od_scaling,
				     Real cloud_fraction_threshold);
#if ADEPT_REAL_TYPE_SIZE == 8
// Instantiate the differentiable function but only in double
// precision
template
void
tcrad::calc_region_properties<true>(int nlev,
				    int nreg,
				    const Array<1,true>& cloud_fraction,
				    const Array<1>& fractional_std,
				    Array<2,true> region_fracs,
				    Array<2,true> od_scaling,
				    Real cloud_fraction_threshold);

#endif
