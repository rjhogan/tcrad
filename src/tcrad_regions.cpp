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

#include "tcrad_regions.hpp"

namespace tcrad {

template <bool IsActive> 
void
calc_region_properties(int nlev,
		       int nreg,
		       int ncol,
		       const aMatrix<IsActive>& cloud_fraction,
		       const Matrix& fractional_std,
		       aArray<3,IsActive> reg_fracs,
		       aArray<3,IsActive> od_scaling,
		       Real cloud_fraction_threshold)
{
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
    reg_fracs(__,__,1) = transpose(cloud_fraction);
    reg_fracs(__,__,0) = 1.0 - reg_fracs(__,__,1);
    od_scaling(__,__,0) = 1.0;
  }
  else {
    // We treat the distribution as a gamma then the 16th
    // percentile is close to the following.  Note that because it
    // becomes vanishingly small for FSD >~ 2, we have a lower
    // limit of 1/40, and for higher FSDs reduce the fractional
    // cover of the denser region - see region_fractions routine
    // below
    for (int jcol = 0; jcol < ncol; ++jcol) {
      for (int jlev = 0; jlev < nlev; ++jlev) {
	if (cloud_fraction(jlev,jcol) < cloud_fraction_threshold) {
	  reg_fracs(jcol,jlev,0)  = 1.0; // Clear region
	  reg_fracs(jcol,jlev,1)  = 0.0; // First cloudy region
	  reg_fracs(jcol,jlev,2)  = 0.0; // Second cloudy region
	  od_scaling(jcol,jlev,0) = 1.0; // First cloudy region
	  od_scaling(jcol,jlev,1) = 1.0; // Second cloudy region
	}
	else {
	  reg_fracs(jcol,jlev,0) = 1.0 - cloud_fraction(jlev,jcol);
	  // Fraction and optical-depth scaling of the lower of the
          // two cloudy regions
	  reg_fracs(jcol,jlev,1) = cloud_fraction(jlev,jcol)
	    * max(MIN_LOWER_FRAC, min(MAX_LOWER_FRAC,
	      LOWER_FRAC_FSD_INTERCEPT + fractional_std(jlev,jcol)*LOWER_FRAC_FSD_GRADIENT));

	  od_scaling(jcol,jlev,1) = MIN_GAMMA_OD_SCALING
	    + (1.0 - MIN_GAMMA_OD_SCALING)
	    * exp(-fractional_std(jlev,jcol)*(1.0 + 0.5*fractional_std(jlev,jcol)
					*(1.0 + 0.5*fractional_std(jlev,jcol))));
	  // Fraction of the upper of the two cloudy regions
	  reg_fracs(jcol,jlev,2) = 1.0 - reg_fracs(jcol,jlev,0) - reg_fracs(jcol,jlev,1);
	  // Ensure conservation of the mean optical depth
	  od_scaling(jcol,jlev,2) = (cloud_fraction(jlev,jcol)
				     -reg_fracs(jcol,jlev,1)*od_scaling(jcol,jlev,1))
	    / reg_fracs(jcol,jlev,2);
	}
      }
    }
  }
}
  
}; // namespace tcrad

template
void
tcrad::calc_region_properties<false>(int nlev,
				     int nreg,
				     int ncol,
				     const aMatrix<false>& cloud_fraction,
				     const Matrix& fractional_std,
				     aArray<3,false> reg_fracs,
				     aArray<3,false> od_scaling,
				     Real cloud_fraction_threshold);
#if ADEPT_REAL_TYPE_SIZE == 8
// Instantiate the differentiable function but only in double
// precision
template
void
tcrad::calc_region_properties<true>(int nlev,
				     int nreg,
				     int ncol,
				     const aMatrix<true>& cloud_fraction,
				     const Matrix& fractional_std,
				     aArray<3,true> reg_fracs,
				     aArray<3,true> od_scaling,
				     Real cloud_fraction_threshold);

#endif
