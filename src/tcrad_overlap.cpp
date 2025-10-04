// tcrad_overlap.F90 - Compute cloud overlap quantities
//
// (C) Copyright 2014- ECMWF.
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
// This is adapted from ecRad's radiation_overlap.F90

#include "tcrad_overlap.hpp"
#include "tcrad_config.hpp"

namespace tcrad {

// ---------------------------------------------------------------------
// Calculate a matrix expressing the overlap of regions in adjacent
// layers, using the Hogan and Illingworth (2000) "alpha" overlap
// parameter, but allowing for the two cloudy regions in the
// Tripleclouds assumption to have different areas. "op" is the
// overlap parameter of (0) cloud boundaries and (1) internal cloud
// heterogeneities, while frac_upper and frac_lower are the area
// fractions of the three regions in the upper and lower layers,
// respectively, each summing to one.
template <bool IsActive>
Matrix33<IsActive>
calc_alpha_overlap_matrix(const Vector2<IsActive>& op,
			  const Vector3<IsActive>& frac_upper,
			  const Vector3<IsActive>& frac_lower,
			  Real cloud_fraction_threshold) {

  typedef typename scalar<IsActive>::type Scalar;

  // Output overlap matrix
  Matrix33<IsActive> overlap_matrix;

  // Cloud fraction is the sum of the two cloudy fractions
  Scalar cf_upper = frac_upper(1)+frac_upper(2);
  Scalar cf_lower = frac_lower(1)+frac_lower(2);

  // Combined cloud cover of the two layers
  Scalar pair_cloud_cover = op(0)*max(cf_upper,cf_lower)
    + (1.0 - op(0)) * (cf_upper + cf_lower - cf_upper*cf_lower);

  // Element of overlap matrix representing the area fraction of the
  // gridbox that is clear in both layers
  overlap_matrix(0,0) = 1.0 - pair_cloud_cover;

  // Clear in upper layer, cloudy in lower layer
  Scalar one_over_cf = 1.0 / max(cf_lower, cloud_fraction_threshold);
  overlap_matrix(1,0) = (pair_cloud_cover - cf_upper)
    * frac_lower(1) * one_over_cf;
  overlap_matrix(2,0) = (pair_cloud_cover - cf_upper)
    * frac_lower(2) * one_over_cf;

  // Clear in lower layer, cloudy in upper
  one_over_cf = 1.0 / max(cf_upper, cloud_fraction_threshold);
  overlap_matrix(0,1) = (pair_cloud_cover - cf_lower)
    * frac_upper(1) * one_over_cf;
  overlap_matrix(0,2) = (pair_cloud_cover - cf_lower)
    * frac_upper(2) * one_over_cf;

  // Cloudy in both layers: frac_both is the fraction of the fraction
  // of the gridbox with cloud in both layers
  Scalar frac_both = cf_upper + cf_lower - pair_cloud_cover;
  // Treat low and high optical-depth regions within frac_both as one
  // treats clear and cloudy skies in the whole domain; redefine the
  // following variables treating the high optical-depth region as the
  // cloud
  cf_upper = frac_upper(2) / max(cf_upper, cloud_fraction_threshold);
  cf_lower = frac_lower(2) / max(cf_lower, cloud_fraction_threshold);
  pair_cloud_cover = op(1) * max(cf_upper,cf_lower)
    + (1.0 - op(1)) * (cf_upper + cf_lower - cf_upper * cf_lower);
  // Assign overlaps for this 2x2 section of the 3x3 matrix as for
  // the 2-region case above, but multiplied by frac_both
  overlap_matrix(1,1) = frac_both * (1.0 - pair_cloud_cover);
  overlap_matrix(2,1) = frac_both * (pair_cloud_cover - cf_upper);
  overlap_matrix(1,2) = frac_both * (pair_cloud_cover - cf_lower);
  overlap_matrix(2,2) = frac_both * (cf_upper + cf_lower - pair_cloud_cover);

  return overlap_matrix;
}
  
//---------------------------------------------------------------------
// Compute the upward and downward overlap matrices u_overlap and
// v_overlap, respectively, where u_overlap is defined as in Hogan et
// al. (JAS, 2016) such that y=u_overlap*x, where x is a vector of
// upwelling fluxes in each region just below an interface, and y is a
// vector of upwelling fluxes in each region just above that
// interface. HOWEVER in translating from Fortran to C++, all
// dimensions were reversed so actually the orientation of the matrix
// is such that y=transpose(u_overlap)*x. For nlev model levels there
// are nlev+1 interfaces including the ground and top-of-atmosphere,
// and so that is one of the dimensions of u_overlap and v_overlap.
template <bool IsActive>
void
calc_overlap_matrices(int nlev,
		      const Array<2,IsActive>& region_fracs,
		      const Array<1,IsActive>& overlap_param,
		      Array<3,IsActive> u_overlap,
		      Array<3,IsActive> v_overlap,
		      Real cloud_fraction_threshold)
{
  // Overlap matrix (non-directional)
  Matrix33<IsActive> overlap_matrix;
  
  // Fraction of the gridbox occupied by each region in the upper and
  // lower layers for an interface
  Vector3<IsActive> frac_upper, frac_lower;

  // Overlap parameter for first two regions
  Vector2<IsActive> op;

  // Outer space is treated as one clear-sky region, so the fractions
  // are assigned as such
  frac_upper(0) = 1.0;
  frac_upper(range(1,end)) = 0.0;
  
  // Overlap parameter is irrelevant when there is only one region
  // in the upper layer
  op = 1.0;

  // Loop down through the atmosphere, where jlev indexes each
  // half-level starting at 1 for the top-of-atmosphere, as well
  // as indexing each level starting at 1 for the top-most level.
  for (int jlev = 0; jlev < nlev+1; ++jlev) {
    // Fraction of each region just below the interface
    if (jlev >= nlev) {
      // We are at the surface: treat as a single clear-sky region
      frac_lower(0) = 1.0;
      frac_lower(range(1,end)) = 0.0;
    }
    else {
      frac_lower = region_fracs(jlev,__);
    }
    
    // Compute the overlap parameter of the interface just below the
    // current full level
    if (jlev == 0 || jlev >= nlev) {
      // We are at the surface or top-of-atmosphere: overlap
      // parameter is irrelevant
      op = 1.0;
    }
    else {
      // We are not at the surface
      op(0) = overlap_param(jlev-1);
      // For cloudy regions, scale the cloud-boundary overlap
      // parameter to obtain the cloud-inhomogeneity overlap
      // parameter as follows
      if (op(0) >= 0.0) {
	op(1) = pow(op(0), 1.0/DECORRELATION_SCALING);
      }
      else {
	op(1) = op(0);
      }
    }
    overlap_matrix = calc_alpha_overlap_matrix(op, frac_upper, frac_lower,
					       cloud_fraction_threshold);
    
    // Convert to directional overlap matrices
    for (int jupper = 0; jupper < NREGIONS; ++jupper) {
      for (int jlower = 0; jlower < NREGIONS; ++jlower) {
	if (frac_lower(jlower) >= cloud_fraction_threshold) {
	  u_overlap(jlev,jlower,jupper) = overlap_matrix(jlower,jupper)
	    / frac_lower(jlower);
	}
	else {
	  u_overlap(jlev,jlower,jupper) = 0.0;
	}
	if (frac_upper(jupper) >= cloud_fraction_threshold) {
	  v_overlap(jlev,jupper,jlower) = overlap_matrix(jlower,jupper)
	    / frac_upper(jupper);
	}
	else {
	  v_overlap(jlev,jupper,jlower) = 0.0;
	}
      }
    }
    frac_upper = frac_lower;
  } // levels
}
  
}; // namespace tcrad

// Instantiate the direct function
template
void
tcrad::calc_overlap_matrices<false>(int nlev,
				    const Array<2,false>& region_fracs,
				    const Array<1,false>& overlap_param,
				    Array<3,false> u_overlap,
				    Array<3,false> v_overlap,
				    Real cloud_fraction_threshold);
#if ADEPT_REAL_TYPE_SIZE == 8
// Instantiate the differentiable function but only in double
// precision
template
void
tcrad::calc_overlap_matrices<true>(int nlev,
				   const Array<2,true>& region_fracs,
				   const Array<1,true>& overlap_param,
				   Array<3,true> u_overlap,
				   Array<3,true> v_overlap,
				   Real cloud_fraction_threshold);
#endif
