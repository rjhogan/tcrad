// radiation_radiance_propagation_lw.cpp - Radiance propagation in longwave
//
// (C) Copyright 2021- ECMWF.
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

#include "tcrad_propagate_radiance_lw.hpp"
#include "tcrad_config.hpp"

namespace tcrad {

// ---------------------------------------------------------------------
// Compute downward radiance profile by solving the Schwarzschild
// radiative transfer equation in each region of each layer assuming
// the source term to vary linearly with optical depth in each
// layer. The overlap matrix is used to translate the radiances exiting
// the base of one layer into radiances in the top of the layer below,
// consistent with the Tripleclouds approximation. 
template <bool IsActive>
void
propagate_radiance_dn_lw(int ng, int nlev,
			 const Array<1,IsActive>& clear_fraction,
			 const Array<3,IsActive>& transmittance,
			 const Array<3,IsActive>& source,
			 const Array<3,IsActive>& v_overlap,
			 Array<2,IsActive> radiance_base) {

  Array<2,IsActive> radiance_top(NREGIONS,ng);
  radiance_top = 0.0;
  
  for (int jlev = 0; jlev < nlev; ++jlev) {
    if (clear_fraction(jlev) == 1.0) {
      // Solution to Schwarzschild equation
      radiance_base(0,__) = transmittance(jlev,0,__)*radiance_top(0,__)
	+ source(jlev,0,__);
      if (jlev < nlev-1) {
	if (clear_fraction(jlev+1) == 1.0) {
	  // Copy over clear-sky radiances
	  radiance_top(0,__) = radiance_base(0,__);
	}
	else {
	  // Clear above, cloudy below; use overlap rules to
	  // distribute radiance from base of clear region above to
	  // top of all three regions below; loop over cloudy regions
	  // of lower layer
	  for (int jreg = 0; jreg < NREGIONS; ++jreg) {
	    radiance_top(jreg,__) = v_overlap(jlev+1,0,jreg)
	      * radiance_base(0,__);
	  }
	}
      }
    }
    else {
      // Solution to Schwazschild equation
      radiance_base = transmittance(jlev,__,__)*radiance_top
	+ source(jlev,__,__);
      if (jlev < nlev-1) {
	// Use overlap rules...
	// Loop over lower layers
	for (int jreg = 0; jreg < NREGIONS; ++jreg) {
	  // Copy over contribution from clear region above
	  radiance_top(jreg,__) = v_overlap(jlev+1,0,jreg)
	    * radiance_base(0,__);
	  // Loop over upper cloudy layers
	  for (int jreg2 = 1; jreg2 < NREGIONS; ++jreg2) {
	    radiance_top(jreg,__) += v_overlap(jlev+1,jreg2,jreg)
	      * radiance_base(jreg2,__);
	  }
	}
      }
    }
  }
  if (clear_fraction(nlev-1) == 1.0) {
    radiance_base(range(1,end),__) = 0.0;
  }
}
		   

// ---------------------------------------------------------------------
// Compute upward radiance profile by solving the Schwarzschild
// radiative transfer equation in each region of each layer assuming
// the source term to vary linearly with optical depth in each
// layer. The overlap matrix is used to translate the radiances exiting
// the top of one layer into radiances in the base of the layer above,
// consistent with the Tripleclouds approximation. 
template <bool IsActive>
void
propagate_radiance_up_lw(int ng, int nlev,
			 const Array<1,IsActive>& clear_fraction,
			 const Array<2,IsActive>& radiance_up_surf,
			 const Array<3,IsActive>& transmittance,
			 const Array<3,IsActive>& source,
			 const Array<3,IsActive>& u_overlap,
			 Array<1,IsActive> radiance_toa) {
  Array<2,IsActive> radiance_base(NREGIONS,ng), radiance_top(NREGIONS,ng);
  radiance_base = radiance_up_surf;
  
  for (int jlev = nlev-1; jlev >= 0; --jlev) {
    if (clear_fraction(jlev) == 1.0) {
      // Solution to Schwarzschild equation
      radiance_top(0,__) = transmittance(jlev,0,__)*radiance_base(0,__)
	+ source(jlev,0,__);
      if (jlev > 0) {
	if (clear_fraction(jlev-1) == 1.0) {
	  // Copy over clear-sky radiances
	  radiance_base(0,__) = radiance_top(0,__);
	}
	else {
	  // Clear below, cloudy above; use overlap rules to
	  // distribute radiance from top of clear region below to
	  // base of all three regions above; loop over cloudy regions
	  // of upper layer
	  for (int jreg = 0; jreg < NREGIONS; ++jreg) {
	    radiance_base(jreg,__) = u_overlap(jlev,0,jreg)
	      * radiance_top(0,__);
	  }
	}
      }
    }
    else {
      // Solution to Schwazschild equation
      radiance_top = transmittance(jlev,__,__)*radiance_base
	+ source(jlev,__,__);
      if (jlev > 0) {
	// Use overlap rules...
	// Loop over upper layers
	for (int jreg = 0; jreg < NREGIONS; ++jreg) {
	  // Copy over contribution from clear region below
	  radiance_base(jreg,__) = u_overlap(jlev,0,jreg) * radiance_top(0,__);
	  // Loop over lower cloudy layers
	  for (int jreg2 = 1; jreg2 < NREGIONS; ++jreg2) {
	    radiance_base(jreg,__) += u_overlap(jlev,jreg2,jreg)
	      * radiance_top(jreg2,__);
	  }
	}
      }
    }
  }
  if (clear_fraction(0) == 1.0) {
    radiance_toa = radiance_top(0,__);
  }
  else {
    radiance_toa = sum(radiance_top,0);
  }
}

};

// Instantiate the direct functions
template void
tcrad::propagate_radiance_dn_lw<false>(int ng, int nlev,
				       const Array<1,false>& clear_fraction,
				       const Array<3,false>& transmittance,
				       const Array<3,false>& source,
				       const Array<3,false>& v_overlap,
				       Array<2,false> radiance_base);
template void
tcrad::propagate_radiance_up_lw<false>(int ng, int nlev,
				       const Array<1,false>& clear_fraction,
				       const Array<2,false>& radiance_up_surf,
				       const Array<3,false>& transmittance,
				       const Array<3,false>& source,
				       const Array<3,false>& u_overlap,
				       Array<1,false> radiance_toa);


#if ADEPT_REAL_TYPE_SIZE == 8
// Instantiate the differentiable functions but only in double
// precision
template void
tcrad::propagate_radiance_dn_lw<true>(int ng, int nlev,
				      const Array<1,true>& clear_fraction,
				      const Array<3,true>& transmittance,
				      const Array<3,true>& source,
				      const Array<3,true>& v_overlap,
				      Array<2,true> radiance_base);
template void
tcrad::propagate_radiance_up_lw<true>(int ng, int nlev,
				      const Array<1,true>& clear_fraction,
				      const Array<2,true>& radiance_up_surf,
				      const Array<3,true>& transmittance,
				      const Array<3,true>& source,
				      const Array<3,true>& u_overlap,
				      Array<1,true> radiance_toa);

#endif
  
