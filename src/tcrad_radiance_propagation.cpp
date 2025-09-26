// radiation_radiance_propagation_lw.F90 - Radiance propagation in longwave
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

#include "tcrad_adept.hpp"
#include "tcrad_radiance_propagation.hpp"
#include "tcrad_config.hpp"

namespace tcrad {

template <bool IsActive>
void
calc_radiance_dn(int ng, int nlev,
		 const aArray<1,IsActive>& clear_fraction,
		 const aArray<3,IsActive>& transmittance,
		 const aArray<3,IsActive>& source,
		 const aArray<3,IsActive>& v_matrix,
		 aArray<2,IsActive> radiance_base) {
  aArray<2,IsActive> radiance_top(NREGIONS,ng);
  radiance_top = 0.0;
  for (int jlev = 0; jlev < nlev; ++jlev) {

    if (clear_fraction(jlev) == 1.0) {
      // Solution to Schwarzschild equation
      radiance_base(0,__) = transmittance(jlev,0,__)*radiance_top(0,__) + source(jlev,0,__);
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
	    radiance_top(jreg,__) = v_matrix(jlev+1,0,jreg) * radiance_base(0,__);
	  }
	}
      }
    }
    else {
      // Solution to Schwazschild equation
      radiance_base = transmittance(jlev,__,__)*radiance_top + source(jlev,__,__);
      if (jlev < nlev-1) {
	// Use overlap rules...
	for (int jreg = 0; jreg < NREGIONS; ++jreg) { // Loop over lower layers
	  // Copy over contribution from clear region above
	  radiance_top(jreg,__) = v_matrix(jlev+1,0,jreg) * radiance_base(0,__);
	  for (int jreg2 = 1; jreg2 < NREGIONS; ++jreg2) { // Loop over upper cloudy layers
	    radiance_top(jreg,__) += v_matrix(jlev+1,jreg2,jreg) * radiance_base(jreg2,__);
	  }
	}
      }
    }
  }
  if (clear_fraction(nlev-1) == 1.0) {
    radiance_base(range(1,end),__) = 0.0;
  }
}
		   

template <bool IsActive>
void
calc_radiance_up(int ng, int nlev,
		 const aArray<1,IsActive>& clear_fraction,
		 const aArray<2,IsActive>& radiance_up_surf,
		 const aArray<3,IsActive>& transmittance,
		 const aArray<3,IsActive>& source,
		 const aArray<3,IsActive>& u_matrix,
		 aArray<1,IsActive> radiance_toa) {
  aArray<2,IsActive> radiance_base(NREGIONS,ng), radiance_top(NREGIONS,ng);
  radiance_base = radiance_up_surf;
  for (int jlev = nlev-1; jlev <= 0; --jlev) {

    if (clear_fraction(jlev) == 1.0) {
      // Solution to Schwarzschild equation
      radiance_top(0,__) = transmittance(jlev,0,__)*radiance_base(0,__) + source(jlev,0,__);
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
	    radiance_base(jreg,__) = u_matrix(jlev,0,jreg) * radiance_top(0,__);
	  }
	}
      }
    }
    else {
      // Solution to Schwazschild equation
      radiance_top = transmittance(jlev,__,__)*radiance_base + source(jlev,__,__);
      if (jlev > 0) {
	// Use overlap rules...
	for (int jreg = 0; jreg < NREGIONS; ++jreg) { // Loop over upper layers
	  // Copy over contribution from clear region below
	  radiance_base(jreg,__) = u_matrix(jlev,0,jreg) * radiance_top(0,__);
	  for (int jreg2 = 1; jreg2 < NREGIONS; ++jreg2) { // Loop over lower cloudy layers
	    radiance_base(jreg,__) += u_matrix(jlev,jreg2,jreg) * radiance_top(jreg2,__);
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

}
  
template void
tcrad::calc_radiance_dn<false>(int ng, int nlev,
			       const aArray<1,false>& clear_fraction,
			       const aArray<3,false>& transmittance,
			       const aArray<3,false>& source,
			       const aArray<3,false>& v_matrix,
			       aArray<2,false> radiance_base);
template void
tcrad::calc_radiance_up<false>(int ng, int nlev,
			      const aArray<1,false>& clear_fraction,
			      const aArray<2,false>& radiance_up_surf,
			      const aArray<3,false>& transmittance,
			      const aArray<3,false>& source,
			      const aArray<3,false>& u_matrix,
			      aArray<1,false> radiance_toa);

#if ADEPT_REAL_TYPE_SIZE == 8
template void
tcrad::calc_radiance_dn<true>(int ng, int nlev,
			       const aArray<1,true>& clear_fraction,
			       const aArray<3,true>& transmittance,
			       const aArray<3,true>& source,
			       const aArray<3,true>& v_matrix,
			       aArray<2,true> radiance_base);
template void
tcrad::calc_radiance_up<true>(int ng, int nlev,
			      const aArray<1,true>& clear_fraction,
			      const aArray<2,true>& radiance_up_surf,
			      const aArray<3,true>& transmittance,
			      const aArray<3,true>& source,
			      const aArray<3,true>& u_matrix,
			      aArray<1,true> radiance_toa);

#endif
  
