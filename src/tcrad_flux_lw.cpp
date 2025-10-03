// tcrad_flux_lw.cpp - Compute longwave fluxes
//
// (C) Copyright 2025- ECMWF.
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

#include "tcrad_flux_lw.hpp"
#include "tcrad_regions.hpp"
#include "tcrad_overlap.hpp"
#include "tcrad_tripleclouds_lw.hpp"
#include "tcrad_layer_solutions_lw.hpp"

namespace tcrad {

// ---------------------------------------------------------------------
// Longwave "Tripleclouds" solver following Shonk & Hogan (2008) that
// treats cloud inhomogeneity by dividing each model level into three
// regions, one clear and two cloudy (with differing optical depth),
// returning a profile of upwelling and downwelling spectral fluxes.
template <bool IsActive>
void calc_tripleclouds_flux_lw(int ng,
			       int nlev,
			       int ncol,
			       const Config& config,
			       const Array<2,IsActive>& surf_emission,
			       const Array<2,IsActive>& surf_albedo,
			       const Array<3,IsActive>& planck_hl,
			       const Array<2,IsActive>& cloud_fraction,
			       const Array<2>& fractional_std,
			       const Array<3,IsActive>& od_clear,
			       const Array<3,IsActive>& od_cloud,
			       const Array<3,IsActive>& ssa_cloud,
			       const Array<3,IsActive>& asymmetry_cloud,
			       const Array<2,IsActive>& overlap_param,
			       Array<3,IsActive> flux_up,
			       Array<3,IsActive> flux_dn)
{

  Array<3,IsActive> region_fracs(ncol,nlev,NREGIONS);
  Array<3,IsActive> od_scaling(ncol,nlev,NREGIONS);
  calc_region_properties(nlev, NREGIONS, ncol, cloud_fraction, fractional_std,
			 region_fracs, od_scaling, config.cloud_fraction_threshold);

  Array<4,IsActive> u_overlap(ncol,nlev+1,NREGIONS,NREGIONS);
  Array<4,IsActive> v_overlap(ncol,nlev+1,NREGIONS,NREGIONS);
  calc_overlap_matrices(nlev, ncol, region_fracs, overlap_param,
			u_overlap, v_overlap, config.cloud_fraction_threshold);

  Array<3,IsActive> flux_up_base(nlev,NREGIONS,ng);
  Array<3,IsActive> flux_dn_base(nlev,NREGIONS,ng);
  Array<3,IsActive> flux_up_top (nlev,NREGIONS,ng);
  Array<3,IsActive> flux_dn_top (nlev,NREGIONS,ng);

  // Merged clear+cloudy optical properties, noting that the asymmetry
  // factor in the cloudy region equals that of the cloud since there
  // is no air/aerosol scattering
  Array<3,IsActive> od(nlev,NREGIONS,ng);
  Array<3,IsActive> ssa(nlev,NREGIONS,ng);

  // Loop over columns
  for (int jcol = 0; jcol < ncol; ++jcol) {
    
    od(__,0,__) = od_clear[jcol];
    for (int jlev = 0; jlev < nlev; ++jlev) {
      if (region_fracs(jcol,jlev,0) < 1.0) {
	// Cloud present in layer: compute combined air+cloud
	// scattering properties
	for (int jreg = 1; jreg < NREGIONS; ++jreg) {
	  od(jlev,jreg,__) = od_clear(jcol,jlev,__)
	    + od_scaling(jcol,jlev,jreg) * od_cloud(jcol,jlev,__);
	  ssa(jlev,jreg,__) = ssa_cloud(jcol,jlev,__)
	    * od_cloud(jcol,jlev,__) * od_scaling(jcol,jlev,jreg)
	    / od(jlev,jreg,__);
	}
      }
    }
    
    // Calculate fluxes
    solver_tripleclouds_lw(ng, nlev, config, region_fracs[jcol],
			   u_overlap[jcol], v_overlap[jcol],
			   od, ssa, asymmetry_cloud[jcol],
			   planck_hl[jcol], surf_emission[jcol], surf_albedo[jcol],
			   flux_up_base, flux_dn_base, flux_up_top, flux_dn_top);

    // Sum over regions and copy to output
    flux_up(jcol,range(0,nlev-1),__) = sum(flux_up_top,1);
    flux_up(jcol,nlev,__)            = sum(flux_up_base(nlev-1,__,__),0);
    flux_dn(jcol,range(0,nlev-1),__) = sum(flux_dn_top,1);
    flux_dn(jcol,nlev,__)            = sum(flux_dn_base(nlev-1,__,__),0);
  }
}

};

template void
tcrad::calc_tripleclouds_flux_lw<false>(int ng,
					int nlev,
					int ncol,
					const Config& config,
					const Array<2,false>& surf_emission,
					const Array<2,false>& surf_albedo,
					const Array<3,false>& planck_hl,
					const Array<2,false>& cloud_fraction,
					const Array<2>& fractional_std,
					const Array<3,false>& od_clear,
					const Array<3,false>& od_cloud,
					const Array<3,false>& ssa_cloud,
					const Array<3,false>& asymmetry_cloud,
					const Array<2,false>& overlap_param,
					Array<3,false> flux_up,
					Array<3,false> flux_dn);

#if ADEPT_REAL_TYPE_SIZE == 8
template void
tcrad::calc_tripleclouds_flux_lw<true>(int ng,
				       int nlev,
				       int ncol,
				       const Config& config,
				       const Array<2,true>& surf_emission,
				       const Array<2,true>& surf_albedo,
				       const Array<3,true>& planck_hl,
				       const Array<2,true>& cloud_fraction,
				       const Array<2>& fractional_std,
				       const Array<3,true>& od_clear,
				       const Array<3,true>& od_cloud,
				       const Array<3,true>& ssa_cloud,
				       const Array<3,true>& asymmetry_cloud,
				       const Array<2,true>& overlap_param,
				       Array<3,true> flux_up,
				       Array<3,true> flux_dn);

#endif

  



