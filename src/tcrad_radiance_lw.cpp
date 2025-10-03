// tcrad_radiance_lw.cpp - Compute longwave radiances
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

#include "tcrad_radiance_lw.hpp"
#include "tcrad_regions.hpp"
#include "tcrad_overlap.hpp"
#include "tcrad_tripleclouds_lw.hpp"
#include "tcrad_propagate_radiance_lw.hpp"
#include "tcrad_layer_solutions_lw.hpp"

namespace tcrad {
  
// ---------------------------------------------------------------------
// Longwave "Tripleclouds" solver following Shonk & Hogan (2008) that
// treats cloud inhomogeneity by dividing each model level into three
// regions, one clear and two cloudy (with differing optical depth),
// returning a top-of-atmosphere spectral upwelling radiances.
template <bool IsActive>
void calc_tripleclouds_radiance_lw(int ng,
				   int nlev,
				   int ncol,
				   const Config& config,
				   const Array<1>& mu, // Cosine of sensor zenith angle
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
				   Array<2,IsActive> radiance)
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

  // Transmittance along the path to the satellite
  Array<3,IsActive> transmittance(nlev,NREGIONS,ng);
  Array<3,IsActive> source_dn(nlev,NREGIONS,ng);
  Array<3,IsActive> source_up(nlev,NREGIONS,ng);

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
    solver_tripleclouds_lw(ng, nlev, config, region_fracs[jcol],
			   u_overlap[jcol], v_overlap[jcol],
			   od, ssa, asymmetry_cloud[jcol],
			   planck_hl[jcol], surf_emission[jcol], surf_albedo[jcol],
			   flux_up_base, flux_dn_base, flux_up_top, flux_dn_top);
    Array<2,IsActive> radiance_up_surf(NREGIONS,ng);
    calc_radiance_trans_source_lw(ng, nlev, config, mu(jcol), region_fracs[jcol],
				  planck_hl[jcol], od, ssa, asymmetry_cloud[jcol],
				  flux_up_base, flux_dn_top,
				  transmittance, source_up, source_dn);
    if (config.do_specular_surface) {
      Array<2,IsActive> radiance_dn_surf(NREGIONS,ng);
      propagate_radiance_dn_lw(ng, nlev, region_fracs(jcol,__,0), transmittance,
			       source_dn, v_overlap[jcol], radiance_dn_surf);
      for (int jreg = 0; jreg < NREGIONS; ++jreg) {
	radiance_up_surf(jreg,__) = region_fracs(jcol,nlev-1,jreg)
	  *surf_emission[jcol]/tcrad::PI
	  + surf_albedo[jcol]*radiance_dn_surf(jreg,__);
      }
    }
    else {
      radiance_up_surf = flux_up_base(nlev-1,__,__)/tcrad::PI;
    }
    propagate_radiance_up_lw(ng, nlev, region_fracs(jcol,__,0), radiance_up_surf,
			     transmittance, source_up,
			     u_overlap[jcol], radiance[jcol]);
  }
}

};

template void
tcrad::calc_tripleclouds_radiance_lw<false>(int ng,
					    int nlev,
					    int ncol,
					    const Config& config,
					    const Array<1>& mu, // Cosine of sensor zenith angle
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
					    Array<2,false> radiance);

#if ADEPT_REAL_TYPE_SIZE == 8
template void
tcrad::calc_tripleclouds_radiance_lw<true>(int ng,
					   int nlev,
					   int ncol,
					   const Config& config,
					   const Array<1>& mu, // Cosine of sensor zenith angle
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
					   Array<2,true> radiance);

#endif

  



