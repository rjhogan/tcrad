// radiation_tripleclouds_lw.F90 - Longwave "Tripleclouds" solver
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
#include "tcrad_tripleclouds_lw.hpp"
#include "tcrad_two_stream.hpp"

namespace tcrad {

template <bool IsActive>
void
solver_tripleclouds_lw(int ng,
		       int nlev,
		       int ncol,
		       const Config& config,
		       const aArray<3,IsActive>& region_fracs, // (ncol,nlev,nreg)
		       const aArray<3,IsActive>& od_scaling,   // (ncol,nlev,nreg-1)
		       const aArray<4,IsActive>& u_matrix,     // (ncol,nlev+1,nreg,nreg)
		       const aArray<4,IsActive>& v_matrix,     // (ncol,nlev+1,nreg,nreg)
		       const aArray<3,IsActive>& od,           // (ncol,nlev,ng)
		       const aArray<3,IsActive>& od_cloud,     // (ncol,nlev,ng)
		       const aArray<3,IsActive>& ssa_cloud,    // (ncol,nlev,ng)
		       const aArray<3,IsActive>& asymmetry_cloud, // (ncol,nlev,ng)
		       const aArray<3,IsActive>& planck_hl,    // (ncol,nlev+1,ng)
		       const aMatrix<IsActive>& surf_emission, // (ncol,ng)
		       const aMatrix<IsActive>& surf_albedo,   // (ncol,ng)
		       /*
		       aMatrix<IsActive> lw_up,                // (ncol,nlev+1)
		       aMatrix<IsActive> lw_dn)                // (ncol,nlev+1)
		       */
		       aArray<3,IsActive> lw_up,                // (ncol,nlev+1,ng)
		       aArray<3,IsActive> lw_dn)                // (ncol,nlev+1,ng)
{

  // Optical depth, single scattering albedo and asymmetry factor in
  // each g-point including gas, aerosol and clouds
  aVector<IsActive> od_total(ng), ssa_total(ng), asymmetry_total(ng);

  // Modified optical depth after Tripleclouds scaling to represent
  // cloud inhomogeneity
  aVector<IsActive> od_cloud_new(ng);

  static const int NREGIONS = 3;

  // In a clear-sky layer this will be 1, otherwise equal to NREGIONS
  //  int nreg;

  // Directional overlap matrices defined at all layer interfaces
  // including top-of-atmosphere and the surface
  aArray<3,IsActive> reflectance(nlev,NREGIONS,ng);
  aArray<3,IsActive> transmittance(nlev,NREGIONS,ng);

  // Emission by a layer into the upwelling or downwelling diffuse
  // streams
  aArray<3,IsActive> source_up(nlev,NREGIONS,ng);
  aArray<3,IsActive> source_dn(nlev,NREGIONS,ng);

  // Total albedo of the atmosphere/surface just above a layer
  // interface with respect to downwelling diffuse radiation at that
  // interface, where level index = 1 corresponds to the
  // top-of-atmosphere
  aArray<3,IsActive> total_albedo(nlev+1,NREGIONS,ng);

  // Upwelling radiation just above a layer interface due to emission
  // below that interface, where level index = 1 corresponds to the
  // top-of-atmosphere
  aArray<3,IsActive> total_source(nlev+1,NREGIONS,ng);

  // Total albedo and source of the atmosphere just below a layer interface
  aMatrix<IsActive> total_albedo_below(NREGIONS,ng);
  aMatrix<IsActive> total_source_below(NREGIONS,ng);

  // Downwelling flux below and above an interface between layers into
  // a plane perpendicular to the direction of the sun
  aMatrix<IsActive> flux_dn(NREGIONS,ng);
  aMatrix<IsActive> flux_dn_below(NREGIONS,ng);
  aMatrix<IsActive> flux_up(NREGIONS,ng);

  aMatrix<IsActive> inv_denom(NREGIONS,ng);

  aVector<IsActive> flux_dn_cloud_top(ng);
  
  // Identify clear-sky layers, with pseudo layers for outer space and
  // below the ground, both treated as single-region clear skies
  boolVector is_clear_sky_layer(nlev+2);

  // Index of highest cloudy layer
  int i_cloud_top;

  // --------------------------------------------------------
  // Section 1: Prepare general variables and arrays
  // --------------------------------------------------------

  // Main loop over columns
  for (int jcol = 0; jcol < ncol; ++jcol) {
  
    // --------------------------------------------------------
    // Section 2: Prepare column-specific variables and arrays
    // --------------------------------------------------------
    is_clear_sky_layer = true;
    i_cloud_top = nlev;
    for (int jlev = 0; jlev < nlev; ++jlev) {
      if (region_fracs(jcol,jlev,0) < 1.0-config.cloud_fraction_threshold) {
	is_clear_sky_layer(jlev) = false;
	if (i_cloud_top > jlev) {
	  i_cloud_top = jlev;
	}
      }
    }
    // --------------------------------------------------------
    // Section 3: Clear-sky calculation
    // --------------------------------------------------------
    calc_no_scattering_transmittance_lw(ng, nlev, od(jcol,__,__),
					planck_hl(jcol,range(0,nlev-1),__), planck_hl(jcol,range(1,nlev),__),
					transmittance(__,0,__), source_up(__,0,__), source_dn(__,0,__));

    reflectance(__,0,__) = 0.0;
    
    // --------------------------------------------------------
    // Section 4: Loop over cloudy layers to compute reflectance and transmittance
    // --------------------------------------------------------
    for (int jlev = i_cloud_top; jlev < nlev; ++jlev) { // Start at cloud top and work down
      if (!is_clear_sky_layer(jlev)) {
	// Cloudy sky
	for (int jreg = 1; jreg < NREGIONS; ++jreg) {
	  // Add scaled cloud optical depth to clear-sky value
	  od_cloud_new = od_cloud(jcol,jlev,__) * od_scaling(jcol,jlev,jreg);
	  od_total = od(jcol,jlev,__) + od_cloud_new;
	  ssa_total = 0.0;
	  asymmetry_total = 0.0;
	  // Note that we are neglecting aerosol scattering
	  ssa_total.where(od_total > 0.0) = ssa_cloud(jcol,jlev,__)*od_cloud_new/od_total;
	  asymmetry_total.where(ssa_total > 0.0 && od_total > 0.0)
	    = (asymmetry_cloud(jcol,jlev,__)*ssa_cloud(jcol,jlev,__)*od_cloud_new)
	    / (ssa_total*od_total);
	  calc_ref_trans_lw(ng, od_total, ssa_total, asymmetry_total,
			    planck_hl(jcol,jlev,__), planck_hl(jcol,jlev+1,__),
			    reflectance(jlev,jreg,__), transmittance(jlev,jreg,__),
			    source_up(jlev,jreg,__), source_dn(jlev,jreg,__));
	}
	for (int jreg = 0; jreg < NREGIONS; ++jreg) {
          // Emission is scaled by the size of each region, including
          // scaling down the clear region
	  source_up(jlev,jreg,__) *= region_fracs(jcol,jlev,jreg);
	  source_dn(jlev,jreg,__) *= region_fracs(jcol,jlev,jreg);
	}
      }
    } // Loop over level

    // --------------------------------------------------------
    // Section 5: Compute total sources and albedos at each half level
    // --------------------------------------------------------

    total_albedo = 0.0;
    total_source = 0.0;

    // Calculate the upwelling radiation emitted by the surface, and
    // copy the surface albedo into total_albedo 
    for (int jreg = 0; jreg < NREGIONS; ++jreg) {
      for (int jg = 0; jg < ng; ++jg) {
	total_source(nlev,jreg,jg) = region_fracs(jcol,nlev-1,jreg)*surf_emission(jcol,jg);
	total_albedo(nlev,jreg,jg) = surf_albedo(jcol,jg);
      }
    }
    // Work up from the surface computing the total albedo of the
    // atmosphere and the total upwelling due to emission below each
    // level below using the adding method
    for (int jlev = nlev-1; jlev >= i_cloud_top; --jlev) {
      if (is_clear_sky_layer(jlev)) {
	for (int jg = 0; jg < ng; ++jg) {
	  inv_denom(0,jg) = 1.0 / (1.0 - total_albedo(jlev+1,0,jg)*reflectance(jlev,0,jg));
	  total_albedo_below(0,jg) = reflectance(jlev,0,jg)
	    + transmittance(jlev,0,jg)*transmittance(jlev,0,jg)*total_albedo(jlev+1,0,jg)
	    * inv_denom(0,jg);
	  total_source_below(0,jg) = source_up(jlev,0,jg)
	    + transmittance(jlev,0,jg)*(total_source(jlev+1,0,jg)
					+total_albedo(jlev+1,0,jg)*source_dn(jlev,0,jg))
	    * inv_denom(0,jg);
	}
	total_albedo_below(range(1,end),__) = 0.0;
	total_source_below(range(1,end),__) = 0.0;
      }
      else {
	inv_denom = 1.0 / (1.0 - total_albedo(jlev+1,__,__)*reflectance(jlev,__,__));
	total_albedo_below = reflectance(jlev,__,__)
	  + transmittance(jlev,__,__)*transmittance(jlev,__,__)*total_albedo(jlev+1,__,__)
	  * inv_denom;
	total_source_below = source_up(jlev,__,__)
	  + transmittance(jlev,__,__)*(total_source(jlev+1,__,__)
				       +total_albedo(jlev+1,__,__)*source_dn(jlev,__,__))
	  * inv_denom;
      }

      // Account for cloud overlap when converting albedo below a
      // layer interface to the equivalent values just above
      if (is_clear_sky_layer(jlev) && is_clear_sky_layer(jlev-1)) {
	total_albedo(jlev,__,__) = total_albedo_below;
	total_source(jlev,__,__) = total_source_below;
      }
      else {
	// Loop over regions in upper layer
	for (int jreg = 0; jreg < NREGIONS; ++jreg) {
	  // Loop over regions in lower layer
	  for (int jreg2 = 0; jreg2 < NREGIONS; ++jreg2) {
	    total_source(jlev,jreg,__) += u_matrix(jcol,jlev,jreg2,jreg) * total_source_below(jreg2,__);
	  }
	}
	// Use overlap matrix and exclude "anomalous" horizontal
	// transport described by Shonk & Hogan (2008).  Therefore,
	// the operation we perform is essentially diag(total_albedo)
	// = matmul(v_matrix, diag(total_albedo_below)).
	for (int jreg = 0; jreg < NREGIONS; ++jreg) {
	  for (int jreg2 = 0; jreg2 < NREGIONS; ++jreg2) {
	    total_albedo(jlev,jreg,__)
	      += total_albedo_below(jreg2,__) * v_matrix(jcol,jlev,jreg,jreg2);
	  }
	}
      }

    } // Reverse loop over layers
    /*
    std::cerr << "total_albedo = " << total_albedo(__,__,0) << "\n";
    std::cerr << "total_source = " << total_source(__,__,0) << "\n";
    */
    // --------------------------------------------------------
    // Section 6: Calculate downwelling fluxes from top-of-atmosphere to cloud top
    // --------------------------------------------------------

    // At top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn(0,__) = 0.0; // Index 0 is region not level
    //lw_dn(jcol,0) = 0.0; // Index 0 is level (TOA)
    lw_dn(jcol,0,__) = 0.0; // Index 0 is level (TOA)
    
    // Work down through the atmosphere computing the downward fluxes
    // at each half-level
    for (int jlev = 0; jlev < i_cloud_top; ++jlev) {
      flux_dn(0,__) = transmittance(jlev,0,__)*flux_dn(0,__) + source_dn(jlev,0,__);
      //lw_dn(jcol,jlev+1) = sum(flux_dn(0,__));
      lw_dn(jcol,jlev+1,__) = flux_dn(0,__);
    }

    // Save spectral fluxes at cloud top
    flux_dn_cloud_top = flux_dn(0,__);
    
    // --------------------------------------------------------
    // Section 7: Compute fluxes up to top-of-atmosphere
    // --------------------------------------------------------

    // Compute the fluxes just above the highest cloud
    flux_up(0,__) = total_source(i_cloud_top,0,__)
      + total_albedo(i_cloud_top,0,__)*flux_dn(0,__);
    //lw_up(jcol,i_cloud_top) = sum(flux_up(0,__));
    lw_up(jcol,i_cloud_top,__) = flux_up(0,__);
    for (int jlev = i_cloud_top-1; jlev >= 0; --jlev) {
      flux_up(0,__) = transmittance(jlev,0,__)*flux_up(0,__) + source_up(jlev,0,__);
      //lw_up(jcol,jlev) = sum(flux_up(0,__));
      lw_up(jcol,jlev,__) = flux_up(0,__);
    }

    // --------------------------------------------------------
    // Section 8: Compute fluxes down to surface
    // --------------------------------------------------------

    // Copy over downwelling spectral fluxes at top of first
    // scattering layer, using overlap matrix to translate to the
    // regions of the first layer of cloud
    for (int jreg = 0; jreg < NREGIONS; ++jreg) {
      flux_dn(jreg,__) = v_matrix(jcol,i_cloud_top,0,jreg) * flux_dn_cloud_top;
    }

    // Final loop back down through the atmosphere to compute fluxes
    for (int jlev = i_cloud_top; jlev < nlev; ++jlev) {
      if (is_clear_sky_layer(jlev)) {
	for (int jg = 0; jg < ng; ++jg) {
	  flux_dn(0,jg) = (transmittance(jlev,0,jg) * flux_dn(0,jg)
			   + reflectance(jlev,0,jg) * total_source(jlev+1,0,jg) + source_dn(jlev,0,jg) )
	    / (1.0 - reflectance(jlev,0,jg) * total_albedo(jlev+1,0,jg));
	  flux_up(0,jg) = total_source(jlev+1,0,jg) + flux_dn(0,jg) * total_albedo(jlev+1,0,jg);
	}
	flux_dn(range(1,end),__)  = 0.0;
	flux_up(range(1,end),__)  = 0.0;
      }
      else {
	flux_dn = (transmittance(jlev,__,__) * flux_dn
		   + reflectance(jlev,__,__) * total_source(jlev+1,__,__) + source_dn(jlev,__,__) )
	  / (1.0 - reflectance(jlev,__,__) * total_albedo(jlev+1,__,__));
	flux_up = total_source(jlev+1,__,__) + flux_dn * total_albedo(jlev+1,__,__);
      }

      if (!(is_clear_sky_layer(jlev) && is_clear_sky_layer(jlev+1))) {
	// Account for overlap rules in translating fluxes just above
        // a layer interface to the values just below;
	// Loop over regions in lower layer
	for (int jreg = 0; jreg < NREGIONS; ++jreg) {
	  flux_dn_below(jreg,__) = v_matrix(jcol,jlev+1,0,jreg) * flux_dn(0,__);
	  // Loop over regions in upper layer
	  for (int jreg2 = 1; jreg2 < NREGIONS; ++jreg2) {
	    flux_dn_below(jreg,__) += v_matrix(jcol,jlev+1,jreg2,jreg) * flux_dn(jreg2,__);
	  }
	}
	flux_dn = flux_dn_below;
      }
      lw_dn(jcol,jlev+1,__) = sum(flux_dn,0);
      lw_up(jcol,jlev+1,__) = sum(flux_up,0);
      
      // Otherwise the fluxes in each region are the same so
      // nothing to do
    }
    
  } // Loop over column
  
}

}; // namespace tcrad

// Instantiate the direct function
template void
tcrad::solver_tripleclouds_lw<false>(int ng,
				     int nlev,
				     int ncol,
				     const Config& config,
				     const aArray<3,false>& region_fracs, // (ncol,nlev,nreg)
				     const aArray<3,false>& od_scaling,   // (ncol,nlev,nreg-1)
				     const aArray<4,false>& u_matrix,     // (ncol,nlev+1,nreg,nreg)
				     const aArray<4,false>& v_matrix,     // (ncol,nlev+1,nreg,nreg)
				     const aArray<3,false>& od,           // (ncol,nlev,ng)
				     const aArray<3,false>& od_cloud,     // (ncol,nlev,ng)
				     const aArray<3,false>& ssa_cloud,    // (ncol,nlev,ng)
				     const aArray<3,false>& asymmetry_cloud, // (ncol,nlev,ng)
				     const aArray<3,false>& planck_hl,    // (ncol,nlev+1,ng)
				     const aMatrix<false>& surf_emission, // (ncol,ng)
				     const aMatrix<false>& surf_albedo,   // (ncol,ng)
				     /*
				     aMatrix<false> lw_up,                // (ncol,nlev+1)
				     aMatrix<false> lw_dn);               // (ncol,nlev+1)
				     */
				     aArray<3,false> lw_up,                // (ncol,nlev+1)
				     aArray<3,false> lw_dn);               // (ncol,nlev+1)

#if ADEPT_REAL_TYPE_SIZE == 8
// Instantiate the differentiable function but only in double
// precision
template void
tcrad::solver_tripleclouds_lw<true>(int ng,
				    int nlev,
				    int ncol,
				    const Config& config,
				    const aArray<3,true>& region_fracs, // (ncol,nlev,nreg)
				    const aArray<3,true>& od_scaling,   // (ncol,nlev,nreg-1)
				    const aArray<4,true>& u_matrix,     // (ncol,nlev+1,nreg,nreg)
				    const aArray<4,true>& v_matrix,     // (ncol,nlev+1,nreg,nreg)
				    const aArray<3,true>& od,           // (ncol,nlev,ng)
				    const aArray<3,true>& od_cloud,     // (ncol,nlev,ng)
				    const aArray<3,true>& ssa_cloud,    // (ncol,nlev,ng)
				    const aArray<3,true>& asymmetry_cloud, // (ncol,nlev,ng)
				    const aArray<3,true>& planck_hl,    // (ncol,nlev+1,ng)
				    const aMatrix<true>& surf_emission, // (ncol,ng)
				    const aMatrix<true>& surf_albedo,   // (ncol,ng)
				    /*
				    aMatrix<true> lw_up,                // (ncol,nlev+1)
				    aMatrix<true> lw_dn);               // (ncol,nlev+1)
				    */
				    aArray<3,true> lw_up,                // (ncol,nlev+1)
				    aArray<3,true> lw_dn);               // (ncol,nlev+1)


#endif
