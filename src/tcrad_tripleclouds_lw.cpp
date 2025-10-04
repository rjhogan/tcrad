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
//
// This file is adapted from ecRad's radiation_tripleclouds_lw.F90

#include "tcrad_tripleclouds_lw.hpp"
#include "tcrad_layer_solutions_lw.hpp"

namespace tcrad {

// ---------------------------------------------------------------------
// Longwave "Tripleclouds" solver following Shonk & Hogan (2008) that
// treats cloud inhomogeneity by dividing each model level into three
// regions, one clear and two cloudy (with differing optical
// depth). This version computes the upwelling and downwelling fluxes
// at the top and bottom of each layer and region, and can used in
// subsequent processing for computing either broadband flux profiles
// or top-of-atmosphere radiances.
template <bool IsActive>
void
solver_tripleclouds_lw(int ng,
		       int nlev,
		       const Config& config,
		       const Array<2,IsActive>& region_fracs, // (nlev,nreg)
		       const Array<3,IsActive>& u_overlap,    // (nlev+1,nreg,nreg)
		       const Array<3,IsActive>& v_overlap,    // (nlev+1,nreg,nreg)
		       const Array<3,IsActive>& od,           // (nlev,nreg,ng)
		       const Array<3,IsActive>& ssa,          // (nlev,nreg,ng)
		       const Array<2,IsActive>& asymmetry,    // (nlev,ng)
		       const Array<2,IsActive>& planck_hl,    // (nlev+1,ng)
		       const Array<1,IsActive>& surf_emission,// (ng)
		       const Array<1,IsActive>& surf_albedo,  // (ng)
		       Array<3,IsActive> flux_up_base,        // (nlev,nreg,ng)
		       Array<3,IsActive> flux_dn_base,        // (nlev,nreg,ng)
		       Array<3,IsActive> flux_up_top,         // (nlev,nreg,ng)
		       Array<3,IsActive> flux_dn_top)         // (nlev,nreg,ng)
{
  // Layer reflectance and transmittance in layers and regions but
  // note that since reflectance is zero in the clear region its
  // region dimension is one less and all indexing is jreg-1.
  Array<3,IsActive> reflectance(nlev,NREGIONS-1,ng);
  Array<3,IsActive> transmittance(nlev,NREGIONS,ng);

  // Emission by a layer into the upwelling or downwelling diffuse
  // streams
  Array<3,IsActive> source_up(nlev,NREGIONS,ng);
  Array<3,IsActive> source_dn(nlev,NREGIONS,ng);

  // Total albedo of the atmosphere/surface just above a layer
  // interface with respect to downwelling diffuse radiation at that
  // interface, where level index = 1 corresponds to the
  // top-of-atmosphere
  Array<3,IsActive> total_albedo(nlev+1,NREGIONS,ng);

  // Upwelling radiation just above a layer interface due to emission
  // below that interface, where level index = 1 corresponds to the
  // top-of-atmosphere
  Array<3,IsActive> total_source(nlev+1,NREGIONS,ng);

  // Total albedo and source of the atmosphere just below a layer interface
  Array<2,IsActive> total_albedo_below(NREGIONS,ng);
  Array<2,IsActive> total_source_below(NREGIONS,ng);

  // Inverse of the denominator in the geometric series representing
  // multiple reflections between layers
  Array<1,IsActive> inv_denom(ng);

  // Identify clear-sky layers, with pseudo layer for below the
  // ground, treated as a single clear-sky region
  boolVector is_clear_sky_layer(nlev+1);

  // Index of highest cloudy layer
  int i_cloud_top;

  // --------------------------------------------------------
  // Section 1: Prepare general variables and arrays (redundant)
  // --------------------------------------------------------

  // --------------------------------------------------------
  // Section 2: Prepare column-specific variables and arrays
  // --------------------------------------------------------
  is_clear_sky_layer = true;
  i_cloud_top = nlev; // If no cloud this is below the surface!
  for (int jlev = 0; jlev < nlev; ++jlev) {
    if (region_fracs(jlev,0) < 1.0-config.cloud_fraction_threshold) {
      is_clear_sky_layer(jlev) = false;
      if (i_cloud_top > jlev) {
	i_cloud_top = jlev;
      }
    }
  }
  
  // --------------------------------------------------------
  // Section 3: Clear-sky calculation
  // --------------------------------------------------------

  // Do we compute the properties of all layers in multiple function
  // calls or a single one? The latter is faster.
#define ALL_LAYERS 1
#ifndef ALL_LAYERS
  calc_no_scat_trans_source_lw(ng, nlev, config, od(__,0,__),
			       planck_hl(range(0,nlev-1),__), planck_hl(range(1,nlev),__),
			       transmittance(__,0,__), source_up(__,0,__), source_dn(__,0,__));
  for (int jlev = i_cloud_top; jlev < nlev; ++jlev) { // Start at cloud top and work down
    if (!is_clear_sky_layer(jlev)) {
      // Cloudy sky
      for (int jreg = 1; jreg < NREGIONS; ++jreg) {
	calc_ref_trans_source_layer_lw(ng, config, od(jlev,jreg,__), ssa(jlev,jreg,__), asymmetry(jlev,__),
				 planck_hl(jlev,__), planck_hl(jlev+1,__),
				 reflectance(jlev,jreg-1,__), transmittance(jlev,jreg,__),
				 source_up(jlev,jreg,__), source_dn(jlev,jreg,__));
      }
      for (int jreg = 0; jreg < NREGIONS; ++jreg) {
	// Emission is scaled by the size of each region, including
	// scaling down the clear region
	source_up(jlev,jreg,__) *= region_fracs(jlev,jreg);
	source_dn(jlev,jreg,__) *= region_fracs(jlev,jreg);
      }
    }
  }
#else
  // --------------------------------------------------------
  // Section 4: Loop over cloudy layers to compute reflectance and transmittance
  // --------------------------------------------------------
  calc_ref_trans_source_lw(ng, nlev, config, region_fracs, od, ssa, asymmetry,
			   planck_hl(range(0,nlev-1),__), planck_hl(range(1,nlev),__),
			   reflectance, transmittance,
			   source_up, source_dn);
#endif
  
  // --------------------------------------------------------
  // Section 5: Compute total sources and albedos at each half level
  // --------------------------------------------------------

  total_albedo = 0.0;
  total_source = 0.0;

  // Calculate the upwelling radiation emitted by the surface, and
  // copy the surface albedo into total_albedo 
  for (int jreg = 0; jreg < NREGIONS; ++jreg) {
    for (int jg = 0; jg < ng; ++jg) {
      total_source(nlev,jreg,jg) = region_fracs(nlev-1,jreg)*surf_emission(jg);
      total_albedo(nlev,jreg,jg) = surf_albedo(jg);
    }
  }
  // Work up from the surface computing the total albedo of the
  // atmosphere and the total upwelling due to emission below each
  // level below using the adding method
  for (int jlev = nlev-1; jlev >= i_cloud_top; --jlev) {
    // Clear region: no reflection from the layer itself so simpler equations
    for (int jg = 0; jg < ng; ++jg) {
      total_albedo_below(0,jg) =
	transmittance(jlev,0,jg)*transmittance(jlev,0,jg)*total_albedo(jlev+1,0,jg);
      total_source_below(0,jg) = source_up(jlev,0,jg)
	+ transmittance(jlev,0,jg)*(total_source(jlev+1,0,jg)
				    +total_albedo(jlev+1,0,jg)*source_dn(jlev,0,jg));
    }
    if (is_clear_sky_layer(jlev)) {
      // Fill unused cloudy regions with zeros
      total_albedo_below(range(1,end),__) = 0.0;
      total_source_below(range(1,end),__) = 0.0;
    }
    else {
      for (int jreg = 1; jreg < NREGIONS; ++jreg) {
	inv_denom = 1.0 / (1.0 - total_albedo(jlev+1,jreg,__)*reflectance(jlev,jreg-1,__));
	total_albedo_below(jreg,__) = reflectance(jlev,jreg-1,__)
	  + transmittance(jlev,jreg,__)*transmittance(jlev,jreg,__)*total_albedo(jlev+1,jreg,__)
	  * inv_denom;
	total_source_below(jreg,__) = source_up(jlev,jreg,__)
	  + transmittance(jlev,jreg,__)*(total_source(jlev+1,jreg,__)
					 +total_albedo(jlev+1,jreg,__)*source_dn(jlev,jreg,__))
	  * inv_denom;
      }
    }

    // Account for cloud overlap when converting albedo below a
    // layer interface to the equivalent values just above
    bool is_clear_above = true;
    if (jlev > 0) {
      is_clear_above = is_clear_sky_layer(jlev-1);
    }
    if (is_clear_sky_layer(jlev) && is_clear_above) {
      total_albedo(jlev,__,__) = total_albedo_below;
      total_source(jlev,__,__) = total_source_below;
    }
    else {
      // Loop over regions in upper layer
      for (int jreg = 0; jreg < NREGIONS; ++jreg) {
	// Loop over regions in lower layer
	for (int jreg2 = 0; jreg2 < NREGIONS; ++jreg2) {
	  total_source(jlev,jreg,__) += u_overlap(jlev,jreg2,jreg) * total_source_below(jreg2,__);
	}
      }
      // Use overlap matrix and exclude "anomalous" horizontal
      // transport described by Shonk & Hogan (2008).  Therefore,
      // the operation we perform is essentially diag(total_albedo)
      // = matmul(v_overlap, diag(total_albedo_below)).
      for (int jreg = 0; jreg < NREGIONS; ++jreg) {
	for (int jreg2 = 0; jreg2 < NREGIONS; ++jreg2) {
	  total_albedo(jlev,jreg,__)
	    += total_albedo_below(jreg2,__) * v_overlap(jlev,jreg,jreg2);
	}
      }
    }

  } // Reverse loop over layers

  
  // --------------------------------------------------------
  // Section 6: Calculate downwelling fluxes from top-of-atmosphere to cloud top
  // --------------------------------------------------------

  // Clear the output flux variables
  flux_up_base = 0.0;
  flux_dn_base = 0.0;
  flux_up_top  = 0.0;
  flux_dn_top  = 0.0;

  // At top-of-atmosphere there is no diffuse downwelling radiation
    
  // Work down through the atmosphere computing the downward fluxes
  // at each half-level
  for (int jlev = 0; jlev < i_cloud_top; ++jlev) {
    flux_dn_base(jlev,0,__) = transmittance(jlev,0,__)*flux_dn_top(jlev,0,__) + source_dn(jlev,0,__);
    // In a completely clear sky, i_cloud_top is not a valid layer
    if (jlev < nlev-1) {
      flux_dn_top(jlev+1,0,__) = flux_dn_base(jlev,0,__);
    }
  }

  // --------------------------------------------------------
  // Section 7: Compute fluxes up to top-of-atmosphere
  // --------------------------------------------------------

  if (i_cloud_top > 0) {
    // Compute the fluxes just above the highest cloud
    flux_up_base(i_cloud_top-1,0,__) = total_source(i_cloud_top,0,__)
      + total_albedo(i_cloud_top,0,__)*flux_dn_base(i_cloud_top-1,0,__);
    flux_up_top(i_cloud_top-1,0,__) = source_up(i_cloud_top-1,0,__)
      + transmittance(i_cloud_top-1,0,__)*flux_up_base(i_cloud_top-1,0,__);
    for (int jlev = i_cloud_top-2; jlev >= 0; --jlev) {
      flux_up_base(jlev,0,__) = flux_up_top(jlev+1,0,__);
      flux_up_top(jlev,0,__)  = source_up(jlev,0,__)
	+ transmittance(jlev,0,__)*flux_up_base(jlev,0,__);
    }
  }

  // --------------------------------------------------------
  // Section 8: Compute fluxes down to surface
  // --------------------------------------------------------

  // Copy over downwelling spectral fluxes at top of first
  // scattering layer, using overlap matrix to translate to the
  // regions of the first layer of cloud
  if (i_cloud_top > 0 && i_cloud_top < nlev) {
    for (int jreg = 0; jreg < NREGIONS; ++jreg) {
      flux_dn_top(i_cloud_top,jreg,__)
	= v_overlap(i_cloud_top,0,jreg) * flux_dn_base(i_cloud_top-1,0,__);
    }
    // else the highest layer is cloudy, in which case
    // flux_dn_top(jlev=0,__,__) is already set to zero
  }

  // Final loop back down through the atmosphere to compute fluxes
  for (int jlev = i_cloud_top; jlev < nlev; ++jlev) {
    // Clear region
    for (int jg = 0; jg < ng; ++jg) {
      flux_dn_base(jlev,0,jg) = transmittance(jlev,0,jg) * flux_dn_top(jlev,0,jg)
	+ source_dn(jlev,0,jg);
      flux_up_base(jlev,0,jg) = total_source(jlev+1,0,jg)
	+ flux_dn_base(jlev,0,jg)*total_albedo(jlev+1,0,jg);
      flux_up_top(jlev,0,jg) = flux_up_base(jlev,0,jg) * transmittance(jlev,0,jg)
	+ source_up(jlev,0,jg);
    }
    if (!is_clear_sky_layer(jlev)) {
      for (int jreg = 1; jreg < NREGIONS; ++jreg) {
	flux_dn_base(jlev,jreg,__) = (transmittance(jlev,jreg,__) * flux_dn_top(jlev,jreg,__)
	      + reflectance(jlev,jreg-1,__) * total_source(jlev+1,jreg,__) + source_dn(jlev,jreg,__) )
	  / (1.0 - reflectance(jlev,jreg-1,__) * total_albedo(jlev+1,jreg,__));
	flux_up_base(jlev,jreg,__) = total_source(jlev+1,jreg,__)
	  + flux_dn_base(jlev,jreg,__) * total_albedo(jlev+1,jreg,__);
	flux_up_top(jlev,jreg,__) = flux_up_base(jlev,jreg,__) * transmittance(jlev,jreg,__)
	  + source_up(jlev,jreg,__) + flux_dn_top(jlev,jreg,__) * reflectance(jlev,jreg-1,__);
      }
    }

    if (jlev < nlev-1) {
      if (!(is_clear_sky_layer(jlev) && is_clear_sky_layer(jlev+1))) {
	// Account for overlap rules in translating fluxes just above
	// a layer interface to the values just below;
	// Loop over regions in lower layer
	for (int jreg = 0; jreg < NREGIONS; ++jreg) {
	  flux_dn_top(jlev+1,jreg,__) = v_overlap(jlev+1,0,jreg) * flux_dn_base(jlev,0,__);
	  // Loop over regions in upper layer
	  for (int jreg2 = 1; jreg2 < NREGIONS; ++jreg2) {
	    flux_dn_top(jlev+1,jreg,__) += v_overlap(jlev+1,jreg2,jreg) * flux_dn_base(jlev,jreg2,__);
	  }
	}
      }
      else {
	flux_dn_top(jlev+1,0,__) = flux_dn_base(jlev,0,__);
      }
    }
      
    // Otherwise the fluxes in each region are the same so
    // nothing to do
  }
  
}

}; // namespace tcrad

// Instantiate the direct function
template void
tcrad::solver_tripleclouds_lw<false>(int ng,
				     int nlev,
				     const Config& config,
				     const Array<2,false>& region_fracs, // (nlev,nreg)
				     const Array<3,false>& u_overlap,    // (nlev+1,nreg,nreg)
				     const Array<3,false>& v_overlap,    // (nlev+1,nreg,nreg)
				     const Array<3,false>& od,           // (nlev,nreg,ng)
				     const Array<3,false>& ssa,          // (nlev,nreg,ng)
				     const Array<2,false>& asymmetry,    // (nlev,ng)
				     const Array<2,false>& planck_hl,    // (nlev+1,ng)
				     const Array<1,false>& surf_emission,// (ng)
				     const Array<1,false>& surf_albedo,  // (ng)
				     Array<3,false> flux_up_base,        // (nlev,nreg,ng)
				     Array<3,false> flux_dn_base,        // (nlev,nreg,ng)
				     Array<3,false> flux_up_top,         // (nlev,nreg,ng)
				     Array<3,false> flux_dn_top);        // (nlev,nreg,ng)

#if ADEPT_REAL_TYPE_SIZE == 8
// Instantiate the differentiable function but only in double
// precision
template void
tcrad::solver_tripleclouds_lw<true>(int ng,
				    int nlev,
				    const Config& config,
				    const Array<2,true>& region_fracs, // (nlev,nreg)
				    const Array<3,true>& u_overlap,    // (nlev+1,nreg,nreg)
				    const Array<3,true>& v_overlap,    // (nlev+1,nreg,nreg)
				    const Array<3,true>& od,           // (nlev,nreg,ng)
				    const Array<3,true>& ssa,          // (nlev,nreg,ng)
				    const Array<2,true>& asymmetry,    // (nlev,ng)
				    const Array<2,true>& planck_hl,    // (nlev+1,ng)
				    const Array<1,true>& surf_emission,// (ng)
				    const Array<1,true>& surf_albedo,  // (ng)
				    Array<3,true> flux_up_base,        // (nlev,nreg,ng)
				    Array<3,true> flux_dn_base,        // (nlev,nreg,ng)
				    Array<3,true> flux_up_top,         // (nlev,nreg,ng)
				    Array<3,true> flux_dn_top);        // (nlev,nreg,ng)

#endif
