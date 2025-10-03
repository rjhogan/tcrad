// tcrad_layer_solutions_lw.cpp - Solutions to longwave two-stream equations in a layer
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

#ifndef TCRAD_LAYER_SOLUTIONS_LW_H
#define TCRAD_LAYER_SOLUTIONS_LW_H 1

#include "tcrad_adept.hpp"
#include "tcrad_config.hpp"

namespace tcrad {

// ---------------------------------------------------------------------
// Compute the longwave transmittance to diffuse radiation in the
// no-scattering (clear-sky) case, as well as the upward flux at the
// top and the downward flux at the base of the layer due to emission
// from within the layer assuming a linear variation of Planck
// function within the layer.
template<bool IsActive>
void
calc_no_scat_trans_source_lw(int ng,   // Number of spectral intervals
			     int nlev, // Number of levels
			     const Config& config,
			     const Array<2,IsActive>& od,
			     const Array<2,IsActive>& planck_top,
			     const Array<2,IsActive>& planck_base,
			     Array<2,IsActive> transmittance,
			     Array<2,IsActive> source_up,
			     Array<2,IsActive> source_dn);

// ---------------------------------------------------------------------
// Compute the longwave reflectance and transmittance to diffuse
// radiation using the Meador & Weaver formulas, as well as the upward
// flux at the top and the downward flux at the base of the layer due
// to emission from within the layer assuming a linear variation of
// Planck function within the layer
template<bool IsActive>
void
calc_ref_trans_source_lw(int ng,
			 const Config& config,
			 const Array<1,IsActive>& od,
			 const Array<1,IsActive>& ssa,
			 const Array<1,IsActive>& asymmetry,
			 const Array<1,IsActive>& planck_top,
			 const Array<1,IsActive>& planck_bot,
			 Array<1,IsActive> reflectance,
			 Array<1,IsActive> transmittance,
			 Array<1,IsActive> source_up,
			 Array<1,IsActive> source_dn);

// ---------------------------------------------------------------------
// Calculate the transmittance and upward/downward sources (if
// allocated) in a particular direction (cosine of zenith angle, mu)
// for a subsequent radiance calculation
template<bool IsActive>
void
calc_radiance_trans_source_lw(int ng, int nlev,
			      const Config& config,
			      Real mu,
			      const Array<2,IsActive>& region_fracs,
			      const Array<2,IsActive>& planck_hl,
			      const Array<3,IsActive>& od,
			      const Array<3,IsActive>& ssa,
			      const Array<2,IsActive>& asymmetry,
			      const Array<3,IsActive>& flux_up_base,
			      const Array<3,IsActive>& flux_dn_top,
			      Array<3,IsActive> transmittance,
			      Array<3,IsActive> source_up,
			      Array<3,IsActive> source_dn);
  
}; // namespace tcrad


#endif
