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
// Return diffusivity factor for a particular two-stream scheme,
// allowing for separate values between clear and cloudy layers
static inline Real
get_lw_diffusivity(TwoStreamScheme scheme)
{
  if (  (scheme == TWO_STREAM_EDDINGTON)
      | (scheme == TWO_STREAM_LEGENDRE) ) {
    return 2.0;
  }
  else {
    return 1.66;
  }
}

static inline Real
get_lw_diffusivity_cloud(TwoStreamScheme scheme)
{
  if (  (scheme == TWO_STREAM_EDDINGTON)
      | (scheme == TWO_STREAM_LEGENDRE)
      | (scheme == TWO_STREAM_HYBRID)
      | (scheme == TWO_STREAM_SCALED_WISCOMBE_GRAMS) ) {
    return 2.0;
  }
  else {
    return 1.66;
  }
}

// ---------------------------------------------------------------------
// Precision-dependent parameters such that
// param<Real>::optical_depth_threshold() always returns the correct
// value
template <typename Type>
struct param {
  static Type optical_depth_threshold() { return 1.0e-7;  }
  static Type min_k_squared()           { return 1.0e-12; }
};

template<>
struct param<float> {
  static float optical_depth_threshold() { return 1.0e-4;  }
  static float min_k_squared()           { return 1.0e-6; }
};
  

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
calc_ref_trans_source_lw(int ng, int nlev,
			 const Config& config,
			 const Array<2,IsActive>& region_fracs,
			 const Array<3,IsActive>& od,
			 const Array<3,IsActive>& ssa,
			 const Array<2,IsActive>& asymmetry,
			 const Array<2,IsActive>& planck_top,
			 const Array<2,IsActive>& planck_bot,
			 Array<3,IsActive> reflectance,
			 Array<3,IsActive> transmittance,
			 Array<3,IsActive> source_up,
			 Array<3,IsActive> source_dn);

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


// ---------------------------------------------------------------------
// Compute the longwave reflectance and transmittance to diffuse
// radiation using the Meador & Weaver formulas, as well as the upward
// flux at the top and the downward flux at the base of the layer due
// to emission from within the layer assuming a linear variation of
// Planck function within the layer. 
template<bool IsActive>
void
calc_ref_trans_source_layer_lw(int ng,
			 const Config& config,
			 const Array<1,IsActive>& od,
			 const Array<1,IsActive>& ssa,
			 const Array<1,IsActive>& asymmetry,
			 const Array<1,IsActive>& planck_top,
			 const Array<1,IsActive>& planck_bot,
			 Array<1,IsActive> reflectance,
			 Array<1,IsActive> transmittance,
			 Array<1,IsActive> source_up,
			 Array<1,IsActive> source_dn)
{
  typedef typename scalar<IsActive>::type Scalar;
  
  Array<1,IsActive> gamma1(ng), gamma2(ng);
  Scalar reftrans_factor, factor;
  Scalar exponential2;
  Scalar coeff, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot;

  const Real lw_diffusivity = get_lw_diffusivity_cloud(config.i_two_stream_scheme);

  if (config.i_two_stream_scheme == TWO_STREAM_EDDINGTON) {
    for (int jg = 0; jg < ng; ++jg) {
      // See Meador & Weaver (1980), Table 1; Toon et al. (1989), Table 1
      gamma1(jg) = 1.75 - ssa(jg) * (1.0 + 0.75*asymmetry(jg));
      gamma2(jg) = ssa(jg) * (1.0 - 0.75*asymmetry(jg) - 0.25);
    }
  }
  else if (config.i_two_stream_scheme == TWO_STREAM_SCALED_WISCOMBE_GRAMS) {
    for (int jg = 0; jg < ng; ++jg) {
      // Wiscombe-Grams backscatter fraction applied to de-scaled
      // asymmety factor
      factor = 0.5 * (1.0 - 0.75*asymmetry(jg)/(1.0-asymmetry(jg)));
      gamma1(jg) = lw_diffusivity * (1.0 - ssa(jg)*(1.0-factor));
      gamma2(jg) = lw_diffusivity * ssa(jg) * factor;
    }
  }
  else { // TWO_STREAM_LEGENDRE or TWO_STREAM_ELSASSER
    for (int jg = 0; jg < ng; ++jg) {
      factor = (lw_diffusivity * 0.5) * ssa(jg);
      gamma1(jg) = lw_diffusivity - factor*(1.0 + asymmetry(jg));
      gamma2(jg) = factor * (1.0 - asymmetry(jg));
    }
  } 
  // Eq 18 of Meador & Weaver (1980)
  Array<1,IsActive> k_exponent = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2),
					 param<Real>::min_k_squared()));
  
  Array<1,IsActive> exponential = exp(-k_exponent*od);
  
  for (int jg = 0; jg < ng; ++jg) {
    if (od(jg) > param<Real>::optical_depth_threshold()) {
      exponential2 = exponential(jg)*exponential(jg);
      reftrans_factor = 1.0 / (k_exponent(jg) + gamma1(jg) + (k_exponent(jg) - gamma1(jg))*exponential2);
      // Meador & Weaver (1980) Eq. 25
      reflectance(jg) = gamma2(jg) * (1.0 - exponential2) * reftrans_factor;
      // Meador & Weaver (1980) Eq. 26
      transmittance(jg) = 2.0 * k_exponent(jg) * exponential(jg) * reftrans_factor;

      // Compute upward and downward emission assuming the Planck
      // function to vary linearly with optical depth within the layer
      // (e.g. Wiscombe , JQSRT 1976).

      // Stackhouse and Stephens (JAS 1991) Eqs 5 & 12
      coeff = (planck_bot(jg)-planck_top(jg)) / (od(jg)*(gamma1(jg)+gamma2(jg)));
      coeff_up_top  =  coeff + planck_top(jg);
      coeff_up_bot  =  coeff + planck_bot(jg);
      coeff_dn_top  = -coeff + planck_top(jg);
      coeff_dn_bot  = -coeff + planck_bot(jg);
      source_up(jg) =  coeff_up_top - reflectance(jg) * coeff_dn_top
	- transmittance(jg) * coeff_up_bot;
      source_dn(jg) =  coeff_dn_bot - reflectance(jg) * coeff_up_bot
	- transmittance(jg) * coeff_dn_top;
    }
    else {
      reflectance(jg) = gamma2(jg) * od(jg);
      transmittance(jg) = (1.0 - k_exponent(jg)*od(jg)) / (1.0 + od(jg)*(gamma1(jg)-k_exponent(jg)));
      source_up(jg) = (1.0 - reflectance(jg) - transmittance(jg))
	* 0.5 * (planck_top(jg) + planck_bot(jg));
      source_dn(jg) = source_up(jg);
    }
  }  
}
  
}; // namespace tcrad


#endif
