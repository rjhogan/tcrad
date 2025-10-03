// tcrad_layer_solutions_lw.cpp - Solutions to longwave two-stream equations in a layer
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
// This file is adapted from ecRad's radiation_two_stream.F90

#include <cmath>
#include <limits>

#include "tcrad_layer_solutions_lw.hpp"
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
calc_no_scat_trans_source_lw(int ng,
			     int nlev,
			     const Config& config,
			     const Array<2,IsActive>& od,
			     const Array<2,IsActive>& planck_top,
			     const Array<2,IsActive>& planck_bot,
			     Array<2,IsActive> transmittance,
			     Array<2,IsActive> source_up,
			     Array<2,IsActive> source_dn)
{
  typedef typename scalar<IsActive>::type Scalar;
  Scalar coeff, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot;

  const Real lw_diffusivity = get_lw_diffusivity(config.i_two_stream_scheme);

  // Outside loop so can vectorize along spectral dimension
  transmittance = exp(-lw_diffusivity*od);
  
  for (int jlev = 0; jlev < nlev; ++jlev) {
    for (int jg = 0; jg < ng; ++jg) {
      if (od(jlev,jg) > 1.0e-3) {
        coeff = (planck_bot(jlev,jg)-planck_top(jlev,jg)) / (lw_diffusivity*od(jlev,jg));
        coeff_up_top  =  coeff + planck_top(jlev,jg);
        coeff_up_bot  =  coeff + planck_bot(jlev,jg);
        coeff_dn_top  = -coeff + planck_top(jlev,jg);
        coeff_dn_bot  = -coeff + planck_bot(jlev,jg);
        source_up(jlev,jg) =  coeff_up_top - transmittance(jlev,jg) * coeff_up_bot;
        source_dn(jlev,jg) =  coeff_dn_bot - transmittance(jlev,jg) * coeff_dn_top;
      }
      else {
        // Linear limit at low optical depth
        source_up(jlev,jg) = lw_diffusivity * od(jlev,jg) * 0.5 * (planck_top(jlev,jg)+planck_bot(jlev,jg));
	source_dn(jlev,jg) = source_up(jlev,jg);
      }
    }
  }
}

// ---------------------------------------------------------------------
// Compute the longwave reflectance and transmittance to diffuse
// radiation using the Meador & Weaver formulas, as well as the upward
// flux at the top and the downward flux at the base of the layer due
// to emission from within the layer assuming a linear variation of
// Planck function within the layer. 
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

// ---------------------------------------------------------------------
// Compute the transmittance and source terms for a radiance
// computation. If source_up or source_dn are empty arrays, they are
// treated as optional arguments and not filled.
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
			      Array<3,IsActive> source_dn) {

  Array<2,IsActive> planck_top (NREGIONS,ng);
  Array<2,IsActive> planck_base(NREGIONS,ng);
  Array<2,IsActive> source_top (NREGIONS,ng);
  Array<2,IsActive> source_base(NREGIONS,ng);

  Array<1,IsActive> gamma1(ng), gamma2(ng), exponential(ng), k_exponent(ng);

  typedef typename scalar<IsActive>::type Scalar;
  Scalar p_same, p_opposite, planck_prime, coeff;
  Scalar rt_factor, c1, c2, factor, scaling1, scaling2, one_minus_kmu;

  const Real lw_diffusivity = get_lw_diffusivity_cloud(config.i_two_stream_scheme);
  const Real minus_secant = -1.0 / std::abs(mu);

  for (int jlev = 0; jlev < nlev; ++jlev) {
    int imaxreg = 0;
    if (region_fracs(jlev,0) < 1.0) {
      // Layer contains cloud
      imaxreg = NREGIONS-1;
      // Transmittance of all regions
      transmittance(jlev,__,__) = exp(minus_secant * od(jlev,__,__));
      // Scale the Planck terms by the region fraction and also by the
      // single-scattering co-albedo
      planck_top(0,__)  = planck_hl(jlev,__)   * region_fracs(jlev,0);
      planck_base(0,__) = planck_hl(jlev+1,__) * region_fracs(jlev,0);
      for (int jreg = 1; jreg < NREGIONS; ++jreg) {
	planck_top(jreg,__)  = planck_hl(jlev,__)   * region_fracs(jlev,jreg);
	planck_base(jreg,__) = planck_hl(jlev+1,__) * region_fracs(jlev,jreg);
      }
    }
    else {
      // Clear layer: transmittance of only the clear region
      transmittance(jlev,0,__) = exp(minus_secant * od(jlev,0,__));
      planck_top(0,__)  = planck_hl(jlev,__);
      planck_base(0,__) = planck_hl(jlev+1,__);
    }

    // Compute the sources of the clear region
    if (!source_up.empty()) {
      for (int jg = 0; jg < ng; ++jg) {
	if (od(jlev,0,jg) > param<Real>::optical_depth_threshold()) {
	  planck_prime = (planck_base(0,jg)-planck_top(0,jg))
	    / od(jlev,0,jg);
	  source_up(jlev,0,jg) = planck_top(0,jg)
	    - planck_base(0,jg) * transmittance(jlev,0,jg)
	    + planck_prime*mu*(1.0 - transmittance(jlev,0,jg));
	}
	else {
	  // At low optical depths the effective Planck function is
          // half the top and bottom values, and we avoid the division
          // by optical depth
	  source_up(jlev,0,jg) = od(jlev,0,jg) * 0.5
	    * (planck_base(0,jg)+planck_top(0,jg)) / mu;
	}
      }
      source_up(jlev,0,__) *= (1.0/tcrad::PI);
    }
    if (!source_dn.empty()) {
      for (int jg = 0; jg < ng; ++jg) {
	if (od(jlev,0,jg) > param<Real>::optical_depth_threshold()) {
	  planck_prime = (planck_base(0,jg)-planck_top(0,jg))
	    / od(jlev,0,jg);
	  source_dn(jlev,0,jg) = planck_base(0,jg)
	    - planck_top(0,jg) * transmittance(jlev,0,jg)
	    - planck_prime*mu*(1.0 - transmittance(jlev,0,jg));
	}
	else {
	  // At low optical depths the effective Planck function is
          // half the top and bottom values, and we avoid the division
          // by optical depth
	  source_dn(jlev,0,jg) = od(jlev,0,jg) * 0.5
	    * (planck_base(0,jg)+planck_top(0,jg)) / mu;
	}
      }
      source_dn(jlev,0,__) *= (1.0/tcrad::PI);
    }
    
    // Scattering from two-stream fluxes: loop over cloudy regions
    for (int jreg = 1; jreg <= imaxreg; ++jreg) {
      if (config.i_two_stream_scheme == TWO_STREAM_EDDINGTON) {
	for (int jg = 0; jg < ng; ++jg) {
	  // See Meador & Weaver (1980), Table 1; Toon et al. (1989), Table 1
	  gamma1(jg) = 1.75 - ssa(jlev,jreg,jg) * (1.0 + 0.75*asymmetry(jlev,jg));
	  gamma2(jg) = ssa(jlev,jreg,jg) * (1.0 - 0.75*asymmetry(jlev,jg) - 0.25);
	}
      }
      else if (config.i_two_stream_scheme == TWO_STREAM_SCALED_WISCOMBE_GRAMS) {
	for (int jg = 0; jg < ng; ++jg) {
	  // Wiscombe-Grams backscatter fraction applied to de-scaled
	  // asymmety factor
	  factor = 0.5 * (1.0 - 0.75*asymmetry(jlev,jg)/(1.0-asymmetry(jlev,jg)));
	  gamma1(jg) = lw_diffusivity * (1.0 - ssa(jlev,jreg,jg)*(1.0-factor));
	  gamma2(jg) = lw_diffusivity * ssa(jlev,jreg,jg) * factor;
	}
      }
      else { // TWO_STREAM_LEGENDRE or TWO_STREAM_ELSASSER
	for (int jg = 0; jg < ng; ++jg) {
	  factor = (lw_diffusivity * 0.5) * ssa(jlev,jreg,jg);
	  gamma1(jg) = lw_diffusivity - factor*(1.0 + asymmetry(jlev,jg));
	  gamma2(jg) = factor * (1.0 - asymmetry(jlev,jg));
	}
      } 
      // Eq 18 of Meador & Weaver (1980)
      k_exponent = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2),
			    param<Real>::min_k_squared()));
      
      exponential = exp(-k_exponent*od(jlev,jreg,__));

      for (int jg = 0; jg < ng; ++jg) {
	if (od(jlev,jreg,jg) > param<Real>::optical_depth_threshold()) {
	  p_same = 3.0 * asymmetry(jlev,jg) * mu / lw_diffusivity;
	  // Phase function from downwelling flux to upwelling
	  // radiance (or up to down)
	  p_opposite = 1.0 - p_same;
	  // Phase functions from upwelling flux to upwelling radiance
	  // (or down to down)
	  p_same += 1.0;

	  planck_prime = (planck_base(jreg,jg)-planck_top(jreg,jg))
	    / od(jlev,jreg,jg);
	  coeff = planck_prime / (gamma1(jg)+gamma2(jg));
	  rt_factor = 1.0 / (k_exponent(jg) + gamma1(jg) + (k_exponent(jg)-gamma1(jg))
			     *exponential(jg)*exponential(jg));
	  factor = exponential(jg) * gamma2(jg) / (gamma1(jg) + k_exponent(jg));
	  c1 = rt_factor * (flux_up_base(jlev,jreg,jg) - factor*flux_dn_top(jlev,jreg,jg)
			    - (planck_base(jreg,jg)+coeff) + factor*(planck_top(jreg,jg)-coeff));
	  c2 = rt_factor * (flux_dn_top(jlev,jreg,jg) - factor*flux_up_base(jlev,jreg,jg)
			    -(planck_top(jreg,jg)-coeff) + factor*(planck_base(jreg,jg)+coeff));

	  // Convert to scaling factors...
	  one_minus_kmu = 1.0 - k_exponent(jg)*mu;
			
	  scaling1 = (exponential(jg) - transmittance(jlev,jreg,jg))
	    / copysign(max(abs(one_minus_kmu),std::numeric_limits<Real>::epsilon()),
		       one_minus_kmu);
	  scaling2 = (1.0-exponential(jg)*transmittance(jlev,jreg,jg))/(1.0+k_exponent(jg)*mu);

	  if (!source_up.empty()) {
	    source_up(jlev,jreg,jg) =
	      // Direct emission plus scattering from the part of the
	      // fluxes due to internal emission and having a linear
	      // structure
	      (1.0 - transmittance(jlev,jreg,jg))
	      * (0.5*ssa(jlev,jreg,jg)*planck_prime*(p_same-p_opposite)
		 /(gamma1(jg)+gamma2(jg)) + planck_prime*mu)
	      + planck_top(jreg,jg)
	      -planck_base(jreg,jg)*transmittance(jlev,jreg,jg)
	      // Scattering from the exponential part of the flux,
	      // whether caused by external or internal sources
	      + 0.5*ssa(jlev,jreg,jg)
	      * (p_same * ((gamma1(jg)+k_exponent(jg))*scaling1*c1
			   + gamma2(jg)*scaling2*c2)
		 + p_opposite * (gamma2(jg)*scaling1*c1
				 + (gamma1(jg)+k_exponent(jg))*scaling2*c2));
	    source_up(jlev,jreg,jg) /= tcrad::PI;
	  }
	  if (!source_dn.empty()) {
	    source_dn(jlev,jreg,jg) = 
	      // Direct emission plus scattering from the part of the
	      // fluxes due to internal emission and having a linear
	      // structure
	      -(1.0 - transmittance(jlev,jreg,jg))
	      * (0.5*ssa(jlev,jreg,jg)*planck_prime*(p_same-p_opposite)
		 /(gamma1(jg)+gamma2(jg)) + planck_prime*mu)
	      + planck_base(jreg,jg)
	      -planck_top(jreg,jg)*transmittance(jlev,jreg,jg)
	      // Scattering from the exponential part of the flux,
	      // whether caused by external or internal sources
	      + 0.5*ssa(jlev,jreg,jg)
	      * (p_opposite * ((gamma1(jg)+k_exponent(jg))*scaling2*c1
			   + gamma2(jg)*scaling1*c2)
		 + p_same * (gamma2(jg)*scaling2*c1
			     + (gamma1(jg)+k_exponent(jg))*scaling1*c2));
	    source_dn(jlev,jreg,jg) /= tcrad::PI;
	  }
	}
	else {
	  if (!source_up.empty()) {
	    source_up(jlev,jreg,jg) = od(jlev,jreg,jg)
	      * 0.5 * (planck_base(jreg,jg)+planck_top(jreg,jg)) / mu;
	    source_up(jlev,jreg,jg) /= tcrad::PI;
	  }
	  if (!source_dn.empty()) {
	    source_dn(jlev,jreg,jg) = od(jlev,jreg,jg)
	      * 0.5 * (planck_base(jreg,jg)+planck_top(jreg,jg)) / mu;
	    source_dn(jlev,jreg,jg) /= tcrad::PI;
	  }
	}
      } // Loop over jg
    } // Loop over jreg
  } // Loop over jlev
}


template
void
calc_no_scat_trans_source_lw<false>(int ng,
				    int nlev,
				    const Config& config,
				    const Array<2,false>& od,
				    const Array<2,false>& planck_top,
				    const Array<2,false>& planck_bot,
				    Array<2,false> transmittance,
				    Array<2,false> source_up,
				    Array<2,false> source_dn);
// Instantiate direct functions
template
void
calc_ref_trans_source_lw<false>(int ng,
				const Config& config,
				const Array<1,false>& od,
				const Array<1,false>& ssa,
				const Array<1,false>& asymmetry,
				const Array<1,false>& planck_top,
				const Array<1,false>& planck_bot,
				Array<1,false> reflectance,
				Array<1,false> transmittance,
				Array<1,false> source_up,
				Array<1,false> source_dn);

template
void
calc_radiance_trans_source_lw<false>(int ng, int nlev,
				     const Config& config,
				     Real mu,
				     const Array<2,false>& region_fracs,
				     const Array<2,false>& planck_hl,
				     const Array<3,false>& od,
				     const Array<3,false>& ssa,
				     const Array<2,false>& asymmetry,
				     const Array<3,false>& flux_up_base,
				     const Array<3,false>& flux_dn_top,
				     Array<3,false> transmittance,
				     Array<3,false> source_up,
				     Array<3,false> source_dn);

// Instantiate active functions
#if ADEPT_REAL_TYPE_SIZE == 8
template
void
calc_no_scat_trans_source_lw<true>(int ng,
				   int nlev,
				   const Config& config,
				   const Array<2,true>& od,
				   const Array<2,true>& planck_top,
				   const Array<2,true>& planck_bot,
				   Array<2,true> transmittance,
				   Array<2,true> source_up,
				   Array<2,true> source_dn);


template
void
calc_ref_trans_source_lw<true>(int ng,
			       const Config& config,
			       const Array<1,true>& od,
			       const Array<1,true>& ssa,
			       const Array<1,true>& asymmetry,
			       const Array<1,true>& planck_top,
			       const Array<1,true>& planck_bot,
			       Array<1,true> reflectance,
			       Array<1,true> transmittance,
			       Array<1,true> source_up,
			       Array<1,true> source_dn);

template
void
calc_radiance_trans_source_lw<true>(int ng, int nlev,
				    const Config& config,
				    Real mu,
				    const Array<2,true>& region_fracs,
				    const Array<2,true>& planck_hl,
				    const Array<3,true>& od,
				    const Array<3,true>& ssa,
				    const Array<2,true>& asymmetry,
				    const Array<3,true>& flux_up_base,
				    const Array<3,true>& flux_dn_top,
				    Array<3,true> transmittance,
				    Array<3,true> source_up,
				    Array<3,true> source_dn);


#endif

};
