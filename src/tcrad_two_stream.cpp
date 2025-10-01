#include <cmath>
#include <limits>

#include "tcrad_two_stream.hpp"
#include "tcrad_config.hpp"

namespace tcrad {

static inline Real
get_lw_diffusivity(TwoStreamScheme scheme)
{
  if (  scheme == TWO_STREAM_EDDINGTON
      | scheme == TWO_STREAM_LEGENDRE) {
    return 2.0;
  }
  else {
    return 1.66;
  }
}

static inline Real
get_lw_diffusivity_cloud(TwoStreamScheme scheme)
{
  if (  scheme == TWO_STREAM_EDDINGTON
      | scheme == TWO_STREAM_LEGENDRE
      | scheme == TWO_STREAM_HYBRID
      | scheme == TWO_STREAM_SCALED_WISCOMBE_GRAMS) {
    return 2.0;
  }
  else {
    return 1.66;
  }
}

// Precision-dependent parameters such that
// param<Real>::OPTICAL_DEPTH_THRESHOLD always returns the correct
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
  

// For clear skies
template<bool IsActive>
void
calc_no_scattering_transmittance_lw(int ng,
				    int nlev,
				    const Config& config,
				    const aMatrix<IsActive>& od,
				    const aMatrix<IsActive>& planck_top,
				    const aMatrix<IsActive>& planck_bot,
  				    aMatrix<IsActive> transmittance,
				    aMatrix<IsActive> source_up,
				    aMatrix<IsActive> source_dn)
{
  typedef typename scalar<IsActive>::type aScalar;
  aScalar coeff, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot;

  const Real lw_diffusivity = get_lw_diffusivity(config.i_two_stream_scheme);
  
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

  // Cloud only
template<bool IsActive>
void
calc_ref_trans_lw(int ng,
		  const Config& config,
		  const aVector<IsActive>& od,
		  const aVector<IsActive>& ssa,
		  const aVector<IsActive>& asymmetry,
		  const aVector<IsActive>& planck_top,
		  const aVector<IsActive>& planck_bot,
		  aVector<IsActive> reflectance,
		  aVector<IsActive> transmittance,
		  aVector<IsActive> source_up,
		  aVector<IsActive> source_dn)
{
  typedef typename scalar<IsActive>::type aScalar;
  aVector<IsActive> gamma1(ng), gamma2(ng);
  aScalar reftrans_factor, factor;
  aScalar exponential2;
  aScalar coeff, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot;

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
  aVector<IsActive> k_exponent = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2),
					  param<Real>::min_k_squared()));
  
  aVector<IsActive> exponential = exp(-k_exponent*od);
  
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

template<bool IsActive>
void
calc_transmittance(int ng, int nlev,
		   Real mu,                      // Cosine of sensor zenith angle
		   const aArray<1,IsActive>& clear_fraction, // (nlev,NREGIONS)
		   const aArray<3,IsActive>& od, // (nlev,NREGIONS,ng)
		   aArray<3,IsActive> transmittance) { // (nlev,NREGIONS,ng)
  Real minus_secant = -1.0 / std::abs(mu);
  for (int jlev = 0; jlev < nlev; ++jlev) {
    if (clear_fraction(jlev) == 1.0) {
      transmittance(jlev,0,__) = exp(minus_secant * od(jlev,0,__));
    }
    else {
      transmittance(jlev,__,__) = exp(minus_secant * od(jlev,__,__));
    }
  }
}

template<bool IsActive>
void
calc_radiance_source(int ng, int nlev,
		     const Config& config,
		     Real mu,
		     const aArray<2,IsActive>& reg_fracs,
		     const aArray<2,IsActive>& planck_hl,
		     const aArray<3,IsActive>& od,
		     const aArray<3,IsActive>& ssa,
		     const aArray<2,IsActive>& asymmetry,
		     const aArray<3,IsActive>& flux_up_base,
		     const aArray<3,IsActive>& flux_dn_top,
		     const aArray<3,IsActive>& transmittance,
		     aArray<3,IsActive> source_up,
		     aArray<3,IsActive> source_dn) {

  aArray<2,IsActive> planck_top (NREGIONS,ng);
  aArray<2,IsActive> planck_base(NREGIONS,ng);
  aArray<2,IsActive> source_top (NREGIONS,ng);
  aArray<2,IsActive> source_base(NREGIONS,ng);

  aArray<1,IsActive> gamma1(ng), gamma2(ng), exponential(ng), k_exponent(ng);

  typedef typename scalar<IsActive>::type aScalar;
  aScalar p_same, p_opposite, y_both, planck_prime, x_up, x_dn, coeff;
  aScalar rt_factor, c1, c2, factor, scaling1, scaling2, one_minus_kmu;

  const Real lw_diffusivity = get_lw_diffusivity_cloud(config.i_two_stream_scheme);
  
  for (int jlev = 0; jlev < nlev; ++jlev) {
    int imaxreg = 0;
    if (reg_fracs(jlev,0) < 1.0) {
      imaxreg = NREGIONS-1;
    }
    if (imaxreg > 1) {
      // Cloudy layer: scale the Planck terms by the region fraction
      // and also by the single-scattering co-albedo
      planck_top(0,__)  = planck_hl(jlev,__)   * reg_fracs(jlev,0);
      planck_base(0,__) = planck_hl(jlev+1,__) * reg_fracs(jlev,0);
      for (int jreg = 1; jreg < NREGIONS; ++jreg) {
	planck_top(jreg,__)  = planck_hl(jlev,__)   * reg_fracs(jlev,jreg);
	planck_base(jreg,__) = planck_hl(jlev+1,__) * reg_fracs(jlev,jreg);
      }
    }
    else {
      // Clear layer
      planck_top(0,__)  = planck_hl(jlev,__);
      planck_base(0,__) = planck_hl(jlev+1,__);
    }

    // Clear region
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
      aVector<IsActive> k_exponent = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2),
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
  
  /*
template<bool IsActive>
void
calc_radiance_trans_source(int ng,
			   int nlev,
			   Real mu,
			   const Config& config,
			   const aArray<2,IsActive>& region_fracs,
			   const aArray<2,IsActive>& planck_hl,    // (nlev+1,ng)
  */

  

template
void
calc_no_scattering_transmittance_lw<false>(int ng,
					   int nlev,
					   const Config& config,
					   const aMatrix<false>& od,
					   const aMatrix<false>& planck_top,
					   const aMatrix<false>& planck_bot,
					   aMatrix<false> transmittance,
					   aMatrix<false> source_up,
					   aMatrix<false> source_dn);

template
void
calc_ref_trans_lw<false>(int ng,
			 const Config& config,
			 const aVector<false>& od,
			 const aVector<false>& ssa,
			 const aVector<false>& asymmetry,
			 const aVector<false>& planck_top,
			 const aVector<false>& planck_bot,
			 aVector<false> reflectance,
			 aVector<false> transmittance,
			 aVector<false> source_up,
			 aVector<false> source_dn);

template
void
calc_transmittance<false>(int ng, int nlev,
			  Real mu,                      // Cosine of sensor zenith angle
			  const aArray<1,false>& clear_fraction, // (nlev,NREGIONS)
			  const aArray<3,false>& od, // (nlev,NREGIONS,ng)
			  aArray<3,false> transmittance); // (nlev,NREGIONS,ng)
  
template
void
calc_radiance_source<false>(int ng, int nlev,
			    const Config& config,
			    Real mu,
			    const aArray<2,false>& reg_fracs,
			    const aArray<2,false>& planck_hl,
			    const aArray<3,false>& od,
			    const aArray<3,false>& ssa,
			    const aArray<2,false>& asymmetry,
			    const aArray<3,false>& flux_up_base,
			    const aArray<3,false>& flux_dn_top,
			    const aArray<3,false>& transmittance,
			    aArray<3,false> source_up,
			    aArray<3,false> source_dn);

#if ADEPT_REAL_TYPE_SIZE == 8
template
void
calc_no_scattering_transmittance_lw<true>(int ng,
					  int nlev,
					  const Config& config,
					  const aMatrix<true>& od,
					  const aMatrix<true>& planck_top,
					  const aMatrix<true>& planck_bot,
					  aMatrix<true> transmittance,
					  aMatrix<true> source_up,
					  aMatrix<true> source_dn);


template
void
calc_ref_trans_lw<true>(int ng,
			const Config& config,
			const aVector<true>& od,
			const aVector<true>& ssa,
			const aVector<true>& asymmetry,
			const aVector<true>& planck_top,
			const aVector<true>& planck_bot,
			aVector<true> reflectance,
			aVector<true> transmittance,
			aVector<true> source_up,
			aVector<true> source_dn);

template
void
calc_transmittance<true>(int ng, int nlev,
			 Real mu,                      // Cosine of sensor zenith angle
			 const aArray<1,true>& clear_fraction, // (nlev,NREGIONS)
			 const aArray<3,true>& od, // (nlev,NREGIONS,ng)
			 aArray<3,true> transmittance); // (nlev,NREGIONS,ng)

template
void
calc_radiance_source<true>(int ng, int nlev,
			   const Config& config,
			   Real mu,
			   const aArray<2,true>& reg_fracs,
			   const aArray<2,true>& planck_hl,
			   const aArray<3,true>& od,
			   const aArray<3,true>& ssa,
			   const aArray<2,true>& asymmetry,
			   const aArray<3,true>& flux_up_base,
			   const aArray<3,true>& flux_dn_top,
			   const aArray<3,true>& transmittance,
			   aArray<3,true> source_up,
			   aArray<3,true> source_dn);


#endif



};
