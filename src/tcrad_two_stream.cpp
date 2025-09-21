#include "tcrad_two_stream.hpp"

namespace tcrad {

template<bool IsActive>
void
calc_no_scattering_transmittance_lw(int ng,
				    int nlev,
				    const aMatrix<IsActive>& od,
				    const aMatrix<IsActive>& planck_top,
				    const aMatrix<IsActive>& planck_bot,
  				    aMatrix<IsActive> transmittance,
				    aMatrix<IsActive> source_up,
				    aMatrix<IsActive> source_dn)
{
  typedef typename scalar<IsActive>::type aScalar;
  aScalar coeff, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot;

  transmittance = exp(-LW_DIFFUSIVITY*od);
  
  for (int jlev = 0; jlev < nlev; ++jlev) {
    for (int jg = 0; jg < ng; ++jg) {
      if (od(jlev,jg) > 1.0e-3) {
        coeff = (planck_bot(jlev,jg)-planck_top(jlev,jg)) / (LW_DIFFUSIVITY*od(jlev,jg));
        coeff_up_top  =  coeff + planck_top(jlev,jg);
        coeff_up_bot  =  coeff + planck_bot(jlev,jg);
        coeff_dn_top  = -coeff + planck_top(jlev,jg);
        coeff_dn_bot  = -coeff + planck_bot(jlev,jg);
        source_up(jlev,jg) =  coeff_up_top - transmittance(jlev,jg) * coeff_up_bot;
        source_dn(jlev,jg) =  coeff_dn_bot - transmittance(jlev,jg) * coeff_dn_top;
      }
      else {
        // Linear limit at low optical depth
        source_up(jlev,jg) = LW_DIFFUSIVITY * od(jlev,jg) * 0.5 * (planck_top(jlev,jg)+planck_bot(jlev,jg));
	source_dn(jlev,jg) = source_up(jlev,jg);
      }
    }
  }
}
  
template<bool IsActive>
void
calc_ref_trans_lw(int ng,
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
  aVector<IsActive> gamma1(ng), gamma2(ng), k_exponent(ng);
  aScalar reftrans_factor, factor;
  aScalar exponential2;
  aScalar coeff, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot;

  for (int jg = 0; jg < ng; ++jg) {
    factor = (LW_DIFFUSIVITY * 0.5) * ssa(jg);
    gamma1(jg) = LW_DIFFUSIVITY - factor*(1.0 + asymmetry(jg));
    gamma2(jg) = factor * (1.0 - asymmetry(jg));
    // Eq 18 of Meador & Weaver (1980)
    k_exponent(jg) = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), 1.0e-12));
  } 

  aVector<IsActive> exponential = exp(-k_exponent*od);
  
  for (int jg = 0; jg < ng; ++jg) {
    if (od(jg) > 1.0e-3) {
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


template
void
calc_no_scattering_transmittance_lw<false>(int ng,
					   int nlev,
					   const aMatrix<false>& od,
					   const aMatrix<false>& planck_top,
					   const aMatrix<false>& planck_bot,
					   aMatrix<false> transmittance,
					   aMatrix<false> source_up,
					   aMatrix<false> source_dn);
template
void
calc_no_scattering_transmittance_lw<true>(int ng,
					  int nlev,
					  const aMatrix<true>& od,
					  const aMatrix<true>& planck_top,
					  const aMatrix<true>& planck_bot,
					  aMatrix<true> transmittance,
					  aMatrix<true> source_up,
					  aMatrix<true> source_dn);

template
void
calc_ref_trans_lw(int ng,
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
calc_ref_trans_lw(int ng,
		  const aVector<true>& od,
		  const aVector<true>& ssa,
		  const aVector<true>& asymmetry,
		  const aVector<true>& planck_top,
		  const aVector<true>& planck_bot,
		  aVector<true> reflectance,
		  aVector<true> transmittance,
		  aVector<true> source_up,
		  aVector<true> source_dn);


};
