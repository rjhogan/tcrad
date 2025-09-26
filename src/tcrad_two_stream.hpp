#ifndef ECRAD_TWO_STREAM_H
#define ECRAD_TWO_STREAM_H 1

#include "tcrad_adept.hpp"
#include "tcrad_config.hpp"

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
				    aMatrix<IsActive> source_dn);

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
		  aVector<IsActive> source_dn);

template<bool IsActive>
void
calc_transmittance(int ng, int nlev,
		   Real mu,                      // Cosine of sensor zenith angle
		   const aArray<1,IsActive>& clear_fraction, // (nlev,NREGIONS)
		   const aArray<3,IsActive>& od, // (nlev,NREGIONS,ng)
		   aArray<3,IsActive> transmittance); // (nlev,NREGIONS,ng)

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
		     aArray<3,IsActive> source_dn);
  
}; // namespace tcrad


#endif
