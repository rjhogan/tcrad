#ifndef ECRAD_TWO_STREAM_H
#define ECRAD_TWO_STREAM_H 1

#include "tcrad_adept.hpp"

namespace tcrad {

static const Real LW_DIFFUSIVITY = 1.66;
  
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

}; // namespace tcrad


#endif
