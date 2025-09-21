#ifndef ECRAD_REGIONS_H
#define ECRAD_REGIONS_H 1

#include "tcrad_adept.hpp"

namespace tcrad {

template <bool IsActive> 
void
calc_region_properties(int nlev,
		       int nreg,
		       int ncol,
		       const aMatrix<IsActive>& cloud_fraction,
		       const Matrix& fractional_std,
		       aArray<3,IsActive> reg_fracs,
		       aArray<3,IsActive> od_scaling,
		       Real cloud_fraction_threshold);

};

#endif
