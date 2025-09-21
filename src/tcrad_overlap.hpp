#ifndef ECRAD_OVERLAP_H
#define ECRAD_OVERLAP_H 1

#include "tcrad_adept.hpp"

namespace tcrad {

template <bool IsActive>
void
calc_overlap_matrices(int nlev,
		      int ncol,
		      const aArray<3,IsActive>& region_fracs,
		      const aMatrix<IsActive>& overlap_param,
		      aArray<4,IsActive> u_matrix,
		      aArray<4,IsActive> v_matrix,
		      Real cloud_fraction_threshold);

};

#endif
