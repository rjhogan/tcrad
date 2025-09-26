// radiation_radiance_propagation_lw.hpp - Radiance propagation in longwave
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

#include "tcrad_adept.hpp"

namespace tcrad {

template <bool IsActive>
void
calc_radiance_dn(int ng, int nlev,
		 const aArray<1,IsActive>& clear_fraction,
		 const aArray<3,IsActive>& transmittance,
		 const aArray<3,IsActive>& source,
		 const aArray<3,IsActive>& v_matrix,
		 aArray<2,IsActive> radiance_base);

template <bool IsActive>
void
calc_radiance_up(int ng, int nlev,
		 const aArray<1,IsActive>& clear_fraction,
		 const aArray<2,IsActive>& radiance_up_surf,
		 const aArray<3,IsActive>& transmittance,
		 const aArray<3,IsActive>& source,
		 const aArray<3,IsActive>& u_matrix,
		 aArray<1,IsActive> radiance_toa);
}
