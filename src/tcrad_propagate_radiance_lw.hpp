// radiation_radiance_propagation_lw.hpp - Radiance propagation in longwave
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

#include "tcrad_adept.hpp"

namespace tcrad {

// ---------------------------------------------------------------------
// Compute downward radiance profile by solving the Schwarzschild
// radiative transfer equation in each region of each layer assuming
// the source term to vary linearly with optical depth in each
// layer. The overlap matrix is used to translate the radiances exiting
// the base of one layer into radiances in the top of the layer below,
// consistent with the Tripleclouds approximation. 
template <bool IsActive>
void
propagate_radiance_dn_lw(int ng, int nlev,
			 const Array<1,IsActive>& clear_fraction,
			 const Array<3,IsActive>& transmittance,
			 const Array<3,IsActive>& source,
			 const Array<3,IsActive>& v_overlap,
			 Array<2,IsActive> radiance_base);

// ---------------------------------------------------------------------
// Compute upward radiance profile by solving the Schwarzschild
// radiative transfer equation in each region of each layer assuming
// the source term to vary linearly with optical depth in each
// layer. The overlap matrix is used to translate the radiances exiting
// the top of one layer into radiances in the base of the layer above,
// consistent with the Tripleclouds approximation. 
template <bool IsActive>
void
propagate_radiance_up_lw(int ng, int nlev,
			 const Array<1,IsActive>& clear_fraction,
			 const Array<2,IsActive>& radiance_up_surf,
			 const Array<3,IsActive>& transmittance,
			 const Array<3,IsActive>& source,
			 const Array<3,IsActive>& u_overlap,
			 Array<1,IsActive> radiance_toa);
}
