// tcrad_flux_lw.hpp - Compute longwave fluxes
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

#ifndef TCRAD_FLUX_LW
#define TCRAD_FLUX_LW 1

#include "tcrad_adept.hpp"
#include "tcrad_config.hpp"

namespace tcrad {
  
// ---------------------------------------------------------------------
// Longwave "Tripleclouds" solver following Shonk & Hogan (2008) that
// treats cloud inhomogeneity by dividing each model level into three
// regions, one clear and two cloudy (with differing optical depth),
// returning a profile of upwelling and downwelling spectral fluxes.
template <bool IsActive>
void calc_tripleclouds_flux_lw(int ng,
			       int nlev,
			       int ncol,
			       const Config& config,
			       const Array<2,IsActive>& surf_emission,
			       const Array<2,IsActive>& surf_albedo,
			       const Array<3,IsActive>& planck_hl,
			       const Array<2,IsActive>& cloud_fraction,
			       const Array<2>& fractional_std,
			       const Array<3,IsActive>& od_clear,
			       const Array<3,IsActive>& od_cloud,
			       const Array<3,IsActive>& ssa_cloud,
			       const Array<3,IsActive>& asymmetry_cloud,
			       const Array<2,IsActive>& overlap_param,
			       Array<3,IsActive> flux_up,
			       Array<3,IsActive> flux_dn);
}
#endif
