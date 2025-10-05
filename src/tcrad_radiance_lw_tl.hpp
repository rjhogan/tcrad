// tcrad_radiance_lw_tl.hpp - Tangent-linear of longwave radiance model
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
// Tangent-linear of longwave radiance model
void calc_tripleclouds_radiance_lw_tl(int ng,
				      int nlev,
				      const Config& config,
				      Real mu, // Cosine of sensor zenith angle
				      const Array<1>& surf_emission,
				      const Array<1>& surf_albedo,
				      const Array<2>& planck_hl,
				      const Array<1>& cloud_fraction,
				      const Array<1>& fractional_std,
				      const Array<2>& od_clear,
				      const Array<2>& od_cloud,
				      const Array<2>& ssa_cloud,
				      const Array<2>& asymmetry_cloud,
				      const Array<1>& overlap_param,
				      Array<1> radiance,
				      const Array<1>& surf_emission_tl,
				      const Array<1>& surf_albedo_tl,
				      const Array<2>& planck_hl_tl,
				      const Array<1>& cloud_fraction_tl,
				      const Array<2>& od_clear_tl,
				      const Array<2>& od_cloud_tl,
				      const Array<2>& ssa_cloud_tl,
				      const Array<2>& asymmetry_cloud_tl,
				      const Array<1>& overlap_param_tl,
				      Array<1> radiance_tl);

};
  



