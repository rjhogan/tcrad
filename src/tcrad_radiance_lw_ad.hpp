// tcrad_radiance_lw_ad.hpp - Adjoint of longwave radiance model
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
// Longwave "Tripleclouds" solver following Shonk & Hogan (2008) that
// treats cloud inhomogeneity by dividing each model level into three
// regions, one clear and two cloudy (with differing optical depth),
// returning a top-of-atmosphere spectral upwelling radiances.
void calc_tripleclouds_radiance_lw_ad(int ng,
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
				      const Array<1>& radiance_ad,
				      Array<1> surf_emission_ad,
				      Array<1> surf_albedo_ad,
				      Array<2> planck_hl_ad,
				      Array<1> cloud_fraction_ad,
				      Array<2> od_clear_ad,
				      Array<2> od_cloud_ad,
				      Array<2> ssa_cloud_ad,
				      Array<2> asymmetry_cloud_ad,
				      Array<1> overlap_param_ad);

};


