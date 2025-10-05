// tcrad_radiance_lw_tl.cpp - Tangent-linear of longwave radiance model
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
#include "tcrad_radiance_lw.hpp"

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
				      Array<1> radiance_tl) {

  // Get the currenty active stack on this thread from global variable
  tcrad_create_autodiff_engine();
  Stack& stack = *ADEPT_ACTIVE_STACK;

  // Create active versions of each input array that link to the
  // existing data; note that this does not conserve constness, so
  // don't change the contents of these arrays.
  aArray<1> surf_emission_active(surf_emission);
  aArray<1> surf_albedo_active(surf_albedo);
  aArray<2> planck_hl_active(planck_hl);
  aArray<1> cloud_fraction_active(cloud_fraction);
  aArray<2> od_clear_active(od_clear);
  aArray<2> od_cloud_active(od_cloud);
  aArray<2> ssa_cloud_active(ssa_cloud);
  aArray<2> asymmetry_cloud_active(asymmetry_cloud);
  aArray<1> overlap_param_active(overlap_param);

  // Start a new recording
  stack.new_recording();

  aArray<1> radiance_active(ng);
  
  // Run the algorithm
  calc_tripleclouds_radiance_lw(ng, nlev, config, mu,
				surf_emission_active, surf_albedo_active,
				planck_hl_active, cloud_fraction_active,
				fractional_std, od_clear_active,
				od_cloud_active, ssa_cloud_active, asymmetry_cloud_active,
				overlap_param_active, radiance_active);

  // Extract inactive array from the result
  radiance = value(radiance_active);

  // Set the tangent-linear inputs
  if (!surf_emission_tl.empty())
    surf_emission_active.set_gradient(surf_emission_tl);
  if (!surf_albedo_tl.empty())
    surf_albedo_active.set_gradient(surf_albedo_tl);
  if (!planck_hl_tl.empty())
    planck_hl_active.set_gradient(planck_hl_tl);
  if (!cloud_fraction_tl.empty())
    cloud_fraction_active.set_gradient(cloud_fraction_tl);
  if (!od_clear_tl.empty()) 
    od_clear_active.set_gradient(od_clear_tl);
  if (!od_cloud_tl.empty())
    od_cloud_active.set_gradient(od_cloud_tl);
  if (!ssa_cloud_tl.empty())
    ssa_cloud_active.set_gradient(ssa_cloud_tl);
  if (!asymmetry_cloud_tl.empty())
    asymmetry_cloud_active.set_gradient(asymmetry_cloud_tl);
  if (!overlap_param_tl.empty())
    overlap_param_active.set_gradient(overlap_param_tl);

  // Run forward through the derivative calculation
  stack.compute_tangent_linear();

  // Extract the tangent-linear of the result
  radiance_tl = radiance_active.get_gradient();
}
  
};
  



