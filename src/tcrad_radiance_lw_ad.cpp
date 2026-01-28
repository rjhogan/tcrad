// tcrad_radiance_lw_ad.cpp - Adjoint of longwave radiance model
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
// Longwave "Tripleclouds" solver following Shonk & Hogan (2008) that
// treats cloud inhomogeneity by dividing each model level into three
// regions, one clear and two cloudy (with differing optical depth),
// returning a top-of-atmosphere spectral upwelling radiances. This is
// the adjoint version.
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
				      Array<1> radiance, // Output radiance
				      const Array<1>& radiance_ad, // Input adjoint
				      Array<1> surf_emission_ad, // Output adjoints...
				      Array<1> surf_albedo_ad,
				      Array<2> planck_hl_ad,
				      Array<1> cloud_fraction_ad,
				      Array<2> od_clear_ad,
				      Array<2> od_cloud_ad,
				      Array<2> ssa_cloud_ad,
				      Array<2> asymmetry_cloud_ad,
				      Array<1> overlap_param_ad) {

  // Create an Adept Stack if one does not exist in this thread
  // already, and get a reference to it (ADEPT_ACTIVE_STACK is a
  // pointer to the thread-local global variable)
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

  // Start recording differential statements
  stack.new_recording();

  aArray<1> radiance_active(ng);

  // Call the radiance model
  calc_tripleclouds_radiance_lw(ng, nlev, config, mu,
				surf_emission_active, surf_albedo_active,
				planck_hl_active, cloud_fraction_active,
				fractional_std, od_clear_active,
				od_cloud_active, ssa_cloud_active, asymmetry_cloud_active,
				overlap_param_active, radiance_active);

  // Extract inactive array from the result
  radiance = value(radiance_active);

  // Set the adjoint of the radiance to the value provided by the user
  radiance_active.set_gradient(radiance_ad);

  // Reverse pass through the differential statements
  stack.compute_adjoint();

  //std::cerr << "Memory: " << stack.memory () << " bytes " << ng << " " << nlev << "\n";
  //std::cerr << stack  << "\n";
  
  // Save the output adjoints requested by the user
  if (!surf_emission_ad.empty())
    surf_emission_ad   = surf_emission_active.get_gradient();
  if (!surf_albedo_ad.empty())
    surf_albedo_ad     = surf_albedo_active.get_gradient();
  if (!planck_hl_ad.empty())
    planck_hl_ad       = planck_hl_active.get_gradient();
  if (!cloud_fraction_ad.empty())
    cloud_fraction_ad  = cloud_fraction_active.get_gradient();
  if (!od_clear_ad.empty())
    od_clear_ad        = od_clear_active.get_gradient();
  if (!od_cloud_ad.empty())
    od_cloud_ad        = od_cloud_active.get_gradient();
  if (!ssa_cloud_ad.empty())
    ssa_cloud_ad       = ssa_cloud_active.get_gradient();
  if (!asymmetry_cloud_ad.empty())
    asymmetry_cloud_ad = asymmetry_cloud_active.get_gradient();
  if (!overlap_param_ad.empty())
    overlap_param_ad   = overlap_param_active.get_gradient();
  
}
  
};
  



