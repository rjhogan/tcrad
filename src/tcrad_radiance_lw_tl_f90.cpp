// tcrad_radiance_lw_tl_f90.cpp - Fortran-90 interface to tangent-linear of radiance computation
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
#include "tcrad_config.hpp"
#include "tcrad_radiance_lw_tl.hpp"

extern "C" {

// ---------------------------------------------------------------------
// Longwave "Tripleclouds" solver following Shonk & Hogan (2008) that
// treats cloud inhomogeneity by dividing each model level into three
// regions, one clear and two cloudy (with differing optical depth),
// returning a top-of-atmosphere spectral upwelling radiances. This
// routine provides the Fortran-90 interface using interoperability
// features provided by the Adept library.
void tcrad_radiance_lw_tl(int ng,
			  int nlev,
			  tcrad::Config* config,
			  Real mu, // Cosine of sensor zenith angle
			  FortranArray* surf_emission_in,
			  FortranArray* surf_albedo_in,
			  FortranArray* planck_hl_in,
			  FortranArray* cloud_fraction_in,
			  FortranArray* fractional_std_in,
			  FortranArray* od_clear_in,
			  FortranArray* od_cloud_in,
			  FortranArray* ssa_cloud_in,
			  FortranArray* asymmetry_cloud_in,
			  FortranArray* overlap_param_in,
			  FortranArray* radiance_out,
			  FortranArray* surf_emission_tl_in,
			  FortranArray* surf_albedo_tl_in,
			  FortranArray* planck_hl_tl_in,
			  FortranArray* cloud_fraction_tl_in,
			  FortranArray* od_clear_tl_in,
			  FortranArray* od_cloud_tl_in,
			  FortranArray* ssa_cloud_tl_in,
			  FortranArray* asymmetry_cloud_tl_in,
			  FortranArray* overlap_param_tl_in,
			  FortranArray* radiance_tl_out)
{
  Array<1> surf_emission, surf_albedo;
  Array<2> planck_hl;
  Array<1> cloud_fraction, fractional_std;
  Array<2> od_clear, od_cloud, ssa_cloud, asymmetry_cloud;
  Array<1> overlap_param;
  Array<1> radiance;
  
  Array<1> surf_emission_tl, surf_albedo_tl;
  Array<2> planck_hl_tl;
  Array<1> cloud_fraction_tl;
  Array<2> od_clear_tl, od_cloud_tl, ssa_cloud_tl, asymmetry_cloud_tl;
  Array<1> overlap_param_tl;
  Array<1> radiance_tl;
  

  surf_emission  >>= surf_emission_in;
  surf_albedo    >>= surf_albedo_in;
  planck_hl      >>= planck_hl_in;
  cloud_fraction >>= cloud_fraction_in;
  fractional_std >>= fractional_std_in;
  od_clear       >>= od_clear_in;
  od_cloud       >>= od_cloud_in;
  ssa_cloud      >>= ssa_cloud_in;
  asymmetry_cloud>>= asymmetry_cloud_in;
  overlap_param  >>= overlap_param_in;
  radiance       >>= radiance_out;

  surf_emission_tl  >>= surf_emission_tl_in;
  surf_albedo_tl    >>= surf_albedo_tl_in;
  planck_hl_tl      >>= planck_hl_tl_in;
  cloud_fraction_tl >>= cloud_fraction_tl_in;
  od_clear_tl       >>= od_clear_tl_in;
  od_cloud_tl       >>= od_cloud_tl_in;
  ssa_cloud_tl      >>= ssa_cloud_tl_in;
  asymmetry_cloud_tl>>= asymmetry_cloud_tl_in;
  overlap_param_tl  >>= overlap_param_tl_in;
  radiance_tl       >>= radiance_tl_out;

  calc_tripleclouds_radiance_lw_tl(ng, nlev, *config, mu,
				   surf_emission, surf_albedo,
				   planck_hl, cloud_fraction,
				   fractional_std, od_clear,
				   od_cloud, ssa_cloud, asymmetry_cloud,
				   overlap_param, radiance,
				   surf_emission_tl, surf_albedo_tl,
				   planck_hl_tl, cloud_fraction_tl,
				   od_clear_tl, od_cloud_tl, ssa_cloud_tl,
				   asymmetry_cloud_tl, overlap_param_tl,
				   radiance_tl);
}
 
};


