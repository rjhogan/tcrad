// tcrad_tripleclouds_lw.hpp - The longwave Tripleclouds flux solver
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

#ifndef TCRAD_TRIPLECLOUDS_LW_H
#define TCRAD_TRIPLECLOUDS_LW_H 1

#include "tcrad_adept.hpp"
#include "tcrad_config.hpp"

namespace tcrad {

// ---------------------------------------------------------------------
// Longwave "Tripleclouds" solver following Shonk & Hogan (2008) that
// treats cloud inhomogeneity by dividing each model level into three
// regions, one clear and two cloudy (with differing optical
// depth). This version computes the upwelling and downwelling fluxes
// at the top and bottom of each layer and region, and can used in
// subsequent processing for computing either broadband flux profiles
// or top-of-atmosphere radiances.
template <bool IsActive>
void
solver_tripleclouds_lw(int ng,
		       int nlev,
		       const Config& config,
		       const Array<2,IsActive>& region_fracs, // (nlev,nreg)
		       const Array<3,IsActive>& u_overlap,    // (nlev+1,nreg,nreg)
		       const Array<3,IsActive>& v_overlap,    // (nlev+1,nreg,nreg)
		       const Array<3,IsActive>& od,           // (nlev,nreg,ng)
		       const Array<3,IsActive>& ssa,          // (nlev,nreg,ng)
		       const Array<2,IsActive>& asymmetry,    // (nlev,ng)
		       const Array<2,IsActive>& planck_hl,    // (nlev+1,ng)
		       const Array<1,IsActive>& surf_emission,// (ng)
		       const Array<1,IsActive>& surf_albedo,  // (ng)
		       Array<3,IsActive> flux_up_base,        // (nlev,nreg,ng)
		       Array<3,IsActive> flux_dn_base,        // (nlev,nreg,ng)
		       Array<3,IsActive> flux_up_top,         // (nlev,nreg,ng)
		       Array<3,IsActive> flux_dn_top);        // (nlev,nreg,ng)

}; // namespace tcrad

#endif
