// tcrad_overlap.cpp - Compute cloud overlap quantities
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

#ifndef TCRAD_OVERLAP_H
#define TCRAD_OVERLAP_H 1

#include "tcrad_adept.hpp"

namespace tcrad {

//---------------------------------------------------------------------
// Compute the upward and downward overlap matrices u_overlap and
// v_overlap, respectively, where u_overlap is defined as in Hogan et
// al. (JAS, 2016) such that y=u_overlap*x, where x is a vector of
// upwelling fluxes in each region just below an interface, and y is a
// vector of upwelling fluxes in each region just above that
// interface. HOWEVER in translating from Fortran to C++, all
// dimensions were reversed so actually the orientation of the matrix
// is such that y=transpose(u_overlap)*x. For nlev model levels there
// are nlev+1 interfaces including the ground and top-of-atmosphere,
// and so that is one of the dimensions of u_overlap and v_overlap.
template <bool IsActive>
void
calc_overlap_matrices(int nlev,
		      int ncol,
		      const Array<3,IsActive>& region_fracs,
		      const Array<2,IsActive>& overlap_param,
		      Array<4,IsActive> u_overlap,
		      Array<4,IsActive> v_overlap,
		      Real cloud_fraction_threshold);

};

#endif
