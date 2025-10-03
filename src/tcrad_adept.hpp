// tcrad_adept.cpp - Provide Adept array & automatic differentiation 
//
// Copyright (C) 20205 ECMWF.
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

#ifndef TCRAD_ADEPT_H
#define TCRAD_ADEPT_H 1

#include <algorithm>

// Adept array functionality and Fortran interoperability
#include <adept_arrays.h>
#include <adept_fortran.h>

using std::min;
using std::max;

// Import basic Adept features
using adept::__;
using adept::end;
using adept::range;
using adept::intVector;
using adept::boolVector;
using adept::FortranArray;

// Working floating-point type, "double" by default but if you compile
// with ADEPT_REAL_TYPE_SIZE=4 then it will be "float", although in
// that case automatic differentiation is not possible.
using adept::Real;

// Create a floating-point array type with the specified number of
// dimensions and whether it is "active" in computing derivatives,
// e.g. Array<2> is an ordinary matrix and Array<1,true> is an active
// vector.
template <int NDims, bool IsActive = false>
using Array = adept::Array<NDims,Real,IsActive>;

// Fixed size arrays
template <bool IsActive>
using Vector2 = adept::FixedArray<Real,IsActive,2>;
template <bool IsActive>
using Vector3 = adept::FixedArray<Real,IsActive,3>;
template <bool IsActive>
using Matrix33 = adept::FixedArray<Real,IsActive,3,3>;

// Select an active versus ordinary scalar depending on bool template
// argument
template <bool IsActive>
struct scalar {
  typedef adept::Real type;
};
template <>
struct scalar<true> {
  typedef adept::aReal type;
};

#endif
