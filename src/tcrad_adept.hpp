#ifndef ECRAD_ADEPT_H
#define ECRAD_ADEPT_H 1

//#include <stdbool.h>
#include <adept_arrays.h>
//#include "/home/parr/git/Adept-2/include/adept_fortran.h"
#include <adept_fortran.h>

using std::min;
using std::max;

using adept::__;
using adept::end;
using adept::Real;
using adept::Vector;
using adept::Matrix;
using adept::Array3;
typedef adept::Array<4,Real> Array4;
using adept::intVector;
using adept::boolVector;
using adept::FortranArray;
using adept::range;

template <bool IsActive>
using aVector = adept::Array<1,Real,IsActive>;

template <bool IsActive>
using aMatrix = adept::Array<2,Real,IsActive>;

template <int NDims>
using Array = adept::Array<NDims,Real,false>;

template <int NDims, bool IsActive>
using aArray = adept::Array<NDims,Real,IsActive>;

template <bool IsActive>
using aVector2 = adept::FixedArray<Real,IsActive,2>;

template <bool IsActive>
using aVector3 = adept::FixedArray<Real,IsActive,3>;

template <bool IsActive>
using aMatrix33 = adept::FixedArray<Real,IsActive,3,3>;

/// Select an active versus static scalar depending on bool template
/// argument
template <bool IsActive>
struct scalar {
  typedef adept::Real type;
};

template <>
struct scalar<true> {
  typedef adept::aReal type;
};

#endif
