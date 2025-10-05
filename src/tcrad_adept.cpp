// tcrad_adept.cpp - Activate Adept's auto-differentiation stack
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

extern "C" {

  // Create an Adept Stack for automatic differentiation in the
  // current thread if one does not exist already; the constructor of
  // the Stack object stores a pointer to itself in a thread-local
  // global variable accessible through ADEPT_ACTIVE_STACK.
  void tcrad_create_autodiff_engine() {
    if (ADEPT_ACTIVE_STACK == 0) {
      Stack* stack = new Stack;
    }
  }
  // Delete an Adept Stack created by "new"; this is likely to fail if
  // the Stack was not created using the create function above.
  void tcrad_delete_autodiff_engine() {
    if (ADEPT_ACTIVE_STACK == 0) {
      delete ADEPT_ACTIVE_STACK;
    }
  }
  // Return whether an Adept automatic differentiation engine exists
  // for this thread.
  bool tcrad_autodiff_engine_exists() {
    if (ADEPT_ACTIVE_STACK != 0) {
      return true;
    }
    else {
      return false;
    }
  }

};
