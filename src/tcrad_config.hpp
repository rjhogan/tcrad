// tcrad_config.cpp - Define configuration structure
//
// Copyright (C) 2025 ECMWF.
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

#ifndef TCRAD_CONFIG_H
#define TCRAD_CONFIG_H 1

#include "tcrad_adept.hpp"

namespace tcrad {

  static const Real PI = acos(-1.0);

  // Triplecluds, the clue is in the name
  static const int NREGIONS = 3;

  // Cloud heterogeneities decorrelate twice as fast as cloud
  // boundaries
  static const Real DECORRELATION_SCALING = 0.5;

  // Various two-stream schemes
  enum TwoStreamScheme : int {
    TWO_STREAM_ELSASSER = 0,
    TWO_STREAM_EDDINGTON = 1,
    TWO_STREAM_LEGENDRE = 2,
    TWO_STREAM_HYBRID = 3,
    TWO_STREAM_SCALED_WISCOMBE_GRAMS = 4
  };

  // Configuration structure
  class Config {
  public:
    // Computing satellite microwave radiances over the sea typically
    // treats the surface as specular, but this is more expensive as
    // both downward and upward radiance calculations are needed
    bool do_specular_surface = false;
    TwoStreamScheme i_two_stream_scheme = TWO_STREAM_ELSASSER;
    Real cloud_fraction_threshold = 1.0e-6;
  };
  
};

#endif
