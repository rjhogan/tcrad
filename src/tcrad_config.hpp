#ifndef ECRAD_CONFIG_H
#define ECRAD_CONFIG_H 1

#include "tcrad_adept.hpp"
#include "tcrad_ecckd.hpp"

namespace tcrad {

  // Triplecluds, the clue is in the name
  static const int NREGIONS = 3;
  static const Real PI = acos(-1.0);
  static const Real LW_DIFFUSIVITY = 1.66;
  static const Real LW_INV_DIFFUSIVITY = (1.0/LW_DIFFUSIVITY);
  // To avoid division by near-zero values use simpler formulae in the
  // low optical depth regime
  static const Real OD_THRESH = 1.0e-3;

  typedef enum {
    TWO_STREAM_ELSASSER = 0,
    TWO_STREAM_EDDINGTON,
    TWO_STREAM_LEGENDRE,
    TWO_STREAM_HYBRID,
    TWO_STREAM_SCALED_WISCOMBE_GRAMS
  } TwoStreamScheme;
    
  class Config {
  public:
    bool do_specular_surface = true;
    TwoStreamScheme i_two_stream_scheme = TWO_STREAM_ELSASSER;
    CkdModel gas_optics_sw;
    CkdModel gas_optics_lw;
    Real cloud_fraction_threshold = 1.0e-6;
  };
  
};

#endif
