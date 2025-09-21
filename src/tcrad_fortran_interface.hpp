#ifndef ECRAD_FORTRAN_INTERFACE_H
#define ECRAD_FORTRAN_INTERFACE_H 1

#include "tcrad_config.hpp"
#include "tcrad_ecckd.hpp"

extern "C" {

void* allocate_tcrad();

void deallocate_tcrad(tcrad::Config* config);

void fill_ckd_model_sw(tcrad::Config* config,
		       int ngas,
		       int* gas_mapping,
		       int npress,
		       int ntemp,
		       Real log_pressure1,
		       Real d_log_pressure,
		       FortranArray* temperature1,
		       Real d_temperature,
		       FortranArray* norm_solar_irradiance,
		       FortranArray* norm_amplitude_solar_irradiance,
		       FortranArray* rayleigh_molar_scat);

void fill_ckd_model_lw(tcrad::Config* config,
		       int ngas,
		       int* gas_mapping,
		       int npress,
		       int ntemp,
		       Real log_pressure1,
		       Real d_log_pressure,
		       FortranArray* temperature1,
		       Real d_temperature,
		       int nplanck,
		       Real temperature1_planck,
		       Real d_temperature_planck,
		       FortranArray* planck_function);

}

#endif
