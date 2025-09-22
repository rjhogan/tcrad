//#include <adept_source.h>
#include "tcrad_adept.hpp"
#include "tcrad_config.hpp"
#include "tcrad_fortran_interface.hpp"
#include "tcrad_ecckd.hpp"
#include "tcrad_regions.hpp"
#include "tcrad_overlap.hpp"
#include "tcrad_tripleclouds_lw.hpp"

#include <fenv.h>
int _feenableexcept_status = feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);

extern "C" {
		
void* allocate_tcrad() {
  return static_cast<void*>(new tcrad::Config());
}

void deallocate_tcrad(tcrad::Config* config) {
  delete config;
}

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
		       FortranArray* rayleigh_molar_scat) {
  tcrad::CkdModel& ckd_model = config->gas_optics_sw;
  ckd_model.ngas = ngas;
  for (int igas = 0; igas < tcrad::MAX_GASES; ++igas) {
    ckd_model.gas_mapping[igas] = gas_mapping[igas];
  }
  ckd_model.npress = npress;
  ckd_model.ntemp  = ntemp;
  ckd_model.log_pressure1 = log_pressure1;
  ckd_model.d_log_pressure = d_log_pressure;
  ckd_model.temperature1 >>= temperature1;
  ckd_model.d_temperature = d_temperature;
  ckd_model.norm_solar_irradiance >>= norm_solar_irradiance;
  ckd_model.norm_amplitude_solar_irradiance >>= norm_amplitude_solar_irradiance;
  ckd_model.rayleigh_molar_scat >>= rayleigh_molar_scat;
}

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
		       FortranArray* planck_function) {
  tcrad::CkdModel& ckd_model = config->gas_optics_lw;
  ckd_model.ngas = ngas;
  for (int igas = 0; igas < tcrad::MAX_GASES; ++igas) {
    ckd_model.gas_mapping[igas] = gas_mapping[igas];
  }
  ckd_model.npress = npress;
  ckd_model.ntemp  = ntemp;
  ckd_model.log_pressure1 = log_pressure1;
  ckd_model.d_log_pressure = d_log_pressure;
  ckd_model.temperature1 >>= temperature1;
  ckd_model.d_temperature = d_temperature;
  ckd_model.nplanck = nplanck;
  ckd_model.temperature1_planck = temperature1_planck;
  ckd_model.planck_function >>= planck_function;
}		    

void calc_tc_flux(int ng,
		  int nlev,
		  int ncol,
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
		  FortranArray* flux_up_out,
		  FortranArray* flux_dn_out)
{

  Matrix surf_emission, surf_albedo;
  Array3 planck_hl;
  Matrix cloud_fraction, fractional_std;
  Array3 od_clear, od_cloud, ssa_cloud, asymmetry_cloud;
  Matrix overlap_param;
  Array3 flux_up, flux_dn;

  Array3 reg_fracs(ncol,nlev,3);
  Array3 od_scaling(ncol,nlev,3);
  Real cloud_fraction_threshold = 1.0e-6;

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
  flux_up        >>= flux_up_out;
  flux_dn        >>= flux_dn_out;
  /*
  std::cerr << "min cloud_frac = " << minval(cloud_fraction) << "\n";
  std::cerr << "max cloud_frac = " << maxval(cloud_fraction) << "\n";
  
  std::cerr << "min fsd = " << minval(fractional_std) << "\n";
  std::cerr << "max fsd = " << maxval(fractional_std) << "\n";
  */
  tcrad::calc_region_properties(nlev, 3, ncol, cloud_fraction, fractional_std,
				reg_fracs, od_scaling, cloud_fraction_threshold);

  Array4 u_matrix(ncol,nlev+1,3,3);
  Array4 v_matrix(ncol,nlev+1,3,3);
  tcrad::calc_overlap_matrices(nlev, ncol, reg_fracs, overlap_param,
			       u_matrix, v_matrix, cloud_fraction_threshold);

  tcrad::Config config;
  Array3 flux_up_base(nlev,3,ng);
  Array3 flux_dn_base(nlev,3,ng);
  Array3 flux_up_top(nlev,3,ng);
  Array3 flux_dn_top(nlev,3,ng);
  for (int jcol = 0; jcol < ncol; ++jcol) {
    tcrad::solver_tripleclouds_lw(ng, nlev, config, reg_fracs[jcol], od_scaling[jcol],
				  u_matrix[jcol], v_matrix[jcol], od_clear[jcol], od_cloud[jcol], ssa_cloud[jcol],
				  asymmetry_cloud[jcol], planck_hl[jcol], surf_emission[jcol], surf_albedo[jcol],
				  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top);
    flux_up(jcol,range(0,nlev-1),__) = sum(flux_up_top,1);
    flux_up(jcol,nlev,__)            = sum(flux_up_base(nlev,__,__),1);
    flux_dn(jcol,range(0,nlev-1),__) = sum(flux_dn_top,1);
    flux_dn(jcol,nlev,__)            = sum(flux_dn_base(nlev,__,__),1);
  }
}
};
