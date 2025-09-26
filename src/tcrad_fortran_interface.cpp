//#include <adept_source.h>
#include "tcrad_adept.hpp"
#include "tcrad_config.hpp"
#include "tcrad_fortran_interface.hpp"
#include "tcrad_ecckd.hpp"
#include "tcrad_regions.hpp"
#include "tcrad_overlap.hpp"
#include "tcrad_tripleclouds_lw.hpp"
#include "tcrad_two_stream.hpp"
#include "tcrad_radiance_propagation.hpp"

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
  using tcrad::NREGIONS;

  Matrix surf_emission, surf_albedo;
  Array3 planck_hl;
  Matrix cloud_fraction, fractional_std;
  Array3 od_clear, od_cloud, ssa_cloud, asymmetry_cloud;
  Matrix overlap_param;
  Array3 flux_up, flux_dn;

  Array3 reg_fracs(ncol,nlev,NREGIONS);
  Array3 od_scaling(ncol,nlev,NREGIONS);
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
  tcrad::calc_region_properties(nlev, NREGIONS, ncol, cloud_fraction, fractional_std,
				reg_fracs, od_scaling, cloud_fraction_threshold);

  Array4 u_matrix(ncol,nlev+1,NREGIONS,NREGIONS);
  Array4 v_matrix(ncol,nlev+1,NREGIONS,NREGIONS);
  tcrad::calc_overlap_matrices(nlev, ncol, reg_fracs, overlap_param,
			       u_matrix, v_matrix, cloud_fraction_threshold);

  tcrad::Config config;
  Array3 flux_up_base(nlev,NREGIONS,ng);
  Array3 flux_dn_base(nlev,NREGIONS,ng);
  Array3 flux_up_top(nlev,NREGIONS,ng);
  Array3 flux_dn_top(nlev,NREGIONS,ng);

  // Merged clear+cloudy optical properties, noting that the asymmetry
  // factor in the cloudy region equals that of the cloud since there
  // is no air/aerosol scattering
  Array3 od(nlev,NREGIONS,ng);
  Array3 ssa(nlev,NREGIONS,ng);
  
  for (int jcol = 0; jcol < ncol; ++jcol) {
    
    od(__,0,__) = od_clear[jcol];
    for (int jlev = 0; jlev < nlev; ++jlev) {
      if (reg_fracs(jcol,jlev,0) < 1.0) {
	// Cloud present
	for (int jreg = 1; jreg < NREGIONS; ++jreg) {
	  od(jlev,jreg,__) = od_clear(jcol,jlev,__) + od_scaling(jcol,jlev,jreg) * od_cloud(jcol,jlev,__);
	  ssa(jlev,jreg,__) = ssa_cloud(jcol,jlev,__) * od_cloud(jcol,jlev,__) * od_scaling(jcol,jlev,jreg)
	    / od(jlev,jreg,__);
	}
      }
    }
    tcrad::solver_tripleclouds_lw(ng, nlev, config, reg_fracs[jcol], u_matrix[jcol], v_matrix[jcol],
				  od, ssa, asymmetry_cloud[jcol],
				  planck_hl[jcol], surf_emission[jcol], surf_albedo[jcol],
				  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top);
    flux_up(jcol,range(0,nlev-1),__) = sum(flux_up_top,1);
    flux_up(jcol,nlev,__)            = sum(flux_up_base(nlev-1,__,__),0);
    flux_dn(jcol,range(0,nlev-1),__) = sum(flux_dn_top,1);
    flux_dn(jcol,nlev,__)            = sum(flux_dn_base(nlev-1,__,__),0);
  }
}

void calc_tc_radiance(int ng,
		      int nlev,
		      int ncol,
		      FortranArray* mu_in, // Cosine of sensor zenith angle
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
		      FortranArray* radiance_out)
{
  using tcrad::NREGIONS;

  Vector mu;
  Matrix surf_emission, surf_albedo;
  Array3 planck_hl;
  Matrix cloud_fraction, fractional_std;
  Array3 od_clear, od_cloud, ssa_cloud, asymmetry_cloud;
  Matrix overlap_param;
  Matrix radiance;

  Array3 reg_fracs(ncol,nlev,NREGIONS);
  Array3 od_scaling(ncol,nlev,NREGIONS);
  Real cloud_fraction_threshold = 1.0e-6;

  mu             >>= mu_in;
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
  /*
  std::cerr << "min cloud_frac = " << minval(cloud_fraction) << "\n";
  std::cerr << "max cloud_frac = " << maxval(cloud_fraction) << "\n";
  
  std::cerr << "min fsd = " << minval(fractional_std) << "\n";
  std::cerr << "max fsd = " << maxval(fractional_std) << "\n";
  */
  tcrad::calc_region_properties(nlev, NREGIONS, ncol, cloud_fraction, fractional_std,
				reg_fracs, od_scaling, cloud_fraction_threshold);

  Array4 u_matrix(ncol,nlev+1,NREGIONS,NREGIONS);
  Array4 v_matrix(ncol,nlev+1,NREGIONS,NREGIONS);
  tcrad::calc_overlap_matrices(nlev, ncol, reg_fracs, overlap_param,
			       u_matrix, v_matrix, cloud_fraction_threshold);

  tcrad::Config config;
  Array3 flux_up_base(nlev,NREGIONS,ng);
  Array3 flux_dn_base(nlev,NREGIONS,ng);
  Array3 flux_up_top(nlev,NREGIONS,ng);
  Array3 flux_dn_top(nlev,NREGIONS,ng);

  // Merged clear+cloudy optical properties, noting that the asymmetry
  // factor in the cloudy region equals that of the cloud since there
  // is no air/aerosol scattering
  Array3 od(nlev,NREGIONS,ng);
  Array3 ssa(nlev,NREGIONS,ng);

  // Transmittance along the path to the satellite
  Array3 transmittance(nlev,NREGIONS,ng);
  Array3 source_dn(nlev,NREGIONS,ng);
  Array3 source_up(nlev,NREGIONS,ng);
  
  for (int jcol = 0; jcol < ncol; ++jcol) {
    
    od(__,0,__) = od_clear[jcol];
    for (int jlev = 0; jlev < nlev; ++jlev) {
      if (reg_fracs(jcol,jlev,0) < 1.0) {
	// Cloud present
	for (int jreg = 1; jreg < NREGIONS; ++jreg) {
	  od(jlev,jreg,__) = od_clear(jcol,jlev,__) + od_scaling(jcol,jlev,jreg) * od_cloud(jcol,jlev,__);
	  ssa(jlev,jreg,__) = ssa_cloud(jcol,jlev,__) * od_cloud(jcol,jlev,__) * od_scaling(jcol,jlev,jreg)
	    / od(jlev,jreg,__);
	}
      }
    }
    tcrad::solver_tripleclouds_lw(ng, nlev, config, reg_fracs[jcol], u_matrix[jcol], v_matrix[jcol],
				  od, ssa, asymmetry_cloud[jcol],
				  planck_hl[jcol], surf_emission[jcol], surf_albedo[jcol],
				  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top);
    tcrad::calc_transmittance(ng, nlev, mu(jcol), reg_fracs(jcol,__,0),
			      od, transmittance);
    Matrix radiance_up_surf(NREGIONS,ng);
    tcrad::calc_radiance_source(ng, nlev, config, mu(jcol), reg_fracs[jcol], planck_hl[jcol],
				od, ssa, asymmetry_cloud[jcol], flux_up_base, flux_dn_top,
				transmittance, source_dn, source_up);
    if (config.do_specular_surface) {
      Matrix radiance_dn_surf(NREGIONS,ng);
      tcrad::calc_radiance_dn(ng, nlev, reg_fracs(jcol,__,0), transmittance,
			      source_dn, v_matrix[jcol], radiance_dn_surf);
      for (int jreg = 0; jreg < NREGIONS; ++jreg) {
	radiance_up_surf(jreg,__) = reg_fracs(jcol,nlev-1,jreg)*surf_emission[jcol]/tcrad::PI
	  + surf_albedo[jcol]*radiance_dn_surf(jreg,__);
      }
    }
    else {
      radiance_up_surf = flux_up_base(nlev-1,__,__)/tcrad::PI;
    }
    tcrad::calc_radiance_up(ng, nlev, reg_fracs(jcol,__,0), radiance_up_surf, transmittance, source_up,
			    u_matrix[jcol], radiance[jcol]);
  }
}

};


