#ifndef ECRAD_TRIPLECLOUDS_LW_H
#define ECRAD_TRIPLECLOUDS_LW_H 1

#include "tcrad_adept.hpp"
#include "tcrad_config.hpp"

namespace tcrad {

template <bool IsActive>
void
solver_tripleclouds_lw(int ng,
		       int nlev,
		       int ncol,
		       const Config& config,
		       const aArray<3,IsActive>& region_fracs, // (ncol,nlev,nreg)
		       const aArray<3,IsActive>& od_scaling,   // (ncol,nlev,nreg-1)
		       const aArray<4,IsActive>& u_matrix,     // (ncol,nlev+1,nreg,nreg)
		       const aArray<4,IsActive>& v_matrix,     // (ncol,nlev+1,nreg,nreg)
		       const aArray<3,IsActive>& od,           // (ncol,nlev,ng)
		       const aArray<3,IsActive>& od_cloud,     // (ncol,nlev,ng)
		       const aArray<3,IsActive>& ssa_cloud,    // (ncol,nlev,ng)
		       const aArray<3,IsActive>& asymmetry_cloud, // (ncol,nlev,ng)
		       const aArray<3,IsActive>& planck_hl,    // (ncol,nlev+1,ng)
		       const aMatrix<IsActive>& surf_emission, // (ncol,ng)
		       const aMatrix<IsActive>& surf_albedo,   // (ncol,ng)
		       /*
		       aMatrix<IsActive> lw_up,                // (ncol,nlev+1)
		       aMatrix<IsActive> lw_dn);               // (ncol,nlev+1)
		       */
		       aArray<3,IsActive> lw_up,                // (ncol,nlev+1)
		       aArray<3,IsActive> lw_dn);               // (ncol,nlev+1)
  

}; // namespace tcrad
#endif
