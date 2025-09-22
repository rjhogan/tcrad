#ifndef ECRAD_TRIPLECLOUDS_LW_H
#define ECRAD_TRIPLECLOUDS_LW_H 1

#include "tcrad_adept.hpp"
#include "tcrad_config.hpp"

namespace tcrad {

template <bool IsActive>
void
solver_tripleclouds_lw(int ng,
		       int nlev,
		       const Config& config,
		       const aArray<2,IsActive>& region_fracs, // (nlev,nreg)
		       const aArray<3,IsActive>& u_matrix,     // (nlev+1,nreg,nreg)
		       const aArray<3,IsActive>& v_matrix,     // (nlev+1,nreg,nreg)
		       const aArray<3,IsActive>& od,           // (nlev,nreg,ng)
		       const aArray<3,IsActive>& ssa,          // (nlev,nreg,ng)
		       const aArray<2,IsActive>& asymmetry,    // (nlev,ng)
		       const aArray<2,IsActive>& planck_hl,    // (nlev+1,ng)
		       const aArray<1,IsActive>& surf_emission,// (ng)
		       const aArray<1,IsActive>& surf_albedo,  // (ng)
		       aArray<3,IsActive> flux_up_base,        // (nlev,nreg,ng)
		       aArray<3,IsActive> flux_dn_base,        // (nlev,nreg,ng)
		       aArray<3,IsActive> flux_up_top,         // (nlev,nreg,ng)
		       aArray<3,IsActive> flux_dn_top);        // (nlev,nreg,ng)

}; // namespace tcrad
#endif
