! -*- f90 -*-
interface

   subroutine calc_tc_flux(ng, nlev, ncol, surf_emission, surf_albedo, &
        &     planck_hl, cloud_fraction, fractional_std, &
        &     od_clear, od_cloud, ssa_cloud, asymmetry_cloud, &
        &     overlap_param, &
        &     flux_up, flux_dn) bind(c)

     use iso_c_binding
  
     !use parkind1, only : jprb
     integer, parameter :: wp = c_double
     integer(c_int), value :: ng, nlev, ncol
     real(kind=wp), intent(in) :: surf_emission(:,:), surf_albedo(:,:)
     real(kind=wp), intent(in) :: planck_hl(:,:,:)
     real(kind=wp), intent(in) :: cloud_fraction(:,:), fractional_std(:,:)
     real(kind=wp), intent(in) :: od_clear(:,:,:), od_cloud(:,:,:)
     real(kind=wp), intent(in) :: ssa_cloud(:,:,:), asymmetry_cloud(:,:,:)
     real(kind=wp), intent(in) :: overlap_param(:,:)
     real(kind=wp), intent(out) :: flux_up(:,:,:), flux_dn(:,:,:)

   end subroutine calc_tc_flux


   subroutine calc_tc_radiance(ng, nlev, ncol, mu, surf_emission, surf_albedo, &
        &     planck_hl, cloud_fraction, fractional_std, &
        &     od_clear, od_cloud, ssa_cloud, asymmetry_cloud, &
        &     overlap_param, &
        &     radiance) bind(c)

     use iso_c_binding
  
     !use parkind1, only : jprb
     integer, parameter :: wp = c_double
     integer(c_int), value :: ng, nlev, ncol
     real(kind=wp), intent(in) :: mu(:)
     real(kind=wp), intent(in) :: surf_emission(:,:), surf_albedo(:,:)
     real(kind=wp), intent(in) :: planck_hl(:,:,:)
     real(kind=wp), intent(in) :: cloud_fraction(:,:), fractional_std(:,:)
     real(kind=wp), intent(in) :: od_clear(:,:,:), od_cloud(:,:,:)
     real(kind=wp), intent(in) :: ssa_cloud(:,:,:), asymmetry_cloud(:,:,:)
     real(kind=wp), intent(in) :: overlap_param(:,:)
     real(kind=wp), intent(out) :: radiance(:,:)

   end subroutine calc_tc_radiance

   
end interface
