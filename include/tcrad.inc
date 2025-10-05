! tcrad.inc - Fortran header file for Tripleclouds radiance model (TCRAD) -*- f90 -*-
!
! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!
! Include this header file with
!   #include "tcrad.inc"
! in a Fortran routine that needs to call any TCRAD functions, but note that a
!   use iso_c_binding
! will need to be included earleir in the file.  

! Enumerators defining the two-stream scheme
enum, bind(c)
  enumerator :: TWO_STREAM_ELSASSER = 0
  enumerator :: TWO_STREAM_EDDINGTON = 1
  enumerator :: TWO_STREAM_LEGENDRE = 2
  enumerator :: TWO_STREAM_HYBRID = 3
  enumerator :: TWO_STREAM_SCALED_WISCOMBE_GRAMS = 4
end enum

! Configuration structure passed to TCRAD routines
type, bind(c) :: tcrad_config
  logical(c_bool) :: do_specular_surface = .false.
  integer(c_int)  :: two_stream_scheme = TWO_STREAM_ELSASSER
  real(c_double)  :: cloud_fraction_threshold = 1.0e-6
end type tcrad_config
  
interface

  ! ---------------------------------------------------------------------
  ! Create an automatic differentiation engine in the current thread,
  ! if one does not exist already, enabling TCRAD tangent-linear and
  ! adjoints to work. This is called automatically by the TL and AD
  ! functions.
  subroutine tcrad_create_autodiff_engine() bind(c)
  end subroutine tcrad_create_autodiff_engine

  ! ---------------------------------------------------------------------
  ! Delete the automatic differentiation engine in the current thread,
  ! if one exists. The engine saves memory for future calls to save
  ! reallocation, so it may be useful to explicitly delete when no
  ! longer needed.
  subroutine tcrad_delete_autodiff_engine() bind(c)
  end subroutine tcrad_delete_autodiff_engine

  ! ---------------------------------------------------------------------
  ! Returns true if an automatic differentiation engine exists in the
  ! current thread.

  function tcrad_autodiff_engine_exists() bind(c)
    use iso_c_binding
    logical(c_bool) :: tcrad_autodiff_engine_exists
  end function tcrad_autodiff_engine_exists
  
  ! ---------------------------------------------------------------------
  ! Longwave "Tripleclouds" solver following Shonk & Hogan (2008) that
  ! treats cloud inhomogeneity by dividing each model level into three
  ! regions, one clear and two cloudy (with differing optical depth),
  ! returning a profile of upwelling and downwelling spectral fluxes.
  subroutine tcrad_flux_lw(ng, nlev, config, &
       &     surf_emission, surf_albedo, &
       &     planck_hl, cloud_fraction, fractional_std, &
       &     od_clear, od_cloud, ssa_cloud, asymmetry_cloud, &
       &     overlap_param, &
       &     flux_up, flux_dn) bind(c)

    use iso_c_binding

    import :: tcrad_config

    integer(c_int),   value      :: ng, nlev
    type(tcrad_config),intent(in):: config
    real(c_double),   intent(in) :: surf_emission(:), surf_albedo(:)
    real(c_double),   intent(in) :: planck_hl(:,:)
    real(c_double),   intent(in) :: cloud_fraction(:), fractional_std(:)
    real(c_double),   intent(in) :: od_clear(:,:), od_cloud(:,:)
    real(c_double),   intent(in) :: ssa_cloud(:,:), asymmetry_cloud(:,:)
    real(c_double),   intent(in) :: overlap_param(:)
    real(c_double),   intent(out):: flux_up(:,:), flux_dn(:,:)
    
  end subroutine tcrad_flux_lw
  
  ! ---------------------------------------------------------------------
  ! Longwave "Tripleclouds" solver following Shonk & Hogan (2008) that
  ! treats cloud inhomogeneity by dividing each model level into three
  ! regions, one clear and two cloudy (with differing optical depth),
  ! returning a top-of-atmosphere spectral upwelling radiances.
  subroutine tcrad_radiance_lw(ng, nlev, config, &
       &     mu, surf_emission, surf_albedo, &
       &     planck_hl, cloud_fraction, fractional_std, &
       &     od_clear, od_cloud, ssa_cloud, asymmetry_cloud, &
       &     overlap_param, &
       &     radiance) bind(c)

    use iso_c_binding

    import :: tcrad_config
    
    integer(c_int), value        :: ng, nlev
    type(tcrad_config),intent(in):: config
    real(c_double), value        :: mu
    real(c_double), intent(in)   :: surf_emission(:), surf_albedo(:)
    real(c_double), intent(in)   :: planck_hl(:,:)
    real(c_double), intent(in)   :: cloud_fraction(:), fractional_std(:)
    real(c_double), intent(in)   :: od_clear(:,:), od_cloud(:,:)
    real(c_double), intent(in)   :: ssa_cloud(:,:), asymmetry_cloud(:,:)
    real(c_double), intent(in)   :: overlap_param(:)
    real(c_double), intent(out)  :: radiance(:)
    
  end subroutine tcrad_radiance_lw
    
  ! ---------------------------------------------------------------------
  ! Longwave "Tripleclouds" solver following Shonk & Hogan (2008) that
  ! treats cloud inhomogeneity by dividing each model level into three
  ! regions, one clear and two cloudy (with differing optical depth),
  ! returning a top-of-atmosphere spectral upwelling radiances. This
  ! version takes as input also the adjoint of the radiance and
  ! returns the adjoint of the active inpus (note that derivatives
  ! with respect to fractional_std are currently not computed).
  subroutine tcrad_radiance_lw_ad(ng, nlev, config, mu, &
       &     surf_emission, surf_albedo, &
       &     planck_hl, cloud_fraction, fractional_std, &
       &     od_clear, od_cloud, ssa_cloud, asymmetry_cloud, &
       &     overlap_param, &
       &     radiance, &
       &     radiance_ad, &
       &     surf_emission_ad, surf_albedo_ad, &
       &     planck_hl_ad, cloud_fraction_ad, &
       &     od_clear_ad, od_cloud_ad, ssa_cloud_ad, asymmetry_cloud_ad, &
       &     overlap_param_ad) bind(c)

    use iso_c_binding

    import :: tcrad_config
    
    integer(c_int), value        :: ng, nlev
    type(tcrad_config),intent(in):: config
    real(c_double), value        :: mu
    real(c_double), intent(in)   :: surf_emission(:), surf_albedo(:)
    real(c_double), intent(in)   :: planck_hl(:,:)
    real(c_double), intent(in)   :: cloud_fraction(:), fractional_std(:)
    real(c_double), intent(in)   :: od_clear(:,:), od_cloud(:,:)
    real(c_double), intent(in)   :: ssa_cloud(:,:), asymmetry_cloud(:,:)
    real(c_double), intent(in)   :: overlap_param(:)
    real(c_double), intent(out)  :: radiance(:)
    real(c_double), intent(in)   :: radiance_ad(:)
    real(c_double), intent(out)  :: surf_emission_ad(:), surf_albedo_ad(:)
    real(c_double), intent(out)  :: planck_hl_ad(:,:)
    real(c_double), intent(out)  :: cloud_fraction_ad(:)
    real(c_double), intent(out)  :: od_clear_ad(:,:), od_cloud_ad(:,:)
    real(c_double), intent(out)  :: ssa_cloud_ad(:,:), asymmetry_cloud_ad(:,:)
    real(c_double), intent(out)  :: overlap_param_ad(:)
    
  end subroutine tcrad_radiance_lw_ad
   
end interface
