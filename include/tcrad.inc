! -*- f90 -*-


enum, bind(c)
  enumerator :: TWO_STREAM_ELSASSER = 0
  enumerator :: TWO_STREAM_EDDINGTON = 1
  enumerator :: TWO_STREAM_LEGENDRE = 2
  enumerator :: TWO_STREAM_HYBRID = 3
  enumerator :: TWO_STREAM_SCALED_WISCOMBE_GRAMS = 4
end enum

type, bind(c) :: tcrad_config
  logical(c_bool)  :: do_specular_surface = .false.
  integer(c_int)   :: two_stream_scheme = TWO_STREAM_ELSASSER
  real(c_double) :: cloud_fraction_threshold = 1.0e-6
end type tcrad_config
  
interface

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

enum, bind(c)
  enumerator :: TWO_STREAM_ELSASSER = 0
  enumerator :: TWO_STREAM_EDDINGTON = 1
  enumerator :: TWO_STREAM_LEGENDRE = 2
  enumerator :: TWO_STREAM_HYBRID = 3
  enumerator :: TWO_STREAM_SCALED_WISCOMBE_GRAMS = 4
end enum

  type, bind(c) :: tcrad_config
    logical(c_bool)  :: do_specular_surface = .false.
    integer(c_int)   :: two_stream_scheme = TWO_STREAM_ELSASSER
    real(c_double) :: cloud_fraction_threshold = 1.0e-6
  end type tcrad_config
  
  
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

 enum, bind(c)
   enumerator :: TWO_STREAM_ELSASSER = 0
   enumerator :: TWO_STREAM_EDDINGTON = 1
   enumerator :: TWO_STREAM_LEGENDRE = 2
   enumerator :: TWO_STREAM_HYBRID = 3
   enumerator :: TWO_STREAM_SCALED_WISCOMBE_GRAMS = 4
 end enum

 type, bind(c) :: tcrad_config
     logical(c_bool)  :: do_specular_surface = .false.
     integer(c_int)   :: two_stream_scheme = TWO_STREAM_ELSASSER
     real(c_double) :: cloud_fraction_threshold = 1.0e-6
   end type tcrad_config
  
  integer(c_int),   value      :: ng, nlev
  type(tcrad_config),intent(in):: config
  real(c_double), value      :: mu
  real(c_double), intent(in) :: surf_emission(:), surf_albedo(:)
  real(c_double), intent(in) :: planck_hl(:,:)
  real(c_double), intent(in) :: cloud_fraction(:), fractional_std(:)
  real(c_double), intent(in) :: od_clear(:,:), od_cloud(:,:)
  real(c_double), intent(in) :: ssa_cloud(:,:), asymmetry_cloud(:,:)
  real(c_double), intent(in) :: overlap_param(:)
  real(c_double), intent(out) :: radiance(:)
  
  end subroutine tcrad_radiance_lw

   
end interface
