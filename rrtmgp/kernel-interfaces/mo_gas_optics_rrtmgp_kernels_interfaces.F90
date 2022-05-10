! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2022-,  Atmospheric and Environmental Research,
!    Regents of the University of Colorado,
!    Trustees of Columbia University in the City of New York
! All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description: Numeric calculations for gas optics.
!  Absorption and Rayleigh optical depths,
!   source functions.

module mo_gas_optics_rrtmgp_kernels_interface
  use mo_rte_kind, only : wp, wl
  implicit none
  private
  public :: interpolation, compute_tau_absorption, compute_tau_rayleigh, compute_Planck_source

  interface
    ! Compute interpolation coefficients
    ! for calculations of major optical depths, minor optical depths, Rayleigh,
    ! and Planck fractions
    module subroutine interpolation(                       &
                  ncol,nlay,ngas,nflav,neta, npres, ntemp, &
                  flavor,                                  &
                  press_ref_log, temp_ref,press_ref_log_delta,    &
                  temp_ref_min,temp_ref_delta,press_ref_trop_log, &
                  vmr_ref,                                        &
                  play,tlay,col_gas,                              &
                  jtemp,fmajor,fminor,col_mix,tropo,jeta,jpress) bind(C, name="rrtmgp_interpolation")
      ! input dimensions
      use mo_rte_kind, only : wp, wl
      integer,                            intent(in) :: ncol,nlay
      integer,                            intent(in) :: ngas,nflav,neta,npres,ntemp
      integer,     dimension(2,nflav),    intent(in) :: flavor
      real(wp),    dimension(npres),      intent(in) :: press_ref_log
      real(wp),    dimension(ntemp),      intent(in) :: temp_ref
      real(wp),                           intent(in) :: press_ref_log_delta, &
                                                        temp_ref_min, temp_ref_delta, &
                                                        press_ref_trop_log
      real(wp),    dimension(2,0:ngas,ntemp), intent(in) :: vmr_ref

      ! inputs from profile or parent function
      real(wp),    dimension(ncol,nlay),        intent(in) :: play, tlay
      real(wp),    dimension(ncol,nlay,0:ngas), intent(in) :: col_gas

      ! outputs
      integer,     dimension(ncol,nlay), intent(out) :: jtemp, jpress
      logical(wl), dimension(ncol,nlay), intent(out) :: tropo
      integer,     dimension(2,    ncol,nlay,nflav), intent(out) :: jeta
      real(wp),    dimension(2,    ncol,nlay,nflav), intent(out) :: col_mix
      real(wp),    dimension(2,2,2,ncol,nlay,nflav), intent(out) :: fmajor
      real(wp),    dimension(2,2,  ncol,nlay,nflav), intent(out) :: fminor
    end subroutine interpolation

    !
    ! Compute minor and major species opitcal depth from pre-computed interpolation coefficients
    !   (jeta,jtemp,jpress)
    !
    module subroutine compute_tau_absorption(         &
                  ncol,nlay,nbnd,ngpt,                &  ! dimensions
                  ngas,nflav,neta,npres,ntemp,        &
                  nminorlower, nminorklower,          & ! number of minor contributors, total num absorption coeffs
                  nminorupper, nminorkupper,          &
                  idx_h2o,                            &
                  gpoint_flavor,                      &
                  band_lims_gpt,                      &
                  kmajor,                             &
                  kminor_lower,                       &
                  kminor_upper,                       &
                  minor_limits_gpt_lower,             &
                  minor_limits_gpt_upper,             &
                  minor_scales_with_density_lower,    &
                  minor_scales_with_density_upper,    &
                  scale_by_complement_lower,          &
                  scale_by_complement_upper,          &
                  idx_minor_lower,                    &
                  idx_minor_upper,                    &
                  idx_minor_scaling_lower,            &
                  idx_minor_scaling_upper,            &
                  kminor_start_lower,                 &
                  kminor_start_upper,                 &
                  tropo,                              &
                  col_mix,fmajor,fminor,              &
                  play,tlay,col_gas,                  &
                  jeta,jtemp,jpress,                  &
                  tau) bind(C, name="rrtmgp_compute_tau_absorption")
      ! ---------------------
      ! input dimensions
      use mo_rte_kind, only : wp, wl
      integer,                                intent(in) :: ncol,nlay,nbnd,ngpt
      integer,                                intent(in) :: ngas,nflav,neta,npres,ntemp
      integer,                                intent(in) :: nminorlower, nminorklower,nminorupper, nminorkupper
      integer,                                intent(in) :: idx_h2o
      ! ---------------------
      ! inputs from object
      integer,     dimension(2,ngpt),                  intent(in) :: gpoint_flavor
      integer,     dimension(2,nbnd),                  intent(in) :: band_lims_gpt
      real(wp),    dimension(ntemp,neta,npres+1,ngpt), intent(in) :: kmajor
      real(wp),    dimension(ntemp,neta,nminorklower), intent(in) :: kminor_lower
      real(wp),    dimension(ntemp,neta,nminorkupper), intent(in) :: kminor_upper
      integer,     dimension(2,nminorlower),           intent(in) :: minor_limits_gpt_lower
      integer,     dimension(2,nminorupper),           intent(in) :: minor_limits_gpt_upper
      logical(wl), dimension(  nminorlower),           intent(in) :: minor_scales_with_density_lower
      logical(wl), dimension(  nminorupper),           intent(in) :: minor_scales_with_density_upper
      logical(wl), dimension(  nminorlower),           intent(in) :: scale_by_complement_lower
      logical(wl), dimension(  nminorupper),           intent(in) :: scale_by_complement_upper
      integer,     dimension(  nminorlower),           intent(in) :: idx_minor_lower
      integer,     dimension(  nminorupper),           intent(in) :: idx_minor_upper
      integer,     dimension(  nminorlower),           intent(in) :: idx_minor_scaling_lower
      integer,     dimension(  nminorupper),           intent(in) :: idx_minor_scaling_upper
      integer,     dimension(  nminorlower),           intent(in) :: kminor_start_lower
      integer,     dimension(  nminorupper),           intent(in) :: kminor_start_upper
      logical(wl), dimension(ncol,nlay),               intent(in) :: tropo
      ! ---------------------
      ! inputs from profile or parent function
      real(wp), dimension(2,    ncol,nlay,nflav       ), intent(in) :: col_mix
      real(wp), dimension(2,2,2,ncol,nlay,nflav       ), intent(in) :: fmajor
      real(wp), dimension(2,2,  ncol,nlay,nflav       ), intent(in) :: fminor
      real(wp), dimension(            ncol,nlay       ), intent(in) :: play, tlay      ! pressure and temperature
      real(wp), dimension(            ncol,nlay,0:ngas), intent(in) :: col_gas
      integer,  dimension(2,    ncol,nlay,nflav       ), intent(in) :: jeta
      integer,  dimension(            ncol,nlay       ), intent(in) :: jtemp
      integer,  dimension(            ncol,nlay       ), intent(in) :: jpress
      ! ---------------------
      ! output - optical depth
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau
    end subroutine compute_tau_absorption

    ! ----------------------------------------------------------
    !
    ! compute Rayleigh scattering optical depths
    !
    module subroutine compute_tau_rayleigh(ncol,nlay,nbnd,ngpt,  &
                                    ngas,nflav,neta,npres,ntemp, &
                                    gpoint_flavor,band_lims_gpt, &
                                    krayl,                       &
                                    idx_h2o, col_dry,col_gas,    &
                                    fminor,jeta,tropo,jtemp,     &
                                    tau_rayleigh) bind(C, name="rrtmgp_compute_tau_rayleigh")
      use mo_rte_kind, only : wp, wl
      integer,                                     intent(in ) :: ncol,nlay,nbnd,ngpt
      integer,                                     intent(in ) :: ngas,nflav,neta,npres,ntemp
      integer,     dimension(2,ngpt),              intent(in ) :: gpoint_flavor
      integer,     dimension(2,nbnd),              intent(in ) :: band_lims_gpt ! start and end g-point for each band
      real(wp),    dimension(ntemp,neta,ngpt,2),   intent(in ) :: krayl
      integer,                                     intent(in ) :: idx_h2o
      real(wp),    dimension(ncol,nlay),           intent(in ) :: col_dry
      real(wp),    dimension(ncol,nlay,0:ngas),    intent(in ) :: col_gas
      real(wp),    dimension(2,2,ncol,nlay,nflav), intent(in ) :: fminor
      integer,     dimension(2,  ncol,nlay,nflav), intent(in ) :: jeta
      logical(wl), dimension(ncol,nlay),           intent(in ) :: tropo
      integer,     dimension(ncol,nlay),           intent(in ) :: jtemp
      ! outputs
      real(wp),    dimension(ncol,nlay,ngpt),      intent(out) :: tau_rayleigh
    end subroutine compute_tau_rayleigh

    ! ----------------------------------------------------------
    module subroutine compute_Planck_source(                 &
                      ncol, nlay, nbnd, ngpt,                &
                      nflav, neta, npres, ntemp, nPlanckTemp,&
                      tlay, tlev, tsfc, sfc_lay,             &
                      fmajor, jeta, tropo, jtemp, jpress,    &
                      gpoint_bands, band_lims_gpt,           &
                      pfracin, temp_ref_min, totplnk_delta, totplnk, gpoint_flavor, &
                      sfc_src, lay_src, lev_src_inc, lev_src_dec, sfc_source_Jac) bind(C, name="rrtmgp_compute_Planck_source")
      use mo_rte_kind, only : wp, wl
      integer,                                    intent(in) :: ncol, nlay, nbnd, ngpt
      integer,                                    intent(in) :: nflav, neta, npres, ntemp, nPlanckTemp
      real(wp),    dimension(ncol,nlay  ),        intent(in) :: tlay
      real(wp),    dimension(ncol,nlay+1),        intent(in) :: tlev
      real(wp),    dimension(ncol       ),        intent(in) :: tsfc
      integer,                                    intent(in) :: sfc_lay
      ! Interpolation variables
      real(wp),    dimension(2,2,2,ncol,nlay,nflav), intent(in) :: fmajor
      integer,     dimension(2,    ncol,nlay,nflav), intent(in) :: jeta
      logical(wl), dimension(            ncol,nlay), intent(in) :: tropo
      integer,     dimension(            ncol,nlay), intent(in) :: jtemp, jpress
      ! Table-specific
      integer, dimension(ngpt),                     intent(in) :: gpoint_bands ! start and end g-point for each band
      integer, dimension(2, nbnd),                  intent(in) :: band_lims_gpt ! start and end g-point for each band
      real(wp),                                     intent(in) :: temp_ref_min, totplnk_delta
      real(wp), dimension(ntemp,neta,npres+1,ngpt), intent(in) :: pfracin
      real(wp), dimension(nPlanckTemp,nbnd),        intent(in) :: totplnk
      integer,  dimension(2,ngpt),                  intent(in) :: gpoint_flavor

      real(wp), dimension(ncol,     ngpt), intent(out) :: sfc_src
      real(wp), dimension(ncol,nlay,ngpt), intent(out) :: lay_src
      real(wp), dimension(ncol,nlay,ngpt), intent(out) :: lev_src_inc, lev_src_dec
      real(wp), dimension(ncol,     ngpt), intent(out) :: sfc_source_Jac
    end subroutine compute_Planck_source
  end interface
end module mo_gas_optics_rrtmgp_kernels_interface
