module micro_pumas_v1
!---------------------------------------------------------------------------------
! Parameterization of Unified Microphysics Across Scales version 1 (PUMASv1)
!
! References:
!
!           Gettelman, A., H. Morrison, T. Eidhammer, K. Thayer-Calder, J. Sun,
!
!           R. Forbes, Z. McGraw, J. Zhu, T. Storelvmo, and J. Dennis (2023):
!
!           Importance of Ice Nucleation and Precipitation on Climate with the
!
!           Parameterization of Unified Microphysics Across Scales version 1
!
!           (PUMASv1). Geosci. Model Dev., 16, 1735-1754.
!
!           https://doi.org/10.5194/gmd-16-1735-2023
!
!
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!---------------------------------------------------------------------------------
!
! NOTE: Modified to allow other microphysics packages (e.g. CARMA) to do ice
! microphysics in cooperation with the MG liquid microphysics. This is
! controlled by the do_cldice variable.
!
! If do_cldice is false, then MG microphysics should not update CLDICE or
! NUMICE; it is assumed that the other microphysics scheme will have updated
! CLDICE and NUMICE. The other microphysics should handle the following
! processes that would have been done by MG:
!   - Detrainment (liquid and ice)
!   - Homogeneous ice nucleation
!   - Heterogeneous ice nucleation
!   - Bergeron process
!   - Melting of ice
!   - Freezing of cloud drops
!   - Autoconversion (ice -> snow)
!   - Growth/Sublimation of ice
!   - Sedimentation of ice
!
! This option has not been updated since the introduction of prognostic
! precipitation, and probably should be adjusted to cover snow as well.
!
!---------------------------------------------------------------------------------
! Version 3.O based on micro_mg2_0.F90 and WRF3.8.1 module_mp_morr_two_moment.F
!---------------------------------------------------------------------------------
! Based on micro_mg (restructuring of former cldwat2m_micro)
! Author: Andrew Gettelman, Hugh Morrison.
! Contributions from: Xiaohong Liu and Steve Ghan
! December 2005-May 2010
! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!---------------------------------------------------------------------------------
! Code comments added by HM, 093011
! General code structure:
!
! Code is divided into two main subroutines:
!   subroutine micro_pumas_init --> initializes microphysics routine, should be called
!                                  once at start of simulation
!   subroutine micro_pumas_tend --> main microphysics routine to be called each time step
!                                this also calls several smaller subroutines to calculate
!                                microphysical processes and other utilities
!
! List of external functions:
!   qsat_water --> for calculating saturation vapor pressure with respect to liquid water
!   qsat_ice --> for calculating saturation vapor pressure with respect to ice
!   gamma   --> standard mathematical gamma function
! .........................................................................
! List of inputs through use statement in fortran90:
! Variable Name                      Description                Units
! .........................................................................
! gravit          acceleration due to gravity                    m s-2
! rair            dry air gas constant for air                  J kg-1 K-1
! tmelt           temperature of melting point for water          K
! cpair           specific heat at constant pressure for dry air J kg-1 K-1
! rh2o            gas constant for water vapor                  J kg-1 K-1
! latvap          latent heat of vaporization                   J kg-1
! latice          latent heat of fusion                         J kg-1
! qsat_water      external function for calculating liquid water
!                 saturation vapor pressure/humidity              -
! qsat_ice        external function for calculating ice
!                 saturation vapor pressure/humidity              pa
! rhmini          relative humidity threshold parameter for
!                 nucleating ice                                  -
! .........................................................................
! NOTE: List of all inputs/outputs passed through the call/subroutine statement
!       for micro_pumas_tend is given below at the start of subroutine micro_pumas_tend.
!---------------------------------------------------------------------------------

! Procedures required:
! 1) An implementation of the gamma function (if not intrinsic).
! 2) saturation vapor pressure and specific humidity over water
! 3) svp over ice

use pumas_kinds, only: r8=>kind_r8

#ifndef HAVE_GAMMA_INTRINSICS
use pumas_gamma_function, only: gamma=>pumas_gamma
#endif

use wv_sat_methods, only: &
     qsat_water => wv_sat_qsat_water_vect, &
     qsat_ice => wv_sat_qsat_ice_vect

! Parameters from the utilities module.
use micro_pumas_utils, only: &
     pi, &
     omsm, &
     qsmall, &
     mincld, &
     rhosn, &
     rhoi, &
     rhow, &
     rhows, &
     ac, bc, &
     ai, bi, &
     aj, bj, &
     ar, br, &
     as, bs, &
     ag, bg, &
     ah, bh, &
     rhog,rhoh, &
     mi0, &
     rising_factorial, &
     VLENS

implicit none
private
save

public :: &
     micro_pumas_init, &
     micro_pumas_get_cols, &
     micro_pumas_tend

! Switches for specification rather than prediction of droplet and crystal number
! note: number will be adjusted as needed to keep mean size within bounds,
! even when specified droplet or ice number is used
!
! If constant cloud ice number is set (nicons = .true.),
! then all microphysical processes except mass transfer due to ice nucleation
! (mnuccd) are based on the fixed cloud ice number. Calculation of
! mnuccd follows from the prognosed ice crystal number ni.

logical :: nccons ! nccons = .true. to specify constant cloud droplet number
logical :: nicons ! nicons = .true. to specify constant cloud ice number
logical :: ngcons ! ngcons = .true. to specify constant graupel number
logical :: nrcons ! constant rain number
logical :: nscons ! constant snow number

! specified ice and droplet number concentrations
! note: these are local in-cloud values, not grid-mean
real(r8) :: ncnst  ! droplet num concentration when nccons=.true. (m-3)
real(r8) :: ninst  ! ice num concentration when nicons=.true. (m-3)
real(r8) :: ngnst   ! graupel num concentration when ngcons=.true. (m-3)
real(r8) :: nrnst
real(r8) :: nsnst

! IFS Switches....
! Switch to turn off evaporation of sedimenting condensate
! Found to interact badly in some models with diagnostic cloud fraction
logical :: evap_sed_off

! Remove RH conditional from ice nucleation
logical :: icenuc_rh_off

! Internally: Meyers Ice Nucleation
logical :: icenuc_use_meyers

! Scale evaporation as IFS does (*0.3)
logical :: evap_scl_ifs

! Evap RH threshold following ifs
logical :: evap_rhthrsh_ifs

! Rain freezing at 0C following ifs

logical :: rainfreeze_ifs

! Snow sedimentation = 1 m/s

logical :: ifs_sed

! Precipitation fall speed, prevent zero velocity if precip above

logical :: precip_fall_corr

!--ag

!=========================================================
! Private module parameters
!=========================================================

!Range of cloudsat reflectivities (dBz) for analytic simulator
real(r8), parameter :: csmin = -30._r8
real(r8), parameter :: csmax = 26._r8
real(r8), parameter :: mindbz = -99._r8
real(r8), parameter :: minrefl = 1.26e-10_r8    ! minrefl = 10._r8**(mindbz/10._r8)

integer, parameter  :: MG_PRECIP_FRAC_INCLOUD = 101
integer, parameter  :: MG_PRECIP_FRAC_OVERLAP = 102

! Reflectivity min for 10cm (Rain) radar reflectivity
real(r8), parameter :: minrefl10 = 1.e-26_r8

! autoconversion size threshold for cloud ice to snow (m)
real(r8) :: dcs

! minimum mass of new crystal due to freezing of cloud droplets done
! externally (kg)
real(r8), parameter :: mi0l_min = 4._r8/3._r8*pi*rhow*(4.e-6_r8)**3

! Ice number sublimation parameter. Assume some decrease in ice number with sublimation if non-zero. Else, no decrease in number with sublimation.
real(r8), parameter :: sublim_factor =0.0_r8      !number sublimation factor.

! Parameters related to GPU computing
integer, parameter :: RQUEUE = 101    ! GPU stream ID for rain
integer, parameter :: SQUEUE = 102    ! GPU stream ID for snow
integer, parameter :: LQUEUE = 103    ! GPU stream ID for liquid
integer, parameter :: IQUEUE = 104    ! GPU stream ID for ice
integer, parameter :: GQUEUE = 105    ! GPU stream ID for hail/graupel

!=========================================================
! Constants set in initialization
!=========================================================

! Set using arguments to micro_pumas_init
real(r8) :: g           ! gravity
real(r8) :: r           ! dry air gas constant
real(r8) :: rv          ! water vapor gas constant
real(r8) :: cpp         ! specific heat of dry air
real(r8) :: tmelt       ! freezing point of water (K)

! latent heats of:
real(r8) :: xxlv        ! vaporization
real(r8) :: xlf         ! freezing
real(r8) :: xxls        ! sublimation

real(r8) :: rhmini      ! Minimum rh for ice cloud fraction > 0.

! flags
logical :: microp_uniform
logical :: do_cldice
logical :: use_hetfrz_classnuc
logical :: do_hail
logical :: do_graupel

real(r8) :: rhosu       ! typical 850mn air density

real(r8) :: icenuct     ! ice nucleation temperature: currently -5 degrees C

real(r8) :: snowmelt    ! what temp to melt all snow: currently 2 degrees C
real(r8) :: rainfrze    ! what temp to freeze all rain: currently -5 degrees C

! additional constants to help speed up code
real(r8) :: gamma_br_plus1
real(r8) :: gamma_br_plus4
real(r8) :: gamma_bs_plus1
real(r8) :: gamma_bs_plus4
real(r8) :: gamma_bi_plus1
real(r8) :: gamma_bi_plus4
real(r8) :: gamma_bj_plus1
real(r8) :: gamma_bj_plus4
real(r8) :: gamma_bg_plus1
real(r8) :: gamma_bg_plus4
real(r8) :: xxlv_squared
real(r8) :: xxls_squared

character(len=16)  :: micro_mg_precip_frac_method  ! type of precipitation fraction method
real(r8)           :: micro_mg_berg_eff_factor     ! berg efficiency factor

real(r8)           :: micro_mg_accre_enhan_fact     ! accretion enhancment factor
real(r8)           :: micro_mg_autocon_fact     ! autoconversion prefactor
real(r8)           :: micro_mg_autocon_nd_exp     ! autoconversion Nd exponent factor
real(r8)           :: micro_mg_autocon_lwp_exp  !autoconversion LWP exponent
real(r8)           :: micro_mg_homog_size ! size of freezing homogeneous ice
real(r8)           :: micro_mg_vtrmi_factor
real(r8)           :: micro_mg_vtrms_factor
real(r8)           :: micro_mg_effi_factor
real(r8)           :: micro_mg_iaccr_factor
real(r8)           :: micro_mg_max_nicons

logical            :: remove_supersat      ! If true, remove supersaturation after sedimentation loop
character(len=16)  :: warm_rain            ! 'tau','emulated','sb2001' or 'kk2000'

!Parameters for Implicit Sedimentation Calculation
real(r8), parameter :: vfactor = 1.0_r8        ! Rain/Snow/Graupel Factor
real(r8), parameter :: vfac_drop = 1.0_r8      ! Cloud Liquid Factor
real(r8), parameter :: vfac_ice  = 1.0_r8      ! Cloud Ice Factor

logical           :: do_implicit_fall !   = .true.

logical           :: accre_sees_auto  != .true.

!$acc declare create (nccons,nicons,ngcons,nrcons,nscons,ncnst,ninst,ngnst,    &
!$acc                 nrnst,nsnst,evap_sed_off,icenuc_rh_off,evap_scl_ifs,     &
!$acc                 icenuc_use_meyers,evap_rhthrsh_ifs,rainfreeze_ifs,       &
!$acc                 ifs_sed,precip_fall_corr,dcs,                            &
!$acc                 g,r,rv,cpp,tmelt,xxlv,xlf,xxls,rhmini,microp_uniform,    &
!$acc                 do_cldice,use_hetfrz_classnuc,do_hail,do_graupel,rhosu,  &
!$acc                 icenuct,snowmelt,rainfrze,xxlv_squared,xxls_squared,     &
!$acc                 gamma_br_plus1,gamma_br_plus4,gamma_bs_plus1,            &
!$acc                 gamma_bs_plus4,gamma_bi_plus1,gamma_bi_plus4,            &
!$acc                 gamma_bj_plus1,gamma_bj_plus4,gamma_bg_plus1,            &
!$acc                 gamma_bg_plus4,micro_mg_berg_eff_factor,                 &
!$acc                 micro_mg_accre_enhan_fact,micro_mg_autocon_fact,         &
!$acc                 micro_mg_autocon_nd_exp,micro_mg_autocon_lwp_exp,        &
!$acc                 micro_mg_homog_size,micro_mg_vtrmi_factor,               &
!$acc                 micro_mg_vtrms_factor,                                   &
!$acc                 micro_mg_effi_factor,micro_mg_iaccr_factor,              &
!$acc                 micro_mg_max_nicons,remove_supersat,do_implicit_fall,    &
!$acc                 accre_sees_auto)

!===============================================================================
contains
!===============================================================================

subroutine micro_pumas_init( &
     kind, gravit, rair, rh2o, cpair,    &
     tmelt_in, latvap, latice,           &
     rhmini_in, micro_mg_dcs,            &
     micro_mg_do_hail_in,micro_mg_do_graupel_in, &
     microp_uniform_in, do_cldice_in, use_hetfrz_classnuc_in, &
     micro_mg_precip_frac_method_in, micro_mg_berg_eff_factor_in, &
     micro_mg_accre_enhan_fact_in, micro_mg_autocon_fact_in, &
     micro_mg_autocon_nd_exp_in, micro_mg_autocon_lwp_exp_in, micro_mg_homog_size_in, &
     micro_mg_vtrmi_factor_in, micro_mg_vtrms_factor_in, micro_mg_effi_factor_in, &
     micro_mg_iaccr_factor_in, micro_mg_max_nicons_in, &
     remove_supersat_in, warm_rain_in, &
     micro_mg_evap_sed_off_in, micro_mg_icenuc_rh_off_in, micro_mg_icenuc_use_meyers_in, &
     micro_mg_evap_scl_ifs_in, micro_mg_evap_rhthrsh_ifs_in, &
     micro_mg_rainfreeze_ifs_in,  micro_mg_ifs_sed_in, micro_mg_precip_fall_corr, &
     micro_mg_accre_sees_auto_in, micro_mg_implicit_fall_in, &
     nccons_in, nicons_in, ncnst_in, ninst_in, ngcons_in, ngnst_in, &
     nrcons_in, nrnst_in, nscons_in, nsnst_in, &
     stochastic_emulated_filename_quantile, stochastic_emulated_filename_input_scale, &
     stochastic_emulated_filename_output_scale, &
     iulog, errstring)

  use micro_pumas_utils,            only: micro_pumas_utils_init
  use pumas_stochastic_collect_tau, only: pumas_stochastic_kernel_init
  use tau_neural_net_quantile,      only: initialize_tau_emulators

  !-----------------------------------------------------------------------
  !
  ! Purpose:
  ! initialize constants for MG microphysics
  !
  ! Author: Andrew Gettelman Dec 2005
  !
  !-----------------------------------------------------------------------

  integer,  intent(in)  :: kind         ! Kind used for reals
  real(r8), intent(in)  :: gravit
  real(r8), intent(in)  :: rair
  real(r8), intent(in)  :: rh2o
  real(r8), intent(in)  :: cpair
  real(r8), intent(in)  :: tmelt_in     ! Freezing point of water (K)
  real(r8), intent(in)  :: latvap
  real(r8), intent(in)  :: latice
  real(r8), intent(in)  :: rhmini_in    ! Minimum rh for ice cloud fraction > 0.
  real(r8), intent(in)  :: micro_mg_dcs

!MG3 dense precipitating ice. Note, only 1 can be true, or both false.
  logical,  intent(in)  :: micro_mg_do_graupel_in    ! .true. = configure with graupel
                                                   ! .false. = no graupel (hail possible)
  logical,  intent(in)  :: micro_mg_do_hail_in    ! .true. = configure with hail
                                                   ! .false. = no hail (graupel possible)
  logical,  intent(in)  :: microp_uniform_in    ! .true. = configure uniform for sub-columns
                                            ! .false. = use w/o sub-columns (standard)
  logical,  intent(in)  :: do_cldice_in     ! .true. = do all processes (standard)
                                            ! .false. = skip all processes affecting
                                            !           cloud ice
  logical,  intent(in)  :: use_hetfrz_classnuc_in ! use heterogeneous freezing

  character(len=16),intent(in)  :: micro_mg_precip_frac_method_in  ! type of precipitation fraction method
  real(r8),         intent(in)  :: micro_mg_berg_eff_factor_in     ! berg efficiency factor
  real(r8),         intent(in)  :: micro_mg_accre_enhan_fact_in     !accretion enhancment factor
  real(r8),         intent(in) ::  micro_mg_autocon_fact_in    !autconversion prefactor
  real(r8),         intent(in) ::  micro_mg_autocon_nd_exp_in !autconversion exponent factor
  real(r8),         intent(in) ::  micro_mg_autocon_lwp_exp_in    !autconversion exponent factor
  real(r8),         intent(in) ::  micro_mg_homog_size_in  ! size of homoegenous freezing ice
  real(r8),         intent(in)  :: micro_mg_vtrmi_factor_in    !factor for ice fall velocity
  real(r8),         intent(in)  :: micro_mg_vtrms_factor_in    !factor for snow fall velocity
  real(r8),         intent(in)  :: micro_mg_effi_factor_in    !factor for ice effective radius
  real(r8),         intent(in)  :: micro_mg_iaccr_factor_in  ! ice accretion factor
  real(r8),         intent(in)  :: micro_mg_max_nicons_in ! maximum number ice crystal allowed

  logical,  intent(in)  ::  remove_supersat_in ! If true, remove supersaturation after sedimentation loop
  character(len=*),  intent(in)  ::  warm_rain_in

! IFS-like Switches

  logical, intent(in) :: micro_mg_evap_sed_off_in ! Turn off evaporation/sublimation based on cloud fraction for sedimenting condensate

  logical, intent(in) :: micro_mg_icenuc_rh_off_in ! Remove RH conditional from ice nucleation
  logical, intent(in) :: micro_mg_icenuc_use_meyers_in ! Internally: Meyers Ice Nucleation
  logical, intent(in) :: micro_mg_evap_scl_ifs_in ! Scale evaporation as IFS does (*0.3)
  logical, intent(in) :: micro_mg_evap_rhthrsh_ifs_in ! Evap RH threshold following ifs
  logical, intent(in) :: micro_mg_rainfreeze_ifs_in ! Rain freezing temp following ifs
  logical, intent(in) :: micro_mg_ifs_sed_in ! snow sedimentation = 1m/s following ifs
  logical, intent(in) :: micro_mg_precip_fall_corr ! ensure rain fall speed non-zero if rain above in column

  logical, intent(in) :: micro_mg_accre_sees_auto_in ! autoconverted rain is passed to accretion

  logical, intent(in) :: micro_mg_implicit_fall_in !Implicit fall speed (sedimentation) calculation for hydrometors



  logical, intent(in)   :: nccons_in
  logical, intent(in)   :: nicons_in
  real(r8), intent(in)  :: ncnst_in
  real(r8), intent(in)  :: ninst_in

  logical, intent(in)   :: ngcons_in
  real(r8), intent(in)  :: ngnst_in
  logical, intent(in)   :: nrcons_in
  real(r8), intent(in)  :: nrnst_in
  logical, intent(in)   :: nscons_in
  real(r8), intent(in)  :: nsnst_in

  character(len=*), intent(in) :: stochastic_emulated_filename_quantile, &
                                  stochastic_emulated_filename_input_scale, &
                                  stochastic_emulated_filename_output_scale ! Files for emulated machine learning

  integer, intent(in) :: iulog                ! Host model log file unit number
  character(128), intent(out) :: errstring    ! Output status (non-blank for error return)

  !-----------------------------------------------------------------------

  ! Initialize ice autoconversion size
  ! to host-model specified value:
  dcs = micro_mg_dcs

  ! Initialize subordinate utilities module.
  call micro_pumas_utils_init(kind, rair, rh2o, cpair, tmelt_in, latvap, latice, &
       dcs, errstring)

  if (trim(errstring) /= "") then !If True then error occured
     return
  end if

  ! declarations for MG code (transforms variable names)

  g= gravit                 ! gravity
  r= rair                   ! dry air gas constant: note units(phys_constants are in J/K/kmol)
  rv= rh2o                  ! water vapor gas constant
  cpp = cpair               ! specific heat of dry air
  tmelt = tmelt_in
  rhmini = rhmini_in
  micro_mg_precip_frac_method = micro_mg_precip_frac_method_in
  micro_mg_berg_eff_factor    = micro_mg_berg_eff_factor_in
  micro_mg_accre_enhan_fact   =  micro_mg_accre_enhan_fact_in
  micro_mg_autocon_fact  = micro_mg_autocon_fact_in
  micro_mg_autocon_nd_exp = micro_mg_autocon_nd_exp_in
  micro_mg_autocon_lwp_exp = micro_mg_autocon_lwp_exp_in
  micro_mg_homog_size   = micro_mg_homog_size_in
  micro_mg_vtrmi_factor = micro_mg_vtrmi_factor_in
  micro_mg_vtrms_factor = micro_mg_vtrms_factor_in
  micro_mg_effi_factor = micro_mg_effi_factor_in
  micro_mg_iaccr_factor = micro_mg_iaccr_factor_in
  micro_mg_max_nicons = micro_mg_max_nicons_in
  remove_supersat          = remove_supersat_in
  warm_rain                = warm_rain_in
  do_implicit_fall   = micro_mg_implicit_fall_in
  accre_sees_auto = micro_mg_accre_sees_auto_in

  nccons = nccons_in
  nicons = nicons_in
  ncnst  = ncnst_in
  ninst  = ninst_in
  ngcons = ngcons_in
  ngnst  = ngnst_in
  nscons = nscons_in
  nsnst  = nsnst_in
  nrcons = nrcons_in
  nrnst  = nrnst_in

  ! latent heats

  xxlv = latvap         ! latent heat vaporization
  xlf  = latice         ! latent heat freezing
  xxls = xxlv + xlf     ! latent heat of sublimation

  ! flags
  microp_uniform = microp_uniform_in
  do_cldice  = do_cldice_in
  use_hetfrz_classnuc = use_hetfrz_classnuc_in
  do_hail = micro_mg_do_hail_in
  do_graupel = micro_mg_do_graupel_in
  evap_sed_off = micro_mg_evap_sed_off_in
  icenuc_rh_off = micro_mg_icenuc_rh_off_in
  icenuc_use_meyers = micro_mg_icenuc_use_meyers_in
  evap_scl_ifs = micro_mg_evap_scl_ifs_in
  evap_rhthrsh_ifs = micro_mg_evap_rhthrsh_ifs_in
  rainfreeze_ifs = micro_mg_rainfreeze_ifs_in
  ifs_sed = micro_mg_ifs_sed_in
  precip_fall_corr = micro_mg_precip_fall_corr
  ! typical air density at 850 mb

  ! Currently "do_hail" and "do_graupel" cannot both
  ! be True, so raise an error if that is the case:
  if (do_hail .and. do_graupel) then
     errstring = "The variables 'do_hail' and 'do_graupel' cannot both be true.  Please set one to '.false.'."
     return
  end if

  rhosu = 85000._r8/(rair * tmelt)

  ! Maximum temperature at which snow is allowed to exist
  snowmelt = tmelt + 2._r8
  ! Minimum temperature at which rain is allowed to exist
   if (rainfreeze_ifs) then
      rainfrze = tmelt
   else
      rainfrze = tmelt - 40._r8
   end if


  ! Ice nucleation temperature
  icenuct  = tmelt - 5._r8

  ! Define constants to help speed up code (this limits calls to gamma function)
  gamma_br_plus1=gamma(1._r8+br)
  gamma_br_plus4=gamma(4._r8+br)
  gamma_bs_plus1=gamma(1._r8+bs)
  gamma_bs_plus4=gamma(4._r8+bs)
  gamma_bi_plus1=gamma(1._r8+bi)
  gamma_bi_plus4=gamma(4._r8+bi)
  gamma_bj_plus1=gamma(1._r8+bj)
  gamma_bj_plus4=gamma(4._r8+bj)
  gamma_bg_plus1=gamma(1._r8)
  gamma_bg_plus4=gamma(4._r8)
  if (do_hail) then
     gamma_bg_plus1 = gamma(1._r8+bh)
     gamma_bg_plus4 = gamma(4._r8+bh)
  end if
  if (do_graupel) then
     gamma_bg_plus1 = gamma(1._r8+bg)
     gamma_bg_plus4 = gamma(4._r8+bg)
  end if

  xxlv_squared=xxlv**2
  xxls_squared=xxls**2

  !$acc update device (nccons,nicons,ngcons,nrcons,nscons,ncnst,ninst,ngnst,   &
  !$acc                nrnst,nsnst,evap_sed_off,icenuc_rh_off,evap_scl_ifs,    &
  !$acc                icenuc_use_meyers,evap_rhthrsh_ifs,rainfreeze_ifs,      &
  !$acc                ifs_sed,precip_fall_corr,dcs,                           &
  !$acc                g,r,rv,cpp,tmelt,xxlv,xlf,xxls,rhmini,microp_uniform,   &
  !$acc                do_cldice,use_hetfrz_classnuc,do_hail,do_graupel,rhosu, &
  !$acc                icenuct,snowmelt,rainfrze,xxlv_squared,xxls_squared,    &
  !$acc                gamma_br_plus1,gamma_br_plus4,gamma_bs_plus1,           &
  !$acc                gamma_bs_plus4,gamma_bi_plus1,gamma_bi_plus4,           &
  !$acc                gamma_bj_plus1,gamma_bj_plus4,gamma_bg_plus1,           &
  !$acc                gamma_bg_plus4,micro_mg_berg_eff_factor,                &
  !$acc                micro_mg_accre_enhan_fact,micro_mg_autocon_fact,        &
  !$acc                micro_mg_autocon_nd_exp,micro_mg_autocon_lwp_exp,       &
  !$acc                micro_mg_homog_size,micro_mg_vtrmi_factor,              &
  !$acc                micro_mg_vtrms_factor,                                  &
  !$acc                micro_mg_effi_factor,micro_mg_iaccr_factor,             &
  !$acc                micro_mg_max_nicons,remove_supersat,do_implicit_fall,   &
  !$acc                accre_sees_auto)

  if (trim(warm_rain) == 'emulated') then
      call initialize_tau_emulators(stochastic_emulated_filename_quantile, stochastic_emulated_filename_input_scale, &
                                    stochastic_emulated_filename_output_scale, iulog, errstring)
  end if

end subroutine micro_pumas_init

!===============================================================================
!microphysics routine for each timestep goes here...

subroutine micro_pumas_tend ( &
     mgncol,             nlev,               deltatin,           &
     t,                            q,                            &
     qcn,                          qin,                          &
     ncn,                          nin,                          &
     qrn,                          qsn,                          &
     nrn,                          nsn,                          &
     qgr,                          ngr,                          &
     relvar,                       accre_enhan,                  &
     p,                            pdel, pint,                   &
     cldn,    liqcldf,        icecldf,       qsatfac,            &
     qcsinksum_rate1ord,                                         &
     naai,                         npccn,                        &
     rndst,                        nacon,                        &
     tlat,                         qvlat,                        &
     qctend,                       qitend,                       &
     nctend,                       nitend,                       &
     qrtend,                       qstend,                       &
     nrtend,                       nstend,                       &
     qgtend,                       ngtend,                       &
     effc,               effc_fn,            effi,               &
     sadice,                       sadsnow,                      &
     prect,                        preci,                        &
     nevapr,                       am_evp_st,                    &
     prain,                                                      &
     cmeout,                       deffi,                        &
     pgamrad,                      lamcrad,                      &
     qsout,                        dsout,                        &
     qgout,     ngout,             dgout,                        &
     lflx,               iflx,                                   &
     gflx,                                                       &
     rflx,               sflx,               qrout,              &
     reff_rain,          reff_snow,          reff_grau,          &
     nrout,                        nsout,                        &
     refl,               arefl,              areflz,             &
     frefl,              csrfl,              acsrfl,             &
     fcsrfl,        refl10cm, reflz10cm,     rercld,             &
     ncai,                         ncal,                         &
     qrout2,                       qsout2,                       &
     nrout2,                       nsout2,                       &
     drout2,                       dsout2,                       &
     qgout2,        ngout2,        dgout2,    freqg,                   &
     freqs,                        freqr,                        &
     nfice,                        qcrat,                        &
     proc_rates,                                                 &
     errstring, &
     tnd_qsnow,          tnd_nsnow,          re_ice,             &
     prer_evap,                                                      &
     frzimm,             frzcnt,             frzdep)

  use pumas_stochastic_collect_tau, only: ncd, pumas_stochastic_collect_tau_tend
  use tau_neural_net_quantile,      only: tau_emulated_cloud_rain_interactions
  use ML_fixer_check,               only: ML_fixer_calc

  ! Constituent properties.
  use micro_pumas_utils, only: &
       mg_liq_props, &
       mg_ice_props, &
       mg_rain_props, &
       mg_graupel_props, &
       mg_hail_props, &
       mg_snow_props

  ! Size calculation functions.
  use micro_pumas_utils, only: &
       size_dist_param_liq, &
       size_dist_param_basic, &
       avg_diameter, &
       avg_diameter_vec

  ! Microphysical processes.
  use micro_pumas_utils, only: &
       ice_deposition_sublimation, &
       sb2001v2_liq_autoconversion,&
       sb2001v2_accre_cld_water_rain,&
       kk2000_liq_autoconversion, &
       ice_autoconversion, &
       immersion_freezing, &
       contact_freezing, &
       snow_self_aggregation, &
       accrete_cloud_water_snow, &
       secondary_ice_production, &
       accrete_rain_snow, &
       heterogeneous_rain_freezing, &
       accrete_cloud_water_rain, &
       self_collection_rain, &
       accrete_cloud_ice_snow, &
       evaporate_sublimate_precip, &
       bergeron_process_snow, &
       graupel_collecting_snow, &
       graupel_collecting_rain, &
       graupel_collecting_cld_water, &
       graupel_riming_liquid_snow, &
       graupel_rain_riming_snow, &
       graupel_rime_splintering, &
       vapor_deposition_onto_snow, &
       evaporate_sublimate_precip_graupel

  use micro_pumas_diags, only: proc_rates_type

  !Authors: Hugh Morrison, Andrew Gettelman, NCAR, Peter Caldwell, LLNL
  ! e-mail: morrison@ucar.edu, andrew@ucar.edu

  ! input arguments
  integer,  intent(in) :: mgncol         ! number of microphysics columns
  integer,  intent(in) :: nlev           ! number of layers
  real(r8), intent(in) :: deltatin       ! time step (s)
  real(r8), intent(in) :: t(mgncol,nlev) ! input temperature (K)
  real(r8), intent(in) :: q(mgncol,nlev) ! input h20 vapor mixing ratio (kg/kg)

  ! note: all input cloud variables are grid-averaged
  real(r8), intent(in) :: qcn(mgncol,nlev)       ! cloud water mixing ratio (kg/kg)
  real(r8), intent(in) :: qin(mgncol,nlev)       ! cloud ice mixing ratio (kg/kg)
  real(r8), intent(in) :: ncn(mgncol,nlev)       ! cloud water number conc (1/kg)
  real(r8), intent(in) :: nin(mgncol,nlev)       ! cloud ice number conc (1/kg)

  real(r8), intent(in) :: qrn(mgncol,nlev)       ! rain mixing ratio (kg/kg)
  real(r8), intent(in) :: qsn(mgncol,nlev)       ! snow mixing ratio (kg/kg)
  real(r8), intent(in) :: nrn(mgncol,nlev)       ! rain number conc (1/kg)
  real(r8), intent(in) :: nsn(mgncol,nlev)       ! snow number conc (1/kg)
  real(r8), intent(in) :: qgr(mgncol,nlev)       ! graupel/hail mixing ratio (kg/kg)
  real(r8), intent(in) :: ngr(mgncol,nlev)       ! graupel/hail number conc (1/kg)

  real(r8), intent(in) :: relvar(mgncol,nlev)      ! cloud water relative variance (-)
  real(r8), intent(in) :: accre_enhan(mgncol,nlev) ! optional accretion
                                             ! enhancement factor (-)

  real(r8), intent(in) :: p(mgncol,nlev)        ! air pressure (pa)
  real(r8), intent(in) :: pdel(mgncol,nlev)     ! pressure difference across level (pa)
  real(r8), intent(in) :: pint(mgncol,nlev+1)   ! pressure at interfaces

  real(r8), intent(in) :: cldn(mgncol,nlev)      ! cloud fraction (no units)
  real(r8), intent(in) :: liqcldf(mgncol,nlev)   ! liquid cloud fraction (no units)
  real(r8), intent(in) :: icecldf(mgncol,nlev)   ! ice cloud fraction (no units)
  real(r8), intent(in) :: qsatfac(mgncol,nlev)   ! subgrid cloud water saturation scaling factor (no units)

  ! used for scavenging
  ! Inputs for aerosol activation
  real(r8), intent(in) :: naai(mgncol,nlev)     ! ice nucleation number (from microp_aero_ts) (1/kg*s)
  real(r8), intent(in) :: npccn(mgncol,nlev)   ! ccn activated number tendency (from microp_aero_ts) (1/kg*s)

  ! Note that for these variables, the dust bin is assumed to be the last index.
  ! (For example, in CAM, the last dimension is always size 4.)
  real(r8), intent(in) :: rndst(:,:,:)  ! radius of each dust bin, for contact freezing (from microp_aero_ts) (m)
  real(r8), intent(in) :: nacon(:,:,:) ! number in each dust bin, for contact freezing  (from microp_aero_ts) (1/m^3)

  ! output arguments

  real(r8), intent(out) :: qcsinksum_rate1ord(mgncol,nlev) ! 1st order rate for
  ! direct cw to precip conversion
  real(r8), intent(out) :: tlat(mgncol,nlev)         ! latent heating rate       (W/kg)
  real(r8), intent(out) :: qvlat(mgncol,nlev)        ! microphysical tendency qv (1/s)
  real(r8), intent(out) :: qctend(mgncol,nlev)       ! microphysical tendency qc (1/s)
  real(r8), intent(out) :: qitend(mgncol,nlev)       ! microphysical tendency qi (1/s)
  real(r8), intent(out) :: nctend(mgncol,nlev)       ! microphysical tendency nc (1/(kg*s))
  real(r8), intent(out) :: nitend(mgncol,nlev)       ! microphysical tendency ni (1/(kg*s))

  real(r8), intent(out) :: qrtend(mgncol,nlev)       ! microphysical tendency qr (1/s)
  real(r8), intent(out) :: qstend(mgncol,nlev)       ! microphysical tendency qs (1/s)
  real(r8), intent(out) :: nrtend(mgncol,nlev)       ! microphysical tendency nr (1/(kg*s))
  real(r8), intent(out) :: nstend(mgncol,nlev)       ! microphysical tendency ns (1/(kg*s))
  real(r8), intent(out) :: qgtend(mgncol,nlev)       ! microphysical tendency qg (1/s)
  real(r8), intent(out) :: ngtend(mgncol,nlev)       ! microphysical tendency ng (1/(kg*s))

  real(r8), intent(out) :: effc(mgncol,nlev)         ! droplet effective radius (micron)
  real(r8), intent(out) :: effc_fn(mgncol,nlev)      ! droplet effective radius, assuming nc = 1.e8 kg-1
  real(r8), intent(out) :: effi(mgncol,nlev)         ! cloud ice effective radius (micron)
  real(r8), intent(out) :: sadice(mgncol,nlev)       ! cloud ice surface area density (cm2/cm3)
  real(r8), intent(out) :: sadsnow(mgncol,nlev)      ! cloud snow surface area density (cm2/cm3)
  real(r8), intent(out) :: prect(mgncol)             ! surface precip rate (m/s)
  real(r8), intent(out) :: preci(mgncol)             ! cloud ice/snow precip rate (m/s)
  real(r8), intent(out) :: nevapr(mgncol,nlev)       ! evaporation rate of rain + snow (kg/kg/s)
  real(r8), intent(out) :: am_evp_st(mgncol,nlev)    ! stratiform evaporation area (frac)
  real(r8), intent(out) :: prain(mgncol,nlev)        ! production of rain + snow (kg/kg/s)
  real(r8), intent(out) :: cmeout(mgncol,nlev)       ! Rate of cond-evap of ice (kg/kg/s)
  real(r8), intent(out) :: deffi(mgncol,nlev)        ! ice effective diameter for optics (radiation) (micron)
  real(r8), intent(out) :: pgamrad(mgncol,nlev)      ! ice gamma parameter for optics (radiation) (no units)
  real(r8), intent(out) :: lamcrad(mgncol,nlev)      ! slope of droplet distribution for optics (radiation) (1/m)
  real(r8), intent(out) :: qsout(mgncol,nlev)        ! snow mixing ratio (kg/kg)
  real(r8), intent(out) :: dsout(mgncol,nlev)        ! snow diameter (m)
  real(r8), intent(out) :: lflx(mgncol,nlev+1)       ! grid-box average liquid condensate flux (kg m^-2 s^-1)
  real(r8), intent(out) :: iflx(mgncol,nlev+1)       ! grid-box average ice condensate flux (kg m^-2 s^-1)
  real(r8), intent(out) :: rflx(mgncol,nlev+1)       ! grid-box average rain flux (kg m^-2 s^-1)
  real(r8), intent(out) :: sflx(mgncol,nlev+1)       ! grid-box average snow flux (kg m^-2 s^-1)
  real(r8), intent(out) :: gflx(mgncol,nlev+1)       ! grid-box average graupel/hail flux (kg m^-2 s^-1)

  real(r8), intent(out) :: qrout(mgncol,nlev)        ! grid-box average rain mixing ratio (kg/kg)
  real(r8), intent(out) :: reff_rain(mgncol,nlev)    ! rain effective radius (micron)
  real(r8), intent(out) :: reff_snow(mgncol,nlev)    ! snow effective radius (micron)
  real(r8), intent(out) :: reff_grau(mgncol,nlev)    ! graupel effective radius (micron)

  real(r8), intent(out) :: nrout(mgncol,nlev)        ! rain number concentration (1/m3)
  real(r8), intent(out) :: nsout(mgncol,nlev)        ! snow number concentration (1/m3)
  real(r8), intent(out) :: refl(mgncol,nlev)         ! analytic radar reflectivity (94GHZ, cloud radar) (dBZ)
  real(r8), intent(out) :: arefl(mgncol,nlev)        ! average reflectivity will zero points outside valid range (dBZ)
  real(r8), intent(out) :: areflz(mgncol,nlev)       ! average reflectivity in z. (mm6 m-3)
  real(r8), intent(out) :: frefl(mgncol,nlev)        ! fractional occurrence of radar reflectivity
  real(r8), intent(out) :: csrfl(mgncol,nlev)        ! cloudsat reflectivity (dBZ)
  real(r8), intent(out) :: acsrfl(mgncol,nlev)       ! cloudsat average (dBZ)
  real(r8), intent(out) :: fcsrfl(mgncol,nlev)       ! cloudsat fractional occurrence of radar reflectivity
  real(r8), intent(out) :: refl10cm(mgncol,nlev)     ! 10cm (rain) analytic radar reflectivity (dBZ)
  real(r8), intent(out) :: reflz10cm(mgncol,nlev)    ! 10cm (rain) analytic radar reflectivity (mm6 m-3)
  real(r8), intent(out) :: rercld(mgncol,nlev)       ! effective radius calculation for rain + cloud
  real(r8), intent(out) :: ncai(mgncol,nlev)         ! output number conc of ice nuclei available (1/m3)
  real(r8), intent(out) :: ncal(mgncol,nlev)         ! output number conc of CCN (1/m3)
  real(r8), intent(out) :: qrout2(mgncol,nlev)       ! copy of qrout as used to compute drout2
  real(r8), intent(out) :: qsout2(mgncol,nlev)       ! copy of qsout as used to compute dsout2
  real(r8), intent(out) :: nrout2(mgncol,nlev)       ! copy of nrout as used to compute drout2
  real(r8), intent(out) :: nsout2(mgncol,nlev)       ! copy of nsout as used to compute dsout2
  real(r8), intent(out) :: drout2(mgncol,nlev)       ! mean rain particle diameter (m)
  real(r8), intent(out) :: dsout2(mgncol,nlev)       ! mean snow particle diameter (m)
  real(r8), intent(out) :: freqs(mgncol,nlev)        ! fractional occurrence of snow
  real(r8), intent(out) :: freqr(mgncol,nlev)        ! fractional occurrence of rain
  real(r8), intent(out) :: nfice(mgncol,nlev)        ! fraction of frozen water to total condensed water
  real(r8), intent(out) :: qcrat(mgncol,nlev)        ! limiter for qc process rates (1=no limit --> 0. no qc)
  real(r8), intent(out) :: qgout(mgncol,nlev)        ! graupel/hail mixing ratio (kg/kg)
  real(r8), intent(out) :: dgout(mgncol,nlev)        ! graupel/hail diameter (m)
  real(r8), intent(out) :: ngout(mgncol,nlev)        ! graupel/hail number concentration (1/m3)
  real(r8), intent(out) :: qgout2(mgncol,nlev)       ! copy of qgout as used to compute dgout2
  real(r8), intent(out) :: ngout2(mgncol,nlev)       ! copy of ngout as used to compute dgout2
  real(r8), intent(out) :: dgout2(mgncol,nlev)       ! mean graupel/hail particle diameter (m)
  real(r8), intent(out) :: freqg(mgncol,nlev)        ! fractional occurrence of graupel

  real(r8), intent(out) :: prer_evap(mgncol,nlev)    ! evaporation rate of rain (kg/kg/s)

  type (proc_rates_type), intent(inout)  :: proc_rates

  character(128),   intent(out) :: errstring  ! output status (non-blank for error return)

  ! Tendencies calculated by external schemes that can replace MG's native
  ! process tendencies.

  ! Used with CARMA cirrus microphysics
  ! (or similar external microphysics model)
  real(r8), intent(in) :: tnd_qsnow(:,:) ! snow mass tendency (kg/kg/s)
  real(r8), intent(in) :: tnd_nsnow(:,:) ! snow number tendency (#/kg/s)
  real(r8), intent(in) :: re_ice(:,:)    ! ice effective radius (m)

  ! From external ice nucleation.
  real(r8), intent(in) :: frzimm(:,:) ! Number tendency due to immersion freezing (1/cm3)
  real(r8), intent(in) :: frzcnt(:,:) ! Number tendency due to contact freezing (1/cm3)
  real(r8), intent(in) :: frzdep(:,:) ! Number tendency due to deposition nucleation (1/cm3)

  ! local workspace
  ! all units mks unless otherwise stated

  ! local copies of input variables
  real(r8) :: cldliq_mass_mixing_ratio_kg_per_kg(mgncol,nlev)      ! cloud liquid mixing ratio (kg/kg)
  real(r8) :: cldice_mass_mixing_ratio_kg_per_kg(mgncol,nlev)      ! cloud ice mixing ratio (kg/kg)
  real(r8) :: cldliq_number_conc_per_kg(mgncol,nlev)      ! cloud liquid number concentration (1/kg)
  real(r8) :: cldice_number_conc_per_kg(mgncol,nlev)      ! cloud liquid number concentration (1/kg)
  real(r8) :: rain_mass_mixing_ratio_kg_per_kg(mgncol,nlev)      ! rain mixing ratio (kg/kg)
  real(r8) :: snow_mass_mixing_ratio_kg_per_kg(mgncol,nlev)      ! snow mixing ratio (kg/kg)
  real(r8) :: rain_number_conc_per_kg(mgncol,nlev)      ! rain number concentration (1/kg)
  real(r8) :: snow_number_conc_per_kg(mgncol,nlev)      ! snow number concentration (1/kg)
  real(r8) :: graupel_mass_mixing_ratio_kg_per_kg(mgncol,nlev)      ! graupel mixing ratio (kg/kg)
  real(r8) :: graupel_number_conc_per_kg(mgncol,nlev)      ! graupel number concentration (1/kg)
  real(r8) :: rhogtmp              ! hail or graupel density (kg m-3)

  ! general purpose variables
  real(r8) :: deltat            ! sub-time step (s)
  real(r8) :: rdeltat           ! reciprocal of sub-time step (1/s)

  ! physical properties of the air at a given point
  real(r8) :: air_density_kg_per_m3(mgncol,nlev)    ! density (kg m-3)
  real(r8) :: dv(mgncol,nlev)     ! diffusivity of water vapor
  real(r8) :: air_dynamic_viscosity_kg_per_m_per_s(mgncol,nlev)     ! viscosity
  real(r8) :: schmidt_number(mgncol,nlev)     ! schmidt number
  real(r8) :: fall_speed_air_density_correction(mgncol,nlev)   ! density correction factor for fallspeed

  ! cloud fractions
  real(r8) :: precipitation_area_fraction(mgncol,nlev) ! precip fraction assuming maximum overlap
  real(r8) :: cloud_area_fraction(mgncol,nlev)   ! cloud fraction
  real(r8) :: ice_cloud_area_fraction(mgncol,nlev)  ! ice cloud fraction
  real(r8) :: liquid_cloud_area_fraction(mgncol,nlev)  ! liq cloud fraction
  real(r8) :: subgrid_saturation_scaling_factor(mgncol,nlev)   ! subgrid cloud water saturation scaling factor

  ! mass mixing ratios
  real(r8) :: incloud_cldliq_mass_mixing_ratio_kg_per_kg(mgncol,nlev)   ! in-cloud cloud liquid
  real(r8) :: incloud_cldice_mass_mixing_ratio_kg_per_kg(mgncol,nlev)   ! in-cloud cloud ice
  real(r8) :: inprecip_snow_mass_mixing_ratio_kg_per_kg(mgncol,nlev)   ! in-precip snow
  real(r8) :: inprecip_rain_mass_mixing_ratio_kg_per_kg(mgncol,nlev)   ! in-precip rain
  real(r8) :: inprecip_graupel_mass_mixing_ratio_kg_per_kg(mgncol,nlev)   ! in-precip graupel/hail

  ! number concentrations
  real(r8) :: incloud_cldliq_number_conc_per_kg(mgncol,nlev)   ! in-cloud droplet
  real(r8) :: incloud_cldice_number_conc_per_kg(mgncol,nlev)   ! in-cloud cloud ice
  real(r8) :: inprecip_snow_number_conc_per_kg(mgncol,nlev)   ! in-precip snow
  real(r8) :: inprecip_rain_number_conc_per_kg(mgncol,nlev)   ! in-precip rain
  real(r8) :: inprecip_graupel_number_conc_per_kg(mgncol,nlev)   ! in-precip graupel/hail

  ! Size distribution parameters for:
  ! cloud ice
  real(r8) :: cldice_size_dist_slope_per_m(mgncol,nlev)   ! slope
  real(r8) :: cldice_size_dist_intercept_per_m4(mgncol,nlev)    ! intercept
  ! cloud liquid
  real(r8) :: cldliq_size_dist_slope_per_m(mgncol,nlev)   ! slope
  real(r8) :: cldliq_size_dist_shape(mgncol,nlev)   ! spectral width parameter
  ! snow
  real(r8) :: snow_size_dist_slope_per_m(mgncol,nlev)   ! slope
  real(r8) :: snow_size_dist_intercept_per_m4(mgncol,nlev)    ! intercept
  ! rain
  real(r8) :: rain_size_dist_slope_per_m(mgncol,nlev)   ! slope
  real(r8) :: rain_size_dist_intercept_per_m4(mgncol,nlev)    ! intercept
  ! graupel/hail
  real(r8) :: graupel_size_dist_slope_per_m(mgncol,nlev)   ! slope
  real(r8) :: graupel_size_dist_intercept_per_m4(mgncol,nlev)    ! intercept
  real(r8) :: bgtmp               ! tmp fall speed parameter

  ! Rates/tendencies due to:

  ! Instantaneous snow melting
  real(r8) :: snow_mass_change_instant_melting_kg_per_kg(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: snow_number_change_instant_melting_per_kg(mgncol,nlev)    ! number concentration
  ! Instantaneous graupel melting
  real(r8) :: graupel_mass_change_instant_melting_kg_per_kg(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: graupel_number_change_instant_melting_per_kg(mgncol,nlev)    ! number concentration

  ! Instantaneous rain freezing
  real(r8) :: rain_mass_change_instant_freezing_kg_per_kg(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: rain_number_change_instant_freezing_per_kg(mgncol,nlev)    ! number concentration

  ! deposition of cloud ice
  real(r8) :: local_cldice_mass_tend_deposition_kg_per_kg_per_s(mgncol,nlev)    ! deposition from vapor to ice PMC 12/3/12
  ! sublimation of cloud ice
  real(r8) :: local_cldice_mass_tend_sublimation_kg_per_kg_per_s(mgncol,nlev) ! sublimation from ice to vapor PMC 12/3/12
  ! vapor deposition onto
  real(r8) :: local_snow_mass_tend_deposition_kg_per_kg_per_s(mgncol,nlev) ! Vapor deposition onto snow.

  ! ice nucleation
  real(r8) :: local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(mgncol,nlev) ! number rate from deposition/cond.-freezing
  real(r8) :: local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(mgncol,nlev) ! mass mixing ratio
  ! freezing of cloud water
  real(r8) :: local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(mgncol,nlev) ! mass mixing ratio
  real(r8) :: local_cldliq_number_tend_immersion_freezing_per_kg_per_s(mgncol,nlev) ! number concentration
  ! contact freezing of cloud water
  real(r8) :: local_cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(mgncol,nlev) ! mass mixing ratio
  real(r8) :: local_cldliq_number_tend_contact_freezing_per_kg_per_s(mgncol,nlev) ! number concentration
  ! deposition nucleation in mixed-phase clouds (from external scheme)
  real(r8) :: local_cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s(mgncol,nlev) ! mass mixing ratio
  real(r8) :: local_cldliq_number_tend_deposition_nucleation_per_kg_per_s(mgncol,nlev) ! number concentration
  ! ice multiplication
  real(r8) :: local_cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s(mgncol,nlev) ! mass mixing ratio
  real(r8) :: local_cldice_number_tend_splintering_by_snow_per_kg_per_s(mgncol,nlev) ! number concentration
  ! autoconversion of cloud droplets
  real(r8) :: local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: local_rain_number_tend_autoconversion_from_cldliq_per_kg_per_s(mgncol,nlev)   ! number concentration (rain)
  real(r8) :: local_cldliq_number_tend_autoconversion_to_rain_per_kg_per_s(mgncol,nlev)  ! number concentration (cloud droplets)
  ! self-aggregation of snow
  real(r8) :: local_snow_number_tend_aggregation_per_kg_per_s(mgncol,nlev)  ! number concentration
  ! self-collection of rain
  real(r8) :: local_rain_number_tend_self_collection_per_kg_per_s(mgncol,nlev)  ! number concentration
  ! collection of droplets by snow
  real(r8) :: local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(mgncol,nlev)     ! mass mixing ratio
  real(r8) :: local_cldliq_number_tend_riming_by_snow_per_kg_per_s(mgncol,nlev)    ! number concentration
  ! collection of rain by snow
  real(r8) :: local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(mgncol,nlev)  ! mass mixing ratio
  real(r8) :: local_rain_number_tend_accretion_by_snow_per_kg_per_s(mgncol,nlev) ! number concentration
  ! freezing of rain
  real(r8) :: local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(mgncol,nlev) ! mass mixing ratio
  real(r8) :: local_rain_number_tend_freezing_to_snow_per_kg_per_s(mgncol,nlev) ! number concentration
  ! freezing of rain to form ice (mg add 4/26/13)
  real(r8) :: local_rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: local_rain_number_tend_freezing_to_cldice_per_kg_per_s(mgncol,nlev)    ! number concentration
  ! accretion of droplets by rain
  real(r8) :: local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: local_cldliq_number_tend_accretion_by_rain_per_kg_per_s(mgncol,nlev)   ! number concentration
  ! autoconversion of cloud ice to snow
  real(r8) :: local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(mgncol,nlev)   ! mass mixing ratio
  real(r8) :: local_cldice_number_tend_autoconversion_to_snow_per_kg_per_s(mgncol,nlev)  ! number concentration
  ! accretion of cloud ice by snow
  real(r8) :: local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(mgncol,nlev)   ! mass mixing ratio
  real(r8) :: local_cldice_number_tend_accretion_by_snow_per_kg_per_s(mgncol,nlev)  ! number concentration
  ! evaporation of rain
  real(r8) :: local_rain_mass_tend_evaporation_kg_per_kg_per_s(mgncol,nlev)    ! mass mixing ratio
  ! sublimation of snow
  real(r8) :: local_snow_mass_tend_sublimation_kg_per_kg_per_s(mgncol,nlev)   ! mass mixing ratio
  ! number evaporation
  real(r8) :: local_cldice_number_tend_sublimation_per_kg_per_s(mgncol,nlev)  ! cloud ice
  real(r8) :: local_cldliq_number_tend_evaporation_per_kg_per_s(mgncol,nlev)  ! droplet
  real(r8) :: local_snow_number_tend_sublimation_per_kg_per_s(mgncol,nlev)  ! snow
  real(r8) :: local_rain_number_tend_evaporation_per_kg_per_s(mgncol,nlev)  ! rain
  ! bergeron process
  real(r8) :: local_cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s(mgncol,nlev)   ! mass mixing ratio (cloud ice)
  real(r8) :: local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(mgncol,nlev)  ! mass mixing ratio (snow)

  !graupel/hail processes
  real(r8) :: local_rain_number_tend_accretion_by_graupel_per_kg_per_s(mgncol,nlev)  ! change n collection rain by graupel  (precipf)
  real(r8) :: local_graupel_number_tend_from_snow_riming_per_kg_per_s(mgncol,nlev)   ! change n conversion to graupel due to collection droplets by snow (lcldm)
  real(r8) :: local_graupel_number_tend_rain_snow_collection_per_kg_per_s(mgncol,nlev)  ! change n conversion to graupel due to collection rain by snow (precipf)
  real(r8) :: local_cldice_number_tend_splintering_by_graupel_per_kg_per_s(mgncol,nlev)  ! ice mult due to acc droplets by graupel  (lcldm)
  real(r8) :: local_cldice_number_tend_splintering_rain_graupel_per_kg_per_s(mgncol,nlev) ! ice mult due to acc rain by graupel  (precipf)
  real(r8) :: local_cldliq_number_tend_riming_by_graupel_per_kg_per_s(mgncol,nlev) ! change n collection droplets by graupel (lcldm)

  real(r8) :: local_snow_mass_tend_collection_by_rain_kg_per_kg_per_s(mgncol,nlev)   ! conversion due to coll of snow by rain (precipf)
  real(r8) :: local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(mgncol,nlev)   ! change in q collection rain by graupel  (precipf)
  real(r8) :: local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(mgncol,nlev)  ! change in q collection droplets by graupel (lcldm)
  real(r8) :: local_graupel_mass_tend_from_snow_riming_kg_per_kg_per_s(mgncol,nlev)  ! conversion q to graupel due to collection droplets by snow  (lcldm)
  real(r8) :: local_graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s(mgncol,nlev)  ! conversion q to graupel due to collection rain by snow (precipf)
  real(r8) :: local_graupel_mass_tend_deposition_kg_per_kg_per_s(mgncol,nlev)    ! dep of graupel (precipf)
  real(r8) :: local_cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s(mgncol,nlev)  ! change q due to ice mult droplets/graupel  (lcldm)
  real(r8) :: local_cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s(mgncol,nlev) ! change q due to ice mult rain/graupel (precipf)


  ! fallspeeds
  ! number-weighted
  real(r8) :: snow_number_weighted_fall_speed_m_per_s(mgncol,nlev)    ! snow
  real(r8) :: rain_number_weighted_fall_speed_m_per_s(mgncol,nlev)    ! rain
  real(r8) :: graupel_number_weighted_fall_speed_m_per_s(mgncol,nlev)    ! graupel/hail

  ! air density corrected fallspeed parameters
  real(r8) :: rain_fall_speed_prefactor(mgncol,nlev)    ! rain
  real(r8) :: snow_fall_speed_prefactor(mgncol,nlev)    ! snow
  real(r8) :: graupel_fall_speed_prefactor(mgncol,nlev)    ! graupel
  real(r8) :: cldliq_fall_speed_prefactor(mgncol,nlev)    ! cloud droplet
  real(r8) :: cldice_fall_speed_prefactor(mgncol,nlev)    ! cloud ice
  real(r8) :: small_cldice_fall_speed_prefactor(mgncol,nlev)    ! cloud small ice

  ! Mass of liquid droplets used with external heterogeneous freezing.
  real(r8) :: mi0l(mgncol,nlev)

  ! saturation vapor pressures
  real(r8) :: saturation_vapor_pressure_over_liquid_Pa(mgncol,nlev)    ! liquid
  real(r8) :: saturation_vapor_pressure_over_ice_Pa(mgncol,nlev)    ! ice
  real(r8) :: esnA(mgncol,nlev)   ! checking for RH after rain evap

  ! saturation vapor mixing ratios
  real(r8) :: saturation_specific_humidity_over_liquid_kg_per_kg(mgncol,nlev)    ! liquid
  real(r8) :: saturation_specific_humidity_over_ice_kg_per_kg(mgncol,nlev)    ! ice
  real(r8) :: qvnA(mgncol,nlev), qvnAI(mgncol,nlev) ! checking for RH after rain evap

  ! relative humidity
  real(r8) :: relative_humidity(mgncol,nlev)

  ! parameters for cloud water and cloud ice sedimentation calculations
  real(r8) :: cldliq_mass_pressure_fall_speed_Pa_per_s(mgncol,nlev)
  real(r8) :: cldliq_number_pressure_fall_speed_Pa_per_s(mgncol,nlev)
  real(r8) :: cldice_mass_pressure_fall_speed_Pa_per_s(mgncol,nlev)
  real(r8) :: cldice_number_pressure_fall_speed_Pa_per_s(mgncol,nlev)
  real(r8) :: graupel_mass_pressure_fall_speed_Pa_per_s(mgncol,nlev)
  real(r8) :: graupel_number_pressure_fall_speed_Pa_per_s(mgncol,nlev)
  real(r8) :: rain_mass_pressure_fall_speed_Pa_per_s(mgncol,nlev)
  real(r8) :: rain_number_pressure_fall_speed_Pa_per_s(mgncol,nlev)
  real(r8) :: snow_mass_pressure_fall_speed_Pa_per_s(mgncol,nlev)
  real(r8) :: snow_number_pressure_fall_speed_Pa_per_s(mgncol,nlev)

  real(r8) :: rthrsh     ! rain rate threshold for reflectivity calculation

  ! dummy variables
  real(r8) :: dum, dum1, dum2, dum3, dum4, qtmp
  real(r8) :: dum1A(mgncol,nlev), dum2A(mgncol,nlev), dum3A(mgncol,nlev)
  real(r8) :: dumni0, dumni0A2D(mgncol,nlev)
  real(r8) :: dumns0, dumns0A2D(mgncol,nlev)
  ! dummies for checking RH
  real(r8) :: ttmpA(mgncol,nlev), qtmpAI(mgncol,nlev)
  ! dummies for conservation check
  real(r8) :: ratio, tmpnr,tmpp
  real(r8) :: tmpfrz
  ! dummies for in-cloud variables
  real(r8) :: dumc(mgncol,nlev)   ! qc
  real(r8) :: dumnc(mgncol,nlev)  ! nc
  real(r8) :: dumi(mgncol,nlev)   ! qi
  real(r8) :: dumni(mgncol,nlev)  ! ni
  real(r8) :: dumr(mgncol,nlev)   ! rain mixing ratio
  real(r8) :: dumnr(mgncol,nlev)  ! rain number concentration
  real(r8) :: dums(mgncol,nlev)   ! snow mixing ratio
  real(r8) :: dumns(mgncol,nlev)  ! snow number concentration
  real(r8) :: dumg(mgncol,nlev)   ! graupel mixing ratio
  real(r8) :: dumng(mgncol,nlev)  ! graupel number concentration
  ! Array dummy variable
  real(r8) :: dum_2D(mgncol,nlev)
  real(r8) :: pdel_inv(mgncol,nlev)

  ! loop array variables
  ! "i" and "k" are column/level iterators for internal (MG) variables
  ! "n" is used for other looping (currently just sedimentation)
  integer i, k, n

  integer mdust
  integer :: precip_frac_method

  ! Varaibles to scale fall velocity between small and regular ice regimes.
  real(r8) :: irad
  real(r8) :: ifrac

  !Variables for accretion seeing autoconverted liquid
  real(r8) :: rtmp(mgncol,nlev) ! dummy for rain + autoconversion
  real(r8) :: ctmp(mgncol,nlev) ! dummy for liq - autoconversion
  real(r8) :: ntmp(mgncol,nlev) ! dummy for liq - autoconversion number

  ! Variables for height calculation (used in Implicit Fall Speed)
  real(r8) :: zint(mgncol,nlev+1) ! interface height
  real(r8) :: H   !Scale height

  ! temporary local variables for asynchronous GPU run
  ! ice
  real(r8) :: prect_i(mgncol)
  real(r8) :: tlat_i(mgncol,nlev)
  real(r8) :: qvlat_i(mgncol,nlev)
  real(r8) :: preci_i(mgncol)
  ! liq
  real(r8) :: prect_l(mgncol)
  real(r8) :: tlat_l(mgncol,nlev)
  real(r8) :: qvlat_l(mgncol,nlev)
  ! rain
  real(r8) :: prect_r(mgncol)
  ! snow
  real(r8) :: prect_s(mgncol)
  real(r8) :: preci_s(mgncol)
  ! graupel
  real(r8) :: prect_g(mgncol)
  real(r8) :: preci_g(mgncol)

  ! number of sub-steps for loops over "n" (for sedimentation)
  ! ice
  integer nstep_i(mgncol)
  real(r8) :: rnstep_i(mgncol)
  ! liq
  integer nstep_l(mgncol)
  real(r8) :: rnstep_l(mgncol)
  ! rain
  integer nstep_r(mgncol)
  real(r8) :: rnstep_r(mgncol)
  ! snow
  integer nstep_s(mgncol)
  real(r8) :: rnstep_s(mgncol)
  ! graupel
  integer nstep_g(mgncol)
  real(r8) :: rnstep_g(mgncol)

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! Initialize scale height (H) for interface height calculation
  ! needed for Implicit Fall Speed
  H=0._r8

  ! Return error message
  errstring = ' '

  ! Process inputs

  ! assign variable deltat to deltatin
  deltat  = deltatin
  rdeltat = 1._r8 / deltat

  if (trim(micro_mg_precip_frac_method) == 'in_cloud') then
     precip_frac_method = MG_PRECIP_FRAC_INCLOUD
  else if(trim(micro_mg_precip_frac_method) == 'max_overlap') then
     precip_frac_method = MG_PRECIP_FRAC_OVERLAP
  endif

  !......................................................................
  !       graupel/hail density set (Hail = 400, Graupel = 500 from M2005)
  bgtmp=0._r8
  rhogtmp=0._r8
  if (do_hail) then
     bgtmp = bh
     rhogtmp = rhoh
  end if
  if (do_graupel) then
     bgtmp = bg
     rhogtmp = rhog
  end if

  ! set mdust as the number of dust bins for use later in contact freezing subroutine
  mdust = size(rndst,3)

  !$acc data copyin  (t,q,qcn,qin,ncn,nin,qrn,qsn,nrn,nsn,qgr,ngr,relvar,     &
  !$acc               accre_enhan,p,pdel,pint,cldn,liqcldf,icecldf,qsatfac,   &
  !$acc               naai,npccn,rndst,nacon,tnd_qsnow,tnd_nsnow,re_ice,      &
  !$acc               frzimm,frzcnt,frzdep,mg_liq_props,mg_ice_props,         &
  !$acc               mg_rain_props,mg_graupel_props,mg_hail_props,           &
  !$acc               mg_snow_props,proc_rates)                               &
  !$acc      copyout (qcsinksum_rate1ord,tlat,qvlat,qctend,qitend,nctend,     &
  !$acc               nitend,qrtend,qstend,nrtend,nstend,qgtend,ngtend,       &
  !$acc               effc,effc_fn,effi,sadice,sadsnow,prect,preci,           &
  !$acc               nevapr,proc_rates%snow_mass_tend_sublimation_kg_per_kg_per_s,am_evp_st,prain,             &
  !$acc               proc_rates%snow_mass_tend_production_kg_per_kg_per_s,cmeout,                             &
  !$acc               deffi,pgamrad,lamcrad,qsout,dsout,lflx,iflx,rflx,       &
  !$acc               sflx,gflx,qrout,reff_rain,reff_snow,reff_grau,          &
  !$acc               proc_rates%cldliq_mass_tend_sedimentation_evaporation_kg_per_kg_per_s,proc_rates%cldice_mass_tend_sedimentation_sublimation_kg_per_kg_per_s,proc_rates%vapor_mass_tend_residual_condensation_kg_per_kg_per_s, &
  !$acc               proc_rates%cldice_mass_tend_net_deposition_kg_per_kg_per_s,proc_rates%cldliq_mass_weighted_fall_speed_m_per_s,proc_rates%cldice_mass_weighted_fall_speed_m_per_s,   &
  !$acc               proc_rates%rain_mass_weighted_fall_speed_m_per_s,proc_rates%snow_mass_weighted_fall_speed_m_per_s,                          &
  !$acc               proc_rates%graupel_mass_weighted_fall_speed_m_per_s,proc_rates%graupel_mass_tend_sedimentation_kg_per_kg_per_s,proc_rates%cldliq_mass_tend_sedimentation_kg_per_kg_per_s, &
  !$acc               proc_rates%cldice_mass_tend_sedimentation_kg_per_kg_per_s,proc_rates%rain_mass_tend_sedimentation_kg_per_kg_per_s,                &
  !$acc               proc_rates%snow_mass_tend_sedimentation_kg_per_kg_per_s,proc_rates%cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s,                  &
  !$acc               proc_rates%cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s,proc_rates%cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s,                 &
  !$acc               proc_rates%cldliq_mass_tend_contact_freezing_kg_per_kg_per_s,proc_rates%cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s,              &
  !$acc               proc_rates%cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s,proc_rates%cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s,               &
  !$acc               proc_rates%snow_mass_tend_deposition_kg_per_kg_per_s,proc_rates%cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s,               &
  !$acc               proc_rates%cldice_mass_tend_melting_kg_per_kg_per_s,proc_rates%snow_mass_tend_melting_kg_per_kg_per_s,                 &
  !$acc               proc_rates%graupel_mass_tend_melting_kg_per_kg_per_s,proc_rates%cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s,               &
  !$acc               proc_rates%cldliq_mass_tend_homogeneous_freezing_kg_per_kg_per_s,                                     &
  !$acc               proc_rates%cldliq_mass_tend_residual_condensation_kg_per_kg_per_s,proc_rates%cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s,                 &
  !$acc               proc_rates%cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s,proc_rates%cldice_mass_tend_residual_deposition_kg_per_kg_per_s,                 &
  !$acc               proc_rates%rain_mass_tend_freezing_to_snow_kg_per_kg_per_s,proc_rates%rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s,             &
  !$acc               proc_rates%rain_mass_tend_accretion_by_snow_kg_per_kg_per_s,proc_rates%air_heating_rate_snow_melting_W_per_kg,              &
  !$acc               proc_rates%air_heating_rate_rain_freezing_W_per_kg,proc_rates%cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s,              &
  !$acc               proc_rates%rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s,proc_rates%cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s,               &
  !$acc               proc_rates%graupel_mass_tend_from_snow_riming_kg_per_kg_per_s,proc_rates%graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s,              &
  !$acc               proc_rates%graupel_mass_tend_deposition_kg_per_kg_per_s,proc_rates%cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s,                &
  !$acc               proc_rates%cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s,proc_rates%snow_mass_tend_collection_by_rain_kg_per_kg_per_s,              &
  !$acc               proc_rates%rain_number_tend_accretion_by_graupel_per_kg_per_s,proc_rates%graupel_number_tend_from_snow_riming_per_kg_per_s,               &
  !$acc               proc_rates%graupel_number_tend_rain_snow_collection_per_kg_per_s,proc_rates%cldice_number_tend_splintering_by_graupel_per_kg_per_s,              &
  !$acc               proc_rates%cldice_number_tend_splintering_rain_graupel_per_kg_per_s,proc_rates%cldliq_number_tend_riming_by_graupel_per_kg_per_s,            &
  !$acc               nrout,nsout,refl,arefl,                                 &
  !$acc               areflz,frefl,csrfl,acsrfl,fcsrfl,refl10cm,reflz10cm,    &
  !$acc               rercld,ncai,ncal,qrout2,qsout2,nrout2,nsout2,drout2,    &
  !$acc               dsout2,freqs,freqr,nfice,qcrat,qgout,dgout,ngout,       &
  !$acc               qgout2,ngout2,dgout2,freqg,prer_evap,                   &
  !$acc               proc_rates%cldliq_number_tend_immersion_freezing_per_kg_per_s,proc_rates%cldliq_number_tend_contact_freezing_per_kg_per_s,              &
  !$acc               proc_rates%cldice_number_tend_nucleation_from_vapor_per_kg_per_s,proc_rates%cldliq_number_tend_deposition_nucleation_per_kg_per_s,              &
  !$acc               proc_rates%cldliq_number_tend_homogeneous_freezing_per_kg_per_s,proc_rates%rain_number_tend_freezing_to_snow_per_kg_per_s,               &
  !$acc               proc_rates%rain_number_tend_freezing_to_cldice_per_kg_per_s,proc_rates%cldice_number_tend_splintering_by_snow_per_kg_per_s,             &
  !$acc               proc_rates%cldliq_number_tend_accretion_by_rain_per_kg_per_s,proc_rates%cldliq_number_tend_riming_by_snow_per_kg_per_s,               &
  !$acc               proc_rates%cldice_number_tend_accretion_by_snow_per_kg_per_s,proc_rates%rain_number_tend_accretion_by_snow_per_kg_per_s,               &
  !$acc               proc_rates%rain_number_tend_autoconversion_from_cldliq_per_kg_per_s,proc_rates%cldice_number_tend_autoconversion_to_snow_per_kg_per_s,                 &
  !$acc               proc_rates%cldliq_number_tend_sedimentation_per_kg_per_s,proc_rates%cldice_number_tend_sedimentation_per_kg_per_s,                &
  !$acc               proc_rates%rain_number_tend_sedimentation_per_kg_per_s,proc_rates%snow_number_tend_sedimentation_per_kg_per_s,                &
  !$acc               proc_rates%graupel_number_tend_sedimentation_per_kg_per_s,proc_rates%cldice_number_tend_melting_per_kg_per_s,                &
  !$acc               proc_rates%snow_number_tend_melting_per_kg_per_s,proc_rates%graupel_number_tend_melting_per_kg_per_s,              &
  !$acc               proc_rates%rain_number_tend_self_collection_per_kg_per_s,proc_rates%scale_qc,                &
  !$acc               proc_rates%scale_nc,proc_rates%scale_qr,                &
  !$acc               proc_rates%scale_nr,proc_rates%amk_c,proc_rates%ank_c,  &
  !$acc               proc_rates%amk_r,proc_rates%ank_r,proc_rates%amk,       &
  !$acc               proc_rates%ank,proc_rates%amk_out,proc_rates%ank_out,   &
  !$acc               proc_rates%qc_out_TAU,proc_rates%nc_out_TAU,            &
  !$acc               proc_rates%qr_out_TAU,proc_rates%nr_out_TAU,            &
  !$acc               proc_rates%qc_in_TAU,proc_rates%nc_in_TAU,              &
  !$acc               proc_rates%qr_in_TAU,proc_rates%nr_in_TAU,              &
  !$acc               proc_rates%cldliq_size_dist_slope_output_per_m,proc_rates%rain_size_dist_slope_output_per_m,                &
  !$acc               proc_rates%cldliq_size_dist_shape_output,proc_rates%rain_size_dist_intercept_output_per_m4,                 &
  !$acc               proc_rates%cldliq_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s,proc_rates%cldliq_number_tend_autoconv_accretion_KK2000_per_kg_per_s,      &
  !$acc               proc_rates%rain_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s,proc_rates%rain_number_tend_autoconv_accretion_KK2000_per_kg_per_s,      &
  !$acc               proc_rates%cldliq_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s,proc_rates%cldliq_number_tend_autoconv_accretion_SB2001_per_kg_per_s,      &
  !$acc               proc_rates%rain_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s,proc_rates%rain_number_tend_autoconv_accretion_SB2001_per_kg_per_s,      &
  !$acc               proc_rates%cldliq_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s,proc_rates%cldliq_number_tend_autoconv_accretion_TAU_per_kg_per_s,            &
  !$acc               proc_rates%rain_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s,proc_rates%rain_number_tend_autoconv_accretion_TAU_per_kg_per_s,            &
  !$acc               proc_rates%gmnnn_lmnnn_TAU)                             &
  !$acc      create  (cldliq_mass_mixing_ratio_kg_per_kg,cldice_mass_mixing_ratio_kg_per_kg,cldliq_number_conc_per_kg,cldice_number_conc_per_kg,rain_mass_mixing_ratio_kg_per_kg,snow_mass_mixing_ratio_kg_per_kg,rain_number_conc_per_kg,snow_number_conc_per_kg,graupel_mass_mixing_ratio_kg_per_kg,graupel_number_conc_per_kg,air_density_kg_per_m3,dv,air_dynamic_viscosity_kg_per_m_per_s,schmidt_number,fall_speed_air_density_correction,        &
  !$acc               precipitation_area_fraction,cloud_area_fraction,ice_cloud_area_fraction,liquid_cloud_area_fraction,subgrid_saturation_scaling_factor,incloud_cldliq_mass_mixing_ratio_kg_per_kg,incloud_cldice_mass_mixing_ratio_kg_per_kg,inprecip_snow_mass_mixing_ratio_kg_per_kg,inprecip_rain_mass_mixing_ratio_kg_per_kg,  &
  !$acc               inprecip_graupel_mass_mixing_ratio_kg_per_kg,incloud_cldliq_number_conc_per_kg,incloud_cldice_number_conc_per_kg,inprecip_snow_number_conc_per_kg,inprecip_rain_number_conc_per_kg,inprecip_graupel_number_conc_per_kg,cldice_size_dist_slope_per_m,cldice_size_dist_intercept_per_m4,cldliq_size_dist_slope_per_m,cldliq_size_dist_shape,snow_size_dist_slope_per_m,  &
  !$acc               snow_size_dist_intercept_per_m4,rain_size_dist_slope_per_m,rain_size_dist_intercept_per_m4,graupel_size_dist_slope_per_m,graupel_size_dist_intercept_per_m4,snow_mass_change_instant_melting_kg_per_kg,snow_number_change_instant_melting_per_kg,graupel_mass_change_instant_melting_kg_per_kg,graupel_number_change_instant_melting_per_kg,  &
  !$acc               rain_mass_change_instant_freezing_kg_per_kg,rain_number_change_instant_freezing_per_kg,local_cldice_mass_tend_deposition_kg_per_kg_per_s,local_cldice_mass_tend_sublimation_kg_per_kg_per_s,local_snow_mass_tend_deposition_kg_per_kg_per_s,local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s,     &
  !$acc               local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s,local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s,local_cldliq_number_tend_immersion_freezing_per_kg_per_s,local_cldliq_mass_tend_contact_freezing_kg_per_kg_per_s,local_cldliq_number_tend_contact_freezing_per_kg_per_s,local_cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s,local_cldliq_number_tend_deposition_nucleation_per_kg_per_s,       &
  !$acc               local_cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s,local_cldice_number_tend_splintering_by_snow_per_kg_per_s,local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s,local_rain_number_tend_autoconversion_from_cldliq_per_kg_per_s,local_cldliq_number_tend_autoconversion_to_rain_per_kg_per_s,local_snow_number_tend_aggregation_per_kg_per_s,local_rain_number_tend_self_collection_per_kg_per_s,local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s,        &
  !$acc               local_cldliq_number_tend_riming_by_snow_per_kg_per_s,local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s,local_rain_number_tend_accretion_by_snow_per_kg_per_s,local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s,local_rain_number_tend_freezing_to_snow_per_kg_per_s,local_rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s,local_rain_number_tend_freezing_to_cldice_per_kg_per_s,local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s, &
  !$acc               local_cldliq_number_tend_accretion_by_rain_per_kg_per_s,local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s,local_cldice_number_tend_autoconversion_to_snow_per_kg_per_s,local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s,local_cldice_number_tend_accretion_by_snow_per_kg_per_s,local_rain_mass_tend_evaporation_kg_per_kg_per_s,local_snow_mass_tend_sublimation_kg_per_kg_per_s,local_cldice_number_tend_sublimation_per_kg_per_s,local_cldliq_number_tend_evaporation_per_kg_per_s,local_snow_number_tend_sublimation_per_kg_per_s,  &
  !$acc               local_rain_number_tend_evaporation_per_kg_per_s,local_cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s,local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s,local_rain_number_tend_accretion_by_graupel_per_kg_per_s,local_graupel_number_tend_from_snow_riming_per_kg_per_s,local_graupel_number_tend_rain_snow_collection_per_kg_per_s,local_cldice_number_tend_splintering_by_graupel_per_kg_per_s,local_cldice_number_tend_splintering_rain_graupel_per_kg_per_s,    &
  !$acc               local_cldliq_number_tend_riming_by_graupel_per_kg_per_s,local_snow_mass_tend_collection_by_rain_kg_per_kg_per_s,local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s,local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s,local_graupel_mass_tend_from_snow_riming_kg_per_kg_per_s,local_graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s,local_graupel_mass_tend_deposition_kg_per_kg_per_s,local_cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s,   &
  !$acc               local_cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s,snow_number_weighted_fall_speed_m_per_s,rain_number_weighted_fall_speed_m_per_s,graupel_number_weighted_fall_speed_m_per_s,rain_fall_speed_prefactor,snow_fall_speed_prefactor,graupel_fall_speed_prefactor,cldliq_fall_speed_prefactor,cldice_fall_speed_prefactor,small_cldice_fall_speed_prefactor,mi0l,saturation_vapor_pressure_over_liquid_Pa,   &
  !$acc               saturation_vapor_pressure_over_ice_Pa,esnA,saturation_specific_humidity_over_liquid_kg_per_kg,saturation_specific_humidity_over_ice_kg_per_kg,qvnA,qvnAI,relative_humidity,cldliq_mass_pressure_fall_speed_Pa_per_s,cldliq_number_pressure_fall_speed_Pa_per_s,cldice_mass_pressure_fall_speed_Pa_per_s,cldice_number_pressure_fall_speed_Pa_per_s,graupel_mass_pressure_fall_speed_Pa_per_s,    &
  !$acc               graupel_number_pressure_fall_speed_Pa_per_s,rain_mass_pressure_fall_speed_Pa_per_s,rain_number_pressure_fall_speed_Pa_per_s,snow_mass_pressure_fall_speed_Pa_per_s,snow_number_pressure_fall_speed_Pa_per_s,dum1A,dum2A,dum3A,dumni0A2D,          &
  !$acc               dumns0A2D,ttmpA,qtmpAI,dumc,dumnc,dumi,dumni,dumr,      &
  !$acc               dumnr,dums,dumns,dumg,dumng,dum_2D,pdel_inv,rtmp,ctmp,  &
  !$acc               ntmp,zint,nstep_i,rnstep_i,nstep_l,rnstep_l,nstep_r,    &
  !$acc               rnstep_r,nstep_s,rnstep_s,nstep_g,rnstep_g,prect_i,     &
  !$acc               tlat_i,qvlat_i,preci_i,prect_l,tlat_l,qvlat_l,prect_r,  &
  !$acc               prect_s,preci_s,prect_g,preci_g)

  ! Copies of input concentrations that may be changed internally.

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k = 1,nlev
     do i = 1,mgncol
        cldliq_mass_mixing_ratio_kg_per_kg(i,k) = qcn(i,k)
        cldliq_number_conc_per_kg(i,k) = ncn(i,k)
        cldice_mass_mixing_ratio_kg_per_kg(i,k) = qin(i,k)
        cldice_number_conc_per_kg(i,k) = nin(i,k)
        rain_mass_mixing_ratio_kg_per_kg(i,k) = qrn(i,k)
        rain_number_conc_per_kg(i,k) = nrn(i,k)
        snow_mass_mixing_ratio_kg_per_kg(i,k) = qsn(i,k)
        snow_number_conc_per_kg(i,k) = nsn(i,k)
        graupel_mass_mixing_ratio_kg_per_kg(i,k) = qgr(i,k)
        graupel_number_conc_per_kg(i,k) = ngr(i,k)
     end do
  end do
  !$acc end parallel

  ! cldn: used to set cldm, unused for subcolumns
  ! liqcldf: used to set lcldm, unused for subcolumns
  ! icecldf: used to set icldm, unused for subcolumns

  if (microp_uniform) then
     ! subcolumns, set cloud fraction variables to one
     ! if cloud water or ice is present, if not present
     ! set to mincld (mincld used instead of zero, to prevent
     ! possible division by zero errors).

     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
       do i=1,mgncol
          if (cldliq_mass_mixing_ratio_kg_per_kg(i,k) >= qsmall) then
             liquid_cloud_area_fraction(i,k) = 1._r8
          else
             liquid_cloud_area_fraction(i,k) = mincld
          end if

          if (cldice_mass_mixing_ratio_kg_per_kg(i,k) >= qsmall) then
             ice_cloud_area_fraction(i,k) = 1._r8
          else
             ice_cloud_area_fraction(i,k) = mincld
          end if

          cloud_area_fraction(i,k) = max(ice_cloud_area_fraction(i,k), liquid_cloud_area_fraction(i,k))
          subgrid_saturation_scaling_factor(i,k) = 1._r8
        end do
     end do
     !$acc end parallel
  else
     ! get cloud fraction, check for minimum

     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
          cloud_area_fraction(i,k) = max(cldn(i,k),mincld)
          liquid_cloud_area_fraction(i,k) = max(liqcldf(i,k),mincld)
          ice_cloud_area_fraction(i,k) = max(icecldf(i,k),mincld)
          subgrid_saturation_scaling_factor(i,k) = qsatfac(i,k)
        end do
     end do
     !$acc end parallel
  end if

  ! Initialize local variables

  ! local physical properties

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        air_density_kg_per_m3(i,k) = p(i,k)/(r*t(i,k))
        dv(i,k) = 8.794E-5_r8 * t(i,k)**1.81_r8 / p(i,k)
        air_dynamic_viscosity_kg_per_m_per_s(i,k) = 1.496E-6_r8 * t(i,k)**1.5_r8 / (t(i,k) + 120._r8)
        schmidt_number(i,k) = air_dynamic_viscosity_kg_per_m_per_s(i,k)/(air_density_kg_per_m3(i,k)*dv(i,k))

        ! air density adjustment for fallspeed parameters
        ! includes air density correction factor to the
        ! power of 0.54 following Heymsfield and Bansemer 2007

        fall_speed_air_density_correction(i,k)=(rhosu/air_density_kg_per_m3(i,k))**0.54_r8

        rain_fall_speed_prefactor(i,k)=ar*fall_speed_air_density_correction(i,k)
        snow_fall_speed_prefactor(i,k)=as*fall_speed_air_density_correction(i,k)
        ! Hail use ah*rhof graupel use ag*rhof
        ! Note that do_hail and do_graupel can't both be true
        if (do_hail) then
           graupel_fall_speed_prefactor(i,k) = ah*fall_speed_air_density_correction(i,k)
        end if
        if (do_graupel) then
           graupel_fall_speed_prefactor(i,k) = ag*fall_speed_air_density_correction(i,k)
        end if
        cldliq_fall_speed_prefactor(i,k)=g*rhow/(18._r8*air_dynamic_viscosity_kg_per_m_per_s(i,k))
        cldice_fall_speed_prefactor(i,k)=ai*(rhosu/air_density_kg_per_m3(i,k))**0.35_r8
        small_cldice_fall_speed_prefactor(i,k)=aj*(rhosu/air_density_kg_per_m3(i,k))**0.35_r8
     end do
  end do
  !$acc end parallel

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Get humidity and saturation vapor pressures

  call qsat_water(t, p, saturation_vapor_pressure_over_liquid_Pa, saturation_specific_humidity_over_liquid_kg_per_kg, mgncol*nlev)
  call qsat_ice(t, p, saturation_vapor_pressure_over_ice_Pa, saturation_specific_humidity_over_ice_kg_per_kg, mgncol*nlev)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        ! make sure when above freezing that esi=esl, not active yet
        if (t(i,k) >= tmelt) then
           saturation_vapor_pressure_over_ice_Pa(i,k)=saturation_vapor_pressure_over_liquid_Pa(i,k)
           saturation_specific_humidity_over_ice_kg_per_kg(i,k)=saturation_specific_humidity_over_liquid_kg_per_kg(i,k)
        else
           ! Scale the water saturation values to reflect subgrid scale
           ! ice cloud fraction, where ice clouds begin forming at a
           ! gridbox average relative humidity of rhmini (not 1).
           !
           ! NOTE: For subcolumns and other non-subgrid clouds, qsfm willi
           ! be 1.
           saturation_specific_humidity_over_ice_kg_per_kg(i,k) = subgrid_saturation_scaling_factor(i,k) * saturation_specific_humidity_over_ice_kg_per_kg(i,k)
           saturation_vapor_pressure_over_ice_Pa(i,k) = subgrid_saturation_scaling_factor(i,k) * saturation_vapor_pressure_over_ice_Pa(i,k)
           saturation_specific_humidity_over_liquid_kg_per_kg(i,k) = subgrid_saturation_scaling_factor(i,k) * saturation_specific_humidity_over_liquid_kg_per_kg(i,k)
           saturation_vapor_pressure_over_liquid_Pa(i,k) = subgrid_saturation_scaling_factor(i,k) * saturation_vapor_pressure_over_liquid_Pa(i,k)
        end if

        relative_humidity(i,k) = q(i,k) / max(saturation_specific_humidity_over_liquid_kg_per_kg(i,k), qsmall)

     end do
  end do
  !$acc end parallel

  ! initialize microphysics output

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        proc_rates%cldliq_mass_tend_sedimentation_evaporation_kg_per_kg_per_s(i,k)            = 0._r8
        proc_rates%cldice_mass_tend_sedimentation_sublimation_kg_per_kg_per_s(i,k)            = 0._r8
        proc_rates%vapor_mass_tend_residual_condensation_kg_per_kg_per_s(i,k)              = 0._r8
        proc_rates%cldice_mass_tend_net_deposition_kg_per_kg_per_s(i,k)            = 0._r8
        proc_rates%cldliq_mass_weighted_fall_speed_m_per_s(i,k)              = 0._r8
        proc_rates%cldice_mass_weighted_fall_speed_m_per_s(i,k)              = 0._r8
        proc_rates%cldliq_mass_tend_sedimentation_kg_per_kg_per_s(i,k)           = 0._r8
        proc_rates%cldice_mass_tend_sedimentation_kg_per_kg_per_s(i,k)           = 0._r8
        proc_rates%rain_mass_tend_sedimentation_kg_per_kg_per_s(i,k)           = 0._r8
        proc_rates%snow_mass_tend_sedimentation_kg_per_kg_per_s(i,k)           = 0._r8
        proc_rates%graupel_mass_tend_sedimentation_kg_per_kg_per_s(i,k)           = 0._r8

        proc_rates%cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)             = 0._r8
        proc_rates%cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k)             = 0._r8
        proc_rates%cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(i,k)          = 0._r8
        proc_rates%cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(i,k)          = 0._r8
        proc_rates%cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s(i,k)          = 0._r8
        proc_rates%cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k)          = 0._r8
        proc_rates%cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k)           = 0._r8
        proc_rates%snow_mass_tend_deposition_kg_per_kg_per_s(i,k)         = 0._r8
        proc_rates%cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s(i,k)            = 0._r8
        proc_rates%cldice_mass_tend_melting_kg_per_kg_per_s(i,k)            = 0._r8

        proc_rates%cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s(i,k)          = 0._r8
        proc_rates%snow_mass_tend_melting_kg_per_kg_per_s(i,k)           = 0._r8
        proc_rates%graupel_mass_tend_melting_kg_per_kg_per_s(i,k)           = 0._r8
        proc_rates%cldliq_mass_tend_homogeneous_freezing_kg_per_kg_per_s(i,k)            = 0._r8
        proc_rates%cldliq_mass_tend_residual_condensation_kg_per_kg_per_s(i,k)           = 0._r8
        proc_rates%cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(i,k)            = 0._r8
        proc_rates%cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)            = 0._r8
        proc_rates%cldice_mass_tend_residual_deposition_kg_per_kg_per_s(i,k)           = 0._r8
        proc_rates%rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k)          = 0._r8
        proc_rates%rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(i,k)         = 0._r8
        proc_rates%rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)           = 0._r8
        proc_rates%air_heating_rate_snow_melting_W_per_kg(i,k)         = 0._r8
        proc_rates%air_heating_rate_rain_freezing_W_per_kg(i,k)          = 0._r8
        proc_rates%cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k)          = 0._r8
        proc_rates%snow_mass_tend_collection_by_rain_kg_per_kg_per_s(i,k)           = 0._r8
        proc_rates%rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(i,k)           = 0._r8
        proc_rates%cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(i,k)          = 0._r8
        proc_rates%graupel_mass_tend_from_snow_riming_kg_per_kg_per_s(i,k)          = 0._r8
        proc_rates%graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s(i,k)          = 0._r8
        proc_rates%graupel_mass_tend_deposition_kg_per_kg_per_s(i,k)            = 0._r8
        proc_rates%cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s(i,k)          = 0._r8
        proc_rates%cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s(i,k)         = 0._r8
        proc_rates%rain_number_tend_accretion_by_graupel_per_kg_per_s(i,k)          = 0._r8
        proc_rates%graupel_number_tend_from_snow_riming_per_kg_per_s(i,k)           = 0._r8
        proc_rates%graupel_number_tend_rain_snow_collection_per_kg_per_s(i,k)          = 0._r8
        proc_rates%cldice_number_tend_splintering_by_graupel_per_kg_per_s(i,k)          = 0._r8
        proc_rates%cldice_number_tend_splintering_rain_graupel_per_kg_per_s(i,k)         = 0._r8
        proc_rates%cldliq_number_tend_riming_by_graupel_per_kg_per_s(i,k)         = 0._r8

        proc_rates%cldliq_number_tend_immersion_freezing_per_kg_per_s(i,k)          = 0._r8
        proc_rates%cldliq_number_tend_contact_freezing_per_kg_per_s(i,k)          = 0._r8
        proc_rates%cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k)          = 0._r8
        proc_rates%cldliq_number_tend_deposition_nucleation_per_kg_per_s(i,k)          = 0._r8
        proc_rates%cldliq_number_tend_homogeneous_freezing_per_kg_per_s(i,k)           = 0._r8
        proc_rates%rain_number_tend_freezing_to_snow_per_kg_per_s(i,k)          = 0._r8
        proc_rates%rain_number_tend_freezing_to_cldice_per_kg_per_s(i,k)         = 0._r8
        proc_rates%cldice_number_tend_splintering_by_snow_per_kg_per_s(i,k)          = 0._r8
        proc_rates%cldliq_number_tend_accretion_by_rain_per_kg_per_s(i,k)            = 0._r8
        proc_rates%cldliq_number_tend_riming_by_snow_per_kg_per_s(i,k)         = 0._r8
        proc_rates%cldice_number_tend_accretion_by_snow_per_kg_per_s(i,k)           = 0._r8
        proc_rates%rain_number_tend_accretion_by_snow_per_kg_per_s(i,k)          = 0._r8
        proc_rates%rain_number_tend_autoconversion_from_cldliq_per_kg_per_s(i,k)            = 0._r8
        proc_rates%rain_number_tend_self_collection_per_kg_per_s(i,k)           = 0._r8
        proc_rates%cldice_number_tend_autoconversion_to_snow_per_kg_per_s(i,k)           = 0._r8
        proc_rates%cldliq_number_tend_sedimentation_per_kg_per_s(i,k)           = 0._r8
        proc_rates%cldice_number_tend_sedimentation_per_kg_per_s(i,k)           = 0._r8
        proc_rates%rain_number_tend_sedimentation_per_kg_per_s(i,k)           = 0._r8
        proc_rates%snow_number_tend_sedimentation_per_kg_per_s(i,k)           = 0._r8
        proc_rates%graupel_number_tend_sedimentation_per_kg_per_s(i,k)           = 0._r8
        proc_rates%cldice_number_tend_melting_per_kg_per_s(i,k)           = 0._r8
        proc_rates%snow_number_tend_melting_per_kg_per_s(i,k)          = 0._r8
        proc_rates%graupel_number_tend_melting_per_kg_per_s(i,k)          = 0._r8

!need to zero these out to be totally switchable (for conservation)
        local_snow_mass_tend_collection_by_rain_kg_per_kg_per_s(i,k)              = 0._r8
        local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(i,k)              = 0._r8
        local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(i,k)             = 0._r8
        local_graupel_mass_tend_from_snow_riming_kg_per_kg_per_s(i,k)             = 0._r8
        local_graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s(i,k)             = 0._r8
        local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k)               = 0._r8
        local_cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s(i,k)             = 0._r8
        local_cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s(i,k)            = 0._r8
        local_rain_number_tend_accretion_by_graupel_per_kg_per_s(i,k)             = 0._r8
        local_graupel_number_tend_from_snow_riming_per_kg_per_s(i,k)              = 0._r8
        local_graupel_number_tend_rain_snow_collection_per_kg_per_s(i,k)             = 0._r8
        local_cldice_number_tend_splintering_by_graupel_per_kg_per_s(i,k)             = 0._r8
        local_cldice_number_tend_splintering_rain_graupel_per_kg_per_s(i,k)            = 0._r8
        local_cldliq_number_tend_riming_by_graupel_per_kg_per_s(i,k)            = 0._r8
        local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k)                = 0._r8
        local_rain_number_tend_autoconversion_from_cldliq_per_kg_per_s(i,k)               = 0._r8
        local_cldliq_number_tend_autoconversion_to_rain_per_kg_per_s(i,k)              = 0._r8
        local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)                = 0._r8
        local_cldliq_number_tend_accretion_by_rain_per_kg_per_s(i,k)               = 0._r8
        cldice_size_dist_intercept_per_m4(i,k)                = 0._r8
        cldice_size_dist_slope_per_m(i,k)               = 0._r8
        snow_size_dist_intercept_per_m4(i,k)                = 0._r8
        snow_size_dist_slope_per_m(i,k)               = 0._r8
        snow_number_change_instant_melting_per_kg(i,k)            = 0._r8
        graupel_number_change_instant_melting_per_kg(i,k)            = 0._r8
        local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k)             = 0._r8
        local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)              = 0._r8
        local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k)              = 0._r8
        local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)               = 0._r8
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev+1
     do i=1,mgncol
        rflx(i,k)               = 0._r8
        sflx(i,k)               = 0._r8
        lflx(i,k)               = 0._r8
        iflx(i,k)               = 0._r8
        gflx(i,k)               = 0._r8
        zint(i,k)               = 0._r8
     end do
  end do
  !$acc end parallel

  ! initialize precip at surface

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,mgncol
     prect(i)                   = 0._r8
     preci(i)                   = 0._r8
     prect_i(i)                 = 0._r8
     preci_i(i)                 = 0._r8
     prect_l(i)                 = 0._r8
     prect_r(i)                 = 0._r8
     prect_s(i)                 = 0._r8
     preci_s(i)                 = 0._r8
     prect_g(i)                 = 0._r8
     preci_g(i)                 = 0._r8
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        ! initialize precip output
        qrout(i,k)              = 0._r8
        qsout(i,k)              = 0._r8
        nrout(i,k)              = 0._r8
        nsout(i,k)              = 0._r8
        qgout(i,k)              = 0._r8
        ngout(i,k)              = 0._r8

        ! initialize rain size
        rercld(i,k)             = 0._r8

        qcsinksum_rate1ord(i,k) = 0._r8

        ! initialize variables for trop_mozart
        nevapr(i,k)             = 0._r8
        prer_evap(i,k)          = 0._r8
        proc_rates%snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)           = 0._r8
        am_evp_st(i,k)          = 0._r8
        prain(i,k)              = 0._r8
        proc_rates%snow_mass_tend_production_kg_per_kg_per_s(i,k)           = 0._r8
        cmeout(i,k)             = 0._r8

        precipitation_area_fraction(i,k)        = mincld
        cldliq_size_dist_slope_per_m(i,k)               = 0._r8
        graupel_size_dist_slope_per_m(i,k)               = 0._r8

        ! Interim variables for accretion
        rtmp(i,k)               = 0._r8
        ctmp(i,k)               = 0._r8
        ntmp(i,k)               = 0._r8

        ! initialize microphysical tendencies
        tlat(i,k)               = 0._r8
        qvlat(i,k)              = 0._r8
        qctend(i,k)             = 0._r8
        qitend(i,k)             = 0._r8
        qstend(i,k)             = 0._r8
        qrtend(i,k)             = 0._r8
        nctend(i,k)             = 0._r8
        nitend(i,k)             = 0._r8
        nrtend(i,k)             = 0._r8
        nstend(i,k)             = 0._r8
        qgtend(i,k)             = 0._r8
        ngtend(i,k)             = 0._r8

        ! initialize in-cloud and in-precip quantities to zero
        incloud_cldliq_mass_mixing_ratio_kg_per_kg(i,k)               = 0._r8
        incloud_cldice_mass_mixing_ratio_kg_per_kg(i,k)               = 0._r8
        inprecip_snow_mass_mixing_ratio_kg_per_kg(i,k)               = 0._r8
        inprecip_rain_mass_mixing_ratio_kg_per_kg(i,k)               = 0._r8
        inprecip_graupel_mass_mixing_ratio_kg_per_kg(i,k)               = 0._r8

        incloud_cldliq_number_conc_per_kg(i,k)               = 0._r8
        incloud_cldice_number_conc_per_kg(i,k)               = 0._r8
        inprecip_snow_number_conc_per_kg(i,k)               = 0._r8
        inprecip_rain_number_conc_per_kg(i,k)               = 0._r8
        inprecip_graupel_number_conc_per_kg(i,k)               = 0._r8
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        ! initialize vapor_deposition
        local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k)            = 0._r8
        local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)           = 0._r8

        ! initialize precip fallspeeds to zero
        proc_rates%snow_mass_weighted_fall_speed_m_per_s(i,k)     = 0._r8
        snow_number_weighted_fall_speed_m_per_s(i,k)                = 0._r8
        proc_rates%rain_mass_weighted_fall_speed_m_per_s(i,k)     = 0._r8
        rain_number_weighted_fall_speed_m_per_s(i,k)                = 0._r8
        proc_rates%graupel_mass_weighted_fall_speed_m_per_s(i,k)     = 0._r8
        graupel_number_weighted_fall_speed_m_per_s(i,k)                = 0._r8

        ! initialize limiter for output
        qcrat(i,k)              = 1._r8

        ! Many outputs have to be initialized here at the top to work around
        ! ifort problems, even if they are always overwritten later.
        effc(i,k)               = 10._r8
        lamcrad(i,k)            = 0._r8
        pgamrad(i,k)            = 0._r8
        effc_fn(i,k)            = 10._r8
        effi(i,k)               = 25._r8
        effi(i,k)               = effi(i,k)*micro_mg_effi_factor
        sadice(i,k)             = 0._r8
        sadsnow(i,k)            = 0._r8
        deffi(i,k)              = 50._r8

        qrout2(i,k)             = 0._r8
        nrout2(i,k)             = 0._r8
        drout2(i,k)             = 0._r8
        qsout2(i,k)             = 0._r8
        nsout2(i,k)             = 0._r8
        dsout(i,k)              = 0._r8
        dsout2(i,k)             = 0._r8
        qgout2(i,k)             = 0._r8
        ngout2(i,k)             = 0._r8
        freqg(i,k)              = 0._r8
        freqr(i,k)              = 0._r8
        freqs(i,k)              = 0._r8

        reff_rain(i,k)          = 0._r8
        reff_snow(i,k)          = 0._r8
        reff_grau(i,k)          = 0._r8

        refl(i,k)               = -9999._r8
        arefl(i,k)              = 0._r8
        areflz(i,k)             = 0._r8
        frefl(i,k)              = 0._r8
        csrfl(i,k)              = 0._r8
        acsrfl(i,k)             = 0._r8
        fcsrfl(i,k)             = 0._r8

        refl10cm(i,k)           = -9999._r8
        reflz10cm(i,k)          =  0._r8

        ncal(i,k)               = 0._r8
        ncai(i,k)               = 0._r8
        nfice(i,k)              = 0._r8

        pdel_inv(i,k)           = 1._r8/pdel(i,k)
        tlat_i(i,k)             = 0._r8
        qvlat_i(i,k)            = 0._r8
        tlat_l(i,k)             = 0._r8
        qvlat_l(i,k)            = 0._r8

        local_cldliq_number_tend_deposition_nucleation_per_kg_per_s(i,k)             = 0._r8
        local_cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s(i,k)             = 0._r8
        local_rain_number_tend_self_collection_per_kg_per_s(i,k)              = 0._r8

        proc_rates%cldliq_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s(i,k) = 0._r8
        proc_rates%cldliq_number_tend_autoconv_accretion_KK2000_per_kg_per_s(i,k) = 0._r8
        proc_rates%rain_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s(i,k) = 0._r8
        proc_rates%rain_number_tend_autoconv_accretion_KK2000_per_kg_per_s(i,k) = 0._r8
        proc_rates%cldliq_size_dist_slope_output_per_m(i,k)      = 0._r8
        proc_rates%rain_size_dist_slope_output_per_m(i,k)      = 0._r8
        proc_rates%cldliq_size_dist_shape_output(i,k)      = 0._r8
        proc_rates%rain_size_dist_intercept_output_per_m4(i,k)       = 0._r8
     end do
  end do
  !$acc end parallel

  if (trim(warm_rain) == 'sb2001') then
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           proc_rates%cldliq_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s(i,k) = 0._r8
           proc_rates%cldliq_number_tend_autoconv_accretion_SB2001_per_kg_per_s(i,k) = 0._r8
           proc_rates%rain_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s(i,k) = 0._r8
           proc_rates%rain_number_tend_autoconv_accretion_SB2001_per_kg_per_s(i,k) = 0._r8
        end do
     end do
     !$acc end parallel
  end if

  if (trim(warm_rain) == 'tau' .or. trim(warm_rain) == 'emulated') then
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           proc_rates%cldliq_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s(i,k) = 0._r8
           proc_rates%cldliq_number_tend_autoconv_accretion_TAU_per_kg_per_s(i,k) = 0._r8
           proc_rates%rain_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s(i,k) = 0._r8
           proc_rates%rain_number_tend_autoconv_accretion_TAU_per_kg_per_s(i,k) = 0._r8
           proc_rates%qc_out_TAU(i,k) = 0._r8
           proc_rates%nc_out_TAU(i,k) = 0._r8
           proc_rates%qr_out_TAU(i,k) = 0._r8
           proc_rates%nr_out_TAU(i,k) = 0._r8
           proc_rates%qc_in_TAU(i,k) = 0._r8
           proc_rates%nc_in_TAU(i,k) = 0._r8
           proc_rates%qr_in_TAU(i,k) = 0._r8
           proc_rates%nr_in_TAU(i,k) = 0._r8
           proc_rates%gmnnn_lmnnn_TAU(i,k) = 0._r8
        end do
     end do
     !$acc end parallel
  end if

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! droplet activation
  ! get provisional droplet number after activation. This is used for
  ! all microphysical process calculations, for consistency with update of
  ! droplet mass before microphysics

  ! calculate potential for droplet activation if cloud water is present
  ! tendency from activation (npccn) is read in from companion routine

  ! output activated liquid and ice (convert from #/kg -> #/m3)
  !--------------------------------------------------

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        if (cldliq_mass_mixing_ratio_kg_per_kg(i,k) >= qsmall) then
           cldliq_number_conc_per_kg(i,k) = max(cldliq_number_conc_per_kg(i,k) + npccn(i,k)*deltat, 0._r8)
           ncal(i,k) = npccn(i,k)
        else
           ncal(i,k) = 0._r8
        end if

        if (t(i,k) < icenuct) then
           ncai(i,k) = naai(i,k)*deltat*air_density_kg_per_m3(i,k)
        else
           ncai(i,k) = 0._r8
        end if

  !===============================================

  ! ice nucleation if activated nuclei exist at t<-5C AND rhmini + 5%
  !
  ! NOTE: If using gridbox average values, condensation will not occur until rh=1,
  ! so the threshold seems like it should be 1.05 and not rhmini + 0.05. For subgrid
  ! clouds (using rhmini and qsfacm), the relhum has already been adjusted, and thus
  ! the nucleation threshold should also be 1.05 and not rhmini + 0.05.
  !-------------------------------------------------------

        if (do_cldice) then
           if (icenuc_rh_off) then
              if (naai(i,k) > 0._r8 .and. t(i,k) < icenuct) then
                 !if NAAI > 0. then set numice = naai (as before)
                 !note: this is gridbox averaged
                 local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k) = naai(i,k)*ice_cloud_area_fraction(i,k)
                 local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k) = max(local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k),0._r8)

                 !Calc mass of new particles using new crystal mass...
                 !also this will be multiplied by mtime as nnuccd is...
                 local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k) = local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k) * mi0
              else
                 local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k) = 0._r8
                 local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k) = 0._r8
              end if
           else
              if (naai(i,k) > 0._r8 .and. t(i,k) < icenuct .and. &
                 relative_humidity(i,k)*saturation_vapor_pressure_over_liquid_Pa(i,k)/saturation_vapor_pressure_over_ice_Pa(i,k) > 1.05_r8) then
                 !if NAAI > 0. then set numice = naai (as before)
                 !note: this is gridbox averaged
                 local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k) = naai(i,k)*ice_cloud_area_fraction(i,k)
                 local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k) = max(local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k),0._r8)

                 !Calc mass of new particles using new crystal mass...
                 !also this will be multiplied by mtime as nnuccd is...
                 local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k) = local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k) * mi0
              else
                 local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k) = 0._r8
                 local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k) = 0._r8
              end if
           end if
        end if
     end do
  end do
  !$acc end parallel

  !=============================================================================

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        ! calculate instantaneous precip processes (melting and homogeneous freezing)
        ! melting of snow at +2 C
        if (t(i,k) > snowmelt) then
           if (snow_mass_mixing_ratio_kg_per_kg(i,k) > 0._r8) then
              ! make sure melting snow doesn't reduce temperature below threshold
              dum = -xlf/cpp*snow_mass_mixing_ratio_kg_per_kg(i,k)
              if (t(i,k)+dum < snowmelt) then
                 dum = (t(i,k)-snowmelt)*cpp/xlf
                 dum = dum/snow_mass_mixing_ratio_kg_per_kg(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              snow_mass_change_instant_melting_kg_per_kg(i,k) = dum*snow_mass_mixing_ratio_kg_per_kg(i,k)
              snow_number_change_instant_melting_per_kg(i,k) = dum*snow_number_conc_per_kg(i,k)

              dum1=-xlf*snow_mass_change_instant_melting_kg_per_kg(i,k)*rdeltat
              tlat(i,k)=tlat(i,k)+dum1
              proc_rates%air_heating_rate_snow_melting_W_per_kg(i,k)=proc_rates%air_heating_rate_snow_melting_W_per_kg(i,k) + dum1
              proc_rates%snow_mass_tend_melting_kg_per_kg_per_s(i,k)=snow_mass_change_instant_melting_kg_per_kg(i,k)*rdeltat

              snow_mass_mixing_ratio_kg_per_kg(i,k) = max(snow_mass_mixing_ratio_kg_per_kg(i,k) - snow_mass_change_instant_melting_kg_per_kg(i,k), 0._r8)
              snow_number_conc_per_kg(i,k) = max(snow_number_conc_per_kg(i,k) - snow_number_change_instant_melting_per_kg(i,k), 0._r8)
              rain_mass_mixing_ratio_kg_per_kg(i,k) = max(rain_mass_mixing_ratio_kg_per_kg(i,k) + snow_mass_change_instant_melting_kg_per_kg(i,k), 0._r8)
              rain_number_conc_per_kg(i,k) = max(rain_number_conc_per_kg(i,k) + snow_number_change_instant_melting_per_kg(i,k), 0._r8)
           end if
        end if

        ! melting of graupel at +2 C

        if (t(i,k) > snowmelt) then
           if (graupel_mass_mixing_ratio_kg_per_kg(i,k) > 0._r8) then

              ! make sure melting graupel doesn't reduce temperature below threshold
              dum = -xlf/cpp*graupel_mass_mixing_ratio_kg_per_kg(i,k)
              if (t(i,k)+dum < snowmelt) then
                 dum = (t(i,k)-snowmelt)*cpp/xlf
                 dum = dum/graupel_mass_mixing_ratio_kg_per_kg(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              graupel_mass_change_instant_melting_kg_per_kg(i,k) = dum*graupel_mass_mixing_ratio_kg_per_kg(i,k)
              graupel_number_change_instant_melting_per_kg(i,k) = dum*graupel_number_conc_per_kg(i,k)

              dum1=-xlf*graupel_mass_change_instant_melting_kg_per_kg(i,k)*rdeltat
              tlat(i,k)=tlat(i,k)+dum1
              proc_rates%air_heating_rate_snow_melting_W_per_kg(i,k)=proc_rates%air_heating_rate_snow_melting_W_per_kg(i,k) + dum1
              proc_rates%graupel_mass_tend_melting_kg_per_kg_per_s(i,k)=graupel_mass_change_instant_melting_kg_per_kg(i,k)*rdeltat

              graupel_mass_mixing_ratio_kg_per_kg(i,k) = max(graupel_mass_mixing_ratio_kg_per_kg(i,k) - graupel_mass_change_instant_melting_kg_per_kg(i,k), 0._r8)
              graupel_number_conc_per_kg(i,k) = max(graupel_number_conc_per_kg(i,k) - graupel_number_change_instant_melting_per_kg(i,k), 0._r8)
              rain_mass_mixing_ratio_kg_per_kg(i,k) = max(rain_mass_mixing_ratio_kg_per_kg(i,k) + graupel_mass_change_instant_melting_kg_per_kg(i,k), 0._r8)
              rain_number_conc_per_kg(i,k) = max(rain_number_conc_per_kg(i,k) + graupel_number_change_instant_melting_per_kg(i,k), 0._r8)
           end if
        end if

        ! freezing of rain at -5 C

        if (t(i,k) < rainfrze) then

           if (rain_mass_mixing_ratio_kg_per_kg(i,k) > 0._r8) then

              ! make sure freezing rain doesn't increase temperature above threshold
              dum = xlf/cpp*rain_mass_mixing_ratio_kg_per_kg(i,k)
              if (t(i,k)+dum > rainfrze) then
                 dum = -(t(i,k)-rainfrze)*cpp/xlf
                 dum = dum/rain_mass_mixing_ratio_kg_per_kg(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              rain_mass_change_instant_freezing_kg_per_kg(i,k) = dum*rain_mass_mixing_ratio_kg_per_kg(i,k)
              rain_number_change_instant_freezing_per_kg(i,k) = dum*rain_number_conc_per_kg(i,k)

              ! heating tendency
              dum1 = xlf*rain_mass_change_instant_freezing_kg_per_kg(i,k)*rdeltat
              tlat(i,k)=tlat(i,k)+dum1
              proc_rates%air_heating_rate_rain_freezing_W_per_kg(i,k)=proc_rates%air_heating_rate_rain_freezing_W_per_kg(i,k) + dum1

              rain_mass_mixing_ratio_kg_per_kg(i,k) = max(rain_mass_mixing_ratio_kg_per_kg(i,k) - rain_mass_change_instant_freezing_kg_per_kg(i,k), 0._r8)
              rain_number_conc_per_kg(i,k) = max(rain_number_conc_per_kg(i,k) - rain_number_change_instant_freezing_per_kg(i,k), 0._r8)

              ! freeze rain to graupel not snow.
              if(do_hail.or.do_graupel) then
                 graupel_mass_mixing_ratio_kg_per_kg(i,k) = max(graupel_mass_mixing_ratio_kg_per_kg(i,k) + rain_mass_change_instant_freezing_kg_per_kg(i,k), 0._r8)
                 graupel_number_conc_per_kg(i,k) = max(graupel_number_conc_per_kg(i,k) + rain_number_change_instant_freezing_per_kg(i,k), 0._r8)
              else
                 snow_mass_mixing_ratio_kg_per_kg(i,k) = max(snow_mass_mixing_ratio_kg_per_kg(i,k) + rain_mass_change_instant_freezing_kg_per_kg(i,k), 0._r8)
                 snow_number_conc_per_kg(i,k) = max(snow_number_conc_per_kg(i,k) + rain_number_change_instant_freezing_per_kg(i,k), 0._r8)
              end if
           end if
        end if
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
    do i=1,mgncol
        ! obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
        !-------------------------------------------------------
        ! for microphysical process calculations
        ! units are kg/kg for mixing ratio, 1/kg for number conc

        if (cldliq_mass_mixing_ratio_kg_per_kg(i,k).ge.qsmall) then
           ! limit in-cloud values to 0.005 kg/kg
           incloud_cldliq_mass_mixing_ratio_kg_per_kg(i,k)=min(cldliq_mass_mixing_ratio_kg_per_kg(i,k)/liquid_cloud_area_fraction(i,k),5.e-3_r8)
           incloud_cldliq_number_conc_per_kg(i,k)=max(cldliq_number_conc_per_kg(i,k)/liquid_cloud_area_fraction(i,k),0._r8)

           ! specify droplet concentration
           if (nccons) then
              incloud_cldliq_number_conc_per_kg(i,k)=ncnst/air_density_kg_per_m3(i,k)
           end if
        else
           incloud_cldliq_mass_mixing_ratio_kg_per_kg(i,k)=0._r8
           incloud_cldliq_number_conc_per_kg(i,k)=0._r8
        end if

        if (cldice_mass_mixing_ratio_kg_per_kg(i,k).ge.qsmall) then
           ! limit in-cloud values to 0.005 kg/kg
           incloud_cldice_mass_mixing_ratio_kg_per_kg(i,k)=min(cldice_mass_mixing_ratio_kg_per_kg(i,k)/ice_cloud_area_fraction(i,k),5.e-3_r8)
           incloud_cldice_number_conc_per_kg(i,k)=max(cldice_number_conc_per_kg(i,k)/ice_cloud_area_fraction(i,k),0._r8)

           ! switch for specification of cloud ice number
           if (nicons) then
              incloud_cldice_number_conc_per_kg(i,k)=ninst/air_density_kg_per_m3(i,k)
           end if
        else
           incloud_cldice_mass_mixing_ratio_kg_per_kg(i,k)=0._r8
           incloud_cldice_number_conc_per_kg(i,k)=0._r8
        end if

  !========================================================================

  ! for sub-columns cldm has already been set to 1 if cloud
  ! water or ice is present, so precip_frac will be correctly set below
  ! and nothing extra needs to be done here

        precipitation_area_fraction(i,k) = cloud_area_fraction(i,k)
     end do
  end do
  !$acc end parallel

  if (precip_frac_method == MG_PRECIP_FRAC_INCLOUD) then
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector
     do i=1,mgncol
        !$acc loop seq
        do k=2,nlev
           if (cldliq_mass_mixing_ratio_kg_per_kg(i,k) < qsmall .and. cldice_mass_mixing_ratio_kg_per_kg(i,k) < qsmall) then
              precipitation_area_fraction(i,k) = precipitation_area_fraction(i,k-1)
           end if
        end do
     end do
     !$acc end parallel
  else if (precip_frac_method == MG_PRECIP_FRAC_OVERLAP) then
     ! calculate precip fraction based on maximum overlap assumption

     ! if rain or snow mix ratios are smaller than threshold,
     ! then leave precip_frac as cloud fraction at current level

     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector
     do i=1,mgncol
        !$acc loop seq
        do k=2,nlev
           if (rain_mass_mixing_ratio_kg_per_kg(i,k-1) >= qsmall .or. snow_mass_mixing_ratio_kg_per_kg(i,k-1) >= qsmall .or. graupel_mass_mixing_ratio_kg_per_kg(i,k-1) >= qsmall) then
              precipitation_area_fraction(i,k)=max(precipitation_area_fraction(i,k-1),precipitation_area_fraction(i,k))
           end if
        end do
     end do
     !$acc end parallel
  end if

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! get size distribution parameters based on in-cloud cloud water
  ! these calculations also ensure consistency between number and mixing ratio
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! cloud liquid
  !-------------------------------------------
  call size_dist_param_liq(mg_liq_props, incloud_cldliq_mass_mixing_ratio_kg_per_kg, incloud_cldliq_number_conc_per_kg, air_density_kg_per_m3, cldliq_size_dist_shape, cldliq_size_dist_slope_per_m, mgncol, nlev)


  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        ! assign qric based on prognostic qr, using assumed precip fraction
        ! note: this could be moved above for consistency with qcic and qiic calculations
        inprecip_rain_mass_mixing_ratio_kg_per_kg(i,k) = rain_mass_mixing_ratio_kg_per_kg(i,k)/precipitation_area_fraction(i,k)
        inprecip_rain_number_conc_per_kg(i,k) = rain_number_conc_per_kg(i,k)/precipitation_area_fraction(i,k)

        ! limit in-precip mixing ratios to 10 g/kg
        inprecip_rain_mass_mixing_ratio_kg_per_kg(i,k)=min(inprecip_rain_mass_mixing_ratio_kg_per_kg(i,k),0.01_r8)

        ! add autoconversion to precip from above to get provisional rain mixing ratio
        ! and number concentration (qric and nric)

        if (inprecip_rain_mass_mixing_ratio_kg_per_kg(i,k).lt.qsmall) then
           inprecip_rain_mass_mixing_ratio_kg_per_kg(i,k)=0._r8
           inprecip_rain_number_conc_per_kg(i,k)=0._r8
        end if

        ! make sure number concentration is a positive number to avoid
        ! taking root of negative later

        inprecip_rain_number_conc_per_kg(i,k)=max(inprecip_rain_number_conc_per_kg(i,k),0._r8)
     end do
  end do
  !$acc end parallel

  ! get size distribution parameters for rain
  !......................................................................

  call size_dist_param_basic(mg_rain_props, inprecip_rain_mass_mixing_ratio_kg_per_kg, inprecip_rain_number_conc_per_kg, rain_size_dist_slope_per_m, mgncol, nlev, n0=rain_size_dist_intercept_per_m4)

  ! Save off size distribution parameters for output
  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        proc_rates%cldliq_size_dist_shape_output(i,k)=cldliq_size_dist_shape(i,k)
        proc_rates%rain_size_dist_intercept_output_per_m4(i,k)=rain_size_dist_intercept_per_m4(i,k)
        proc_rates%cldliq_size_dist_slope_output_per_m(i,k)=cldliq_size_dist_slope_per_m(i,k)
        proc_rates%rain_size_dist_slope_output_per_m(i,k)=rain_size_dist_slope_per_m(i,k)
     end do
  end do
  !$acc end parallel

  !========================================================================
  ! autoconversion of cloud liquid water to rain
  ! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
  ! minimum qc of 1 x 10^-8 prevents floating point error

  call kk2000_liq_autoconversion(microp_uniform, incloud_cldliq_mass_mixing_ratio_kg_per_kg, incloud_cldliq_number_conc_per_kg, air_density_kg_per_m3, relvar, &
                      proc_rates%cldliq_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s, proc_rates%rain_number_tend_autoconv_accretion_KK2000_per_kg_per_s, &
                      proc_rates%cldliq_number_tend_autoconv_accretion_KK2000_per_kg_per_s, micro_mg_autocon_fact, &
                      micro_mg_autocon_nd_exp, micro_mg_autocon_lwp_exp, &
                      mgncol*nlev)

  ! Write to pumas tendency arrays if kk2000 is active, otherwise just record diagnostics

  if ( trim(warm_rain) == 'kk2000' ) then
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k)=proc_rates%cldliq_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s(i,k)
           local_cldliq_number_tend_autoconversion_to_rain_per_kg_per_s(i,k)=proc_rates%cldliq_number_tend_autoconv_accretion_KK2000_per_kg_per_s(i,k)
           local_rain_number_tend_autoconversion_from_cldliq_per_kg_per_s(i,k)=proc_rates%rain_number_tend_autoconv_accretion_KK2000_per_kg_per_s(i,k)
           proc_rates%rain_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s(i,k)=-proc_rates%cldliq_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s(i,k)
        end do
     end do
     !$acc end parallel
  end if

  if ( trim(warm_rain) == 'tau' ) then
     call pumas_stochastic_collect_tau_tend(deltatin,t,air_density_kg_per_m3,qcn,qrn,incloud_cldliq_mass_mixing_ratio_kg_per_kg,incloud_cldliq_number_conc_per_kg,inprecip_rain_mass_mixing_ratio_kg_per_kg, &
                                      inprecip_rain_number_conc_per_kg,liquid_cloud_area_fraction,precipitation_area_fraction,cldliq_size_dist_shape,cldliq_size_dist_slope_per_m,rain_size_dist_intercept_per_m4,rain_size_dist_slope_per_m,   &
                                      proc_rates%qc_out_TAU,proc_rates%nc_out_TAU, &
                                      proc_rates%qr_out_TAU,proc_rates%nr_out_TAU, &
                                      qctend,nctend,qrtend,nrtend, &
                                      proc_rates%cldliq_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s,proc_rates%cldliq_number_tend_autoconv_accretion_TAU_per_kg_per_s, &
                                      proc_rates%rain_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s,proc_rates%rain_number_tend_autoconv_accretion_TAU_per_kg_per_s, &
                                      proc_rates%scale_qc,proc_rates%scale_nc, &
                                      proc_rates%scale_qr,proc_rates%scale_nr, &
                                      proc_rates%amk_c,proc_rates%ank_c, &
                                      proc_rates%amk_r,proc_rates%ank_r, &
                                      proc_rates%amk, proc_rates%ank, &
                                      proc_rates%amk_out, proc_rates%ank_out, &
                                      proc_rates%gmnnn_lmnnn_TAU,mgncol,nlev)

     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           proc_rates%qc_in_TAU(i,k)=incloud_cldliq_mass_mixing_ratio_kg_per_kg(i,k)
           proc_rates%nc_in_TAU(i,k)=incloud_cldliq_number_conc_per_kg(i,k)
           proc_rates%qr_in_TAU(i,k)=inprecip_rain_mass_mixing_ratio_kg_per_kg(i,k)
           proc_rates%nr_in_TAU(i,k)=inprecip_rain_number_conc_per_kg(i,k)
           ! PUMAS expects prc and nprc1 (cloud rates) are positive
           local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k)= -proc_rates%cldliq_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s(i,k)
           local_cldliq_number_tend_autoconversion_to_rain_per_kg_per_s(i,k)= -proc_rates%cldliq_number_tend_autoconv_accretion_TAU_per_kg_per_s(i,k)

           ! PUMAS expects nprc to be positive. Negative nrtend_TAU is from self collection, so put it into nragg
           if ( proc_rates%rain_number_tend_autoconv_accretion_TAU_per_kg_per_s(i,k) > 0._r8 ) then
              local_rain_number_tend_autoconversion_from_cldliq_per_kg_per_s(i,k)= proc_rates%rain_number_tend_autoconv_accretion_TAU_per_kg_per_s(i,k)
           else
              local_rain_number_tend_self_collection_per_kg_per_s(i,k) = proc_rates%rain_number_tend_autoconv_accretion_TAU_per_kg_per_s(i,k)
           end if
        end do
     end do
     !$acc end parallel

  else if (trim(warm_rain) == 'emulated') then
     ! JS - 08/22/2023: this code block only works on CPU

     !$acc update self(incloud_cldliq_mass_mixing_ratio_kg_per_kg,incloud_cldliq_number_conc_per_kg,inprecip_rain_mass_mixing_ratio_kg_per_kg,inprecip_rain_number_conc_per_kg,air_density_kg_per_m3,liquid_cloud_area_fraction,precipitation_area_fraction, &
     !$acc             proc_rates%cldliq_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s,proc_rates%rain_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s, &
     !$acc             proc_rates%cldliq_number_tend_autoconv_accretion_TAU_per_kg_per_s,proc_rates%rain_number_tend_autoconv_accretion_TAU_per_kg_per_s, &
     !$acc             cldliq_mass_mixing_ratio_kg_per_kg,cldliq_number_conc_per_kg,rain_mass_mixing_ratio_kg_per_kg,rain_number_conc_per_kg,local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s,local_cldliq_number_tend_autoconversion_to_rain_per_kg_per_s,local_rain_number_tend_autoconversion_from_cldliq_per_kg_per_s,local_rain_number_tend_self_collection_per_kg_per_s)

     do k=1,nlev
        call tau_emulated_cloud_rain_interactions(incloud_cldliq_mass_mixing_ratio_kg_per_kg(1:mgncol,k), incloud_cldliq_number_conc_per_kg(1:mgncol,k), &
                                                  inprecip_rain_mass_mixing_ratio_kg_per_kg(1:mgncol,k), inprecip_rain_number_conc_per_kg(1:mgncol,k), &
                                                  air_density_kg_per_m3(1:mgncol,k), liquid_cloud_area_fraction(1:mgncol,k), &
                                                  precipitation_area_fraction(1:mgncol,k), mgncol, qsmall, &
                                                  proc_rates%cldliq_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s(1:mgncol,k), &
                                                  proc_rates%rain_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s(1:mgncol,k), &
                                                  proc_rates%cldliq_number_tend_autoconv_accretion_TAU_per_kg_per_s(1:mgncol,k), &
                                                  proc_rates%rain_number_tend_autoconv_accretion_TAU_per_kg_per_s(1:mgncol,k))

        call ML_fixer_calc(mgncol, deltatin, cldliq_mass_mixing_ratio_kg_per_kg(1:mgncol,k), cldliq_number_conc_per_kg(1:mgncol,k), &
                           rain_mass_mixing_ratio_kg_per_kg(1:mgncol,k), rain_number_conc_per_kg(1:mgncol,k), &
                           proc_rates%cldliq_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s(1:mgncol,k),&
                           proc_rates%cldliq_number_tend_autoconv_accretion_TAU_per_kg_per_s(1:mgncol,k), &
                           proc_rates%rain_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s(1:mgncol,k), &
                           proc_rates%rain_number_tend_autoconv_accretion_TAU_per_kg_per_s(1:mgncol,k), &
                           proc_rates%ML_fixer(1:mgncol,k), &
                           proc_rates%QC_fixer(1:mgncol,k), &
                           proc_rates%NC_fixer(1:mgncol,k), &
                           proc_rates%QR_fixer(1:mgncol,k), &
                           proc_rates%NR_fixer(1:mgncol,k))

        ! PUMAS expects prc and nprc1 (cloud rates) are positive
        local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(1:mgncol,k)= -proc_rates%cldliq_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s(1:mgncol,k)
        local_cldliq_number_tend_autoconversion_to_rain_per_kg_per_s(1:mgncol,k)= -proc_rates%cldliq_number_tend_autoconv_accretion_TAU_per_kg_per_s(1:mgncol,k)

        ! PUMAS expects nprc to be positive. Negative nrtend_TAU is from self
        ! collection, so put it into nragg
        do i=1,mgncol
           if (proc_rates%rain_number_tend_autoconv_accretion_TAU_per_kg_per_s(i,k).gt.0._r8) then
              local_rain_number_tend_autoconversion_from_cldliq_per_kg_per_s(i,k)= proc_rates%rain_number_tend_autoconv_accretion_TAU_per_kg_per_s(i,k)
           else
              local_rain_number_tend_self_collection_per_kg_per_s(i,k)= proc_rates%rain_number_tend_autoconv_accretion_TAU_per_kg_per_s(i,k)
           end if
        end do

     end do

     !$acc update device(proc_rates%cldliq_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s,proc_rates%rain_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s, &
     !$acc               proc_rates%cldliq_number_tend_autoconv_accretion_TAU_per_kg_per_s,proc_rates%rain_number_tend_autoconv_accretion_TAU_per_kg_per_s, &
     !$acc               local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s,local_cldliq_number_tend_autoconversion_to_rain_per_kg_per_s,local_rain_number_tend_autoconversion_from_cldliq_per_kg_per_s,local_rain_number_tend_self_collection_per_kg_per_s)

  end if

  ! Alternative autoconversion
  if (trim(warm_rain) == 'sb2001') then
     call sb2001v2_liq_autoconversion(cldliq_size_dist_shape,incloud_cldliq_mass_mixing_ratio_kg_per_kg,incloud_cldliq_number_conc_per_kg,inprecip_rain_mass_mixing_ratio_kg_per_kg,air_density_kg_per_m3,relvar, &
                                      proc_rates%cldliq_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s, &
                                      proc_rates%rain_number_tend_autoconv_accretion_SB2001_per_kg_per_s, &
                                      proc_rates%cldliq_number_tend_autoconv_accretion_SB2001_per_kg_per_s, &
                                      mgncol*nlev)

     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k)=proc_rates%cldliq_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s(i,k)
           local_rain_number_tend_autoconversion_from_cldliq_per_kg_per_s(i,k)=proc_rates%rain_number_tend_autoconv_accretion_SB2001_per_kg_per_s(i,k)
           local_cldliq_number_tend_autoconversion_to_rain_per_kg_per_s(i,k)=proc_rates%cldliq_number_tend_autoconv_accretion_SB2001_per_kg_per_s(i,k)
           proc_rates%rain_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s(i,k)= -proc_rates%cldliq_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s(i,k)
        end do
     end do
     !$acc end parallel
  end if

  ! Get size distribution parameters for cloud ice
  call size_dist_param_basic(mg_ice_props, incloud_cldice_mass_mixing_ratio_kg_per_kg, incloud_cldice_number_conc_per_kg, cldice_size_dist_slope_per_m, mgncol, nlev, n0=cldice_size_dist_intercept_per_m4)

  !.......................................................................
  ! Autoconversion of cloud ice to snow
  ! similar to Ferrier (1994)
  if (do_cldice) then
     call ice_autoconversion(t, incloud_cldice_mass_mixing_ratio_kg_per_kg, cldice_size_dist_slope_per_m, cldice_size_dist_intercept_per_m4, dcs, local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s, local_cldice_number_tend_autoconversion_to_snow_per_kg_per_s, mgncol*nlev)
  else
     ! Add in the particles that we have already converted to snow, and
     ! don't do any further autoconversion of ice.

     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(i,k)  = tnd_qsnow(i,k) / cloud_area_fraction(i,k)
           local_cldice_number_tend_autoconversion_to_snow_per_kg_per_s(i,k) = tnd_nsnow(i,k) / cloud_area_fraction(i,k)
        end do
     end do
     !$acc end parallel
  end if

  ! note, currently we don't have this
  ! inside the do_cldice block, should be changed later
  ! assign qsic based on prognostic qs, using assumed precip fraction

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        inprecip_snow_mass_mixing_ratio_kg_per_kg(i,k) = snow_mass_mixing_ratio_kg_per_kg(i,k)/precipitation_area_fraction(i,k)
        inprecip_snow_number_conc_per_kg(i,k) = snow_number_conc_per_kg(i,k)/precipitation_area_fraction(i,k)

        ! limit in-precip mixing ratios to 10 g/kg
        inprecip_snow_mass_mixing_ratio_kg_per_kg(i,k)=min(inprecip_snow_mass_mixing_ratio_kg_per_kg(i,k),0.01_r8)

        ! if precip mix ratio is zero so should number concentration
        if (inprecip_snow_mass_mixing_ratio_kg_per_kg(i,k) < qsmall) then
           inprecip_snow_mass_mixing_ratio_kg_per_kg(i,k)=0._r8
           inprecip_snow_number_conc_per_kg(i,k)=0._r8
        end if

        ! make sure number concentration is a positive number to avoid
        ! taking root of negative later
        inprecip_snow_number_conc_per_kg(i,k)=max(inprecip_snow_number_conc_per_kg(i,k),0._r8)

        ! also do this for graupel, which is assumed to be 'precip_frac'
        inprecip_graupel_mass_mixing_ratio_kg_per_kg(i,k) = graupel_mass_mixing_ratio_kg_per_kg(i,k)/precipitation_area_fraction(i,k)
        inprecip_graupel_number_conc_per_kg(i,k) = graupel_number_conc_per_kg(i,k)/precipitation_area_fraction(i,k)

        ! limit in-precip mixing ratios to 10 g/kg
        inprecip_graupel_mass_mixing_ratio_kg_per_kg(i,k)=min(inprecip_graupel_mass_mixing_ratio_kg_per_kg(i,k),0.01_r8)

        ! if precip mix ratio is zero so should number concentration
        if (inprecip_graupel_mass_mixing_ratio_kg_per_kg(i,k) < qsmall) then
           inprecip_graupel_mass_mixing_ratio_kg_per_kg(i,k)=0._r8
           inprecip_graupel_number_conc_per_kg(i,k)=0._r8
        end if

        ! make sure number concentration is a positive number to avoid
        ! taking root of negative later
        inprecip_graupel_number_conc_per_kg(i,k)=max(inprecip_graupel_number_conc_per_kg(i,k),0._r8)
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        if (rain_size_dist_slope_per_m(i,k) >= qsmall) then
           dum_2D(i,k)= rain_size_dist_slope_per_m(i,k)**br
           ! provisional rain number and mass weighted mean fallspeed (m/s)
           rain_number_weighted_fall_speed_m_per_s(i,k) = min(rain_fall_speed_prefactor(i,k)*gamma_br_plus1/dum_2D(i,k),9.1_r8*fall_speed_air_density_correction(i,k))
           proc_rates%rain_mass_weighted_fall_speed_m_per_s(i,k) = min(rain_fall_speed_prefactor(i,k)*gamma_br_plus4/(6._r8*dum_2D(i,k)),9.1_r8*fall_speed_air_density_correction(i,k))
        else
           proc_rates%rain_mass_weighted_fall_speed_m_per_s(i,k) = 0._r8
           rain_number_weighted_fall_speed_m_per_s(i,k) = 0._r8
        end if
     end do
  end do
  !$acc end parallel

  !......................................................................
  ! snow
  call size_dist_param_basic(mg_snow_props, inprecip_snow_mass_mixing_ratio_kg_per_kg, inprecip_snow_number_conc_per_kg, snow_size_dist_slope_per_m, mgncol, nlev, n0=snow_size_dist_intercept_per_m4)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        if (ifs_sed) then
           if (snow_size_dist_slope_per_m(i,k) > 0._r8) then
              proc_rates%snow_mass_weighted_fall_speed_m_per_s(i,k) = 1._r8
              snow_number_weighted_fall_speed_m_per_s(i,k) = 1._r8
           else
              proc_rates%snow_mass_weighted_fall_speed_m_per_s(i,k) = 0._r8
              snow_number_weighted_fall_speed_m_per_s(i,k) = 0._r8
           end if
        else
           if (snow_size_dist_slope_per_m(i,k) > 0._r8) then
              dum_2D(i,k) = snow_size_dist_slope_per_m(i,k)**bs
              ! provisional snow number and mass weighted mean fallspeed (m/s)
              proc_rates%snow_mass_weighted_fall_speed_m_per_s(i,k) = min(snow_fall_speed_prefactor(i,k)*gamma_bs_plus4/(6._r8*dum_2D(i,k)),1.2_r8*fall_speed_air_density_correction(i,k))
              proc_rates%snow_mass_weighted_fall_speed_m_per_s(i,k) = proc_rates%snow_mass_weighted_fall_speed_m_per_s(i,k)*micro_mg_vtrms_factor
              snow_number_weighted_fall_speed_m_per_s(i,k) = min(snow_fall_speed_prefactor(i,k)*gamma_bs_plus1/dum_2D(i,k),1.2_r8*fall_speed_air_density_correction(i,k))
           else
              proc_rates%snow_mass_weighted_fall_speed_m_per_s(i,k) = 0._r8
              snow_number_weighted_fall_speed_m_per_s(i,k) = 0._r8
           end if
        end if
     end do
  end do
  !$acc end parallel

  !  graupel/hail size distributions and properties

  if (do_hail) then
     call size_dist_param_basic(mg_hail_props, inprecip_graupel_mass_mixing_ratio_kg_per_kg, inprecip_graupel_number_conc_per_kg, graupel_size_dist_slope_per_m, mgncol, nlev, n0=graupel_size_dist_intercept_per_m4)
  end if
  if (do_graupel) then
     call size_dist_param_basic(mg_graupel_props, inprecip_graupel_mass_mixing_ratio_kg_per_kg, inprecip_graupel_number_conc_per_kg, graupel_size_dist_slope_per_m, mgncol, nlev, n0=graupel_size_dist_intercept_per_m4)
  end if

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        if (graupel_size_dist_slope_per_m(i,k) > 0._r8) then
           dum_2D(i,k) = graupel_size_dist_slope_per_m(i,k)**bgtmp
           ! provisional graupel/hail number and mass weighted mean fallspeed (m/s)
           proc_rates%graupel_mass_weighted_fall_speed_m_per_s(i,k) = min(graupel_fall_speed_prefactor(i,k)*gamma_bg_plus4/(6._r8*dum_2D(i,k)),20._r8*fall_speed_air_density_correction(i,k))
           graupel_number_weighted_fall_speed_m_per_s(i,k) = min(graupel_fall_speed_prefactor(i,k)*gamma_bg_plus1/dum_2D(i,k),20._r8*fall_speed_air_density_correction(i,k))
        else
           proc_rates%graupel_mass_weighted_fall_speed_m_per_s(i,k) = 0._r8
           graupel_number_weighted_fall_speed_m_per_s(i,k) = 0._r8
        end if
     end do
  end do
  !$acc end parallel

  if (do_cldice) then
     if (.not. use_hetfrz_classnuc) then
        ! heterogeneous freezing of cloud water via Bigg, 1953
        !----------------------------------------------
        call immersion_freezing(microp_uniform, t, cldliq_size_dist_shape, cldliq_size_dist_slope_per_m, incloud_cldliq_mass_mixing_ratio_kg_per_kg, incloud_cldliq_number_conc_per_kg, relvar, local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s, local_cldliq_number_tend_immersion_freezing_per_kg_per_s, mgncol*nlev)

        ! make sure number of droplets frozen does not exceed available ice nuclei concentration
        ! this prevents 'runaway' droplet freezing

        !$acc parallel vector_length(VLENS) default(present)
        !$acc loop gang vector collapse(2)
        do k=1,nlev
           do i=1,mgncol
              if (incloud_cldliq_mass_mixing_ratio_kg_per_kg(i,k).ge.qsmall .and. t(i,k).lt.269.15_r8 .and. &
                   local_cldliq_number_tend_immersion_freezing_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k).gt.local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k)) then
                 ! scale mixing ratio of droplet freezing with limit
                 local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(i,k)=local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(i,k)*(local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k)/(local_cldliq_number_tend_immersion_freezing_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)))
                 local_cldliq_number_tend_immersion_freezing_per_kg_per_s(i,k)=local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k)/liquid_cloud_area_fraction(i,k)
              end if
           end do
        end do
        !$acc end parallel

        call contact_freezing(microp_uniform, t, p, rndst, nacon, cldliq_size_dist_shape, cldliq_size_dist_slope_per_m, incloud_cldliq_mass_mixing_ratio_kg_per_kg, incloud_cldliq_number_conc_per_kg, &
                              relvar, local_cldliq_mass_tend_contact_freezing_kg_per_kg_per_s, local_cldliq_number_tend_contact_freezing_per_kg_per_s, mgncol*nlev, mdust)
     else
        ! Mass of droplets frozen is the average droplet mass, except
        ! with two limiters: concentration must be at least 1/cm^3, and
        ! mass must be at least the minimum defined above.

        !$acc parallel vector_length(VLENS) default(present)
        !$acc loop gang vector collapse(2)
        do k=1,nlev
           do i=1,mgncol
              mi0l(i,k) = incloud_cldliq_mass_mixing_ratio_kg_per_kg(i,k)/max(incloud_cldliq_number_conc_per_kg(i,k), 1.0e6_r8/air_density_kg_per_m3(i,k))
              mi0l(i,k) = max(mi0l_min, mi0l(i,k))
              if (incloud_cldliq_mass_mixing_ratio_kg_per_kg(i,k) >= qsmall) then
                 local_cldliq_number_tend_immersion_freezing_per_kg_per_s(i,k) = frzimm(i,k)*1.0e6_r8/air_density_kg_per_m3(i,k)
                 local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(i,k) = local_cldliq_number_tend_immersion_freezing_per_kg_per_s(i,k)*mi0l(i,k)
                 local_cldliq_number_tend_contact_freezing_per_kg_per_s(i,k) = frzcnt(i,k)*1.0e6_r8/air_density_kg_per_m3(i,k)
                 local_cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(i,k) = local_cldliq_number_tend_contact_freezing_per_kg_per_s(i,k)*mi0l(i,k)
                 local_cldliq_number_tend_deposition_nucleation_per_kg_per_s(i,k) = frzdep(i,k)*1.0e6_r8/air_density_kg_per_m3(i,k)
                 local_cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s(i,k) = local_cldliq_number_tend_deposition_nucleation_per_kg_per_s(i,k)*mi0
              else
                 local_cldliq_number_tend_immersion_freezing_per_kg_per_s(i,k) = 0._r8
                 local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(i,k) = 0._r8
                 local_cldliq_number_tend_contact_freezing_per_kg_per_s(i,k) = 0._r8
                 local_cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(i,k) = 0._r8
                 local_cldliq_number_tend_deposition_nucleation_per_kg_per_s(i,k) = 0._r8
                 local_cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s(i,k) = 0._r8
              end if
           end do
        end do
        !$acc end parallel
     end if
  else
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(i,k)=0._r8
           local_cldliq_number_tend_immersion_freezing_per_kg_per_s(i,k)=0._r8
           local_cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(i,k)=0._r8
           local_cldliq_number_tend_contact_freezing_per_kg_per_s(i,k)=0._r8
           local_cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s(i,k)=0._r8
           local_cldliq_number_tend_deposition_nucleation_per_kg_per_s(i,k)=0._r8
        end do
     end do
     !$acc end parallel
  end if


  call snow_self_aggregation(t, air_density_kg_per_m3, snow_fall_speed_prefactor, rhosn, inprecip_snow_mass_mixing_ratio_kg_per_kg, inprecip_snow_number_conc_per_kg, local_snow_number_tend_aggregation_per_kg_per_s, mgncol*nlev)

  call accrete_cloud_water_snow(t, air_density_kg_per_m3, snow_fall_speed_prefactor, snow_number_weighted_fall_speed_m_per_s, air_dynamic_viscosity_kg_per_m_per_s, incloud_cldliq_mass_mixing_ratio_kg_per_kg, incloud_cldliq_number_conc_per_kg, inprecip_snow_mass_mixing_ratio_kg_per_kg, cldliq_size_dist_shape, &
                                cldliq_size_dist_slope_per_m, snow_size_dist_slope_per_m, snow_size_dist_intercept_per_m4, local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s, local_cldliq_number_tend_riming_by_snow_per_kg_per_s, mgncol*nlev)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k) = local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k)*micro_mg_iaccr_factor
        local_cldliq_number_tend_riming_by_snow_per_kg_per_s(i,k) = local_cldliq_number_tend_riming_by_snow_per_kg_per_s(i,k)*micro_mg_iaccr_factor
     end do
  end do
  !$acc end parallel

  if (do_cldice) then
     call secondary_ice_production(t, local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s, local_cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s, local_cldice_number_tend_splintering_by_snow_per_kg_per_s, mgncol*nlev)
  else
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           local_cldice_number_tend_splintering_by_snow_per_kg_per_s(i,k) = 0.0_r8
           local_cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s(i,k) = 0.0_r8
        end do
     end do
     !$acc end parallel
  end if

  call accrete_rain_snow(t, air_density_kg_per_m3, proc_rates%rain_mass_weighted_fall_speed_m_per_s, proc_rates%snow_mass_weighted_fall_speed_m_per_s, rain_number_weighted_fall_speed_m_per_s, snow_number_weighted_fall_speed_m_per_s, inprecip_rain_mass_mixing_ratio_kg_per_kg, inprecip_snow_mass_mixing_ratio_kg_per_kg, rain_size_dist_slope_per_m, &
                         rain_size_dist_intercept_per_m4, snow_size_dist_slope_per_m, snow_size_dist_intercept_per_m4, local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s, local_rain_number_tend_accretion_by_snow_per_kg_per_s, mgncol*nlev)

  call heterogeneous_rain_freezing(t, inprecip_rain_mass_mixing_ratio_kg_per_kg, inprecip_rain_number_conc_per_kg, rain_size_dist_slope_per_m, local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s, local_rain_number_tend_freezing_to_snow_per_kg_per_s, mgncol*nlev)

  if (trim(warm_rain) == 'sb2001') then
     call sb2001v2_accre_cld_water_rain(incloud_cldliq_mass_mixing_ratio_kg_per_kg, incloud_cldliq_number_conc_per_kg, inprecip_rain_mass_mixing_ratio_kg_per_kg, air_density_kg_per_m3, relvar, local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s, local_cldliq_number_tend_accretion_by_rain_per_kg_per_s, mgncol*nlev)
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           proc_rates%cldliq_number_tend_autoconv_accretion_SB2001_per_kg_per_s(i,k)=proc_rates%cldliq_number_tend_autoconv_accretion_SB2001_per_kg_per_s(i,k)+local_cldliq_number_tend_accretion_by_rain_per_kg_per_s(i,k)
           proc_rates%cldliq_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s(i,k)=proc_rates%cldliq_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)
           proc_rates%rain_number_tend_autoconv_accretion_SB2001_per_kg_per_s(i,k)=proc_rates%rain_number_tend_autoconv_accretion_SB2001_per_kg_per_s(i,k)+local_cldliq_number_tend_accretion_by_rain_per_kg_per_s(i,k)  !Sign should be same as prc?
           proc_rates%rain_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s(i,k)=proc_rates%rain_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s(i,k)-local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)
        end do
     end do
     !$acc end parallel
  end if

  if (trim(warm_rain) == 'kk2000') then
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k = 1,nlev
        do i = 1,mgncol
           rtmp(i,k) = inprecip_rain_mass_mixing_ratio_kg_per_kg(i,k)
           ctmp(i,k) = incloud_cldliq_mass_mixing_ratio_kg_per_kg(i,k)
           ntmp(i,k) = incloud_cldliq_number_conc_per_kg(i,k)

           !Option: include recently autoconverted rain (prc, nprc) in accretion
           if (accre_sees_auto) then
              rtmp(i,k) = rtmp(i,k) + local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k)*deltat
              ctmp(i,k) = ctmp(i,k) - local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k)*deltat
              ntmp(i,k) = ntmp(i,k) - local_rain_number_tend_autoconversion_from_cldliq_per_kg_per_s(i,k)*deltat
           endif
        end do
     end do
     !$acc end parallel

     call accrete_cloud_water_rain(microp_uniform,rtmp,ctmp,ntmp,relvar, &
                                   accre_enhan,local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s,local_cldliq_number_tend_accretion_by_rain_per_kg_per_s,mgncol*nlev)

     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k = 1,nlev
        do i = 1,mgncol
           proc_rates%cldliq_number_tend_autoconv_accretion_KK2000_per_kg_per_s(i,k)=proc_rates%cldliq_number_tend_autoconv_accretion_KK2000_per_kg_per_s(i,k)+local_cldliq_number_tend_accretion_by_rain_per_kg_per_s(i,k)
           proc_rates%cldliq_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s(i,k)=proc_rates%cldliq_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)
           proc_rates%rain_number_tend_autoconv_accretion_KK2000_per_kg_per_s(i,k)=proc_rates%rain_number_tend_autoconv_accretion_KK2000_per_kg_per_s(i,k)+local_cldliq_number_tend_accretion_by_rain_per_kg_per_s(i,k)  !Sign consistent with prc,nprc
           proc_rates%rain_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s(i,k)=proc_rates%rain_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s(i,k)-local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)
        end do
     end do
     !$acc end parallel
  endif

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k) = local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)*micro_mg_accre_enhan_fact
        local_cldliq_number_tend_accretion_by_rain_per_kg_per_s(i,k) = local_cldliq_number_tend_accretion_by_rain_per_kg_per_s(i,k)*micro_mg_accre_enhan_fact
     end do
  end do
  !$acc end parallel

  if (trim(warm_rain) == 'kk2000' .or. trim(warm_rain) == 'sb2001') then
     call self_collection_rain(air_density_kg_per_m3, inprecip_rain_mass_mixing_ratio_kg_per_kg, inprecip_rain_number_conc_per_kg, local_rain_number_tend_self_collection_per_kg_per_s, mgncol*nlev)
  end if

  if (do_cldice) then
     call accrete_cloud_ice_snow(t, air_density_kg_per_m3, snow_fall_speed_prefactor, incloud_cldice_mass_mixing_ratio_kg_per_kg, incloud_cldice_number_conc_per_kg, inprecip_snow_mass_mixing_ratio_kg_per_kg, snow_size_dist_slope_per_m, snow_size_dist_intercept_per_m4, local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s, local_cldice_number_tend_accretion_by_snow_per_kg_per_s, mgncol*nlev)
  else
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k) = 0._r8
           local_cldice_number_tend_accretion_by_snow_per_kg_per_s(i,k) = 0._r8
        end do
     end do
     !$acc end parallel
  end if

  call bergeron_process_snow(t, air_density_kg_per_m3, dv, air_dynamic_viscosity_kg_per_m_per_s, schmidt_number, saturation_specific_humidity_over_liquid_kg_per_kg, saturation_specific_humidity_over_ice_kg_per_kg, snow_fall_speed_prefactor, incloud_cldliq_mass_mixing_ratio_kg_per_kg, inprecip_snow_mass_mixing_ratio_kg_per_kg, snow_size_dist_slope_per_m, snow_size_dist_intercept_per_m4, local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s, mgncol*nlev)
  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k)=local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k)*micro_mg_berg_eff_factor
     end do
  end do
  !$acc end parallel


  call vapor_deposition_onto_snow(t, q, snow_mass_mixing_ratio_kg_per_kg, snow_number_conc_per_kg, precipitation_area_fraction, air_density_kg_per_m3, dv, saturation_specific_humidity_over_liquid_kg_per_kg, &
      saturation_specific_humidity_over_ice_kg_per_kg, snow_fall_speed_prefactor, air_dynamic_viscosity_kg_per_m_per_s, schmidt_number, local_snow_mass_tend_deposition_kg_per_kg_per_s, mgncol*nlev)


  if (do_cldice) then
     call ice_deposition_sublimation(t, q, cldice_mass_mixing_ratio_kg_per_kg, cldice_number_conc_per_kg, ice_cloud_area_fraction, air_density_kg_per_m3, dv, saturation_specific_humidity_over_liquid_kg_per_kg, saturation_specific_humidity_over_ice_kg_per_kg, &
                                     local_cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s, local_cldice_mass_tend_deposition_kg_per_kg_per_s, local_cldice_mass_tend_sublimation_kg_per_kg_per_s, mgncol*nlev)
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           local_cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s(i,k)=local_cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s(i,k)*micro_mg_berg_eff_factor
           if (local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k) < 0._r8 .and. cldice_mass_mixing_ratio_kg_per_kg(i,k) > qsmall .and. ice_cloud_area_fraction(i,k) > mincld) then
              local_cldice_number_tend_sublimation_per_kg_per_s(i,k) = sublim_factor*local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k) / cldice_mass_mixing_ratio_kg_per_kg(i,k) * cldice_number_conc_per_kg(i,k) / ice_cloud_area_fraction(i,k)
           else
              local_cldice_number_tend_sublimation_per_kg_per_s(i,k) = 0._r8
           end if

           ! bergeron process should not reduce nc unless
           ! all ql is removed (which is handled elsewhere)
           !in fact, nothing in this entire file makes nsubc nonzero.
           local_cldliq_number_tend_evaporation_per_kg_per_s(i,k) = 0._r8

        end do
     end do
     !$acc end parallel
  end if !do_cldice

! Process rate calls for graupel
!===================================================================

  if (do_hail.or.do_graupel) then
     call graupel_collecting_snow(inprecip_snow_mass_mixing_ratio_kg_per_kg, inprecip_rain_mass_mixing_ratio_kg_per_kg, proc_rates%rain_mass_weighted_fall_speed_m_per_s, proc_rates%snow_mass_weighted_fall_speed_m_per_s, air_density_kg_per_m3, &
                                  rain_size_dist_slope_per_m, rain_size_dist_intercept_per_m4, snow_size_dist_slope_per_m, snow_size_dist_intercept_per_m4, local_snow_mass_tend_collection_by_rain_kg_per_kg_per_s, mgncol*nlev)

     call graupel_collecting_cld_water(inprecip_graupel_mass_mixing_ratio_kg_per_kg, incloud_cldliq_mass_mixing_ratio_kg_per_kg, incloud_cldliq_number_conc_per_kg, air_density_kg_per_m3, graupel_size_dist_intercept_per_m4, graupel_size_dist_slope_per_m, bgtmp, graupel_fall_speed_prefactor, local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s, local_cldliq_number_tend_riming_by_graupel_per_kg_per_s, mgncol*nlev)

     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(i,k) = local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(i,k)*micro_mg_iaccr_factor
           local_cldliq_number_tend_riming_by_graupel_per_kg_per_s(i,k) = local_cldliq_number_tend_riming_by_graupel_per_kg_per_s(i,k)*micro_mg_iaccr_factor
        end do
     end do
     !$acc end parallel

     call graupel_riming_liquid_snow(local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s, inprecip_snow_mass_mixing_ratio_kg_per_kg, incloud_cldliq_mass_mixing_ratio_kg_per_kg, inprecip_snow_number_conc_per_kg, air_density_kg_per_m3, rhosn, rhogtmp, snow_fall_speed_prefactor, &
                                     snow_size_dist_slope_per_m, snow_size_dist_intercept_per_m4, deltat, local_graupel_mass_tend_from_snow_riming_kg_per_kg_per_s, local_graupel_number_tend_from_snow_riming_per_kg_per_s, mgncol*nlev)

     call graupel_collecting_rain(inprecip_rain_mass_mixing_ratio_kg_per_kg, inprecip_graupel_mass_mixing_ratio_kg_per_kg, proc_rates%graupel_mass_weighted_fall_speed_m_per_s, proc_rates%rain_mass_weighted_fall_speed_m_per_s, graupel_number_weighted_fall_speed_m_per_s, rain_number_weighted_fall_speed_m_per_s, air_density_kg_per_m3, rain_size_dist_intercept_per_m4, &
                                  rain_size_dist_slope_per_m, graupel_size_dist_intercept_per_m4, graupel_size_dist_slope_per_m, local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s, local_rain_number_tend_accretion_by_graupel_per_kg_per_s, mgncol*nlev)

     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(i,k) = local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(i,k)*micro_mg_iaccr_factor
           local_rain_number_tend_accretion_by_graupel_per_kg_per_s(i,k) = local_rain_number_tend_accretion_by_graupel_per_kg_per_s(i,k)*micro_mg_iaccr_factor
        end do
     end do
     !$acc end parallel

!AG note: Graupel rain riming snow changes
!    pracs, npracs, (accretion of rain by snow)  psacr (collection of snow by rain)

     call graupel_rain_riming_snow(local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s, local_rain_number_tend_accretion_by_snow_per_kg_per_s, local_snow_mass_tend_collection_by_rain_kg_per_kg_per_s, inprecip_snow_mass_mixing_ratio_kg_per_kg, inprecip_rain_mass_mixing_ratio_kg_per_kg, inprecip_rain_number_conc_per_kg, inprecip_snow_number_conc_per_kg, &
                                   snow_size_dist_intercept_per_m4, snow_size_dist_slope_per_m, rain_size_dist_intercept_per_m4, rain_size_dist_slope_per_m, deltat, local_graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s, local_graupel_number_tend_rain_snow_collection_per_kg_per_s, mgncol*nlev)

     call graupel_rime_splintering(t, incloud_cldliq_mass_mixing_ratio_kg_per_kg, inprecip_rain_mass_mixing_ratio_kg_per_kg, inprecip_graupel_mass_mixing_ratio_kg_per_kg, local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s, local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s, local_cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s, local_cldice_number_tend_splintering_by_graupel_per_kg_per_s, local_cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s, local_cldice_number_tend_splintering_rain_graupel_per_kg_per_s,mgncol*nlev)


     call evaporate_sublimate_precip_graupel(t, air_density_kg_per_m3, dv, air_dynamic_viscosity_kg_per_m_per_s, schmidt_number, q, saturation_specific_humidity_over_liquid_kg_per_kg, saturation_specific_humidity_over_ice_kg_per_kg, liquid_cloud_area_fraction, precipitation_area_fraction, rain_fall_speed_prefactor, snow_fall_speed_prefactor, graupel_fall_speed_prefactor, &
                                             bgtmp, incloud_cldliq_mass_mixing_ratio_kg_per_kg, incloud_cldice_mass_mixing_ratio_kg_per_kg, inprecip_rain_mass_mixing_ratio_kg_per_kg, inprecip_snow_mass_mixing_ratio_kg_per_kg, inprecip_graupel_mass_mixing_ratio_kg_per_kg, rain_size_dist_slope_per_m, rain_size_dist_intercept_per_m4, snow_size_dist_slope_per_m, snow_size_dist_intercept_per_m4, graupel_size_dist_slope_per_m, graupel_size_dist_intercept_per_m4, &
                                             local_rain_mass_tend_evaporation_kg_per_kg_per_s, local_snow_mass_tend_sublimation_kg_per_kg_per_s, local_graupel_mass_tend_deposition_kg_per_kg_per_s, am_evp_st, mgncol*nlev, evap_rhthrsh_ifs)
  else
     ! Routine without Graupel (original)
     call evaporate_sublimate_precip(t, air_density_kg_per_m3, dv, air_dynamic_viscosity_kg_per_m_per_s, schmidt_number, q, saturation_specific_humidity_over_liquid_kg_per_kg, saturation_specific_humidity_over_ice_kg_per_kg, liquid_cloud_area_fraction, precipitation_area_fraction, rain_fall_speed_prefactor, snow_fall_speed_prefactor, incloud_cldliq_mass_mixing_ratio_kg_per_kg, incloud_cldice_mass_mixing_ratio_kg_per_kg, &
                                     inprecip_rain_mass_mixing_ratio_kg_per_kg, inprecip_snow_mass_mixing_ratio_kg_per_kg, rain_size_dist_slope_per_m, rain_size_dist_intercept_per_m4, snow_size_dist_slope_per_m, snow_size_dist_intercept_per_m4, local_rain_mass_tend_evaporation_kg_per_kg_per_s, local_snow_mass_tend_sublimation_kg_per_kg_per_s, am_evp_st, mgncol*nlev, evap_rhthrsh_ifs)
  end if ! end do_graupel/hail loop

! scale precip evaporation to match IFS 'new' version (option 2)
  if (evap_scl_ifs) then
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)= 0.15_r8 * local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)
        end do
     end do
     !$acc end parallel
  end if

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        ! conservation to ensure no negative values of cloud water/precipitation
        ! in case microphysical process rates are large
        !===================================================================

        ! note: for check on conservation, processes are multiplied by omsm
        ! to prevent problems due to round off error

        ! conservation of qc
        !-------------------------------------------------------------------
        dum = ((local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s(i,k)+ &
             local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_from_snow_riming_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)+ &
             local_cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s(i,k))*deltat
        if (dum.gt.cldliq_mass_mixing_ratio_kg_per_kg(i,k)) then
           ratio = cldliq_mass_mixing_ratio_kg_per_kg(i,k)*rdeltat/((local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(i,k)+ &
                local_cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_from_snow_riming_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)+&
                local_cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s(i,k))*omsm
           local_cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s(i,k)=local_cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s(i,k)*ratio
           local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(i,k)=local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(i,k)*ratio
           local_graupel_mass_tend_from_snow_riming_kg_per_kg_per_s(i,k)=local_graupel_mass_tend_from_snow_riming_kg_per_kg_per_s(i,k)*ratio
           local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k) = local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k)*ratio
           local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k) = local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)*ratio
           local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(i,k) = local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(i,k)*ratio
           local_cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(i,k) = local_cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(i,k)*ratio
           local_cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s(i,k) = local_cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s(i,k)*ratio
           local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k) = local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k)*ratio
           local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k) = local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k)*ratio
           local_cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s(i,k) = local_cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s(i,k)*ratio
           qcrat(i,k) = ratio
        else
           qcrat(i,k) = 1._r8
        end if
        !PMC 12/3/12: ratio is also frac of step w/ liquid.
        !thus we apply berg for "ratio" of timestep and vapor
        !deposition for the remaining frac of the timestep.
        if (cldliq_mass_mixing_ratio_kg_per_kg(i,k) >= qsmall) then
           local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k) = local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k)*(1._r8-qcrat(i,k))
           local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k) = local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)*(1._r8-qcrat(i,k))
        end if

        !=================================================================
        ! apply limiter to ensure that ice/snow sublimation and rain evap
        ! don't push conditions into supersaturation, and ice deposition/nucleation don't
        ! push conditions into sub-saturation
        ! note this is done after qc conservation since we don't know how large
        ! vap_dep is before then
        ! estimates are only approximate since other process terms haven't been limited
        ! for conservation yet

        ! first limit ice deposition/nucleation vap_dep + mnuccd + vap_deps
        local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k) = max(0._r8,local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k))
        local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k) = max(0._r8,local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k))
        local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k) = max(0._r8,local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k))

        dum1 = local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k) + local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k) + local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)
        if (dum1 > 1.e-20_r8) then
           dum = (q(i,k)-saturation_specific_humidity_over_ice_kg_per_kg(i,k))/(1._r8 + xxls_squared*saturation_specific_humidity_over_ice_kg_per_kg(i,k)/(cpp*rv*t(i,k)**2))*rdeltat
           dum = max(dum,0._r8)
           if (dum1 > dum) then
              ! Allocate the limited "dum" tendency to mnuccd and vap_dep
              ! processes. Don't divide by cloud fraction; these are grid-
              ! mean rates.
              local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k) = dum*local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k)/dum1
              local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k) = dum*local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k)/dum1
              local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k) = dum*local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)/dum1

           end if
        end if
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        !===================================================================
        ! conservation of nc
        !-------------------------------------------------------------------
        dum = (local_cldliq_number_tend_autoconversion_to_rain_per_kg_per_s(i,k)+local_cldliq_number_tend_accretion_by_rain_per_kg_per_s(i,k)+local_cldliq_number_tend_immersion_freezing_per_kg_per_s(i,k)+local_cldliq_number_tend_contact_freezing_per_kg_per_s(i,k)+ &
               local_cldliq_number_tend_riming_by_snow_per_kg_per_s(i,k)-local_cldliq_number_tend_evaporation_per_kg_per_s(i,k)+local_cldliq_number_tend_riming_by_graupel_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)*deltat

        if (dum.gt.cldliq_number_conc_per_kg(i,k)) then
           ratio = cldliq_number_conc_per_kg(i,k)*rdeltat/((local_cldliq_number_tend_autoconversion_to_rain_per_kg_per_s(i,k)+local_cldliq_number_tend_accretion_by_rain_per_kg_per_s(i,k)+local_cldliq_number_tend_immersion_freezing_per_kg_per_s(i,k)+local_cldliq_number_tend_contact_freezing_per_kg_per_s(i,k)+&
                   local_cldliq_number_tend_riming_by_snow_per_kg_per_s(i,k)-local_cldliq_number_tend_evaporation_per_kg_per_s(i,k)+local_cldliq_number_tend_riming_by_graupel_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k))*omsm
           local_cldliq_number_tend_riming_by_graupel_per_kg_per_s(i,k) = local_cldliq_number_tend_riming_by_graupel_per_kg_per_s(i,k)*ratio
           local_cldliq_number_tend_autoconversion_to_rain_per_kg_per_s(i,k)   = local_cldliq_number_tend_autoconversion_to_rain_per_kg_per_s(i,k)*ratio
           local_cldliq_number_tend_accretion_by_rain_per_kg_per_s(i,k)    = local_cldliq_number_tend_accretion_by_rain_per_kg_per_s(i,k)*ratio
           local_cldliq_number_tend_immersion_freezing_per_kg_per_s(i,k)  = local_cldliq_number_tend_immersion_freezing_per_kg_per_s(i,k)*ratio
           local_cldliq_number_tend_contact_freezing_per_kg_per_s(i,k)  = local_cldliq_number_tend_contact_freezing_per_kg_per_s(i,k)*ratio
           local_cldliq_number_tend_riming_by_snow_per_kg_per_s(i,k) = local_cldliq_number_tend_riming_by_snow_per_kg_per_s(i,k)*ratio
           local_cldliq_number_tend_evaporation_per_kg_per_s(i,k)   = local_cldliq_number_tend_evaporation_per_kg_per_s(i,k)*ratio
        end if
        local_rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(i,k)=0._r8
        local_rain_number_tend_freezing_to_cldice_per_kg_per_s(i,k)=0._r8

        if (do_cldice) then
           ! freezing of rain to produce ice if mean rain size is smaller than Dcs
           if (rain_size_dist_slope_per_m(i,k) > qsmall) then
              if (1._r8/rain_size_dist_slope_per_m(i,k) < Dcs) then
                 local_rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(i,k)=local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k)
                 local_rain_number_tend_freezing_to_cldice_per_kg_per_s(i,k)=local_rain_number_tend_freezing_to_snow_per_kg_per_s(i,k)
                 local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k)=0._r8
                 local_rain_number_tend_freezing_to_snow_per_kg_per_s(i,k)=0._r8
              end if
           end if
        end if
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        ! conservation of rain mixing ratio
        !-------------------------------------------------------------------
        dum = ((-local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)+local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)+local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k)+local_rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(i,k) &
             +local_cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s(i,k)+local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)- &
             (local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k))*deltat
        ! note that qrtend is included below because of instantaneous freezing/melt
        if (dum.gt.rain_mass_mixing_ratio_kg_per_kg(i,k).and. &
             (-local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)+local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)+local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k)+local_rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s(i,k)+local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s(i,k)).ge.qsmall) then
           ratio = (rain_mass_mixing_ratio_kg_per_kg(i,k)*rdeltat+(local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k))/   &
                precipitation_area_fraction(i,k)/(-local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)+local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)+local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k)+local_rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(i,k) &
                +local_cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s(i,k)+local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s(i,k))*omsm
           local_cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s(i,k)= local_cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s(i,k)*ratio
           local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(i,k)=local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(i,k)*ratio
           local_graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s(i,k)=local_graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s(i,k)*ratio
           local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)=local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)*ratio
           local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)=local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)*ratio
           local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k)=local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k)*ratio
           local_rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(i,k)=local_rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(i,k)*ratio
        end if

        ! conservation of rain number
        !-------------------------------------------------------------------
        ! Add evaporation of rain number.
        if (local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k) < 0._r8) then
           local_rain_number_tend_evaporation_per_kg_per_s(i,k) = local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)*rain_number_conc_per_kg(i,k)/rain_mass_mixing_ratio_kg_per_kg(i,k)
        else
           local_rain_number_tend_evaporation_per_kg_per_s(i,k) = 0._r8
        end if

        dum = ((-local_rain_number_tend_evaporation_per_kg_per_s(i,k)+local_rain_number_tend_accretion_by_snow_per_kg_per_s(i,k)+local_rain_number_tend_freezing_to_snow_per_kg_per_s(i,k)+local_rain_number_tend_freezing_to_cldice_per_kg_per_s(i,k)-local_rain_number_tend_self_collection_per_kg_per_s(i,k)+local_rain_number_tend_accretion_by_graupel_per_kg_per_s(i,k)+local_graupel_number_tend_rain_snow_collection_per_kg_per_s(i,k)) &
             *precipitation_area_fraction(i,k)- local_rain_number_tend_autoconversion_from_cldliq_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k))*deltat

        ! Added a check to trap for division by zero errors

        tmpnr = -local_rain_number_tend_evaporation_per_kg_per_s(i,k)+local_rain_number_tend_accretion_by_snow_per_kg_per_s(i,k)+local_rain_number_tend_freezing_to_snow_per_kg_per_s(i,k)+local_rain_number_tend_freezing_to_cldice_per_kg_per_s(i,k)-local_rain_number_tend_self_collection_per_kg_per_s(i,k)+local_rain_number_tend_accretion_by_graupel_per_kg_per_s(i,k)+local_graupel_number_tend_rain_snow_collection_per_kg_per_s(i,k)
        tmpp  = rain_number_conc_per_kg(i,k)*rdeltat + local_rain_number_tend_autoconversion_from_cldliq_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)

        if (dum.gt.rain_number_conc_per_kg(i,k) .and. tmpnr.gt.0._r8 .and. tmpp.gt.0._r8 .and. precipitation_area_fraction(i,k).gt.0._r8) then
           ratio = (rain_number_conc_per_kg(i,k)*rdeltat+local_rain_number_tend_autoconversion_from_cldliq_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k))/precipitation_area_fraction(i,k)/ &
                (-local_rain_number_tend_evaporation_per_kg_per_s(i,k)+local_rain_number_tend_accretion_by_snow_per_kg_per_s(i,k)+local_rain_number_tend_freezing_to_snow_per_kg_per_s(i,k)+local_rain_number_tend_freezing_to_cldice_per_kg_per_s(i,k)-local_rain_number_tend_self_collection_per_kg_per_s(i,k)+local_rain_number_tend_accretion_by_graupel_per_kg_per_s(i,k)+local_graupel_number_tend_rain_snow_collection_per_kg_per_s(i,k))*omsm

           local_rain_number_tend_accretion_by_graupel_per_kg_per_s(i,k)=local_rain_number_tend_accretion_by_graupel_per_kg_per_s(i,k)*ratio
           local_graupel_number_tend_rain_snow_collection_per_kg_per_s(i,k)=local_graupel_number_tend_rain_snow_collection_per_kg_per_s(i,k)*ratio
           local_rain_number_tend_self_collection_per_kg_per_s(i,k)=local_rain_number_tend_self_collection_per_kg_per_s(i,k)*ratio
           local_rain_number_tend_accretion_by_snow_per_kg_per_s(i,k)=local_rain_number_tend_accretion_by_snow_per_kg_per_s(i,k)*ratio
           local_rain_number_tend_freezing_to_snow_per_kg_per_s(i,k)=local_rain_number_tend_freezing_to_snow_per_kg_per_s(i,k)*ratio
           local_rain_number_tend_evaporation_per_kg_per_s(i,k)=local_rain_number_tend_evaporation_per_kg_per_s(i,k)*ratio
           local_rain_number_tend_freezing_to_cldice_per_kg_per_s(i,k)=local_rain_number_tend_freezing_to_cldice_per_kg_per_s(i,k)*ratio
        end if
     end do
  end do
  !$acc end parallel

  if (do_cldice) then
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           ! conservation of qi
           !-------------------------------------------------------------------
           dum = ((-local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(i,k)-local_cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(i,k)-local_cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s(i,k)-local_cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s(i,k)-local_cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)+(local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(i,k)+ &
                local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k))*ice_cloud_area_fraction(i,k)+(-local_cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s(i,k)-local_rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k) &
                -local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k)-local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k)-local_cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s(i,k)-local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k))*deltat
           if (dum.gt.cldice_mass_mixing_ratio_kg_per_kg(i,k)) then
              ratio = (cldice_mass_mixing_ratio_kg_per_kg(i,k)*rdeltat+local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k)+ &
                   (local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)+ &
                   (local_cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s(i,k)+local_rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k))/ &
                   ((local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k))*ice_cloud_area_fraction(i,k)-local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k))*omsm
              local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(i,k) = local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(i,k)*ratio
              local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k) = local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)*ratio
              local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k) = local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k)*ratio
           end if

           ! conservation of ni
           !-------------------------------------------------------------------
           if (use_hetfrz_classnuc) then
              tmpfrz = local_cldliq_number_tend_immersion_freezing_per_kg_per_s(i,k)
           else
              tmpfrz = 0._r8
           end if
           dum = ((-local_cldliq_number_tend_contact_freezing_per_kg_per_s(i,k)-tmpfrz-local_cldliq_number_tend_deposition_nucleation_per_kg_per_s(i,k)-local_cldice_number_tend_splintering_by_snow_per_kg_per_s(i,k)-local_cldice_number_tend_splintering_by_graupel_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)+(local_cldice_number_tend_autoconversion_to_snow_per_kg_per_s(i,k)+ &
                local_cldice_number_tend_accretion_by_snow_per_kg_per_s(i,k)-local_cldice_number_tend_sublimation_per_kg_per_s(i,k))*ice_cloud_area_fraction(i,k)+(-local_cldice_number_tend_splintering_rain_graupel_per_kg_per_s(i,k)-local_rain_number_tend_freezing_to_cldice_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)- &
                local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k))*deltat
           if (dum.gt.cldice_number_conc_per_kg(i,k)) then
              ratio = (cldice_number_conc_per_kg(i,k)*rdeltat+local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k)+ &
                 (local_cldliq_number_tend_contact_freezing_per_kg_per_s(i,k)+tmpfrz+local_cldliq_number_tend_deposition_nucleation_per_kg_per_s(i,k)+local_cldice_number_tend_splintering_by_snow_per_kg_per_s(i,k)+local_cldice_number_tend_splintering_by_graupel_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)+ &
                 (local_rain_number_tend_freezing_to_cldice_per_kg_per_s(i,k)+local_cldice_number_tend_splintering_rain_graupel_per_kg_per_s(i,k))*precipitation_area_fraction(i,k))/ &
                 ((local_cldice_number_tend_autoconversion_to_snow_per_kg_per_s(i,k)+local_cldice_number_tend_accretion_by_snow_per_kg_per_s(i,k)-local_cldice_number_tend_sublimation_per_kg_per_s(i,k))*ice_cloud_area_fraction(i,k))*omsm
              local_cldice_number_tend_autoconversion_to_snow_per_kg_per_s(i,k) = local_cldice_number_tend_autoconversion_to_snow_per_kg_per_s(i,k)*ratio
              local_cldice_number_tend_accretion_by_snow_per_kg_per_s(i,k) = local_cldice_number_tend_accretion_by_snow_per_kg_per_s(i,k)*ratio
              local_cldice_number_tend_sublimation_per_kg_per_s(i,k) = local_cldice_number_tend_sublimation_per_kg_per_s(i,k)*ratio
           end if
        end do
     end do
     !$acc end parallel
  end if

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        ! conservation of snow mixing ratio
        !-------------------------------------------------------------------
        if (do_hail .or. do_graupel) then
        ! NOTE: mnuccr is moved to graupel when active
        ! psacr is a positive value, but a loss for snow
        !HM: psacr is positive in dum (two negatives)
           dum = (-(local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)+local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)-local_snow_mass_tend_collection_by_rain_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)-(local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(i,k))*ice_cloud_area_fraction(i,k) &
             -(local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k) - local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k))*deltat
        else
           dum = (-(local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)+local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)+local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)-(local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(i,k))*ice_cloud_area_fraction(i,k) &
             -(local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k) - local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k))*deltat
        end if
        if (dum.gt.snow_mass_mixing_ratio_kg_per_kg(i,k).and.(local_snow_mass_tend_collection_by_rain_kg_per_kg_per_s(i,k)-local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)).ge.qsmall) then
           if (do_hail .or. do_graupel) then
              ratio = (snow_mass_mixing_ratio_kg_per_kg(i,k)*rdeltat+(local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(i,k))*ice_cloud_area_fraction(i,k)+ &
                   (local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)+local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)+local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)*precipitation_area_fraction(i,k))/ &
                   precipitation_area_fraction(i,k)/(local_snow_mass_tend_collection_by_rain_kg_per_kg_per_s(i,k)-local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k))*omsm
              local_snow_mass_tend_collection_by_rain_kg_per_kg_per_s(i,k)=local_snow_mass_tend_collection_by_rain_kg_per_kg_per_s(i,k)*ratio
           else
              ratio = (snow_mass_mixing_ratio_kg_per_kg(i,k)*rdeltat+(local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(i,k))*ice_cloud_area_fraction(i,k)+ &
                   (local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)+local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)+(local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)+local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k))/ &
                   precipitation_area_fraction(i,k)/(-local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k))*omsm
           end if
           local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)=local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)*ratio
        end if

        ! conservation of snow number
        !-------------------------------------------------------------------
        ! calculate loss of number due to sublimation
        ! for now neglect sublimation of ns
        local_snow_number_tend_sublimation_per_kg_per_s(i,k)=0._r8
        if (do_hail .or. do_graupel) then
           dum = ((-local_snow_number_tend_aggregation_per_kg_per_s(i,k)-local_snow_number_tend_sublimation_per_kg_per_s(i,k)+local_graupel_number_tend_rain_snow_collection_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)-local_cldice_number_tend_autoconversion_to_snow_per_kg_per_s(i,k)*ice_cloud_area_fraction(i,k)+local_graupel_number_tend_from_snow_riming_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k))*deltat
        else
           dum = ((-local_snow_number_tend_aggregation_per_kg_per_s(i,k)-local_snow_number_tend_sublimation_per_kg_per_s(i,k)-local_rain_number_tend_freezing_to_snow_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)-local_cldice_number_tend_autoconversion_to_snow_per_kg_per_s(i,k)*ice_cloud_area_fraction(i,k))*deltat
        end if
        if (dum.gt.snow_number_conc_per_kg(i,k)) then
           if (do_hail .or. do_graupel) then
              ratio = (snow_number_conc_per_kg(i,k)*rdeltat+local_cldice_number_tend_autoconversion_to_snow_per_kg_per_s(i,k)*ice_cloud_area_fraction(i,k))/precipitation_area_fraction(i,k)/ &
                   (-local_snow_number_tend_sublimation_per_kg_per_s(i,k)-local_snow_number_tend_aggregation_per_kg_per_s(i,k)+local_graupel_number_tend_rain_snow_collection_per_kg_per_s(i,k)+liquid_cloud_area_fraction(i,k)/precipitation_area_fraction(i,k)*local_graupel_number_tend_from_snow_riming_per_kg_per_s(i,k))*omsm
              local_graupel_number_tend_from_snow_riming_per_kg_per_s(i,k)=local_graupel_number_tend_from_snow_riming_per_kg_per_s(i,k)*ratio
              local_graupel_number_tend_rain_snow_collection_per_kg_per_s(i,k)=local_graupel_number_tend_rain_snow_collection_per_kg_per_s(i,k)*ratio
           else
              ratio = (snow_number_conc_per_kg(i,k)*rdeltat+local_rain_number_tend_freezing_to_snow_per_kg_per_s(i,k)* &
                   precipitation_area_fraction(i,k)+local_cldice_number_tend_autoconversion_to_snow_per_kg_per_s(i,k)*ice_cloud_area_fraction(i,k))/precipitation_area_fraction(i,k)/ &
                   (-local_snow_number_tend_sublimation_per_kg_per_s(i,k)-local_snow_number_tend_aggregation_per_kg_per_s(i,k))*omsm
           endif
           local_snow_number_tend_sublimation_per_kg_per_s(i,k)=local_snow_number_tend_sublimation_per_kg_per_s(i,k)*ratio
           local_snow_number_tend_aggregation_per_kg_per_s(i,k)=local_snow_number_tend_aggregation_per_kg_per_s(i,k)*ratio
        end if
     end do
  end do
  !$acc end parallel

! Graupel Conservation Checks
!-------------------------------------------------------------------

  if (do_hail.or.do_graupel) then
     ! conservation of graupel mass
     !-------------------------------------------------------------------
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           dum= ((-local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(i,k)-local_graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s(i,k)-local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k)-local_snow_mass_tend_collection_by_rain_kg_per_kg_per_s(i,k)-local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k) &
                + (-local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(i,k)-local_graupel_mass_tend_from_snow_riming_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k))*deltat
           if (dum.gt.graupel_mass_mixing_ratio_kg_per_kg(i,k)) then
              ! note: prdg is always negative (like prds), so it needs to be subtracted in ratio
              ratio = (graupel_mass_mixing_ratio_kg_per_kg(i,k)*rdeltat + (local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s(i,k)+local_snow_mass_tend_collection_by_rain_kg_per_kg_per_s(i,k)+local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k) &
                       + (local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_from_snow_riming_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)) / ((-local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)) * omsm
              local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k)= local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k)*ratio
           end if
        end do
     end do
     !$acc end parallel
     ! conservation of graupel number: not needed, no sinks
     !-------------------------------------------------------------------
  end if

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        ! next limit ice and snow sublimation and rain evaporation
        ! get estimate of q and t at end of time step
        ! don't include other microphysical processes since they haven't
        ! been limited via conservation checks yet
        qtmpAI(i,k)=q(i,k)-(local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k)+local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)+ &
                (local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)+local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k))*deltat
        ttmpA(i,k)=t(i,k)+((local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)*precipitation_area_fraction(i,k))*xxlv+ &
             ((local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)+local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k)+local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k))*xxls)*deltat/cpp
     end do
  end do
  !$acc end parallel

  ! use rhw to allow ice supersaturation
  call qsat_water(ttmpA, p, esnA, qvnAI, mgncol*nlev)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        if ((local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)+local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)+local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k) < -1.e-20_r8) then
           ! modify ice/precip evaporation rate if q > qsat
           if (qtmpAI(i,k) > qvnAI(i,k)) then
              dum1A(i,k)=local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)/((local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)+local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)+local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k))
              dum2A(i,k)=local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)/((local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)+local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)+local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k))
              dum3A(i,k)=local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)/((local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)+local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)+local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k))
              ! recalculate q and t after vap_dep and mnuccd but without evap or sublim
              ttmpA(i,k)=t(i,k)+((local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k)+local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k))*xxls)*deltat/cpp
              dum_2D(i,k)=q(i,k)-(local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k)+local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k))*deltat
           end if
        end if
     end do
  end do
  !$acc end parallel

  ! use rhw to allow ice supersaturation
  call qsat_water(ttmpA, p, esnA, qvnA, mgncol*nlev)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        if ((local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)+local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)+local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k) < -1.e-20_r8) then
           ! modify ice/precip evaporation rate if q > qsat
           if (qtmpAI(i,k) > qvnAI(i,k)) then
              dum=(dum_2D(i,k)-qvnA(i,k))/(1._r8 + xxlv_squared*qvnA(i,k)/(cpp*rv*ttmpA(i,k)**2))
              dum=min(dum,0._r8)
              ! modify rates if needed, divide by precip_frac to get local (in-precip) value
              local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)=dum*dum1A(i,k)*rdeltat/precipitation_area_fraction(i,k)
           end if
        end if
     end do
  end do
  !$acc end parallel

  ! do separately using RHI for prds and ice_sublim
  call qsat_ice(ttmpA, p, esnA, qvnA, mgncol*nlev)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        if ((local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)+local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)+local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k) < -1.e-20_r8) then
           ! modify ice/precip evaporation rate if q > qsat
           if (qtmpAI(i,k) > qvnAI(i,k)) then
              dum=(dum_2D(i,k)-qvnA(i,k))/(1._r8 + xxls_squared*qvnA(i,k)/(cpp*rv*ttmpA(i,k)**2))
              dum=min(dum,0._r8)
              ! modify rates if needed, divide by precip_frac to get local (in-precip) value
              local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k) = dum*dum2A(i,k)*rdeltat/precipitation_area_fraction(i,k)
              local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k) = dum*dum3A(i,k)*rdeltat/precipitation_area_fraction(i,k)
              ! don't divide ice_sublim by cloud fraction since it is grid-averaged
              dum1A(i,k) = (1._r8-dum1A(i,k)-dum2A(i,k)-dum3A(i,k))
              local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k) = dum*dum1A(i,k)*rdeltat
           end if
        end if

        ! get tendencies due to microphysical conversion processes
        !==========================================================
        ! note: tendencies are multiplied by appropriate cloud/precip
        ! fraction to get grid-scale values
        ! note: vap_dep is already grid-average values

        ! The net tendencies need to be added to rather than overwritten,
        ! because they may have a value already set for instantaneous
        ! melting/freezing.
        qvlat(i,k) = qvlat(i,k)-(local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)+local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)-&
             local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k)-local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)-local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k)-local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k)-local_cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k) &
             -local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        tlat(i,k) = tlat(i,k)+((local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)*precipitation_area_fraction(i,k))*xxlv+ &
             ((local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)+local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k)+local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k)+ &
             local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k))*xxls+ &
             ((local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(i,k)+ &
             local_cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_from_snow_riming_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)+ &
             (local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k)+local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)+local_rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(i,k)+local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)+ &
             local_cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s(i,k))*xlf)
        qctend(i,k) = qctend(i,k)+ &
             (-local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)-local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k)-local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(i,k)-local_cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(i,k)-local_cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s(i,k)- &
             local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k)-local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k)-local_cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s(i,k)-local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(i,k)-local_graupel_mass_tend_from_snow_riming_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)-local_cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s(i,k)

        if (do_cldice) then
           qitend(i,k) = qitend(i,k)+ &
                (local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)+(-local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(i,k)- &
                local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k))*ice_cloud_area_fraction(i,k)+local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k)+ &
                local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k)+(local_rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)
        end if

        qrtend(i,k) = qrtend(i,k)+ &
             (local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)+(local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)-local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)- &
             local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k)-local_rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(i,k)-local_cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s(i,k)-local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(i,k)-local_graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)

        if (do_hail.or.do_graupel) then
           qgtend(i,k) = qgtend(i,k) + (local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k)+local_snow_mass_tend_collection_by_rain_kg_per_kg_per_s(i,k)+local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k) &
                + (local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_from_snow_riming_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)
           qstend(i,k) = qstend(i,k)+ &
                (local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(i,k))*ice_cloud_area_fraction(i,k)+(local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)+(local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)+ &
                local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)-local_snow_mass_tend_collection_by_rain_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)+local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)
        else
           !necessary since mnuccr moved to graupel
           qstend(i,k) = qstend(i,k)+ &
                (local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(i,k))*ice_cloud_area_fraction(i,k)+(local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)+(local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)+ &
                local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)+local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)+local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)
        end if
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        cmeout(i,k) = local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k) + local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k) + local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k) + local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)
        ! add output for cmei (accumulate)
        proc_rates%cldice_mass_tend_net_deposition_kg_per_kg_per_s(i,k) = local_cldice_mass_tend_deposition_kg_per_kg_per_s(i,k) + local_cldice_mass_tend_sublimation_kg_per_kg_per_s(i,k) + local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k) + local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)
        !-------------------------------------------------------------------
        ! evaporation/sublimation is stored here as positive term
        ! Add to evapsnow via prdg
        proc_rates%snow_mass_tend_sublimation_kg_per_kg_per_s(i,k) = (-local_snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)-local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)
        nevapr(i,k) = -local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        prer_evap(i,k) = -local_rain_mass_tend_evaporation_kg_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        ! change to make sure prain is positive: do not remove snow from
        ! prain used for wet deposition
        prain(i,k) = (local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)+(-local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)- &
             local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k)-local_rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)
        if (do_hail .or. do_graupel) then
           proc_rates%snow_mass_tend_production_kg_per_kg_per_s(i,k) = (local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(i,k))*ice_cloud_area_fraction(i,k)+(local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)+(&
                local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)+local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)
        else
           proc_rates%snow_mass_tend_production_kg_per_kg_per_s(i,k) = (local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)+local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(i,k))*ice_cloud_area_fraction(i,k)+(local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)+(&
                local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)+local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)+local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)
        end if
        ! following are used to calculate 1st order conversion rate of cloud water
        !    to rain and snow (1/s), for later use in aerosol wet removal routine
        ! previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may be smaller than the qc
        !    used to calculate pra, prc, ... in this routine
        ! qcsinksum_rate1ord = { rate of direct transfer of cloud water to rain & snow }
        !                      (no cloud ice or bergeron terms)
        qcsinksum_rate1ord(i,k) = (local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k)+local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(i,k)+local_graupel_mass_tend_from_snow_riming_kg_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)
        ! Avoid zero/near-zero division.
        qcsinksum_rate1ord(i,k) = qcsinksum_rate1ord(i,k) / &
             max(cldliq_mass_mixing_ratio_kg_per_kg(i,k),1.0e-30_r8)
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        ! microphysics output, note this is grid-averaged
        proc_rates%cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)     = local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k)     = local_cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(i,k)  = local_cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s(i,k)  = local_cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(i,k)  = local_cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s(i,k)  = local_cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k)  = local_cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k)   = local_cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%snow_mass_tend_deposition_kg_per_kg_per_s(i,k) = local_snow_mass_tend_deposition_kg_per_kg_per_s(i,k)
        proc_rates%cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s(i,k)    = local_cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s(i,k)
        proc_rates%cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(i,k)    = local_cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(i,k)*ice_cloud_area_fraction(i,k)
        proc_rates%cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)    = local_cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)*ice_cloud_area_fraction(i,k)
        proc_rates%cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k)  = local_cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(i,k)*ice_cloud_area_fraction(i,k)
        proc_rates%rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)   = local_rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        proc_rates%rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k)  = local_rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        proc_rates%rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(i,k) = local_rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        proc_rates%snow_mass_tend_collection_by_rain_kg_per_kg_per_s(i,k)   = local_snow_mass_tend_collection_by_rain_kg_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        proc_rates%rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(i,k)   = local_rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        proc_rates%cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(i,k)  = local_cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%graupel_mass_tend_from_snow_riming_kg_per_kg_per_s(i,k)  = local_graupel_mass_tend_from_snow_riming_kg_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s(i,k)  = local_graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        proc_rates%graupel_mass_tend_deposition_kg_per_kg_per_s(i,k)    = local_graupel_mass_tend_deposition_kg_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        proc_rates%cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s(i,k)  = local_cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s(i,k) = local_cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        proc_rates%rain_number_tend_accretion_by_graupel_per_kg_per_s(i,k)  = local_rain_number_tend_accretion_by_graupel_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        proc_rates%graupel_number_tend_from_snow_riming_per_kg_per_s(i,k)   = local_graupel_number_tend_from_snow_riming_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%graupel_number_tend_rain_snow_collection_per_kg_per_s(i,k)  = local_graupel_number_tend_rain_snow_collection_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        proc_rates%cldice_number_tend_splintering_by_graupel_per_kg_per_s(i,k)  = local_cldice_number_tend_splintering_by_graupel_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%cldice_number_tend_splintering_rain_graupel_per_kg_per_s(i,k) = local_cldice_number_tend_splintering_rain_graupel_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        proc_rates%cldliq_number_tend_riming_by_graupel_per_kg_per_s(i,k) = local_cldliq_number_tend_riming_by_graupel_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)

        proc_rates%cldliq_number_tend_immersion_freezing_per_kg_per_s(i,k) = local_cldliq_number_tend_immersion_freezing_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%cldliq_number_tend_contact_freezing_per_kg_per_s(i,k) = local_cldliq_number_tend_contact_freezing_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k) = local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k)*ice_cloud_area_fraction(i,k)
        proc_rates%cldliq_number_tend_deposition_nucleation_per_kg_per_s(i,k) = local_cldliq_number_tend_deposition_nucleation_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%rain_number_tend_freezing_to_snow_per_kg_per_s(i,k) = local_rain_number_tend_freezing_to_snow_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        proc_rates%rain_number_tend_freezing_to_cldice_per_kg_per_s(i,k) = local_rain_number_tend_freezing_to_cldice_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        proc_rates%cldice_number_tend_splintering_by_snow_per_kg_per_s(i,k) = local_cldice_number_tend_splintering_by_snow_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%cldliq_number_tend_accretion_by_rain_per_kg_per_s(i,k) = local_cldliq_number_tend_accretion_by_rain_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%cldliq_number_tend_riming_by_snow_per_kg_per_s(i,k) = local_cldliq_number_tend_riming_by_snow_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%cldice_number_tend_accretion_by_snow_per_kg_per_s(i,k) = local_cldice_number_tend_accretion_by_snow_per_kg_per_s(i,k)*ice_cloud_area_fraction(i,k)
        proc_rates%rain_number_tend_accretion_by_snow_per_kg_per_s(i,k) = local_rain_number_tend_accretion_by_snow_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        proc_rates%rain_number_tend_autoconversion_from_cldliq_per_kg_per_s(i,k) = local_rain_number_tend_autoconversion_from_cldliq_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
        proc_rates%rain_number_tend_self_collection_per_kg_per_s(i,k) = local_rain_number_tend_self_collection_per_kg_per_s(i,k)*precipitation_area_fraction(i,k)
        proc_rates%cldice_number_tend_autoconversion_to_snow_per_kg_per_s(i,k) = local_cldice_number_tend_autoconversion_to_snow_per_kg_per_s(i,k)*ice_cloud_area_fraction(i,k)
        proc_rates%snow_number_tend_melting_per_kg_per_s(i,k) = snow_number_change_instant_melting_per_kg(i,k)/deltat
        proc_rates%graupel_number_tend_melting_per_kg_per_s(i,k) = graupel_number_change_instant_melting_per_kg(i,k)/deltat
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        nctend(i,k) = nctend(i,k)+&
           (-local_cldliq_number_tend_immersion_freezing_per_kg_per_s(i,k)-local_cldliq_number_tend_contact_freezing_per_kg_per_s(i,k)-local_cldliq_number_tend_riming_by_snow_per_kg_per_s(i,k)+local_cldliq_number_tend_evaporation_per_kg_per_s(i,k) &
           -local_cldliq_number_tend_accretion_by_rain_per_kg_per_s(i,k)-local_cldliq_number_tend_autoconversion_to_rain_per_kg_per_s(i,k)-local_cldliq_number_tend_riming_by_graupel_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)

        if (do_cldice) then
           if (use_hetfrz_classnuc) then
              tmpfrz = local_cldliq_number_tend_immersion_freezing_per_kg_per_s(i,k)
           else
              tmpfrz = 0._r8
           end if
           nitend(i,k) = nitend(i,k)+ local_cldice_number_tend_nucleation_from_vapor_per_kg_per_s(i,k)+ &
                (local_cldliq_number_tend_contact_freezing_per_kg_per_s(i,k)+tmpfrz+local_cldliq_number_tend_deposition_nucleation_per_kg_per_s(i,k)+local_cldice_number_tend_splintering_by_snow_per_kg_per_s(i,k)+local_cldice_number_tend_splintering_by_graupel_per_kg_per_s(i,k))*liquid_cloud_area_fraction(i,k)+(local_cldice_number_tend_sublimation_per_kg_per_s(i,k)-local_cldice_number_tend_autoconversion_to_snow_per_kg_per_s(i,k)- &
                local_cldice_number_tend_accretion_by_snow_per_kg_per_s(i,k))*ice_cloud_area_fraction(i,k)+(local_rain_number_tend_freezing_to_cldice_per_kg_per_s(i,k)+local_cldice_number_tend_splintering_rain_graupel_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)
        end if

        if(do_graupel.or.do_hail) then
           nstend(i,k) = nstend(i,k)+(local_snow_number_tend_sublimation_per_kg_per_s(i,k)+ &
                local_snow_number_tend_aggregation_per_kg_per_s(i,k)-local_graupel_number_tend_rain_snow_collection_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)+local_cldice_number_tend_autoconversion_to_snow_per_kg_per_s(i,k)*ice_cloud_area_fraction(i,k)-local_graupel_number_tend_from_snow_riming_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)
           ngtend(i,k) = ngtend(i,k)+local_graupel_number_tend_from_snow_riming_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)+(local_graupel_number_tend_rain_snow_collection_per_kg_per_s(i,k)+local_rain_number_tend_freezing_to_snow_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)
        else
           !necessary since mnuccr moved to graupel
           nstend(i,k) = nstend(i,k)+(local_snow_number_tend_sublimation_per_kg_per_s(i,k)+ &
                local_snow_number_tend_aggregation_per_kg_per_s(i,k)+local_rain_number_tend_freezing_to_snow_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)+local_cldice_number_tend_autoconversion_to_snow_per_kg_per_s(i,k)*ice_cloud_area_fraction(i,k)
        end if

        nrtend(i,k) = nrtend(i,k)+ &
             local_rain_number_tend_autoconversion_from_cldliq_per_kg_per_s(i,k)*liquid_cloud_area_fraction(i,k)+(local_rain_number_tend_evaporation_per_kg_per_s(i,k)-local_rain_number_tend_accretion_by_snow_per_kg_per_s(i,k)-local_rain_number_tend_freezing_to_snow_per_kg_per_s(i,k) &
             -local_rain_number_tend_freezing_to_cldice_per_kg_per_s(i,k)+local_rain_number_tend_self_collection_per_kg_per_s(i,k)-local_rain_number_tend_accretion_by_graupel_per_kg_per_s(i,k)-local_graupel_number_tend_rain_snow_collection_per_kg_per_s(i,k))*precipitation_area_fraction(i,k)

        !-----------------------------------------------------
        ! convert rain/snow q and N for output to history, note,
        ! output is for gridbox average

        qrout(i,k) = rain_mass_mixing_ratio_kg_per_kg(i,k)
        nrout(i,k) = rain_number_conc_per_kg(i,k) * air_density_kg_per_m3(i,k)
        qsout(i,k) = snow_mass_mixing_ratio_kg_per_kg(i,k)
        nsout(i,k) = snow_number_conc_per_kg(i,k) * air_density_kg_per_m3(i,k)
        qgout(i,k) = graupel_mass_mixing_ratio_kg_per_kg(i,k)
        ngout(i,k) = graupel_number_conc_per_kg(i,k) * air_density_kg_per_m3(i,k)
     end do
  end do
  !$acc end parallel

  ! calculate n0r and lamr from rain mass and number
  ! divide by precip fraction to get in-precip (local) values of
  ! rain mass and number, divide by rhow to get rain number in kg^-1
  call size_dist_param_basic(mg_rain_props, inprecip_rain_mass_mixing_ratio_kg_per_kg, inprecip_rain_number_conc_per_kg, rain_size_dist_slope_per_m, mgncol, nlev, n0=rain_size_dist_intercept_per_m4)

  ! Calculate rercld
  ! calculate mean size of combined rain and cloud water
  call calc_rercld(rain_size_dist_slope_per_m, rain_size_dist_intercept_per_m4, cldliq_size_dist_slope_per_m, cldliq_size_dist_shape, inprecip_rain_mass_mixing_ratio_kg_per_kg, incloud_cldliq_mass_mixing_ratio_kg_per_kg, incloud_cldliq_number_conc_per_kg, rercld, mgncol*nlev)

  ! Assign variables back to start-of-timestep values
  ! Some state variables are changed before the main microphysics loop
  ! to make "instantaneous" adjustments. Afterward, we must move those changes
  ! back into the tendencies.
  ! These processes:
  !  - Droplet activation (npccn, impacts nc)
  !  - Instantaneous snow melting  (minstsm/ninstsm, impacts qr/qs/nr/ns)
  !  - Instantaneous rain freezing (minstfr/ninstrf, impacts qr/qs/nr/ns)
  !================================================================================
  ! Re-apply droplet activation tendency

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        cldliq_number_conc_per_kg(i,k) = ncn(i,k)
        nctend(i,k) = nctend(i,k) + npccn(i,k)
        ! Re-apply rain freezing and snow melting.
        dum_2D(i,k) = snow_mass_mixing_ratio_kg_per_kg(i,k)
        snow_mass_mixing_ratio_kg_per_kg(i,k)     = qsn(i,k)
        qstend(i,k) = qstend(i,k) + (dum_2D(i,k)-snow_mass_mixing_ratio_kg_per_kg(i,k))*rdeltat

        dum_2D(i,k) = snow_number_conc_per_kg(i,k)
        snow_number_conc_per_kg(i,k)     = nsn(i,k)
        nstend(i,k) = nstend(i,k) + (dum_2D(i,k)-snow_number_conc_per_kg(i,k))*rdeltat

        dum_2D(i,k) = rain_mass_mixing_ratio_kg_per_kg(i,k)
        rain_mass_mixing_ratio_kg_per_kg(i,k)     = qrn(i,k)
        qrtend(i,k) = qrtend(i,k) + (dum_2D(i,k)-rain_mass_mixing_ratio_kg_per_kg(i,k))*rdeltat

        dum_2D(i,k) = rain_number_conc_per_kg(i,k)
        rain_number_conc_per_kg(i,k)     = nrn(i,k)
        nrtend(i,k) = nrtend(i,k) + (dum_2D(i,k)-rain_number_conc_per_kg(i,k))*rdeltat

        ! Re-apply graupel freezing/melting
        dum_2D(i,k) = graupel_mass_mixing_ratio_kg_per_kg(i,k)
        graupel_mass_mixing_ratio_kg_per_kg(i,k)     = qgr(i,k)
        qgtend(i,k) = qgtend(i,k) + (dum_2D(i,k)-graupel_mass_mixing_ratio_kg_per_kg(i,k))*rdeltat

        dum_2D(i,k) = graupel_number_conc_per_kg(i,k)
        graupel_number_conc_per_kg(i,k)     = ngr(i,k)
        ngtend(i,k) = ngtend(i,k) + (dum_2D(i,k)-graupel_number_conc_per_kg(i,k))*rdeltat
        !.............................................................................
        !================================================================================
        ! modify to include snow. in prain & evap (diagnostic here: for wet dep)
        nevapr(i,k) = nevapr(i,k) + proc_rates%snow_mass_tend_sublimation_kg_per_kg_per_s(i,k)
        prain(i,k) = prain(i,k) + proc_rates%snow_mass_tend_production_kg_per_kg_per_s(i,k)
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        ! calculate sedimentation for cloud water and ice
        ! and Graupel (mg3)
        !================================================================================
        ! update in-cloud cloud mixing ratio and number concentration
        ! with microphysical tendencies to calculate sedimentation, assign to dummy vars
        ! note: these are in-cloud values***, hence we divide by cloud fraction
        dumc(i,k)  = (cldliq_mass_mixing_ratio_kg_per_kg(i,k)+qctend(i,k)*deltat)/liquid_cloud_area_fraction(i,k)
        dumi(i,k)  = (cldice_mass_mixing_ratio_kg_per_kg(i,k)+qitend(i,k)*deltat)/ice_cloud_area_fraction(i,k)
        dumnc(i,k) = max((cldliq_number_conc_per_kg(i,k)+nctend(i,k)*deltat)/liquid_cloud_area_fraction(i,k),0._r8)
        dumni(i,k) = max((cldice_number_conc_per_kg(i,k)+nitend(i,k)*deltat)/ice_cloud_area_fraction(i,k),0._r8)

        dumr(i,k)  = (rain_mass_mixing_ratio_kg_per_kg(i,k)+qrtend(i,k)*deltat)/precipitation_area_fraction(i,k)
        dumnr(i,k) = max((rain_number_conc_per_kg(i,k)+nrtend(i,k)*deltat)/precipitation_area_fraction(i,k),0._r8)
        dums(i,k)  = (snow_mass_mixing_ratio_kg_per_kg(i,k)+qstend(i,k)*deltat)/precipitation_area_fraction(i,k)
        dumns(i,k) = max((snow_number_conc_per_kg(i,k)+nstend(i,k)*deltat)/precipitation_area_fraction(i,k),0._r8)

        dumg(i,k)  = (graupel_mass_mixing_ratio_kg_per_kg(i,k)+qgtend(i,k)*deltat)/precipitation_area_fraction(i,k)
        dumng(i,k) = max((graupel_number_conc_per_kg(i,k)+ngtend(i,k)*deltat)/precipitation_area_fraction(i,k),0._r8)

        ! switch for specification of droplet and crystal number
        if (ngcons) then
           dumng(i,k)=ngnst/air_density_kg_per_m3(i,k)
        end if

        ! switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(i,k)=ncnst/air_density_kg_per_m3(i,k)
        end if

        ! switch for specification of cloud ice number
        if (nicons) then
           dumni(i,k)=ninst/air_density_kg_per_m3(i,k)
        end if

        ! switch for specification of constant number
        if (nscons) then
            dumns(i,k)=nsnst/air_density_kg_per_m3(i,k)
        end if

        ! switch for specification of constant number
        if (nrcons) then
            dumnr(i,k)=nrnst/air_density_kg_per_m3(i,k)
        end if
     end do
  end do
  !$acc end parallel

  ! obtain new slope parameter to avoid possible singularity
  call size_dist_param_basic(mg_ice_props, dumi, dumni, cldice_size_dist_slope_per_m, mgncol, nlev)
  call size_dist_param_liq(mg_liq_props, dumc, dumnc, air_density_kg_per_m3, cldliq_size_dist_shape, cldliq_size_dist_slope_per_m, mgncol, nlev)

  ! fallspeed for rain
  call size_dist_param_basic(mg_rain_props, dumr, dumnr, rain_size_dist_slope_per_m, mgncol, nlev)
  ! fallspeed for snow
  call size_dist_param_basic(mg_snow_props, dums, dumns, snow_size_dist_slope_per_m, mgncol, nlev)
  ! fallspeed for graupel
  if (do_hail) then
     call size_dist_param_basic(mg_hail_props, dumg, dumng, graupel_size_dist_slope_per_m, mgncol, nlev)
  end if
  if (do_graupel) then
     call size_dist_param_basic(mg_graupel_props, dumg, dumng, graupel_size_dist_slope_per_m, mgncol, nlev)
  end if

  if ( do_implicit_fall ) then
!    calculate interface height for implicit sedimentation
!    uses Hypsometric equation

     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector
     do i=1,mgncol
        zint(i,nlev+1)=0._r8
        !$acc loop seq
        do k = nlev,1,-1
           H = r*t(i,k)/g*log(pint(i,k+1)/pint(i,k))
           zint(i,k)=zint(i,k+1)+H
        enddo
     enddo
     !$acc end parallel
  end if

  !$acc parallel vector_length(VLENS) default(present) async(LQUEUE)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        ! calculate number and mass weighted fall velocity for droplets and cloud ice
        !-------------------------------------------------------------------
        if (dumc(i,k).ge.qsmall) then
           dum1 = 4._r8+bc+cldliq_size_dist_shape(i,k)
           dum2 = cldliq_size_dist_shape(i,k)+4._r8
           proc_rates%cldliq_mass_weighted_fall_speed_m_per_s(i,k)=cldliq_fall_speed_prefactor(i,k)*gamma(dum1)/(cldliq_size_dist_slope_per_m(i,k)**bc*gamma(dum2))
           ! Following ifs, no condensate sedimentation
           if (ifs_sed) then
              cldliq_mass_pressure_fall_speed_Pa_per_s(i,k)  = 0._r8
              cldliq_number_pressure_fall_speed_Pa_per_s(i,k) = 0._r8
           else
              dum3     = 1._r8+bc+cldliq_size_dist_shape(i,k)
              dum4     = cldliq_size_dist_shape(i,k)+1._r8
              cldliq_mass_pressure_fall_speed_Pa_per_s(i,k)  = g*air_density_kg_per_m3(i,k)*proc_rates%cldliq_mass_weighted_fall_speed_m_per_s(i,k)
              cldliq_number_pressure_fall_speed_Pa_per_s(i,k) = g*air_density_kg_per_m3(i,k)* &
                   cldliq_fall_speed_prefactor(i,k)*gamma(dum3)/ &
                   (cldliq_size_dist_slope_per_m(i,k)**bc*gamma(dum4))
           end if
        else
           cldliq_mass_pressure_fall_speed_Pa_per_s(i,k) = 0._r8
           cldliq_number_pressure_fall_speed_Pa_per_s(i,k)= 0._r8
        end if

        ! redefine dummy variables - sedimentation is calculated over grid-scale
        ! quantities to ensure conservation
        dumc(i,k) = (cldliq_mass_mixing_ratio_kg_per_kg(i,k)+qctend(i,k)*deltat)
        dumnc(i,k) = max((cldliq_number_conc_per_kg(i,k)+nctend(i,k)*deltat),0._r8)
        if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present) async(IQUEUE)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        ! calculate number and mass weighted fall velocity for cloud ice
        if (dumi(i,k).ge.qsmall) then
           proc_rates%cldice_mass_weighted_fall_speed_m_per_s(i,k)=min(cldice_fall_speed_prefactor(i,k)*gamma_bi_plus4/(6._r8*cldice_size_dist_slope_per_m(i,k)**bi), &
                1.2_r8*fall_speed_air_density_correction(i,k))
           proc_rates%cldice_mass_weighted_fall_speed_m_per_s(i,k)=proc_rates%cldice_mass_weighted_fall_speed_m_per_s(i,k)*micro_mg_vtrmi_factor

           cldice_mass_pressure_fall_speed_Pa_per_s(i,k) = g*air_density_kg_per_m3(i,k)*proc_rates%cldice_mass_weighted_fall_speed_m_per_s(i,k)
           cldice_number_pressure_fall_speed_Pa_per_s(i,k) = g*air_density_kg_per_m3(i,k)* &
                min(cldice_fall_speed_prefactor(i,k)*gamma_bi_plus1/cldice_size_dist_slope_per_m(i,k)**bi,1.2_r8*fall_speed_air_density_correction(i,k))

           ! adjust the ice fall velocity for smaller (r < 20 um) ice
           ! particles (blend over 8-20 um)
           irad = 1.5_r8 / cldice_size_dist_slope_per_m(i,k) * 1e6_r8
           ifrac = min(1._r8, max(0._r8, (irad - 18._r8) / 2._r8))

           if (ifrac .lt. 1._r8) then
              proc_rates%cldice_mass_weighted_fall_speed_m_per_s(i,k) = ifrac * proc_rates%cldice_mass_weighted_fall_speed_m_per_s(i,k) + &
                 (1._r8 - ifrac) * &
                 min(small_cldice_fall_speed_prefactor(i,k)*gamma_bj_plus4/(6._r8*cldice_size_dist_slope_per_m(i,k)**bj), &
                 1.2_r8*fall_speed_air_density_correction(i,k))
              proc_rates%cldice_mass_weighted_fall_speed_m_per_s(i,k)=proc_rates%cldice_mass_weighted_fall_speed_m_per_s(i,k)*micro_mg_vtrmi_factor

              cldice_mass_pressure_fall_speed_Pa_per_s(i,k)  = g*air_density_kg_per_m3(i,k)*proc_rates%cldice_mass_weighted_fall_speed_m_per_s(i,k)
              cldice_number_pressure_fall_speed_Pa_per_s(i,k) = ifrac * cldice_number_pressure_fall_speed_Pa_per_s(i,k) + &
                 (1._r8 - ifrac) * &
                 g*air_density_kg_per_m3(i,k)* &
                 min(small_cldice_fall_speed_prefactor(i,k)*gamma_bj_plus1/cldice_size_dist_slope_per_m(i,k)**bj,1.2_r8*fall_speed_air_density_correction(i,k))
           end if

           ! Fix ice fall speed following IFS microphysics
           if (ifs_sed) then
              cldice_mass_pressure_fall_speed_Pa_per_s(i,k)=g*air_density_kg_per_m3(i,k)*0.1_r8
              cldice_number_pressure_fall_speed_Pa_per_s(i,k)=g*air_density_kg_per_m3(i,k)*0.1_r8
           end if
        else
           cldice_mass_pressure_fall_speed_Pa_per_s(i,k) = 0._r8
           cldice_number_pressure_fall_speed_Pa_per_s(i,k)= 0._r8
        end if

        ! redefine dummy variables - sedimentation is calculated over grid-scale
        ! quantities to ensure conservation
        dumi(i,k) = (cldice_mass_mixing_ratio_kg_per_kg(i,k)+qitend(i,k)*deltat)
        dumni(i,k) = max((cldice_number_conc_per_kg(i,k)+nitend(i,k)*deltat),0._r8)
        if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present) async(RQUEUE)
  !$acc loop gang vector
  do i=1,mgncol
     !$acc loop seq
     do k=1,nlev
        if (rain_size_dist_slope_per_m(i,k).ge.qsmall) then
           qtmp = rain_size_dist_slope_per_m(i,k)**br
           ! 'final' values of number and mass weighted mean fallspeed for rain (m/s)
           rain_number_weighted_fall_speed_m_per_s(i,k) = min(rain_fall_speed_prefactor(i,k)*gamma_br_plus1/qtmp,9.1_r8*fall_speed_air_density_correction(i,k))
           rain_number_pressure_fall_speed_Pa_per_s(i,k) = g*air_density_kg_per_m3(i,k)*rain_number_weighted_fall_speed_m_per_s(i,k)
           proc_rates%rain_mass_weighted_fall_speed_m_per_s(i,k) = min(rain_fall_speed_prefactor(i,k)*gamma_br_plus4/(6._r8*qtmp),9.1_r8*fall_speed_air_density_correction(i,k))
           rain_mass_pressure_fall_speed_Pa_per_s(i,k) = g*air_density_kg_per_m3(i,k)*proc_rates%rain_mass_weighted_fall_speed_m_per_s(i,k)
        else
           rain_mass_pressure_fall_speed_Pa_per_s(i,k)=0._r8
           rain_number_pressure_fall_speed_Pa_per_s(i,k)=0._r8
        end if

        ! Fallspeed correction to ensure non-zero if rain in the column
        ! from updated Morrison (WRFv3.3) and P3 schemes
        ! If fallspeed exists at a higher level, apply it below to eliminate
        if (precip_fall_corr) then
           if (k.gt.2) then
              if (rain_mass_pressure_fall_speed_Pa_per_s(i,k).lt.1.e-10_r8) then
                 rain_mass_pressure_fall_speed_Pa_per_s(i,k)=rain_mass_pressure_fall_speed_Pa_per_s(i,k-1)
                 rain_number_pressure_fall_speed_Pa_per_s(i,k)=rain_number_pressure_fall_speed_Pa_per_s(i,k-1)
              end if
           end if
        end if

        ! redefine dummy variables - sedimentation is calculated over grid-scale
        ! quantities to ensure conservation
        dumr(i,k) = (rain_mass_mixing_ratio_kg_per_kg(i,k)+qrtend(i,k)*deltat)
        dumnr(i,k) = max((rain_number_conc_per_kg(i,k)+nrtend(i,k)*deltat),0._r8)
        if (dumr(i,k).lt.qsmall) dumnr(i,k)=0._r8
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present) async(SQUEUE)
  !$acc loop gang vector
  do i=1,mgncol
     !$acc loop seq
     do k=1,nlev
        if (snow_size_dist_slope_per_m(i,k).ge.qsmall) then
           qtmp = snow_size_dist_slope_per_m(i,k)**bs
           ! 'final' values of number and mass weighted mean fallspeed for snow (m/s)
           proc_rates%snow_mass_weighted_fall_speed_m_per_s(i,k) = min(snow_fall_speed_prefactor(i,k)*gamma_bs_plus4/(6._r8*qtmp),1.2_r8*fall_speed_air_density_correction(i,k))
           proc_rates%snow_mass_weighted_fall_speed_m_per_s(i,k) = proc_rates%snow_mass_weighted_fall_speed_m_per_s(i,k)*micro_mg_vtrms_factor

           snow_mass_pressure_fall_speed_Pa_per_s(i,k)  = g*air_density_kg_per_m3(i,k)*proc_rates%snow_mass_weighted_fall_speed_m_per_s(i,k)
           snow_number_weighted_fall_speed_m_per_s(i,k) = min(snow_fall_speed_prefactor(i,k)*gamma_bs_plus1/qtmp,1.2_r8*fall_speed_air_density_correction(i,k))
           snow_number_pressure_fall_speed_Pa_per_s(i,k) = g*air_density_kg_per_m3(i,k)*snow_number_weighted_fall_speed_m_per_s(i,k)
           ! Fix fallspeed for snow
           if (ifs_sed) then
              proc_rates%snow_mass_weighted_fall_speed_m_per_s(i,k) = 1._r8
              snow_number_weighted_fall_speed_m_per_s(i,k) = 1._r8
            end if
        else
           snow_mass_pressure_fall_speed_Pa_per_s(i,k)=0._r8
           snow_number_pressure_fall_speed_Pa_per_s(i,k)=0._r8
        end if

        if (precip_fall_corr) then
           if (k.gt.2) then
              if (snow_mass_pressure_fall_speed_Pa_per_s(i,k).lt.1.e-10_r8) then
                 snow_mass_pressure_fall_speed_Pa_per_s(i,k)=snow_mass_pressure_fall_speed_Pa_per_s(i,k-1)
                 snow_number_pressure_fall_speed_Pa_per_s(i,k)=snow_number_pressure_fall_speed_Pa_per_s(i,k-1)
              end if
           end if
        end if

        ! redefine dummy variables - sedimentation is calculated over grid-scale
        ! quantities to ensure conservation
        dums(i,k) = (snow_mass_mixing_ratio_kg_per_kg(i,k)+qstend(i,k)*deltat)
        dumns(i,k) = max((snow_number_conc_per_kg(i,k)+nstend(i,k)*deltat),0._r8)
        if (dums(i,k).lt.qsmall) dumns(i,k)=0._r8
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present) async(GQUEUE)
  !$acc loop gang vector
  do i=1,mgncol
     !$acc loop seq
     do k=1,nlev
        if (graupel_size_dist_slope_per_m(i,k).ge.qsmall) then
           qtmp = graupel_size_dist_slope_per_m(i,k)**bgtmp
           ! 'final' values of number and mass weighted mean fallspeed for graupel (m/s)
           proc_rates%graupel_mass_weighted_fall_speed_m_per_s(i,k) = min(graupel_fall_speed_prefactor(i,k)*gamma_bg_plus4/(6._r8*qtmp),20._r8*fall_speed_air_density_correction(i,k))
           graupel_mass_pressure_fall_speed_Pa_per_s(i,k) = g*air_density_kg_per_m3(i,k)*proc_rates%graupel_mass_weighted_fall_speed_m_per_s(i,k)
           graupel_number_weighted_fall_speed_m_per_s(i,k) = min(graupel_fall_speed_prefactor(i,k)*gamma_bg_plus1/qtmp,20._r8*fall_speed_air_density_correction(i,k))
           graupel_number_pressure_fall_speed_Pa_per_s(i,k) = g*air_density_kg_per_m3(i,k)*graupel_number_weighted_fall_speed_m_per_s(i,k)
        else
           graupel_mass_pressure_fall_speed_Pa_per_s(i,k)=0._r8
           graupel_number_pressure_fall_speed_Pa_per_s(i,k)=0._r8
        end if

        if (precip_fall_corr) then
           if (k.gt.2) then
              if (graupel_mass_pressure_fall_speed_Pa_per_s(i,k).lt.1.e-10_r8) then
                 graupel_mass_pressure_fall_speed_Pa_per_s(i,k)=graupel_mass_pressure_fall_speed_Pa_per_s(i,k-1)
                 graupel_number_pressure_fall_speed_Pa_per_s(i,k)=graupel_number_pressure_fall_speed_Pa_per_s(i,k-1)
              end if
           end if
        end if

        ! redefine dummy variables - sedimentation is calculated over grid-scale
        ! quantities to ensure conservation
        dumg(i,k) = (graupel_mass_mixing_ratio_kg_per_kg(i,k)+qgtend(i,k)*deltat)
        dumng(i,k) = max((graupel_number_conc_per_kg(i,k)+ngtend(i,k)*deltat),0._r8)
        if (dumg(i,k).lt.qsmall) dumng(i,k)=0._r8
     end do
  end do
  !$acc end parallel

! ----------------------------------------------
! Sedimentation
! ----------------------------------------------

if ( do_implicit_fall ) then

! Implicit Sedimentation calculation: from Guo et al, 2021, GFDL version.

  !$acc parallel vector_length(VLENS) default(present) async(LQUEUE)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        cldliq_mass_pressure_fall_speed_Pa_per_s(i,k)  = vfac_drop * cldliq_mass_pressure_fall_speed_Pa_per_s(i,k)/g/air_density_kg_per_m3(i,k)
        cldliq_number_pressure_fall_speed_Pa_per_s(i,k) = vfac_drop * cldliq_number_pressure_fall_speed_Pa_per_s(i,k)/g/air_density_kg_per_m3(i,k)
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present) async(IQUEUE)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        cldice_mass_pressure_fall_speed_Pa_per_s(i,k)  = vfac_ice  * cldice_mass_pressure_fall_speed_Pa_per_s(i,k)/g/air_density_kg_per_m3(i,k)
        cldice_number_pressure_fall_speed_Pa_per_s(i,k) = vfac_ice  * cldice_number_pressure_fall_speed_Pa_per_s(i,k)/g/air_density_kg_per_m3(i,k)
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present) async(RQUEUE)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        rain_mass_pressure_fall_speed_Pa_per_s(i,k)  = vfactor * rain_mass_pressure_fall_speed_Pa_per_s(i,k)/g/air_density_kg_per_m3(i,k)
        rain_number_pressure_fall_speed_Pa_per_s(i,k) = vfactor * rain_number_pressure_fall_speed_Pa_per_s(i,k)/g/air_density_kg_per_m3(i,k)
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present) async(SQUEUE)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        snow_mass_pressure_fall_speed_Pa_per_s(i,k)  = vfactor * snow_mass_pressure_fall_speed_Pa_per_s(i,k)/g/air_density_kg_per_m3(i,k)
        snow_number_pressure_fall_speed_Pa_per_s(i,k) = vfactor * snow_number_pressure_fall_speed_Pa_per_s(i,k)/g/air_density_kg_per_m3(i,k)
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present) async(GQUEUE)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        graupel_mass_pressure_fall_speed_Pa_per_s(i,k)  = vfactor * graupel_mass_pressure_fall_speed_Pa_per_s(i,k)/g/air_density_kg_per_m3(i,k)
        graupel_number_pressure_fall_speed_Pa_per_s(i,k) = vfactor * graupel_number_pressure_fall_speed_Pa_per_s(i,k)/g/air_density_kg_per_m3(i,k)
     end do
  end do
  !$acc end parallel

  ! cloud water mass sedimentation

  call Sedimentation_implicit(mgncol,nlev,deltat,zint,pdel,dumc,cldliq_mass_pressure_fall_speed_Pa_per_s,.FALSE.,qctend, &
                              LQUEUE,xflx=lflx,qxsedten=proc_rates%cldliq_mass_tend_sedimentation_kg_per_kg_per_s,prect=prect_l)

  ! cloud water number sedimentation
  call Sedimentation_implicit(mgncol,nlev,deltat,zint,pdel,dumnc,cldliq_number_pressure_fall_speed_Pa_per_s,.FALSE.,nctend, &
                              LQUEUE,qxsedten=proc_rates%cldliq_number_tend_sedimentation_per_kg_per_s)

  ! cloud ice mass sedimentation

   call Sedimentation_implicit(mgncol,nlev,deltat,zint,pdel,dumi,cldice_mass_pressure_fall_speed_Pa_per_s,.FALSE.,qitend, &
                               IQUEUE,xflx=iflx,qxsedten=proc_rates%cldice_mass_tend_sedimentation_kg_per_kg_per_s,prect=prect_i,preci=preci_i)

  ! cloud ice number sedimentation

  call Sedimentation_implicit(mgncol,nlev,deltat,zint,pdel,dumni,cldice_number_pressure_fall_speed_Pa_per_s,.FALSE.,nitend, &
                              IQUEUE,qxsedten=proc_rates%cldice_number_tend_sedimentation_per_kg_per_s)

  ! rain water mass sedimentation

  call Sedimentation_implicit(mgncol,nlev,deltat,zint,pdel,dumr,rain_mass_pressure_fall_speed_Pa_per_s,.TRUE.,qrtend, &
                              RQUEUE,xflx=rflx,qxsedten=proc_rates%rain_mass_tend_sedimentation_kg_per_kg_per_s,prect=prect_r)

  ! rain water number sedimentation

  call Sedimentation_implicit(mgncol,nlev,deltat,zint,pdel,dumnr,rain_number_pressure_fall_speed_Pa_per_s,.TRUE.,nrtend, &
                              RQUEUE,qxsedten=proc_rates%rain_number_tend_sedimentation_per_kg_per_s)

  ! snow water mass sedimentation

  call Sedimentation_implicit(mgncol,nlev,deltat,zint,pdel,dums,snow_mass_pressure_fall_speed_Pa_per_s,.TRUE.,qstend, &
                              SQUEUE,xflx=sflx,qxsedten=proc_rates%snow_mass_tend_sedimentation_kg_per_kg_per_s,prect=prect_s,preci=preci_s)

  ! snow water number sedimentation

  call Sedimentation_implicit(mgncol,nlev,deltat,zint,pdel,dumns,snow_number_pressure_fall_speed_Pa_per_s,.TRUE.,nstend, &
                              SQUEUE,qxsedten=proc_rates%snow_number_tend_sedimentation_per_kg_per_s)

  ! graupel mass sedimentation

   call Sedimentation_implicit(mgncol,nlev,deltat,zint,pdel,dumg,graupel_mass_pressure_fall_speed_Pa_per_s,.TRUE.,qgtend, &
                               GQUEUE,xflx=gflx,qxsedten=proc_rates%graupel_mass_tend_sedimentation_kg_per_kg_per_s,prect=prect_g,preci=preci_g)

  ! graupel number sedimentation

  call Sedimentation_implicit(mgncol,nlev,deltat,zint,pdel,dumng,graupel_number_pressure_fall_speed_Pa_per_s,.TRUE.,ngtend, &
                              GQUEUE,qxsedten=proc_rates%graupel_number_tend_sedimentation_per_kg_per_s)

else

! Explicit Sedimentation calculation

  !$acc parallel vector_length(VLENS) default(present) async(IQUEUE)
  !$acc loop gang vector
  do i = 1, mgncol
     nstep_i(i) = 1 + int( max( maxval( cldice_mass_pressure_fall_speed_Pa_per_s(i,:)*pdel_inv(i,:) ), maxval( cldice_number_pressure_fall_speed_Pa_per_s(i,:)*pdel_inv(i,:) ) ) * deltat )
     rnstep_i(i) = 1._r8/real(nstep_i(i))
  end do
  !$acc end parallel

  ! ice mass sediment
  call Sedimentation(mgncol,nlev,do_cldice,deltat,nstep_i,rnstep_i,cldice_mass_pressure_fall_speed_Pa_per_s,dumi,pdel_inv, &
                     qitend,IQUEUE,qxsedten=proc_rates%cldice_mass_tend_sedimentation_kg_per_kg_per_s,prect=prect_i,xflx=iflx,xxlx=xxls, &
                     qxsevap=proc_rates%cldice_mass_tend_sedimentation_sublimation_kg_per_kg_per_s,tlat=tlat_i,qvlat=qvlat_i,xcldm=ice_cloud_area_fraction,preci=preci_i)

  ! ice number sediment
  call Sedimentation(mgncol,nlev,do_cldice,deltat,nstep_i,rnstep_i,cldice_number_pressure_fall_speed_Pa_per_s,dumni,pdel_inv, &
                     nitend,IQUEUE,xcldm=ice_cloud_area_fraction,qxsedten=proc_rates%cldice_number_tend_sedimentation_per_kg_per_s)

  !$acc parallel vector_length(VLENS) default(present) async(LQUEUE)
  !$acc loop gang vector
  do i = 1, mgncol
     nstep_l(i) = 1 + int( max( maxval( cldliq_mass_pressure_fall_speed_Pa_per_s(i,:)*pdel_inv(i,:) ), maxval( cldliq_number_pressure_fall_speed_Pa_per_s(i,:)*pdel_inv(i,:) ) ) * deltat )
     rnstep_l(i) = 1._r8/real(nstep_l(i))
  end do
  !$acc end parallel

  ! liq mass sediment
  call Sedimentation(mgncol,nlev,.TRUE.,deltat,nstep_l,rnstep_l,cldliq_mass_pressure_fall_speed_Pa_per_s,dumc,pdel_inv, &
                     qctend,LQUEUE,qxsedten=proc_rates%cldliq_mass_tend_sedimentation_kg_per_kg_per_s,prect=prect_l,xflx=lflx,xxlx=xxlv, &
                     qxsevap=proc_rates%cldliq_mass_tend_sedimentation_evaporation_kg_per_kg_per_s,tlat=tlat_l,qvlat=qvlat_l,xcldm=liquid_cloud_area_fraction)

  ! liq number sediment
  call Sedimentation(mgncol,nlev,.TRUE.,deltat,nstep_l,rnstep_l,cldliq_number_pressure_fall_speed_Pa_per_s,dumnc,pdel_inv, &
                     nctend,LQUEUE,xcldm=liquid_cloud_area_fraction,qxsedten=proc_rates%cldliq_number_tend_sedimentation_per_kg_per_s)

  !$acc parallel vector_length(VLENS) default(present) async(RQUEUE)
  !$acc loop gang vector
  do i = 1, mgncol
     nstep_r(i) = 1 + int( max( maxval( rain_mass_pressure_fall_speed_Pa_per_s(i,:)*pdel_inv(i,:) ), maxval( rain_number_pressure_fall_speed_Pa_per_s(i,:)*pdel_inv(i,:) ) ) * deltat )
     rnstep_r(i) = 1._r8/real(nstep_r(i))
  end do
  !$acc end parallel

  ! rain mass sediment
  call Sedimentation(mgncol,nlev,.TRUE.,deltat,nstep_r,rnstep_r,rain_mass_pressure_fall_speed_Pa_per_s,dumr,pdel_inv, &
                     qrtend,RQUEUE,qxsedten=proc_rates%rain_mass_tend_sedimentation_kg_per_kg_per_s,prect=prect_r,xflx=rflx)

  ! rain number sediment
  call Sedimentation(mgncol,nlev,.TRUE.,deltat,nstep_r,rnstep_r,rain_number_pressure_fall_speed_Pa_per_s,dumnr,pdel_inv, &
                     nrtend,RQUEUE,qxsedten=proc_rates%rain_number_tend_sedimentation_per_kg_per_s)

  !$acc parallel vector_length(VLENS) default(present) async(SQUEUE)
  !$acc loop gang vector
  do i = 1, mgncol
     nstep_s(i) = 1 + int( max( maxval( snow_mass_pressure_fall_speed_Pa_per_s(i,:)*pdel_inv(i,:) ), maxval( snow_number_pressure_fall_speed_Pa_per_s(i,:)*pdel_inv(i,:) ) ) * deltat )
     rnstep_s(i) = 1._r8/real(nstep_s(i))
  end do
  !$acc end parallel

  ! snow mass sediment
  call Sedimentation(mgncol,nlev,.TRUE.,deltat,nstep_s,rnstep_s,snow_mass_pressure_fall_speed_Pa_per_s,dums,pdel_inv, &
                     qstend,SQUEUE,qxsedten=proc_rates%snow_mass_tend_sedimentation_kg_per_kg_per_s,prect=prect_s,xflx=sflx,preci=preci_s)

  ! snow number sediment
  call Sedimentation(mgncol,nlev,.TRUE.,deltat,nstep_s,rnstep_s,snow_number_pressure_fall_speed_Pa_per_s,dumns,pdel_inv, &
                     nstend,SQUEUE,qxsedten=proc_rates%snow_number_tend_sedimentation_per_kg_per_s)

  !$acc parallel vector_length(VLENS) default(present) async(GQUEUE)
  !$acc loop gang vector
  do i = 1, mgncol
     nstep_g(i) = 1 + int( max( maxval( graupel_mass_pressure_fall_speed_Pa_per_s(i,:)*pdel_inv(i,:) ), maxval( graupel_number_pressure_fall_speed_Pa_per_s(i,:)*pdel_inv(i,:) ) ) * deltat )
     rnstep_g(i) = 1._r8/real(nstep_g(i))
  end do
  !$acc end parallel

  ! graupel mass sediment
  call Sedimentation(mgncol,nlev,.TRUE.,deltat,nstep_g,rnstep_g,graupel_mass_pressure_fall_speed_Pa_per_s,dumg,pdel_inv, &
                     qgtend,GQUEUE,qxsedten=proc_rates%graupel_mass_tend_sedimentation_kg_per_kg_per_s,prect=prect_g,xflx=gflx,preci=preci_g)

  ! graupel number sediment
  call Sedimentation(mgncol,nlev,.TRUE.,deltat,nstep_g,rnstep_g,graupel_number_pressure_fall_speed_Pa_per_s,dumng,pdel_inv, &
                     ngtend,GQUEUE,qxsedten=proc_rates%graupel_number_tend_sedimentation_per_kg_per_s)

end if
! ----------------------------------------------
! End Sedimentation
! ----------------------------------------------

  ! sum up the changes due to sedimentation process for different hydrometeors

  !$acc parallel vector_length(VLENS) default(present) wait(IQUEUE,LQUEUE)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        tlat(i,k)  = tlat(i,k) + tlat_i(i,k) + tlat_l(i,k)
        qvlat(i,k) = qvlat(i,k) + qvlat_i(i,k) + qvlat_l(i,k)
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) wait(RQUEUE,SQUEUE,GQUEUE)
  !$acc loop gang vector
  do i=1,mgncol
     prect(i)  = prect(i) + prect_i(i) + prect_l(i) + prect_r(i) + prect_s(i) + prect_g(i)
     preci(i)  = preci(i) + preci_i(i) + preci_s(i) + preci_g(i)
  end do
  !$acc end parallel

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! get new update for variables that includes sedimentation tendency
  ! note : here dum variables are grid-average, NOT in-cloud

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        dumc(i,k)  = max(cldliq_mass_mixing_ratio_kg_per_kg(i,k)+qctend(i,k)*deltat,0._r8)
        dumi(i,k)  = max(cldice_mass_mixing_ratio_kg_per_kg(i,k)+qitend(i,k)*deltat,0._r8)
        dumnc(i,k) = max(cldliq_number_conc_per_kg(i,k)+nctend(i,k)*deltat,0._r8)
        dumni(i,k) = max(cldice_number_conc_per_kg(i,k)+nitend(i,k)*deltat,0._r8)

        dumr(i,k)  = max(rain_mass_mixing_ratio_kg_per_kg(i,k)+qrtend(i,k)*deltat,0._r8)
        dumnr(i,k) = max(rain_number_conc_per_kg(i,k)+nrtend(i,k)*deltat,0._r8)
        dums(i,k)  = max(snow_mass_mixing_ratio_kg_per_kg(i,k)+qstend(i,k)*deltat,0._r8)
        dumns(i,k) = max(snow_number_conc_per_kg(i,k)+nstend(i,k)*deltat,0._r8)
        dumg(i,k)  = max(graupel_mass_mixing_ratio_kg_per_kg(i,k)+qgtend(i,k)*deltat,0._r8)
        dumng(i,k) = max(graupel_number_conc_per_kg(i,k)+ngtend(i,k)*deltat,0._r8)

        ! switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(i,k)=ncnst/air_density_kg_per_m3(i,k)*liquid_cloud_area_fraction(i,k)
        end if

        ! switch for specification of cloud ice number
        if (nicons) then
           dumni(i,k)=ninst/air_density_kg_per_m3(i,k)*ice_cloud_area_fraction(i,k)
        end if

        ! switch for specification of graupel number
        if (ngcons) then
           dumng(i,k)=ngnst/air_density_kg_per_m3(i,k)*precipitation_area_fraction(i,k)
        end if

        ! switch for specification of constant snow number
        if (nscons) then
            dumns(i,k)=nsnst/air_density_kg_per_m3(i,k)
        end if

        ! switch for specification of constant rain number
        if (nrcons) then
            dumnr(i,k)=nrnst/air_density_kg_per_m3(i,k)
        end if

        if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
        if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8
        if (dumr(i,k).lt.qsmall) dumnr(i,k)=0._r8
        if (dums(i,k).lt.qsmall) dumns(i,k)=0._r8
        if (dumg(i,k).lt.qsmall) dumng(i,k)=0._r8

  ! calculate instantaneous processes (melting, homogeneous freezing)
  !====================================================================
  ! melting of snow at +2 C

        if (t(i,k)+tlat(i,k)/cpp*deltat > snowmelt) then
           if (dums(i,k) > 0._r8) then
              ! make sure melting snow doesn't reduce temperature below threshold
              dum = -xlf/cpp*dums(i,k)
              if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt. snowmelt) then
                 dum = (t(i,k)+tlat(i,k)/cpp*deltat-snowmelt)*cpp/xlf
                 dum = dum/dums(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              qstend(i,k)=qstend(i,k)-dum*dums(i,k)*rdeltat
              nstend(i,k)=nstend(i,k)-dum*dumns(i,k)*rdeltat
              qrtend(i,k)=qrtend(i,k)+dum*dums(i,k)*rdeltat
              nrtend(i,k)=nrtend(i,k)+dum*dumns(i,k)*rdeltat

              dum1=-xlf*dum*dums(i,k)*rdeltat
              tlat(i,k)=tlat(i,k)+dum1
              proc_rates%air_heating_rate_snow_melting_W_per_kg(i,k)=proc_rates%air_heating_rate_snow_melting_W_per_kg(i,k) + dum1

!STOPPED FIX FOR SNOW NUMBER
!ensure that snow... number does not go negative with constant number set
!necessary because dumng is updated above.
              if (nscons .and. ((snow_number_conc_per_kg(i,k)+nstend(i,k)*deltat) .lt. 0._r8)) then
                 nstend(i,k)=-snow_number_conc_per_kg(i,k)*rdeltat
              end if
           end if
        end if

  ! melting of graupel at +2 C

        if (t(i,k)+tlat(i,k)/cpp*deltat > snowmelt) then
           if (dumg(i,k) > 0._r8) then
              ! make sure melting graupel doesn't reduce temperature below threshold
              dum = -xlf/cpp*dumg(i,k)
              if (t(i,k)+tlat(i,k)/cpp*deltat+dum .lt. snowmelt) then
                 dum = (t(i,k)+tlat(i,k)/cpp*deltat-snowmelt)*cpp/xlf
                 dum = dum/dumg(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              qgtend(i,k)=qgtend(i,k)-dum*dumg(i,k)*rdeltat
              ngtend(i,k)=ngtend(i,k)-dum*dumng(i,k)*rdeltat
              qrtend(i,k)=qrtend(i,k)+dum*dumg(i,k)*rdeltat
              nrtend(i,k)=nrtend(i,k)+dum*dumng(i,k)*rdeltat

              dum1=-xlf*dum*dumg(i,k)*rdeltat
              tlat(i,k)=tlat(i,k)+dum1
              proc_rates%air_heating_rate_snow_melting_W_per_kg(i,k)=proc_rates%air_heating_rate_snow_melting_W_per_kg(i,k) + dum1

!ensure that graupel number does not go negative with constant number set
!necessary because dumng is updated above.
              if (ngcons .and. ((graupel_number_conc_per_kg(i,k)+ngtend(i,k)*deltat) .lt. 0._r8)) then
                 ngtend(i,k)=-graupel_number_conc_per_kg(i,k)*rdeltat
              end if
           end if
        end if
     end do
  end do
  !$acc end parallel

  ! get mean size of rain = 1/lamr, add frozen rain to either snow or cloud ice
  ! depending on mean rain size
  ! add to graupel if using that option....
  call size_dist_param_basic(mg_rain_props, dumr, dumnr, rain_size_dist_slope_per_m, mgncol, nlev)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        ! freezing of rain at -5 C
        if (t(i,k)+tlat(i,k)/cpp*deltat < rainfrze) then
           if (dumr(i,k) > 0._r8) then
              ! make sure freezing rain doesn't increase temperature above threshold
              dum = xlf/cpp*dumr(i,k)
              if (t(i,k)+tlat(i,k)/cpp*deltat+dum.gt.rainfrze) then
                 dum = -(t(i,k)+tlat(i,k)/cpp*deltat-rainfrze)*cpp/xlf
                 dum = dum/dumr(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              qrtend(i,k)=qrtend(i,k)-dum*dumr(i,k)*rdeltat
              nrtend(i,k)=nrtend(i,k)-dum*dumnr(i,k)*rdeltat

              if (rain_size_dist_slope_per_m(i,k) < 1._r8/Dcs) then
                 if (do_hail.or.do_graupel) then
                    qgtend(i,k)=qgtend(i,k)+dum*dumr(i,k)*rdeltat
                    ngtend(i,k)=ngtend(i,k)+dum*dumnr(i,k)*rdeltat
                 else
                    qstend(i,k)=qstend(i,k)+dum*dumr(i,k)*rdeltat
                    nstend(i,k)=nstend(i,k)+dum*dumnr(i,k)*rdeltat
                 end if
              else
                 qitend(i,k)=qitend(i,k)+dum*dumr(i,k)*rdeltat
                 nitend(i,k)=nitend(i,k)+dum*dumnr(i,k)*rdeltat
              end if

              ! heating tendency
              dum1 = xlf*dum*dumr(i,k)*rdeltat
              proc_rates%air_heating_rate_rain_freezing_W_per_kg(i,k)=proc_rates%air_heating_rate_rain_freezing_W_per_kg(i,k) + dum1
              tlat(i,k)=tlat(i,k)+dum1
           end if
        end if
      end do
   end do
   !$acc end parallel

   if (do_cldice) then
      !$acc parallel vector_length(VLENS) default(present)
      !$acc loop gang vector collapse(2)
      do k=1,nlev
        do i=1,mgncol
           if (t(i,k)+tlat(i,k)/cpp*deltat > tmelt) then
              if (dumi(i,k) > 0._r8) then
                 ! limit so that melting does not push temperature below freezing
                 !-----------------------------------------------------------------
                 dum = -dumi(i,k)*xlf/cpp
                 if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt.tmelt) then
                    dum = (t(i,k)+tlat(i,k)/cpp*deltat-tmelt)*cpp/xlf
                    dum = dum/dumi(i,k)
                    dum = max(0._r8,dum)
                    dum = min(1._r8,dum)
                 else
                    dum = 1._r8
                 end if

                 qctend(i,k)=qctend(i,k)+dum*dumi(i,k)*rdeltat

                 ! for output
                 proc_rates%cldice_mass_tend_melting_kg_per_kg_per_s(i,k)=dum*dumi(i,k)*rdeltat

                 ! assume melting ice produces droplet
                 ! mean volume radius of 8 micron

                 proc_rates%cldice_number_tend_melting_per_kg_per_s(i,k)=3._r8*dum*dumi(i,k)*rdeltat/ &
                      (4._r8*pi*5.12e-16_r8*rhow)
                 nctend(i,k)=nctend(i,k)+proc_rates%cldice_number_tend_melting_per_kg_per_s(i,k)

                 qitend(i,k)=((1._r8-dum)*dumi(i,k)-cldice_mass_mixing_ratio_kg_per_kg(i,k))*rdeltat
                 nitend(i,k)=((1._r8-dum)*dumni(i,k)-cldice_number_conc_per_kg(i,k))*rdeltat
                 tlat(i,k)=tlat(i,k)-xlf*dum*dumi(i,k)*rdeltat
              end if
           end if

           ! homogeneously freeze droplets at -40 C
           !-----------------------------------------------------------------

           if (t(i,k)+tlat(i,k)/cpp*deltat < 233.15_r8) then
              if (dumc(i,k) > 0._r8) then
                 ! limit so that freezing does not push temperature above threshold
                 dum = dumc(i,k)*xlf/cpp
                 if (t(i,k)+tlat(i,k)/cpp*deltat+dum.gt.233.15_r8) then
                    dum = -(t(i,k)+tlat(i,k)/cpp*deltat-233.15_r8)*cpp/xlf
                    dum = dum/dumc(i,k)
                    dum = max(0._r8,dum)
                    dum = min(1._r8,dum)
                 else
                    dum = 1._r8
                 end if

                 qitend(i,k)=qitend(i,k)+dum*dumc(i,k)*rdeltat
                 ! for output
                 proc_rates%cldliq_mass_tend_homogeneous_freezing_kg_per_kg_per_s(i,k)=dum*dumc(i,k)*rdeltat

                 ! assume 25 micron mean volume radius of homogeneously frozen droplets
                 ! consistent with size of detrained ice in stratiform.F90
                 proc_rates%cldliq_number_tend_homogeneous_freezing_per_kg_per_s(i,k)=dum*3._r8*dumc(i,k)/(4._r8*3.14_r8*micro_mg_homog_size**3._r8*500._r8)*rdeltat
                 nitend(i,k)=nitend(i,k)+proc_rates%cldliq_number_tend_homogeneous_freezing_per_kg_per_s(i,k)

                 qctend(i,k)=((1._r8-dum)*dumc(i,k)-cldliq_mass_mixing_ratio_kg_per_kg(i,k))*rdeltat
                 nctend(i,k)=((1._r8-dum)*dumnc(i,k)-cldliq_number_conc_per_kg(i,k))*rdeltat
                 tlat(i,k)=tlat(i,k)+xlf*dum*dumc(i,k)*rdeltat
              end if
           end if

           ! ice number limiter
           if (do_cldice .and. nitend(i,k).gt.0._r8.and.cldice_number_conc_per_kg(i,k)+nitend(i,k)*deltat.gt.micro_mg_max_nicons*ice_cloud_area_fraction(i,k)/air_density_kg_per_m3(i,k)) then
              nitend(i,k)=max(0._r8,(micro_mg_max_nicons*ice_cloud_area_fraction(i,k)/air_density_kg_per_m3(i,k)-cldice_number_conc_per_kg(i,k))/deltat)
           end if

     ! remove any excess over-saturation, which is possible due to non-linearity when adding
     ! together all microphysical processes
     !-----------------------------------------------------------------
     ! follow code similar to old CAM scheme

           dum_2D(i,k)=q(i,k)+qvlat(i,k)*deltat
           ttmpA(i,k)=t(i,k)+tlat(i,k)/cpp*deltat
        end do
     end do
     !$acc end parallel

     ! use rhw to allow ice supersaturation
     call qsat_water(ttmpA, p, esnA, qvnA, mgncol*nlev)

     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           if (dum_2D(i,k) > qvnA(i,k) .and. qvnA(i,k) > 0 .and. remove_supersat) then
              ! expression below is approximate since there may be ice deposition
              dum = (dum_2D(i,k)-qvnA(i,k))/(1._r8+xxlv_squared*qvnA(i,k)/(cpp*rv*ttmpA(i,k)**2))*rdeltat
              ! add to output cme
              cmeout(i,k) = cmeout(i,k)+dum
              ! now add to tendencies, partition between liquid and ice based on temperature
              if (ttmpA(i,k) > 268.15_r8) then
                 dum1=0.0_r8
                 ! now add to tendencies, partition between liquid and ice based on te
                 !-------------------------------------------------------
              else if (ttmpA(i,k) < 238.15_r8) then
                 dum1=1.0_r8
              else
                 dum1=(268.15_r8-ttmpA(i,k))/30._r8
              end if
              dum = (dum_2D(i,k)-qvnA(i,k))/(1._r8+(xxls*dum1+xxlv*(1._r8-dum1))**2 &
                    *qvnA(i,k)/(cpp*rv*ttmpA(i,k)**2))*rdeltat
              qctend(i,k)=qctend(i,k)+dum*(1._r8-dum1)
              ! for output
              proc_rates%cldliq_mass_tend_residual_condensation_kg_per_kg_per_s(i,k)=dum*(1._r8-dum1)
              qitend(i,k)=qitend(i,k)+dum*dum1
              proc_rates%cldice_mass_tend_residual_deposition_kg_per_kg_per_s(i,k)=dum*dum1
              qvlat(i,k)=qvlat(i,k)-dum
              ! for output
              proc_rates%vapor_mass_tend_residual_condensation_kg_per_kg_per_s(i,k)=-dum
              tlat(i,k)=tlat(i,k)+dum*(1._r8-dum1)*xxlv+dum*dum1*xxls
           end if
        end do
     end do
     !$acc end parallel
  end if

  ! calculate effective radius for pass to radiation code
  !=========================================================
  ! if no cloud water, default value is 10 micron for droplets,
  ! 25 micron for cloud ice

  ! update cloud variables after instantaneous processes to get effective radius
  ! variables are in-cloud to calculate size dist parameters

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        dumc(i,k) = max(cldliq_mass_mixing_ratio_kg_per_kg(i,k)+qctend(i,k)*deltat,0._r8)/liquid_cloud_area_fraction(i,k)
        dumi(i,k) = max(cldice_mass_mixing_ratio_kg_per_kg(i,k)+qitend(i,k)*deltat,0._r8)/ice_cloud_area_fraction(i,k)
        dumnc(i,k) = max(cldliq_number_conc_per_kg(i,k)+nctend(i,k)*deltat,0._r8)/liquid_cloud_area_fraction(i,k)
        dumni(i,k) = max(cldice_number_conc_per_kg(i,k)+nitend(i,k)*deltat,0._r8)/ice_cloud_area_fraction(i,k)

        dumr(i,k) = max(rain_mass_mixing_ratio_kg_per_kg(i,k)+qrtend(i,k)*deltat,0._r8)/precipitation_area_fraction(i,k)
        dumnr(i,k) = max(rain_number_conc_per_kg(i,k)+nrtend(i,k)*deltat,0._r8)/precipitation_area_fraction(i,k)
        dums(i,k) = max(snow_mass_mixing_ratio_kg_per_kg(i,k)+qstend(i,k)*deltat,0._r8)/precipitation_area_fraction(i,k)
        dumns(i,k) = max(snow_number_conc_per_kg(i,k)+nstend(i,k)*deltat,0._r8)/precipitation_area_fraction(i,k)
        dumg(i,k) = max(graupel_mass_mixing_ratio_kg_per_kg(i,k)+qgtend(i,k)*deltat,0._r8)
        dumng(i,k) = max(graupel_number_conc_per_kg(i,k)+ngtend(i,k)*deltat,0._r8)

        ! switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(i,k)=ncnst/air_density_kg_per_m3(i,k)
        end if

        ! switch for specification of cloud ice number
        if (nicons) then
           dumni(i,k)=ninst/air_density_kg_per_m3(i,k)
        end if

        ! switch for specification of graupel number
        if (ngcons) then
           dumng(i,k)=ngnst/air_density_kg_per_m3(i,k)*precipitation_area_fraction(i,k)
        end if

        ! switch for specification of constant snow number
        if (nscons) then
            dumns(i,k)=nsnst/air_density_kg_per_m3(i,k)
        end if

        ! switch for specification of constant rain number
        if (nrcons) then
            dumnr(i,k)=nrnst/air_density_kg_per_m3(i,k)
        end if

        ! limit in-cloud mixing ratio to reasonable value of 5 g kg-1
        dumc(i,k)=min(dumc(i,k),5.e-3_r8)
        dumi(i,k)=min(dumi(i,k),5.e-3_r8)
        ! limit in-precip mixing ratios
        dumr(i,k)=min(dumr(i,k),10.e-3_r8)
        dums(i,k)=min(dums(i,k),10.e-3_r8)
        dumg(i,k)=min(dumg(i,k),10.e-3_r8)
     end do
  end do
  !$acc end parallel

  ! cloud ice effective radius
  !-----------------------------------------------------------------
  if (do_cldice) then
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           dum_2D(i,k) = dumni(i,k)
        end do
     end do
     !$acc end parallel

     call size_dist_param_basic(mg_ice_props, dumi, dumni, cldice_size_dist_slope_per_m, mgncol, nlev, n0=dumni0A2D)

     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           if (dumi(i,k).ge.qsmall) then
              if (dumni(i,k) /=dum_2D(i,k)) then
                 ! adjust number conc if needed to keep mean size in reasonable range
                 nitend(i,k)=(dumni(i,k)*ice_cloud_area_fraction(i,k)-cldice_number_conc_per_kg(i,k))*rdeltat
              end if
              effi(i,k)   = 1.5_r8/cldice_size_dist_slope_per_m(i,k)*1.e6_r8
              effi(i,k)   = effi(i,k)*micro_mg_effi_factor

              sadice(i,k) = 2._r8*pi*(cldice_size_dist_slope_per_m(i,k)**(-3))*dumni0A2D(i,k)*air_density_kg_per_m3(i,k)*1.e-2_r8  ! m2/m3 -> cm2/cm3
           else
              effi(i,k)   = 25._r8
              effi(i,k)   = effi(i,k)*micro_mg_effi_factor

              sadice(i,k) = 0._r8
           end if
           ! ice effective diameter for david mitchell's optics
           deffi(i,k)=effi(i,k)*rhoi/rhows*2._r8
        end do
     end do
     !$acc end parallel
  else
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           ! NOTE: If CARMA is doing the ice microphysics, then the ice effective
           ! radius has already been determined from the size distribution.
           effi(i,k)   = re_ice(i,k) * 1.e6_r8      ! m -> um
           effi(i,k)   = effi(i,k)*micro_mg_effi_factor

           deffi(i,k)  = effi(i,k) * 2._r8
           sadice(i,k) = 4._r8*pi*(effi(i,k)**2)*cldice_number_conc_per_kg(i,k)*air_density_kg_per_m3(i,k)*1e-2_r8
        end do
     end do
     !$acc end parallel
  end if

  ! cloud droplet effective radius
  !-----------------------------------------------------------------

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        dum_2D(i,k) = dumnc(i,k)
     end do
  end do
  !$acc end parallel

  call size_dist_param_liq(mg_liq_props, dumc, dumnc, air_density_kg_per_m3, cldliq_size_dist_shape, cldliq_size_dist_slope_per_m, mgncol, nlev)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        if (dumc(i,k).ge.qsmall) then
           ! switch for specification of droplet and crystal number
           if (nccons) then
              ! make sure nc is consistence with the constant N by adjusting tendency, need
              ! to multiply by cloud fraction
              ! note that nctend may be further adjusted below if mean droplet size is
              ! out of bounds
              nctend(i,k)=(ncnst/air_density_kg_per_m3(i,k)*liquid_cloud_area_fraction(i,k)-cldliq_number_conc_per_kg(i,k))*rdeltat
           end if
           if (dum_2D(i,k) /= dumnc(i,k)) then
              ! adjust number conc if needed to keep mean size in reasonable range
              nctend(i,k)=(dumnc(i,k)*liquid_cloud_area_fraction(i,k)-cldliq_number_conc_per_kg(i,k))*rdeltat
           end if

           effc(i,k) = (cldliq_size_dist_shape(i,k)+3._r8)/cldliq_size_dist_slope_per_m(i,k)/2._r8*1.e6_r8
           !assign output fields for shape here
           lamcrad(i,k)=cldliq_size_dist_slope_per_m(i,k)
           pgamrad(i,k)=cldliq_size_dist_shape(i,k)

           ! recalculate effective radius for constant number, in order to separate
           ! first and second indirect effects
           !======================================
           ! assume constant number of 10^8 kg-1
           dumnc(i,k)=1.e8_r8
        end if
     end do
  end do
  !$acc end parallel

  ! Pass in "false" adjust flag to prevent number from being changed within
  ! size distribution subroutine.
  call size_dist_param_liq(mg_liq_props, dumc, dumnc, air_density_kg_per_m3, cldliq_size_dist_shape, cldliq_size_dist_slope_per_m, mgncol, nlev)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k =1,nlev
     do i=1,mgncol
        if (dumc(i,k).ge.qsmall) then
           effc_fn(i,k) = (cldliq_size_dist_shape(i,k)+3._r8)/cldliq_size_dist_slope_per_m(i,k)/2._r8*1.e6_r8
        else
           effc(i,k) = 10._r8
           lamcrad(i,k)=0._r8
           pgamrad(i,k)=0._r8
           effc_fn(i,k) = 10._r8
        end if

        ! recalculate 'final' rain size distribution parameters
        ! to ensure that rain size is in bounds, adjust rain number if needed
        dum_2D(i,k) = dumnr(i,k)
     end do
  end do
  !$acc end parallel

  call size_dist_param_basic(mg_rain_props, dumr, dumnr, rain_size_dist_slope_per_m, mgncol, nlev, n0=rain_size_dist_intercept_per_m4)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        if (dumr(i,k).ge.qsmall) then
           if (dum_2D(i,k) /= dumnr(i,k)) then
              ! adjust number conc if needed to keep mean size in reasonable range
              nrtend(i,k)=(dumnr(i,k)*precipitation_area_fraction(i,k)-rain_number_conc_per_kg(i,k))*rdeltat
           end if

        end if

        ! recalculate 'final' snow size distribution parameters
        ! to ensure that snow size is in bounds, adjust snow number if needed
        dum_2D(i,k) = dumns(i,k)
     end do
  end do
  !$acc end parallel

  call size_dist_param_basic(mg_snow_props, dums, dumns, snow_size_dist_slope_per_m, mgncol, nlev, n0=dumns0A2D)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        if (dums(i,k).ge.qsmall) then

           if (dum_2D(i,k) /= dumns(i,k)) then
              ! adjust number conc if needed to keep mean size in reasonable range
              nstend(i,k)=(dumns(i,k)*precipitation_area_fraction(i,k)-snow_number_conc_per_kg(i,k))*rdeltat
           end if

           sadsnow(i,k) = 2._r8*pi*(snow_size_dist_slope_per_m(i,k)**(-3))*dumns0A2D(i,k)*air_density_kg_per_m3(i,k)*1.e-2_r8  ! m2/m3 -> cm2/cm3

        end if

        ! recalculate 'final' graupel size distribution parameters
        ! to ensure that  size is in bounds, addjust number if needed
        dum_2D(i,k) = dumng(i,k)
     end do
  end do
  !$acc end parallel

  if (do_hail) then
     call size_dist_param_basic(mg_hail_props, dumg, dumng, graupel_size_dist_slope_per_m, mgncol, nlev)
  end if
  if (do_graupel) then
     call size_dist_param_basic(mg_graupel_props, dumg, dumng, graupel_size_dist_slope_per_m, mgncol, nlev)
  end if

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        if (dumg(i,k).ge.qsmall) then
           if (dum_2D(i,k) /= dumng(i,k)) then
              ! adjust number conc if needed to keep mean size in reasonable range
              ngtend(i,k)=(dumng(i,k)*precipitation_area_fraction(i,k)-graupel_number_conc_per_kg(i,k))*rdeltat
           end if
        end if

        ! if updated q (after microphysics) is zero, then ensure updated n is also zero
        !=================================================================================
        if (cldliq_mass_mixing_ratio_kg_per_kg(i,k)+qctend(i,k)*deltat.lt.qsmall) nctend(i,k)=-cldliq_number_conc_per_kg(i,k)*rdeltat
        if (do_cldice .and. cldice_mass_mixing_ratio_kg_per_kg(i,k)+qitend(i,k)*deltat.lt.qsmall) nitend(i,k)=-cldice_number_conc_per_kg(i,k)*rdeltat
        if (rain_mass_mixing_ratio_kg_per_kg(i,k)+qrtend(i,k)*deltat.lt.qsmall) nrtend(i,k)=-rain_number_conc_per_kg(i,k)*rdeltat
        if (snow_mass_mixing_ratio_kg_per_kg(i,k)+qstend(i,k)*deltat.lt.qsmall) nstend(i,k)=-snow_number_conc_per_kg(i,k)*rdeltat
        if (graupel_mass_mixing_ratio_kg_per_kg(i,k)+qgtend(i,k)*deltat.lt.qsmall) ngtend(i,k)=-graupel_number_conc_per_kg(i,k)*rdeltat

  ! DO STUFF FOR OUTPUT:
  !==================================================
  ! qc and qi are only used for output calculations past here,
  ! so add qctend and qitend back in one more time

        cldliq_mass_mixing_ratio_kg_per_kg(i,k) = cldliq_mass_mixing_ratio_kg_per_kg(i,k) + qctend(i,k)*deltat
        cldice_mass_mixing_ratio_kg_per_kg(i,k) = cldice_mass_mixing_ratio_kg_per_kg(i,k) + qitend(i,k)*deltat

  ! averaging for snow and rain number and diameter
  !--------------------------------------------------
  ! drout2/dsout2:
  ! diameter of rain and snow
  ! dsout:
  ! scaled diameter of snow (passed to radiation in CAM)
  ! reff_rain/reff_snow:
  ! calculate effective radius of rain and snow in microns for COSP using Eq. 9 of COSP v1.3 manual

        ! avoid divide by zero in avg_diameter_vec
        if (nrout(i,k) .eq. 0._r8) nrout(i,k)=1.e-34_r8
     end do
  end do
  !$acc end parallel

  ! The avg_diameter_vec call does the actual calculation; other diameter
  ! outputs are just drout2 times constants.
  call avg_diameter_vec(qrout,nrout,air_density_kg_per_m3,rhow,drout2,mgncol*nlev)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        if (qrout(i,k) .gt. 1.e-7_r8 .and. nrout(i,k) .gt. 0._r8) then
           qrout2(i,k) = qrout(i,k) * precipitation_area_fraction(i,k)
           nrout2(i,k) = nrout(i,k) * precipitation_area_fraction(i,k)
           freqr(i,k) = precipitation_area_fraction(i,k)
           reff_rain(i,k)=1.5_r8*drout2(i,k)*1.e6_r8
        else
           qrout2(i,k) = 0._r8
           nrout2(i,k) = 0._r8
           drout2(i,k) = 0._r8
           freqr(i,k) = 0._r8
           reff_rain(i,k) = 0._r8
        end if

        ! avoid divide by zero in avg_diameter_vec
        if (nsout(i,k) .eq. 0._r8) nsout(i,k) = 1.e-34_r8
     end do
  end do
  !$acc end parallel

  ! The avg_diameter_vec call does the actual calculation; other diameter
  ! outputs are just dsout2 times constants.
  call avg_diameter_vec(qsout, nsout, air_density_kg_per_m3, rhosn,dsout2,mgncol*nlev)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        if (qsout(i,k) .gt. 1.e-7_r8 .and. nsout(i,k) .gt. 0._r8) then
           qsout2(i,k) = qsout(i,k) * precipitation_area_fraction(i,k)
           nsout2(i,k) = nsout(i,k) * precipitation_area_fraction(i,k)
           freqs(i,k) = precipitation_area_fraction(i,k)
           dsout(i,k)=3._r8*rhosn/rhows*dsout2(i,k)
           reff_snow(i,k)=1.5_r8*dsout2(i,k)*1.e6_r8
        else
           dsout(i,k)  = 0._r8
           qsout2(i,k) = 0._r8
           nsout2(i,k) = 0._r8
           dsout2(i,k) = 0._r8
           freqs(i,k)  = 0._r8
           reff_snow(i,k)=0._r8
        end if

        ! avoid divide by zero in avg_diameter_vec
        if (ngout(i,k) .eq. 0._r8) ngout(i,k) = 1.e-34_r8
     end do
  end do
  !$acc end parallel

  ! The avg_diameter_vec call does the actual calculation; other diameter
  ! outputs are just dgout2 times constants.
  if (do_hail .or. do_graupel) then
     call avg_diameter_vec(qgout, ngout, air_density_kg_per_m3, rhogtmp, dgout2, mgncol*nlev)
  else
     ! need this if statement for MG2, where rhogtmp = 0

     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector collapse(2)
     do k=1,nlev
        do i=1,mgncol
           dgout2(i,k) = 0._r8
        end do
     end do
     !$acc end parallel
  end if

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        if (qgout(i,k) .gt. 1.e-7_r8 .and. ngout(i,k) .gt. 0._r8) then
           qgout2(i,k) = qgout(i,k) * precipitation_area_fraction(i,k)
           ngout2(i,k) = ngout(i,k) * precipitation_area_fraction(i,k)
           freqg(i,k) = precipitation_area_fraction(i,k)
           dgout(i,k)=3._r8*rhogtmp/rhows*dgout2(i,k)
           reff_grau(i,k)=1.5_r8*dgout2(i,k)*1.e6_r8
        else
           dgout(i,k)  = 0._r8
           qgout2(i,k) = 0._r8
           ngout2(i,k) = 0._r8
           dgout2(i,k) = 0._r8
           freqg(i,k)  = 0._r8
           reff_grau(i,k)=0._r8
        end if
     end do
  end do
  !$acc end parallel

  ! analytic radar reflectivity
  !--------------------------------------------------
  ! formulas from Matthew Shupe, NOAA/CERES
  ! *****note: radar reflectivity is local (in-precip average)
  ! units of mm^6/m^3

  ! Min rain rate of 0.1 mm/hr
  rthrsh=0.0001_r8/3600._r8

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        if (cldliq_mass_mixing_ratio_kg_per_kg(i,k).ge.qsmall .and. (cldliq_number_conc_per_kg(i,k)+nctend(i,k)*deltat).gt.10._r8) then
           dum=(cldliq_mass_mixing_ratio_kg_per_kg(i,k)/liquid_cloud_area_fraction(i,k)*air_density_kg_per_m3(i,k)*1000._r8)**2 &
                /(0.109_r8*(cldliq_number_conc_per_kg(i,k)+nctend(i,k)*deltat)/liquid_cloud_area_fraction(i,k)*air_density_kg_per_m3(i,k)/1.e6_r8)*liquid_cloud_area_fraction(i,k)/precipitation_area_fraction(i,k)
        else
           dum=0._r8
        end if
        if (cldice_mass_mixing_ratio_kg_per_kg(i,k).ge.qsmall) then
           dum1=(cldice_mass_mixing_ratio_kg_per_kg(i,k)*air_density_kg_per_m3(i,k)/ice_cloud_area_fraction(i,k)*1000._r8/0.1_r8)**(1._r8/0.63_r8)*ice_cloud_area_fraction(i,k)/precipitation_area_fraction(i,k)
        else
           dum1=0._r8
        end if
        if (qsout(i,k).ge.qsmall) then
           dum1=dum1+(qsout(i,k)*air_density_kg_per_m3(i,k)*1000._r8/0.1_r8)**(1._r8/0.63_r8)
        end if
        refl(i,k)=dum+dum1

        ! add rain to reflectivity (rain rate in mm/hr)
        ! reflectivity (dum) is in DBz
        ! New version Aircraft cloud values
        !Z=a*R^b (R in mm/hr) from Comstock et al 2004

        if (rflx(i,k+1).ge.rthrsh) then
           dum=32._r8*(rflx(i,k+1)*3600._r8)**1.4_r8
        else
           ! don't include rain rate in R calculation for values less than 0.001 mm/hr
           dum=0._r8
        end if

        ! add to refl
        refl(i,k)=refl(i,k)+dum

        !output reflectivity in Z.
        areflz(i,k)=refl(i,k) * precipitation_area_fraction(i,k)

        ! convert back to DBz
        if (refl(i,k).gt.minrefl) then
           refl(i,k)=10._r8*dlog10(refl(i,k))
        else
           refl(i,k)=-9999._r8
        end if

        !set averaging flag
        if (refl(i,k).gt.mindbz) then
           arefl(i,k)=refl(i,k) * precipitation_area_fraction(i,k)
           frefl(i,k)=precipitation_area_fraction(i,k)
        else
           arefl(i,k)=0._r8
           areflz(i,k)=0._r8
           frefl(i,k)=0._r8
        end if

        ! bound cloudsat reflectivity
        csrfl(i,k)=min(csmax,refl(i,k))

        !set averaging flag
        if (csrfl(i,k).gt.csmin) then
           acsrfl(i,k)=refl(i,k) * precipitation_area_fraction(i,k)
           fcsrfl(i,k)=precipitation_area_fraction(i,k)
        else
           acsrfl(i,k)=0._r8
           fcsrfl(i,k)=0._r8
        end if
     end do
  end do
  !$acc end parallel

  ! 10cm analytic radar reflectivity (rain radar)
  !--------------------------------------------------
  ! Formula from Hugh Morrison
  ! Ice dielectric correction from Smith 1984, Equation  10 and Snow correction from Smith 1984 Equation 14
  ! Smith, Paul L. “Equivalent Radar Reflectivity Factors for Snow and Ice Particles.
  !                ” Journal of Climate and Applied Meteorology 23, no. 8 (1984): 1258–60.
  !                DOI:  10.1175/1520-0450(1984)023<1258:ERRFFS>2.0.CO;2

  ! *****note: radar reflectivity is local (in-precip average)
  ! units of mm^6/m^3

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol

        dum1  = minrefl10
        dum2  = minrefl10
        dum3  = minrefl10
        dum4  = minrefl10
        dum   = minrefl10

!     Rain
        if (rain_size_dist_slope_per_m(i,k) > 0._r8) then
           dum1 = air_density_kg_per_m3(i,k)*rain_size_dist_intercept_per_m4(i,k)*720._r8/rain_size_dist_slope_per_m(i,k)**3/rain_size_dist_slope_per_m(i,k)**3/rain_size_dist_slope_per_m(i,k)
           dum1 = max(dum1,minrefl10)
        end if

!     Ice
        !  Add diaelectric factor from Smith 1984 equation 10
        if (cldice_size_dist_slope_per_m(i,k) > 0._r8) then
           dum2= air_density_kg_per_m3(i,k)*(0.176_r8/0.93_r8) * 720._r8*dumni0A2D(i,k)*(rhoi/900._r8)**2/cldice_size_dist_slope_per_m(i,k)**7
           dum2 = max(dum2,minrefl10)
        endif

!     Snow
        if (snow_size_dist_slope_per_m(i,k) > 0._r8) then
           dum3= air_density_kg_per_m3(i,k)*(0.176_r8/0.93_r8) * 720._r8*dumns0A2D(i,k)*(rhosn/900._r8)**2/snow_size_dist_slope_per_m(i,k)**7._r8
           dum3 = max(dum3,minrefl10)
        endif

!     Graupel
        if (do_hail .or. do_graupel .and. graupel_size_dist_slope_per_m(i,k) > 0._r8) then
          dum4= air_density_kg_per_m3(i,k)*(0.176_r8/0.93_r8) * 720._r8*graupel_size_dist_intercept_per_m4(i,k)*(rhogtmp/900._r8)**2/graupel_size_dist_slope_per_m(i,k)**7._r8
          dum4 =max(dum4,minrefl10)
        end if

        reflz10cm(i,k) = (dum1+dum2+dum3+dum4) * precipitation_area_fraction(i,k)

  ! Convert to dBz....

        dum = reflz10cm(i,k)*1.e18_r8
        refl10cm(i,k) = 10._r8*dlog10(dum)

  !redefine fice here....

        dum_2D(i,k) = qsout(i,k) + qrout(i,k) + cldliq_mass_mixing_ratio_kg_per_kg(i,k) + cldice_mass_mixing_ratio_kg_per_kg(i,k)
        dumi(i,k)   = qsout(i,k) + cldice_mass_mixing_ratio_kg_per_kg(i,k)
        if (dumi(i,k) .gt. qsmall .and. dum_2D(i,k) .gt. qsmall) then
           nfice(i,k) = min(dumi(i,k)/dum_2D(i,k),1._r8)
        else
           nfice(i,k) = 0._r8
        end if

     end do
  end do
  !$acc end parallel

  !$acc end data

end subroutine micro_pumas_tend

!========================================================================
!OUTPUT CALCULATIONS
!========================================================================

subroutine calc_rercld(rain_size_dist_slope_per_m, rain_size_dist_intercept_per_m4, cldliq_size_dist_slope_per_m, cldliq_size_dist_shape, inprecip_rain_mass_mixing_ratio_kg_per_kg, incloud_cldliq_mass_mixing_ratio_kg_per_kg, incloud_cldliq_number_conc_per_kg, rercld, vlen)
  integer,                   intent(in) :: vlen
  real(r8), dimension(vlen), intent(in) :: rain_size_dist_slope_per_m          ! rain size parameter (slope)
  real(r8), dimension(vlen), intent(in) :: rain_size_dist_intercept_per_m4           ! rain size parameter (intercept)
  real(r8), dimension(vlen), intent(in) :: cldliq_size_dist_slope_per_m          ! size distribution parameter (slope)
  real(r8), dimension(vlen), intent(in) :: cldliq_size_dist_shape          ! droplet size parameter
  real(r8), dimension(vlen), intent(in) :: inprecip_rain_mass_mixing_ratio_kg_per_kg          ! in-cloud rain mass mixing ratio
  real(r8), dimension(vlen), intent(in) :: incloud_cldliq_mass_mixing_ratio_kg_per_kg          ! in-cloud cloud liquid
  real(r8), dimension(vlen), intent(in) :: incloud_cldliq_number_conc_per_kg          ! in-cloud droplet number concentration

  real(r8), dimension(vlen), intent(inout) :: rercld     ! effective radius calculation for rain + cloud

  ! combined size of precip & cloud drops
  real(r8) :: Atmp,tmp(vlen), pgamp1(vlen)

  integer :: i

  !$acc data create (tmp,pgamp1)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     pgamp1(i) = cldliq_size_dist_shape(i)+1._r8
  end do
  !$acc end parallel

  call rising_factorial(pgamp1, 2, tmp, vlen)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     ! Rain drops
     if (rain_size_dist_slope_per_m(i) > 0._r8) then
        Atmp = rain_size_dist_intercept_per_m4(i) * pi / (2._r8 * rain_size_dist_slope_per_m(i)**3._r8)
     else
        Atmp = 0._r8
     end if
     ! Add cloud drops
     if (cldliq_size_dist_slope_per_m(i) > 0._r8) then
        Atmp = Atmp + &
             incloud_cldliq_number_conc_per_kg(i) * pi * tmp(i) / (4._r8 * cldliq_size_dist_slope_per_m(i)**2._r8)
     end if
     if (Atmp > 0._r8) then
        rercld(i) = rercld(i) + 3._r8 *(inprecip_rain_mass_mixing_ratio_kg_per_kg(i) + incloud_cldliq_mass_mixing_ratio_kg_per_kg(i)) / (4._r8 * rhow * Atmp)
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine calc_rercld

!========================================================================
!2020-09-15: Follow John Dennis's version to generate a new interface
!            to update tendency in the sedimentation loop;
!2021-10-19: Separate the mass and ice sediment for each class;
!========================================================================
subroutine Sedimentation(mgncol,nlev,do_cldice,deltat,nstep,rnstep,fx,dumx,pdel_inv, &
                         xxtend,queue,qxsedten,prect,xflx,xxlx,qxsevap,tlat,qvlat,xcldm,preci)

   integer, intent(in)               :: mgncol,nlev
   logical, intent(in)               :: do_cldice
   real(r8),intent(in)               :: deltat
   integer, intent(in)               :: nstep(mgncol)
   real(r8), intent(in)              :: rnstep(mgncol)
   real(r8), intent(in)              :: fx(mgncol,nlev)
   real(r8), intent(inout)           :: dumx(mgncol,nlev)
   real(r8), intent(in)              :: pdel_inv(mgncol,nlev)
   real(r8), intent(inout)           :: xxtend(mgncol,nlev)
   integer, intent(in)               :: queue
   real(r8), intent(inout), optional :: qxsedten(mgncol,nlev)
   real(r8), intent(inout), optional :: prect(mgncol)
   real(r8), intent(inout), optional :: xflx(mgncol,nlev+1)
   real(r8), intent(in)   , optional :: xxlx
   real(r8), intent(inout), optional :: qxsevap(mgncol,nlev)
   real(r8), intent(in)   , optional :: xcldm(mgncol,nlev)
   real(r8), intent(inout), optional :: tlat(mgncol,nlev)
   real(r8), intent(inout), optional :: qvlat(mgncol,nlev)
   real(r8), intent(inout), optional :: preci(mgncol)

   ! local variables
   integer  :: i,k,n,nstepmax
   real(r8) :: faltndx,rnstepmax,faltndqxe
   real(r8) :: dum1(mgncol,nlev),faloutx(mgncol,0:nlev)
   logical  :: present_tlat, present_qvlat, present_xcldm, present_qxsevap, &
               present_prect, present_preci, present_qxsedten, present_xflx

   present_tlat     = present(tlat)
   present_qvlat    = present(qvlat)
   present_xcldm    = present(xcldm)
   present_qxsevap  = present(qxsevap)
   present_preci    = present(preci)
   present_prect    = present(prect)
   present_qxsedten = present(qxsedten)
   present_xflx     = present(xflx)

   ! loop over sedimentation sub-time step to ensure stability
   !==============================================================

   !$acc enter data create (faloutx,dum1) async(queue)

   !$acc parallel vector_length(VLENS) default(present) async(queue)
   !$acc loop gang vector
   do i = 1,mgncol
      nstepmax = nstep(i)
      rnstepmax = rnstep(i)

      dum1(i,1) = 0._r8
      if (present_xcldm) then
         do k = 2,nlev
            dum1(i,k) = xcldm(i,k)/xcldm(i,k-1)
            dum1(i,k) = min(dum1(i,k),1._r8)
         end do
      else
         do k=2,nlev
            dum1(i,k) = 1._r8
         end do
      end if

      !$acc loop seq
      do n = 1,nstepmax
         faloutx(i,0)  = 0._r8
         if (do_cldice) then
            do k=1,nlev
               faloutx(i,k)  = fx(i,k)  * dumx(i,k)
            end do
         else
            do k=1,nlev
               faloutx(i,k)  = 0._r8
            end do
         end if

         do k = 1,nlev
            ! for cloud liquid and ice, if cloud fraction increases with height
            ! then add flux from above to both vapor and cloud water of current level
            ! this means that flux entering clear portion of cell from above evaporates
            ! instantly
            ! note: this is not an issue with precip, since we assume max overlap
            faltndx = (faloutx(i,k) - dum1(i,k) * faloutx(i,k-1)) * pdel_inv(i,k)
            ! add fallout terms to eulerian tendencies
            xxtend(i,k) = xxtend(i,k) - faltndx * rnstepmax
            ! sedimentation tendency for output
            if (present_qxsedten) qxsedten(i,k) = qxsedten(i,k)-faltndx*rnstepmax
            ! add terms to to evap/sub of cloud water
            dumx(i,k) = dumx(i,k) - faltndx*deltat*rnstepmax

            if (k>1) then
               if (present_qxsevap .or. present_qvlat .or. present_tlat) then
                  faltndqxe = (faloutx(i,k)-faloutx(i,k-1))*pdel_inv(i,k)
                  ! for output
                  if (present_qxsevap) qxsevap(i,k) = qxsevap(i,k) - (faltndqxe-faltndx)*rnstepmax
                  if (present_qvlat) qvlat(i,k) = qvlat(i,k) - (faltndqxe-faltndx)*rnstepmax
                  if (present_tlat) tlat(i,k) = tlat(i,k) + (faltndqxe-faltndx)*xxlx*rnstepmax
               end if
            end if

            if (present_xflx) xflx(i,k+1) = xflx(i,k+1) + faloutx(i,k) / g * rnstepmax
         end do

         ! units below are m/s
         ! sedimentation flux at surface is added to precip flux at surface
         ! to get total precip (cloud + precip water) rate
         if (present_prect) prect(i) = prect(i) + faloutx(i,nlev) / g * rnstepmax / 1000._r8
         if (present_preci) preci(i) = preci(i) + faloutx(i,nlev) / g * rnstepmax / 1000._r8
      end do  ! n loop of 1, nstep
   end do  ! i loop of 1, mgncol
   !$acc end parallel

   !$acc exit data delete(faloutx,dum1) async(queue)
end subroutine Sedimentation

!========================================================================
!2021-10-19: Add a new interface for the implicit sedimentation calculation;
!            Separate number/mass sediment for each class;
!========================================================================
subroutine Sedimentation_implicit(mgncol,nlev,deltat,zint,pdel,dumx,fx,check_qsmall, &
                                  xxtend,queue,xflx,qxsedten,prect,preci)

   integer,  intent(in)              :: mgncol,nlev
   real(r8), intent(in)              :: deltat
   real(r8), intent(in)              :: zint(mgncol,nlev+1)
   real(r8), intent(in)              :: pdel(mgncol,nlev)
   real(r8), intent(in)              :: dumx(mgncol,nlev)
   real(r8), intent(in)              :: fx(mgncol,nlev)
   logical,  intent(in)              :: check_qsmall
   real(r8), intent(inout)           :: xxtend(mgncol,nlev)
   integer,  intent(in)              :: queue
   real(r8), intent(inout), optional :: xflx(mgncol,nlev+1)
   real(r8), intent(inout), optional :: qxsedten(mgncol,nlev)
   real(r8), intent(inout), optional :: prect(mgncol)
   real(r8), intent(inout), optional :: preci(mgncol)

   ! Local variables
   integer  :: i,k
   real(r8) :: dum_2D(mgncol,nlev),flx(mgncol,nlev),precip(mgncol)
   logical  :: present_preci, present_xflx, present_qxsedten, present_prect

   present_preci = present(preci)
   present_xflx = present(xflx)
   present_qxsedten = present(qxsedten)
   present_prect = present(prect)

   !$acc enter data create (flx,dum_2D,precip) async(queue)

   !$acc parallel vector_length(VLENS) default(present) async(queue)
   !$acc loop gang vector collapse(2)
   do k=1,nlev
      do i=1,mgncol
         dum_2D(i,k) = dumx(i,k)
      enddo
   enddo
   !$acc end parallel

   call implicit_fall ( deltat, mgncol, 1, nlev, zint, fx, pdel, dum_2D, precip, flx, queue)

   !$acc parallel vector_length(VLENS) default(present) async(queue)
   !$acc loop gang vector collapse(2)
   do k=1,nlev
      do i=1,mgncol
         if ( check_qsmall ) then
            !h1g, 2019-11-26, ensure numerical stability
            if ( flx(i,k) .ge. qsmall .and. present_xflx ) xflx(i,k+1) = xflx(i,k+1) + flx(i,k) / g / deltat
         else
            if ( present_xflx ) xflx(i,k+1) = xflx(i,k+1) + flx(i,k) / g / deltat
         end if
         if ( present_qxsedten) qxsedten(i,k) = qxsedten(i,k) + (dum_2D(i,k) - dumx(i,k)) / deltat
         xxtend(i,k) = xxtend(i,k) + (dum_2D(i,k) - dumx(i,k)) / deltat
      enddo
   enddo
   !$acc end parallel

   !$acc parallel vector_length(VLENS) default(present) async(queue)
   !$acc loop gang vector
   do i=1,mgncol
      if ( precip(i) .ge. 0.0_r8 ) then        !h1g, 2019-11-26, ensure numerical stability
         if ( present_prect ) prect(i) = prect(i) + precip(i) / g / deltat / 1000._r8
         if ( present_preci ) preci(i) = preci(i) + precip(i) / g / deltat / 1000._r8
      endif
   enddo
   !$acc end parallel

   !$acc exit data delete(flx,dum_2D,precip) async(queue)

end subroutine Sedimentation_implicit

!========================================================================
!UTILITIES
!========================================================================


pure subroutine micro_pumas_get_cols(ncol, nlev, top_lev, mgncol, mgcols, &
     qcn, qin, qrn, qsn, qgr)

  ! Determines which columns microphysics should operate over by
  ! checking for non-zero cloud water/ice.

  integer, intent(in) :: ncol      ! Number of columns with meaningful data
  integer, intent(in) :: nlev      ! Number of levels to use
  integer, intent(in) :: top_lev   ! Top level for microphysics
  integer, intent(out) :: mgncol   ! Number of columns MG will use
  integer, allocatable, intent(out) :: mgcols(:) ! column indices

  real(r8), intent(in) :: qcn(:,:) ! cloud water mixing ratio (kg/kg)
  real(r8), intent(in) :: qin(:,:) ! cloud ice mixing ratio (kg/kg)
  real(r8), intent(in) :: qrn(:,:) ! rain mixing ratio (kg/kg)
  real(r8), intent(in) :: qsn(:,:) ! snow mixing ratio (kg/kg)
  real(r8), optional, intent(in) :: qgr(:,:) ! graupel mixing ratio (kg/kg)

  integer :: lev_offset  ! top_lev - 1 (defined here for consistency)
  logical :: ltrue(ncol) ! store tests for each column

  integer :: i, ii ! column indices

  if (allocated(mgcols)) deallocate(mgcols)

  lev_offset = top_lev - 1

  ! Using "any" along dimension 2 collapses across levels, but
  ! not columns, so we know if water is present at any level
  ! in each column.

  ltrue = any(qcn(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)
  ltrue = ltrue .or. any(qin(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)
  ltrue = ltrue .or. any(qrn(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)
  ltrue = ltrue .or. any(qsn(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)

  if(present(qgr)) ltrue = ltrue .or. any(qgr(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)

  ! Scan for true values to get a usable list of indices.

  mgncol = count(ltrue)
  allocate(mgcols(mgncol))
  i = 0
  do ii = 1,ncol
     if (ltrue(ii)) then
        i = i + 1
        mgcols(i) = ii
     end if
  end do

end subroutine micro_pumas_get_cols

! =======================================================================
! time - implicit monotonic scheme
! developed by sj lin, 2016
! =======================================================================

subroutine implicit_fall (dt, mgncol, ktop, kbot, ze, vt, dp, q, precip, m1, queue)

    implicit none

    integer, intent (in) :: mgncol                                   ! Number of columns in MG
    integer, intent (in) :: ktop,kbot                                ! Level range (top to bottom)
    real(r8), intent (in) :: dt                                      ! Time step
    real(r8), intent (in), dimension (mgncol,ktop:kbot+1) :: ze      ! Interface height (m)
    real(r8), intent (in), dimension (mgncol,ktop:kbot) :: vt, dp    ! fall speed and pressure difference across level
    real(r8), intent (inout), dimension (mgncol,ktop:kbot) :: q      ! mass
    real(r8), intent (out), dimension (mgncol,ktop:kbot) :: m1       ! Surface Flux
    real(r8), intent (out), dimension (mgncol) :: precip             ! Surface Precipitation
    integer, intent (in) :: queue                                    ! Stream ID for GPU asynchronous run

    ! Local variables
    real(r8), dimension (mgncol,ktop:kbot) :: dz, qm, dd
    integer :: i,k

    !$acc enter data create (dz,qm,dd) async(queue)

    !$acc parallel vector_length(VLENS) default(present) async(queue)
    !$acc loop gang vector collapse(2)
    do i = 1, mgncol
       do k = ktop, kbot
          dz (i,k) = ze (i,k) - ze (i,k + 1)
          dd (i,k) = dt * vt (i,k)
          q (i,k) = q (i,k) * dp (i,k)
       end do
    end do
    !$acc end parallel

    ! -----------------------------------------------------------------------
    ! sedimentation: non - vectorizable loop
    ! -----------------------------------------------------------------------

    !$acc parallel vector_length(VLENS) default(present) async(queue)
    !$acc loop gang vector
    do i = 1, mgncol
       qm (i,ktop) = q (i,ktop) / (dz (i,ktop) + dd (i,ktop))

       !$acc loop seq
       do k = ktop + 1, kbot
          qm (i,k) = (q (i,k) + dd (i,k - 1) * qm (i,k - 1)) / (dz (i,k) + dd (i,k))
       end do
    end do
    !$acc end parallel

    ! -----------------------------------------------------------------------
    ! qm is density at this stage
    ! -----------------------------------------------------------------------

    !$acc parallel vector_length(VLENS) default(present) async(queue)
    !$acc loop gang vector collapse(2)
    do i = 1, mgncol
       do k = ktop, kbot
          qm (i,k) = qm (i,k) * dz (i,k)
       end do
    end do
    !$acc end parallel

    ! -----------------------------------------------------------------------
    ! output mass fluxes: non - vectorizable loop
    ! -----------------------------------------------------------------------

    !$acc parallel vector_length(VLENS) default(present) async(queue)
    !$acc loop gang vector
    do i = 1, mgncol
       m1 (i,ktop) = q (i,ktop) - qm (i,ktop)

       !$acc loop seq
       do k = ktop + 1, kbot
          m1 (i,k) = m1 (i,k - 1) + q (i,k) - qm (i,k)
       end do
    end do
    !$acc end parallel

    !$acc parallel vector_length(VLENS) default(present) async(queue)
    !$acc loop gang vector
    do i = 1, mgncol
       precip(i) = m1 (i,kbot)
    end do
    !$acc end parallel

    ! -----------------------------------------------------------------------
    ! update:
    ! -----------------------------------------------------------------------

    !$acc parallel vector_length(VLENS) default(present) async(queue)
    !$acc loop gang vector collapse(2)
    do i = 1, mgncol
       do k = ktop, kbot
          q (i,k) = qm (i,k) / dp (i,k)
       end do
    end do
    !$acc end parallel

    !$acc exit data delete (dz,qm,dd) async(queue)

end subroutine implicit_fall


end module micro_pumas_v1
