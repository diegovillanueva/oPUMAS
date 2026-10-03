! 🧪📦 Standalone column driver for micro_pumas_tend (oPUMAS Julia twin, milestone M1).
! 🔁🧊 Reads one inputs NetCDF and one settings namelist, calls micro_pumas_init and
! micro_pumas_tend exactly as CAM does (interface/cesm215/micro_mg_cam.F90), and writes
! every output and every allocated proc_rates component to one NetCDF file.
! ⏱️📈 An optional fourth argument nsteps repeats the call, applying the tendencies to the
! state with forward Euler between calls (host forcing held fixed); the outputs are those
! of the last call. The Julia twin is time_loop! in julia/src/time_loop.jl.
! 🛠️📐 Build and run: tests/column/build.sh, tests/column/run.sh (see tests/column/README.md).
program pumas_column_driver
  use pumas_kinds, only: r8 => kind_r8
  use netcdf
  use wv_sat_methods, only: wv_sat_methods_init
  use micro_pumas_v1, only: micro_pumas_init, micro_pumas_tend
  use micro_pumas_diags, only: proc_rates_type
  use pumas_stochastic_collect_tau, only: ncd, pumas_stochastic_kernel_init
  use pumas_column_io, only: write_outputs
#ifdef CAPTURE
  use pumas_column_capture, only: capture_close
#endif
  use micro_pumas_sppe, only: micro_pumas_sppe_set
  implicit none

  ! 🌍📏 CESM 2.1.5 constants (cime/src/share/util/shr_const_mod.F90, cam physconst.F90).
  real(r8), parameter :: shr_boltz = 1.38065e-23_r8, shr_avogad = 6.02214e26_r8
  real(r8), parameter :: shr_rgas = shr_avogad*shr_boltz
  real(r8), parameter :: shr_mwdair = 28.966_r8, shr_mwwv = 18.016_r8
  real(r8), parameter :: gravit = 9.80616_r8
  real(r8), parameter :: rair = shr_rgas/shr_mwdair, rh2o = shr_rgas/shr_mwwv
  real(r8), parameter :: cpair = 1.00464e3_r8, tmelt = 273.15_r8, h2otrip = 273.16_r8
  real(r8), parameter :: latvap = 2.501e6_r8, latice = 3.337e5_r8
  real(r8), parameter :: tboil = 373.16_r8, ttrice = 20.00_r8
  real(r8), parameter :: epsilo = shr_mwwv/shr_mwdair

  ! 🎛️📄 Settings namelist; defaults are the values of case cesm215_pumas.F2000climo.bfb_dev.
  real(r8) :: rhmini = 0.80_r8, micro_mg_dcs = 500.e-6_r8
  logical :: micro_mg_do_hail = .false., micro_mg_do_graupel = .true.
  logical :: microp_uniform = .false., do_cldice = .true., use_hetfrz_classnuc = .true.
  character(len=16) :: micro_mg_precip_frac_method = 'in_cloud'
  real(r8) :: micro_mg_berg_eff_factor = 1._r8, micro_mg_accre_enhan_fact = 0.731_r8
  real(r8) :: micro_mg_autocon_fact = 0.01_r8, micro_mg_autocon_nd_exp = -1.1_r8
  real(r8) :: micro_mg_autocon_lwp_exp = 2.37_r8, micro_mg_homog_size = 25.e-6_r8
  real(r8) :: micro_mg_vtrmi_factor = 1.39_r8, micro_mg_vtrms_factor = 1._r8
  real(r8) :: micro_mg_effi_factor = 1._r8, micro_mg_iaccr_factor = 1._r8
  real(r8) :: micro_mg_max_nicons = 1.e8_r8
  logical :: allow_sed_supersat = .false.
  character(len=16) :: micro_pumas_warm_rain = 'kk2000'
  logical :: micro_mg_evap_sed_off = .false., micro_mg_icenuc_rh_off = .false.
  logical :: micro_mg_icenuc_use_meyers = .false., micro_mg_evap_scl_ifs = .false.
  logical :: micro_mg_evap_rhthrsh_ifs = .false., micro_mg_rainfreeze_ifs = .false.
  logical :: micro_mg_ifs_sed = .false., micro_mg_precip_fall_corr = .true.
  logical :: micro_mg_accre_sees_auto = .true., micro_mg_implicit_fall = .true.
  logical :: micro_mg_nccons = .false., micro_mg_nicons = .false.
  real(r8) :: micro_mg_ncnst = 100.e6_r8, micro_mg_ninst = 0.1e6_r8
  logical :: micro_mg_ngcons = .false.
  real(r8) :: micro_mg_ngnst = 0.0005e6_r8
  character(len=256) :: micro_pumas_kernel_file = ' '
  namelist /pumas_settings/ rhmini, micro_mg_dcs, micro_mg_do_hail, micro_mg_do_graupel, &
       microp_uniform, do_cldice, use_hetfrz_classnuc, micro_mg_precip_frac_method, &
       micro_mg_berg_eff_factor, micro_mg_accre_enhan_fact, micro_mg_autocon_fact, &
       micro_mg_autocon_nd_exp, micro_mg_autocon_lwp_exp, micro_mg_homog_size, &
       micro_mg_vtrmi_factor, micro_mg_vtrms_factor, micro_mg_effi_factor, &
       micro_mg_iaccr_factor, micro_mg_max_nicons, allow_sed_supersat, &
       micro_pumas_warm_rain, micro_mg_evap_sed_off, micro_mg_icenuc_rh_off, &
       micro_mg_icenuc_use_meyers, micro_mg_evap_scl_ifs, micro_mg_evap_rhthrsh_ifs, &
       micro_mg_rainfreeze_ifs, micro_mg_ifs_sed, micro_mg_precip_fall_corr, &
       micro_mg_accre_sees_auto, micro_mg_implicit_fall, micro_mg_nccons, micro_mg_nicons, &
       micro_mg_ncnst, micro_mg_ninst, micro_mg_ngcons, micro_mg_ngnst, &
       micro_pumas_kernel_file
  ! 🧊🎛️ Sppe process controls (micro_pumas_sppe.F90); defaults are upstream PUMAS.
  logical :: naai_het_also_in_mpc = .true., enable_nimax = .true., detrainment_ramp_liq = .false.
  logical :: onoff_contact_freezing = .true., onoff_sip = .true., rain_freeze = .true.
  real(r8) :: auto_liq_factor = 1._r8, riming_factor = 1._r8, factor_HM = 1._r8
  real(r8) :: micro_mg_ice_berg_ice_factor = 1._r8, micro_mg_ice_berg_snow_factor = 1._r8
  namelist /pumas_settings/ naai_het_also_in_mpc, enable_nimax, detrainment_ramp_liq, &
       onoff_contact_freezing, onoff_sip, rain_freeze, auto_liq_factor, riming_factor, &
       micro_mg_ice_berg_ice_factor, micro_mg_ice_berg_snow_factor, factor_HM

  character(len=1024) :: settings_file, inputs_file, output_file
  character(len=128) :: errstring
  character(len=32) :: nsteps_arg
  integer :: ncid, ncol, nlev, ndust, unit, ios, istep, nsteps
  real(r8) :: deltat
  type(proc_rates_type) :: pr

  real(r8), allocatable :: t(:,:)
  real(r8), allocatable :: q(:,:)
  real(r8), allocatable :: qcn(:,:)
  real(r8), allocatable :: qin(:,:)
  real(r8), allocatable :: ncn(:,:)
  real(r8), allocatable :: nin(:,:)
  real(r8), allocatable :: qrn(:,:)
  real(r8), allocatable :: qsn(:,:)
  real(r8), allocatable :: nrn(:,:)
  real(r8), allocatable :: nsn(:,:)
  real(r8), allocatable :: qgr(:,:)
  real(r8), allocatable :: ngr(:,:)
  real(r8), allocatable :: relvar(:,:)
  real(r8), allocatable :: accre_enhan(:,:)
  real(r8), allocatable :: p(:,:)
  real(r8), allocatable :: pdel(:,:)
  real(r8), allocatable :: cldn(:,:)
  real(r8), allocatable :: liqcldf(:,:)
  real(r8), allocatable :: icecldf(:,:)
  real(r8), allocatable :: qsatfac(:,:)
  real(r8), allocatable :: naai(:,:)
  real(r8), allocatable :: npccn(:,:)
  real(r8), allocatable :: tnd_qsnow(:,:)
  real(r8), allocatable :: tnd_nsnow(:,:)
  real(r8), allocatable :: re_ice(:,:)
  real(r8), allocatable :: frzimm(:,:)
  real(r8), allocatable :: frzcnt(:,:)
  real(r8), allocatable :: frzdep(:,:)
  real(r8), allocatable :: pint(:,:)
  real(r8), allocatable :: rndst(:,:,:), nacon(:,:,:)
  real(r8), allocatable :: tlat(:,:)
  real(r8), allocatable :: qvlat(:,:)
  real(r8), allocatable :: qctend(:,:)
  real(r8), allocatable :: qitend(:,:)
  real(r8), allocatable :: nctend(:,:)
  real(r8), allocatable :: nitend(:,:)
  real(r8), allocatable :: qrtend(:,:)
  real(r8), allocatable :: qstend(:,:)
  real(r8), allocatable :: nrtend(:,:)
  real(r8), allocatable :: nstend(:,:)
  real(r8), allocatable :: qgtend(:,:)
  real(r8), allocatable :: ngtend(:,:)
  real(r8), allocatable :: effc(:,:)
  real(r8), allocatable :: effc_fn(:,:)
  real(r8), allocatable :: effi(:,:)
  real(r8), allocatable :: sadice(:,:)
  real(r8), allocatable :: sadsnow(:,:)
  real(r8), allocatable :: nevapr(:,:)
  real(r8), allocatable :: am_evp_st(:,:)
  real(r8), allocatable :: prain(:,:)
  real(r8), allocatable :: cmeout(:,:)
  real(r8), allocatable :: deffi(:,:)
  real(r8), allocatable :: pgamrad(:,:)
  real(r8), allocatable :: lamcrad(:,:)
  real(r8), allocatable :: qsout(:,:)
  real(r8), allocatable :: dsout(:,:)
  real(r8), allocatable :: qgout(:,:)
  real(r8), allocatable :: ngout(:,:)
  real(r8), allocatable :: dgout(:,:)
  real(r8), allocatable :: qrout(:,:)
  real(r8), allocatable :: reff_rain(:,:)
  real(r8), allocatable :: reff_snow(:,:)
  real(r8), allocatable :: reff_grau(:,:)
  real(r8), allocatable :: nrout(:,:)
  real(r8), allocatable :: nsout(:,:)
  real(r8), allocatable :: refl(:,:)
  real(r8), allocatable :: arefl(:,:)
  real(r8), allocatable :: areflz(:,:)
  real(r8), allocatable :: frefl(:,:)
  real(r8), allocatable :: csrfl(:,:)
  real(r8), allocatable :: acsrfl(:,:)
  real(r8), allocatable :: fcsrfl(:,:)
  real(r8), allocatable :: refl10cm(:,:)
  real(r8), allocatable :: reflz10cm(:,:)
  real(r8), allocatable :: rercld(:,:)
  real(r8), allocatable :: ncai(:,:)
  real(r8), allocatable :: ncal(:,:)
  real(r8), allocatable :: qrout2(:,:)
  real(r8), allocatable :: qsout2(:,:)
  real(r8), allocatable :: nrout2(:,:)
  real(r8), allocatable :: nsout2(:,:)
  real(r8), allocatable :: drout2(:,:)
  real(r8), allocatable :: dsout2(:,:)
  real(r8), allocatable :: qgout2(:,:)
  real(r8), allocatable :: ngout2(:,:)
  real(r8), allocatable :: dgout2(:,:)
  real(r8), allocatable :: freqg(:,:)
  real(r8), allocatable :: freqs(:,:)
  real(r8), allocatable :: freqr(:,:)
  real(r8), allocatable :: nfice(:,:)
  real(r8), allocatable :: qcrat(:,:)
  real(r8), allocatable :: prer_evap(:,:)
  real(r8), allocatable :: qcsinksum_rate1ord(:,:)
  real(r8), allocatable :: lflx(:,:)
  real(r8), allocatable :: iflx(:,:)
  real(r8), allocatable :: gflx(:,:)
  real(r8), allocatable :: rflx(:,:)
  real(r8), allocatable :: sflx(:,:)
  real(r8), allocatable :: prect(:), preci(:)

  call get_command_argument(1, settings_file)
  call get_command_argument(2, inputs_file)
  call get_command_argument(3, output_file)
  if (len_trim(output_file) == 0) then
     write(*,*) 'usage: pumas_column_driver settings.nml inputs.nc output.nc [nsteps]'
     stop 2
  end if
  nsteps = 1
  if (command_argument_count() >= 4) then
     call get_command_argument(4, nsteps_arg)
     read(nsteps_arg, *) nsteps
  end if

  open(newunit=unit, file=trim(settings_file), status='old', action='read', iostat=ios)
  if (ios /= 0) call die('cannot open '//trim(settings_file))
  read(unit, nml=pumas_settings, iostat=ios)
  if (ios /= 0) call die('cannot read namelist pumas_settings from '//trim(settings_file))
  close(unit)

  call wv_sat_methods_init(r8, tmelt, h2otrip, tboil, ttrice, epsilo, errstring)
  call check_err(errstring, 'wv_sat_methods_init')

  ! 🔌📞 Same argument list as micro_mg_cam.F90 (nrcons, nscons off; no emulator files).
  call micro_pumas_init( &
       r8, gravit, rair, rh2o, cpair, &
       tmelt, latvap, latice, rhmini, &
       micro_mg_dcs, micro_mg_do_hail, micro_mg_do_graupel, &
       microp_uniform, do_cldice, use_hetfrz_classnuc, &
       micro_mg_precip_frac_method, micro_mg_berg_eff_factor, &
       micro_mg_accre_enhan_fact, micro_mg_autocon_fact, &
       micro_mg_autocon_nd_exp, micro_mg_autocon_lwp_exp, micro_mg_homog_size, &
       micro_mg_vtrmi_factor, micro_mg_vtrms_factor, micro_mg_effi_factor, &
       micro_mg_iaccr_factor, micro_mg_max_nicons, &
       allow_sed_supersat, micro_pumas_warm_rain, &
       micro_mg_evap_sed_off, micro_mg_icenuc_rh_off, micro_mg_icenuc_use_meyers, &
       micro_mg_evap_scl_ifs, micro_mg_evap_rhthrsh_ifs, &
       micro_mg_rainfreeze_ifs, micro_mg_ifs_sed, micro_mg_precip_fall_corr, &
       micro_mg_accre_sees_auto, micro_mg_implicit_fall, &
       micro_mg_nccons, micro_mg_nicons, micro_mg_ncnst, micro_mg_ninst, &
       micro_mg_ngcons, micro_mg_ngnst, .false., 0._r8, .false., 0._r8, &
       '', '', '', 6, errstring)
  call check_err(errstring, 'micro_pumas_init')
  ! 🧊🎛️ As micro_mg_cam_init: the controls are set once, after micro_pumas_init.
  call micro_pumas_sppe_set(naai_het_also_in_mpc, enable_nimax, detrainment_ramp_liq, &
       onoff_contact_freezing, onoff_sip, rain_freeze, auto_liq_factor, riming_factor, &
       micro_mg_ice_berg_ice_factor, micro_mg_ice_berg_snow_factor, factor_HM)
  if (trim(micro_pumas_warm_rain) == 'tau') &
       call pumas_stochastic_kernel_init(trim(micro_pumas_kernel_file))

  ! 📥🌧️ Inputs.
  call nc(nf90_open(trim(inputs_file), nf90_nowrite, ncid), 'open '//trim(inputs_file))
  ncol = dimlen(ncid, 'col'); nlev = dimlen(ncid, 'lev'); ndust = dimlen(ncid, 'dust_bin')
  call nc(nf90_get_att(ncid, nf90_global, 'deltat', deltat), 'attribute deltat')
  call get2(ncid, 't', t, ncol, nlev)
  call get2(ncid, 'q', q, ncol, nlev)
  call get2(ncid, 'qcn', qcn, ncol, nlev)
  call get2(ncid, 'qin', qin, ncol, nlev)
  call get2(ncid, 'ncn', ncn, ncol, nlev)
  call get2(ncid, 'nin', nin, ncol, nlev)
  call get2(ncid, 'qrn', qrn, ncol, nlev)
  call get2(ncid, 'qsn', qsn, ncol, nlev)
  call get2(ncid, 'nrn', nrn, ncol, nlev)
  call get2(ncid, 'nsn', nsn, ncol, nlev)
  call get2(ncid, 'qgr', qgr, ncol, nlev)
  call get2(ncid, 'ngr', ngr, ncol, nlev)
  call get2(ncid, 'relvar', relvar, ncol, nlev)
  call get2(ncid, 'accre_enhan', accre_enhan, ncol, nlev)
  call get2(ncid, 'p', p, ncol, nlev)
  call get2(ncid, 'pdel', pdel, ncol, nlev)
  call get2(ncid, 'cldn', cldn, ncol, nlev)
  call get2(ncid, 'liqcldf', liqcldf, ncol, nlev)
  call get2(ncid, 'icecldf', icecldf, ncol, nlev)
  call get2(ncid, 'qsatfac', qsatfac, ncol, nlev)
  call get2(ncid, 'naai', naai, ncol, nlev)
  call get2(ncid, 'npccn', npccn, ncol, nlev)
  call get2(ncid, 'tnd_qsnow', tnd_qsnow, ncol, nlev)
  call get2(ncid, 'tnd_nsnow', tnd_nsnow, ncol, nlev)
  call get2(ncid, 're_ice', re_ice, ncol, nlev)
  call get2(ncid, 'frzimm', frzimm, ncol, nlev)
  call get2(ncid, 'frzcnt', frzcnt, ncol, nlev)
  call get2(ncid, 'frzdep', frzdep, ncol, nlev)
  allocate(pint(ncol, nlev+1)); call get_any(ncid, 'pint', pint)
  allocate(rndst(ncol, nlev, ndust)); call get_any(ncid, 'rndst', rndst)
  allocate(nacon(ncol, nlev, ndust)); call get_any(ncid, 'nacon', nacon)
  call nc(nf90_close(ncid), 'close inputs')

  allocate(tlat(ncol, nlev)); tlat = 0._r8
  allocate(qvlat(ncol, nlev)); qvlat = 0._r8
  allocate(qctend(ncol, nlev)); qctend = 0._r8
  allocate(qitend(ncol, nlev)); qitend = 0._r8
  allocate(nctend(ncol, nlev)); nctend = 0._r8
  allocate(nitend(ncol, nlev)); nitend = 0._r8
  allocate(qrtend(ncol, nlev)); qrtend = 0._r8
  allocate(qstend(ncol, nlev)); qstend = 0._r8
  allocate(nrtend(ncol, nlev)); nrtend = 0._r8
  allocate(nstend(ncol, nlev)); nstend = 0._r8
  allocate(qgtend(ncol, nlev)); qgtend = 0._r8
  allocate(ngtend(ncol, nlev)); ngtend = 0._r8
  allocate(effc(ncol, nlev)); effc = 0._r8
  allocate(effc_fn(ncol, nlev)); effc_fn = 0._r8
  allocate(effi(ncol, nlev)); effi = 0._r8
  allocate(sadice(ncol, nlev)); sadice = 0._r8
  allocate(sadsnow(ncol, nlev)); sadsnow = 0._r8
  allocate(nevapr(ncol, nlev)); nevapr = 0._r8
  allocate(am_evp_st(ncol, nlev)); am_evp_st = 0._r8
  allocate(prain(ncol, nlev)); prain = 0._r8
  allocate(cmeout(ncol, nlev)); cmeout = 0._r8
  allocate(deffi(ncol, nlev)); deffi = 0._r8
  allocate(pgamrad(ncol, nlev)); pgamrad = 0._r8
  allocate(lamcrad(ncol, nlev)); lamcrad = 0._r8
  allocate(qsout(ncol, nlev)); qsout = 0._r8
  allocate(dsout(ncol, nlev)); dsout = 0._r8
  allocate(qgout(ncol, nlev)); qgout = 0._r8
  allocate(ngout(ncol, nlev)); ngout = 0._r8
  allocate(dgout(ncol, nlev)); dgout = 0._r8
  allocate(qrout(ncol, nlev)); qrout = 0._r8
  allocate(reff_rain(ncol, nlev)); reff_rain = 0._r8
  allocate(reff_snow(ncol, nlev)); reff_snow = 0._r8
  allocate(reff_grau(ncol, nlev)); reff_grau = 0._r8
  allocate(nrout(ncol, nlev)); nrout = 0._r8
  allocate(nsout(ncol, nlev)); nsout = 0._r8
  allocate(refl(ncol, nlev)); refl = 0._r8
  allocate(arefl(ncol, nlev)); arefl = 0._r8
  allocate(areflz(ncol, nlev)); areflz = 0._r8
  allocate(frefl(ncol, nlev)); frefl = 0._r8
  allocate(csrfl(ncol, nlev)); csrfl = 0._r8
  allocate(acsrfl(ncol, nlev)); acsrfl = 0._r8
  allocate(fcsrfl(ncol, nlev)); fcsrfl = 0._r8
  allocate(refl10cm(ncol, nlev)); refl10cm = 0._r8
  allocate(reflz10cm(ncol, nlev)); reflz10cm = 0._r8
  allocate(rercld(ncol, nlev)); rercld = 0._r8
  allocate(ncai(ncol, nlev)); ncai = 0._r8
  allocate(ncal(ncol, nlev)); ncal = 0._r8
  allocate(qrout2(ncol, nlev)); qrout2 = 0._r8
  allocate(qsout2(ncol, nlev)); qsout2 = 0._r8
  allocate(nrout2(ncol, nlev)); nrout2 = 0._r8
  allocate(nsout2(ncol, nlev)); nsout2 = 0._r8
  allocate(drout2(ncol, nlev)); drout2 = 0._r8
  allocate(dsout2(ncol, nlev)); dsout2 = 0._r8
  allocate(qgout2(ncol, nlev)); qgout2 = 0._r8
  allocate(ngout2(ncol, nlev)); ngout2 = 0._r8
  allocate(dgout2(ncol, nlev)); dgout2 = 0._r8
  allocate(freqg(ncol, nlev)); freqg = 0._r8
  allocate(freqs(ncol, nlev)); freqs = 0._r8
  allocate(freqr(ncol, nlev)); freqr = 0._r8
  allocate(nfice(ncol, nlev)); nfice = 0._r8
  allocate(qcrat(ncol, nlev)); qcrat = 0._r8
  allocate(prer_evap(ncol, nlev)); prer_evap = 0._r8
  allocate(qcsinksum_rate1ord(ncol, nlev)); qcsinksum_rate1ord = 0._r8
  allocate(lflx(ncol, nlev+1)); lflx = 0._r8
  allocate(iflx(ncol, nlev+1)); iflx = 0._r8
  allocate(gflx(ncol, nlev+1)); gflx = 0._r8
  allocate(rflx(ncol, nlev+1)); rflx = 0._r8
  allocate(sflx(ncol, nlev+1)); sflx = 0._r8
  allocate(prect(ncol), preci(ncol)); prect = 0._r8; preci = 0._r8

  call pr%allocate(ncol, nlev, ncd, micro_pumas_warm_rain, errstring)
  call check_err(errstring, 'proc_rates%allocate')

  do istep = 1, nsteps
  if (istep > 1) call advance_state()
  call micro_pumas_tend( &
       ncol, nlev, deltat, &
       t, q, qcn, qin, ncn, nin, qrn, qsn, nrn, nsn, qgr, ngr, &
       relvar, accre_enhan, p, pdel, pint, &
       cldn, liqcldf, icecldf, qsatfac, &
       qcsinksum_rate1ord, naai, npccn, rndst, nacon, &
       tlat, qvlat, qctend, qitend, nctend, nitend, &
       qrtend, qstend, nrtend, nstend, qgtend, ngtend, &
       effc, effc_fn, effi, sadice, sadsnow, prect, preci, &
       nevapr, am_evp_st, prain, cmeout, deffi, pgamrad, lamcrad, &
       qsout, dsout, qgout, ngout, dgout, &
       lflx, iflx, gflx, rflx, sflx, qrout, &
       reff_rain, reff_snow, reff_grau, nrout, nsout, &
       refl, arefl, areflz, frefl, csrfl, acsrfl, &
       fcsrfl, refl10cm, reflz10cm, rercld, ncai, ncal, &
       qrout2, qsout2, nrout2, nsout2, drout2, dsout2, &
       qgout2, ngout2, dgout2, freqg, freqs, freqr, &
       nfice, qcrat, pr, errstring, &
       tnd_qsnow, tnd_nsnow, re_ice, prer_evap, &
       frzimm, frzcnt, frzdep)
  call check_err(errstring, 'micro_pumas_tend')
  end do
#ifdef CAPTURE
  call capture_close()
#endif

  ! 📤🧾 Outputs, written by the same module as the CAM dump (pumas_column_io.F90).
  call write_outputs(output_file, deltat, qcsinksum_rate1ord, tlat, qvlat, qctend, qitend, &
       nctend, nitend, qrtend, qstend, nrtend, nstend, qgtend, ngtend, effc, effc_fn, effi, &
       sadice, sadsnow, nevapr, am_evp_st, prain, cmeout, deffi, pgamrad, lamcrad, qsout, &
       dsout, qgout, ngout, dgout, qrout, reff_rain, reff_snow, reff_grau, nrout, nsout, &
       refl, arefl, areflz, frefl, csrfl, acsrfl, fcsrfl, refl10cm, reflz10cm, rercld, ncai, &
       ncal, qrout2, qsout2, nrout2, nsout2, drout2, dsout2, qgout2, ngout2, dgout2, freqg, &
       freqs, freqr, nfice, qcrat, prer_evap, lflx, iflx, gflx, rflx, sflx, prect, preci, pr)
  write(*,'(a)') trim(output_file)

contains

  ! ⏱️📈 Forward Euler with the tendencies of the last call; negative q and N become zero.
  subroutine advance_state()
    t = t + tlat/cpair*deltat
    call euler(q, qvlat)
    call euler(qcn, qctend); call euler(ncn, nctend)
    call euler(qin, qitend); call euler(nin, nitend)
    call euler(qrn, qrtend); call euler(nrn, nrtend)
    call euler(qsn, qstend); call euler(nsn, nstend)
    call euler(qgr, qgtend); call euler(ngr, ngtend)
  end subroutine advance_state

  subroutine euler(x, tendency)
    real(r8), intent(inout) :: x(:,:)
    real(r8), intent(in) :: tendency(:,:)
    x = x + tendency*deltat
    where (x < 0._r8) x = 0._r8
  end subroutine euler

  subroutine die(msg)
    character(len=*), intent(in) :: msg
    write(0,'(a)') 'pumas_column_driver: '//trim(msg)
    stop 1
  end subroutine die

  subroutine check_err(err, where)
    character(len=*), intent(in) :: err, where
    if (len_trim(err) > 0) call die(trim(where)//': '//trim(err))
  end subroutine check_err

  subroutine nc(status, what)
    integer, intent(in) :: status
    character(len=*), intent(in) :: what
    if (status /= nf90_noerr) call die(trim(what)//': '//trim(nf90_strerror(status)))
  end subroutine nc

  integer function dimlen(ncid, name)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer :: dimid
    call nc(nf90_inq_dimid(ncid, name, dimid), 'dimension '//name)
    call nc(nf90_inquire_dimension(ncid, dimid, len=dimlen), 'dimension '//name)
  end function dimlen

  subroutine get2(ncid, name, arr, ncol, nlev)
    integer, intent(in) :: ncid, ncol, nlev
    character(len=*), intent(in) :: name
    real(r8), allocatable, intent(out) :: arr(:,:)
    allocate(arr(ncol, nlev))
    call get_any(ncid, name, arr)
  end subroutine get2

  subroutine get_any(ncid, name, arr)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(r8), intent(out) :: arr(..)
    integer :: varid
    call nc(nf90_inq_varid(ncid, name, varid), 'variable '//name)
    select rank (arr)
    rank (2)
       call nc(nf90_get_var(ncid, varid, arr), 'read '//name)
    rank (3)
       call nc(nf90_get_var(ncid, varid, arr), 'read '//name)
    end select
  end subroutine get_any

end program pumas_column_driver
