! 🧪📦 NetCDF I/O shared by the CAM dump patch (tests/column/dump/) and the standalone
! column driver (tests/column/driver.F90), so both write inputs and outputs of
! micro_pumas_tend in one format: dimensions (col, lev), (col, ilev), (col, lev, dust_bin)
! and (col); variables named like the arguments of micro_pumas_tend; proc_rates
! components prefixed with proc_rates_.
! ⚙️🔁 The argument lists follow micro_pumas_tend; the proc_rates list holds every 2D
! component of proc_rates_type in micro_pumas_diags.F90.
module pumas_column_io
  use pumas_kinds, only: r8 => kind_r8
  use netcdf
  use micro_pumas_diags, only: proc_rates_type
  implicit none
  private
  public :: write_inputs, write_outputs, read_dump_control

contains

  ! 📥🧾 All inputs of one micro_pumas_tend call; lat and lon (degrees) are optional labels.
  subroutine write_inputs(file, deltat, t, q, qcn, qin, ncn, nin, qrn, qsn, nrn, nsn, qgr, ngr, relvar, accre_enhan, &
       p, pdel, pint, cldn, liqcldf, icecldf, qsatfac, naai, npccn, rndst, nacon, &
       tnd_qsnow, tnd_nsnow, re_ice, frzimm, frzcnt, frzdep, &
       lat, lon, nstep, lchnk)
    character(len=*), intent(in) :: file
    real(r8), intent(in) :: deltat
    real(r8), intent(in) :: t(:,:)
    real(r8), intent(in) :: q(:,:)
    real(r8), intent(in) :: qcn(:,:)
    real(r8), intent(in) :: qin(:,:)
    real(r8), intent(in) :: ncn(:,:)
    real(r8), intent(in) :: nin(:,:)
    real(r8), intent(in) :: qrn(:,:)
    real(r8), intent(in) :: qsn(:,:)
    real(r8), intent(in) :: nrn(:,:)
    real(r8), intent(in) :: nsn(:,:)
    real(r8), intent(in) :: qgr(:,:)
    real(r8), intent(in) :: ngr(:,:)
    real(r8), intent(in) :: relvar(:,:)
    real(r8), intent(in) :: accre_enhan(:,:)
    real(r8), intent(in) :: p(:,:)
    real(r8), intent(in) :: pdel(:,:)
    real(r8), intent(in) :: pint(:,:)
    real(r8), intent(in) :: cldn(:,:)
    real(r8), intent(in) :: liqcldf(:,:)
    real(r8), intent(in) :: icecldf(:,:)
    real(r8), intent(in) :: qsatfac(:,:)
    real(r8), intent(in) :: naai(:,:)
    real(r8), intent(in) :: npccn(:,:)
    real(r8), intent(in) :: tnd_qsnow(:,:)
    real(r8), intent(in) :: tnd_nsnow(:,:)
    real(r8), intent(in) :: re_ice(:,:)
    real(r8), intent(in) :: frzimm(:,:)
    real(r8), intent(in) :: frzcnt(:,:)
    real(r8), intent(in) :: frzdep(:,:)
    real(r8), intent(in) :: rndst(:,:,:), nacon(:,:,:)
    real(r8), intent(in), optional :: lat(:), lon(:)
    integer, intent(in), optional :: nstep, lchnk
    integer :: ncid
    ncid = create(file, deltat, size(t, 2))
    if (present(nstep)) call nc(nf90_put_att(ncid, nf90_global, 'nstep', nstep), 'nstep')
    if (present(lchnk)) call nc(nf90_put_att(ncid, nf90_global, 'lchnk', lchnk), 'lchnk')
    call put2(ncid, 't', t)
    call put2(ncid, 'q', q)
    call put2(ncid, 'qcn', qcn)
    call put2(ncid, 'qin', qin)
    call put2(ncid, 'ncn', ncn)
    call put2(ncid, 'nin', nin)
    call put2(ncid, 'qrn', qrn)
    call put2(ncid, 'qsn', qsn)
    call put2(ncid, 'nrn', nrn)
    call put2(ncid, 'nsn', nsn)
    call put2(ncid, 'qgr', qgr)
    call put2(ncid, 'ngr', ngr)
    call put2(ncid, 'relvar', relvar)
    call put2(ncid, 'accre_enhan', accre_enhan)
    call put2(ncid, 'p', p)
    call put2(ncid, 'pdel', pdel)
    call put2(ncid, 'pint', pint)
    call put2(ncid, 'cldn', cldn)
    call put2(ncid, 'liqcldf', liqcldf)
    call put2(ncid, 'icecldf', icecldf)
    call put2(ncid, 'qsatfac', qsatfac)
    call put2(ncid, 'naai', naai)
    call put2(ncid, 'npccn', npccn)
    call put2(ncid, 'tnd_qsnow', tnd_qsnow)
    call put2(ncid, 'tnd_nsnow', tnd_nsnow)
    call put2(ncid, 're_ice', re_ice)
    call put2(ncid, 'frzimm', frzimm)
    call put2(ncid, 'frzcnt', frzcnt)
    call put2(ncid, 'frzdep', frzdep)
    call put3(ncid, 'rndst', rndst)
    call put3(ncid, 'nacon', nacon)
    if (present(lat)) call put1(ncid, 'lat', lat)
    if (present(lon)) call put1(ncid, 'lon', lon)
    call nc(nf90_close(ncid), 'close '//trim(file))
  end subroutine write_inputs

  ! 📤🧾 All outputs of one micro_pumas_tend call and every allocated proc_rates field.
  subroutine write_outputs(file, deltat, qcsinksum_rate1ord, tlat, qvlat, qctend, qitend, nctend, nitend, qrtend, &
       qstend, nrtend, nstend, qgtend, ngtend, effc, effc_fn, effi, sadice, sadsnow, &
       nevapr, am_evp_st, prain, cmeout, deffi, pgamrad, lamcrad, qsout, dsout, &
       qgout, ngout, dgout, qrout, reff_rain, reff_snow, reff_grau, nrout, nsout, &
       refl, arefl, areflz, frefl, csrfl, acsrfl, fcsrfl, refl10cm, reflz10cm, &
       rercld, ncai, ncal, qrout2, qsout2, nrout2, nsout2, drout2, dsout2, qgout2, &
       ngout2, dgout2, freqg, freqs, freqr, nfice, qcrat, prer_evap, lflx, iflx, &
       gflx, rflx, sflx, prect, preci, &
       pr)
    character(len=*), intent(in) :: file
    real(r8), intent(in) :: deltat
    real(r8), intent(in) :: qcsinksum_rate1ord(:,:)
    real(r8), intent(in) :: tlat(:,:)
    real(r8), intent(in) :: qvlat(:,:)
    real(r8), intent(in) :: qctend(:,:)
    real(r8), intent(in) :: qitend(:,:)
    real(r8), intent(in) :: nctend(:,:)
    real(r8), intent(in) :: nitend(:,:)
    real(r8), intent(in) :: qrtend(:,:)
    real(r8), intent(in) :: qstend(:,:)
    real(r8), intent(in) :: nrtend(:,:)
    real(r8), intent(in) :: nstend(:,:)
    real(r8), intent(in) :: qgtend(:,:)
    real(r8), intent(in) :: ngtend(:,:)
    real(r8), intent(in) :: effc(:,:)
    real(r8), intent(in) :: effc_fn(:,:)
    real(r8), intent(in) :: effi(:,:)
    real(r8), intent(in) :: sadice(:,:)
    real(r8), intent(in) :: sadsnow(:,:)
    real(r8), intent(in) :: nevapr(:,:)
    real(r8), intent(in) :: am_evp_st(:,:)
    real(r8), intent(in) :: prain(:,:)
    real(r8), intent(in) :: cmeout(:,:)
    real(r8), intent(in) :: deffi(:,:)
    real(r8), intent(in) :: pgamrad(:,:)
    real(r8), intent(in) :: lamcrad(:,:)
    real(r8), intent(in) :: qsout(:,:)
    real(r8), intent(in) :: dsout(:,:)
    real(r8), intent(in) :: qgout(:,:)
    real(r8), intent(in) :: ngout(:,:)
    real(r8), intent(in) :: dgout(:,:)
    real(r8), intent(in) :: qrout(:,:)
    real(r8), intent(in) :: reff_rain(:,:)
    real(r8), intent(in) :: reff_snow(:,:)
    real(r8), intent(in) :: reff_grau(:,:)
    real(r8), intent(in) :: nrout(:,:)
    real(r8), intent(in) :: nsout(:,:)
    real(r8), intent(in) :: refl(:,:)
    real(r8), intent(in) :: arefl(:,:)
    real(r8), intent(in) :: areflz(:,:)
    real(r8), intent(in) :: frefl(:,:)
    real(r8), intent(in) :: csrfl(:,:)
    real(r8), intent(in) :: acsrfl(:,:)
    real(r8), intent(in) :: fcsrfl(:,:)
    real(r8), intent(in) :: refl10cm(:,:)
    real(r8), intent(in) :: reflz10cm(:,:)
    real(r8), intent(in) :: rercld(:,:)
    real(r8), intent(in) :: ncai(:,:)
    real(r8), intent(in) :: ncal(:,:)
    real(r8), intent(in) :: qrout2(:,:)
    real(r8), intent(in) :: qsout2(:,:)
    real(r8), intent(in) :: nrout2(:,:)
    real(r8), intent(in) :: nsout2(:,:)
    real(r8), intent(in) :: drout2(:,:)
    real(r8), intent(in) :: dsout2(:,:)
    real(r8), intent(in) :: qgout2(:,:)
    real(r8), intent(in) :: ngout2(:,:)
    real(r8), intent(in) :: dgout2(:,:)
    real(r8), intent(in) :: freqg(:,:)
    real(r8), intent(in) :: freqs(:,:)
    real(r8), intent(in) :: freqr(:,:)
    real(r8), intent(in) :: nfice(:,:)
    real(r8), intent(in) :: qcrat(:,:)
    real(r8), intent(in) :: prer_evap(:,:)
    real(r8), intent(in) :: lflx(:,:)
    real(r8), intent(in) :: iflx(:,:)
    real(r8), intent(in) :: gflx(:,:)
    real(r8), intent(in) :: rflx(:,:)
    real(r8), intent(in) :: sflx(:,:)
    real(r8), intent(in) :: prect(:), preci(:)
    type(proc_rates_type), intent(in) :: pr
    integer :: ncid
    ncid = create(file, deltat, size(tlat, 2))
    call put2(ncid, 'qcsinksum_rate1ord', qcsinksum_rate1ord)
    call put2(ncid, 'tlat', tlat)
    call put2(ncid, 'qvlat', qvlat)
    call put2(ncid, 'qctend', qctend)
    call put2(ncid, 'qitend', qitend)
    call put2(ncid, 'nctend', nctend)
    call put2(ncid, 'nitend', nitend)
    call put2(ncid, 'qrtend', qrtend)
    call put2(ncid, 'qstend', qstend)
    call put2(ncid, 'nrtend', nrtend)
    call put2(ncid, 'nstend', nstend)
    call put2(ncid, 'qgtend', qgtend)
    call put2(ncid, 'ngtend', ngtend)
    call put2(ncid, 'effc', effc)
    call put2(ncid, 'effc_fn', effc_fn)
    call put2(ncid, 'effi', effi)
    call put2(ncid, 'sadice', sadice)
    call put2(ncid, 'sadsnow', sadsnow)
    call put2(ncid, 'nevapr', nevapr)
    call put2(ncid, 'am_evp_st', am_evp_st)
    call put2(ncid, 'prain', prain)
    call put2(ncid, 'cmeout', cmeout)
    call put2(ncid, 'deffi', deffi)
    call put2(ncid, 'pgamrad', pgamrad)
    call put2(ncid, 'lamcrad', lamcrad)
    call put2(ncid, 'qsout', qsout)
    call put2(ncid, 'dsout', dsout)
    call put2(ncid, 'qgout', qgout)
    call put2(ncid, 'ngout', ngout)
    call put2(ncid, 'dgout', dgout)
    call put2(ncid, 'qrout', qrout)
    call put2(ncid, 'reff_rain', reff_rain)
    call put2(ncid, 'reff_snow', reff_snow)
    call put2(ncid, 'reff_grau', reff_grau)
    call put2(ncid, 'nrout', nrout)
    call put2(ncid, 'nsout', nsout)
    call put2(ncid, 'refl', refl)
    call put2(ncid, 'arefl', arefl)
    call put2(ncid, 'areflz', areflz)
    call put2(ncid, 'frefl', frefl)
    call put2(ncid, 'csrfl', csrfl)
    call put2(ncid, 'acsrfl', acsrfl)
    call put2(ncid, 'fcsrfl', fcsrfl)
    call put2(ncid, 'refl10cm', refl10cm)
    call put2(ncid, 'reflz10cm', reflz10cm)
    call put2(ncid, 'rercld', rercld)
    call put2(ncid, 'ncai', ncai)
    call put2(ncid, 'ncal', ncal)
    call put2(ncid, 'qrout2', qrout2)
    call put2(ncid, 'qsout2', qsout2)
    call put2(ncid, 'nrout2', nrout2)
    call put2(ncid, 'nsout2', nsout2)
    call put2(ncid, 'drout2', drout2)
    call put2(ncid, 'dsout2', dsout2)
    call put2(ncid, 'qgout2', qgout2)
    call put2(ncid, 'ngout2', ngout2)
    call put2(ncid, 'dgout2', dgout2)
    call put2(ncid, 'freqg', freqg)
    call put2(ncid, 'freqs', freqs)
    call put2(ncid, 'freqr', freqr)
    call put2(ncid, 'nfice', nfice)
    call put2(ncid, 'qcrat', qcrat)
    call put2(ncid, 'prer_evap', prer_evap)
    call put2(ncid, 'lflx', lflx)
    call put2(ncid, 'iflx', iflx)
    call put2(ncid, 'gflx', gflx)
    call put2(ncid, 'rflx', rflx)
    call put2(ncid, 'sflx', sflx)
    call put1(ncid, 'prect', prect)
    call put1(ncid, 'preci', preci)
    if (allocated(pr%prodsnow)) call put2(ncid, 'proc_rates_prodsnow', pr%prodsnow)
    if (allocated(pr%evapsnow)) call put2(ncid, 'proc_rates_evapsnow', pr%evapsnow)
    if (allocated(pr%qcsevap)) call put2(ncid, 'proc_rates_qcsevap', pr%qcsevap)
    if (allocated(pr%qisevap)) call put2(ncid, 'proc_rates_qisevap', pr%qisevap)
    if (allocated(pr%qvres)) call put2(ncid, 'proc_rates_qvres', pr%qvres)
    if (allocated(pr%cmeitot)) call put2(ncid, 'proc_rates_cmeitot', pr%cmeitot)
    if (allocated(pr%vtrmc)) call put2(ncid, 'proc_rates_vtrmc', pr%vtrmc)
    if (allocated(pr%vtrmi)) call put2(ncid, 'proc_rates_vtrmi', pr%vtrmi)
    if (allocated(pr%umr)) call put2(ncid, 'proc_rates_umr', pr%umr)
    if (allocated(pr%ums)) call put2(ncid, 'proc_rates_ums', pr%ums)
    if (allocated(pr%umg)) call put2(ncid, 'proc_rates_umg', pr%umg)
    if (allocated(pr%qgsedten)) call put2(ncid, 'proc_rates_qgsedten', pr%qgsedten)
    if (allocated(pr%qcsedten)) call put2(ncid, 'proc_rates_qcsedten', pr%qcsedten)
    if (allocated(pr%qisedten)) call put2(ncid, 'proc_rates_qisedten', pr%qisedten)
    if (allocated(pr%qrsedten)) call put2(ncid, 'proc_rates_qrsedten', pr%qrsedten)
    if (allocated(pr%qssedten)) call put2(ncid, 'proc_rates_qssedten', pr%qssedten)
    if (allocated(pr%pratot)) call put2(ncid, 'proc_rates_pratot', pr%pratot)
    if (allocated(pr%prctot)) call put2(ncid, 'proc_rates_prctot', pr%prctot)
    if (allocated(pr%mnuccctot)) call put2(ncid, 'proc_rates_mnuccctot', pr%mnuccctot)
    if (allocated(pr%mnuccttot)) call put2(ncid, 'proc_rates_mnuccttot', pr%mnuccttot)
    if (allocated(pr%msacwitot)) call put2(ncid, 'proc_rates_msacwitot', pr%msacwitot)
    if (allocated(pr%psacwstot)) call put2(ncid, 'proc_rates_psacwstot', pr%psacwstot)
    if (allocated(pr%bergstot)) call put2(ncid, 'proc_rates_bergstot', pr%bergstot)
    if (allocated(pr%vapdepstot)) call put2(ncid, 'proc_rates_vapdepstot', pr%vapdepstot)
    if (allocated(pr%bergtot)) call put2(ncid, 'proc_rates_bergtot', pr%bergtot)
    if (allocated(pr%melttot)) call put2(ncid, 'proc_rates_melttot', pr%melttot)
    if (allocated(pr%meltstot)) call put2(ncid, 'proc_rates_meltstot', pr%meltstot)
    if (allocated(pr%meltgtot)) call put2(ncid, 'proc_rates_meltgtot', pr%meltgtot)
    if (allocated(pr%homotot)) call put2(ncid, 'proc_rates_homotot', pr%homotot)
    if (allocated(pr%qcrestot)) call put2(ncid, 'proc_rates_qcrestot', pr%qcrestot)
    if (allocated(pr%prcitot)) call put2(ncid, 'proc_rates_prcitot', pr%prcitot)
    if (allocated(pr%praitot)) call put2(ncid, 'proc_rates_praitot', pr%praitot)
    if (allocated(pr%qirestot)) call put2(ncid, 'proc_rates_qirestot', pr%qirestot)
    if (allocated(pr%mnuccrtot)) call put2(ncid, 'proc_rates_mnuccrtot', pr%mnuccrtot)
    if (allocated(pr%mnudeptot)) call put2(ncid, 'proc_rates_mnudeptot', pr%mnudeptot)
    if (allocated(pr%mnuccritot)) call put2(ncid, 'proc_rates_mnuccritot', pr%mnuccritot)
    if (allocated(pr%pracstot)) call put2(ncid, 'proc_rates_pracstot', pr%pracstot)
    if (allocated(pr%meltsdttot)) call put2(ncid, 'proc_rates_meltsdttot', pr%meltsdttot)
    if (allocated(pr%frzrdttot)) call put2(ncid, 'proc_rates_frzrdttot', pr%frzrdttot)
    if (allocated(pr%mnuccdtot)) call put2(ncid, 'proc_rates_mnuccdtot', pr%mnuccdtot)
    if (allocated(pr%pracgtot)) call put2(ncid, 'proc_rates_pracgtot', pr%pracgtot)
    if (allocated(pr%psacwgtot)) call put2(ncid, 'proc_rates_psacwgtot', pr%psacwgtot)
    if (allocated(pr%pgsacwtot)) call put2(ncid, 'proc_rates_pgsacwtot', pr%pgsacwtot)
    if (allocated(pr%pgracstot)) call put2(ncid, 'proc_rates_pgracstot', pr%pgracstot)
    if (allocated(pr%prdgtot)) call put2(ncid, 'proc_rates_prdgtot', pr%prdgtot)
    if (allocated(pr%qmultgtot)) call put2(ncid, 'proc_rates_qmultgtot', pr%qmultgtot)
    if (allocated(pr%qmultrgtot)) call put2(ncid, 'proc_rates_qmultrgtot', pr%qmultrgtot)
    if (allocated(pr%psacrtot)) call put2(ncid, 'proc_rates_psacrtot', pr%psacrtot)
    if (allocated(pr%npracgtot)) call put2(ncid, 'proc_rates_npracgtot', pr%npracgtot)
    if (allocated(pr%nscngtot)) call put2(ncid, 'proc_rates_nscngtot', pr%nscngtot)
    if (allocated(pr%ngracstot)) call put2(ncid, 'proc_rates_ngracstot', pr%ngracstot)
    if (allocated(pr%nmultgtot)) call put2(ncid, 'proc_rates_nmultgtot', pr%nmultgtot)
    if (allocated(pr%nmultrgtot)) call put2(ncid, 'proc_rates_nmultrgtot', pr%nmultrgtot)
    if (allocated(pr%npsacwgtot)) call put2(ncid, 'proc_rates_npsacwgtot', pr%npsacwgtot)
    if (allocated(pr%nnuccctot)) call put2(ncid, 'proc_rates_nnuccctot', pr%nnuccctot)
    if (allocated(pr%nnuccttot)) call put2(ncid, 'proc_rates_nnuccttot', pr%nnuccttot)
    if (allocated(pr%nnuccdtot)) call put2(ncid, 'proc_rates_nnuccdtot', pr%nnuccdtot)
    if (allocated(pr%nnudeptot)) call put2(ncid, 'proc_rates_nnudeptot', pr%nnudeptot)
    if (allocated(pr%nhomotot)) call put2(ncid, 'proc_rates_nhomotot', pr%nhomotot)
    if (allocated(pr%nnuccrtot)) call put2(ncid, 'proc_rates_nnuccrtot', pr%nnuccrtot)
    if (allocated(pr%nnuccritot)) call put2(ncid, 'proc_rates_nnuccritot', pr%nnuccritot)
    if (allocated(pr%nsacwitot)) call put2(ncid, 'proc_rates_nsacwitot', pr%nsacwitot)
    if (allocated(pr%npratot)) call put2(ncid, 'proc_rates_npratot', pr%npratot)
    if (allocated(pr%npsacwstot)) call put2(ncid, 'proc_rates_npsacwstot', pr%npsacwstot)
    if (allocated(pr%npraitot)) call put2(ncid, 'proc_rates_npraitot', pr%npraitot)
    if (allocated(pr%npracstot)) call put2(ncid, 'proc_rates_npracstot', pr%npracstot)
    if (allocated(pr%nprctot)) call put2(ncid, 'proc_rates_nprctot', pr%nprctot)
    if (allocated(pr%nprcitot)) call put2(ncid, 'proc_rates_nprcitot', pr%nprcitot)
    if (allocated(pr%ncsedten)) call put2(ncid, 'proc_rates_ncsedten', pr%ncsedten)
    if (allocated(pr%nisedten)) call put2(ncid, 'proc_rates_nisedten', pr%nisedten)
    if (allocated(pr%nrsedten)) call put2(ncid, 'proc_rates_nrsedten', pr%nrsedten)
    if (allocated(pr%nssedten)) call put2(ncid, 'proc_rates_nssedten', pr%nssedten)
    if (allocated(pr%ngsedten)) call put2(ncid, 'proc_rates_ngsedten', pr%ngsedten)
    if (allocated(pr%nmelttot)) call put2(ncid, 'proc_rates_nmelttot', pr%nmelttot)
    if (allocated(pr%nmeltstot)) call put2(ncid, 'proc_rates_nmeltstot', pr%nmeltstot)
    if (allocated(pr%nmeltgtot)) call put2(ncid, 'proc_rates_nmeltgtot', pr%nmeltgtot)
    if (allocated(pr%nraggtot)) call put2(ncid, 'proc_rates_nraggtot', pr%nraggtot)
    if (allocated(pr%pgam_out)) call put2(ncid, 'proc_rates_pgam_out', pr%pgam_out)
    if (allocated(pr%lamc_out)) call put2(ncid, 'proc_rates_lamc_out', pr%lamc_out)
    if (allocated(pr%lamr_out)) call put2(ncid, 'proc_rates_lamr_out', pr%lamr_out)
    if (allocated(pr%n0r_out)) call put2(ncid, 'proc_rates_n0r_out', pr%n0r_out)
    if (allocated(pr%scale_qc)) call put2(ncid, 'proc_rates_scale_qc', pr%scale_qc)
    if (allocated(pr%scale_nc)) call put2(ncid, 'proc_rates_scale_nc', pr%scale_nc)
    if (allocated(pr%scale_qr)) call put2(ncid, 'proc_rates_scale_qr', pr%scale_qr)
    if (allocated(pr%scale_nr)) call put2(ncid, 'proc_rates_scale_nr', pr%scale_nr)
    if (allocated(pr%qc_out_TAU)) call put2(ncid, 'proc_rates_qc_out_TAU', pr%qc_out_TAU)
    if (allocated(pr%nc_out_TAU)) call put2(ncid, 'proc_rates_nc_out_TAU', pr%nc_out_TAU)
    if (allocated(pr%qr_out_TAU)) call put2(ncid, 'proc_rates_qr_out_TAU', pr%qr_out_TAU)
    if (allocated(pr%nr_out_TAU)) call put2(ncid, 'proc_rates_nr_out_TAU', pr%nr_out_TAU)
    if (allocated(pr%qc_in_TAU)) call put2(ncid, 'proc_rates_qc_in_TAU', pr%qc_in_TAU)
    if (allocated(pr%nc_in_TAU)) call put2(ncid, 'proc_rates_nc_in_TAU', pr%nc_in_TAU)
    if (allocated(pr%qr_in_TAU)) call put2(ncid, 'proc_rates_qr_in_TAU', pr%qr_in_TAU)
    if (allocated(pr%nr_in_TAU)) call put2(ncid, 'proc_rates_nr_in_TAU', pr%nr_in_TAU)
    if (allocated(pr%qctend_KK2000)) call put2(ncid, 'proc_rates_qctend_KK2000', pr%qctend_KK2000)
    if (allocated(pr%nctend_KK2000)) call put2(ncid, 'proc_rates_nctend_KK2000', pr%nctend_KK2000)
    if (allocated(pr%qrtend_KK2000)) call put2(ncid, 'proc_rates_qrtend_KK2000', pr%qrtend_KK2000)
    if (allocated(pr%nrtend_KK2000)) call put2(ncid, 'proc_rates_nrtend_KK2000', pr%nrtend_KK2000)
    if (allocated(pr%qctend_SB2001)) call put2(ncid, 'proc_rates_qctend_SB2001', pr%qctend_SB2001)
    if (allocated(pr%nctend_SB2001)) call put2(ncid, 'proc_rates_nctend_SB2001', pr%nctend_SB2001)
    if (allocated(pr%qrtend_SB2001)) call put2(ncid, 'proc_rates_qrtend_SB2001', pr%qrtend_SB2001)
    if (allocated(pr%nrtend_SB2001)) call put2(ncid, 'proc_rates_nrtend_SB2001', pr%nrtend_SB2001)
    if (allocated(pr%qctend_TAU)) call put2(ncid, 'proc_rates_qctend_TAU', pr%qctend_TAU)
    if (allocated(pr%nctend_TAU)) call put2(ncid, 'proc_rates_nctend_TAU', pr%nctend_TAU)
    if (allocated(pr%qrtend_TAU)) call put2(ncid, 'proc_rates_qrtend_TAU', pr%qrtend_TAU)
    if (allocated(pr%nrtend_TAU)) call put2(ncid, 'proc_rates_nrtend_TAU', pr%nrtend_TAU)
    if (allocated(pr%gmnnn_lmnnn_TAU)) call put2(ncid, 'proc_rates_gmnnn_lmnnn_TAU', pr%gmnnn_lmnnn_TAU)
    if (allocated(pr%ML_fixer)) call put2(ncid, 'proc_rates_ML_fixer', pr%ML_fixer)
    if (allocated(pr%QC_fixer)) call put2(ncid, 'proc_rates_QC_fixer', pr%QC_fixer)
    if (allocated(pr%NC_fixer)) call put2(ncid, 'proc_rates_NC_fixer', pr%NC_fixer)
    if (allocated(pr%QR_fixer)) call put2(ncid, 'proc_rates_QR_fixer', pr%QR_fixer)
    if (allocated(pr%NR_fixer)) call put2(ncid, 'proc_rates_NR_fixer', pr%NR_fixer)
    call nc(nf90_close(ncid), 'close '//trim(file))
  end subroutine write_outputs

  ! 🎛️📄 Dump control: active only when the file pumas_dump.nml exists in the run directory.
  subroutine read_dump_control(active, dump_dir, dump_nsteps, chunk_stride)
    ! 🗂️⚡ Reads pumas_dump.nml once per process and caches the answer.
    logical, intent(out) :: active
    character(len=512), intent(out) :: dump_dir
    integer, intent(out) :: dump_nsteps(16), chunk_stride
    logical, save :: cached = .false., c_active = .false.
    character(len=512), save :: c_dir = '.'
    integer, save :: c_nsteps(16) = -1, c_stride = 1
    integer :: unit, ios
    namelist /pumas_dump/ dump_dir, dump_nsteps, chunk_stride
    if (.not. cached) then
      dump_dir = '.'
      dump_nsteps = -1
      chunk_stride = 1
      inquire(file='pumas_dump.nml', exist=c_active)
      if (c_active) then
        open(newunit=unit, file='pumas_dump.nml', status='old', action='read', iostat=ios)
        if (ios == 0) read(unit, nml=pumas_dump, iostat=ios)
        close(unit)
        c_active = (ios == 0)
      end if
      c_dir = dump_dir
      c_nsteps = dump_nsteps
      c_stride = max(chunk_stride, 1)
      cached = .true.
    end if
    active = c_active
    dump_dir = c_dir
    dump_nsteps = c_nsteps
    chunk_stride = c_stride
  end subroutine read_dump_control

  integer function create(file, deltat, nlev) result(ncid)
    character(len=*), intent(in) :: file
    real(r8), intent(in) :: deltat
    integer, intent(in) :: nlev
    call nc(nf90_create(trim(file), ior(nf90_clobber, nf90_netcdf4), ncid), 'create '//trim(file))
    call nc(nf90_put_att(ncid, nf90_global, 'deltat', deltat), 'deltat')
    call nc(nf90_put_att(ncid, nf90_global, 'nlev', nlev), 'nlev')
  end function create

  subroutine nc(status, what)
    integer, intent(in) :: status
    character(len=*), intent(in) :: what
    if (status /= nf90_noerr) then
       write(0,'(a)') 'pumas_column_io: '//trim(what)//': '//trim(nf90_strerror(status))
       error stop 1
    end if
  end subroutine nc

  integer function def_dim(ncid, name, n)
    integer, intent(in) :: ncid, n
    character(len=*), intent(in) :: name
    if (nf90_inq_dimid(ncid, name, def_dim) == nf90_noerr) return
    call nc(nf90_redef(ncid), 'redef')
    call nc(nf90_def_dim(ncid, name, n, def_dim), 'define dimension '//name)
    call nc(nf90_enddef(ncid), 'enddef')
  end function def_dim

  integer function level_dim(ncid, n)
    integer, intent(in) :: ncid, n
    integer :: nlev
    call nc(nf90_get_att(ncid, nf90_global, 'nlev', nlev), 'nlev')
    if (n == nlev) then
       level_dim = def_dim(ncid, 'lev', n)
    else
       level_dim = def_dim(ncid, 'ilev', n)
    end if
  end function level_dim

  integer function def_var(ncid, name, dims) result(varid)
    integer, intent(in) :: ncid, dims(:)
    character(len=*), intent(in) :: name
    call nc(nf90_redef(ncid), 'redef')
    call nc(nf90_def_var(ncid, name, nf90_double, dims, varid), 'define '//name)
    call nc(nf90_enddef(ncid), 'enddef')
  end function def_var

  subroutine put1(ncid, name, arr)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: arr(:)
    integer :: varid
    varid = def_var(ncid, name, [def_dim(ncid, 'col', size(arr))])
    call nc(nf90_put_var(ncid, varid, arr), 'write '//name)
  end subroutine put1

  subroutine put2(ncid, name, arr)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: arr(:,:)
    integer :: varid
    varid = def_var(ncid, name, [def_dim(ncid, 'col', size(arr, 1)), level_dim(ncid, size(arr, 2))])
    call nc(nf90_put_var(ncid, varid, arr), 'write '//name)
  end subroutine put2

  subroutine put3(ncid, name, arr)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: arr(:,:,:)
    integer :: varid
    varid = def_var(ncid, name, [def_dim(ncid, 'col', size(arr, 1)), &
         level_dim(ncid, size(arr, 2)), def_dim(ncid, 'dust_bin', size(arr, 3))])
    call nc(nf90_put_var(ncid, varid, arr), 'write '//name)
  end subroutine put3

end module pumas_column_io
