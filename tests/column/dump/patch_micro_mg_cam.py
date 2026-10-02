#!/usr/bin/env python3
"""🧪📤 Insert the column dump into a copy of interface/cesm215/micro_mg_cam.F90.

🎛️🔕 The dump is active only when the run directory holds pumas_dump.nml; the patched file
lives only in the dump source directory, never in interface/cesm215 (golden rule).
Usage: patch_micro_mg_cam.py <in micro_mg_cam.F90> <out micro_mg_cam.F90>
"""

import sys
from pathlib import Path

USES = """
   ! >>> oPUMAS column dump (tests/column/dump) >>>
   use time_manager,    only: get_nstep
   use pumas_column_io, only: write_inputs, write_outputs, read_dump_control
   ! <<< oPUMAS column dump <<<
"""

DECLS = """
   ! >>> oPUMAS column dump (tests/column/dump) >>>
   logical :: dump_on
   character(len=512) :: dump_dir
   character(len=600) :: dump_file
   integer :: dump_nsteps(16), dump_stride
   ! <<< oPUMAS column dump <<<
"""

CALL = """
            ! >>> oPUMAS column dump (tests/column/dump) >>>
            call read_dump_control(dump_on, dump_dir, dump_nsteps, dump_stride)
            if (dump_on .and. it == 1 .and. any(dump_nsteps == get_nstep()) .and. &
                 mod(state%lchnk, dump_stride) == 0) then
               write(dump_file, '(a,"/inputs.n",i5.5,".c",i6.6,".nc")') &
                    trim(dump_dir), get_nstep(), state%lchnk
               call write_inputs(dump_file, dtime/num_steps, &
                    packed_t, packed_q, packed_qc, packed_qi, packed_nc, packed_ni, &
                    packed_qr, packed_qs, packed_nr, packed_ns, packed_qg, packed_ng, &
                    packed_relvar, packed_accre_enhan, packed_p, packed_pdel, packed_pint, &
                    packed_cldn, packed_liqcldf, packed_icecldf, packed_qsatfac, &
                    packed_naai, packed_npccn, packed_rndst, packed_nacon, &
                    packed_tnd_qsnow, packed_tnd_nsnow, packed_re_ice, &
                    packed_frzimm, packed_frzcnt, packed_frzdep, &
                    lat=state%lat(mgcols(:mgncol))*180._r8/3.14159265358979323846_r8, &
                    lon=state%lon(mgcols(:mgncol))*180._r8/3.14159265358979323846_r8, &
                    nstep=get_nstep(), lchnk=state%lchnk)
               write(dump_file, '(a,"/outputs.n",i5.5,".c",i6.6,".nc")') &
                    trim(dump_dir), get_nstep(), state%lchnk
               call write_outputs(dump_file, dtime/num_steps, &
                    packed_rate1ord_cw2pr_st, packed_tlat, packed_qvlat, &
                    packed_qctend, packed_qitend, packed_nctend, packed_nitend, &
                    packed_qrtend, packed_qstend, packed_nrtend, packed_nstend, &
                    packed_qgtend, packed_ngtend, packed_rel, rel_fn_dum, packed_rei, &
                    packed_sadice, packed_sadsnow, packed_nevapr, packed_am_evp_st, &
                    packed_prain, packed_cmeout, packed_dei, packed_mu, packed_lambdac, &
                    packed_qsout, packed_des, pumas_dum(:,:,3), pumas_dum(:,:,4), &
                    pumas_dum(:,:,5), packed_qrout, reff_rain_dum, reff_snow_dum, &
                    pumas_dum(:,:,6), packed_nrout, packed_nsout, packed_refl, &
                    packed_arefl, packed_areflz, packed_frefl, packed_csrfl, &
                    packed_acsrfl, packed_fcsrfl, pumas_dum(:,:,7), pumas_dum(:,:,8), &
                    packed_rercld, packed_ncai, packed_ncal, packed_qrout2, &
                    packed_qsout2, packed_nrout2, packed_nsout2, drout_dum, dsout2_dum, &
                    packed_qgout2, packed_ngout2, pumas_dum(:,:,9), pumas_dum(:,:,10), &
                    packed_freqs, packed_freqr, packed_nfice, packed_qcrat, &
                    packed_prer_evap, packed_cflx, packed_iflx, packed_gflx, &
                    packed_rflx, packed_sflx, packed_prect, packed_preci, proc_rates)
            end if
            ! <<< oPUMAS column dump <<<
"""

SETTINGS_CALL = """
   ! >>> oPUMAS column dump (tests/column/dump) >>>
   if (micro_mg_version == 3) call pumas_dump_settings()
   ! <<< oPUMAS column dump <<<
"""

SETTINGS_SUB = """
! >>> oPUMAS column dump (tests/column/dump) >>>
! 🎛️📄 Write the micro_pumas_init arguments as namelist pumas_settings for the column driver.
subroutine pumas_dump_settings()
   use pumas_column_io, only: read_dump_control
   logical :: dump_on
   character(len=512) :: dump_dir
   integer :: dump_nsteps(16), dump_stride, unit
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
   if (.not. masterproc) return
   call read_dump_control(dump_on, dump_dir, dump_nsteps, dump_stride)
   if (.not. dump_on) return
   open(newunit=unit, file=trim(dump_dir)//'/settings.nml', status='replace', action='write')
   write(unit, nml=pumas_settings)
   close(unit)
end subroutine pumas_dump_settings
! <<< oPUMAS column dump <<<
"""


def insert_after(text: str, anchor: str, block: str, start: int = 0) -> str:
    i = text.index(anchor, start)
    j = text.index("\n", i + len(anchor) - 1) + 1
    return text[:j] + block + text[j:]


def insert_before(text: str, anchor: str, block: str, start: int = 0) -> str:
    i = text.index(anchor, start)
    i = text.rindex("\n", 0, i) + 1
    return text[:i] + block.lstrip("\n") + text[i:]


def main() -> None:
    src, dst = Path(sys.argv[1]), Path(sys.argv[2])
    s = src.read_text()
    tp = s.index("subroutine micro_mg_cam_tend_pack(")
    s = insert_after(s, "use pumas_stochastic_collect_tau, only: ncd", USES, tp)
    tp = s.index("subroutine micro_mg_cam_tend_pack(")
    # 🧾 Declarations go right before the first local declaration after the uses.
    s = insert_before(s, "type(physics_state),         intent(in)    :: state", DECLS,
                      s.index(USES, tp))
    s = insert_before(s, "! Copy PUMAS process rates into the MG2 output arrays", CALL)
    s = insert_after(s, 'call handle_errmsg(errstring, subname="micro_mg_init")', SETTINGS_CALL)
    s = insert_before(s, "end module micro_mg_cam", SETTINGS_SUB)
    dst.write_text(s)
    print(dst.resolve())


if __name__ == "__main__":
    main()
