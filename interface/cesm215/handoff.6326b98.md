# 🤝🧊 Handoff 6326b98

- **Commit**: 6326b98 ✨ feat: 🧊🔌 interface/cesm215 to run oPUMAS in CESM 2.1.5 CAM (micro_mg_version=3)
- **What was done**: 🧊🔌 PUMAS runs in CESM 2.1.5 through `-usr_src` of this dir. Related commits: CAM worktree `dv_pumas_rel60` 81285e66 (namelist XML), config_and_bash 383f787 + db43f9d (`CAM_USR_SRC`, `--srcroot`), cesm_wrapper 95d5a8e (README).
- **Current state**: ✅🏗️ Case `cesm215_pumas.F2000climo.pumas_test` builds; 1-day test submitted as job 15714706 (result not yet checked).
- **Next steps**: 🔍📈 check the 1-day run (`cesm.log`, `atm.log` for "PUMAS namelist"), then a 3-month run vs `cesm215.F2000climo.Sppe_ndg` (PRECT, CLDLIQ, CLDICE, TGCLDLWP, SWCF, LWCF). 🔄 Port to `dv_diag` by applying `./make_interface.sh diff` to that branch's `micro_mg_cam.F90`.
- **Key decisions**: 📄 full override files instead of generated patches (base is an immutable tag); graupel off; `emulated` warm rain not wired; MG2 `micro_mg_utils_init` still called for CAM post-processing.

## 🧮✅ Standalone compile check (herdr-compute, never a login node)

```bash
E=<any cesm215 EXEROOT>; O=$E/atm/obj; B=<scratch dir>; cd $B
LINE=$(zcat $(ls -t $E/atm.bldlog.*.gz | head -1) | grep mpif90 | grep micro_mg_cam.F90 | head -1); LINE=${LINE% *}
LINE=$(echo "$LINE" | sed -e "s# -I\.\./# -I$O/../#g" -e "s# -I\. # -I$O #g" -e "s#^mpif90  -c#mpif90 -c -I$B -J$B -fallow-argument-mismatch#")
for f in pumas_kinds pumas_gamma_function micro_pumas_utils micro_pumas_diags pumas_stochastic_collect_tau \
         module_neural_net tau_neural_net_quantile ML_fixer_check wv_sat_methods micro_pumas_v1 micro_mg_cam; do
  $LINE <interface dir>/$f.F90 -o $B/$f.o; done
```
