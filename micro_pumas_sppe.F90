! 🧊🎛️ Sppe process controls for PUMAS (oPUMAS addition, not in ESCOMP/PUMAS).
! 🔕✅ Every default reproduces upstream PUMAS bit-for-bit; a host changes them once,
!    after micro_pumas_init, through micro_pumas_sppe_set.
! 📘🔗 Origin: github.com/diegovillanueva/Sppe_revision_public (CESM 2.1.5 MG2 patches).

module micro_pumas_sppe

   use pumas_kinds, only: r8 => kind_r8

   implicit none
   private
   save

   public :: micro_pumas_sppe_set

   ! 🧊🚧 .false.: deposition nucleation to naai only below 233.15 K (no cirrus leak into mixed-phase clouds)
   logical,  public :: naai_het_also_in_mpc = .true.
   ! 🧊🔓 .false.: remove the micro_mg_max_nicons cap on the ice number tendency
   logical,  public :: enable_nimax = .true.
   ! ☁️💧 .true.: condensate from removed supersaturation goes to liquid in the 238-268 K ramp
   logical,  public :: detrainment_ramp_liq = .false.
   ! 💧🛑 .false.: no contact freezing of cloud droplets (utils contact_freezing)
   logical,  public :: onoff_contact_freezing = .true.
   ! 🧊🛑 .false.: no Hallett-Mossop splinters from snow riming
   logical,  public :: onoff_sip = .true.
   ! 🌧️🛑 .false.: no Bigg heterogeneous rain freezing
   logical,  public :: rain_freeze = .true.
   ! 🎛️🧮 Process multipliers
   real(r8), public :: auto_liq_factor = 1._r8               ! liquid autoconversion (prc, nprc, nprc1), kk2000 and sb2001
   real(r8), public :: riming_factor = 1._r8                 ! snow riming (psacws, npsacws), before Hallett-Mossop
   real(r8), public :: micro_mg_ice_berg_ice_factor = 1._r8  ! Bergeron deposition onto cloud ice (berg)
   real(r8), public :: micro_mg_ice_berg_snow_factor = 1._r8 ! Bergeron deposition onto snow (bergs)
   real(r8), public :: factor_HM = 1._r8                     ! Hallett-Mossop splinter number

   !$acc declare create (naai_het_also_in_mpc, enable_nimax, detrainment_ramp_liq, &
   !$acc                 onoff_contact_freezing, onoff_sip, rain_freeze,           &
   !$acc                 auto_liq_factor, riming_factor,                           &
   !$acc                 micro_mg_ice_berg_ice_factor, micro_mg_ice_berg_snow_factor, factor_HM)

contains

subroutine micro_pumas_sppe_set(naai_het_also_in_mpc_in, enable_nimax_in, detrainment_ramp_liq_in, &
     onoff_contact_freezing_in, onoff_sip_in, rain_freeze_in, auto_liq_factor_in, riming_factor_in, &
     micro_mg_ice_berg_ice_factor_in, micro_mg_ice_berg_snow_factor_in, factor_HM_in)

   ! 🎛️📥 Set any subset of the controls; arguments left out keep their current value.
   logical,  intent(in), optional :: naai_het_also_in_mpc_in, enable_nimax_in, detrainment_ramp_liq_in
   logical,  intent(in), optional :: onoff_contact_freezing_in, onoff_sip_in, rain_freeze_in
   real(r8), intent(in), optional :: auto_liq_factor_in, riming_factor_in
   real(r8), intent(in), optional :: micro_mg_ice_berg_ice_factor_in, micro_mg_ice_berg_snow_factor_in, factor_HM_in

   if (present(naai_het_also_in_mpc_in)) naai_het_also_in_mpc = naai_het_also_in_mpc_in
   if (present(enable_nimax_in)) enable_nimax = enable_nimax_in
   if (present(detrainment_ramp_liq_in)) detrainment_ramp_liq = detrainment_ramp_liq_in
   if (present(onoff_contact_freezing_in)) onoff_contact_freezing = onoff_contact_freezing_in
   if (present(onoff_sip_in)) onoff_sip = onoff_sip_in
   if (present(rain_freeze_in)) rain_freeze = rain_freeze_in
   if (present(auto_liq_factor_in)) auto_liq_factor = auto_liq_factor_in
   if (present(riming_factor_in)) riming_factor = riming_factor_in
   if (present(micro_mg_ice_berg_ice_factor_in)) micro_mg_ice_berg_ice_factor = micro_mg_ice_berg_ice_factor_in
   if (present(micro_mg_ice_berg_snow_factor_in)) micro_mg_ice_berg_snow_factor = micro_mg_ice_berg_snow_factor_in
   if (present(factor_HM_in)) factor_HM = factor_HM_in

   !$acc update device (naai_het_also_in_mpc, enable_nimax, detrainment_ramp_liq, &
   !$acc                onoff_contact_freezing, onoff_sip, rain_freeze,           &
   !$acc                auto_liq_factor, riming_factor,                           &
   !$acc                micro_mg_ice_berg_ice_factor, micro_mg_ice_berg_snow_factor, factor_HM)

end subroutine micro_pumas_sppe_set

end module micro_pumas_sppe
