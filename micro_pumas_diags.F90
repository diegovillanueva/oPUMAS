module micro_pumas_diags

!----------------------------------------
! PUMAS diagnostics support package
!----------------------------------------

use pumas_kinds,   only: r8=>kind_r8

!> \section arg_table_proc_rates_type  Argument Table
!! \htmlinclude proc_rates_type.html
  type, public :: proc_rates_type

    real(r8), allocatable :: snow_mass_tend_production_kg_per_kg_per_s(:,:)     ! production of snow (kg/kg/s)
    real(r8), allocatable :: snow_mass_tend_sublimation_kg_per_kg_per_s(:,:)     ! sublimation rate of snow (1/s)
    real(r8), allocatable :: cldliq_mass_tend_sedimentation_evaporation_kg_per_kg_per_s(:,:)      ! cloud water evaporation due to sedimentation (1/s)
    real(r8), allocatable :: cldice_mass_tend_sedimentation_sublimation_kg_per_kg_per_s(:,:)      ! cloud ice sublimation due to sublimation (1/s)
    real(r8), allocatable :: vapor_mass_tend_residual_condensation_kg_per_kg_per_s(:,:)        ! residual condensation term to ensure RH < 100% (1/s)
    real(r8), allocatable :: cldice_mass_tend_net_deposition_kg_per_kg_per_s(:,:)      ! grid-mean cloud ice sub/dep (1/s)
    real(r8), allocatable :: cldliq_mass_weighted_fall_speed_m_per_s(:,:)        ! mass-weighted cloud water fallspeed (m/s)
    real(r8), allocatable :: cldice_mass_weighted_fall_speed_m_per_s(:,:)        ! mass-weighted cloud ice fallspeed (m/s)
    real(r8), allocatable :: rain_mass_weighted_fall_speed_m_per_s(:,:)          ! mass weighted rain fallspeed (m/s)
    real(r8), allocatable :: snow_mass_weighted_fall_speed_m_per_s(:,:)          ! mass weighted snow fallspeed (m/s)
    real(r8), allocatable :: graupel_mass_weighted_fall_speed_m_per_s(:,:)          ! mass weighted graupel/hail fallspeed (m/s)
    real(r8), allocatable :: graupel_mass_tend_sedimentation_kg_per_kg_per_s(:,:)     ! qg sedimentation tendency (1/s)
    real(r8), allocatable :: cldliq_mass_tend_sedimentation_kg_per_kg_per_s(:,:)     ! qc sedimentation tendency (1/s)
    real(r8), allocatable :: cldice_mass_tend_sedimentation_kg_per_kg_per_s(:,:)     ! qi sedimentation tendency (1/s)
    real(r8), allocatable :: rain_mass_tend_sedimentation_kg_per_kg_per_s(:,:)     ! qr sedimentation tendency (1/s)
    real(r8), allocatable :: snow_mass_tend_sedimentation_kg_per_kg_per_s(:,:)     ! qs sedimentation tendency (1/s)

    real(r8), allocatable :: cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: snow_mass_tend_deposition_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: cldice_mass_tend_melting_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: snow_mass_tend_melting_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: graupel_mass_tend_melting_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: cldliq_mass_tend_homogeneous_freezing_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: cldliq_mass_tend_residual_condensation_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: cldice_mass_tend_residual_deposition_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: air_heating_rate_snow_melting_W_per_kg(:,:)
    real(r8), allocatable :: air_heating_rate_rain_freezing_W_per_kg(:,:)
    real(r8), allocatable :: cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: graupel_mass_tend_from_snow_riming_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: graupel_mass_tend_deposition_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: snow_mass_tend_collection_by_rain_kg_per_kg_per_s(:,:)
    real(r8), allocatable :: rain_number_tend_accretion_by_graupel_per_kg_per_s(:,:)
    real(r8), allocatable :: graupel_number_tend_from_snow_riming_per_kg_per_s(:,:)
    real(r8), allocatable :: graupel_number_tend_rain_snow_collection_per_kg_per_s(:,:)
    real(r8), allocatable :: cldice_number_tend_splintering_by_graupel_per_kg_per_s(:,:)
    real(r8), allocatable :: cldice_number_tend_splintering_rain_graupel_per_kg_per_s(:,:)
    real(r8), allocatable :: cldliq_number_tend_riming_by_graupel_per_kg_per_s(:,:)

  real(r8), allocatable :: cldliq_number_tend_immersion_freezing_per_kg_per_s(:,:)        ! change n  due to Immersion freezing of cloud water
  real(r8), allocatable :: cldliq_number_tend_contact_freezing_per_kg_per_s(:,:)        ! change n  due to Contact freezing of cloud water
  real(r8), allocatable :: cldice_number_tend_nucleation_from_vapor_per_kg_per_s(:,:)        ! change n  due to Ice nucleation
  real(r8), allocatable :: cldliq_number_tend_deposition_nucleation_per_kg_per_s(:,:)        ! change n  due to Deposition Nucleation
  real(r8), allocatable :: cldliq_number_tend_homogeneous_freezing_per_kg_per_s(:,:)         ! change n  due to Homogeneous freezing of cloud water
  real(r8), allocatable :: rain_number_tend_freezing_to_snow_per_kg_per_s(:,:)        ! change n  due to heterogeneous freezing of rain to snow (1/s)
  real(r8), allocatable :: rain_number_tend_freezing_to_cldice_per_kg_per_s(:,:)       ! change n  due to Heterogeneous freezing of rain to ice
  real(r8), allocatable :: cldice_number_tend_splintering_by_snow_per_kg_per_s(:,:)        ! change n  due to Conversion of cloud water [to cloud ice]
                                                         !                  from rime-splintering
  real(r8), allocatable :: cldliq_number_tend_accretion_by_rain_per_kg_per_s(:,:)          ! change n  due to Accretion of cloud water by rain
  real(r8), allocatable :: cldliq_number_tend_riming_by_snow_per_kg_per_s(:,:)       ! change n  due to Accretion of cloud water by snow
  real(r8), allocatable :: cldice_number_tend_accretion_by_snow_per_kg_per_s(:,:)         ! change n  due to Accretion of cloud ice to snow
  real(r8), allocatable :: rain_number_tend_accretion_by_snow_per_kg_per_s(:,:)        ! change n  due to Accretion of rain by snow
  real(r8), allocatable :: rain_number_tend_autoconversion_from_cldliq_per_kg_per_s(:,:)          ! change nr  due to Autoconversion of cloud water [to rain]
  real(r8), allocatable :: cldice_number_tend_autoconversion_to_snow_per_kg_per_s(:,:)         ! change n  due to Autoconversion of cloud ice to snow
  real(r8), allocatable :: cldliq_number_tend_sedimentation_per_kg_per_s(:,:)         ! change n  due to cloud liquid sedimentation
  real(r8), allocatable :: cldice_number_tend_sedimentation_per_kg_per_s(:,:)         ! change n  due to cloud ice sedimentation
  real(r8), allocatable :: rain_number_tend_sedimentation_per_kg_per_s(:,:)         ! change n  due to rain sedimentation
  real(r8), allocatable :: snow_number_tend_sedimentation_per_kg_per_s(:,:)         ! change n  due to snow sedimentation
  real(r8), allocatable :: graupel_number_tend_sedimentation_per_kg_per_s(:,:)         ! change n  due to graupel sedimentation
  real(r8), allocatable :: cldice_number_tend_melting_per_kg_per_s(:,:)         ! change n  due to Melting of cloud ice
  real(r8), allocatable :: snow_number_tend_melting_per_kg_per_s(:,:)        ! change n  due to Melting of snow
  real(r8), allocatable :: graupel_number_tend_melting_per_kg_per_s(:,:)        ! change n  due to Melting of graupel

  ! TAU diagnostic variables
  real(r8), allocatable :: rain_number_tend_self_collection_per_kg_per_s(:,:)          ! change nr  due to self collection of rain


  real(r8), allocatable :: cldliq_size_dist_shape_output(:,:)      ! Liquid Size distribution parameter Mu for output
  real(r8), allocatable :: cldliq_size_dist_slope_output_per_m(:,:)      ! Liquid Size distribution parameter Lambda for output
  real(r8), allocatable :: rain_size_dist_slope_output_per_m(:,:)      ! Rain Size distribution parameter Lambda for output
  real(r8), allocatable :: rain_size_dist_intercept_output_per_m4(:,:)       ! Size distribution parameter n0 for output
  real(r8), allocatable :: scale_qc(:,:)        !TAU scaling factor for liquid mass to ensure conservation
  real(r8), allocatable :: scale_nc(:,:)        !TAU scaling factor for liquid number to ensure conservation
  real(r8), allocatable :: scale_qr(:,:)        !TAU scaling factor for rain mass to ensure conservation
  real(r8), allocatable :: scale_nr(:,:)        !TAU scaling factor for rain_number to ensure conservation
  real(r8), allocatable :: amk_c(:,:,:)         !TAU cloud liquid mass from bins
  real(r8), allocatable :: ank_c(:,:,:)         !TAU cloud liquid number from bins
  real(r8), allocatable :: amk_r(:,:,:)         !TAU cloud rain mass from bins
  real(r8), allocatable :: ank_r(:,:,:)         !TAU cloud rain number from bins
  real(r8), allocatable :: amk(:,:,:)           !TAU all liquid mass from bins
  real(r8), allocatable :: ank(:,:,:)           !TAU all liquid number from bins
  real(r8), allocatable :: amk_out(:,:,:)       !TAU all liquid number from bins output
  real(r8), allocatable :: ank_out(:,:,:)       !TAU all liquid mass from bins output
  real(r8), allocatable :: qc_out_TAU(:,:)      !TAU: output total cloud liquid mass
  real(r8), allocatable :: nc_out_TAU(:,:)      !TAU: output total cloud liquid number
  real(r8), allocatable :: qr_out_TAU(:,:)      !TAU: output total rain mass
  real(r8), allocatable :: nr_out_TAU(:,:)      !TAU: output total cloud rain number
  real(r8), allocatable :: qc_in_TAU(:,:)       !TAU: input total cloud liquid mass
  real(r8), allocatable :: nc_in_TAU(:,:)       !TAU: input total cloud liquid number
  real(r8), allocatable :: qr_in_TAU(:,:)       !TAU: input total rain mass
  real(r8), allocatable :: nr_in_TAU(:,:)       !TAU: input total cloud rain number
  real(r8), allocatable :: cldliq_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s(:,:)   !cloud liquid mass tendency due to autoconversion  & accretion from KK2000
  real(r8), allocatable :: cldliq_number_tend_autoconv_accretion_KK2000_per_kg_per_s(:,:)   !cloud liquid number tendency due to autoconversion  & accretion from KK2000
  real(r8), allocatable :: rain_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s(:,:)   !rain mass tendency due to autoconversion  & accretion from KK2000
  real(r8), allocatable :: rain_number_tend_autoconv_accretion_KK2000_per_kg_per_s(:,:)   !rain number tendency due to autoconversion   & accretion from KK2000
  real(r8), allocatable :: cldliq_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s(:,:)   !cloud liquid mass tendency due to autoconversion  & accretion from SB2001
  real(r8), allocatable :: cldliq_number_tend_autoconv_accretion_SB2001_per_kg_per_s(:,:)   !cloud liquid number tendency due to autoconversion  & accretion from SB2001
  real(r8), allocatable :: rain_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s(:,:)   !rain mass tendency due to autoconversion  & accretion from SB2001
  real(r8), allocatable :: rain_number_tend_autoconv_accretion_SB2001_per_kg_per_s(:,:)   !rain number tendency due to autoconversion  & accretion from SB2001
  real(r8), allocatable :: cldliq_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s(:,:)      !cloud liquid mass tendency due to autoconversion & accretion from TAU or Emulator code
  real(r8), allocatable :: cldliq_number_tend_autoconv_accretion_TAU_per_kg_per_s(:,:)      !cloud liquid number tendency due to autoconversion & accretion from TAU or Emulator code
  real(r8), allocatable :: rain_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s(:,:)      !rain mass tendency due to autoconversion & accretion from TAU or Emulator code
  real(r8), allocatable :: rain_number_tend_autoconv_accretion_TAU_per_kg_per_s(:,:)      !rain number tendency due to autoconversion & accretion from TAU or Emulatorcode
  real(r8), allocatable :: gmnnn_lmnnn_TAU(:,:) ! TAU sum of mass gain and loss from bin code
  real(r8), allocatable :: ML_fixer(:,:)     !Emulated: frequency of ML fixer is activated
  real(r8), allocatable :: QC_fixer(:,:)     !Emulated: change in cloud liquid mass due to ML fixer
  real(r8), allocatable :: NC_fixer(:,:)     !Emulated: change in cloud number number due to ML fixer
  real(r8), allocatable :: QR_fixer(:,:)     !Emulated: change in rain mass due to ML fixer
  real(r8), allocatable :: NR_fixer(:,:)     !Emulated: change in rain number due to ML fixer

    contains
      procedure :: allocate => proc_rates_allocate
      procedure :: deallocate => proc_rates_deallocate
  end type proc_rates_type

contains

   subroutine proc_rates_allocate(this, psetcols, nlev, ncd, warm_rain, errstring)
   !--------------------------------------------------------------
   ! Routine to allocate the elements of the proc_rates DDT
   !--------------------------------------------------------------

      implicit none

      class(proc_rates_type) :: this

      integer,           intent(in) :: psetcols, nlev
      integer,           intent(in) :: ncd
      character(len=16), intent(in) :: warm_rain            ! 'tau','emulated','sb2001' or 'kk2000'
      character(128),   intent(out) :: errstring

      integer :: ierr

      errstring=' '

      allocate(this%snow_mass_tend_production_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%prodsnow'
      end if
      allocate(this%snow_mass_tend_sublimation_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%evapsnow'
      end if
      allocate(this%cldliq_mass_tend_sedimentation_evaporation_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qcsevap'
      end if
      allocate(this%cldice_mass_tend_sedimentation_sublimation_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qisevap'
      end if
      allocate(this%vapor_mass_tend_residual_condensation_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qvres'
      end if
      allocate(this%cldice_mass_tend_net_deposition_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%cmeitot'
      end if
      allocate(this%cldliq_mass_weighted_fall_speed_m_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%vtrmc'
      end if
      allocate(this%cldice_mass_weighted_fall_speed_m_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%vtrmi'
      end if
      allocate(this%rain_mass_weighted_fall_speed_m_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%umr'
      end if
      allocate(this%snow_mass_weighted_fall_speed_m_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%ums'
      end if
      allocate(this%graupel_mass_weighted_fall_speed_m_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%umg'
      end if
      allocate(this%graupel_mass_tend_sedimentation_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qgsedten'
      end if
      allocate(this%cldliq_mass_tend_sedimentation_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qcsedten'
      end if
      allocate(this%cldice_mass_tend_sedimentation_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qisedten'
      end if
      allocate(this%rain_mass_tend_sedimentation_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qrsedten'
      end if
      allocate(this%snow_mass_tend_sedimentation_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qssedten'
      end if
      allocate(this%cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%pratot'
      end if
      allocate(this%cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%prctot'
      end if
      allocate(this%cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%mnuccctot'
      end if
      allocate(this%cldliq_mass_tend_contact_freezing_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%mnuccttot'
      end if
      allocate(this%cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%msacwitot'
      end if
      allocate(this%cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%psacwstot'
      end if
      allocate(this%cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%bergstot'
      end if
      allocate(this%snow_mass_tend_deposition_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%vapdepstot'
      end if
      allocate(this%cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%bergtot'
      end if
      allocate(this%cldice_mass_tend_melting_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%melttot'
      end if
      allocate(this%snow_mass_tend_melting_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%meltstot'
      end if
      allocate(this%graupel_mass_tend_melting_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%meltgtot'
      end if
      allocate(this%cldliq_mass_tend_homogeneous_freezing_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%homotot'
      end if
      allocate(this%cldliq_mass_tend_residual_condensation_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qcrestot'
      end if
      allocate(this%cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%prcitot'
      end if
      allocate(this%cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%praitot'
      end if
      allocate(this%cldice_mass_tend_residual_deposition_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qirestot'
      end if
      allocate(this%rain_mass_tend_freezing_to_snow_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%mnuccrtot'
      end if
      allocate(this%cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%mnudeptot'
      end if
      allocate(this%rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%mnuccritot'
      end if
      allocate(this%rain_mass_tend_accretion_by_snow_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%pracstot'
      end if
      allocate(this%air_heating_rate_snow_melting_W_per_kg(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%meltsdttot'
      end if
      allocate(this%air_heating_rate_rain_freezing_W_per_kg(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%frzrdttot'
      end if
      allocate(this%cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%mnuccdtot'
      end if
      allocate(this%rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%pracgtot'
      end if
      allocate(this%cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%psacwgtot'
      end if
      allocate(this%graupel_mass_tend_from_snow_riming_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%pgsacwtot'
      end if
      allocate(this%graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%pgracstot'
      end if
      allocate(this%graupel_mass_tend_deposition_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%prdgtot'
      end if
      allocate(this%cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qmultgtot'
      end if
      allocate(this%cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qmultrgtot'
      end if
      allocate(this%snow_mass_tend_collection_by_rain_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%psacrtot'
      end if
      allocate(this%rain_number_tend_accretion_by_graupel_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%npracgtot'
      end if
      allocate(this%graupel_number_tend_from_snow_riming_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nscngtot'
      end if
      allocate(this%graupel_number_tend_rain_snow_collection_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%ngracstot'
      end if
      allocate(this%cldice_number_tend_splintering_by_graupel_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nmultgtot'
      end if
      allocate(this%cldice_number_tend_splintering_rain_graupel_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nmultrgtot'
      end if
      allocate(this%cldliq_number_tend_riming_by_graupel_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%npsacwgtot'
      end if
      allocate(this%cldliq_number_tend_immersion_freezing_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nnuccctot'
      end if
      allocate(this%cldliq_number_tend_contact_freezing_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nnuccttot'
      end if
      allocate(this%cldice_number_tend_nucleation_from_vapor_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nnuccdtot'
      end if
      allocate(this%cldliq_number_tend_deposition_nucleation_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nnudeptot'
      end if
      allocate(this%cldliq_number_tend_homogeneous_freezing_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nhomotot'
      end if
      allocate(this%rain_number_tend_freezing_to_snow_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nnuccrtot'
      end if
      allocate(this%rain_number_tend_freezing_to_cldice_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nnuccritot'
      end if
      allocate(this%cldice_number_tend_splintering_by_snow_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nsacwitot'
      end if
      allocate(this%cldliq_number_tend_accretion_by_rain_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%npratot'
      end if
      allocate(this%cldliq_number_tend_riming_by_snow_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%npsacwstot'
      end if
      allocate(this%cldice_number_tend_accretion_by_snow_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%npraitot'
      end if
      allocate(this%rain_number_tend_accretion_by_snow_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%npracstot'
      end if
      allocate(this%rain_number_tend_autoconversion_from_cldliq_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nprctot'
      end if
      allocate(this%rain_number_tend_self_collection_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nraggtot'
      end if
      allocate(this%cldice_number_tend_autoconversion_to_snow_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nprcitot'
      end if
      allocate(this%cldliq_number_tend_sedimentation_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%ncsedten'
      end if
      allocate(this%cldice_number_tend_sedimentation_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nisedten'
      end if
      allocate(this%rain_number_tend_sedimentation_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nrsedten'
      end if
      allocate(this%snow_number_tend_sedimentation_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nssedten'
      end if
      allocate(this%graupel_number_tend_sedimentation_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%ngsedten'
      end if
      allocate(this%cldice_number_tend_melting_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nmelttot'
      end if
      allocate(this%snow_number_tend_melting_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nmeltstot'
      end if
      allocate(this%graupel_number_tend_melting_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nmeltgtot'
      end if
      allocate(this%cldliq_size_dist_slope_output_per_m(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
         errstring='Error allocating this%lamc_out'
      end if
      allocate(this%rain_size_dist_slope_output_per_m(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
         errstring='Error allocating this%lamr_out'
      end if
      allocate(this%cldliq_size_dist_shape_output(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
         errstring='Error allocating this%pgam_out'
      end if
      allocate(this%rain_size_dist_intercept_output_per_m4(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
         errstring='Error allocating this%n0r_out'
      end if

      ! Only allocate these variables if machine learning turned on

      if (trim(warm_rain) == 'tau' .or. trim(warm_rain) == 'emulated') then
         allocate(this%scale_qc(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%scale_qc'
         end if
         allocate(this%scale_nc(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%scale_nc'
         end if
         allocate(this%scale_qr(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%scale_qr'
         end if
         allocate(this%scale_nr(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%scale_nr'
         end if
         allocate(this%amk_c(psetcols,nlev,ncd), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%amk_c'
         end if
         allocate(this%ank_c(psetcols,nlev,ncd), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%ank_c'
         end if
         allocate(this%amk_r(psetcols,nlev,ncd), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%amk_r'
         end if
         allocate(this%ank_r(psetcols,nlev,ncd), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%ank_r'
         end if
         allocate(this%amk(psetcols,nlev,ncd), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%amk'
         end if
         allocate(this%ank(psetcols,nlev,ncd), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%ank'
         end if
         allocate(this%amk_out(psetcols,nlev,ncd), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%amk_out'
         end if
         allocate(this%ank_out(psetcols,nlev,ncd), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%ank_out'
         end if
         allocate(this%qc_out_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%qc_out_TAU'
         end if
         allocate(this%nc_out_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%nc_out_TAU'
         end if
         allocate(this%qr_out_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%qr_out_TAU'
         end if
         allocate(this%nr_out_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%nr_out_TAU'
         end if
         allocate(this%qc_in_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%qc_in_TAU'
         end if
         allocate(this%nc_in_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%nc_in_TAU'
         end if
         allocate(this%qr_in_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%qr_in_TAU'
         end if
         allocate(this%nr_in_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%nr_in_TAU'
         end if
         allocate(this%cldliq_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%qctend_TAU'
         end if
         allocate(this%cldliq_number_tend_autoconv_accretion_TAU_per_kg_per_s(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%nctend_TAU'
         end if
         allocate(this%rain_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%qrtend_TAU'
         end if
         allocate(this%rain_number_tend_autoconv_accretion_TAU_per_kg_per_s(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%nrtend_TAU'
         end if
         allocate(this%gmnnn_lmnnn_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%gmnnn_lmnnn_TAU'
         end if
         allocate(this%ML_fixer(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%ML_fixer'
         end if
         allocate(this%QC_fixer(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%QC_fixer'
         end if
         allocate(this%NC_fixer(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%NC_fixer'
         end if
         allocate(this%QR_fixer(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%QR_fixer'
         end if
         allocate(this%NR_fixer(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%NR_fixer'
         end if
      else if (warm_rain == 'sb2001') then
         ! Classic default (non-ML) microphysics
         allocate(this%cldliq_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%qctend_SB2001'
         end if
         allocate(this%cldliq_number_tend_autoconv_accretion_SB2001_per_kg_per_s(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%nctend_SB2001'
         end if
         allocate(this%rain_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%artend_SB2001'
         end if
         allocate(this%rain_number_tend_autoconv_accretion_SB2001_per_kg_per_s(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%nrtend_SB2001'
         end if
      end if

      ! Variables which are needed by all code (Machine Learning and non-ML)
      allocate(this%cldliq_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qctend_KK2000'
      end if
      allocate(this%cldliq_number_tend_autoconv_accretion_KK2000_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nctend_KK2000'
      end if
      allocate(this%rain_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%artend_KK2000'
      end if
      allocate(this%rain_number_tend_autoconv_accretion_KK2000_per_kg_per_s(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nrtend_KK2000'
      end if

   end subroutine proc_rates_allocate

   subroutine proc_rates_deallocate(this, warm_rain)
   !--------------------------------------------------------------
   ! Routine to deallocate the elements of the proc_rates DDT
   !--------------------------------------------------------------

      class(proc_rates_type) :: this
      character(len=16), intent(in) :: warm_rain            ! 'tau','emulated','sb2001' or 'kk2000'

      deallocate(this%snow_mass_tend_production_kg_per_kg_per_s)
      deallocate(this%snow_mass_tend_sublimation_kg_per_kg_per_s)
      deallocate(this%cldliq_mass_tend_sedimentation_evaporation_kg_per_kg_per_s)
      deallocate(this%cldice_mass_tend_sedimentation_sublimation_kg_per_kg_per_s)
      deallocate(this%vapor_mass_tend_residual_condensation_kg_per_kg_per_s)
      deallocate(this%cldice_mass_tend_net_deposition_kg_per_kg_per_s)
      deallocate(this%cldliq_mass_weighted_fall_speed_m_per_s)
      deallocate(this%cldice_mass_weighted_fall_speed_m_per_s)
      deallocate(this%rain_mass_weighted_fall_speed_m_per_s)
      deallocate(this%snow_mass_weighted_fall_speed_m_per_s)
      deallocate(this%graupel_mass_weighted_fall_speed_m_per_s)
      deallocate(this%graupel_mass_tend_sedimentation_kg_per_kg_per_s)
      deallocate(this%cldliq_mass_tend_sedimentation_kg_per_kg_per_s)
      deallocate(this%cldice_mass_tend_sedimentation_kg_per_kg_per_s)
      deallocate(this%rain_mass_tend_sedimentation_kg_per_kg_per_s)
      deallocate(this%snow_mass_tend_sedimentation_kg_per_kg_per_s)
      deallocate(this%cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s)
      deallocate(this%cldliq_mass_tend_autoconversion_to_rain_kg_per_kg_per_s)
      deallocate(this%cldliq_mass_tend_immersion_freezing_kg_per_kg_per_s)
      deallocate(this%cldliq_mass_tend_contact_freezing_kg_per_kg_per_s)
      deallocate(this%cldice_mass_tend_splintering_by_snow_kg_per_kg_per_s)
      deallocate(this%cldliq_mass_tend_riming_by_snow_kg_per_kg_per_s)
      deallocate(this%cldliq_mass_tend_bergeron_to_snow_kg_per_kg_per_s)
      deallocate(this%snow_mass_tend_deposition_kg_per_kg_per_s)
      deallocate(this%cldliq_mass_tend_bergeron_to_cldice_kg_per_kg_per_s)
      deallocate(this%cldice_mass_tend_melting_kg_per_kg_per_s)
      deallocate(this%snow_mass_tend_melting_kg_per_kg_per_s)
      deallocate(this%graupel_mass_tend_melting_kg_per_kg_per_s)
      deallocate(this%cldliq_mass_tend_homogeneous_freezing_kg_per_kg_per_s)
      deallocate(this%cldliq_mass_tend_residual_condensation_kg_per_kg_per_s)
      deallocate(this%cldice_mass_tend_autoconversion_to_snow_kg_per_kg_per_s)
      deallocate(this%cldice_mass_tend_accretion_by_snow_kg_per_kg_per_s)
      deallocate(this%cldice_mass_tend_residual_deposition_kg_per_kg_per_s)
      deallocate(this%rain_mass_tend_freezing_to_snow_kg_per_kg_per_s)
      deallocate(this%cldliq_mass_tend_deposition_nucleation_kg_per_kg_per_s)
      deallocate(this%rain_mass_tend_freezing_to_cldice_kg_per_kg_per_s)
      deallocate(this%rain_mass_tend_accretion_by_snow_kg_per_kg_per_s)
      deallocate(this%air_heating_rate_snow_melting_W_per_kg)
      deallocate(this%air_heating_rate_rain_freezing_W_per_kg)
      deallocate(this%cldice_mass_tend_nucleation_from_vapor_kg_per_kg_per_s)
      deallocate(this%rain_mass_tend_accretion_by_graupel_kg_per_kg_per_s)
      deallocate(this%cldliq_mass_tend_riming_by_graupel_kg_per_kg_per_s)
      deallocate(this%graupel_mass_tend_from_snow_riming_kg_per_kg_per_s)
      deallocate(this%graupel_mass_tend_rain_snow_collection_kg_per_kg_per_s)
      deallocate(this%graupel_mass_tend_deposition_kg_per_kg_per_s)
      deallocate(this%cldice_mass_tend_splintering_by_graupel_kg_per_kg_per_s)
      deallocate(this%cldice_mass_tend_splintering_rain_graupel_kg_per_kg_per_s)
      deallocate(this%snow_mass_tend_collection_by_rain_kg_per_kg_per_s)
      deallocate(this%rain_number_tend_accretion_by_graupel_per_kg_per_s)
      deallocate(this%graupel_number_tend_from_snow_riming_per_kg_per_s)
      deallocate(this%graupel_number_tend_rain_snow_collection_per_kg_per_s)
      deallocate(this%cldice_number_tend_splintering_by_graupel_per_kg_per_s)
      deallocate(this%cldice_number_tend_splintering_rain_graupel_per_kg_per_s)
      deallocate(this%cldliq_number_tend_riming_by_graupel_per_kg_per_s)
      deallocate(this%cldliq_number_tend_immersion_freezing_per_kg_per_s)
      deallocate(this%cldliq_number_tend_contact_freezing_per_kg_per_s)
      deallocate(this%cldice_number_tend_nucleation_from_vapor_per_kg_per_s)
      deallocate(this%cldliq_number_tend_deposition_nucleation_per_kg_per_s)
      deallocate(this%cldliq_number_tend_homogeneous_freezing_per_kg_per_s)
      deallocate(this%rain_number_tend_freezing_to_snow_per_kg_per_s)
      deallocate(this%rain_number_tend_freezing_to_cldice_per_kg_per_s)
      deallocate(this%cldice_number_tend_splintering_by_snow_per_kg_per_s)
      deallocate(this%cldliq_number_tend_accretion_by_rain_per_kg_per_s)
      deallocate(this%cldliq_number_tend_riming_by_snow_per_kg_per_s)
      deallocate(this%cldice_number_tend_accretion_by_snow_per_kg_per_s)
      deallocate(this%rain_number_tend_accretion_by_snow_per_kg_per_s)
      deallocate(this%rain_number_tend_autoconversion_from_cldliq_per_kg_per_s)
      deallocate(this%rain_number_tend_self_collection_per_kg_per_s)
      deallocate(this%cldice_number_tend_autoconversion_to_snow_per_kg_per_s)
      deallocate(this%cldliq_number_tend_sedimentation_per_kg_per_s)
      deallocate(this%cldice_number_tend_sedimentation_per_kg_per_s)
      deallocate(this%rain_number_tend_sedimentation_per_kg_per_s)
      deallocate(this%snow_number_tend_sedimentation_per_kg_per_s)
      deallocate(this%graupel_number_tend_sedimentation_per_kg_per_s)
      deallocate(this%cldice_number_tend_melting_per_kg_per_s)
      deallocate(this%snow_number_tend_melting_per_kg_per_s)
      deallocate(this%graupel_number_tend_melting_per_kg_per_s)

      deallocate(this%cldliq_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s)
      deallocate(this%cldliq_number_tend_autoconv_accretion_KK2000_per_kg_per_s)
      deallocate(this%rain_mass_tend_autoconv_accretion_KK2000_kg_per_kg_per_s)
      deallocate(this%rain_number_tend_autoconv_accretion_KK2000_per_kg_per_s)

      deallocate(this%cldliq_size_dist_slope_output_per_m)
      deallocate(this%rain_size_dist_slope_output_per_m)
      deallocate(this%cldliq_size_dist_shape_output)
      deallocate(this%rain_size_dist_intercept_output_per_m4)

      if (trim(warm_rain) == 'tau' .or. trim(warm_rain) == 'emulated') then
         deallocate(this%scale_qc)
         deallocate(this%scale_nc)
         deallocate(this%scale_qr)
         deallocate(this%scale_nr)
         deallocate(this%amk_c)
         deallocate(this%ank_c)
         deallocate(this%amk_r)
         deallocate(this%ank_r)
         deallocate(this%amk)
         deallocate(this%ank)
         deallocate(this%amk_out)
         deallocate(this%ank_out)
         deallocate(this%qc_out_TAU)
         deallocate(this%nc_out_TAU)
         deallocate(this%qr_out_TAU)
         deallocate(this%nr_out_TAU)
         deallocate(this%qc_in_TAU)
         deallocate(this%nc_in_TAU)
         deallocate(this%qr_in_TAU)
         deallocate(this%nr_in_TAU)
         deallocate(this%cldliq_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s)
         deallocate(this%cldliq_number_tend_autoconv_accretion_TAU_per_kg_per_s)
         deallocate(this%rain_mass_tend_autoconv_accretion_TAU_kg_per_kg_per_s)
         deallocate(this%rain_number_tend_autoconv_accretion_TAU_per_kg_per_s)
         deallocate(this%gmnnn_lmnnn_TAU)
         deallocate(this%ML_fixer)
         deallocate(this%QC_fixer)
         deallocate(this%NC_fixer)
         deallocate(this%QR_fixer)
         deallocate(this%NR_fixer)
      else if (trim(warm_rain) == 'sb2001') then
         deallocate(this%cldliq_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s)
         deallocate(this%cldliq_number_tend_autoconv_accretion_SB2001_per_kg_per_s)
         deallocate(this%rain_mass_tend_autoconv_accretion_SB2001_kg_per_kg_per_s)
         deallocate(this%rain_number_tend_autoconv_accretion_SB2001_per_kg_per_s)
      end if

   end subroutine proc_rates_deallocate

end module micro_pumas_diags
