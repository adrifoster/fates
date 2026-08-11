program FatesCanopyLevelPhoto
  !
  ! DESCRIPTION:
  ! Canopy-level photosynthesis sensitivity sweep: five independent sweeps (PAR,
  ! CO2, VPD, leaf temperature, and soil water content), each varying one
  ! driver variable while holding everything else at a fixed default,
  ! evaluated for a single PFT 
  !
  ! The Kumarathunge et al. (2019) temperature-acclimation model needs two running-mean
  ! reference temperatures (t_growth, t_home) in addition to the instantaneous leaf
  ! temperature being swept. Both are held fixed at the same default leaf temperature
  ! for every sweep and every point.
  !
  ! The soil water content sweep does not vary btran directly. It sweeps a soil
  ! water content fraction in [0, 1], maps it onto a soil matric potential, since 
  ! these drivers have no soil texture from which a real retention curve could be built), 
  ! and derives btran from that via BtranFromSMP, which is production's own 
  ! EDBtranMod.F90::btran_ed formula specialized to the single unfrozen layer at 
  ! root fraction 1 these drivers assume

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesConstantsMod,           only : fates_unset_r8
  use FatesConstantsMod,           only : wm2_to_umolm2s
  use FatesConstantsMod,           only : t_water_freeze_k_1atm
  use EDParamsMod,                 only : nlevleaf
  use EDParamsMod,                 only : GetNVegLayers
  use EDPftvarcon,                 only : EDPftvarcon_inst
  use PRTParametersMod,            only : prt_params
  use PRTInitParamsFatesMod,       only : PRTDerivedParams
  use FatesParameterDerivedMod,    only : param_derived
  use FatesUnitTestParamReaderMod, only : ReadParameters
  use FatesArgumentUtils,          only : command_line_arg
  use FatesFactoryMod,             only : InitializeGlobals
  use FatesGlobals,                only : FatesGlobalsInit
  use FatesInterfaceTypesMod,      only : numpft
  use FatesTwoStreamUtilsMod,      only : TransferRadParams
  use LeafBiophysicsMod,           only : lb_params
  use LeafBiophysicsMod,           only : FvCB1980, medlyn_model, net_assim_model
  use LeafBiophysicsMod,           only : photosynth_acclim_model_kumarathunge_etal_2019
  use LeafBiophysicsMod,           only : QSat
  use LeafBiophysicsMod,           only : GetCanopyGasParameters
  use FatesTestEnvironmentMod,     only : environment_type
  use FatesTestEnvironmentMod,     only : BtranFromSMP, SoilMatricPotential
  use FatesTestSiteMod,            only : constant_vpd
  use FatesTestLightEnvMod,        only : light_env_type
  use FatesTestLightEnvMod,        only : direct_frac, diffuse_frac
  use FatesTestLeafPhotoMod,       only : leaf_capacity_type
  use FatesTestLeafPhotoMod,       only : LeafLayerCapacity
  use FatesTestLeafPhotoMod,       only : LeafLayerSunShade
  use FatesTestLeafPhotoMod,       only : LeafNitrogenContent
  use FatesTestLeafPhotoMod,       only : LeafLayerNitrogenScaling
  use FatesUnitTestIOMod,          only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod,          only : WriteVar, RegisterVarAtts, EndNCDef
  use FatesUnitTestIOMod,          only : RegisterFillValue
  use FatesUnitTestIOMod,          only : type_double, type_int

  implicit none

  ! LOCALS:
  character(len=:), allocatable :: param_file ! input parameter file
  type(environment_type)        :: env        ! prescribed atmospheric boundary conditions
  type(light_env_type)          :: light_env  ! prescribed light environment for the current LAI
  real(r8)                      :: lnc_top    ! leaf N content at the canopy top [gN/m2 leaf]
  real(r8)                      :: vcmax25top ! top-of-canopy carboxylation rate at 25degC [umol/m2/s]
  real(r8)                      :: jmax25top  ! top-of-canopy electron transport rate at 25degC [umol/m2/s]
  real(r8)                      :: kp25top    ! top-of-canopy initial slope of C4 CO2 response at 25degC [umol/m2/s]
  real(r8)                      :: smpsc      ! soil matric potential at full stomatal closure [mm, negative]
  real(r8)                      :: smpso      ! soil matric potential at full stomatal opening [mm, negative]

  ! per-leaf-layer light profile at each prescribed LAI, (nlevleaf, n_lai)
  real(r8), allocatable :: parsun_z_out(:,:)  ! absorbed PAR, sunlit [W/m2 ground]
  real(r8), allocatable :: parsha_z_out(:,:)  ! absorbed PAR, shaded [W/m2 ground]
  real(r8), allocatable :: laisun_z_out(:,:)  ! sunlit leaf area index [m2/m2]
  real(r8), allocatable :: laisha_z_out(:,:)  ! shaded leaf area index [m2/m2]
  real(r8), allocatable :: nscaler_z_out(:,:) ! per-layer nitrogen-scaling factor [0-1]
  real(r8), allocatable :: anet_z_out(:,:)    ! per-layer area-weighted net photosynthesis [umolC/m2 leaf/s]
  integer,  allocatable :: nv_out(:)          ! number of occupied leaf layers at each prescribed LAI
  integer,  allocatable :: layer_index(:)     ! leaf-layer index coordinate [1..nlevleaf]

  ! canopy-integrated output at the reference condition, per prescribed LAI
  real(r8), allocatable :: canopy_anet(:)   ! canopy net photosynthesis [umolC/m2 ground/s]
  real(r8), allocatable :: canopy_agross(:) ! canopy gross photosynthesis [umolC/m2 ground/s]

  ! per-layer working arrays, sized to the current LAI's nv
  real(r8), allocatable :: nscaler_z(:) ! per-leaf-layer nitrogen-scaling factor [0-1]

  ! the shared default vapor-pressure state at the standard reference conditions
  real(r8) :: default_veg_esat   ! saturation vapor pressure at default veg tempK [Pa]
  real(r8) :: default_can_vpress ! canopy air vapor pressure at the default VPD [Pa]
  real(r8) :: qs_dummy           ! saturation specific humidity output from QSat (unused here)

  ! the per-layer photosynthesis working variables are locals of CanopyNetAssim
  ! rather than program-scope names it reaches by host association - a contained
  ! subroutine silently sharing the program's scratch space (its loop index in
  ! particular) is a hazard the moment a second loop exists

  integer :: ilai ! prescribed-LAI looping index
  integer :: iv   ! leaf-layer looping index (this program's own layer_index loop; CanopyNetAssim declares its own)
  integer :: i    ! sweep-point looping index

  ! swept-variable value arrays - same swept variables and ranges as
  ! test_LeafLevelPhoto.F90, so the leaf and canopy panels are the same
  ! experiment at two scales
  real(r8), allocatable :: par_vals(:)      ! swept incident PPFD at the top of the canopy [umol/m2/s]
  real(r8), allocatable :: co2_vals(:)      ! swept CO2 partial pressure values [Pa]
  real(r8), allocatable :: vpd_vals(:)      ! swept leaf-to-air VPD values [Pa]
  real(r8), allocatable :: temp_vals(:)     ! swept leaf temperature values [K]
  real(r8), allocatable :: soilfrac_vals(:) ! swept soil water content, fraction of saturation [0-1]

  ! derived vapor-pressure/btran diagnostics, dimensioned like their own sweep
  ! (only temperature and VPD change vapor pressure, only soil water changes btran)
  real(r8), allocatable :: veg_esat_bytemp(:)   ! saturation vapor pressure at each swept leaf temperature [Pa]
  real(r8), allocatable :: can_vpress_bytemp(:) ! canopy air vapor pressure at each swept leaf temperature, fixed default VPD [Pa]
  real(r8), allocatable :: veg_esat_byvpd(:)    ! saturation vapor pressure at the fixed default leaf temperature [Pa]
  real(r8), allocatable :: can_vpress_byvpd(:)  ! canopy air vapor pressure at each swept VPD (= veg_esat - vpd) [Pa]
  real(r8), allocatable :: btran_bysoilfrac(:)  ! btran derived from each swept soil water content fraction [0-1]

  ! sweep output, each dimensioned (n_sweep_points, n_lai)
  real(r8), allocatable :: canopy_anet_bypar(:,:), canopy_agross_bypar(:,:)
  real(r8), allocatable :: canopy_anet_byco2(:,:), canopy_agross_byco2(:,:)
  real(r8), allocatable :: canopy_anet_byvpd(:,:), canopy_agross_byvpd(:,:)
  real(r8), allocatable :: canopy_anet_bytemp(:,:), canopy_agross_bytemp(:,:)
  real(r8), allocatable :: canopy_anet_bysoilfrac(:,:), canopy_agross_bysoilfrac(:,:)

  integer :: n_par, n_co2, n_vpd, n_temp, n_soilfrac ! sweep array sizes

  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'canopy_level_photo_out.nc' ! output file
  integer,  parameter :: target_pft = 6 ! PFT index to evaluate (1-based)

  ! prescribed canopy
  integer,  parameter :: n_lai = 3
  real(r8), parameter :: lai_vals(n_lai) = [1.0_r8, 3.0_r8, 7.0_r8] ! prescribed leaf area index [m2 leaf/m2 ground]
  real(r8), parameter :: canopy_sai    = 0.0_r8  ! prescribed stem area index [m2 stem/m2 ground]
  real(r8), parameter :: canopy_height = 20.0_r8 ! prescribed canopy height [m]

  ! illumination geometry
  real(r8), parameter :: diagnostic_coszen = 1.0_r8 ! cosine of solar zenith angle (sun directly overhead)

  ! incident PAR
  real(r8), parameter :: reference_par = 1500.0_r8/wm2_to_umolm2s ! [W/m2]

  ! Rogers' standard reference conditions - identical to
  real(r8), parameter :: default_veg_tempk = 25.0_r8 + t_water_freeze_k_1atm ! [K]
  ! the reference VPD is FatesTestSiteMod's constant_vpd rather than a local
  ! constant, so that this sweep and single_cohort's constant-VPD boundary
  ! condition are the same number by construction and cannot drift apart
  real(r8), parameter :: default_ppfd      = 1500.0_r8 ! [umol/m2/s] incident PPFD at the top of the canopy (= reference_par, in the units the PAR sweep is expressed in)

  ! sweep ranges
  real(r8), parameter :: min_temp = 5.0_r8,    max_temp = 40.0_r8,    temp_inc = 0.5_r8   ! [degC]
  real(r8), parameter :: min_par  = 0.0_r8,    max_par  = 1600.0_r8,  par_inc  = 5.0_r8   ! [umol/m2/s]
  real(r8), parameter :: min_vpd  = 500.0_r8,  max_vpd  = 2500.0_r8,  vpd_inc  = 20.0_r8  ! [Pa] (0.5-2.5 kPa)
  real(r8), parameter :: min_co2  = 250.0_r8,  max_co2  = 1000.0_r8,  co2_inc  = 5.0_r8   ! [umol/mol]
  real(r8), parameter :: soilfrac_inc = 0.02_r8 ! [0-1]

  ! read in parameter file name from command line
  param_file = command_line_arg(1)

  ! read in parameter file
  call ReadParameters(param_file)

  ! step_size is unused by anything this test exercises (no time-stepping
  ! occurs), but InitializeGlobals requires a value
  call InitializeGlobals(86400.0_r8)

  ! initialize FATES logging (the two-stream module's debug/error paths write
  ! to this unit and it is otherwise never set in a standalone driver) and the
  ! two-stream radiation parameters (leaf/stem optical properties, per pft)
  numpft = size(prt_params%wood_density, dim=1)
  call FatesGlobalsInit(6, .false.)
  call TransferRadParams()

  ! derive the organ_id -> parameter-file-index reverse lookup map
  ! (prt_params%organ_param_id) - normally done by the host model's own
  ! interface setup (FatesInterfaceMod.F90), which this standalone driver
  ! bypasses entirely
  call PRTDerivedParams()

  ! host-model-namelist-controlled leaf biophysics switches - matches
  ! test_LeafLevelPhoto.F90/test_SingleCohort.F90's choices for consistency
  ! across this test suite. Not consumed at this light-profile-only stage, but
  ! set here so the setup sequence is complete and correct as photosynthesis
  ! is added on top
  lb_params%electron_transport_model = FvCB1980
  lb_params%stomatal_model           = medlyn_model
  lb_params%stomatal_assim_model     = net_assim_model
  lb_params%photo_tempsens_model     = photosynth_acclim_model_kumarathunge_etal_2019

  ! leaf N content and reference (25C, canopy-top) photosynthetic capacity for
  ! target_pft
  lnc_top        = LeafNitrogenContent(target_pft)
  vcmax25top_pft = EDPftvarcon_inst%vcmax25top(target_pft,1)
  jmax25top_pft  = param_derived%jmax25top(target_pft,1)
  kp25top_pft    = param_derived%kp25top(target_pft,1)

  ! target_pft's soil matric potential thresholds - used only by the soil
  ! water content sweep
  smpsc_pft = EDPftvarcon_inst%smpsc(target_pft)
  smpso_pft = EDPftvarcon_inst%smpso(target_pft)

  ! prescribed atmospheric boundary conditions (canopy pressure, background
  ! CO2/O2 partial pressure, leaf boundary-layer conductance, dayl_factor,
  ! btran - see FatesTestEnvironmentMod's shared reference-atmosphere defaults)
  call env%Init()

  ! ---------------------------------------------------------------------------------------
  ! resolve the within-canopy light profile at each prescribed LAI
  ! ---------------------------------------------------------------------------------------

  allocate(parsun_z_out(nlevleaf, n_lai), parsha_z_out(nlevleaf, n_lai))
  allocate(laisun_z_out(nlevleaf, n_lai), laisha_z_out(nlevleaf, n_lai))
  allocate(nscaler_z_out(nlevleaf, n_lai), anet_z_out(nlevleaf, n_lai))
  allocate(nv_out(n_lai), layer_index(nlevleaf))
  allocate(canopy_anet(n_lai), canopy_agross(n_lai))

  ! layers above a given LAI's actual nv are never written, and keep this fill
  ! value (registered as each variable's _FillValue below)
  parsun_z_out(:,:)  = fates_unset_r8
  parsha_z_out(:,:)  = fates_unset_r8
  laisun_z_out(:,:)  = fates_unset_r8
  laisha_z_out(:,:)  = fates_unset_r8
  nscaler_z_out(:,:) = fates_unset_r8
  anet_z_out(:,:)    = fates_unset_r8

  do iv = 1, nlevleaf
    layer_index(iv) = iv
  end do

  ! the standard reference condition's vapor-pressure state - leaf
  ! temperature and VPD are both fixed at this stage, so this is computed
  ! once for the whole run
  call QSat(default_veg_tempk, env%can_press, qs_dummy, default_veg_esat)
  default_can_vpress = default_veg_esat - constant_vpd

  ! ---------------------------------------------------------------------------------------
  ! build the swept-value arrays and the diagnostics derived from them
  ! ---------------------------------------------------------------------------------------

  n_par      = int((max_par - min_par)/par_inc) + 1
  n_co2      = int((max_co2 - min_co2)/co2_inc) + 1
  n_vpd      = int((max_vpd - min_vpd)/vpd_inc) + 1
  n_temp     = int((max_temp - min_temp)/temp_inc) + 1
  n_soilfrac = int((1.0_r8 - 0.0_r8)/soilfrac_inc) + 1

  allocate(par_vals(n_par), co2_vals(n_co2), vpd_vals(n_vpd))
  allocate(temp_vals(n_temp), soilfrac_vals(n_soilfrac))
  allocate(veg_esat_bytemp(n_temp), can_vpress_bytemp(n_temp))
  allocate(veg_esat_byvpd(n_vpd), can_vpress_byvpd(n_vpd))
  allocate(btran_bysoilfrac(n_soilfrac))

  allocate(canopy_anet_bypar(n_par, n_lai), canopy_agross_bypar(n_par, n_lai))
  allocate(canopy_anet_byco2(n_co2, n_lai), canopy_agross_byco2(n_co2, n_lai))
  allocate(canopy_anet_byvpd(n_vpd, n_lai), canopy_agross_byvpd(n_vpd, n_lai))
  allocate(canopy_anet_bytemp(n_temp, n_lai), canopy_agross_bytemp(n_temp, n_lai))
  allocate(canopy_anet_bysoilfrac(n_soilfrac, n_lai), canopy_agross_bysoilfrac(n_soilfrac, n_lai))

  do i = 1, n_par
    par_vals(i) = min_par + par_inc*real(i-1, r8)
  end do
  do i = 1, n_co2
    ! ppm -> Pa, at the default (Init()-prescribed) canopy pressure
    co2_vals(i) = ((min_co2 + co2_inc*real(i-1, r8))/1.0e6_r8) * env%can_press
  end do
  do i = 1, n_vpd
    ! leaf temperature is fixed across this sweep, so veg_esat stays at its
    ! default and the swept VPD moves canopy air vapor pressure alone
    vpd_vals(i) = min_vpd + vpd_inc*real(i-1, r8)
    veg_esat_byvpd(i) = default_veg_esat
    can_vpress_byvpd(i) = default_veg_esat - vpd_vals(i)
  end do
  do i = 1, n_temp
    ! VPD is held at the default as leaf temperature varies, so veg_esat is
    ! recomputed per point and canopy air vapor pressure follows it
    temp_vals(i) = (min_temp + temp_inc*real(i-1, r8)) + t_water_freeze_k_1atm
    call QSat(temp_vals(i), env%can_press, qs_dummy, veg_esat_bytemp(i))
    can_vpress_bytemp(i) = veg_esat_bytemp(i) - constant_vpd
  end do
  do i = 1, n_soilfrac
    ! soil water content -> matric potential -> btran, via the same shared
    ! routines test_LeafLevelPhoto.F90's soil water sweep uses
    soilfrac_vals(i) = 0.0_r8 + soilfrac_inc*real(i-1, r8)
    btran_bysoilfrac(i) = BtranFromSMP(                                        &
      SoilMatricPotential(soilfrac_vals(i), smpsc_pft), smpsc_pft, smpso_pft)
  end do

  do ilai = 1, n_lai

    ! build the prescribed canopy's two-stream light environment - no cohort:
    ! Init takes canopy structure as scalars precisely so a prescribed-LAI
    ! driver like this one can use it
    call light_env%Init(lai_vals(ilai), canopy_sai, canopy_height, target_pft)

    nv_out(ilai) = GetNVegLayers(lai_vals(ilai) + canopy_sai)

    ! this canopy's per-layer nitrogen scaling - the vertical decay of
    ! photosynthetic capacity with cumulative leaf area above each layer,
    ! from the same shared routine test_SingleCohort.F90 drives off a real
    ! cohort's allometry
    if (allocated(nscaler_z)) deallocate(nscaler_z)
    allocate(nscaler_z(nv_out(ilai)))
    call LeafLayerNitrogenScaling(lai_vals(ilai), canopy_sai, canopy_height,   &
      nv_out(ilai), target_pft, vcmax25top_pft, nscaler_z)

    ! the reference case: this driver's prescribed clear-sky beam fraction.
    ! This is also the only call that stores the per-leaf-layer profile - the
    ! sweeps below would otherwise write nlevleaf x (sum of sweep lengths) x
    ! n_lai layer profiles, which is a great deal of output for a diagnostic
    ! whose point is the canopy-integrated response
    call CanopyNetAssim(nv_out(ilai), nscaler_z, reference_par, direct_frac,   &
      default_veg_tempk, default_veg_tempk, default_veg_tempk,                 &
      default_veg_esat, default_can_vpress, env%can_co2_ppress, env%btran,     &
      vcmax25top_pft, canopy_anet(ilai), canopy_agross(ilai),                  &
      store_profile=.true., ilai_store=ilai)

    ! every sweep point, at this LAI. t_growth/t_home are held at the default
    ! leaf temperature for every sweep and every point - including the
    ! temperature sweep

    ! PAR sweep - incident PPFD converted to the W/m2 the light environment works in
    do i = 1, n_par
      call CanopyNetAssim(nv_out(ilai), nscaler_z, par_vals(i)/wm2_to_umolm2s, &
        direct_frac, default_veg_tempk, default_veg_tempk, default_veg_tempk,  &
        default_veg_esat, default_can_vpress, env%can_co2_ppress, env%btran,   &
        vcmax25top_pft, canopy_anet_bypar(i,ilai), canopy_agross_bypar(i,ilai))
    end do

    ! CO2 sweep
    do i = 1, n_co2
      call CanopyNetAssim(nv_out(ilai), nscaler_z, reference_par, direct_frac, &
        default_veg_tempk, default_veg_tempk, default_veg_tempk,               &
        default_veg_esat, default_can_vpress, co2_vals(i), env%btran,          &
        vcmax25top_pft, canopy_anet_byco2(i,ilai), canopy_agross_byco2(i,ilai))
    end do

    ! VPD sweep
    do i = 1, n_vpd
      call CanopyNetAssim(nv_out(ilai), nscaler_z, reference_par, direct_frac, &
        default_veg_tempk, default_veg_tempk, default_veg_tempk,               &
        veg_esat_byvpd(i), can_vpress_byvpd(i), env%can_co2_ppress,            &
        env%btran, vcmax25top_pft, canopy_anet_byvpd(i,ilai),                  &
        canopy_agross_byvpd(i,ilai))
    end do

    ! leaf temperature sweep
    do i = 1, n_temp
      call CanopyNetAssim(nv_out(ilai), nscaler_z, reference_par, direct_frac, &
        temp_vals(i), default_veg_tempk, default_veg_tempk,                    &
        veg_esat_bytemp(i), can_vpress_bytemp(i), env%can_co2_ppress,          &
        env%btran, vcmax25top_pft, canopy_anet_bytemp(i,ilai),                 &
        canopy_agross_bytemp(i,ilai))
    end do

    ! soil water content sweep
    do i = 1, n_soilfrac
      call CanopyNetAssim(nv_out(ilai), nscaler_z, reference_par, direct_frac, &
        default_veg_tempk, default_veg_tempk, default_veg_tempk,               &
        default_veg_esat, default_can_vpress, env%can_co2_ppress,              &
        btran_bysoilfrac(i), vcmax25top_pft,                                   &
        canopy_anet_bysoilfrac(i,ilai), canopy_agross_bysoilfrac(i,ilai))
    end do

    call light_env%Free()

  end do

  call WriteOutput()

contains


  ! ==========================================================================

  subroutine CanopyNetAssim(nv, nscaler_z, par_toc, beam_frac, veg_tempk,      &
    t_growth, t_home, veg_esat, can_vpress, can_co2_ppress, btran, vcmax25top, &
    canopy_anet_out, canopy_agross_out, store_profile, ilai_store)
    !
    ! DESCRIPTION:
    ! Integrates leaf photosynthesis down a canopy
    ! (light_env, holding this LAI's structure) to a canopy net assimilation
    ! per unit ground area - see the program header for the method.
    !
    ! Writes nothing except canopy_anet_out unless store_profile is set, in
    ! which case the per-layer diagnostics for column ilai_store of the
    ! output arrays are filled as well.

    ! ARGUMENTS:
    integer,  intent(in)  :: nv              ! number of occupied leaf layers
    real(r8), intent(in)  :: nscaler_z(:)    ! per-leaf-layer nitrogen-scaling factor [0-1]
    real(r8), intent(in)  :: par_toc         ! incident PAR at the top of the canopy [W/m2]
    real(r8), intent(in)  :: beam_frac       ! fraction of par_toc arriving as direct beam [0-1]
    real(r8), intent(in)  :: veg_tempk       ! instantaneous leaf temperature [K]
    real(r8), intent(in)  :: t_growth        ! 10-day running-mean growth temperature [K] - held at the reference temperature by every sweep, see the call site
    real(r8), intent(in)  :: t_home          ! long-term running-mean home temperature [K] - held at the reference temperature by every sweep, see the call site
    real(r8), intent(in)  :: veg_esat        ! saturation vapor pressure at veg_tempk [Pa]
    real(r8), intent(in)  :: can_vpress      ! canopy air vapor pressure [Pa]
    real(r8), intent(in)  :: can_co2_ppress  ! CO2 partial pressure at the leaf surface [Pa]
    real(r8), intent(in)  :: btran           ! soil moisture stress factor [0-1]
    real(r8), intent(in)  :: vcmax25top      ! reference (25C, canopy-top) maximum carboxylation rate [umol/m2/s]
    real(r8), intent(out) :: canopy_anet_out   ! canopy net photosynthesis [umolC/m2 ground/s]
    real(r8), intent(out) :: canopy_agross_out ! canopy gross photosynthesis [umolC/m2 ground/s]
    logical,  intent(in), optional :: store_profile ! also fill the per-layer output arrays (default .false.)
    integer,  intent(in), optional :: ilai_store    ! output-array column to fill when store_profile is set

    ! LOCALS:
    logical :: do_store ! resolved store_profile
    integer :: iv       ! leaf-layer looping index
    type(leaf_capacity_type) :: cap ! this layer's capacity/dark respiration at the swept conditions
    real(r8) :: mm_kco2, mm_ko2, co2_cpoint ! Michaelis-Menten constants for CO2/O2, CO2 compensation point at veg_tempk [Pa]
    real(r8) :: agross_layer ! area-weighted gross photosynthesis for this layer [umolC/m2 leaf/s]
    real(r8) :: anet_layer   ! area-weighted net photosynthesis for this layer [umolC/m2 leaf/s]
    real(r8) :: lai_layer    ! this layer's total leaf area index [m2 leaf/m2 ground]

    do_store = .false.
    if (present(store_profile)) do_store = store_profile

    ! attenuate the incident PAR through this canopy at the prescribed
    ! overhead-sun geometry. AttenuateCanopy (rather than Profile) is called
    ! directly because Profile derives its incident PAR from a light fraction
    ! and a day/hour solar cycle, which this static, prescribed-Q protocol
    ! has no equivalent of
    call light_env%AttenuateCanopy(beam_frac*par_toc,                          &
      (1.0_r8 - beam_frac)*par_toc, diagnostic_coszen,                         &
      light_env%parsun_z, light_env%parsha_z, light_env%laisun_z,              &
      light_env%laisha_z)

    ! the Michaelis-Menten constants/CO2 compensation point depend only on
    ! can_press/can_o2_ppress/veg_tempk, none of which vary down the canopy at a
    ! single swept point, so GetCanopyGasParameters is hoisted out of the layer
    ! loop rather than recomputed per photosynthesis call
    call GetCanopyGasParameters(env%can_press, env%can_o2_ppress, veg_tempk,   &
      mm_kco2, mm_ko2, co2_cpoint)

    canopy_anet_out = 0.0_r8
    canopy_agross_out = 0.0_r8
    do iv = 1, nv

      ! this layer's capacity at the swept conditions - recomputed per point
      ! rather than held, since several of the sweeps move something it
      ! depends on (leaf temperature, btran)
      call LeafLayerCapacity(target_pft, veg_tempk, t_growth, t_home,          &
        nscaler_z(iv), env%dayl_factor, btran, vcmax25top, jmax25top_pft,      &
        kp25top_pft, lnc_top, cap)

      ! sunlit and shaded leaves in this layer, area-weighted into one rate
      call LeafLayerSunShade(target_pft, cap, light_env%parsun_z(iv),          &
        light_env%parsha_z(iv), light_env%laisun_z(iv),                        &
        light_env%laisha_z(iv), veg_tempk, env%can_press, can_co2_ppress,      &
        env%can_o2_ppress, veg_esat, can_vpress, env%gb, mm_kco2, mm_ko2,      &
        co2_cpoint, agross_layer, anet_layer, lai_layer)

      ! scaled by this layer's leaf area index: [umolC/m2 leaf/s] *
      ! [m2 leaf/m2 ground] -> [umolC/m2 ground/s], accumulated down the canopy.
      ! Gross is accumulated alongside net at no extra cost (LeafLayerSunShade
      ! returns both) - net is the canopy flux Rogers' panels report, gross
      ! separates the leaf dark respiration internal to it
      canopy_anet_out = canopy_anet_out + anet_layer * lai_layer
      canopy_agross_out = canopy_agross_out + agross_layer * lai_layer

      if (do_store) then
        parsun_z_out(iv, ilai_store)  = light_env%parsun_z(iv)
        parsha_z_out(iv, ilai_store)  = light_env%parsha_z(iv)
        laisun_z_out(iv, ilai_store)  = light_env%laisun_z(iv)
        laisha_z_out(iv, ilai_store)  = light_env%laisha_z(iv)
        nscaler_z_out(iv, ilai_store) = nscaler_z(iv)
        anet_z_out(iv, ilai_store)    = anet_layer
      end if

    end do

  end subroutine CanopyNetAssim

  ! ==========================================================================

  subroutine WriteOutput()
    !
    ! DESCRIPTION:
    ! Writes the prescribed LAI values, the per-leaf-layer light profile
    ! resolved at each of them (at the reference condition only - the sweeps
    ! would otherwise write nlevleaf x every sweep point x n_lai layer
    ! profiles), and every sweep's swept values, derived vapor-pressure/btran
    ! diagnostics, and canopy-integrated gross/net photosynthesis per LAI.

    ! LOCALS:
    integer           :: ncid          ! netcdf file id
    character(len=20) :: dim_names(7)  ! dimension names
    integer           :: dimIDs(7)     ! dimension IDs
    integer           :: laiID, layerID, nvID, canopyanetID, canopyagrossID
    integer           :: parsunID, parshaID, laisunID, laishaID
    integer           :: nscalerID, anetzID
    integer           :: parID, co2ID, vpdID, tempID, soilfracID
    integer           :: vegesatbytempID, canvpressbytempID
    integer           :: vegesatbyvpdID, canvpressbyvpdID, btranbysoilfracID
    integer           :: anetbyparID, agrossbyparID
    integer           :: anetbyco2ID, agrossbyco2ID
    integer           :: anetbyvpdID, agrossbyvpdID
    integer           :: anetbytempID, agrossbytempID
    integer           :: anetbysoilfracID, agrossbysoilfracID

    dim_names = [character(len=20) :: 'layer', 'lai', 'par', 'co2', 'vpd',     &
      'temp', 'soilfrac']

    call OpenNCFile(trim(out_file), ncid, 'readwrite')
    call RegisterNCDims(ncid, dim_names, (/nlevleaf, n_lai, n_par, n_co2,      &
      n_vpd, n_temp, n_soilfrac/), 7, dimIDs)

    call RegisterVarAtts(ncid, 'layer', dimIDs(1:1), type_int, '-',                     &
      'leaf layer index, 1 = top of canopy', layerID)
    call RegisterVarAtts(ncid, 'lai', dimIDs(2:2), type_double, 'm2 m-2',               &
      'prescribed canopy leaf area index', laiID)
    call RegisterVarAtts(ncid, 'nv', dimIDs(2:2), type_int, '-',                        &
      'number of occupied leaf layers at each prescribed LAI', nvID)
    call RegisterVarAtts(ncid, 'canopy_anet', dimIDs(2:2), type_double,                 &
      'umolC m-2 s-1',                                                                  &
      'canopy net photosynthesis per unit ground area, at the reference condition',     &
      canopyanetID)
    call RegisterVarAtts(ncid, 'canopy_agross', dimIDs(2:2), type_double,               &
      'umolC m-2 s-1',                                                                  &
      'canopy gross photosynthesis per unit ground area, at the reference condition',   &
      canopyagrossID)

    call RegisterVarAtts(ncid, 'parsun_z', (/dimIDs(1), dimIDs(2)/), type_double,       &
      'W m-2', 'absorbed PAR per unit ground area, sunlit leaves', parsunID,            &
      coordinates='layer lai')
    call RegisterFillValue(ncid, parsunID, fates_unset_r8)
    call RegisterVarAtts(ncid, 'parsha_z', (/dimIDs(1), dimIDs(2)/), type_double,       &
      'W m-2', 'absorbed PAR per unit ground area, shaded leaves', parshaID,            &
      coordinates='layer lai')
    call RegisterFillValue(ncid, parshaID, fates_unset_r8)
    call RegisterVarAtts(ncid, 'laisun_z', (/dimIDs(1), dimIDs(2)/), type_double,       &
      'm2 m-2', 'sunlit leaf area index per layer', laisunID, coordinates='layer lai')
    call RegisterFillValue(ncid, laisunID, fates_unset_r8)
    call RegisterVarAtts(ncid, 'laisha_z', (/dimIDs(1), dimIDs(2)/), type_double,       &
      'm2 m-2', 'shaded leaf area index per layer', laishaID, coordinates='layer lai')
    call RegisterFillValue(ncid, laishaID, fates_unset_r8)
    call RegisterVarAtts(ncid, 'nscaler_z', (/dimIDs(1), dimIDs(2)/), type_double, '-', &
      'nitrogen-scaling factor per layer', nscalerID, coordinates='layer lai')
    call RegisterFillValue(ncid, nscalerID, fates_unset_r8)
    call RegisterVarAtts(ncid, 'anet_z', (/dimIDs(1), dimIDs(2)/), type_double,         &
      'umolC m-2 s-1',                                                                  &
      'area-weighted net photosynthesis per unit leaf area, per layer', anetzID,        &
      coordinates='layer lai')
    call RegisterFillValue(ncid, anetzID, fates_unset_r8)

    ! swept coordinates
    call RegisterVarAtts(ncid, 'par', dimIDs(3:3), type_double,                &
      'umol m-2 s-1', 'swept incident PPFD at the top of the canopy', parID)
    call RegisterVarAtts(ncid, 'co2', dimIDs(4:4), type_double, 'Pa',          &
      'swept CO2 partial pressure', co2ID)
    call RegisterVarAtts(ncid, 'vpd', dimIDs(5:5), type_double, 'Pa',          &
      'swept leaf-to-air vapor pressure deficit', vpdID)
    call RegisterVarAtts(ncid, 'temp', dimIDs(6:6), type_double, 'K',          &
      'swept leaf temperature', tempID)
    call RegisterVarAtts(ncid, 'soilfrac', dimIDs(7:7), type_double, '-',      &
      'swept soil water content, fraction of saturation', soilfracID)

    ! diagnostics derived from the swept values
    call RegisterVarAtts(ncid, 'veg_esat_bytemp', dimIDs(6:6), type_double,    &
      'Pa', 'saturation vapor pressure at each swept leaf temperature',        &
      vegesatbytempID)
    call RegisterVarAtts(ncid, 'can_vpress_bytemp', dimIDs(6:6), type_double,  &
      'Pa', 'canopy air vapor pressure at each swept leaf temperature',        &
      canvpressbytempID)
    call RegisterVarAtts(ncid, 'veg_esat_byvpd', dimIDs(5:5), type_double,     &
      'Pa', 'saturation vapor pressure at the fixed default leaf temperature', &
      vegesatbyvpdID)
    call RegisterVarAtts(ncid, 'can_vpress_byvpd', dimIDs(5:5), type_double,   &
      'Pa', 'canopy air vapor pressure at each swept VPD', canvpressbyvpdID)
    call RegisterVarAtts(ncid, 'btran_bysoilfrac', dimIDs(7:7), type_double,   &
      '-', 'btran derived from each swept soil water content fraction',        &
      btranbysoilfracID)

    ! canopy-integrated photosynthesis, per sweep and prescribed LAI
    call RegisterVarAtts(ncid, 'canopy_anet_bypar', (/dimIDs(3), dimIDs(2)/),  &
      type_double, 'umolC m-2 s-1', 'canopy net photosynthesis vs. PAR',       &
      anetbyparID, coordinates='par lai')
    call RegisterVarAtts(ncid, 'canopy_agross_bypar', (/dimIDs(3), dimIDs(2)/),&
      type_double, 'umolC m-2 s-1', 'canopy gross photosynthesis vs. PAR',     &
      agrossbyparID, coordinates='par lai')
    call RegisterVarAtts(ncid, 'canopy_anet_byco2', (/dimIDs(4), dimIDs(2)/),  &
      type_double, 'umolC m-2 s-1', 'canopy net photosynthesis vs. CO2',       &
      anetbyco2ID, coordinates='co2 lai')
    call RegisterVarAtts(ncid, 'canopy_agross_byco2', (/dimIDs(4), dimIDs(2)/),&
      type_double, 'umolC m-2 s-1', 'canopy gross photosynthesis vs. CO2',     &
      agrossbyco2ID, coordinates='co2 lai')
    call RegisterVarAtts(ncid, 'canopy_anet_byvpd', (/dimIDs(5), dimIDs(2)/),  &
      type_double, 'umolC m-2 s-1', 'canopy net photosynthesis vs. VPD',       &
      anetbyvpdID, coordinates='vpd lai')
    call RegisterVarAtts(ncid, 'canopy_agross_byvpd', (/dimIDs(5), dimIDs(2)/),&
      type_double, 'umolC m-2 s-1', 'canopy gross photosynthesis vs. VPD',     &
      agrossbyvpdID, coordinates='vpd lai')
    call RegisterVarAtts(ncid, 'canopy_anet_bytemp', (/dimIDs(6), dimIDs(2)/), &
      type_double, 'umolC m-2 s-1',                                            &
      'canopy net photosynthesis vs. leaf temperature', anetbytempID,          &
      coordinates='temp lai')
    call RegisterVarAtts(ncid, 'canopy_agross_bytemp',                         &
      (/dimIDs(6), dimIDs(2)/), type_double, 'umolC m-2 s-1',                  &
      'canopy gross photosynthesis vs. leaf temperature', agrossbytempID,      &
      coordinates='temp lai')
    call RegisterVarAtts(ncid, 'canopy_anet_bysoilfrac',                       &
      (/dimIDs(7), dimIDs(2)/), type_double, 'umolC m-2 s-1',                  &
      'canopy net photosynthesis vs. soil water content', anetbysoilfracID,    &
      coordinates='soilfrac lai')
    call RegisterVarAtts(ncid, 'canopy_agross_bysoilfrac',                     &
      (/dimIDs(7), dimIDs(2)/), type_double, 'umolC m-2 s-1',                  &
      'canopy gross photosynthesis vs. soil water content',                    &
      agrossbysoilfracID, coordinates='soilfrac lai')

    call EndNCDef(ncid)
    call WriteVar(ncid, layerID, layer_index(:))
    call WriteVar(ncid, laiID, lai_vals(:))
    call WriteVar(ncid, nvID, nv_out(:))
    call WriteVar(ncid, canopyanetID, canopy_anet(:))
    call WriteVar(ncid, canopyagrossID, canopy_agross(:))
    call WriteVar(ncid, parsunID, parsun_z_out(:,:))
    call WriteVar(ncid, parshaID, parsha_z_out(:,:))
    call WriteVar(ncid, laisunID, laisun_z_out(:,:))
    call WriteVar(ncid, laishaID, laisha_z_out(:,:))
    call WriteVar(ncid, nscalerID, nscaler_z_out(:,:))
    call WriteVar(ncid, anetzID, anet_z_out(:,:))

    call WriteVar(ncid, parID, par_vals(:))
    call WriteVar(ncid, co2ID, co2_vals(:))
    call WriteVar(ncid, vpdID, vpd_vals(:))
    call WriteVar(ncid, tempID, temp_vals(:))
    call WriteVar(ncid, soilfracID, soilfrac_vals(:))

    call WriteVar(ncid, vegesatbytempID, veg_esat_bytemp(:))
    call WriteVar(ncid, canvpressbytempID, can_vpress_bytemp(:))
    call WriteVar(ncid, vegesatbyvpdID, veg_esat_byvpd(:))
    call WriteVar(ncid, canvpressbyvpdID, can_vpress_byvpd(:))
    call WriteVar(ncid, btranbysoilfracID, btran_bysoilfrac(:))

    call WriteVar(ncid, anetbyparID, canopy_anet_bypar(:,:))
    call WriteVar(ncid, agrossbyparID, canopy_agross_bypar(:,:))
    call WriteVar(ncid, anetbyco2ID, canopy_anet_byco2(:,:))
    call WriteVar(ncid, agrossbyco2ID, canopy_agross_byco2(:,:))
    call WriteVar(ncid, anetbyvpdID, canopy_anet_byvpd(:,:))
    call WriteVar(ncid, agrossbyvpdID, canopy_agross_byvpd(:,:))
    call WriteVar(ncid, anetbytempID, canopy_anet_bytemp(:,:))
    call WriteVar(ncid, agrossbytempID, canopy_agross_bytemp(:,:))
    call WriteVar(ncid, anetbysoilfracID, canopy_anet_bysoilfrac(:,:))
    call WriteVar(ncid, agrossbysoilfracID, canopy_agross_bysoilfrac(:,:))

    call CloseNCFile(ncid)

  end subroutine WriteOutput

end program FatesCanopyLevelPhoto
