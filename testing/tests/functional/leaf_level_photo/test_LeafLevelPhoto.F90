program FatesTestLeafLevelPhoto
  !
  ! DESCRIPTION:
  ! Leaf-level photosynthesis sensitivity sweep: five independent sweeps (PAR,
  ! CO2, VPD, leaf temperature, and soil water content), each varying ONE
  ! driver variable while holding everything else at a fixed default,
  ! evaluated for a single PFT (target_pft, 1-based, into the parameter file)
  ! via FatesTestLeafPhotoMod's EvaluateLeafPhotosynthesis (see that module's
  ! header comment for the full current production call sequence).
  !
  ! DEFAULT REFERENCE CONDITIONS: this test evaluates isolated, static points,
  ! so "default" here means a deliberately simple, site-independent reference state
  ! (25 C leaf temperature, 1 kPa leaf-to-air VPD, 380 umol/mol CO2, 210 mmol/mol O2,
  ! ~full-sun PAR, non-limiting nscaler/btran)
  !
  ! The Kumarathunge et al. (2019) temperature-acclimation model needs two running-mean
  ! reference temperatures (t_growth, t_home) in addition to the instantaneous leaf
  ! temperature being swept - both are held fixed at the same default leaf temperature
  ! for every sweep and every point, which isolates the instantaneous response of an
  ! already-acclimated leaf from any shift in the acclimation state itself
  !
  ! The soil water content sweep does NOT vary btran directly. FATES's real
  ! btran (biogeophys/EDBtranMod.F90::btran_ed) is computed from soil matric
  ! potential (smp) via, for a single unfrozen soil layer with root fraction 1
  ! (this test has no soil column, so those simplifications stand in for the
  ! multi-layer/root-weighted production sum):
  !     smp_node = max(smpsc, smp)
  !     btran    = min((smp_node - smpsc)/(smpso - smpsc), 1.0)
  ! where smpsc/smpso are the real PFT parameters (soil matric potential at
  ! full stomatal closure/opening, read from the parameter file for
  ! target_pft). This test sweeps a soil water content fraction in [0, 1] and
  ! maps it linearly onto smp between smpsc (fraction=0) and saturation,
  ! smp=0 (fraction=1): smp(frac) = smpsc*(1-frac). Feeding that into the
  ! production ramp above reproduces the flat-btran-at-1 / declining-toward-0
  ! shape most soil-moisture-stress formulations show, using only real FATES
  ! PFT parameters - no soil texture (which this test has no basis to invent)
  ! is involved.

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesConstantsMod,           only : t_water_freeze_k_1atm
  use EDPftvarcon,                 only : EDPftvarcon_inst
  use PRTParametersMod,            only : prt_params
  use PRTInitParamsFatesMod,       only : PRTDerivedParams
  use FatesParameterDerivedMod,    only : param_derived
  use FatesUnitTestParamReaderMod, only : ReadParameters
  use FatesArgumentUtils,          only : command_line_arg
  use FatesFactoryMod,             only : InitializeGlobals
  use FatesGlobals,                only : FatesGlobalsInit
  use FatesInterfaceTypesMod,      only : numpft
  use LeafBiophysicsMod,           only : lb_params
  use LeafBiophysicsMod,           only : FvCB1980, medlyn_model, net_assim_model
  use LeafBiophysicsMod,           only : photosynth_acclim_model_kumarathunge_etal_2019
  use LeafBiophysicsMod,           only : QSat
  use FatesTestEnvironmentMod,     only : environment_type
  use FatesTestLeafPhotoMod,       only : EvaluateLeafPhotosynthesis, LeafNitrogenContent
  use FatesUnitTestIOMod,          only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod,          only : WriteVar, RegisterVar, EndNCDef
  use FatesUnitTestIOMod,          only : type_double, type_int

  implicit none

  ! LOCALS:
  character(len=:), allocatable :: param_file ! input parameter file
  type(environment_type)        :: env        ! prescribed atmospheric boundary conditions
  real(r8), allocatable         :: lnc_top(:)  ! leaf N content at the canopy top, per pft [gN/m2 leaf]
  real(r8) :: vcmax25top_pft ! reference (25C, canopy-top) maximum carboxylation rate, target_pft's flat PFT default [umol/m2/s]
  real(r8) :: jmax25top_pft  ! reference (25C, canopy-top) maximum electron transport rate, target_pft's flat PFT default [umol/m2/s]
  real(r8) :: kp25top_pft    ! reference (25C, canopy-top) initial slope of C4 CO2 response, target_pft's flat PFT default [umol/m2/s]

  ! swept-variable value arrays
  real(r8), allocatable :: par_vals(:)      ! swept PAR values [umol/m2/s]
  real(r8), allocatable :: co2_vals(:)      ! swept CO2 partial pressure values [Pa]
  real(r8), allocatable :: vpd_vals(:)      ! swept leaf-to-air VPD values [Pa]
  real(r8), allocatable :: temp_vals(:)     ! swept leaf temperature values [K]
  real(r8), allocatable :: soilfrac_vals(:) ! swept soil water content, fraction of saturation [0-1]

  ! derived vapor-pressure diagnostics, dimensioned like their sweep (only
  ! temperature and VPD change vapor pressure at all)
  real(r8), allocatable :: veg_esat_bytemp(:)   ! saturation vapor pressure at each swept leaf temperature [Pa]
  real(r8), allocatable :: can_vpress_bytemp(:) ! canopy air vapor pressure at each swept leaf temperature, fixed default VPD [Pa]
  real(r8), allocatable :: veg_esat_byvpd(:)    ! saturation vapor pressure at the fixed default leaf temperature [Pa] (constant across this sweep - included for symmetry with veg_esat_bytemp)
  real(r8), allocatable :: can_vpress_byvpd(:)  ! canopy air vapor pressure at each swept VPD (= veg_esat - vpd), fixed default leaf temperature [Pa]
  real(r8), allocatable :: btran_bysoilfrac(:)  ! btran derived from each swept soil water content fraction (see header comment) [0-1]

  ! sweep output, each dimensioned (n_sweep_points, numpft)
  real(r8), allocatable :: anet_bypar(:,:), agross_bypar(:,:), gs_bypar(:,:), ci_bypar(:,:)
  real(r8), allocatable :: anet_byco2(:,:), agross_byco2(:,:), gs_byco2(:,:), ci_byco2(:,:)
  real(r8), allocatable :: anet_byvpd(:,:), agross_byvpd(:,:), gs_byvpd(:,:), ci_byvpd(:,:)
  real(r8), allocatable :: anet_bytemp(:,:), agross_bytemp(:,:), gs_bytemp(:,:), ci_bytemp(:,:)
  real(r8), allocatable :: anet_bysoilfrac(:,:), agross_bysoilfrac(:,:), gs_bysoilfrac(:,:), ci_bysoilfrac(:,:)

  integer :: n_par, n_co2, n_vpd, n_temp, n_soilfrac ! sweep array sizes
  integer :: i ! looping index

  ! target_pft's soil matric potential thresholds (mm, negative), read from
  ! the parameter file below and used only by the soil water content sweep
  real(r8) :: smpsc_pft ! soil matric potential at full stomatal closure [mm]
  real(r8) :: smpso_pft ! soil matric potential at full stomatal opening [mm]

  ! CONSTANTS:
  ! PFT 6 is broadleaf_colddecid_extratrop_tree, matching Rogers et al.
  ! (2017)'s stated "generic temperate broad leaved deciduous tree".
  ! test_CanopyLevelPhoto.F90 uses this same PFT so the leaf and canopy
  ! panels describe the same plant and can be shown as a pair - see that
  ! file's header for why the choice matters much more at canopy scale
  integer, parameter :: target_pft = 6 ! PFT index to evaluate (1-based)

  ! default reference conditions - CO2/O2/dayl_factor/btran come from
  ! FatesTestEnvironmentMod's shared reference-atmosphere defaults (env%Init,
  ! below) rather than being duplicated here; this test's own defaults are
  ! only the ones genuinely specific to it (leaf temperature/VPD, both fixed
  ! independently of env's site-climatology/RH-driven state, and nscaler,
  ! which env has no equivalent of)
  real(r8), parameter :: default_veg_tempk   = 25.0_r8 + t_water_freeze_k_1atm ! [K]
  real(r8), parameter :: default_vpd         = 1000.0_r8  ! [Pa] leaf-to-air VPD, esat(Tleaf) - eair
  real(r8), parameter :: default_par         = 1500.0_r8  ! [umol/m2/s]
  real(r8), parameter :: default_nscaler     = 1.0_r8     ! [0-1]

  ! sweep ranges
  real(r8), parameter :: min_temp = 5.0_r8,    max_temp = 40.0_r8,    temp_inc = 0.5_r8   ! [degC]
  real(r8), parameter :: min_par  = 0.0_r8,    max_par  = 1600.0_r8,  par_inc  = 5.0_r8   ! [umol/m2/s]
  real(r8), parameter :: min_vpd  = 500.0_r8,  max_vpd  = 2500.0_r8,  vpd_inc  = 20.0_r8  ! [Pa] (0.5-2.5 kPa)
  real(r8), parameter :: min_co2  = 250.0_r8,  max_co2  = 1000.0_r8,  co2_inc  = 5.0_r8   ! [umol/mol]
  real(r8), parameter :: soilfrac_inc = 0.02_r8 ! [0-1]

  character(len=*), parameter :: out_file = 'leaf_level_photo_out.nc' ! output file

  param_file = command_line_arg(1)
  call ReadParameters(param_file)

  smpsc_pft = EDPftvarcon_inst%smpsc(target_pft)
  smpso_pft = EDPftvarcon_inst%smpso(target_pft)

  ! step_size is unused by anything this test exercises (no time-stepping occurs),
  ! but InitializeGlobals requires a value
  call InitializeGlobals(86400.0_r8)
  numpft = size(prt_params%wood_density, dim=1)
  call FatesGlobalsInit(6, .false.)
  call PRTDerivedParams()

  ! host-model-namelist-controlled leaf biophysics switches - matches
  ! test_SingleCohort.F90's choices for consistency across this test suite
  lb_params%electron_transport_model = FvCB1980
  lb_params%stomatal_model           = medlyn_model
  lb_params%stomatal_assim_model     = net_assim_model
  lb_params%photo_tempsens_model     = photosynth_acclim_model_kumarathunge_etal_2019

  ! leaf N content for target_pft - constant for the whole run
  allocate(lnc_top(1))
  lnc_top(1) = LeafNitrogenContent(target_pft)

  ! reference (25C, canopy-top) photosynthetic capacity for target_pft - the
  ! flat PFT default (no cohort/acclimation state exists in this driver) -
  ! constant for the whole run
  vcmax25top_pft = EDPftvarcon_inst%vcmax25top(target_pft,1)
  jmax25top_pft  = param_derived%jmax25top(target_pft,1)
  kp25top_pft    = param_derived%kp25top(target_pft,1)

  ! atmospheric constants not swept by this test (canopy pressure, background
  ! CO2/O2 partial pressure, leaf boundary-layer conductance, dayl_factor,
  ! btran - see FatesTestEnvironmentMod's shared reference-atmosphere defaults)
  call env%Init()

  ! ---------------------------------------------------------------------------------------
  ! build the swept-value arrays
  ! ---------------------------------------------------------------------------------------

  n_par = int((max_par - min_par)/par_inc) + 1
  n_co2 = int((max_co2 - min_co2)/co2_inc) + 1
  n_vpd  = int((max_vpd - min_vpd)/vpd_inc) + 1
  n_temp = int((max_temp - min_temp)/temp_inc) + 1
  n_soilfrac = int((1.0_r8 - 0.0_r8)/soilfrac_inc) + 1

  allocate(par_vals(n_par))
  allocate(co2_vals(n_co2))
  allocate(vpd_vals(n_vpd))
  allocate(temp_vals(n_temp))
  allocate(soilfrac_vals(n_soilfrac))

  do i = 1, n_par
    par_vals(i) = min_par + par_inc*real(i-1, r8)
  end do
  do i = 1, n_co2
    ! ppm -> Pa, at the default (Init()-prescribed) canopy pressure
    co2_vals(i) = ((min_co2 + co2_inc*real(i-1, r8))/1.0e6_r8) * env%can_press
  end do
  do i = 1, n_vpd
    vpd_vals(i) = min_vpd + vpd_inc*real(i-1, r8)
  end do
  do i = 1, n_temp
    temp_vals(i) = (min_temp + temp_inc*real(i-1, r8)) + t_water_freeze_k_1atm
  end do
  do i = 1, n_soilfrac
    soilfrac_vals(i) = 0.0_r8 + soilfrac_inc*real(i-1, r8)
  end do

  allocate(veg_esat_bytemp(n_temp), can_vpress_bytemp(n_temp))
  allocate(veg_esat_byvpd(n_vpd), can_vpress_byvpd(n_vpd))
  allocate(btran_bysoilfrac(n_soilfrac))

  allocate(anet_bypar(n_par, 1), agross_bypar(n_par, 1), gs_bypar(n_par, 1), ci_bypar(n_par, 1))
  allocate(anet_byco2(n_co2, 1), agross_byco2(n_co2, 1), gs_byco2(n_co2, 1), ci_byco2(n_co2, 1))
  allocate(anet_byvpd(n_vpd, 1), agross_byvpd(n_vpd, 1), gs_byvpd(n_vpd, 1), ci_byvpd(n_vpd, 1))
  allocate(anet_bytemp(n_temp, 1), agross_bytemp(n_temp, 1), gs_bytemp(n_temp, 1), ci_bytemp(n_temp, 1))
  allocate(anet_bysoilfrac(n_soilfrac, 1), agross_bysoilfrac(n_soilfrac, 1), gs_bysoilfrac(n_soilfrac, 1), ci_bysoilfrac(n_soilfrac, 1))

  ! a shared default vapor-pressure state (leaf temperature and VPD both at
  ! their defaults) - reused as the fixed background for every sweep except
  ! the temperature and VPD sweeps themselves, which recompute it per point
  block
    real(r8) :: default_veg_esat, default_can_vpress
    real(r8) :: qs_dummy ! saturation specific humidity output from QSat
    real(r8) :: smp, smp_node ! soil matric potential at a swept soil water content fraction, and its full-closure-clamped value [mm]

    call QSat(default_veg_tempk, env%can_press, qs_dummy, default_veg_esat)
    default_can_vpress = default_veg_esat - default_vpd

    ! ---------------------------------------------------------------------
    ! PAR sweep
    ! ---------------------------------------------------------------------
    print *, '----------------------------------------------------------------'
    print *, 'Exercising leaf photosynthesis for PAR'
    do i = 1, n_par
      call EvaluateLeafPhotosynthesis(target_pft, par_vals(i), default_veg_tempk, &
        default_veg_tempk, default_veg_tempk, env%can_press, env%can_co2_ppress, &
        env%can_o2_ppress, default_veg_esat, default_can_vpress, env%gb,         &
        default_nscaler, env%dayl_factor, env%btran, vcmax25top_pft,             &
        jmax25top_pft, kp25top_pft, lnc_top(1),                                  &
        agross_bypar(i,1), anet_bypar(i,1), gs_bypar(i,1), ci_bypar(i,1))
    end do

    ! ---------------------------------------------------------------------
    ! CO2 sweep
    ! ---------------------------------------------------------------------
    print *, '----------------------------------------------------------------'
    print *, 'Exercising leaf photosynthesis for CO2'
    do i = 1, n_co2
      call EvaluateLeafPhotosynthesis(target_pft, default_par, default_veg_tempk, &
        default_veg_tempk, default_veg_tempk, env%can_press, co2_vals(i),        &
        env%can_o2_ppress, default_veg_esat, default_can_vpress, env%gb,         &
        default_nscaler, env%dayl_factor, env%btran, vcmax25top_pft,             &
        jmax25top_pft, kp25top_pft, lnc_top(1),                                  &
        agross_byco2(i,1), anet_byco2(i,1), gs_byco2(i,1), ci_byco2(i,1))
    end do

    ! ---------------------------------------------------------------------
    ! VPD sweep - leaf temperature fixed at the default, so veg_esat is
    ! constant and can_vpress is derived directly from the swept VPD
    ! ---------------------------------------------------------------------
    print *, '----------------------------------------------------------------'
    print *, 'Exercising leaf photosynthesis for VPD'
    do i = 1, n_vpd
      veg_esat_byvpd(i) = default_veg_esat
      can_vpress_byvpd(i) = default_veg_esat - vpd_vals(i)
      call EvaluateLeafPhotosynthesis(target_pft, default_par, default_veg_tempk, &
        default_veg_tempk, default_veg_tempk, env%can_press, env%can_co2_ppress, &
        env%can_o2_ppress, veg_esat_byvpd(i), can_vpress_byvpd(i), env%gb,       &
        default_nscaler, env%dayl_factor, env%btran, vcmax25top_pft,             &
        jmax25top_pft, kp25top_pft, lnc_top(1),                                  &
        agross_byvpd(i,1), anet_byvpd(i,1), gs_byvpd(i,1), ci_byvpd(i,1))
    end do

    ! ---------------------------------------------------------------------
    ! leaf temperature sweep - t_growth/t_home held at the default leaf
    ! temperature throughout; VPD held fixed at the default as
    ! leaf temperature varies
    ! ---------------------------------------------------------------------
    print *, '----------------------------------------------------------------'
    print *, 'Exercising leaf photosynthesis for leaf temperature'
    do i = 1, n_temp
      call QSat(temp_vals(i), env%can_press, qs_dummy, veg_esat_bytemp(i))
      can_vpress_bytemp(i) = veg_esat_bytemp(i) - default_vpd
      call EvaluateLeafPhotosynthesis(target_pft, default_par, temp_vals(i),     &
        default_veg_tempk, default_veg_tempk, env%can_press, env%can_co2_ppress, &
        env%can_o2_ppress, veg_esat_bytemp(i), can_vpress_bytemp(i), env%gb,     &
        default_nscaler, env%dayl_factor, env%btran, vcmax25top_pft,             &
        jmax25top_pft, kp25top_pft, lnc_top(1),                                  &
        agross_bytemp(i,1), anet_bytemp(i,1), gs_bytemp(i,1), ci_bytemp(i,1))
    end do

    ! ---------------------------------------------------------------------
    ! soil water content sweep - btran derived from the real smpsc/smpso
    ! ramp at each swept fraction 
    ! ---------------------------------------------------------------------
    print *, '----------------------------------------------------------------'
    print *, 'Exercising leaf photosynthesis for soil water content'
    do i = 1, n_soilfrac
      smp = smpsc_pft * (1.0_r8 - soilfrac_vals(i))
      smp_node = max(smpsc_pft, smp)
      btran_bysoilfrac(i) = min((smp_node - smpsc_pft) / (smpso_pft - smpsc_pft), 1.0_r8)
      call EvaluateLeafPhotosynthesis(target_pft, default_par, default_veg_tempk, &
        default_veg_tempk, default_veg_tempk, env%can_press, env%can_co2_ppress, &
        env%can_o2_ppress, default_veg_esat, default_can_vpress, env%gb,         &
        default_nscaler, env%dayl_factor, btran_bysoilfrac(i), vcmax25top_pft,   &
        jmax25top_pft, kp25top_pft, lnc_top(1),                                  &
        agross_bysoilfrac(i,1), anet_bysoilfrac(i,1), gs_bysoilfrac(i,1), ci_bysoilfrac(i,1))
    end do
  end block

  ! ---------------------------------------------------------------------------------------
  ! write output
  ! ---------------------------------------------------------------------------------------

  call WriteOutput()

contains

  ! ==========================================================================

  subroutine WriteOutput()
    !
    ! DESCRIPTION:
    ! Writes every sweep's swept values, derived vapor-pressure/btran
    ! diagnostics (temperature/VPD/soil-water-content sweeps only), and
    ! per-pft agross/anet/gs/ci to netCDF.

    ! LOCALS:
    integer, allocatable :: pft_indices(:) ! pft index coordinate values
    integer              :: ncid           ! netcdf file id
    character(len=20)    :: dim_names(6)   ! dimension names
    integer              :: dimIDs(6)      ! dimension IDs
    integer              :: parID, co2ID, vpdID, tempID, soilfracID, pftID
    integer              :: vegesatbytempID, canvpressbytempID
    integer              :: vegesatbyvpdID, canvpressbyvpdID
    integer              :: btranbysoilfracID
    integer              :: agrossbyparID, anetbyparID, gsbyparID, cibyparID
    integer              :: agrossbyco2ID, anetbyco2ID, gsbyco2ID, cibyco2ID
    integer              :: agrossbyvpdID, anetbyvpdID, gsbyvpdID, cibyvpdID
    integer              :: agrossbytempID, anetbytempID, gsbytempID, cibytempID
    integer              :: agrossbysoilfracID, anetbysoilfracID, gsbysoilfracID, cibysoilfracID

    allocate(pft_indices(1))
    pft_indices(1) = target_pft

    dim_names = [character(len=20) :: 'par', 'co2', 'vpd', 'temp', 'soilfrac', 'pft']

    call OpenNCFile(trim(out_file), ncid, 'readwrite')

    call RegisterNCDims(ncid, dim_names, (/n_par, n_co2, n_vpd, n_temp, n_soilfrac, 1/), 6, dimIDs)

    call RegisterVar(ncid, dim_names(1), dimIDs(1:1), type_double,               &
      [character(len=20) :: 'units', 'long_name'],                              &
      [character(len=150) :: 'umol m-2 s-1', 'swept incident PAR'], 2, parID)
    call RegisterVar(ncid, dim_names(2), dimIDs(2:2), type_double,               &
      [character(len=20) :: 'units', 'long_name'],                              &
      [character(len=150) :: 'Pa', 'swept CO2 partial pressure'], 2, co2ID)
    call RegisterVar(ncid, dim_names(3), dimIDs(3:3), type_double,               &
      [character(len=20) :: 'units', 'long_name'],                              &
      [character(len=150) :: 'Pa', 'swept leaf-to-air vapor pressure deficit'], 2, vpdID)
    call RegisterVar(ncid, dim_names(4), dimIDs(4:4), type_double,               &
      [character(len=20) :: 'units', 'long_name'],                              &
      [character(len=150) :: 'K', 'swept leaf temperature'], 2, tempID)
    call RegisterVar(ncid, dim_names(5), dimIDs(5:5), type_double,               &
      [character(len=20) :: 'units', 'long_name'],                              &
      [character(len=150) :: '-', 'swept soil water content, fraction of saturation'], 2, soilfracID)
    call RegisterVar(ncid, dim_names(6), dimIDs(6:6), type_int,                  &
      [character(len=20) :: 'units', 'long_name'],                              &
      [character(len=150) :: '-', 'plant functional type'], 2, pftID)

    call RegisterVar(ncid, 'veg_esat_bytemp', dimIDs(4:4), type_double,          &
      [character(len=20) :: 'units', 'long_name'],                              &
      [character(len=150) :: 'Pa', 'saturation vapor pressure at each swept leaf temperature'], 2, vegesatbytempID)
    call RegisterVar(ncid, 'can_vpress_bytemp', dimIDs(4:4), type_double,        &
      [character(len=20) :: 'units', 'long_name'],                              &
      [character(len=150) :: 'Pa', 'canopy air vapor pressure at each swept leaf temperature, fixed default VPD'], 2, canvpressbytempID)
    call RegisterVar(ncid, 'veg_esat_byvpd', dimIDs(3:3), type_double,            &
      [character(len=20) :: 'units', 'long_name'],                              &
      [character(len=150) :: 'Pa', 'saturation vapor pressure at the fixed default leaf temperature'], 2, vegesatbyvpdID)
    call RegisterVar(ncid, 'can_vpress_byvpd', dimIDs(3:3), type_double,          &
      [character(len=20) :: 'units', 'long_name'],                              &
      [character(len=150) :: 'Pa', 'canopy air vapor pressure at each swept VPD, fixed default leaf temperature'], 2, canvpressbyvpdID)
    call RegisterVar(ncid, 'btran_bysoilfrac', dimIDs(5:5), type_double,          &
      [character(len=20) :: 'units', 'long_name'],                              &
      [character(len=150) :: '-', 'btran derived from each swept soil water content fraction'], 2, btranbysoilfracID)

    call RegisterVar(ncid, 'agross_bypar', (/dimIDs(1), dimIDs(6)/), type_double, &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'par pft', 'umolC m-2 s-1', 'gross photosynthesis vs. PAR'], 3, agrossbyparID)
    call RegisterVar(ncid, 'anet_bypar', (/dimIDs(1), dimIDs(6)/), type_double,  &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'par pft', 'umolC m-2 s-1', 'net photosynthesis vs. PAR'], 3, anetbyparID)
    call RegisterVar(ncid, 'gs_bypar', (/dimIDs(1), dimIDs(6)/), type_double,    &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'par pft', 'umol H2O m-2 s-1', 'stomatal conductance vs. PAR'], 3, gsbyparID)
    call RegisterVar(ncid, 'ci_bypar', (/dimIDs(1), dimIDs(6)/), type_double,    &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'par pft', 'Pa', 'intracellular CO2 vs. PAR'], 3, cibyparID)

    call RegisterVar(ncid, 'agross_byco2', (/dimIDs(2), dimIDs(6)/), type_double, &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'co2 pft', 'umolC m-2 s-1', 'gross photosynthesis vs. CO2'], 3, agrossbyco2ID)
    call RegisterVar(ncid, 'anet_byco2', (/dimIDs(2), dimIDs(6)/), type_double,  &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'co2 pft', 'umolC m-2 s-1', 'net photosynthesis vs. CO2'], 3, anetbyco2ID)
    call RegisterVar(ncid, 'gs_byco2', (/dimIDs(2), dimIDs(6)/), type_double,    &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'co2 pft', 'umol H2O m-2 s-1', 'stomatal conductance vs. CO2'], 3, gsbyco2ID)
    call RegisterVar(ncid, 'ci_byco2', (/dimIDs(2), dimIDs(6)/), type_double,    &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'co2 pft', 'Pa', 'intracellular CO2 vs. CO2'], 3, cibyco2ID)

    call RegisterVar(ncid, 'agross_byvpd', (/dimIDs(3), dimIDs(6)/), type_double, &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'vpd pft', 'umolC m-2 s-1', 'gross photosynthesis vs. VPD'], 3, agrossbyvpdID)
    call RegisterVar(ncid, 'anet_byvpd', (/dimIDs(3), dimIDs(6)/), type_double,   &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'vpd pft', 'umolC m-2 s-1', 'net photosynthesis vs. VPD'], 3, anetbyvpdID)
    call RegisterVar(ncid, 'gs_byvpd', (/dimIDs(3), dimIDs(6)/), type_double,     &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'vpd pft', 'umol H2O m-2 s-1', 'stomatal conductance vs. VPD'], 3, gsbyvpdID)
    call RegisterVar(ncid, 'ci_byvpd', (/dimIDs(3), dimIDs(6)/), type_double,     &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'vpd pft', 'Pa', 'intracellular CO2 vs. VPD'], 3, cibyvpdID)

    call RegisterVar(ncid, 'agross_bytemp', (/dimIDs(4), dimIDs(6)/), type_double, &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'temperature pft', 'umolC m-2 s-1', 'gross photosynthesis vs. leaf temperature'], 3, agrossbytempID)
    call RegisterVar(ncid, 'anet_bytemp', (/dimIDs(4), dimIDs(6)/), type_double, &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'temperature pft', 'umolC m-2 s-1', 'net photosynthesis vs. leaf temperature'], 3, anetbytempID)
    call RegisterVar(ncid, 'gs_bytemp', (/dimIDs(4), dimIDs(6)/), type_double,   &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'temperature pft', 'umol H2O m-2 s-1', 'stomatal conductance vs. leaf temperature'], 3, gsbytempID)
    call RegisterVar(ncid, 'ci_bytemp', (/dimIDs(4), dimIDs(6)/), type_double,   &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'temperature pft', 'Pa', 'intracellular CO2 vs. leaf temperature'], 3, cibytempID)

    call RegisterVar(ncid, 'agross_bysoilfrac', (/dimIDs(5), dimIDs(6)/), type_double, &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'soilfrac pft', 'umolC m-2 s-1', 'gross photosynthesis vs. soil water content'], 3, agrossbysoilfracID)
    call RegisterVar(ncid, 'anet_bysoilfrac', (/dimIDs(5), dimIDs(6)/), type_double, &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'soilfrac pft', 'umolC m-2 s-1', 'net photosynthesis vs. soil water content'], 3, anetbysoilfracID)
    call RegisterVar(ncid, 'gs_bysoilfrac', (/dimIDs(5), dimIDs(6)/), type_double, &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'soilfrac pft', 'umol H2O m-2 s-1', 'stomatal conductance vs. soil water content'], 3, gsbysoilfracID)
    call RegisterVar(ncid, 'ci_bysoilfrac', (/dimIDs(5), dimIDs(6)/), type_double, &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],                &
      [character(len=150) :: 'soilfrac pft', 'Pa', 'intracellular CO2 vs. soil water content'], 3, cibysoilfracID)

    call EndNCDef(ncid)

    call WriteVar(ncid, parID, par_vals(:))
    call WriteVar(ncid, co2ID, co2_vals(:))
    call WriteVar(ncid, vpdID, vpd_vals(:))
    call WriteVar(ncid, tempID, temp_vals(:))
    call WriteVar(ncid, soilfracID, soilfrac_vals(:))
    call WriteVar(ncid, pftID, pft_indices(:))

    call WriteVar(ncid, vegesatbytempID, veg_esat_bytemp(:))
    call WriteVar(ncid, canvpressbytempID, can_vpress_bytemp(:))
    call WriteVar(ncid, vegesatbyvpdID, veg_esat_byvpd(:))
    call WriteVar(ncid, canvpressbyvpdID, can_vpress_byvpd(:))
    call WriteVar(ncid, btranbysoilfracID, btran_bysoilfrac(:))

    call WriteVar(ncid, agrossbyparID, agross_bypar(:,:))
    call WriteVar(ncid, anetbyparID, anet_bypar(:,:))
    call WriteVar(ncid, gsbyparID, gs_bypar(:,:))
    call WriteVar(ncid, cibyparID, ci_bypar(:,:))

    call WriteVar(ncid, agrossbyco2ID, agross_byco2(:,:))
    call WriteVar(ncid, anetbyco2ID, anet_byco2(:,:))
    call WriteVar(ncid, gsbyco2ID, gs_byco2(:,:))
    call WriteVar(ncid, cibyco2ID, ci_byco2(:,:))

    call WriteVar(ncid, agrossbyvpdID, agross_byvpd(:,:))
    call WriteVar(ncid, anetbyvpdID, anet_byvpd(:,:))
    call WriteVar(ncid, gsbyvpdID, gs_byvpd(:,:))
    call WriteVar(ncid, cibyvpdID, ci_byvpd(:,:))

    call WriteVar(ncid, agrossbytempID, agross_bytemp(:,:))
    call WriteVar(ncid, anetbytempID, anet_bytemp(:,:))
    call WriteVar(ncid, gsbytempID, gs_bytemp(:,:))
    call WriteVar(ncid, cibytempID, ci_bytemp(:,:))

    call WriteVar(ncid, agrossbysoilfracID, agross_bysoilfrac(:,:))
    call WriteVar(ncid, anetbysoilfracID, anet_bysoilfrac(:,:))
    call WriteVar(ncid, gsbysoilfracID, gs_bysoilfrac(:,:))
    call WriteVar(ncid, cibysoilfracID, ci_bysoilfrac(:,:))

    call CloseNCFile(ncid)

  end subroutine WriteOutput

end program FatesTestLeafLevelPhoto
