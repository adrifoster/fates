program FatesSingleCohort
  !
  ! DESCRIPTION:
  ! Fixed-light, single-cohort driver: a sweep over prescribed incident light levels.
  ! For each light level, a single cohort (no patch or site) is built at recruitment
  ! size and stepped forward through a year/day/sub-daily loop nest (RunOneLightLevel,
  ! below). Each light level is an independent trajectory - the cohort is created
  ! fresh at recruitment size and destroyed within the light-level loop, so levels
  ! cannot interact.
  !
  ! The driver's supporting physics/bookkeeping is split across a few small modules,
  ! so that this file can stay focused on the call order of one simulated day:
  !   - FatesTestEnvironmentMod - the prescribed atmospheric and soil boundary
  !     conditions (temperature, pressure, CO2/O2, humidity, boundary-layer
  !     conductance, soil moisture/rooting), held fixed for the entire run.
  !   - FatesTestLightEnvMod - the prescribed light environment (reference full-sun
  !     PAR, direct/diffuse split, diurnal shape), attenuated through the cohort's
  !     own leaf layers via FATES's two-stream radiation solver.
  !   - FatesTestCohortPhysMod - the per-leaf-layer photosynthetic capacity/dark-
  !     respiration working arrays, the once-per-day setup that refreshes them, and
  !     the per-substep carbon uptake step (leaf photosynthesis, scaled up to a
  !     per-individual carbon flux mirroring FatesPlantRespPhotosynthMod's
  !     ScaleLeafLayerFluxToCohort - private to that module, so reimplemented here
  !     rather than reused - plus non-leaf maintenance respiration via
  !     NonleafMaintenanceRespiration).
  !   - FatesTestHistoryMod - accumulates the output time series/light-profile
  !     snapshots in memory and writes them to netCDF once at the end.
  !
  ! NSC/PRT allocation runs once per day, after the sub-daily loop, in
  ! RunOneLightLevel's call to DailyGrowthAndMortality (below): growth respiration
  ! is taxed off the day's net carbon (matching EDMainMod.F90's resp_g_acc_hold
  ! formula), the result is handed to PARTEH via cohort%npp_acc (the real
  ! carbon_balance boundary condition PARTEH reads - see FatesCohortMod.F90's
  ! InitPRTBoundaryConditions), and PRTMaintTurnover/DailyPRT(phase=1,2,3) run in the
  ! same order as EDMainMod.F90's daily dynamics sequence. Crown area (cohort%c_area)
  ! is refreshed daily via carea_allom, alongside treelai/treesai/nv via
  ! UpdateCohortLAI - production only ever recomputes crown area as part of
  ! patch-level canopy structure (EDCanopyStructureMod), which this patch-less driver
  ! has no equivalent of, so it is called directly here instead.
  !
  ! Carbon starvation mortality is wired in via the cmort term from
  ! EDMortalityFunctionsMod.F90's mortality_rates (linear model, matching CTSM's
  ! default fates_cstarvation_model='linear' namelist setting) - the only mortality
  ! term implemented here, since the others (background, hydraulic, freezing,
  ! senescence, damage) are either PFT-constant or depend on site/patch state this
  ! driver doesn't have. cmort is applied as a continuous daily proportional
  ! reduction of cohort%n (an Euler step on the annual rate), matching
  ! EDMainMod.F90's cohort%n update - not a deterministic or stochastic kill, since
  ! FATES cohorts represent a population, not an individual. There is no
  ! disturbance mechanism in this patch-less driver, so (unlike production's
  ! canopy-layer-1 cohorts) the full rate always goes to %n rather than partly
  ! being rerouted to spawn a new patch.
  !
  ! OUTPUT: two groups of variables are accumulated in memory (FatesTestHistoryMod)
  ! and written once at the end:
  !   - Daily whole-cohort time series, dimensioned (time, light_level), where time is
  !     the day index within each light level's independent trajectory (1 to
  !     nyears*days_per_year): dbh, height, treelai, crown area, the PARTEH carbon
  !     pools (leaf/fine root/sapwood/structure/storage), daily net carbon, daily GPP,
  !     the carbon balance handed to PARTEH (npp_acc, net of growth respiration),
  !     frac_store (storage as a fraction of target leaf carbon), the carbon
  !     starvation mortality rate (cmort, indiv/year), and cohort number density (n).
  !   - A per-leaf-layer light profile (parsun_z/parsha_z/laisun_z/laisha_z),
  !     dimensioned (nlevleaf, year, light_level), captured once per year (the first
  !     day of each year, at the sub-daily step nearest solar noon - matching
  !     single_cohort_test.py's established noon convention - since capturing this
  !     daily would make the output file awkwardly large: nlevleaf x n_days_total x
  !     n_light_levels x 4 variables x 8 bytes is ~90MB, vs. ~0.25MB captured
  !     annually). nlevleaf (EDParamsMod) is a compile-time maximum leaf+stem layer
  !     count; layers above a given day's cohort%nv are left at the fates_unset_r8
  !     fill value (registered as each variable's _FillValue attribute), so the array
  !     survives nv changing over time and differing across light levels.
  !

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesConstantsMod,           only : nearzero
  use FatesConstantsMod,           only : sec_per_day
  use FatesConstantsMod,           only : cstarvation_model_lin, cstarvation_model_exp
  use EDParamsMod,                 only : nclmax
  use EDParamsMod,                 only : nlevleaf
  use EDParamsMod,                 only : GetNVegLayers
  use EDPftvarcon,                 only : EDPftvarcon_inst
  use EDTypesMod,                  only : init_spread_inventory
  use FatesAllometryMod,           only : h2d_allom
  use FatesAllometryMod,           only : h_allom
  use FatesAllometryMod,           only : carea_allom
  use FatesCohortMod,              only : fates_cohort_type
  use EDCanopyStructureMod,        only : UpdateCohortLAI
  use EDCohortDynamicsMod,         only : EvaluateAndCorrectDBH
  use PRTLossFluxesMod,            only : PRTMaintTurnover
  use FatesUnitTestParamReaderMod, only : ReadParameters
  use FatesArgumentUtils,          only : command_line_arg
  use FatesFactoryMod,             only : InitializeGlobals, CohortFactory
  use FatesGlobals,                only : FatesGlobalsInit
  use FatesInterfaceTypesMod,      only : numpft, hlm_mort_cstarvation_model
  use PRTParametersMod,            only : prt_params
  use PRTGenericMod,                only : leaf_organ, store_organ, carbon12_element
  use PRTInitParamsFatesMod,       only : PRTDerivedParams
  use FatesTwoStreamUtilsMod,      only : TransferRadParams
  use FatesTestEnvironmentMod,     only : environment_type
  use FatesTestLightEnvMod,        only : light_env_type
  use FatesTestCohortPhysMod,      only : cohort_phys_type
  use FatesTestHistoryMod,         only : history_type
  use LeafBiophysicsMod,           only : GetCanopyGasParameters
  use LeafBiophysicsMod,           only : lb_params
  use LeafBiophysicsMod,           only : FvCB1980, medlyn_model
  use LeafBiophysicsMod,           only : net_assim_model, photosynth_acclim_model_kumarathunge_etal_2019

  implicit none

  ! LOCALS:
  character(len=:), allocatable    :: param_file    ! input parameter file
  type(environment_type)           :: env           ! prescribed atmospheric/soil boundary conditions, fixed for the whole run
  type(history_type)               :: hist          ! output accumulator/writer, shared across the light sweep
  real(r8)                         :: dbh_recruit   ! recruitment-size dbh [cm]
  real(r8), allocatable            :: light_frac(:) ! swept incident light fractions [0-1]
  integer                          :: ilight        ! light-level looping index

  ! canopy gas parameters and leaf N content - constant for the whole run since the
  ! atmospheric boundary conditions never change
  real(r8) :: mm_kco2, mm_ko2, co2_cpoint ! Michaelis-Menten constants for CO2/O2, CO2 compensation point [Pa]
  real(r8) :: lnc_top                     ! leaf N content at the canopy top [gN/m2 leaf]

  ! Ci-solver diagnostics: LeafLayerPhotosynthesis's Newton loop bails out to the
  ! much more expensive CiBisection after newton_max_iters - default parameters
  ! converge fast, but a Latin Hypercube parameter sweep may land on corners of
  ! parameter space that don't, so wall time won't necessarily scale linearly from
  ! a timing run against the default parameter file. Tracked per light level (in
  ! RunOneLightLevel) and overall (here) so the fallback rate is visible while
  ! prototyping.
  integer :: n_photo_calls_total, n_bisection_calls_total, max_solve_iter_total ! whole run

  ! CONSTANTS:
  integer,  parameter :: pft = 1                             ! plant functional type to simulate
  real(r8), dimension(nclmax), parameter :: can_tlai = 0.0_r8 ! canopy-layer LAI above the cohort; kept at zero for the whole run - shading comes only from the prescribed light fraction below, not a fictitious overstory
  real(r8), parameter :: patch_area = 1.0e4_r8               ! reference ground area the cohort occupies [m2]
  real(r8), parameter :: n_indiv = 1.0_r8                    ! number of individuals in the cohort

  ! light sweep
  integer,  parameter :: n_light_levels = 25       ! number of incident light levels to sweep
  real(r8), parameter :: light_frac_min = 0.005_r8 ! lowest incident light fraction swept [fraction of full sun]
  real(r8), parameter :: light_frac_max = 1.0_r8   ! highest incident light fraction swept [fraction of full sun]

  ! time stepping
  integer,  parameter :: n_substeps_per_day = 24                     ! sub-daily steps per day (hourly)
  real(r8), parameter :: step_size = 86400.0_r8/n_substeps_per_day   ! model time step [s]
  integer,  parameter :: days_per_year = 365                         ! days per simulated year
  integer,  parameter :: nyears = 10                                 ! number of years to simulate
  integer,  parameter :: n_days_total = nyears * days_per_year       ! total days per light level's trajectory
  integer,  parameter :: noon_substep = n_substeps_per_day/2 + 1     ! sub-daily index nearest solar noon (hour 12.5 of 24 hourly substeps) - matches single_cohort_test.py's _NOON_HOUR_IDX convention

  ! output file
  character(len=*), parameter :: out_file = 'single_cohort_out.nc' ! output file

  ! read in parameter file name from command line
  param_file = command_line_arg(1)

  ! read in parameter file
  call ReadParameters(param_file)

  ! initialize global PRT/allometry data needed by the cohort machinery
  call InitializeGlobals(step_size)

  ! initialize FATES logging (the two-stream module's debug/error paths write to this
  ! unit and it is otherwise never set in a standalone driver) and the two-stream
  ! radiation parameters (leaf/stem optical properties, per pft)
  numpft = size(prt_params%wood_density, dim=1)
  call FatesGlobalsInit(6, .false.)
  call TransferRadParams()

  ! derive the organ_id -> parameter-file-index reverse lookup map
  ! (prt_params%organ_param_id) - normally done by the host model's own interface
  ! setup (FatesInterfaceMod.F90), which this standalone driver bypasses entirely
  call PRTDerivedParams()

  ! host-model-namelist-controlled leaf biophysics switches
  lb_params%electron_transport_model = FvCB1980 ! Farquhar-von Caemmerer-Berry (1980)
  lb_params%stomatal_model           = medlyn_model
  lb_params%stomatal_assim_model     = net_assim_model
  lb_params%photo_tempsens_model     = photosynth_acclim_model_kumarathunge_etal_2019

  ! host-model-namelist-controlled carbon-starvation mortality model - matches
  ! CTSM's fates_cstarvation_model default ('linear', bld/namelist_files/
  ! namelist_defaults_ctsm.xml)
  hlm_mort_cstarvation_model = cstarvation_model_lin

  ! prescribed atmospheric and soil boundary conditions - held fixed for the entire
  ! run (see FatesTestEnvironmentMod)
  call env%Init()

  ! canopy gas parameters and leaf N content - constant for the whole run
  call GetCanopyGasParameters(env%can_press, env%can_o2_ppress, env%tempk,       &
    mm_kco2, mm_ko2, co2_cpoint)
  lnc_top = prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(leaf_organ)) / &
    prt_params%slatop(pft)

  ! recruitment-size initialization: start every light level's cohort at the diameter
  ! implied by this PFT's minimum (sapling) recruitment height
  call h2d_allom(EDPftvarcon_inst%hgt_min(pft), pft, dbh_recruit)

  ! build the log-spaced incident light fractions to sweep
  allocate(light_frac(n_light_levels))
  do ilight = 1, n_light_levels
    light_frac(ilight) = light_frac_min * (light_frac_max/light_frac_min) **      &
      (real(ilight - 1, r8)/real(n_light_levels - 1, r8))
  end do

  n_photo_calls_total = 0
  n_bisection_calls_total = 0
  max_solve_iter_total = 0

  ! main light-level sweep: each level is an independent trajectory from recruitment size
  do ilight = 1, n_light_levels
    call RunOneLightLevel(light_frac(ilight), ilight)
  end do

  write(*,*) 'TOTAL photosynthesis solves = ', n_photo_calls_total,                 &
    ', bisection fallback = ', n_bisection_calls_total, ', max solve_iter = ',      &
    max_solve_iter_total

  ! write out the daily whole-cohort time series and the annual light-profile
  ! snapshot, both across the light sweep
  call hist%Write(out_file, light_frac)

contains

  ! ==========================================================================

  subroutine RunOneLightLevel(light_frac_val, ilight)
    !
    ! DESCRIPTION:
    ! Run one light level's independent year/day/sub-daily trajectory, from a
    ! freshly built recruitment-size cohort, recording daily output and the annual
    ! light profile into the driver's shared history object (hist).

    ! ARGUMENTS:
    real(r8), intent(in) :: light_frac_val ! this light level's incident light fraction [0-1]
    integer,  intent(in) :: ilight         ! light-level index

    ! LOCALS:
    type(fates_cohort_type), pointer :: cohort        ! cohort for this light level
    type(light_env_type)             :: light_env     ! prescribed light environment for this light level
    type(cohort_phys_type)           :: phys          ! per-leaf-layer photosynthetic capacity/leaf physics
    integer  :: iyear, iday   ! year/day looping indices
    integer  :: iday_all      ! day index within this light level's trajectory (1..n_days_total)
    integer  :: isubday       ! sub-daily looping index
    real(r8) :: hour_of_day   ! hour of day at the midpoint of the current substep [0-24]
    real(r8) :: frac_store    ! ratio of storage carbon to target_leaf_c [-], from today's DailySetup - used both for output and (pre-growth) by DailyGrowthAndMortality's mortality term
    real(r8) :: npp_acc_to_prt ! carbon_balance handed to PARTEH via cohort%npp_acc, net of growth respiration
    real(r8) :: cmort          ! carbon starvation mortality rate [indiv/year]
    real(r8) :: gpp_tstep, rdark_tstep, nonleaf_mr_tstep ! this substep's GPP/leaf dark resp/non-leaf MR [kgC/indiv/s]
    real(r8) :: daily_net_c    ! GPP - leaf dark resp - non-leaf MR, integrated over the day [kgC/indiv/day] (reset each day, printed at year end, and fed to PARTEH via DailyGrowthAndMortality)
    real(r8) :: daily_gpp      ! GPP alone, integrated over the day [kgC/indiv/day] (reset each day)
    real(r8) :: storage_c      ! current storage carbon [kgC/indiv], recomputed at year end for the print statement below
    integer  :: n_photo_calls, n_bisection_calls, max_solve_iter ! this light level's Ci-solver diagnostics (see the header comment above)

    n_photo_calls = 0
    n_bisection_calls = 0
    max_solve_iter = 0

    ! build a fresh cohort at recruitment size for this light level - CohortFactory
    ! already sets treelai/treesai/c_area correctly for the recruitment-size cohort
    ! (via its own internal tree_lai_sai/carea_allom calls), so no extra refresh is
    ! needed before day 1
    allocate(cohort)
    call CohortFactory(cohort, pft, can_tlai, dbh=dbh_recruit, number=n_indiv,      &
      patch_area=patch_area)
    cohort%nv = GetNVegLayers(cohort%treelai + cohort%treesai)

    call light_env%Init(cohort, pft)

    ! allocate the history arrays once, from the first light level
    if (ilight == 1) then
      call hist%Init(n_days_total, nlevleaf, n_light_levels, nyears)
    end if

    write(*,*) 'Light level ', ilight, ' (fraction of full sun = ', light_frac_val, '):'
    write(*,*) '  initial dbh = ', cohort%dbh, ' cm, height = ', cohort%height, ' m'

    do iyear = 1, nyears
      do iday = 1, days_per_year

        iday_all = (iyear - 1)*days_per_year + iday

        daily_net_c = 0.0_r8
        daily_gpp = 0.0_r8

        ! once-per-day setup: MR throttle, sapwood/fine-root N, and the per-layer
        ! capacity/dark-respiration profiles (see FatesTestCohortPhysMod)
        call phys%DailySetup(cohort, pft, env, lnc_top, frac_store)

        do isubday = 1, n_substeps_per_day

          ! prescribed incident PAR at the top of the crown, attenuated through the
          ! cohort's own leaf layers via the two-stream solver, to get
          ! parsun_z/parsha_z and laisun_z/laisha_z per leaf layer
          hour_of_day = real(isubday, r8) - 0.5_r8
          call light_env%Profile(light_frac_val, hour_of_day)

          ! per-year light-profile snapshot: first day of each year, substep
          ! nearest solar noon - see the header comment above for why this isn't
          ! captured daily. cohort%nv <= nlevleaf always (nlevleaf is a
          ! compile-time maximum), so layers nv+1:nlevleaf are simply never
          ! touched here and keep the fill value pre-set by hist%Init
          if (iday == 1 .and. isubday == noon_substep) then
            call hist%RecordLightProfile(iyear, ilight, cohort, light_env)
          end if

          ! a single sub-daily carbon uptake step: leaf photosynthesis -> GPP and
          ! leaf dark respiration, plus non-leaf maintenance respiration (see
          ! FatesTestCohortPhysMod)
          call phys%SubdailyStep(cohort, pft, env, light_env, mm_kco2, mm_ko2,     &
            co2_cpoint, step_size, n_photo_calls, n_bisection_calls,               &
            max_solve_iter, gpp_tstep, rdark_tstep, nonleaf_mr_tstep)

          ! daily net carbon = GPP - leaf dark resp - non-leaf MR, integrated over
          ! the sub-daily steps [kgC/indiv/day]. Net of growth respiration below,
          ! this is what feeds PARTEH via cohort%npp_acc.
          daily_net_c = daily_net_c + (gpp_tstep - rdark_tstep - nonleaf_mr_tstep) * step_size
          daily_gpp = daily_gpp + gpp_tstep * step_size

        end do

        ! the daily growth sequence (NSC/PRT allocation) and carbon starvation
        ! mortality - kept as a single flat subroutine; see its header comment
        call DailyGrowthAndMortality(cohort, daily_net_c, frac_store, npp_acc_to_prt, cmort)

        ! refresh crown area, treelai/treesai/nv, and the light environment's
        ! canopy structure to reflect today's growth - ready for tomorrow's
        ! photosynthesis, and needed here so today's output row (below) reflects
        ! today's growth consistently. carea_allom is otherwise only ever called
        ! once, at cohort creation (FatesFactoryMod's CohortFactory) - production
        ! only recomputes it as part of patch-level canopy structure
        ! (EDCanopyStructureMod), which this patch-less driver has no equivalent
        ! of. Without this call, cohort%c_area - and the leaf-area scaling of
        ! gpp_tstep that depends on it - would silently stay frozen at the
        ! recruitment-size crown area for the entire multi-year run
        call carea_allom(cohort%dbh, cohort%n, init_spread_inventory, pft,         &
          cohort%crowndamage, cohort%c_area)
        call UpdateCohortLAI(cohort, can_tlai, patch_area)
        call light_env%Refresh(cohort)

        ! capture today's daily time series row
        call hist%RecordDay(iday_all, ilight, cohort, daily_net_c, daily_gpp,      &
          npp_acc_to_prt, frac_store, cmort)

      end do

      storage_c = cohort%prt%GetState(store_organ, carbon12_element)

      write(*,*) '  year ', iyear, ': daily net C (GPP - leaf Rd - nonleaf MR) = ', &
        daily_net_c, ' kgC/indiv/day, dbh = ', cohort%dbh, ' cm, storage C = ',     &
        storage_c, ' kgC'

    end do

    write(*,*) '  photosynthesis solves = ', n_photo_calls, ', bisection fallback = ', &
      n_bisection_calls, ', max solve_iter = ', max_solve_iter

    n_photo_calls_total = n_photo_calls_total + n_photo_calls
    n_bisection_calls_total = n_bisection_calls_total + n_bisection_calls
    max_solve_iter_total = max(max_solve_iter_total, max_solve_iter)

    ! tear down the cohort and its light environment before moving to the next
    ! light level so each level starts over independently from recruitment size
    call light_env%Free()
    call cohort%FreeMemory()
    deallocate(cohort)

  end subroutine RunOneLightLevel

  ! ==========================================================================

  subroutine DailyGrowthAndMortality(cohort, daily_net_c, frac_store, npp_acc_to_prt, cmort)
    !
    ! DESCRIPTION:
    ! The daily growth sequence (NSC update / PRT allocation, once per day after
    ! the sub-daily loop, matching EDMainMod.F90's daily dynamics sequence) followed
    ! by carbon starvation mortality (lifted from EDMortalityFunctionsMod.F90's
    ! mortality_rates). Deliberately kept as a single flat subroutine, not
    ! decomposed further or hidden behind a type: the call order below is the thing
    ! this driver exists to make explicit.

    ! ARGUMENTS:
    type(fates_cohort_type), intent(inout) :: cohort         ! cohort to grow and apply mortality to
    real(r8),                 intent(in)    :: daily_net_c    ! today's net carbon (GPP - leaf Rd - nonleaf MR) [kgC/indiv/day]
    real(r8),                 intent(in)    :: frac_store     ! storage carbon / target leaf carbon, computed pre-growth at the top of today's loop [-]
    real(r8),                 intent(out)   :: npp_acc_to_prt ! carbon balance handed to PARTEH (net of growth respiration) [kgC/indiv/day]
    real(r8),                 intent(out)   :: cmort          ! carbon starvation mortality rate [indiv/year]

    ! LOCALS:
    real(r8) :: resp_g_acc              ! growth respiration: a tax on the positive part of daily_net_c [kgC/indiv/day], matches EDMainMod.F90's resp_g_acc_hold (without that field's day/year unit round-trip, since it cancels out)
    real(r8) :: delta_dbh, delta_height ! unused outputs of EvaluateAndCorrectDBH (it corrects cohort%dbh/height in place; these report how much it moved them)

    ! ---------------------------------------------------------------------
    ! NSC update / PRT allocation - once per day, after the sub-daily loop,
    ! matching EDMainMod.F90's daily dynamics sequence
    ! ---------------------------------------------------------------------

    ! growth respiration: a tax on the positive part of today's net carbon
    ! (matches EDMainMod.F90's resp_g_acc_hold, which additionally round-trips
    ! through an annual-equivalent "hold" value for I/O - that cancels out
    ! algebraically and is skipped here since this driver has no such
    ! diagnostic)
    resp_g_acc = prt_params%grperc(pft) * max(0.0_r8, daily_net_c)

    ! the actual carbon_balance boundary condition PARTEH reads (see
    ! FatesCohortMod.F90's InitPRTBoundaryConditions, which registers
    ! cohort%npp_acc as ac_bc_inout_id_netdc) - net of growth respiration.
    ! Captured into npp_acc_to_prt for output before DailyPRT decrements
    ! cohort%npp_acc in place as it allocates
    cohort%npp_acc = daily_net_c - resp_g_acc
    npp_acc_to_prt = cohort%npp_acc

    ! maintenance turnover moves senesced tissue mass to storage/litter
    ! before allocation replaces it (phase 1, below). is_drought = .false.
    ! always: this driver's evergreen PFT with non-limiting soil moisture
    ! (btran=1) never triggers phenological drought status
    call PRTMaintTurnover(cohort%prt, pft, cohort%canopy_layer, .false.)
    call cohort%prt%AgeLeaves(pft, cohort%canopy_layer, sec_per_day)

    ! correct dbh if structural carbon has drifted from what's allometrically
    ! consistent (integration error correction) - this driver never damages
    ! the crown, so EDMainMod.F90's "newly_recovered" skip never applies here
    call EvaluateAndCorrectDBH(cohort, delta_dbh, delta_height)

    ! phase 1: replace turnover losses (from storage if npp_acc < 0); phase 2:
    ! push toward allometric targets; phase 3: grow stature with whatever
    ! carbon is left over. cohort%npp_acc is decremented in place by PARTEH as
    ! it allocates (it is a bc_inout, not bc_in, boundary condition)
    call cohort%prt%DailyPRT(phase=1)
    call cohort%prt%DailyPRT(phase=2)
    call cohort%prt%DailyPRT(phase=3)

    ! leaf-age-class mass fractions changed - refresh the biophysical rates
    ! (vcmax25top/jmax25top/tpu25top/kp25top) that tomorrow's photosynthesis
    ! loop will use
    call cohort%UpdateCohortBioPhysRates()

    ! dbh may have grown in phase 3 - update height to match
    call h_allom(cohort%dbh, pft, cohort%height)

    ! carbon has been spent - zero the accumulator for tomorrow, matching
    ! EDMainMod.F90's end-of-dynamics reset
    cohort%npp_acc = 0.0_r8

    ! ---------------------------------------------------------------------
    ! carbon starvation mortality - lifted from EDMortalityFunctionsMod.F90's
    ! mortality_rates (the cmort block only; background/hydraulic/freezing/
    ! senescence/damage mortality are not implemented here). Deliberately
    ! uses frac_store as computed at the top of today's loop (pre-growth),
    ! matching production's Mortality_Derivative, which runs before
    ! PRTMaintTurnover/DailyPRT using the previous day's storage state
    ! ---------------------------------------------------------------------
    select case (hlm_mort_cstarvation_model)
    case (cstarvation_model_lin)
      ! zero once frac_store reaches mort_upthresh_cstarvation (storage at or
      ! above target leaf carbon), rising linearly to mort_scalar_cstarvation
      ! (the max rate, /year) as frac_store -> 0
      cmort = EDPftvarcon_inst%mort_scalar_cstarvation(pft) * max(0.0_r8,       &
        (EDPftvarcon_inst%mort_upthresh_cstarvation(pft) - frac_store) /        &
        EDPftvarcon_inst%mort_upthresh_cstarvation(pft))
    case (cstarvation_model_exp)
      ! mort_scalar_cstarvation (the max rate, /year) at frac_store=0, decaying
      ! exponentially as frac_store rises; mort_upthresh_cstarvation sets the
      ! e-folding scale
      cmort = EDPftvarcon_inst%mort_scalar_cstarvation(pft) *                   &
        exp(-frac_store / EDPftvarcon_inst%mort_upthresh_cstarvation(pft))
    end select
    if (cmort <= nearzero) cmort = 0.0_r8

    ! Euler step on cohort number density: dn/dt = -cmort*n [indiv/year],
    ! stepped by 1/days_per_year - matches EDMainMod.F90's
    ! cohort%n = max(0, n + dndt*hlm_freq_day). No disturbance mechanism
    ! exists in this patch-less driver, so the full rate is applied directly
    ! to n (see the header comment above)
    cohort%n = max(0.0_r8, cohort%n * (1.0_r8 - cmort/real(days_per_year, r8)))

  end subroutine DailyGrowthAndMortality

end program FatesSingleCohort
