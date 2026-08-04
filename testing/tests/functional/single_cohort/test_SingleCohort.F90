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
  !     conditions (pressure, CO2/O2, boundary-layer conductance, soil
  !     moisture/rooting - held fixed for the entire run; temperature follows a
  !     real annual/diurnal cycle fit to BCI DATM data, with genuine running-mean
  !     T_growth/T_home acclimation temperatures).
  !   - FatesTestLightEnvMod - the prescribed light environment (reference full-sun
  !     PAR, direct/diffuse split, and a real annual/diurnal solar cycle - solar
  !     declination from day of year and coszen(hour) from a prescribed site
  !     latitude and declination), attenuated through the cohort's own leaf
  !     layers via FATES's two-stream radiation solver.
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
  ! Prescribed cold-deciduous phenology (a fixed leaf-on/leaf-off day of year, see
  ! DailyPhenology below) runs once per day, before the growth sequence, matching
  ! EDMainMod.F90's relative ordering (phenology, then PRTMaintTurnover/DailyPRT).
  ! It reuses production's actual flush/tissue-drop carbon mechanics
  ! (PRTPhenologyFlush/PRTDeciduousTurnover), substituting a direct day-of-year
  ! comparison for the site-level GDD/chilling-day tracking that would otherwise
  ! decide the leaf-on/off transition (EDPhysiologyMod.F90's phenology_leafonoff)
  ! - this driver's default PFT is evergreen, so phenology is inert unless pointed
  ! at an ihard_season_decid PFT.
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
  ! OUTPUT: three groups of variables are accumulated in memory (FatesTestHistoryMod)
  ! and written once at the end:
  !   - Daily whole-cohort time series, dimensioned (time, light_level), where time is
  !     the day index within each light level's independent trajectory (1 to
  !     nyears*days_per_year): dbh, height, treelai, treesai, nv, crown area, the
  !     PARTEH carbon pools (leaf/fine root/sapwood/structure/storage/reproduction),
  !     frac_store (storage as a fraction of target leaf carbon), cohort number
  !     density (n), and the daily flux terms - kept disaggregated, not netted,
  !     since decomposing the compensation point into leaf/architectural/support-
  !     respiration contributions is impossible from a netted total: daily GPP,
  !     leaf dark respiration, live stem/live coarse root/fine root maintenance
  !     respiration, growth respiration, per-organ turnover loss (leaf/fine root/
  !     sapwood/structure), npp_acc as handed to PARTEH (net of growth
  !     respiration), and the carbon starvation mortality rate (cmort,
  !     indiv/year). daily_net_c is also kept as a convenience total.
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
  !   - An instantaneous whole-plant light-response diagnostic (gross_assim/
  !     total_resp, via LightResponseSweep), dimensioned (ppfd, year, light_level),
  !     captured at the same once-per-year cadence as the light profile above -
  !     see LightResponseSweep's header comment for the method.
  !

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesConstantsMod,           only : nearzero
  use FatesConstantsMod,           only : sec_per_day
  use FatesConstantsMod,           only : cstarvation_model_lin, cstarvation_model_exp
  use FatesConstantsMod,           only : itrue, ifalse
  use FatesConstantsMod,           only : wm2_to_umolm2s
  use FatesConstantsMod,           only : leaves_on, leaves_off
  use FatesConstantsMod,           only : ihard_season_decid
  use EDParamsMod,                 only : nclmax
  use EDParamsMod,                 only : nlevleaf
  use EDParamsMod,                 only : GetNVegLayers
  use EDPftvarcon,                 only : EDPftvarcon_inst
  use FatesAllometryMod,           only : h2d_allom
  use FatesAllometryMod,           only : h_allom
  use FatesAllometryMod,           only : carea_allom
  use FatesAllometryMod,           only : bleaf
  use FatesAllometryMod,           only : bfineroot
  use FatesAllometryMod,           only : bsap_allom
  use FatesAllometryMod,           only : bagw_allom
  use FatesAllometryMod,           only : bbgw_allom
  use FatesAllometryMod,           only : bdead_allom
  use FatesCohortMod,              only : fates_cohort_type
  use EDCanopyStructureMod,        only : UpdateCohortLAI
  use EDCohortDynamicsMod,         only : EvaluateAndCorrectDBH
  use PRTLossFluxesMod,            only : PRTMaintTurnover
  use PRTLossFluxesMod,            only : PRTPhenologyFlush
  use PRTLossFluxesMod,            only : PRTDeciduousTurnover
  use FatesUnitTestParamReaderMod, only : ReadParameters
  use FatesArgumentUtils,          only : command_line_arg
  use FatesFactoryMod,             only : InitializeGlobals, CohortFactory
  use FatesGlobals,                only : FatesGlobalsInit
  use FatesGlobals,                only : fates_log
  use FatesInterfaceTypesMod,      only : numpft, hlm_mort_cstarvation_model
  use PRTParametersMod,            only : prt_params
  use PRTGenericMod,               only : leaf_organ, fnrt_organ, sapw_organ, struct_organ
  use PRTGenericMod,               only : store_organ, carbon12_element
  use PRTInitParamsFatesMod,       only : PRTDerivedParams
  use FatesTwoStreamUtilsMod,      only : TransferRadParams
  use FatesTestEnvironmentMod,     only : environment_type
  use FatesTestLightEnvMod,        only : light_env_type
  use FatesTestCohortPhysMod,      only : cohort_phys_type
  use FatesTestHistoryMod,         only : history_type
  use LeafBiophysicsMod,           only : lb_params
  use LeafBiophysicsMod,           only : FvCB1980, medlyn_model
  use LeafBiophysicsMod,           only : net_assim_model, photosynth_acclim_model_kumarathunge_etal_2019
  use LeafBiophysicsMod,           only : LowstorageMainRespReduction
  use EDParamsMod               , only : dinc_vai

  implicit none

  ! LOCALS:
  character(len=:), allocatable    :: param_file    ! input parameter file
  type(environment_type)           :: env           ! prescribed atmospheric/soil boundary conditions - re-Init'd fresh per light level (see RunOneLightLevel): the annual/diurnal temperature cycle is a deterministic function of day/hour, but its T_growth/T_home running-mean bookkeeping needs to restart at each light level's day 1
  type(history_type)               :: hist          ! output accumulator/writer, shared across the light sweep
  real(r8)                         :: dbh_recruit   ! recruitment-size dbh [cm]
  real(r8), allocatable            :: light_frac(:) ! swept incident light fractions [0-1]
  real(r8), allocatable            :: diagnostic_ppfd(:) ! swept PPFD values for LightResponseSweep [umol/m2/s]
  integer                          :: ilight        ! light-level looping index

  ! leaf N content - constant for the whole run since the atmospheric boundary
  ! conditions never change
  real(r8) :: lnc_top ! leaf N content at the canopy top [gN/m2 leaf]

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
  real(r8), parameter :: coh_age = 0.0_r8                    ! cohort age
  real(r8), parameter :: site_spread = 1.0_r8

  ! phenology - prescribed leaf-on/leaf-off day of year (cold-deciduous PFTs only;
  ! see DailyPhenology). Illustrative Northern-Hemisphere mid-latitude timing, not
  ! fit to any real site - this driver's PFT 1 default is evergreen, so these are
  ! inert unless pft above is pointed at an ihard_season_decid PFT
  integer, parameter :: leaf_on_doy  = 60  ! day of year leaves flush [1-365]
  integer, parameter :: leaf_off_doy = 305 ! day of year leaves shed [1-365]

  ! light sweep
  integer,  parameter :: n_light_levels = 25       ! number of incident light levels to sweep
  real(r8), parameter :: light_frac_min = 0.005_r8 ! lowest incident light fraction swept [fraction of full sun]
  real(r8), parameter :: light_frac_max = 1.0_r8   ! highest incident light fraction swept [fraction of full sun]

  ! instantaneous light-response diagnostic (LightResponseSweep) - one snapshot
  ! per year (first day, solar noon), log-spaced PPFD from ppfd_diagnostic_min
  ! to ppfd_diagnostic_max (~full sun). ppfd_diagnostic_min is set well below
  ! the whole-plant compensation point at the largest simulated cohort (the
  ! highest light level's final day): empirically, that compensation point
  ! (pure-diffuse illumination, coszen=1) falls between 31 and 47 umol/m2/s -
  ! a fixed range chosen from early, small-cohort behavior would be far too
  ! narrow by the run's end, since architectural respiration (and so the
  ! compensation point) grows substantially as the plant grows
  logical,  parameter :: enable_light_response_diagnostic = .true. ! set false only to regression-test that the diagnostic perturbs no state (see RunOneLightLevel)
  integer,  parameter :: n_ppfd_diagnostic = 40      ! number of PPFD values to sweep
  real(r8), parameter :: ppfd_diagnostic_min = 0.01_r8   ! lowest swept PPFD [umol/m2/s]
  real(r8), parameter :: ppfd_diagnostic_max = 2000.0_r8 ! highest swept PPFD [umol/m2/s] (~full sun)

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
  
  print *, dinc_vai, nlevleaf

  ! ! derive the organ_id -> parameter-file-index reverse lookup map
  ! ! (prt_params%organ_param_id) - normally done by the host model's own interface
  ! ! setup (FatesInterfaceMod.F90), which this standalone driver bypasses entirely
  ! call PRTDerivedParams()

  ! ! host-model-namelist-controlled leaf biophysics switches
  ! lb_params%electron_transport_model = FvCB1980 ! Farquhar-von Caemmerer-Berry (1980)
  ! lb_params%stomatal_model           = medlyn_model
  ! lb_params%stomatal_assim_model     = net_assim_model
  ! lb_params%photo_tempsens_model     = photosynth_acclim_model_kumarathunge_etal_2019

  ! ! host-model-namelist-controlled carbon-starvation mortality model - matches
  ! ! CTSM's fates_cstarvation_model default ('linear', bld/namelist_files/
  ! ! namelist_defaults_ctsm.xml)
  ! hlm_mort_cstarvation_model = cstarvation_model_lin

  ! ! leaf N content - constant for the whole run
  ! lnc_top = prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(leaf_organ)) / &
  !   prt_params%slatop(pft)

  ! ! recruitment-size initialization: start every light level's cohort at the diameter
  ! ! implied by this PFT's minimum (sapling) recruitment height
  ! call h2d_allom(EDPftvarcon_inst%hgt_min(pft), pft, dbh_recruit)

  ! ! build the log-spaced incident light fractions to sweep
  ! allocate(light_frac(n_light_levels))
  ! do ilight = 1, n_light_levels
  !   light_frac(ilight) = light_frac_min * (light_frac_max/light_frac_min) **      &
  !     (real(ilight - 1, r8)/real(n_light_levels - 1, r8))
  ! end do

  ! ! build the log-spaced diagnostic PPFD sweep (see the parameter block above
  ! ! for why ppfd_diagnostic_min is set where it is)
  ! block
  !   integer :: ippfd
  !   allocate(diagnostic_ppfd(n_ppfd_diagnostic))
  !   do ippfd = 1, n_ppfd_diagnostic
  !     diagnostic_ppfd(ippfd) = ppfd_diagnostic_min *                             &
  !       (ppfd_diagnostic_max/ppfd_diagnostic_min) **                             &
  !       (real(ippfd - 1, r8)/real(n_ppfd_diagnostic - 1, r8))
  !   end do
  ! end block

  ! n_photo_calls_total = 0
  ! n_bisection_calls_total = 0
  ! max_solve_iter_total = 0

  ! ! main light-level sweep: each level is an independent trajectory from recruitment size
  ! do ilight = 1, n_light_levels
  !   call RunOneLightLevel(light_frac(ilight), ilight)
  ! end do

  ! ! write out the daily whole-cohort time series and the annual light-profile
  ! ! snapshot, both across the light sweep
  ! call hist%Write(out_file, light_frac, diagnostic_ppfd)

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
    real(r8) :: daily_rdark        ! leaf dark respiration, integrated over the day [kgC/indiv/day] (reset each day)
    real(r8) :: daily_livestem_mr  ! live stem MR, integrated over the day [kgC/indiv/day] (reset each day)
    real(r8) :: daily_livecroot_mr ! live coarse root MR, integrated over the day [kgC/indiv/day] (reset each day)
    real(r8) :: daily_froot_mr     ! fine root MR, integrated over the day [kgC/indiv/day] (reset each day)
    real(r8) :: growth_resp        ! today's growth respiration [kgC/indiv/day], from DailyGrowthAndMortality
    real(r8) :: leaf_turnover, fnrt_turnover, sapw_turnover, struct_turnover ! today's per-organ turnover loss [kgC/indiv/day], from DailyGrowthAndMortality
    real(r8) :: storage_c      ! current storage carbon [kgC/indiv], recomputed at year end for the print statement below
    integer  :: n_photo_calls, n_bisection_calls, max_solve_iter, sum_solve_iter ! this light level's Ci-solver diagnostics (see the header comment above)
    real(r8) :: mean_solve_iter ! this light level's mean Ci-solver iteration count (sum_solve_iter/n_photo_calls), computed once at the end of this trajectory
    real(r8) :: diagnostic_gross_assim(n_ppfd_diagnostic) ! this year's LightResponseSweep gross-assimilation output [kgC/indiv/s]
    real(r8) :: diagnostic_total_resp(n_ppfd_diagnostic)  ! this year's LightResponseSweep total-respiration output [kgC/indiv/s]
    real(r8) :: par_toc               ! this substep's incident PAR at the top of the crown [W/m2], from light_env%Profile
    real(r8) :: daily_absorbed_par    ! whole-plant absorbed PAR (parsun_z+parsha_z summed over layers), integrated over the day [J/m2 ground] (reset each day)
    real(r8) :: daily_incident_par    ! incident PAR at the top of the crown, integrated over the day [J/m2 ground] (reset each day)
    real(r8) :: light_intercept_eff   ! today's light interception efficiency: daily_absorbed_par/daily_incident_par [-] (Sterck et al. 2013) - crown_area cancels out of both the numerator and denominator of "whole-plant absorbed PAR / (incident PAR at top of crown x crown area)", since parsun_z/parsha_z/par_toc are all already per-unit-ground-area and crown_area would multiply through both sides identically
    real(r8) :: maintresp_reduction_factor ! today's storage-based maintenance-respiration throttle [0-1], recomputed here via the same LowstorageMainRespReduction(frac_store, pft, ...) call phys%DailySetup already made internally (not read off phys, which keeps it private) - deterministic given frac_store, so bit-identical

    n_photo_calls = 0
    n_bisection_calls = 0
    max_solve_iter = 0
    sum_solve_iter = 0

    ! build a fresh cohort at recruitment size for this light level - CohortFactory
    ! already sets treelai/treesai/c_area correctly for the recruitment-size cohort
    ! (via its own internal tree_lai_sai/carea_allom calls), so no extra refresh is
    ! needed before day 1
    allocate(cohort)
    call CohortFactory(cohort, pft, can_tlai, dbh=dbh_recruit, number=n_indiv,      &
      patch_area=patch_area, age=coh_age, site_spread=site_spread)
    cohort%nv = GetNVegLayers(cohort%treelai + cohort%treesai)

    call light_env%Init(cohort, pft)

    ! prescribed atmospheric/soil boundary conditions for this light level's
    ! trajectory (see FatesTestEnvironmentMod) - re-Init'd fresh per light level so
    ! the T_growth/T_home running-mean bookkeeping restarts at day 1 rather than
    ! carrying over the previous light level's ending state
    call env%Init()

    ! allocate the history arrays once, from the first light level
    if (ilight == 1) then
      call hist%Init(n_days_total, nlevleaf, n_light_levels, nyears, n_ppfd_diagnostic)
    end if

    do iyear = 1, nyears
      do iday = 1, days_per_year

        iday_all = (iyear - 1)*days_per_year + iday

        daily_net_c = 0.0_r8
        daily_gpp = 0.0_r8
        daily_rdark = 0.0_r8
        daily_livestem_mr = 0.0_r8
        daily_livecroot_mr = 0.0_r8
        daily_froot_mr = 0.0_r8
        daily_absorbed_par = 0.0_r8
        daily_incident_par = 0.0_r8

        ! once-per-day setup: MR throttle, sapwood/fine-root N, and the per-layer
        ! nitrogen-scaling factor (see FatesTestCohortPhysMod)
        call phys%DailySetup(cohort, pft, frac_store)

        ! today's storage-based maintenance-respiration throttle, for output -
        ! recomputes what phys%DailySetup already derived internally (same
        ! production call, same frac_store input, so bit-identical)
        call LowstorageMainRespReduction(frac_store, pft, maintresp_reduction_factor)

        do isubday = 1, n_substeps_per_day

          ! prescribed incident PAR at the top of the crown, attenuated through the
          ! cohort's own leaf layers via the two-stream solver, to get
          ! parsun_z/parsha_z and laisun_z/laisha_z per leaf layer
          hour_of_day = real(isubday, r8) - 0.5_r8
          call light_env%Profile(light_frac_val, iday, hour_of_day,              &
            par_toc_out=par_toc)

          ! prescribed temperature (and everything derived from it) for this
          ! substep, from the annual/diurnal cycle fit to BCI data (see
          ! FatesTestEnvironmentMod); accumulates into today's running sum for
          ! UpdateDailyMeans to consume once the sub-daily loop finishes
          call env%SetHour(iday, hour_of_day)

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
          call phys%SubdailyStep(cohort, pft, env, light_env, lnc_top, step_size,  &
            n_photo_calls, n_bisection_calls, max_solve_iter, sum_solve_iter,      &
            gpp_tstep, rdark_tstep, nonleaf_mr_tstep)

          ! light interception efficiency (Sterck et al. 2013): whole-plant
          ! absorbed PAR vs. incident PAR at the top of the crown, integrated
          ! over the day and taken as a ratio once the day is done (below) -
          ! avoids a 0/0 division at night, when both sums are simply 0
          daily_absorbed_par = daily_absorbed_par +                              &
            sum(light_env%parsun_z(:) + light_env%parsha_z(:)) * step_size
          daily_incident_par = daily_incident_par + par_toc * step_size

          ! instantaneous whole-plant light-response diagnostic: one snapshot per
          ! year (first day, solar noon - matching RecordLightProfile's existing
          ! annual cadence), gated by enable_light_response_diagnostic so its
          ! on/off state can be regression-tested for exact bit-identical daily
          ! output (see the header comment above)
          if (enable_light_response_diagnostic .and. iday == 1 .and.             &
            isubday == noon_substep) then
            call LightResponseSweep(cohort, pft, light_env, phys, env,           &
              diagnostic_ppfd, light_frac_val, iday, hour_of_day,                &
              diagnostic_gross_assim, diagnostic_total_resp)
            call hist%RecordLightResponse(iyear, ilight, diagnostic_gross_assim, &
              diagnostic_total_resp)
          end if

          ! daily net carbon = GPP - leaf dark resp - non-leaf MR, integrated over
          ! the sub-daily steps [kgC/indiv/day]. Net of growth respiration below,
          ! this is what feeds PARTEH via cohort%npp_acc. The four terms making up
          ! nonleaf_mr_tstep/rdark_tstep are also integrated separately below (not
          ! netted) - cohort%livestem_mr/livecroot_mr/froot_mr were just set as a
          ! side effect of SubdailyStep, for this same substep
          daily_net_c = daily_net_c + (gpp_tstep - rdark_tstep - nonleaf_mr_tstep) * step_size
          daily_gpp = daily_gpp + gpp_tstep * step_size
          daily_rdark = daily_rdark + rdark_tstep * step_size
          daily_livestem_mr = daily_livestem_mr + cohort%livestem_mr * step_size
          daily_livecroot_mr = daily_livecroot_mr + cohort%livecroot_mr * step_size
          daily_froot_mr = daily_froot_mr + cohort%froot_mr * step_size

        end do

        ! today's light interception efficiency - see the header comment above
        ! the sub-daily loop. daily_incident_par is 0 only if there is no
        ! daylight at all today, which never happens at this driver's latitude
        if (daily_incident_par > nearzero) then
          light_intercept_eff = daily_absorbed_par / daily_incident_par
        else
          light_intercept_eff = 0.0_r8
        end if

        ! roll today's sub-daily tempk samples into the T_growth/T_home running
        ! means (see FatesTestEnvironmentMod) - ready for tomorrow's photosynthetic
        ! acclimation
        call env%UpdateDailyMeans()

        ! prescribed leaf-on/leaf-off phenology (cold-deciduous PFTs only) - runs
        ! before the growth sequence below, matching EDMainMod.F90's relative
        ! ordering (phenology, then PRTMaintTurnover/DailyPRT)
        call DailyPhenology(cohort, pft, iday)

        ! the daily growth sequence (NSC/PRT allocation) and carbon starvation
        ! mortality - kept as a single flat subroutine; see its header comment
        call DailyGrowthAndMortality(cohort, daily_net_c, frac_store,            &
          npp_acc_to_prt, cmort, growth_resp, leaf_turnover, fnrt_turnover,      &
          sapw_turnover, struct_turnover)

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
        call carea_allom(cohort%dbh, cohort%n, site_spread, pft,         &
          cohort%crowndamage, cohort%c_area)
        call UpdateCohortLAI(cohort, can_tlai, patch_area)
        call light_env%Refresh(cohort)

        ! capture today's daily time series row
        call hist%RecordDay(iday_all, ilight, cohort, daily_net_c, daily_gpp,      &
          daily_rdark, daily_livestem_mr, daily_livecroot_mr, daily_froot_mr,      &
          growth_resp, leaf_turnover, fnrt_turnover, sapw_turnover,                &
          struct_turnover, npp_acc_to_prt, frac_store, cmort,                      &
          light_intercept_eff, maintresp_reduction_factor)

      end do

      storage_c = cohort%prt%GetState(store_organ, carbon12_element)

    end do

    n_photo_calls_total = n_photo_calls_total + n_photo_calls
    n_bisection_calls_total = n_bisection_calls_total + n_bisection_calls
    max_solve_iter_total = max(max_solve_iter_total, max_solve_iter)

    ! this light level's whole-trajectory Ci-solver summary - mean iteration
    ! count and bisection-fallback count, over every LeafLayerPhotosynthesis
    ! call across all n_days_total days (not a daily quantity: with ~2 calls
    ! per leaf layer per substep, a per-day mean would be noisy and not
    ! particularly informative on its own)
    mean_solve_iter = real(sum_solve_iter, r8) / real(n_photo_calls, r8)
    call hist%RecordLightLevelSummary(ilight, mean_solve_iter, n_bisection_calls)

    ! tear down the cohort and its light environment before moving to the next
    ! light level so each level starts over independently from recruitment size
    call light_env%Free()
    call cohort%FreeMemory()
    deallocate(cohort)

  end subroutine RunOneLightLevel

  ! ==========================================================================

  subroutine DailyPhenology(cohort, pft, day_of_year)
    !
    ! DESCRIPTION:
    ! Prescribed cold-deciduous (ihard_season_decid) leaf-on/leaf-off phenology,
    ! reusing production's actual carbon mechanics - allometric targets (bleaf/
    ! bfineroot/bsap_allom/bagw_allom/bbgw_allom/bdead_allom) and the storage
    ! flush/tissue drop physics (PRTPhenologyFlush/PRTDeciduousTurnover) - lifted
    ! from EDPhysiologyMod.F90's phenology_leafonoff (its per-cohort body, which
    ! needs no patch/site state beyond the site-level elongation factor and cold
    ! status this driver doesn't have). Those two site-level inputs are replaced
    ! here by a direct day-of-year comparison against prescribed leaf_on_doy/
    ! leaf_off_doy, which is what buys independence from that routine's GDD/
    ! chilling-day tracking (currentSite%cstatus in the production version).
    ! Evergreen and drought-deciduous PFTs (ihard_stress_decid/isemi_stress_decid
    ! - meaningless here anyway with btran permanently non-limiting) are
    ! untouched, matching their existing always-fully-flushed behavior.

    ! ARGUMENTS:
    type(fates_cohort_type), intent(inout) :: cohort      ! cohort to update
    integer,                  intent(in)    :: pft         ! plant functional type index
    integer,                  intent(in)    :: day_of_year ! day of year [1-365]

    ! LOCALS:
    logical  :: is_leaf_on_season ! true on days within [leaf_on_doy, leaf_off_doy)
    logical  :: is_flushing_time  ! true only on the day leaves transition off -> on
    logical  :: is_shedding_time  ! true only on the day leaves transition on -> off
    real(r8) :: elong_factor_today ! today's leaf elongation factor [0 or 1]
    real(r8) :: leaf_c, fnrt_c, sapw_c, struct_c, store_c ! current tissue carbon [kgC]
    real(r8) :: leaf_deficit_c, fnrt_deficit_c, sapw_deficit_c, struct_deficit_c, total_deficit_c ! flush-time carbon deficits relative to target [kgC]
    real(r8) :: eff_leaf_drop_fraction, eff_fnrt_drop_fraction ! shed-time drop fractions [0-1]
    real(r8) :: eff_sapw_drop_fraction, eff_struct_drop_fraction ! shed-time drop fractions [0-1]
    real(r8) :: target_leaf_c, target_fnrt_c, target_sapw_c ! target tissue carbon at today's elongation factor [kgC]
    real(r8) :: target_agw_c, target_bgw_c, target_struct_c ! target tissue carbon at today's elongation factor [kgC]
    real(r8) :: sapw_area          ! unused diagnostic output of bsap_allom [m2]
    real(r8) :: store_c_transfer_frac ! fraction of storage carbon transferred to flush tissues [-]
    real(r8), parameter :: carbon_store_buffer = 0.10_r8 ! matches phenology_leafonoff's identical local constant - leaves this fraction of storage untouched by a flush, to avoid triggering carbon-starvation mortality from a single flush event

    if (prt_params%phen_leaf_habit(pft) /= ihard_season_decid) return

    is_leaf_on_season = (day_of_year >= leaf_on_doy) .and. (day_of_year < leaf_off_doy)
    elong_factor_today = merge(1.0_r8, 0.0_r8, is_leaf_on_season)

    is_flushing_time = is_leaf_on_season .and. (cohort%status_coh == leaves_off)
    is_shedding_time = (.not. is_leaf_on_season) .and. (cohort%status_coh == leaves_on) .and. &
      (cohort%dbh > EDPftvarcon_inst%phen_cold_size_threshold(pft) .or.          &
       prt_params%woody(pft) == itrue)

    ! elongation factors always track today's phenological season, independent of
    ! whether today is a transition day - matches phenology_leafonoff's
    ! unconditional efleaf_coh/effnrt_coh/efstem_coh update
    cohort%efleaf_coh = elong_factor_today
    cohort%effnrt_coh = 1.0_r8 - (1.0_r8 - cohort%efleaf_coh) * prt_params%phen_fnrt_drop_fraction(pft)
    cohort%efstem_coh = 1.0_r8 - (1.0_r8 - cohort%efleaf_coh) * prt_params%phen_stem_drop_fraction(pft)

    if (.not. (is_flushing_time .or. is_shedding_time)) return

    ! target tissue carbon at today's elongation factor - only actually consumed
    ! below, so (unlike phenology_leafonoff) only computed on a transition day
    call bleaf(cohort%dbh, pft, cohort%crowndamage, cohort%canopy_trim,          &
      cohort%efleaf_coh, target_leaf_c)
    call bfineroot(cohort%dbh, pft, cohort%canopy_trim, cohort%l2fr,            &
      cohort%effnrt_coh, target_fnrt_c)
    call bsap_allom(cohort%dbh, pft, cohort%crowndamage, cohort%canopy_trim,     &
      cohort%efstem_coh, sapw_area, target_sapw_c)
    call bagw_allom(cohort%dbh, pft, cohort%crowndamage, cohort%efstem_coh,      &
      target_agw_c)
    call bbgw_allom(cohort%dbh, pft, cohort%efstem_coh, target_bgw_c)
    call bdead_allom(target_agw_c, target_bgw_c, target_sapw_c, pft, target_struct_c)

    if (is_flushing_time) then
      cohort%status_coh = leaves_on

      store_c = cohort%prt%GetState(store_organ, carbon12_element)
      if (store_c > nearzero) then
        leaf_c   = cohort%prt%GetState(leaf_organ,   carbon12_element)
        fnrt_c   = cohort%prt%GetState(fnrt_organ,   carbon12_element)
        sapw_c   = cohort%prt%GetState(sapw_organ,   carbon12_element)
        struct_c = cohort%prt%GetState(struct_organ, carbon12_element)

        leaf_deficit_c   = max(0.0_r8, target_leaf_c   - leaf_c)
        fnrt_deficit_c   = max(0.0_r8, target_fnrt_c   - fnrt_c)
        sapw_deficit_c   = max(0.0_r8, target_sapw_c   - sapw_c)
        struct_deficit_c = max(0.0_r8, target_struct_c - struct_c)
        total_deficit_c  = leaf_deficit_c + fnrt_deficit_c + sapw_deficit_c + struct_deficit_c

        if (total_deficit_c > nearzero) then
          store_c_transfer_frac = min(EDPftvarcon_inst%phenflush_fraction(pft) *  &
            total_deficit_c / store_c, 1.0_r8 - carbon_store_buffer)

          call PRTPhenologyFlush(cohort%prt, pft, leaf_organ,                    &
            store_c_transfer_frac*leaf_deficit_c/total_deficit_c)
          call PRTPhenologyFlush(cohort%prt, pft, fnrt_organ,                    &
            store_c_transfer_frac*fnrt_deficit_c/total_deficit_c)
          if (prt_params%woody(pft) == ifalse) then
            call PRTPhenologyFlush(cohort%prt, pft, sapw_organ,                  &
              store_c_transfer_frac*sapw_deficit_c/total_deficit_c)
            call PRTPhenologyFlush(cohort%prt, pft, struct_organ,                &
              store_c_transfer_frac*struct_deficit_c/total_deficit_c)
          end if
        end if
      end if
    end if

    if (is_shedding_time) then
      cohort%status_coh = leaves_off

      leaf_c   = cohort%prt%GetState(leaf_organ,   carbon12_element)
      fnrt_c   = cohort%prt%GetState(fnrt_organ,   carbon12_element)
      sapw_c   = cohort%prt%GetState(sapw_organ,   carbon12_element)
      struct_c = cohort%prt%GetState(struct_organ, carbon12_element)

      eff_leaf_drop_fraction   = max(0.0_r8, min(1.0_r8, 1.0_r8 - target_leaf_c/max(leaf_c, nearzero)))
      eff_fnrt_drop_fraction   = max(0.0_r8, min(1.0_r8, 1.0_r8 - target_fnrt_c/max(fnrt_c, nearzero)))
      eff_sapw_drop_fraction   = max(0.0_r8, min(1.0_r8, 1.0_r8 - target_sapw_c/max(sapw_c, nearzero)))
      eff_struct_drop_fraction = max(0.0_r8, min(1.0_r8, 1.0_r8 - target_struct_c/max(struct_c, nearzero)))

      call PRTDeciduousTurnover(cohort%prt, pft, leaf_organ, eff_leaf_drop_fraction)
      call PRTDeciduousTurnover(cohort%prt, pft, fnrt_organ, eff_fnrt_drop_fraction)
      if (prt_params%woody(pft) == ifalse) then
        call PRTDeciduousTurnover(cohort%prt, pft, sapw_organ,   eff_sapw_drop_fraction)
        call PRTDeciduousTurnover(cohort%prt, pft, struct_organ, eff_struct_drop_fraction)
      end if
    end if

  end subroutine DailyPhenology

  ! ==========================================================================

  subroutine LightResponseSweep(cohort, pft, light_env, phys, env, ppfd_values, &
    light_frac_val, day_of_year, hour_of_day, gross_assim, total_resp,          &
    use_beam_illumination)
    !
    ! DESCRIPTION:
    ! Instantaneous whole-plant light-response diagnostic: freezes the cohort's
    ! current state (today's already-computed leaf capacity in phys, and
    ! cohort's current dbh/c_area/nv/treelai/treesai) and sweeps a caller-
    ! supplied list of incident PPFD values through the same two-stream
    ! self-shading physics as the real light environment
    ! (light_env%AttenuateCanopy), recording whole-plant gross assimilation and
    ! total respiration at each swept PPFD via phys%GrossAssimAndResp - without
    ! advancing any state. This is a separate calculation from the daily/
    ! sub-daily loop (which only produces daily integrals): unlike Profile/
    ! SubdailyStep, it writes nothing to cohort, phys, or env.
    !
    ! Illumination defaults to pure diffuse (all incident PPFD as diffuse, no
    ! direct beam) - matching Sterck et al. (2013)'s definition of light
    ! interception efficiency, the comparison target for this diagnostic.
    ! use_beam_illumination switches to the opposite extreme (all incident
    ! PPFD as direct beam at zenith) so the two can be compared. Sun-angle
    ! geometry is fixed at coszen=1 (sun directly overhead) regardless of the
    ! illumination mode or the real solar geometry for day_of_year/hour_of_day
    ! - with pure diffuse illumination coszen is nearly inert anyway (it only
    ! enters via the two-stream diffuse-radiation streams' own internal
    ! integration over sky angles, not a cosine projection the way beam
    ! radiation uses it), but fixing it regardless decouples the response
    ! curve's shape from whatever declination/hour a particular diagnostic day
    ! happens to fall on, so curves from different diagnostic days are
    ! directly comparable and differ only via the cohort's own physiological/
    ! canopy state. A realistic-seasonal-sun-angle option is deliberately not
    ! offered - it would reintroduce exactly the confound this frozen-cohort,
    ! fixed-geometry design exists to remove.
    !
    ! AttenuateCanopy is called into local parsun_z/parsha_z/laisun_z/laisha_z
    ! here, never into light_env's own persistent per-substep state, so that is
    ! never perturbed. light_env%twostr's internal solved state is still
    ! touched by each sweep step (ZenithPrep/Solve are not side-effect-free) -
    ! restored at the end by calling light_env%Profile with the real
    ! (light_frac_val, day_of_year, hour_of_day), which unconditionally
    ! re-solves it fresh regardless of what it was left at.

    ! ARGUMENTS:
    type(fates_cohort_type), intent(in)    :: cohort         ! cohort to diagnose (frozen - read-only)
    integer,                  intent(in)    :: pft             ! plant functional type index
    type(light_env_type),    intent(inout) :: light_env       ! this light level's light environment (mutated only transiently; restored before return)
    type(cohort_phys_type),  intent(in)    :: phys            ! today's already-computed leaf physiology (frozen - read-only)
    type(environment_type),  intent(in)    :: env             ! prescribed atmospheric/soil boundary conditions (frozen - read-only)
    real(r8),                 intent(in)    :: ppfd_values(:)  ! incident PPFD values to sweep, at the top of the crown [umol/m2/s]
    real(r8),                 intent(in)    :: light_frac_val  ! this light level's real incident light fraction [0-1] - used only to restore light_env afterward
    integer,                  intent(in)    :: day_of_year     ! real day of year [1-365] - used only to restore light_env afterward
    real(r8),                 intent(in)    :: hour_of_day     ! real hour of day [0-24] - used only to restore light_env afterward
    real(r8),                 intent(out)   :: gross_assim(:)  ! whole-plant gross assimilation at each swept PPFD [kgC/indiv/s], size(ppfd_values)
    real(r8),                 intent(out)   :: total_resp(:)   ! whole-plant total respiration at each swept PPFD [kgC/indiv/s], size(ppfd_values)
    logical,                  intent(in), optional :: use_beam_illumination ! if true, sweep with all incident PPFD as direct beam at zenith instead of the default pure diffuse (Sterck et al. 2013)

    ! LOCALS:
    real(r8), parameter :: sweep_coszen = 1.0_r8 ! fixed sun-angle geometry for the sweep (sun directly overhead) - see header comment
    logical  :: use_beam ! resolved value of use_beam_illumination (defaults to false: pure diffuse)
    real(r8), allocatable :: parsun_z(:), parsha_z(:), laisun_z(:), laisha_z(:) ! this sweep step's attenuated PAR profile - local, never light_env%parsun_z etc.
    real(r8) :: par_toc  ! this sweep step's total incident PAR at the top of the crown [W/m2]
    real(r8) :: par_beam ! this sweep step's direct-beam incident PAR at the top of the crown [W/m2]
    real(r8) :: par_diff ! this sweep step's diffuse incident PAR at the top of the crown [W/m2]
    integer  :: ippfd   ! PPFD-sweep looping index

    use_beam = .false.
    if (present(use_beam_illumination)) use_beam = use_beam_illumination

    allocate(parsun_z(cohort%nv), parsha_z(cohort%nv))
    allocate(laisun_z(cohort%nv), laisha_z(cohort%nv))

    do ippfd = 1, size(ppfd_values)
      par_toc = ppfd_values(ippfd) / wm2_to_umolm2s
      if (use_beam) then
        par_beam = par_toc
        par_diff = 0.0_r8
      else
        par_beam = 0.0_r8
        par_diff = par_toc
      end if
      call light_env%AttenuateCanopy(par_beam, par_diff, sweep_coszen,          &
        parsun_z, parsha_z, laisun_z, laisha_z)
      call phys%GrossAssimAndResp(cohort, pft, env, parsun_z, parsha_z,        &
        laisun_z, laisha_z, step_size, gross_assim(ippfd), total_resp(ippfd))
    end do

    deallocate(parsun_z, parsha_z, laisun_z, laisha_z)

    ! restore light_env's real per-substep state (see header comment)
    call light_env%Profile(light_frac_val, day_of_year, hour_of_day)

  end subroutine LightResponseSweep

  ! ==========================================================================

  subroutine DailyGrowthAndMortality(cohort, daily_net_c, frac_store,          &
    npp_acc_to_prt, cmort, growth_resp, leaf_turnover, fnrt_turnover,          &
    sapw_turnover, struct_turnover)
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
    real(r8),                 intent(out)   :: growth_resp    ! today's growth respiration [kgC/indiv/day]
    real(r8),                 intent(out)   :: leaf_turnover   ! today's leaf turnover loss, from PRTMaintTurnover [kgC/indiv/day]
    real(r8),                 intent(out)   :: fnrt_turnover   ! today's fine root turnover loss, from PRTMaintTurnover [kgC/indiv/day]
    real(r8),                 intent(out)   :: sapw_turnover   ! today's sapwood turnover loss, from PRTMaintTurnover [kgC/indiv/day]
    real(r8),                 intent(out)   :: struct_turnover ! today's structural turnover loss, from PRTMaintTurnover [kgC/indiv/day]

    ! LOCALS:
    real(r8) :: resp_g_acc              ! growth respiration: a tax on the positive part of daily_net_c [kgC/indiv/day], matches EDMainMod.F90's resp_g_acc_hold (without that field's day/year unit round-trip, since it cancels out)
    real(r8) :: delta_dbh, delta_height ! unused outputs of EvaluateAndCorrectDBH (it corrects cohort%dbh/height in place; these report how much it moved them)
    real(r8) :: leaf_c_before, fnrt_c_before, sapw_c_before, struct_c_before ! tissue carbon immediately before PRTMaintTurnover [kgC/indiv] - PRTMaintTurnover has no flux-reporting output of its own, so turnover loss is captured as the state change it causes

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
    growth_resp = resp_g_acc

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
    ! (btran=1) never triggers phenological drought status. PRTMaintTurnover
    ! reports no flux of its own, so per-organ turnover loss is captured here
    ! as the pool-state change it causes
    leaf_c_before = cohort%prt%GetState(leaf_organ, carbon12_element)
    fnrt_c_before = cohort%prt%GetState(fnrt_organ, carbon12_element)
    sapw_c_before = cohort%prt%GetState(sapw_organ, carbon12_element)
    struct_c_before = cohort%prt%GetState(struct_organ, carbon12_element)
    call PRTMaintTurnover(cohort%prt, pft, cohort%canopy_layer, .false.)
    leaf_turnover = leaf_c_before - cohort%prt%GetState(leaf_organ, carbon12_element)
    fnrt_turnover = fnrt_c_before - cohort%prt%GetState(fnrt_organ, carbon12_element)
    sapw_turnover = sapw_c_before - cohort%prt%GetState(sapw_organ, carbon12_element)
    struct_turnover = struct_c_before - cohort%prt%GetState(struct_organ, carbon12_element)
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
