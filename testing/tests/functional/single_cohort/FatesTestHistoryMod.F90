module FatesTestHistoryMod
  !
  ! DESCRIPTION:
  ! Output accumulation and netCDF writing for standalone, patch-less/site-less
  ! cohort test drivers that sweep a set of independent trajectories (e.g. a light
  ! level sweep). Two groups of variables are accumulated in memory (via RecordDay/
  ! RecordLightProfile) and written once at the end (via Write):
  !   - Daily whole-cohort time series, dimensioned (time, light_level), where time
  !     is the day index within each light level's independent trajectory: dbh,
  !     height, treelai, treesai, nv (occupied leaf+stem layers), crown area, the
  !     PARTEH carbon pools (leaf/fine root/sapwood/structure/storage/
  !     reproduction), frac_store (storage as a fraction of target leaf carbon),
  !     cohort number density (n), and the daily flux terms kept disaggregated
  !     rather than netted - GPP, leaf dark respiration, live stem/live coarse
  !     root/fine root maintenance respiration, growth respiration, per-organ
  !     turnover loss (leaf/fine root/sapwood/structure), npp_acc as handed to
  !     PARTEH (net of growth respiration), and the carbon starvation mortality
  !     rate (cmort). The disaggregation is deliberate: decomposing the
  !     compensation point into leaf, architectural (stem/coarse root), and
  !     fine-root support-respiration contributions is impossible from a netted
  !     total. daily_net_c is also kept, as a convenience total equal to
  !     daily_gpp - daily_rdark - daily_livestem_mr - daily_livecroot_mr -
  !     daily_froot_mr by construction (same underlying per-substep terms).
  !     daily_absorbed_par_indiv (Onoda et al. 2013's Phi, whole-plant absorbed
  !     PAR per individual) is also included here, left undivided by any
  !     particular leaf-area or biomass normalization - see
  !     test_SingleCohort.F90's header comment for how LIE_LA/LIE_M/LUE are
  !     built from it in post-processing.
  !   - The daily environmental forcing the trajectory was actually driven with
  !     (daily/midday temperature, saturation vapor pressure and canopy vapor
  !     pressure, the T_growth/T_home acclimation state, and a count of substeps
  !     at which the vapor pressure had to be constrained), dimensioned (time)
  !     alone. These are deliberately NOT on the light_level axis: the prescribed
  !     forcing is a function of day and hour only (see FatesTestEnvironmentMod's
  !     SetHour), so it is identical across light levels by construction and
  !     carrying n_light identical copies would be misleading as well as
  !     wasteful. They exist so that a finished run can be interpreted from its
  !     output file alone - which humidity_mode was in effect, whether a site
  !     swap pushed the vapor pressure against a bound - rather than by
  !     re-deriving the climatology by hand from the namelist. Saturation vapor
  !     pressure is carried alongside vapor pressure rather than a precomputed
  !     VPD or RH so post-processing can form either exactly. Like the Ci-solver
  !     summary below, these are kept regardless of reduced_output: one value per
  !     day is negligible next to the (time, light_level) arrays, and knowing
  !     what a draw was forced with is worth more in a large sweep, not less.
  !   - A per-leaf-layer light profile (parsun_z/parsha_z/laisun_z/laisha_z),
  !     dimensioned (nlevleaf, year, light_level), captured once per year. nlevleaf
  !     is a compile-time maximum leaf+stem layer count; layers above a given
  !     year's cohort%nv are left at the fates_unset_r8 fill value (registered as
  !     each variable's _FillValue attribute), so the array survives nv changing
  !     over time and differing across light levels.
  !   - An instantaneous whole-plant light-response diagnostic (gross_assim/
  !     total_resp, via RecordLightResponse), dimensioned (ppfd, year,
  !     light_level), captured at the same once-per-year cadence as the light
  !     profile above (see test_SingleCohort.F90's LightResponseSweep).
  !   - A companion leaf-level light-response diagnostic (leaf_anet, via
  !     RecordLeafNetAssim), dimensioned (ppfd, year, light_level), at the same
  !     cadence (see test_SingleCohort.F90's LeafNetAssimSweep) - a single
  !     canopy layer's net photosynthesis swept with no canopy attenuation,
  !     unlike the whole-plant diagnostic above.
  !
  ! REDUCED OUTPUT: Init's optional reduced_output argument (default .false.,
  ! full output) switches off everything except the variables actually needed
  ! to compute the shade-tolerance metric set single_cohort_test.py builds:
  ! dbh, treelai, crown_area, n, leaf_c, fnrt_c, sapw_c, struct_c, storage_c,
  ! light_intercept_eff and daily_absorbed_par_indiv (all time x light_level),
  ! plus gross_assim/total_resp/leaf_anet (ppfd x year x light_level). The
  ! per-leaf-layer light profile (parsun_z/parsha_z/laisun_z/laisha_z,
  ! nlevleaf x year x light_level - the single biggest array group) and every
  ! other daily diagnostic (height, treesai, nv, repro_c, the disaggregated
  ! daily flux terms, per-organ turnover, npp_acc, frac_store, cmort,
  ! maintresp_reduction_factor) are skipped entirely: not allocated, not
  ! recorded, not written. This only changes what gets recorded/written -
  ! every one of these quantities is still computed exactly as before
  ! wherever it feeds real cohort dynamics (e.g. daily_net_c still drives
  ! cohort%npp_acc via DailyGrowthAndMortality regardless of this flag); this
  ! module only ever discards values after the physics that needed them has
  ! already run. Intended for the eventual large parameter sweep (Morris/LHC)
  ! over the tolerance-inversion project, where per-run file size and netCDF
  ! write time compound across thousands of draws - mean_ci_solve_iter/
  ! n_bisection_fallbacks (light_level only, negligible size) are kept
  ! regardless, since they are a cheap numerical-health check worth having on
  ! every draw.
  !

  use FatesConstantsMod,  only : r8 => fates_r8
  use FatesConstantsMod,  only : fates_unset_r8
  use FatesCohortMod,     only : fates_cohort_type
  use PRTGenericMod,      only : leaf_organ, fnrt_organ, sapw_organ, struct_organ, store_organ
  use PRTGenericMod,      only : repro_organ
  use PRTGenericMod,      only : carbon12_element
  use FatesTestLightEnvMod, only : light_env_type
  use FatesTestEnvironmentMod, only : environment_type
  use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod, only : WriteVar, RegisterVarAtts, RegisterFillValue, EndNCDef
  use FatesUnitTestIOMod, only : type_double, type_int

  implicit none
  private

  type, public :: history_type

     private

     integer :: n_time  ! days per light level's trajectory (n_days_total)
     integer :: n_leaf   ! compile-time max leaf+stem layers (nlevleaf)
     integer :: n_light  ! number of light levels
     integer :: n_year   ! number of simulated years
     integer :: n_ppfd   ! number of swept PPFD values in the light-response diagnostic
     logical :: reduced_output ! if true, skip allocating/recording/writing everything but the shade-tolerance metric set (see the module header comment)

     ! daily whole-cohort time series, dimensioned (time, light_level)
     real(r8), allocatable :: dbh(:,:)         ! dbh [cm]
     real(r8), allocatable :: height(:,:)      ! height [m]
     real(r8), allocatable :: treelai(:,:)     ! total leaf area index [m2/m2]
     real(r8), allocatable :: treesai(:,:)     ! total stem area index [m2/m2]
     integer,  allocatable :: nv(:,:)          ! number of occupied leaf+stem layers [-]
     real(r8), allocatable :: crown_area(:,:)  ! crown area [m2]
     real(r8), allocatable :: leaf_c(:,:)      ! leaf carbon [kgC/indiv]
     real(r8), allocatable :: fnrt_c(:,:)      ! fine root carbon [kgC/indiv]
     real(r8), allocatable :: sapw_c(:,:)      ! sapwood carbon [kgC/indiv]
     real(r8), allocatable :: struct_c(:,:)    ! structural carbon [kgC/indiv]
     real(r8), allocatable :: storage_c(:,:)   ! storage carbon [kgC/indiv]
     real(r8), allocatable :: repro_c(:,:)     ! reproductive carbon [kgC/indiv]
     real(r8), allocatable :: daily_net_c(:,:) ! daily net C (GPP - leaf Rd - nonleaf MR) [kgC/indiv/day]
     real(r8), allocatable :: daily_gpp(:,:)   ! daily GPP [kgC/indiv/day]
     real(r8), allocatable :: daily_rdark(:,:)      ! daily leaf dark respiration [kgC/indiv/day]
     real(r8), allocatable :: daily_livestem_mr(:,:)  ! daily live stem (aboveground sapwood) maintenance respiration [kgC/indiv/day]
     real(r8), allocatable :: daily_livecroot_mr(:,:) ! daily live coarse root (belowground sapwood) maintenance respiration [kgC/indiv/day]
     real(r8), allocatable :: daily_froot_mr(:,:)     ! daily fine root maintenance respiration [kgC/indiv/day]
     real(r8), allocatable :: daily_growth_resp(:,:)  ! daily growth respiration [kgC/indiv/day]
     real(r8), allocatable :: leaf_turnover(:,:)   ! daily leaf turnover loss [kgC/indiv/day]
     real(r8), allocatable :: fnrt_turnover(:,:)   ! daily fine root turnover loss [kgC/indiv/day]
     real(r8), allocatable :: sapw_turnover(:,:)   ! daily sapwood turnover loss [kgC/indiv/day]
     real(r8), allocatable :: struct_turnover(:,:) ! daily structural turnover loss [kgC/indiv/day]
     real(r8), allocatable :: npp_acc(:,:)     ! carbon balance handed to PARTEH (npp_acc, net of growth respiration) [kgC/indiv/day]
     real(r8), allocatable :: frac_store(:,:)  ! storage carbon as a fraction of target leaf carbon [-]
     real(r8), allocatable :: cmort(:,:)       ! carbon starvation mortality rate [indiv/year]
     real(r8), allocatable :: n(:,:)           ! cohort number density, surviving fraction of the original recruitment cohort [indiv]
     real(r8), allocatable :: light_intercept_eff(:,:)       ! light interception efficiency: whole-plant absorbed PAR / absorbed PAR of an equal-leaf-area, zero-self-shading reference surface, energy-weighted over the day (Sterck et al. 2013) [-]
     real(r8), allocatable :: maintresp_reduction_factor(:,:) ! storage-based maintenance-respiration throttle [0-1]
     real(r8), allocatable :: daily_absorbed_par_indiv(:,:)  ! whole-plant absorbed PAR per individual (Onoda et al. 2013's Phi), integrated over the day [J indiv-1 day-1]

     ! daily environmental forcing, dimensioned (time) only - identical across
     ! light levels by construction, see the module header comment
     real(r8), allocatable :: daily_temp(:)         ! daily-mean vegetation temperature [K]
     real(r8), allocatable :: daily_veg_esat(:)     ! daily-mean saturation vapor pressure [Pa]
     real(r8), allocatable :: daily_can_vpress(:)   ! daily-mean canopy air vapor pressure [Pa]
     real(r8), allocatable :: midday_temp(:)        ! vegetation temperature at the substep nearest solar noon [K]
     real(r8), allocatable :: midday_veg_esat(:)    ! saturation vapor pressure at that substep [Pa]
     real(r8), allocatable :: midday_can_vpress(:)  ! canopy air vapor pressure at that substep [Pa]
     real(r8), allocatable :: t_growth(:)           ! 10-day running-mean growth temperature [K]
     real(r8), allocatable :: t_home(:)             ! long-term running-mean home temperature [K]
     integer,  allocatable :: n_vpress_constrained(:) ! substeps this day at which the requested vapor pressure hit a bound [-]

     ! whole-trajectory Ci-solver summary, dimensioned (light_level) - see
     ! RecordLightLevelSummary
     real(r8), allocatable :: mean_ci_solve_iter(:)   ! mean Ci-solver iteration count over every LeafLayerPhotosynthesis call in this trajectory [-]
     integer,  allocatable :: n_bisection_fallbacks(:) ! count of those calls that fell back to CiBisection [-]

     ! per-leaf-layer light profile, dimensioned (nlevleaf, year, light_level)
     real(r8), allocatable :: parsun_z(:,:,:) ! absorbed PAR, sunlit leaves [W/m2 crown footprint]
     real(r8), allocatable :: parsha_z(:,:,:) ! absorbed PAR, shaded leaves [W/m2 crown footprint]
     real(r8), allocatable :: laisun_z(:,:,:) ! sunlit LAI [m2/m2]
     real(r8), allocatable :: laisha_z(:,:,:) ! shaded LAI [m2/m2]

     ! instantaneous light-response diagnostic, dimensioned (ppfd, year, light_level)
     real(r8), allocatable :: gross_assim(:,:,:) ! whole-plant gross assimilation [kgC/indiv/s]
     real(r8), allocatable :: total_resp(:,:,:)  ! whole-plant total respiration (leaf dark + non-leaf maintenance) [kgC/indiv/s]
     real(r8), allocatable :: leaf_anet(:,:,:)   ! leaf-level net photosynthesis (Aarea), single canopy layer, no canopy attenuation (Sterck et al. 2013's LCPleaf diagnostic) [umolC/m2 leaf/s]

   contains

     procedure, public :: Init
     procedure, public :: RecordDay
     procedure, public :: RecordLightProfile
     procedure, public :: RecordLightResponse
     procedure, public :: RecordLeafNetAssim
     procedure, public :: RecordLightLevelSummary
     procedure, public :: Write

  end type history_type

contains

  ! ==========================================================================

  subroutine Init(this, n_time, n_leaf, n_light, n_year, n_ppfd, reduced_output)
    !
    ! DESCRIPTION:
    ! Allocate the time-series arrays and pre-fill the light-profile arrays with
    ! the fill value - only entries at leaf layers 1:nv, for each (year, light
    ! level) actually recorded, ever get overwritten by RecordLightProfile.
    ! reduced_output (see the module header comment) skips allocating
    ! everything except the shade-tolerance metric set - RecordDay/
    ! RecordLightProfile/Write check this%reduced_output themselves before
    ! touching the skipped arrays, so nothing here ever writes into memory
    ! that wasn't allocated.

    ! ARGUMENTS:
    class(history_type), intent(inout) :: this    ! history object
    integer,              intent(in)    :: n_time  ! days per light level's trajectory
    integer,              intent(in)    :: n_leaf   ! compile-time max leaf+stem layers (nlevleaf)
    integer,              intent(in)    :: n_light  ! number of light levels
    integer,              intent(in)    :: n_year   ! number of simulated years
    integer,              intent(in)    :: n_ppfd   ! number of swept PPFD values in the light-response diagnostic
    logical,              intent(in), optional :: reduced_output ! if true, skip the non-essential output groups (default .false.: full output)

    this%n_time  = n_time
    this%n_leaf  = n_leaf
    this%n_light = n_light
    this%n_year  = n_year
    this%n_ppfd  = n_ppfd

    this%reduced_output = .false.
    if (present(reduced_output)) this%reduced_output = reduced_output

    ! shade-tolerance metric set - always allocated
    allocate(this%dbh(n_time, n_light))
    allocate(this%treelai(n_time, n_light))
    allocate(this%crown_area(n_time, n_light))
    allocate(this%leaf_c(n_time, n_light))
    allocate(this%fnrt_c(n_time, n_light))
    allocate(this%sapw_c(n_time, n_light))
    allocate(this%struct_c(n_time, n_light))
    allocate(this%storage_c(n_time, n_light))
    allocate(this%n(n_time, n_light))
    allocate(this%light_intercept_eff(n_time, n_light))
    allocate(this%daily_absorbed_par_indiv(n_time, n_light))

    allocate(this%gross_assim(n_ppfd, n_year, n_light))
    allocate(this%total_resp(n_ppfd, n_year, n_light))
    allocate(this%leaf_anet(n_ppfd, n_year, n_light))

    ! cheap (light_level only) numerical-health check - always allocated
    allocate(this%mean_ci_solve_iter(n_light))
    allocate(this%n_bisection_fallbacks(n_light))

    ! daily environmental forcing (time only) - always allocated, see the
    ! module header comment
    allocate(this%daily_temp(n_time))
    allocate(this%daily_veg_esat(n_time))
    allocate(this%daily_can_vpress(n_time))
    allocate(this%midday_temp(n_time))
    allocate(this%midday_veg_esat(n_time))
    allocate(this%midday_can_vpress(n_time))
    allocate(this%t_growth(n_time))
    allocate(this%t_home(n_time))
    allocate(this%n_vpress_constrained(n_time))

    ! everything else - skipped entirely in reduced_output mode
    if (.not. this%reduced_output) then
      allocate(this%height(n_time, n_light))
      allocate(this%treesai(n_time, n_light))
      allocate(this%nv(n_time, n_light))
      allocate(this%repro_c(n_time, n_light))
      allocate(this%daily_net_c(n_time, n_light))
      allocate(this%daily_gpp(n_time, n_light))
      allocate(this%daily_rdark(n_time, n_light))
      allocate(this%daily_livestem_mr(n_time, n_light))
      allocate(this%daily_livecroot_mr(n_time, n_light))
      allocate(this%daily_froot_mr(n_time, n_light))
      allocate(this%daily_growth_resp(n_time, n_light))
      allocate(this%leaf_turnover(n_time, n_light))
      allocate(this%fnrt_turnover(n_time, n_light))
      allocate(this%sapw_turnover(n_time, n_light))
      allocate(this%struct_turnover(n_time, n_light))
      allocate(this%npp_acc(n_time, n_light))
      allocate(this%frac_store(n_time, n_light))
      allocate(this%cmort(n_time, n_light))
      allocate(this%maintresp_reduction_factor(n_time, n_light))

      allocate(this%parsun_z(n_leaf, n_year, n_light))
      allocate(this%parsha_z(n_leaf, n_year, n_light))
      allocate(this%laisun_z(n_leaf, n_year, n_light))
      allocate(this%laisha_z(n_leaf, n_year, n_light))
      this%parsun_z(:,:,:) = fates_unset_r8
      this%parsha_z(:,:,:) = fates_unset_r8
      this%laisun_z(:,:,:) = fates_unset_r8
      this%laisha_z(:,:,:) = fates_unset_r8
    end if

  end subroutine Init

  ! ==========================================================================

  subroutine RecordDay(this, iday_all, ilight, cohort, daily_net_c, daily_gpp, &
    daily_rdark, daily_livestem_mr, daily_livecroot_mr, daily_froot_mr,        &
    daily_growth_resp, leaf_turnover, fnrt_turnover, sapw_turnover,            &
    struct_turnover, npp_acc, frac_store, cmort, light_intercept_eff,          &
    maintresp_reduction_factor, daily_absorbed_par_indiv, env)
    !
    ! DESCRIPTION:
    ! Capture one day's daily time series row. dbh/height/treelai/treesai/nv/
    ! crown_area and the PARTEH carbon pools (including repro_c) are read
    ! directly off cohort, since they are simply today's state; the daily flux
    ! terms (daily_net_c and its disaggregated components, growth respiration,
    ! per-organ turnover loss) are passed in, since they are accumulated (or
    ! computed once) over the course of today's sub-daily loop/growth sequence
    ! in RunOneLightLevel/DailyGrowthAndMortality, not recoverable from cohort's
    ! instantaneous state alone.
    !
    ! The environmental forcing is read off env, whose daily diagnostics are
    ! finalized by UpdateDailyMeans - so this must be called after it, which is
    ! also where the daily flux terms above become complete. Those diagnostics
    ! carry no light_level axis (the forcing is identical across light levels,
    ! see the module header), so every light level simply writes the same row;
    ! this is redundant rather than wrong, and keeps the caller from having to
    ! special-case one light level.

    ! ARGUMENTS:
    class(history_type),     intent(inout) :: this        ! history object
    integer,                  intent(in)    :: iday_all    ! day index within this light level's trajectory (1..n_time)
    integer,                  intent(in)    :: ilight      ! light-level index
    type(fates_cohort_type), intent(in)    :: cohort      ! cohort to record
    real(r8),                 intent(in)    :: daily_net_c ! daily net C (GPP - leaf Rd - nonleaf MR) [kgC/indiv/day]
    real(r8),                 intent(in)    :: daily_gpp   ! daily GPP [kgC/indiv/day]
    real(r8),                 intent(in)    :: daily_rdark        ! daily leaf dark respiration [kgC/indiv/day]
    real(r8),                 intent(in)    :: daily_livestem_mr  ! daily live stem maintenance respiration [kgC/indiv/day]
    real(r8),                 intent(in)    :: daily_livecroot_mr ! daily live coarse root maintenance respiration [kgC/indiv/day]
    real(r8),                 intent(in)    :: daily_froot_mr     ! daily fine root maintenance respiration [kgC/indiv/day]
    real(r8),                 intent(in)    :: daily_growth_resp  ! daily growth respiration [kgC/indiv/day]
    real(r8),                 intent(in)    :: leaf_turnover      ! daily leaf turnover loss [kgC/indiv/day]
    real(r8),                 intent(in)    :: fnrt_turnover      ! daily fine root turnover loss [kgC/indiv/day]
    real(r8),                 intent(in)    :: sapw_turnover      ! daily sapwood turnover loss [kgC/indiv/day]
    real(r8),                 intent(in)    :: struct_turnover    ! daily structural turnover loss [kgC/indiv/day]
    real(r8),                 intent(in)    :: npp_acc     ! carbon balance handed to PARTEH (net of growth respiration) [kgC/indiv/day]
    real(r8),                 intent(in)    :: frac_store  ! storage carbon as a fraction of target leaf carbon [-]
    real(r8),                 intent(in)    :: cmort       ! carbon starvation mortality rate [indiv/year]
    real(r8),                 intent(in)    :: light_intercept_eff       ! light interception efficiency [-]
    real(r8),                 intent(in)    :: maintresp_reduction_factor ! storage-based maintenance-respiration throttle [0-1]
    real(r8),                 intent(in)    :: daily_absorbed_par_indiv  ! whole-plant absorbed PAR per individual (Onoda et al. 2013's Phi) [J indiv-1 day-1]
    type(environment_type),  intent(in)    :: env         ! environment holding today's finalized forcing diagnostics

    ! daily environmental forcing - always recorded (see the module header)
    this%daily_temp(iday_all) = env%daily_temp
    this%daily_veg_esat(iday_all) = env%daily_veg_esat
    this%daily_can_vpress(iday_all) = env%daily_can_vpress
    this%midday_temp(iday_all) = env%midday_temp
    this%midday_veg_esat(iday_all) = env%midday_veg_esat
    this%midday_can_vpress(iday_all) = env%midday_can_vpress
    this%t_growth(iday_all) = env%t_growth
    this%t_home(iday_all) = env%t_home
    this%n_vpress_constrained(iday_all) = env%n_vpress_constrained

    ! shade-tolerance metric set - always recorded
    this%dbh(iday_all, ilight) = cohort%dbh
    this%treelai(iday_all, ilight) = cohort%treelai
    this%crown_area(iday_all, ilight) = cohort%c_area
    this%leaf_c(iday_all, ilight) = cohort%prt%GetState(leaf_organ, carbon12_element)
    this%fnrt_c(iday_all, ilight) = cohort%prt%GetState(fnrt_organ, carbon12_element)
    this%sapw_c(iday_all, ilight) = cohort%prt%GetState(sapw_organ, carbon12_element)
    this%struct_c(iday_all, ilight) = cohort%prt%GetState(struct_organ, carbon12_element)
    this%storage_c(iday_all, ilight) = cohort%prt%GetState(store_organ, carbon12_element)
    this%n(iday_all, ilight) = cohort%n
    this%light_intercept_eff(iday_all, ilight) = light_intercept_eff
    this%daily_absorbed_par_indiv(iday_all, ilight) = daily_absorbed_par_indiv

    ! everything else - skipped entirely in reduced_output mode (see Init)
    if (.not. this%reduced_output) then
      this%height(iday_all, ilight) = cohort%height
      this%treesai(iday_all, ilight) = cohort%treesai
      this%nv(iday_all, ilight) = cohort%nv
      this%repro_c(iday_all, ilight) = cohort%prt%GetState(repro_organ, carbon12_element)
      this%daily_net_c(iday_all, ilight) = daily_net_c
      this%daily_rdark(iday_all, ilight) = daily_rdark
      this%daily_livestem_mr(iday_all, ilight) = daily_livestem_mr
      this%daily_livecroot_mr(iday_all, ilight) = daily_livecroot_mr
      this%daily_froot_mr(iday_all, ilight) = daily_froot_mr
      this%daily_growth_resp(iday_all, ilight) = daily_growth_resp
      this%leaf_turnover(iday_all, ilight) = leaf_turnover
      this%fnrt_turnover(iday_all, ilight) = fnrt_turnover
      this%sapw_turnover(iday_all, ilight) = sapw_turnover
      this%struct_turnover(iday_all, ilight) = struct_turnover
      this%daily_gpp(iday_all, ilight) = daily_gpp
      this%npp_acc(iday_all, ilight) = npp_acc
      this%frac_store(iday_all, ilight) = frac_store
      this%cmort(iday_all, ilight) = cmort
      this%maintresp_reduction_factor(iday_all, ilight) = maintresp_reduction_factor
    end if

  end subroutine RecordDay

  ! ==========================================================================

  subroutine RecordLightProfile(this, iyear, ilight, cohort, light_env)
    !
    ! DESCRIPTION:
    ! Capture the per-leaf-layer light profile snapshot for the given year and
    ! light level. cohort%nv <= n_leaf always (n_leaf is a compile-time maximum),
    ! so layers nv+1:n_leaf are simply never touched here and keep the fill value
    ! set by Init. A no-op in reduced_output mode (see the module header
    ! comment) - the arrays this writes into are not even allocated then, so
    ! this must return before touching them, not just skip writing them out

    ! ARGUMENTS:
    class(history_type),     intent(inout) :: this      ! history object
    integer,                  intent(in)    :: iyear     ! simulated year index (1..n_year)
    integer,                  intent(in)    :: ilight    ! light-level index
    type(fates_cohort_type), intent(in)    :: cohort    ! cohort this light profile belongs to
    type(light_env_type),    intent(in)    :: light_env ! light environment holding this substep's profile

    if (this%reduced_output) return

    this%parsun_z(1:cohort%nv, iyear, ilight) = light_env%parsun_z(:)
    this%parsha_z(1:cohort%nv, iyear, ilight) = light_env%parsha_z(:)
    this%laisun_z(1:cohort%nv, iyear, ilight) = light_env%laisun_z(:)
    this%laisha_z(1:cohort%nv, iyear, ilight) = light_env%laisha_z(:)

  end subroutine RecordLightProfile

  ! ==========================================================================

  subroutine RecordLightResponse(this, iyear, ilight, gross_assim, total_resp)
    !
    ! DESCRIPTION:
    ! Capture the instantaneous whole-plant light-response diagnostic (see
    ! test_SingleCohort.F90's LightResponseSweep) for the given year and light
    ! level.

    ! ARGUMENTS:
    class(history_type), intent(inout) :: this            ! history object
    integer,               intent(in)    :: iyear           ! simulated year index (1..n_year)
    integer,               intent(in)    :: ilight          ! light-level index
    real(r8),              intent(in)    :: gross_assim(:)  ! whole-plant gross assimilation at each swept PPFD [kgC/indiv/s], size(n_ppfd)
    real(r8),              intent(in)    :: total_resp(:)   ! whole-plant total respiration at each swept PPFD [kgC/indiv/s], size(n_ppfd)

    this%gross_assim(:, iyear, ilight) = gross_assim(:)
    this%total_resp(:, iyear, ilight) = total_resp(:)

  end subroutine RecordLightResponse

  ! ==========================================================================

  subroutine RecordLeafNetAssim(this, iyear, ilight, leaf_anet)
    !
    ! DESCRIPTION:
    ! Capture the instantaneous leaf-level light-response diagnostic (see
    ! test_SingleCohort.F90's LeafNetAssimSweep) for the given year and light
    ! level.

    ! ARGUMENTS:
    class(history_type), intent(inout) :: this          ! history object
    integer,               intent(in)    :: iyear         ! simulated year index (1..n_year)
    integer,               intent(in)    :: ilight        ! light-level index
    real(r8),              intent(in)    :: leaf_anet(:)  ! leaf-level net photosynthesis at each swept PPFD [umolC/m2 leaf/s], size(n_ppfd)

    this%leaf_anet(:, iyear, ilight) = leaf_anet(:)

  end subroutine RecordLeafNetAssim

  ! ==========================================================================

  subroutine RecordLightLevelSummary(this, ilight, mean_solve_iter, n_bisection_calls)
    !
    ! DESCRIPTION:
    ! Capture the whole-trajectory Ci-solver summary for the given light level
    ! (see test_SingleCohort.F90's RunOneLightLevel) - a single value per light
    ! level, not a daily quantity.

    ! ARGUMENTS:
    class(history_type), intent(inout) :: this              ! history object
    integer,               intent(in)    :: ilight            ! light-level index
    real(r8),              intent(in)    :: mean_solve_iter   ! mean Ci-solver iteration count over this trajectory [-]
    integer,               intent(in)    :: n_bisection_calls ! count of calls that fell back to CiBisection over this trajectory [-]

    this%mean_ci_solve_iter(ilight) = mean_solve_iter
    this%n_bisection_fallbacks(ilight) = n_bisection_calls

  end subroutine RecordLightLevelSummary

  ! ==========================================================================

  subroutine Write(this, out_file, light_frac, ppfd_values)
    !
    ! DESCRIPTION:
    ! Writes out the daily whole-cohort time series (see RecordDay), the daily
    ! environmental forcing (dimensioned (time) alone), the instantaneous
    ! whole-plant/leaf-level light-response diagnostics (gross_assim/total_resp/
    ! leaf_anet, dimensioned (ppfd, year, light_level)), and the whole-trajectory
    ! Ci-solver summary (mean_ci_solve_iter/n_bisection_fallbacks, dimensioned
    ! (light_level)) - these are the shade-tolerance metric set plus the two
    ! cheap always-on diagnostic groups (see the module header comment), and are
    ! always registered/written. Everything else (the per-leaf-layer light profile,
    ! dimensioned (nlevleaf, year, light_level), and the remaining daily
    ! diagnostics) is registered/written only if this%reduced_output is false -
    ! grouped into one block below, rather than guarding each RegisterVarAtts/
    ! WriteVar call individually, since Fortran has no conditional variable
    ! declaration and the two groups are otherwise interleaved in a natural
    ! reading order.

    ! ARGUMENTS:
    class(history_type), intent(in) :: this           ! history object
    character(len=*),     intent(in) :: out_file       ! output file name
    real(r8),              intent(in) :: light_frac(:)  ! swept incident light fractions [0-1]
    real(r8),              intent(in) :: ppfd_values(:) ! swept PPFD values in the light-response diagnostic [umol/m2/s]

    ! LOCALS:
    integer, allocatable :: time_idx(:)   ! day index, 1 = first simulated day
    integer, allocatable :: leaf_layer(:) ! leaf-layer index, 1 = top of crown
    integer, allocatable :: year_idx(:)   ! year index, 1 = first simulated year
    integer               :: i            ! looping index
    integer               :: ncid         ! netcdf file id
    character(len=20)     :: dim_names(5) ! dimension names
    integer               :: dimIDs(5)    ! dimension IDs
    integer               :: timeID, leaflayerID, lightfracID, yearID, ppfdID
    integer               :: dbhID, treelaiID, crownareaID
    integer               :: leafcID, fnrtcID, sapwcID, structcID, storagecID
    integer               :: nID, lightintercepteffID, dailyabsorbedparindivID
    integer               :: meancisolveiterID, nbisectionfallbacksID
    integer               :: dailytempID, dailyvegesatID, dailycanvpressID
    integer               :: middaytempID, middayvegesatID, middaycanvpressID
    integer               :: tgrowthID, thomeID, nvpressconstrainedID
    integer               :: grossassimID, totalrespID, leafanetID
    integer               :: heightID, treesaiID, nvID, reprocID
    integer               :: dailynetcID, dailygppID, nppaccID, fracstoreID
    integer               :: dailyrdarkID, dailylivestemmrID, dailylivecrootmrID
    integer               :: dailyfrootmrID, dailygrowthrespID
    integer               :: leafturnoverID, fnrtturnoverID, sapwturnoverID, structturnoverID
    integer               :: cmortID, maintrespreductionfactorID
    integer               :: parsunID, parshaID, laisunID, laishaID

    ! create day, leaf layer, and year indices
    allocate(time_idx(this%n_time))
    do i = 1, this%n_time
      time_idx(i) = i
    end do
    allocate(leaf_layer(this%n_leaf))
    do i = 1, this%n_leaf
      leaf_layer(i) = i
    end do
    allocate(year_idx(this%n_year))
    do i = 1, this%n_year
      year_idx(i) = i
    end do

    ! dimension names - nlevleaf is registered regardless of reduced_output
    ! (its own cost is negligible; only the 4 big arrays it dimensions are
    ! actually skipped)
    dim_names = [character(len=20) :: 'time', 'nlevleaf', 'light_level', 'year', 'ppfd']

    ! open file
    call OpenNCFile(trim(out_file), ncid, 'readwrite')

    ! register dimensions
    call RegisterNCDims(ncid, dim_names,                                       &
      (/this%n_time, this%n_leaf, this%n_light, this%n_year, this%n_ppfd/), 5, dimIDs)

    ! register day index
    call RegisterVarAtts(ncid, dim_names(1), dimIDs(1:1), type_int, '',                 &
      'day index within this light level''s trajectory, 1 = first simulated day',       &
      timeID)

    ! register leaf layer
    call RegisterVarAtts(ncid, dim_names(2), dimIDs(2:2), type_int, '',                 &
      'leaf layer index, 1 = top of crown (compile-time max nlevleaf; unoccupied layers filled)', &
      leaflayerID)

    ! register light level
    call RegisterVarAtts(ncid, dim_names(3), dimIDs(3:3), type_double, 'fraction',      &
      'incident light fraction, relative to full sun', lightfracID)

    ! register year
    call RegisterVarAtts(ncid, dim_names(4), dimIDs(4:4), type_int, '',                 &
      'simulated year index, 1 = first year', yearID)

    ! register ppfd
    call RegisterVarAtts(ncid, dim_names(5), dimIDs(5:5), type_double, 'umol m-2 s-1',  &
      'incident PPFD swept by the instantaneous light-response diagnostic', ppfdID)

    ! ------------------------------------------------------------------
    ! shade-tolerance metric set - always registered (see module header)
    ! ------------------------------------------------------------------

    call RegisterVarAtts(ncid, 'dbh', (/dimIDs(1), dimIDs(3)/), type_double, 'cm',      &
      'dbh', dbhID, coordinates='time light_level')

    call RegisterVarAtts(ncid, 'treelai', (/dimIDs(1), dimIDs(3)/), type_double,        &
      'm2 m-2', 'total leaf area index', treelaiID, coordinates='time light_level')

    call RegisterVarAtts(ncid, 'crown_area', (/dimIDs(1), dimIDs(3)/), type_double,     &
      'm2', 'crown area', crownareaID, coordinates='time light_level')

    call RegisterVarAtts(ncid, 'leaf_c', (/dimIDs(1), dimIDs(3)/), type_double,         &
      'kgC indiv-1', 'leaf carbon', leafcID, coordinates='time light_level')

    call RegisterVarAtts(ncid, 'fnrt_c', (/dimIDs(1), dimIDs(3)/), type_double,         &
      'kgC indiv-1', 'fine root carbon', fnrtcID, coordinates='time light_level')

    call RegisterVarAtts(ncid, 'sapw_c', (/dimIDs(1), dimIDs(3)/), type_double,         &
      'kgC indiv-1', 'sapwood carbon', sapwcID, coordinates='time light_level')

    call RegisterVarAtts(ncid, 'struct_c', (/dimIDs(1), dimIDs(3)/), type_double,       &
      'kgC indiv-1', 'structural carbon', structcID, coordinates='time light_level')

    call RegisterVarAtts(ncid, 'storage_c', (/dimIDs(1), dimIDs(3)/), type_double,      &
      'kgC indiv-1', 'storage carbon', storagecID, coordinates='time light_level')

    call RegisterVarAtts(ncid, 'n', (/dimIDs(1), dimIDs(3)/), type_double, 'indiv',     &
      'cohort number density (surviving fraction of the original recruitment cohort)',  &
      nID, coordinates='time light_level')

    call RegisterVarAtts(ncid, 'light_intercept_eff', (/dimIDs(1), dimIDs(3)/),         &
      type_double, '-',                                                                 &
      'light interception efficiency: whole-plant absorbed PAR / absorbed PAR of an equal-leaf-area, zero-self-shading reference surface, energy-weighted over the day (Sterck et al. 2013)', &
      lightintercepteffID, coordinates='time light_level')

    call RegisterVarAtts(ncid, 'daily_absorbed_par_indiv', (/dimIDs(1), dimIDs(3)/),    &
      type_double, 'J indiv-1 day-1',                                                   &
      'whole-plant absorbed PAR per individual (Onoda et al. 2013''s Phi) - divide by treelai*crown_area/n for Phi/LA (LIE_LA), by a chosen biomass-pool total for Phi/M (LIE_M), or use its running sum against a biomass increment for LUE', &
      dailyabsorbedparindivID, coordinates='time light_level')

    ! daily environmental forcing, dimensioned (time) only - the forcing is a
    ! function of day and hour alone, so it does not vary across light levels
    ! (see the module header comment). Kept regardless of reduced_output.
    call RegisterVarAtts(ncid, 'daily_temp', dimIDs(1:1), type_double, 'K',             &
      'daily-mean vegetation temperature (mean over the day''s sub-daily substeps)',    &
      dailytempID, coordinates='time')

    call RegisterVarAtts(ncid, 'daily_veg_esat', dimIDs(1:1), type_double, 'Pa',        &
      'daily-mean saturation vapor pressure at the vegetation temperature',             &
      dailyvegesatID, coordinates='time')

    call RegisterVarAtts(ncid, 'daily_can_vpress', dimIDs(1:1), type_double, 'Pa',      &
      'daily-mean canopy air vapor pressure - with daily_veg_esat gives daily-mean VPD (veg_esat - can_vpress) and RH (can_vpress/veg_esat)', &
      dailycanvpressID, coordinates='time')

    call RegisterVarAtts(ncid, 'midday_temp', dimIDs(1:1), type_double, 'K',            &
      'vegetation temperature at the sub-daily substep nearest solar noon',             &
      middaytempID, coordinates='time')

    call RegisterVarAtts(ncid, 'midday_veg_esat', dimIDs(1:1), type_double, 'Pa',       &
      'saturation vapor pressure at the substep nearest solar noon',                    &
      middayvegesatID, coordinates='time')

    call RegisterVarAtts(ncid, 'midday_can_vpress', dimIDs(1:1), type_double, 'Pa',     &
      'canopy air vapor pressure at the substep nearest solar noon - with midday_veg_esat gives the midday VPD that drives stomatal conductance when PAR is highest', &
      middaycanvpressID, coordinates='time')

    call RegisterVarAtts(ncid, 't_growth', dimIDs(1:1), type_double, 'K',               &
      '10-day running-mean growth temperature (photosynthetic acclimation boundary condition)', &
      tgrowthID, coordinates='time')

    call RegisterVarAtts(ncid, 't_home', dimIDs(1:1), type_double, 'K',                 &
      'long-term running-mean home temperature (photosynthetic acclimation boundary condition)', &
      thomeID, coordinates='time')

    call RegisterVarAtts(ncid, 'n_vpress_constrained', dimIDs(1:1), type_int, '-',      &
      'number of sub-daily substeps this day at which the prescribed canopy vapor pressure was altered by the saturation/minimum constraint - nonzero means the run was not forced with the humidity boundary condition it asked for', &
      nvpressconstrainedID, coordinates='time')

    ! whole-trajectory Ci-solver summary, dimensioned (light_level) - cheap,
    ! kept regardless of reduced_output (see module header comment)
    call RegisterVarAtts(ncid, 'mean_ci_solve_iter', dimIDs(3:3), type_double, '-',     &
      'mean Ci-solver iteration count over every LeafLayerPhotosynthesis call in this trajectory', &
      meancisolveiterID, coordinates='light_level')

    call RegisterVarAtts(ncid, 'n_bisection_fallbacks', dimIDs(3:3), type_int, '-',     &
      'count of LeafLayerPhotosynthesis calls that fell back to CiBisection in this trajectory', &
      nbisectionfallbacksID, coordinates='light_level')

    ! instantaneous light-response diagnostics, dimensioned (ppfd, year,
    ! light_level) - see test_SingleCohort.F90's LightResponseSweep/
    ! LeafNetAssimSweep
    call RegisterVarAtts(ncid, 'gross_assim', (/dimIDs(5), dimIDs(4), dimIDs(3)/),      &
      type_double, 'kgC indiv-1 s-1',                                                   &
      'whole-plant gross assimilation at each swept PPFD (first day of year, pure diffuse illumination, coszen=1)', &
      grossassimID, coordinates='ppfd year light_level')

    call RegisterVarAtts(ncid, 'total_resp', (/dimIDs(5), dimIDs(4), dimIDs(3)/),       &
      type_double, 'kgC indiv-1 s-1',                                                   &
      'whole-plant total respiration (leaf dark + non-leaf maintenance) at each swept PPFD (first day of year)', &
      totalrespID, coordinates='ppfd year light_level')

    ! unlike gross_assim/total_resp above (whole-plant, kgC indiv-1 s-1), this
    ! is an area-based leaf-level rate (umol m-2 s-1), swept with no canopy
    ! attenuation, at a single fixed canopy layer (leaf_lcp_layer)
    call RegisterVarAtts(ncid, 'leaf_anet', (/dimIDs(5), dimIDs(4), dimIDs(3)/),        &
      type_double, 'umol m-2 s-1',                                                      &
      'leaf-level net photosynthesis (Aarea) at each swept PPFD, applied directly to a single canopy layer with no canopy attenuation or self-shading (first day of year) (Sterck et al. 2013)', &
      leafanetID, coordinates='ppfd year light_level')

    ! ------------------------------------------------------------------
    ! everything else - skipped entirely in reduced_output mode (see Init)
    ! ------------------------------------------------------------------

    if (.not. this%reduced_output) then

      call RegisterVarAtts(ncid, 'height', (/dimIDs(1), dimIDs(3)/), type_double, 'm',  &
        'height', heightID, coordinates='time light_level')

      call RegisterVarAtts(ncid, 'treesai', (/dimIDs(1), dimIDs(3)/), type_double,      &
        'm2 m-2', 'total stem area index', treesaiID, coordinates='time light_level')

      call RegisterVarAtts(ncid, 'nv', (/dimIDs(1), dimIDs(3)/), type_int, '-',         &
        'number of occupied leaf+stem layers', nvID, coordinates='time light_level')

      call RegisterVarAtts(ncid, 'repro_c', (/dimIDs(1), dimIDs(3)/), type_double,      &
        'kgC indiv-1', 'reproductive carbon', reprocID, coordinates='time light_level')

      call RegisterVarAtts(ncid, 'daily_net_c', (/dimIDs(1), dimIDs(3)/), type_double,  &
        'kgC indiv-1 day-1',                                                            &
        'daily net carbon (GPP - leaf dark resp - nonleaf MR) - equal to daily_gpp - daily_rdark - daily_livestem_mr - daily_livecroot_mr - daily_froot_mr by construction', &
        dailynetcID, coordinates='time light_level')

      call RegisterVarAtts(ncid, 'daily_gpp', (/dimIDs(1), dimIDs(3)/), type_double,    &
        'kgC indiv-1 day-1', 'daily GPP', dailygppID, coordinates='time light_level')

      call RegisterVarAtts(ncid, 'daily_rdark', (/dimIDs(1), dimIDs(3)/), type_double,  &
        'kgC indiv-1 day-1', 'daily leaf dark respiration', dailyrdarkID,               &
        coordinates='time light_level')

      call RegisterVarAtts(ncid, 'daily_livestem_mr', (/dimIDs(1), dimIDs(3)/),         &
        type_double, 'kgC indiv-1 day-1',                                               &
        'daily live stem (aboveground sapwood) maintenance respiration',                &
        dailylivestemmrID, coordinates='time light_level')

      call RegisterVarAtts(ncid, 'daily_livecroot_mr', (/dimIDs(1), dimIDs(3)/),        &
        type_double, 'kgC indiv-1 day-1',                                               &
        'daily live coarse root (belowground sapwood) maintenance respiration',         &
        dailylivecrootmrID, coordinates='time light_level')

      call RegisterVarAtts(ncid, 'daily_froot_mr', (/dimIDs(1), dimIDs(3)/),            &
        type_double, 'kgC indiv-1 day-1', 'daily fine root maintenance respiration',    &
        dailyfrootmrID, coordinates='time light_level')

      call RegisterVarAtts(ncid, 'daily_growth_resp', (/dimIDs(1), dimIDs(3)/),         &
        type_double, 'kgC indiv-1 day-1', 'daily growth respiration',                   &
        dailygrowthrespID, coordinates='time light_level')

      call RegisterVarAtts(ncid, 'leaf_turnover', (/dimIDs(1), dimIDs(3)/),             &
        type_double, 'kgC indiv-1 day-1', 'daily leaf turnover loss', leafturnoverID,   &
        coordinates='time light_level')

      call RegisterVarAtts(ncid, 'fnrt_turnover', (/dimIDs(1), dimIDs(3)/),             &
        type_double, 'kgC indiv-1 day-1', 'daily fine root turnover loss',              &
        fnrtturnoverID, coordinates='time light_level')

      call RegisterVarAtts(ncid, 'sapw_turnover', (/dimIDs(1), dimIDs(3)/),             &
        type_double, 'kgC indiv-1 day-1', 'daily sapwood turnover loss',                &
        sapwturnoverID, coordinates='time light_level')

      call RegisterVarAtts(ncid, 'struct_turnover', (/dimIDs(1), dimIDs(3)/),           &
        type_double, 'kgC indiv-1 day-1', 'daily structural turnover loss',             &
        structturnoverID, coordinates='time light_level')

      call RegisterVarAtts(ncid, 'npp_acc', (/dimIDs(1), dimIDs(3)/), type_double,      &
        'kgC indiv-1 day-1',                                                            &
        'carbon balance handed to PARTEH (net of growth respiration)', nppaccID,        &
        coordinates='time light_level')

      call RegisterVarAtts(ncid, 'frac_store', (/dimIDs(1), dimIDs(3)/), type_double,   &
        '-', 'storage carbon as a fraction of target leaf carbon', fracstoreID,         &
        coordinates='time light_level')

      call RegisterVarAtts(ncid, 'cmort', (/dimIDs(1), dimIDs(3)/), type_double,        &
        'indiv yr-1', 'carbon starvation mortality rate', cmortID,                      &
        coordinates='time light_level')

      call RegisterVarAtts(ncid, 'maintresp_reduction_factor',                          &
        (/dimIDs(1), dimIDs(3)/), type_double, '-',                                     &
        'storage-based maintenance-respiration throttle', maintrespreductionfactorID,   &
        coordinates='time light_level')

      ! per-leaf-layer light profile, dimensioned (nlevleaf, year, light_level),
      ! with unoccupied layers above that year's nv filled with fates_unset_r8
      call RegisterVarAtts(ncid, 'parsun_z', (/dimIDs(2), dimIDs(4), dimIDs(3)/),       &
        type_double, 'W m-2',                                                           &
        'absorbed PAR, sunlit leaves, per unit crown footprint area (first day of year, solar noon)', &
        parsunID, coordinates='nlevleaf year light_level')
      call RegisterFillValue(ncid, parsunID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'parsha_z', (/dimIDs(2), dimIDs(4), dimIDs(3)/),       &
        type_double, 'W m-2',                                                           &
        'absorbed PAR, shaded leaves, per unit crown footprint area (first day of year, solar noon)', &
        parshaID, coordinates='nlevleaf year light_level')
      call RegisterFillValue(ncid, parshaID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'laisun_z', (/dimIDs(2), dimIDs(4), dimIDs(3)/),       &
        type_double, 'm2 m-2', 'sunlit LAI (first day of year, solar noon)', laisunID,  &
        coordinates='nlevleaf year light_level')
      call RegisterFillValue(ncid, laisunID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'laisha_z', (/dimIDs(2), dimIDs(4), dimIDs(3)/),       &
        type_double, 'm2 m-2', 'shaded LAI (first day of year, solar noon)', laishaID,  &
        coordinates='nlevleaf year light_level')
      call RegisterFillValue(ncid, laishaID, fates_unset_r8)

    end if

    ! finish defining variables
    call EndNCDef(ncid)

    ! write out data - shade-tolerance metric set, always written
    call WriteVar(ncid, timeID, time_idx(:))
    call WriteVar(ncid, leaflayerID, leaf_layer(:))
    call WriteVar(ncid, lightfracID, light_frac(:))
    call WriteVar(ncid, yearID, year_idx(:))
    call WriteVar(ncid, dbhID, this%dbh(:,:))
    call WriteVar(ncid, treelaiID, this%treelai(:,:))
    call WriteVar(ncid, crownareaID, this%crown_area(:,:))
    call WriteVar(ncid, leafcID, this%leaf_c(:,:))
    call WriteVar(ncid, fnrtcID, this%fnrt_c(:,:))
    call WriteVar(ncid, sapwcID, this%sapw_c(:,:))
    call WriteVar(ncid, structcID, this%struct_c(:,:))
    call WriteVar(ncid, storagecID, this%storage_c(:,:))
    call WriteVar(ncid, nID, this%n(:,:))
    call WriteVar(ncid, lightintercepteffID, this%light_intercept_eff(:,:))
    call WriteVar(ncid, dailyabsorbedparindivID, this%daily_absorbed_par_indiv(:,:))
    call WriteVar(ncid, dailytempID, this%daily_temp(:))
    call WriteVar(ncid, dailyvegesatID, this%daily_veg_esat(:))
    call WriteVar(ncid, dailycanvpressID, this%daily_can_vpress(:))
    call WriteVar(ncid, middaytempID, this%midday_temp(:))
    call WriteVar(ncid, middayvegesatID, this%midday_veg_esat(:))
    call WriteVar(ncid, middaycanvpressID, this%midday_can_vpress(:))
    call WriteVar(ncid, tgrowthID, this%t_growth(:))
    call WriteVar(ncid, thomeID, this%t_home(:))
    call WriteVar(ncid, nvpressconstrainedID, this%n_vpress_constrained(:))
    call WriteVar(ncid, meancisolveiterID, this%mean_ci_solve_iter(:))
    call WriteVar(ncid, nbisectionfallbacksID, this%n_bisection_fallbacks(:))
    call WriteVar(ncid, ppfdID, ppfd_values(:))
    call WriteVar(ncid, grossassimID, this%gross_assim(:,:,:))
    call WriteVar(ncid, totalrespID, this%total_resp(:,:,:))
    call WriteVar(ncid, leafanetID, this%leaf_anet(:,:,:))

    ! everything else - skipped entirely in reduced_output mode
    if (.not. this%reduced_output) then
      call WriteVar(ncid, heightID, this%height(:,:))
      call WriteVar(ncid, treesaiID, this%treesai(:,:))
      call WriteVar(ncid, nvID, this%nv(:,:))
      call WriteVar(ncid, reprocID, this%repro_c(:,:))
      call WriteVar(ncid, dailynetcID, this%daily_net_c(:,:))
      call WriteVar(ncid, dailygppID, this%daily_gpp(:,:))
      call WriteVar(ncid, dailyrdarkID, this%daily_rdark(:,:))
      call WriteVar(ncid, dailylivestemmrID, this%daily_livestem_mr(:,:))
      call WriteVar(ncid, dailylivecrootmrID, this%daily_livecroot_mr(:,:))
      call WriteVar(ncid, dailyfrootmrID, this%daily_froot_mr(:,:))
      call WriteVar(ncid, dailygrowthrespID, this%daily_growth_resp(:,:))
      call WriteVar(ncid, leafturnoverID, this%leaf_turnover(:,:))
      call WriteVar(ncid, fnrtturnoverID, this%fnrt_turnover(:,:))
      call WriteVar(ncid, sapwturnoverID, this%sapw_turnover(:,:))
      call WriteVar(ncid, structturnoverID, this%struct_turnover(:,:))
      call WriteVar(ncid, nppaccID, this%npp_acc(:,:))
      call WriteVar(ncid, fracstoreID, this%frac_store(:,:))
      call WriteVar(ncid, cmortID, this%cmort(:,:))
      call WriteVar(ncid, maintrespreductionfactorID, this%maintresp_reduction_factor(:,:))
      call WriteVar(ncid, parsunID, this%parsun_z(:,:,:))
      call WriteVar(ncid, parshaID, this%parsha_z(:,:,:))
      call WriteVar(ncid, laisunID, this%laisun_z(:,:,:))
      call WriteVar(ncid, laishaID, this%laisha_z(:,:,:))
    end if

    ! close the file
    call CloseNCFile(ncid)

  end subroutine Write

end module FatesTestHistoryMod
