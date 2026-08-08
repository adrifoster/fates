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
  use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod, only : WriteVar, RegisterVar, RegisterFillValue, EndNCDef
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

     ! whole-trajectory Ci-solver summary, dimensioned (light_level) - see
     ! RecordLightLevelSummary
     real(r8), allocatable :: mean_ci_solve_iter(:)   ! mean Ci-solver iteration count over every LeafLayerPhotosynthesis call in this trajectory [-]
     integer,  allocatable :: n_bisection_fallbacks(:) ! count of those calls that fell back to CiBisection [-]

     ! per-leaf-layer light profile, dimensioned (nlevleaf, year, light_level)
     real(r8), allocatable :: parsun_z(:,:,:) ! absorbed PAR, sunlit leaves [W/m2 ground]
     real(r8), allocatable :: parsha_z(:,:,:) ! absorbed PAR, shaded leaves [W/m2 ground]
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
    maintresp_reduction_factor, daily_absorbed_par_indiv)
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
    ! Writes out the daily whole-cohort time series (see RecordDay), the
    ! instantaneous whole-plant/leaf-level light-response diagnostics
    ! (gross_assim/total_resp/leaf_anet, dimensioned (ppfd, year, light_level)),
    ! and the whole-trajectory Ci-solver summary (mean_ci_solve_iter/
    ! n_bisection_fallbacks, dimensioned (light_level)) - these are the
    ! shade-tolerance metric set (see the module header comment) and are always
    ! registered/written. Everything else (the per-leaf-layer light profile,
    ! dimensioned (nlevleaf, year, light_level), and the remaining daily
    ! diagnostics) is registered/written only if this%reduced_output is false -
    ! grouped into one block below, rather than guarding each RegisterVar/
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
    call RegisterVar(ncid, dim_names(1), dimIDs(1:1), type_int,                &
      [character(len=20)  :: 'units', 'long_name'],                           &
      [character(len=150) :: '', 'day index within this light level''s trajectory, 1 = first simulated day'], &
      2, timeID)

    ! register leaf layer
    call RegisterVar(ncid, dim_names(2), dimIDs(2:2), type_int,                &
      [character(len=20)  :: 'units', 'long_name'],                           &
      [character(len=150) :: '', 'leaf layer index, 1 = top of crown (compile-time max nlevleaf; unoccupied layers filled)'], &
      2, leaflayerID)

    ! register light level
    call RegisterVar(ncid, dim_names(3), dimIDs(3:3), type_double,             &
      [character(len=20)  :: 'units', 'long_name'],                           &
      [character(len=150) :: 'fraction', 'incident light fraction, relative to full sun'], &
      2, lightfracID)

    ! register year
    call RegisterVar(ncid, dim_names(4), dimIDs(4:4), type_int,                &
      [character(len=20)  :: 'units', 'long_name'],                           &
      [character(len=150) :: '', 'simulated year index, 1 = first year'], 2, yearID)

    ! register ppfd
    call RegisterVar(ncid, dim_names(5), dimIDs(5:5), type_double,             &
      [character(len=20)  :: 'units', 'long_name'],                           &
      [character(len=150) :: 'umol m-2 s-1',                                  &
      'incident PPFD swept by the instantaneous light-response diagnostic'],   &
      2, ppfdID)

    ! ------------------------------------------------------------------
    ! shade-tolerance metric set - always registered (see module header)
    ! ------------------------------------------------------------------

    call RegisterVar(ncid, 'dbh', (/dimIDs(1), dimIDs(3)/), type_double,       &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'cm', 'dbh'], 3, dbhID)

    call RegisterVar(ncid, 'treelai', (/dimIDs(1), dimIDs(3)/), type_double,   &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'm2 m-2', 'total leaf area index'], &
      3, treelaiID)

    call RegisterVar(ncid, 'crown_area', (/dimIDs(1), dimIDs(3)/), type_double, &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'm2', 'crown area'], 3, crownareaID)

    call RegisterVar(ncid, 'leaf_c', (/dimIDs(1), dimIDs(3)/), type_double,    &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'kgC indiv-1', 'leaf carbon'], &
      3, leafcID)

    call RegisterVar(ncid, 'fnrt_c', (/dimIDs(1), dimIDs(3)/), type_double,    &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'kgC indiv-1', 'fine root carbon'], &
      3, fnrtcID)

    call RegisterVar(ncid, 'sapw_c', (/dimIDs(1), dimIDs(3)/), type_double,    &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'kgC indiv-1', 'sapwood carbon'], &
      3, sapwcID)

    call RegisterVar(ncid, 'struct_c', (/dimIDs(1), dimIDs(3)/), type_double,  &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'kgC indiv-1', 'structural carbon'], &
      3, structcID)

    call RegisterVar(ncid, 'storage_c', (/dimIDs(1), dimIDs(3)/), type_double, &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'kgC indiv-1', 'storage carbon'], &
      3, storagecID)

    call RegisterVar(ncid, 'n', (/dimIDs(1), dimIDs(3)/), type_double,         &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'indiv',                    &
      'cohort number density (surviving fraction of the original recruitment cohort)'], &
      3, nID)

    call RegisterVar(ncid, 'light_intercept_eff', (/dimIDs(1), dimIDs(3)/), type_double, &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', '-',                        &
      'light interception efficiency: whole-plant absorbed PAR / absorbed PAR of an equal-leaf-area, zero-self-shading reference surface, energy-weighted over the day (Sterck et al. 2013)'], &
      3, lightintercepteffID)

    call RegisterVar(ncid, 'daily_absorbed_par_indiv', (/dimIDs(1), dimIDs(3)/), type_double, &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'J indiv-1 day-1',           &
      'whole-plant absorbed PAR per individual (Onoda et al. 2013''s Phi) - divide by treelai*crown_area/n for Phi/LA (LIE_LA), by a chosen biomass-pool total for Phi/M (LIE_M), or use its running sum against a biomass increment for LUE'], &
      3, dailyabsorbedparindivID)

    ! whole-trajectory Ci-solver summary, dimensioned (light_level) - cheap,
    ! kept regardless of reduced_output (see module header comment)
    call RegisterVar(ncid, 'mean_ci_solve_iter', dimIDs(3:3), type_double,     &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'light_level', '-',                             &
      'mean Ci-solver iteration count over every LeafLayerPhotosynthesis call in this trajectory'], &
      3, meancisolveiterID)

    call RegisterVar(ncid, 'n_bisection_fallbacks', dimIDs(3:3), type_int,     &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'light_level', '-',                             &
      'count of LeafLayerPhotosynthesis calls that fell back to CiBisection in this trajectory'], &
      3, nbisectionfallbacksID)

    ! instantaneous light-response diagnostics, dimensioned (ppfd, year,
    ! light_level) - see test_SingleCohort.F90's LightResponseSweep/
    ! LeafNetAssimSweep
    call RegisterVar(ncid, 'gross_assim', (/dimIDs(5), dimIDs(4), dimIDs(3)/), type_double, &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'ppfd year light_level', 'kgC indiv-1 s-1',      &
      'whole-plant gross assimilation at each swept PPFD (first day of year, pure diffuse illumination, coszen=1)'], &
      3, grossassimID)

    call RegisterVar(ncid, 'total_resp', (/dimIDs(5), dimIDs(4), dimIDs(3)/), type_double, &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'ppfd year light_level', 'kgC indiv-1 s-1',      &
      'whole-plant total respiration (leaf dark + non-leaf maintenance) at each swept PPFD (first day of year)'], &
      3, totalrespID)

    ! unlike gross_assim/total_resp above (whole-plant, kgC indiv-1 s-1), this
    ! is an area-based leaf-level rate (umol m-2 s-1), swept with no canopy
    ! attenuation, at a single fixed canopy layer (leaf_lcp_layer)
    call RegisterVar(ncid, 'leaf_anet', (/dimIDs(5), dimIDs(4), dimIDs(3)/), type_double, &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'ppfd year light_level', 'umol m-2 s-1',         &
      'leaf-level net photosynthesis (Aarea) at each swept PPFD, applied directly to a single canopy layer with no canopy attenuation or self-shading (first day of year) (Sterck et al. 2013)'], &
      3, leafanetID)

    ! ------------------------------------------------------------------
    ! everything else - skipped entirely in reduced_output mode (see Init)
    ! ------------------------------------------------------------------

    if (.not. this%reduced_output) then

      call RegisterVar(ncid, 'height', (/dimIDs(1), dimIDs(3)/), type_double,    &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', 'm', 'height'], 3, heightID)

      call RegisterVar(ncid, 'treesai', (/dimIDs(1), dimIDs(3)/), type_double,   &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', 'm2 m-2', 'total stem area index'], &
        3, treesaiID)

      call RegisterVar(ncid, 'nv', (/dimIDs(1), dimIDs(3)/), type_int,           &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', '-', 'number of occupied leaf+stem layers'], &
        3, nvID)

      call RegisterVar(ncid, 'repro_c', (/dimIDs(1), dimIDs(3)/), type_double,   &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', 'kgC indiv-1', 'reproductive carbon'], &
        3, reprocID)

      call RegisterVar(ncid, 'daily_net_c', (/dimIDs(1), dimIDs(3)/), type_double, &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', 'kgC indiv-1 day-1',          &
        'daily net carbon (GPP - leaf dark resp - nonleaf MR) - equal to daily_gpp - daily_rdark - daily_livestem_mr - daily_livecroot_mr - daily_froot_mr by construction'], &
        3, dailynetcID)

      call RegisterVar(ncid, 'daily_gpp', (/dimIDs(1), dimIDs(3)/), type_double, &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', 'kgC indiv-1 day-1', 'daily GPP'], &
        3, dailygppID)

      call RegisterVar(ncid, 'daily_rdark', (/dimIDs(1), dimIDs(3)/), type_double, &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', 'kgC indiv-1 day-1', 'daily leaf dark respiration'], &
        3, dailyrdarkID)

      call RegisterVar(ncid, 'daily_livestem_mr', (/dimIDs(1), dimIDs(3)/), type_double, &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', 'kgC indiv-1 day-1',          &
        'daily live stem (aboveground sapwood) maintenance respiration'], 3, dailylivestemmrID)

      call RegisterVar(ncid, 'daily_livecroot_mr', (/dimIDs(1), dimIDs(3)/), type_double, &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', 'kgC indiv-1 day-1',          &
        'daily live coarse root (belowground sapwood) maintenance respiration'], 3, dailylivecrootmrID)

      call RegisterVar(ncid, 'daily_froot_mr', (/dimIDs(1), dimIDs(3)/), type_double, &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', 'kgC indiv-1 day-1', 'daily fine root maintenance respiration'], &
        3, dailyfrootmrID)

      call RegisterVar(ncid, 'daily_growth_resp', (/dimIDs(1), dimIDs(3)/), type_double, &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', 'kgC indiv-1 day-1', 'daily growth respiration'], &
        3, dailygrowthrespID)

      call RegisterVar(ncid, 'leaf_turnover', (/dimIDs(1), dimIDs(3)/), type_double, &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', 'kgC indiv-1 day-1', 'daily leaf turnover loss'], &
        3, leafturnoverID)

      call RegisterVar(ncid, 'fnrt_turnover', (/dimIDs(1), dimIDs(3)/), type_double, &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', 'kgC indiv-1 day-1', 'daily fine root turnover loss'], &
        3, fnrtturnoverID)

      call RegisterVar(ncid, 'sapw_turnover', (/dimIDs(1), dimIDs(3)/), type_double, &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', 'kgC indiv-1 day-1', 'daily sapwood turnover loss'], &
        3, sapwturnoverID)

      call RegisterVar(ncid, 'struct_turnover', (/dimIDs(1), dimIDs(3)/), type_double, &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', 'kgC indiv-1 day-1', 'daily structural turnover loss'], &
        3, structturnoverID)

      call RegisterVar(ncid, 'npp_acc', (/dimIDs(1), dimIDs(3)/), type_double,   &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', 'kgC indiv-1 day-1',          &
        'carbon balance handed to PARTEH (net of growth respiration)'], 3, nppaccID)

      call RegisterVar(ncid, 'frac_store', (/dimIDs(1), dimIDs(3)/), type_double, &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', '-',                        &
        'storage carbon as a fraction of target leaf carbon'], 3, fracstoreID)

      call RegisterVar(ncid, 'cmort', (/dimIDs(1), dimIDs(3)/), type_double,     &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', 'indiv yr-1',                &
        'carbon starvation mortality rate'], 3, cmortID)

      call RegisterVar(ncid, 'maintresp_reduction_factor', (/dimIDs(1), dimIDs(3)/), type_double, &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'time light_level', '-',                        &
        'storage-based maintenance-respiration throttle'], 3, maintrespreductionfactorID)

      ! per-leaf-layer light profile, dimensioned (nlevleaf, year, light_level),
      ! with unoccupied layers above that year's nv filled with fates_unset_r8
      call RegisterVar(ncid, 'parsun_z', (/dimIDs(2), dimIDs(4), dimIDs(3)/), type_double, &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'nlevleaf year light_level', 'W m-2',            &
        'absorbed PAR, sunlit leaves, per unit ground area (first day of year, solar noon)'], &
        3, parsunID)
      call RegisterFillValue(ncid, parsunID, fates_unset_r8)

      call RegisterVar(ncid, 'parsha_z', (/dimIDs(2), dimIDs(4), dimIDs(3)/), type_double, &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'nlevleaf year light_level', 'W m-2',            &
        'absorbed PAR, shaded leaves, per unit ground area (first day of year, solar noon)'], &
        3, parshaID)
      call RegisterFillValue(ncid, parshaID, fates_unset_r8)

      call RegisterVar(ncid, 'laisun_z', (/dimIDs(2), dimIDs(4), dimIDs(3)/), type_double, &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'nlevleaf year light_level', 'm2 m-2',           &
        'sunlit LAI (first day of year, solar noon)'], 3, laisunID)
      call RegisterFillValue(ncid, laisunID, fates_unset_r8)

      call RegisterVar(ncid, 'laisha_z', (/dimIDs(2), dimIDs(4), dimIDs(3)/), type_double, &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
        [character(len=150) :: 'nlevleaf year light_level', 'm2 m-2',           &
        'shaded LAI (first day of year, solar noon)'], 3, laishaID)
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
