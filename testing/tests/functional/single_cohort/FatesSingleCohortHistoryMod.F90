module FatesSingleCohortHistoryMod
  !
  ! DESCRIPTION:
  ! Outputs for single_cohort test:
  !
  ! - Daily cohort time series, dimensioned (time, light_level), where time is the
  ! day index within each light level's trajectory
  ! - Daily environmental forcing the test was driven with, dimensioned (time)
  ! - Per-leaf-layer light profile (PAR sun/shade, LAI sun/shade), dimensioned 
  ! (nlevleaf, year, light_level). This is captured once per year
  ! - Instantaneous whole-plant light response diagnostics, dimensioned 
  ! (ppfd, year, light_level). Also captured once per year (same day as light profile)
  ! - Instantaneous leaf-level light-response diagnostics, dimensioned 
  ! (ppfd, year, light_level). Also captured once per year on the same day as above.
  ! This outputs a single canopy layer's (default top) net photosynthesis.
  ! 
  ! The test will optionally output ONLY a reduced_output set of variables.

  use FatesConstantsMod,       only : r8 => fates_r8
  use FatesConstantsMod,       only : fates_unset_r8
  use FatesCohortMod,          only : fates_cohort_type
  use PRTGenericMod,           only : leaf_organ, fnrt_organ, sapw_organ
  use PRTGenericMod,           only : struct_organ, store_organ
  use PRTGenericMod,           only : repro_organ
  use PRTGenericMod,           only : carbon12_element
  use FatesTestLightEnvMod,    only : light_env_type
  use FatesTestEnvironmentMod, only : environment_type
  use FatesUnitTestIOMod,      only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod,      only : WriteVar, RegisterVarAtts, RegisterFillValue
  use FatesUnitTestIOMod,      only : EndNCDef
  use FatesUnitTestIOMod,      only : type_double, type_int

  implicit none
  private

  type, public :: history_type

  private
  
  ! daily whole-cohort time series, dimensioned (time, light_level)
  real(r8), allocatable :: dbh(:,:)                        ! dbh [cm]
  real(r8), allocatable :: height(:,:)                     ! height [m]
  real(r8), allocatable :: treelai(:,:)                    ! total leaf area index [m2/m2]
  real(r8), allocatable :: treesai(:,:)                    ! total stem area index [m2/m2]
  real(r8), allocatable :: crown_area(:,:)                 ! crown area [m2]
  real(r8), allocatable :: leaf_c(:,:)                     ! leaf carbon [kgC/indiv]
  real(r8), allocatable :: fnrt_c(:,:)                     ! fine root carbon [kgC/indiv]
  real(r8), allocatable :: sapw_c(:,:)                     ! sapwood carbon [kgC/indiv]
  real(r8), allocatable :: struct_c(:,:)                   ! structural carbon [kgC/indiv]
  real(r8), allocatable :: storage_c(:,:)                  ! storage carbon [kgC/indiv]
  real(r8), allocatable :: repro_c(:,:)                    ! reproductive carbon [kgC/indiv]
  real(r8), allocatable :: daily_net_c(:,:)                ! daily net C (GPP - leaf Rd - nonleaf MR) [kgC/indiv/day]
  real(r8), allocatable :: daily_gpp(:,:)                  ! daily GPP [kgC/indiv/day]
  real(r8), allocatable :: daily_rdark(:,:)                ! daily leaf dark respiration [kgC/indiv/day]
  real(r8), allocatable :: daily_livestem_mr(:,:)          ! daily live stem (aboveground sapwood) maintenance respiration [kgC/indiv/day]
  real(r8), allocatable :: daily_livecroot_mr(:,:)         ! daily live coarse root (belowground sapwood) maintenance respiration [kgC/indiv/day]
  real(r8), allocatable :: daily_froot_mr(:,:)             ! daily fine root maintenance respiration [kgC/indiv/day]
  real(r8), allocatable :: daily_growth_resp(:,:)          ! daily growth respiration [kgC/indiv/day]
  real(r8), allocatable :: leaf_turnover(:,:)              ! daily leaf turnover loss [kgC/indiv/day]
  real(r8), allocatable :: fnrt_turnover(:,:)              ! daily fine root turnover loss [kgC/indiv/day]
  real(r8), allocatable :: sapw_turnover(:,:)              ! daily sapwood turnover loss [kgC/indiv/day]
  real(r8), allocatable :: struct_turnover(:,:)            ! daily structural turnover loss [kgC/indiv/day]
  real(r8), allocatable :: npp_acc(:,:)                    ! carbon balance handed to PARTEH (npp_acc, net of growth respiration) [kgC/indiv/day]
  real(r8), allocatable :: frac_store(:,:)                 ! storage carbon as a fraction of target leaf carbon [-]
  real(r8), allocatable :: cmort(:,:)                      ! carbon starvation mortality rate [indiv/year]
  real(r8), allocatable :: canopy_trim(:,:)                ! fraction of the maximum leaf biomass targeted [0-1]
  real(r8), allocatable :: leaf_cost(:,:)                  ! cost of maintaining leaves, bottom leaf layer [kgC/m2 leaf/year]
  real(r8), allocatable :: n(:,:)                          ! cohort number density, surviving fraction of the original recruitment cohort [0-1]
  real(r8), allocatable :: light_intercept_eff(:,:)        ! light interception efficiency [-]
  real(r8), allocatable :: maintresp_reduction_factor(:,:) ! storage-based maintenance-respiration throttle [0-1]
  real(r8), allocatable :: daily_absorbed_par_area(:,:)    ! whole-plant absorbed PAR per unit crown footprint, integrated over the day [J m-2 crown footprint day-1]
  real(r8), allocatable :: daily_absorbed_par_indiv(:,:)   ! whole-plant absorbed PAR per individual, integrated over the day [J indiv-1 day-1]
  real(r8), allocatable :: daily_incident_par(:,:)         ! incident PAR at the top of the crown, integrated over the day [J m-2 crown footprint day-1]
  integer,  allocatable :: nv(:,:)                         ! number of occupied leaf+stem layers [-]

  ! daily environmental forcing, dimensioned (time), identical across light levels
  real(r8), allocatable :: daily_temp(:)           ! daily-mean vegetation temperature [K]
  real(r8), allocatable :: daily_veg_esat(:)       ! daily-mean saturation vapor pressure [Pa]
  real(r8), allocatable :: daily_can_vpress(:)     ! daily-mean canopy air vapor pressure [Pa]
  real(r8), allocatable :: midday_temp(:)          ! vegetation temperature at the substep nearest solar noon [K]
  real(r8), allocatable :: midday_veg_esat(:)      ! saturation vapor pressure at that substep [Pa]
  real(r8), allocatable :: midday_can_vpress(:)    ! canopy air vapor pressure at that substep [Pa]
  real(r8), allocatable :: t_growth(:)             ! 10-day running-mean growth temperature [K]
  real(r8), allocatable :: t_home(:)               ! long-term running-mean home temperature [K]
  integer,  allocatable :: n_vpress_constrained(:) ! substeps this day at which the requested vapor pressure hit a bound [-]

  ! whole-trajectory Ci-solver summary, dimensioned (light_level)
  real(r8), allocatable :: mean_ci_solve_iter(:)    ! mean Ci-solver iteration count over every LeafLayerPhotosynthesis call in this trajectory [-]
  integer,  allocatable :: n_bisection_fallbacks(:) ! count of those calls that fell back to CiBisection [-]
  
  ! termination year and reason
  integer,  allocatable :: term_year(:)   ! year the cohort met a termination reason, 0 = never [-]
  integer,  allocatable :: term_reason(:) ! which reason: 0 none, 1 storage, 2 live pools, 3 negative biomass, 4 number density [-]

  ! per-leaf-layer light profile, dimensioned (nlevleaf, year, light_level)
  real(r8), allocatable :: parsun_z(:,:,:) ! absorbed PAR, sunlit leaves [W/m2 crown footprint]
  real(r8), allocatable :: parsha_z(:,:,:) ! absorbed PAR, shaded leaves [W/m2 crown footprint]
  real(r8), allocatable :: laisun_z(:,:,:) ! sunlit LAI [m2/m2]
  real(r8), allocatable :: laisha_z(:,:,:) ! shaded LAI [m2/m2]

  ! instantaneous light-response diagnostic, dimensioned (ppfd, year, light_level)
  real(r8), allocatable :: gross_assim(:,:,:) ! whole-plant gross assimilation [kgC/indiv/s]
  real(r8), allocatable :: leaf_resp(:,:,:)  ! whole-plant leaf dark respiration [kgC/indiv/s]
  real(r8), allocatable :: total_resp(:,:,:)  ! whole-plant total respiration (leaf dark + non-leaf maintenance) [kgC/indiv/s]
  real(r8), allocatable :: leaf_anet(:,:,:)   ! leaf-level net photosynthesis (Aarea), single canopy layer, no canopy attenuation (Sterck et al. 2013's LCPleaf diagnostic) [umolC/m2 leaf/s]
  
  ! dimensions
  integer :: n_time         ! days per light level's trajectory (n_days_total)
  integer :: n_leaf         ! compile-time max leaf+stem layers (nlevleaf)
  integer :: n_light        ! number of light levels
  integer :: n_year         ! number of simulated years
  integer :: n_ppfd         ! number of swept PPFD values in the light-response diagnostic
  logical :: reduced_output ! if true, only write reduced output

  contains

    procedure, public :: Init
    procedure, public :: RecordDay
    procedure, public :: RecordLightProfile
    procedure, public :: RecordLightResponse
    procedure, public :: RecordLeafNetAssim
    procedure, public :: RecordLightLevelSummary
    procedure, public :: RecordTermination
    procedure, public :: WriteVals

  end type history_type

contains

  ! ==========================================================================

  subroutine Init(this, n_time, n_leaf, n_light, n_year, n_ppfd, reduced_output)
    !
    ! DESCRIPTION:
    ! Allocate the time-series arrays and pre-fill the light-profile arrays with
    ! the fill value. Some arrays are not allocated depending whethere we are using 
    ! reduced_ouutput mode

    ! ARGUMENTS:
    class(history_type),  intent(inout) :: this    ! history object
    integer,              intent(in)    :: n_time  ! days per light level's trajectory
    integer,              intent(in)    :: n_leaf  ! compile-time max leaf+stem layers (nlevleaf)
    integer,              intent(in)    :: n_light ! number of light levels
    integer,              intent(in)    :: n_year  ! number of simulated years
    integer,              intent(in)    :: n_ppfd  ! number of swept PPFD values in the light-response diagnostic
    logical,              intent(in), optional :: reduced_output ! if true, skip the non-essential output groups (default .false.: full output)

    this%n_time = n_time
    this%n_leaf = n_leaf
    this%n_light = n_light
    this%n_year = n_year
    this%n_ppfd = n_ppfd

    this%reduced_output = .false.
    if (present(reduced_output)) this%reduced_output = reduced_output

    ! shade-tolerance metric set
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
    allocate(this%daily_absorbed_par_area(n_time, n_light))
    allocate(this%daily_absorbed_par_indiv(n_time, n_light))
    allocate(this%daily_incident_par(n_time, n_light))
    allocate(this%frac_store(n_time, n_light))
    allocate(this%cmort(n_time, n_light))
    allocate(this%canopy_trim(n_time, n_light))
    allocate(this%leaf_cost(n_time, n_light))
    allocate(this%gross_assim(n_ppfd, n_year, n_light))
    allocate(this%leaf_resp(n_ppfd, n_year, n_light))
    allocate(this%total_resp(n_ppfd, n_year, n_light))
    allocate(this%leaf_anet(n_ppfd, n_year, n_light))

    ! cheap (light_level only) numerical-health check
    allocate(this%mean_ci_solve_iter(n_light))
    allocate(this%n_bisection_fallbacks(n_light))
    
    allocate(this%term_year(n_light))
    allocate(this%term_reason(n_light))
    this%term_year(:) = 0
    this%term_reason(:) = 0

    ! daily environmental forcing (time only)
    allocate(this%daily_temp(n_time))
    allocate(this%daily_veg_esat(n_time))
    allocate(this%daily_can_vpress(n_time))
    allocate(this%midday_temp(n_time))
    allocate(this%midday_veg_esat(n_time))
    allocate(this%midday_can_vpress(n_time))
    allocate(this%t_growth(n_time))
    allocate(this%t_home(n_time))
    allocate(this%n_vpress_constrained(n_time))
    
    ! fill with unset
    this%dbh(:,:) = fates_unset_r8
    this%treelai(:,:) = fates_unset_r8
    this%crown_area(:,:) = fates_unset_r8
    this%leaf_c(:,:) = fates_unset_r8
    this%fnrt_c(:,:) = fates_unset_r8
    this%sapw_c(:,:) = fates_unset_r8
    this%struct_c(:,:) = fates_unset_r8
    this%storage_c(:,:) = fates_unset_r8
    this%n(:,:) = fates_unset_r8
    this%light_intercept_eff(:,:) = fates_unset_r8
    this%daily_absorbed_par_area(:,:) = fates_unset_r8
    this%daily_absorbed_par_indiv(:,:) = fates_unset_r8
    this%daily_incident_par(:,:) = fates_unset_r8
    this%frac_store(:,:) = fates_unset_r8
    this%cmort(:,:) = fates_unset_r8
    this%canopy_trim(:,:) = fates_unset_r8
    this%leaf_cost(:,:) = fates_unset_r8
    this%gross_assim(:,:,:) = fates_unset_r8
    this%leaf_resp(:,:,:) = fates_unset_r8
    this%total_resp(:,:,:) = fates_unset_r8
    this%leaf_anet(:,:,:) = fates_unset_r8
    this%mean_ci_solve_iter(:) = fates_unset_r8
    this%n_bisection_fallbacks(:) = 0
    
    this%daily_temp(:) = fates_unset_r8
    this%daily_veg_esat(:) = fates_unset_r8
    this%daily_can_vpress(:) = fates_unset_r8
    this%midday_temp(:) = fates_unset_r8
    this%midday_veg_esat(:) = fates_unset_r8
    this%midday_can_vpress(:) = fates_unset_r8
    this%t_growth(:) = fates_unset_r8
    this%t_home(:) = fates_unset_r8
    this%n_vpress_constrained(:) = 0
    

    ! only allocate in reduced_output mode
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
      allocate(this%maintresp_reduction_factor(n_time, n_light))
      allocate(this%parsun_z(n_leaf, n_year, n_light))
      allocate(this%parsha_z(n_leaf, n_year, n_light))
      allocate(this%laisun_z(n_leaf, n_year, n_light))
      allocate(this%laisha_z(n_leaf, n_year, n_light))
           
      this%height(:, :) = fates_unset_r8
      this%treesai(:, :) = fates_unset_r8
      this%nv(:, :) = 0
      this%repro_c(:, :) = fates_unset_r8
      this%daily_net_c(:, :) = fates_unset_r8
      this%daily_gpp(:, :) = fates_unset_r8
      this%daily_rdark(:, :) = fates_unset_r8
      this%daily_livestem_mr(:, :) = fates_unset_r8
      this%daily_livecroot_mr(:, :) = fates_unset_r8
      this%daily_froot_mr(:, :) = fates_unset_r8
      this%daily_growth_resp(:, :) = fates_unset_r8
      this%leaf_turnover(:, :) = fates_unset_r8
      this%fnrt_turnover(:, :) = fates_unset_r8
      this%sapw_turnover(:, :) = fates_unset_r8
      this%struct_turnover(:, :) = fates_unset_r8
      this%npp_acc(:, :) = fates_unset_r8
      this%maintresp_reduction_factor(:, :) = fates_unset_r8
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
    maintresp_reduction_factor, daily_incident_par, daily_absorbed_par_area,   &
    daily_absorbed_par_indiv, env)
    !
    ! DESCRIPTION:
    ! Capture one day's daily time series, optionally skip some if we are doing
    ! reduced_output

    ! ARGUMENTS:
    class(history_type),     intent(inout) :: this                       ! history object
    type(fates_cohort_type), intent(in)    :: cohort                     ! cohort to record
    type(environment_type),  intent(in)    :: env                        ! environment holding today's finalized forcing diagnostics
    integer,                 intent(in)    :: iday_all                   ! day index within this light level's trajectory (1..n_time)
    integer,                 intent(in)    :: ilight                     ! light-level index
    real(r8),                intent(in)    :: daily_net_c                ! daily net C (GPP - leaf Rd - nonleaf MR) [kgC/indiv/day]
    real(r8),                intent(in)    :: daily_gpp                  ! daily GPP [kgC/indiv/day]
    real(r8),                intent(in)    :: daily_rdark                ! daily leaf dark respiration [kgC/indiv/day]
    real(r8),                intent(in)    :: daily_livestem_mr          ! daily live stem maintenance respiration [kgC/indiv/day]
    real(r8),                intent(in)    :: daily_livecroot_mr         ! daily live coarse root maintenance respiration [kgC/indiv/day]
    real(r8),                intent(in)    :: daily_froot_mr             ! daily fine root maintenance respiration [kgC/indiv/day]
    real(r8),                intent(in)    :: daily_growth_resp          ! daily growth respiration [kgC/indiv/day]
    real(r8),                intent(in)    :: leaf_turnover              ! daily leaf turnover loss [kgC/indiv/day]
    real(r8),                intent(in)    :: fnrt_turnover              ! daily fine root turnover loss [kgC/indiv/day]
    real(r8),                intent(in)    :: sapw_turnover              ! daily sapwood turnover loss [kgC/indiv/day]
    real(r8),                intent(in)    :: struct_turnover            ! daily structural turnover loss [kgC/indiv/day]
    real(r8),                intent(in)    :: npp_acc                    ! carbon balance handed to PARTEH (net of growth respiration) [kgC/indiv/day]
    real(r8),                intent(in)    :: frac_store                 ! storage carbon as a fraction of target leaf carbon [-]
    real(r8),                intent(in)    :: cmort                      ! carbon starvation mortality rate [indiv/year]
    real(r8),                intent(in)    :: light_intercept_eff        ! light interception efficiency [-]
    real(r8),                intent(in)    :: maintresp_reduction_factor ! storage-based maintenance-respiration throttle [0-1]
    real(r8),                intent(in)    :: daily_incident_par         ! incident PAR [J m-2 crown footprint day-1] 
    real(r8),                intent(in)    :: daily_absorbed_par_area    ! whole-plant absorbed PAR per unit crown footprint [J m-2 crown footprint day-1]
    real(r8),                intent(in)    :: daily_absorbed_par_indiv   ! whole-plant absorbed PAR per individual [J indiv-1 day-1]

    ! daily environmental forcing
    this%daily_temp(iday_all) = env%daily_temp
    this%daily_veg_esat(iday_all) = env%daily_veg_esat
    this%daily_can_vpress(iday_all) = env%daily_can_vpress
    this%midday_temp(iday_all) = env%midday_temp
    this%midday_veg_esat(iday_all) = env%midday_veg_esat
    this%midday_can_vpress(iday_all) = env%midday_can_vpress
    this%t_growth(iday_all) = env%t_growth
    this%t_home(iday_all) = env%t_home
    this%n_vpress_constrained(iday_all) = env%n_vpress_constrained

    ! shade-tolerance metric set
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
    this%daily_incident_par(iday_all, ilight) = daily_incident_par
    this%daily_absorbed_par_area(iday_all, ilight) = daily_absorbed_par_area
    this%daily_absorbed_par_indiv(iday_all, ilight) = daily_absorbed_par_indiv
    this%frac_store(iday_all, ilight) = frac_store
    this%cmort(iday_all, ilight) = cmort
    this%canopy_trim(iday_all, ilight) = cohort%canopy_trim
    this%leaf_cost(iday_all, ilight) = cohort%leaf_cost

    ! only save if we are not doing reduced_output
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
      this%maintresp_reduction_factor(iday_all, ilight) = maintresp_reduction_factor
    end if

  end subroutine RecordDay

  ! ==========================================================================

  subroutine RecordLightProfile(this, iyear, ilight, cohort, light_env)
    !
    ! DESCRIPTION:
    ! Capture the per-leaf-layer light profile snapshot for the given year and
    ! light level. cohort%nv <= n_leaf always (n_leaf is a compile-time maximum)
    ! Returns if we are doing reduced_output

    ! ARGUMENTS:
    class(history_type),     intent(inout) :: this      ! history object
    type(fates_cohort_type), intent(in)    :: cohort    ! cohort this light profile belongs to
    type(light_env_type),    intent(in)    :: light_env ! light environment holding this substep's profile
    integer,                 intent(in)    :: iyear     ! simulated year index (1..n_year)
    integer,                 intent(in)    :: ilight    ! light-level index

    if (this%reduced_output) return

    this%parsun_z(1:cohort%nv, iyear, ilight) = light_env%parsun_z(:)
    this%parsha_z(1:cohort%nv, iyear, ilight) = light_env%parsha_z(:)
    this%laisun_z(1:cohort%nv, iyear, ilight) = light_env%laisun_z(:)
    this%laisha_z(1:cohort%nv, iyear, ilight) = light_env%laisha_z(:)

  end subroutine RecordLightProfile

  ! ==========================================================================

  subroutine RecordLightResponse(this, iyear, ilight, gross_assim, leaf_resp, total_resp)
    !
    ! DESCRIPTION:
    ! Capture the instantaneous whole-plant light-response diagnostic for the given year 
    ! and light level.

    ! ARGUMENTS:
    class(history_type), intent(inout) :: this           ! history object
    integer,             intent(in)    :: iyear          ! simulated year index (1..n_year)
    integer,             intent(in)    :: ilight         ! light-level index
    real(r8),            intent(in)    :: gross_assim(:) ! whole-plant gross assimilation at each swept PPFD [kgC/indiv/s], size(n_ppfd)
    real(r8),            intent(in)    :: leaf_resp(:)   ! whole-plant leaf dark respiration at each swept PPFD [kgC/indiv/s], size(n_ppfd)
    real(r8),            intent(in)    :: total_resp(:)  ! whole-plant total respiration at each swept PPFD [kgC/indiv/s], size(n_ppfd)

    this%gross_assim(:, iyear, ilight) = gross_assim(:)
    this%leaf_resp(:, iyear, ilight) = leaf_resp(:)
    this%total_resp(:, iyear, ilight) = total_resp(:)

  end subroutine RecordLightResponse

  ! ==========================================================================

  subroutine RecordLeafNetAssim(this, iyear, ilight, leaf_anet)
    !
    ! DESCRIPTION:
    ! Capture the instantaneous leaf-level light-response diagnostic for the given year 
    ! and light level.

    ! ARGUMENTS:
    class(history_type), intent(inout) :: this         ! history object
    real(r8),            intent(in)    :: leaf_anet(:) ! leaf-level net photosynthesis at each swept PPFD [umolC/m2 leaf/s], size(n_ppfd)
    integer,             intent(in)    :: iyear        ! simulated year index (1..n_year)
    integer,             intent(in)    :: ilight       ! light-level index
    
    this%leaf_anet(:, iyear, ilight) = leaf_anet(:)

  end subroutine RecordLeafNetAssim

  ! ==========================================================================

  subroutine RecordLightLevelSummary(this, ilight, mean_solve_iter, n_bisection_calls)
    !
    ! DESCRIPTION:
    ! Capture the whole-trajectory Ci-solver summary for the given light level
    ! a single value per light level, not a daily quantity.

    ! ARGUMENTS:
    class(history_type), intent(inout) :: this              ! history object
    integer,             intent(in)    :: ilight            ! light-level index
    real(r8),            intent(in)    :: mean_solve_iter   ! mean Ci-solver iteration count over this trajectory [-]
    integer,             intent(in)    :: n_bisection_calls ! count of calls that fell back to CiBisection over this trajectory [-]

    this%mean_ci_solve_iter(ilight) = mean_solve_iter
    this%n_bisection_fallbacks(ilight) = n_bisection_calls

  end subroutine RecordLightLevelSummary
  
  ! ==========================================================================
  
  subroutine RecordTermination(this, ilight, term_year, term_reason)
    !
    ! DESCRIPTION:
    ! Record whether and why this light level's cohort was terminated
    !

    ! ARGUMENTS:
    class(history_type), intent(inout) :: this        ! history object
    integer,             intent(in)    :: ilight      ! light-level index
    integer,             intent(in)    :: term_year   ! year the cohort terminated, 0 = never [-]
    integer,             intent(in)    :: term_reason ! termination criterion met, 0 = none [-]

    this%term_year(ilight) = term_year
    this%term_reason(ilight) = term_reason

  end subroutine RecordTermination

  ! ==========================================================================

  subroutine WriteVals(this, out_file, light_frac, ppfd_values)
    !
    ! DESCRIPTION:
    ! Write out the data

    ! ARGUMENTS:
    class(history_type), intent(in) :: this           ! history object
    character(len=*),    intent(in) :: out_file       ! output file name
    real(r8),            intent(in) :: light_frac(:)  ! swept incident light fractions [0-1]
    real(r8),            intent(in) :: ppfd_values(:) ! swept PPFD values in the light-response diagnostic [umol/m2/s]

    ! LOCALS:
    integer, allocatable :: time_idx(:)   ! day index, 1 = first simulated day
    integer, allocatable :: leaf_layer(:) ! leaf-layer index, 1 = top of crown
    integer, allocatable :: year_idx(:)   ! year index, 1 = first simulated year
    integer              :: i             ! looping index
    integer              :: ncid          ! netcdf file id
    character(len=20)    :: dim_names(5)  ! dimension names
    integer              :: dimIDs(5)     ! dimension IDs
    
    ! variable IDs
    integer              :: timeID, leaflayerID, lightfracID, yearID, ppfdID
    integer              :: dbhID, treelaiID, crownareaID
    integer              :: leafcID, fnrtcID, sapwcID, structcID, storagecID
    integer              :: nID, lightintercepteffID, dailyabsorbedparindivID
    integer              :: dailyincidentparID, dailyabsorbedparareaID
    integer              :: meancisolveiterID, nbisectionfallbacksID
    integer              :: dailytempID, dailyvegesatID, dailycanvpressID
    integer              :: middaytempID, middayvegesatID, middaycanvpressID
    integer              :: tgrowthID, thomeID, nvpressconstrainedID
    integer              :: grossassimID, totalrespID, leafrespID, leafanetID
    integer              :: heightID, treesaiID, nvID, reprocID
    integer              :: dailynetcID, dailygppID, nppaccID, fracstoreID
    integer              :: dailyrdarkID, dailylivestemmrID, dailylivecrootmrID
    integer              :: dailyfrootmrID, dailygrowthrespID
    integer              :: leafturnoverID, fnrtturnoverID, sapwturnoverID, structturnoverID
    integer              :: cmortID, maintrespreductionfactorID
    integer              :: canopytrimID, leafcostID
    integer              :: parsunID, parshaID, laisunID, laishaID
    integer              :: termyearID, termreasonID

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

    call RegisterVarAtts(ncid, 'dbh', (/dimIDs(1), dimIDs(3)/), type_double, 'cm',      &
      'dbh', dbhID, coordinates='time light_level')
    call RegisterFillValue(ncid, dbhID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'treelai', (/dimIDs(1), dimIDs(3)/), type_double,        &
      'm2 m-2', 'total leaf area index', treelaiID, coordinates='time light_level')
    call RegisterFillValue(ncid, treelaiID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'crown_area', (/dimIDs(1), dimIDs(3)/), type_double,     &
      'm2', 'crown area', crownareaID, coordinates='time light_level')
    call RegisterFillValue(ncid, crownareaID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'leaf_c', (/dimIDs(1), dimIDs(3)/), type_double,         &
      'kgC indiv-1', 'leaf carbon', leafcID, coordinates='time light_level')
    call RegisterFillValue(ncid, leafcID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'fnrt_c', (/dimIDs(1), dimIDs(3)/), type_double,         &
      'kgC indiv-1', 'fine root carbon', fnrtcID, coordinates='time light_level')
    call RegisterFillValue(ncid, fnrtcID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'sapw_c', (/dimIDs(1), dimIDs(3)/), type_double,         &
      'kgC indiv-1', 'sapwood carbon', sapwcID, coordinates='time light_level')
    call RegisterFillValue(ncid, sapwcID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'struct_c', (/dimIDs(1), dimIDs(3)/), type_double,       &
      'kgC indiv-1', 'structural carbon', structcID, coordinates='time light_level')
    call RegisterFillValue(ncid, structcID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'storage_c', (/dimIDs(1), dimIDs(3)/), type_double,      &
      'kgC indiv-1', 'storage carbon', storagecID, coordinates='time light_level')
    call RegisterFillValue(ncid, storagecID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'n', (/dimIDs(1), dimIDs(3)/), type_double, 'indiv',     &
      'cohort number density (surviving fraction of the original recruitment cohort)',  &
      nID, coordinates='time light_level')
    call RegisterFillValue(ncid, nID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'light_intercept_eff', (/dimIDs(1), dimIDs(3)/),         &
      type_double, '-',                                                                 &
      'light interception efficiency: whole-plant absorbed PAR / absorbed PAR of an equal-leaf-area, zero-self-shading reference surface, energy-weighted over the day (Sterck et al. 2013)', &
      lightintercepteffID, coordinates='time light_level')
    call RegisterFillValue(ncid, lightintercepteffID, fates_unset_r8)
      
    call RegisterVarAtts(ncid, 'daily_incident_par', (/dimIDs(1), dimIDs(3)/),    &
      type_double, 'J m-2 crown footprint day-1',                                 &
      'incident PAR at the top of the crown, integrated over the day', &
      dailyincidentparID, coordinates='time light_level')
    call RegisterFillValue(ncid, dailyincidentparID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'daily_absorbed_par_area', (/dimIDs(1), dimIDs(3)/),     &
      type_double, 'J m-2 crown footprint day-1',                                       &
      'absorbed PAR per unit crown footprint, integrated over the day - same basis as daily_incident_par, so their ratio is the absorbed fraction', &
      dailyabsorbedparareaID, coordinates='time light_level')
    call RegisterFillValue(ncid, dailyabsorbedparareaID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'daily_absorbed_par_indiv', (/dimIDs(1), dimIDs(3)/),    &
      type_double, 'J indiv-1 day-1',                                                   &
      'whole-plant absorbed PAR per individual (Onoda et al. 2013''s Phi) - divide by treelai*crown_area/n for Phi/LA (LIE_LA), by a chosen biomass-pool total for Phi/M (LIE_M), or use its running sum against a biomass increment for LUE', &
      dailyabsorbedparindivID, coordinates='time light_level')
    call RegisterFillValue(ncid, dailyabsorbedparindivID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'daily_temp', dimIDs(1:1), type_double, 'K',             &
      'daily-mean vegetation temperature (mean over the day''s sub-daily substeps)',    &
      dailytempID, coordinates='time')
    call RegisterFillValue(ncid, dailytempID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'daily_veg_esat', dimIDs(1:1), type_double, 'Pa',        &
      'daily-mean saturation vapor pressure at the vegetation temperature',             &
      dailyvegesatID, coordinates='time')
    call RegisterFillValue(ncid, dailyvegesatID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'daily_can_vpress', dimIDs(1:1), type_double, 'Pa',      &
      'daily-mean canopy air vapor pressure - with daily_veg_esat gives daily-mean VPD (veg_esat - can_vpress) and RH (can_vpress/veg_esat)', &
      dailycanvpressID, coordinates='time')
    call RegisterFillValue(ncid, dailycanvpressID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'midday_temp', dimIDs(1:1), type_double, 'K',            &
      'vegetation temperature at the sub-daily substep nearest solar noon',             &
      middaytempID, coordinates='time')
    call RegisterFillValue(ncid, middaytempID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'midday_veg_esat', dimIDs(1:1), type_double, 'Pa',       &
      'saturation vapor pressure at the substep nearest solar noon',                    &
      middayvegesatID, coordinates='time')
    call RegisterFillValue(ncid, middayvegesatID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'midday_can_vpress', dimIDs(1:1), type_double, 'Pa',     &
      'canopy air vapor pressure at the substep nearest solar noon - with midday_veg_esat gives the midday VPD that drives stomatal conductance when PAR is highest', &
      middaycanvpressID, coordinates='time')
    call RegisterFillValue(ncid, middaycanvpressID, fates_unset_r8)

    call RegisterVarAtts(ncid, 't_growth', dimIDs(1:1), type_double, 'K',               &
      '10-day running-mean growth temperature (photosynthetic acclimation boundary condition)', &
      tgrowthID, coordinates='time')
    call RegisterFillValue(ncid, tgrowthID, fates_unset_r8)

    call RegisterVarAtts(ncid, 't_home', dimIDs(1:1), type_double, 'K',                 &
      'long-term running-mean home temperature (photosynthetic acclimation boundary condition)', &
      thomeID, coordinates='time')
    call RegisterFillValue(ncid, thomeID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'n_vpress_constrained', dimIDs(1:1), type_int, '-',      &
      'number of sub-daily substeps this day at which the prescribed canopy vapor pressure was altered by the saturation/minimum constraint - nonzero means the run was not forced with the humidity boundary condition it asked for', &
      nvpressconstrainedID, coordinates='time')

    call RegisterVarAtts(ncid, 'mean_ci_solve_iter', dimIDs(3:3), type_double, '-',     &
      'mean Ci-solver iteration count over every LeafLayerPhotosynthesis call in this trajectory', &
      meancisolveiterID, coordinates='light_level')
    call RegisterFillValue(ncid, meancisolveiterID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'n_bisection_fallbacks', dimIDs(3:3), type_int, '-',     &
      'count of LeafLayerPhotosynthesis calls that fell back to CiBisection in this trajectory', &
      nbisectionfallbacksID, coordinates='light_level')
      
    call RegisterVarAtts(ncid, 'term_year', dimIDs(3:3), type_int, '-',                &
      'year the cohort met a termination criterion, 0 = ran the full trajectory',      &
      termyearID, coordinates='light_level')

    call RegisterVarAtts(ncid, 'term_reason', dimIDs(3:3), type_int, '-',              &
      'termination criterion met: 0 none, 1 storage depleted, 2 live pools depleted, ' // &
      '3 negative total biomass, 4 number density below numerical safety',             &
      termreasonID, coordinates='light_level')

    call RegisterVarAtts(ncid, 'gross_assim', (/dimIDs(5), dimIDs(4), dimIDs(3)/),      &
      type_double, 'kgC indiv-1 s-1',                                                   &
      'whole-plant gross assimilation at each swept PPFD (first day of year, pure diffuse illumination, coszen=1)', &
      grossassimID, coordinates='ppfd year light_level')
    call RegisterFillValue(ncid, grossassimID, fates_unset_r8)
      
    call RegisterVarAtts(ncid, 'leaf_resp', (/dimIDs(5), dimIDs(4), dimIDs(3)/),       &
      type_double, 'kgC indiv-1 s-1',                                                   &
      'whole-plant leaf dark respiration at each swept PPFD (first day of year)', &
      leafrespID, coordinates='ppfd year light_level')
    call RegisterFillValue(ncid, leafrespID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'total_resp', (/dimIDs(5), dimIDs(4), dimIDs(3)/),       &
      type_double, 'kgC indiv-1 s-1',                                                   &
      'whole-plant total respiration (leaf dark + non-leaf maintenance) at each swept PPFD (first day of year)', &
      totalrespID, coordinates='ppfd year light_level')
    call RegisterFillValue(ncid, totalrespID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'leaf_anet', (/dimIDs(5), dimIDs(4), dimIDs(3)/),        &
      type_double, 'umol m-2 s-1',                                                      &
      'leaf-level net photosynthesis (Aarea) at each swept PPFD, applied directly to a single canopy layer with no canopy attenuation or self-shading (first day of year) (Sterck et al. 2013)', &
      leafanetID, coordinates='ppfd year light_level')
    call RegisterFillValue(ncid, leafanetID, fates_unset_r8)
      
    call RegisterVarAtts(ncid, 'frac_store', (/dimIDs(1), dimIDs(3)/), type_double,   &
        '-', 'storage carbon as a fraction of target leaf carbon', fracstoreID,         &
        coordinates='time light_level')
    call RegisterFillValue(ncid, fracstoreID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'cmort', (/dimIDs(1), dimIDs(3)/), type_double,        &
        'indiv yr-1', 'carbon starvation mortality rate', cmortID,                      &
        coordinates='time light_level')
    call RegisterFillValue(ncid, cmortID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'canopy_trim', (/dimIDs(1), dimIDs(3)/), type_double,  &
        '-',                                                                          &
        'fraction of the maximum leaf biomass targeted - updated once a year by the canopy trim', &
        canopytrimID, coordinates='time light_level')
    call RegisterFillValue(ncid, canopytrimID, fates_unset_r8)

    call RegisterVarAtts(ncid, 'leaf_cost', (/dimIDs(1), dimIDs(3)/), type_double,    &
        'kgC m-2 leaf yr-1',                                                          &
        'cost of maintaining leaves in the bottom leaf layer, the quantity year_net_uptake is tested against by the canopy trim - held at the previous trim''s value between trims', &
        leafcostID, coordinates='time light_level')
    call RegisterFillValue(ncid, leafcostID, fates_unset_r8)

    if (.not. this%reduced_output) then

      call RegisterVarAtts(ncid, 'height', (/dimIDs(1), dimIDs(3)/), type_double, 'm',  &
        'height', heightID, coordinates='time light_level')
      call RegisterFillValue(ncid, heightID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'treesai', (/dimIDs(1), dimIDs(3)/), type_double,      &
        'm2 m-2', 'total stem area index', treesaiID, coordinates='time light_level')
      call RegisterFillValue(ncid, treesaiID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'nv', (/dimIDs(1), dimIDs(3)/), type_int, '-',         &
        'number of occupied leaf+stem layers', nvID, coordinates='time light_level')

      call RegisterVarAtts(ncid, 'repro_c', (/dimIDs(1), dimIDs(3)/), type_double,      &
        'kgC indiv-1', 'reproductive carbon', reprocID, coordinates='time light_level')
      call RegisterFillValue(ncid, reprocID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'daily_net_c', (/dimIDs(1), dimIDs(3)/), type_double,  &
        'kgC indiv-1 day-1',                                                            &
        'daily net carbon (GPP - leaf dark resp - nonleaf MR) - equal to daily_gpp - daily_rdark - daily_livestem_mr - daily_livecroot_mr - daily_froot_mr by construction', &
        dailynetcID, coordinates='time light_level')
      call RegisterFillValue(ncid, dailynetcID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'daily_gpp', (/dimIDs(1), dimIDs(3)/), type_double,    &
        'kgC indiv-1 day-1', 'daily GPP', dailygppID, coordinates='time light_level')
      call RegisterFillValue(ncid, dailygppID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'daily_rdark', (/dimIDs(1), dimIDs(3)/), type_double,  &
        'kgC indiv-1 day-1', 'daily leaf dark respiration', dailyrdarkID,               &
        coordinates='time light_level')
      call RegisterFillValue(ncid, dailyrdarkID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'daily_livestem_mr', (/dimIDs(1), dimIDs(3)/),         &
        type_double, 'kgC indiv-1 day-1',                                               &
        'daily live stem (aboveground sapwood) maintenance respiration',                &
        dailylivestemmrID, coordinates='time light_level')
      call RegisterFillValue(ncid, dailylivestemmrID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'daily_livecroot_mr', (/dimIDs(1), dimIDs(3)/),        &
        type_double, 'kgC indiv-1 day-1',                                               &
        'daily live coarse root (belowground sapwood) maintenance respiration',         &
        dailylivecrootmrID, coordinates='time light_level')
      call RegisterFillValue(ncid, dailylivecrootmrID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'daily_froot_mr', (/dimIDs(1), dimIDs(3)/),            &
        type_double, 'kgC indiv-1 day-1', 'daily fine root maintenance respiration',    &
        dailyfrootmrID, coordinates='time light_level')
      call RegisterFillValue(ncid, dailyfrootmrID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'daily_growth_resp', (/dimIDs(1), dimIDs(3)/),         &
        type_double, 'kgC indiv-1 day-1', 'daily growth respiration',                   &
        dailygrowthrespID, coordinates='time light_level')
      call RegisterFillValue(ncid, dailygrowthrespID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'leaf_turnover', (/dimIDs(1), dimIDs(3)/),             &
        type_double, 'kgC indiv-1 day-1', 'daily leaf turnover loss', leafturnoverID,   &
        coordinates='time light_level')
      call RegisterFillValue(ncid, leafturnoverID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'fnrt_turnover', (/dimIDs(1), dimIDs(3)/),             &
        type_double, 'kgC indiv-1 day-1', 'daily fine root turnover loss',              &
        fnrtturnoverID, coordinates='time light_level')
      call RegisterFillValue(ncid, fnrtturnoverID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'sapw_turnover', (/dimIDs(1), dimIDs(3)/),             &
        type_double, 'kgC indiv-1 day-1', 'daily sapwood turnover loss',                &
        sapwturnoverID, coordinates='time light_level')
      call RegisterFillValue(ncid, sapwturnoverID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'struct_turnover', (/dimIDs(1), dimIDs(3)/),           &
        type_double, 'kgC indiv-1 day-1', 'daily structural turnover loss',             &
        structturnoverID, coordinates='time light_level')
      call RegisterFillValue(ncid, structturnoverID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'npp_acc', (/dimIDs(1), dimIDs(3)/), type_double,      &
        'kgC indiv-1 day-1',                                                            &
        'carbon balance handed to PARTEH (net of growth respiration)', nppaccID,        &
        coordinates='time light_level')
      call RegisterFillValue(ncid, nppaccID, fates_unset_r8)

      call RegisterVarAtts(ncid, 'maintresp_reduction_factor',                          &
        (/dimIDs(1), dimIDs(3)/), type_double, '-',                                     &
        'storage-based maintenance-respiration throttle', maintrespreductionfactorID,   &
        coordinates='time light_level')
      call RegisterFillValue(ncid, maintrespreductionfactorID, fates_unset_r8)

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
    call WriteVar(ncid, dailyincidentparID, this%daily_incident_par(:,:))
    call WriteVar(ncid, dailyabsorbedparareaID, this%daily_absorbed_par_area(:,:))
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
    call WriteVar(ncid, termyearID, this%term_year(:))
    call WriteVar(ncid, termreasonID, this%term_reason(:))
    call WriteVar(ncid, ppfdID, ppfd_values(:))
    call WriteVar(ncid, grossassimID, this%gross_assim(:,:,:))
    call WriteVar(ncid, leafrespID, this%leaf_resp(:,:,:))
    call WriteVar(ncid, totalrespID, this%total_resp(:,:,:))
    call WriteVar(ncid, leafanetID, this%leaf_anet(:,:,:))
    call WriteVar(ncid, fracstoreID, this%frac_store(:,:))
    call WriteVar(ncid, cmortID, this%cmort(:,:))
    call WriteVar(ncid, canopytrimID, this%canopy_trim(:,:))
    call WriteVar(ncid, leafcostID, this%leaf_cost(:,:))

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
      call WriteVar(ncid, maintrespreductionfactorID, this%maintresp_reduction_factor(:,:))
      call WriteVar(ncid, parsunID, this%parsun_z(:,:,:))
      call WriteVar(ncid, parshaID, this%parsha_z(:,:,:))
      call WriteVar(ncid, laisunID, this%laisun_z(:,:,:))
      call WriteVar(ncid, laishaID, this%laisha_z(:,:,:))
    end if

    call CloseNCFile(ncid)

  end subroutine WriteVals

end module FatesSingleCohortHistoryMod
