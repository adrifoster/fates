program FatesSingleCohort
  !
  ! DESCRIPTION:
  ! A test of a single cohort across fractional light levels
  ! For each light level, a single cohort is created at recruitment size and simulated
  ! for a specified number of years, including radiation, photosynthesis, and daily
  ! allocation, phenology, and growth. Each light level is an independent trajectory. 
  !
  ! Only cold-deciduous phenology is available currently for non-evergreen trees (i.e.,
  ! no drought deciduous)
  !
  ! Only cabon starvation mortality is calculated and decreases cohort%n from the initial 
  ! 1.0 value

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesConstantsMod,           only : nearzero
  use FatesConstantsMod,           only : sec_per_day
  use FatesConstantsMod,           only : cstarvation_model_lin, cstarvation_model_exp
  use FatesConstantsMod,           only : itrue, ifalse
  use FatesConstantsMod,           only : wm2_to_umolm2s
  use FatesConstantsMod,           only : leaves_on, leaves_off
  use FatesConstantsMod,           only : ihard_season_decid
  use EDTypesMod,                  only : min_n_safemath
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
  use FatesTestSiteMod,            only : ReadSiteNamelist
  use FatesTestEnvironmentMod,     only : environment_type
  use FatesTestLightEnvMod,        only : light_env_type
  use FatesTestCohortPhysMod,      only : cohort_phys_type
  use FatesSingleCohortHistoryMod, only : history_type
  use FatesTestLeafPhotoMod,       only : LeafNitrogenContent
  use LeafBiophysicsMod,           only : lb_params
  use LeafBiophysicsMod,           only : FvCB1980, medlyn_model
  use LeafBiophysicsMod,           only : net_assim_model, photosynth_acclim_model_kumarathunge_etal_2019
  use LeafBiophysicsMod,           only : LowstorageMainRespReduction
  use EDParamsMod,                 only : dinc_vai
  use FatesRadiationMemMod,        only : ipar

  implicit none

  ! LOCALS:
  type(environment_type)        :: env                     ! prescribed atmospheric/soil boundary conditions 
  type(history_type)            :: hist                    ! output accumulator/writer, shared across the light sweep
  real(r8), allocatable         :: light_frac(:)           ! swept incident light fractions [0-1]
  real(r8), allocatable         :: diagnostic_ppfd(:)      ! swept PPFD values for LightResponseSweep [umol/m2/s]
  real(r8)                      :: dbh_recruit             ! recruitment-size dbh [cm]
  real(r8)                      :: lnc_top                 ! leaf N content at the canopy top [gN/m2 leaf]
  real(r8)                      :: par_absorptance         ! PFT fractional PAR absorptance, 1 - rhol(ipar) - taul(ipar) [-]
  integer                       :: ilight                  ! fractional light level looping index
  integer                       :: ippfd                   ! PPFD looping index
  integer                       :: n_photo_calls_total     !
  integer                       :: n_bisection_calls_total !
  integer                       :: max_solve_iter_total    ! whole run
  logical                       :: reduced_output          ! optional 2nd command-line arg
  character(len=:), allocatable :: param_file              ! input parameter file
  character(len=:), allocatable :: site_namelist_file      ! optional 4th command-line arg
  character(len=:), allocatable :: out_file                ! output file

  ! CONSTANTS:
  integer,                     parameter :: pft = 1                         ! plant functional type to simulate
  real(r8), dimension(nclmax), parameter :: can_tlai = 0.0_r8               ! canopy-layer LAI above the cohort
  real(r8),                    parameter :: patch_area = 1.0e4_r8           ! reference ground area the cohort occupies [m2]
  real(r8),                    parameter :: n_indiv = 1.0_r8                ! number of individuals in the cohort
  real(r8),                    parameter :: coh_age = 0.0_r8                ! cohort age
  real(r8),                    parameter :: site_spread = 0.0_r8            ! site spread index
  integer,                     parameter :: leaf_on_doy = 60                ! prescribed leaf on day of year (for cold deciduous PFTs)
  integer,                     parameter :: leaf_off_doy = 305.             ! prescribed leaf off day of year (for cold deciduous PFTs)
  integer,                     parameter :: n_light_levels = 25             ! number of incident light levels to sweep
  real(r8),                    parameter :: light_frac_min = 0.005_r8       ! lowest incident light fraction swept [fraction of full sun]
  real(r8),                    parameter :: light_frac_max = 1.0_r8         ! highest incident light fraction swept [fraction of full sun]
  integer,                     parameter :: n_ppfd_diagnostic = 40          ! number of PPFD values to sweep
  real(r8),                    parameter :: ppfd_diagnostic_min = 0.01_r8   ! lowest swept PPFD [umol/m2/s]
  real(r8),                    parameter :: ppfd_diagnostic_max = 2000.0_r8 ! highest swept PPFD [umol/m2/s] (~full sun)
  integer,                     parameter :: leaf_lcp_layer = 1              ! canopy layer swept for the leaf-level (LCPleaf) diagnostic
  integer,                     parameter :: days_per_year = 365             ! days per simulated year
  integer,                     parameter :: n_substeps_per_day = 48         ! sub-daily steps per day (half-hourly), must be even
  integer,                     parameter :: nyears = 50                     ! number of years to simulate
  
  ! timing calculations
  real(r8), parameter :: step_size = 86400.0_r8/n_substeps_per_day ! model time step [s]
  real(r8), parameter :: hrs_per_substep = 24.0_r8/n_substeps_per_day ! sub-daily step length [h]
  integer,  parameter :: n_days_total = nyears * days_per_year     ! total days per light level's trajectory
  integer,  parameter :: noon_substep = n_substeps_per_day/2 + 1   ! sub-daily index nearest solar noon (first sample after noon)
  
  ! cohort termination
  ! this driver does not kill the cohorts but we will just stop running them
  integer, parameter :: term_none        = 0  ! still alive
  integer, parameter :: term_storage     = 1  ! storage carbon terminally depleted
  integer, parameter :: term_live_pools  = 2  ! sapwood+leaf+fineroot terminally depleted
  integer, parameter :: term_neg_biomass = 3  ! total cohort biomass negative
  integer, parameter :: term_num_dens    = 4  ! number density below numerical safety

  ! read in parameter file name from command line
  param_file = command_line_arg(1)

  ! optional 2nd command-line argument: switches the history writer to reduced_output mode
  reduced_output = .false.
  if (command_argument_count() >= 2) then
    reduced_output = (trim(command_line_arg(2)) == 'reduced_output')
  end if

  ! output file name, depends on either input arg3 or if 'reduced_output' is used
  if (command_argument_count() >= 3) then
    out_file = trim(command_line_arg(3))
  else if (reduced_output) then
    out_file = 'single_cohort_out_reduced.nc'
  else
    out_file = 'single_cohort_out.nc'
  end if
  
  ! optional 4th command-line argument: a &fates_test_site namelist file for 
  ! overriding default environmental/climate conditions
  if (command_argument_count() >= 4) then
    site_namelist_file = trim(command_line_arg(4))
    call ReadSiteNamelist(site_namelist_file)
  end if

  ! read in parameter file
  call ReadParameters(param_file)

  ! initialize global PRT/allometry data needed by the cohort machinery
  call InitializeGlobals(step_size)

  ! initialize FATES logging (the two-stream module's debug/error paths write to this
  ! unit and it is otherwise never set in a standalone driver) and the two-stream
  ! radiation parameters (leaf/stem optical properties, per pft)
  call FatesGlobalsInit(6, .false.)
  call TransferRadParams()

  ! derive the organ_id -> parameter-file-index reverse lookup map
  call PRTDerivedParams()

  ! host-model-namelist-controlled leaf biophysics switches
  lb_params%electron_transport_model = FvCB1980 ! Farquhar-von Caemmerer-Berry (1980)
  lb_params%stomatal_model           = medlyn_model
  lb_params%stomatal_assim_model     = net_assim_model
  lb_params%photo_tempsens_model     = photosynth_acclim_model_kumarathunge_etal_2019

  ! host-model-namelist-controlled carbon-starvation mortality model
  hlm_mort_cstarvation_model = cstarvation_model_lin

  ! leaf N content - constant for the whole run
  lnc_top = LeafNitrogenContent(pft)

  ! PFT PAR absorptance
  ! used to normalize light_intercept_eff (see RunOneLightLevel) to a 
  ! zero-self-shading asymptote of 1.0
  par_absorptance = 1.0_r8 - EDPftvarcon_inst%rhol(pft,ipar) - EDPftvarcon_inst%taul(pft,ipar)

  ! recruitment-size initialization
  call h2d_allom(EDPftvarcon_inst%hgt_min(pft), pft, dbh_recruit)

  ! build the log-spaced incident light fractions to sweep
  allocate(light_frac(n_light_levels))
  do ilight = 1, n_light_levels
    light_frac(ilight) = light_frac_min * (light_frac_max/light_frac_min) **      &
      (real(ilight - 1, r8)/real(n_light_levels - 1, r8))
  end do

  ! build the log-spaced diagnostic PPFD sweep
  allocate(diagnostic_ppfd(n_ppfd_diagnostic))
  do ippfd = 1, n_ppfd_diagnostic
    diagnostic_ppfd(ippfd) = ppfd_diagnostic_min *   &
      (ppfd_diagnostic_max/ppfd_diagnostic_min) **   &
      (real(ippfd - 1, r8)/real(n_ppfd_diagnostic - 1, r8))
  end do

  ! initialize these to 0
  n_photo_calls_total = 0
  n_bisection_calls_total = 0
  max_solve_iter_total = 0

  ! main light-level sweep
  ! each level is an independent trajectory from recruitment size
  do ilight = 1, n_light_levels
    call RunOneLightLevel(light_frac(ilight), ilight)
  end do

  ! write out the daily whole-cohort time series and the annual light-profile
  ! snapshot, both across the light sweep
  call hist%WriteVals(out_file, light_frac, diagnostic_ppfd)

contains

  ! ==========================================================================

  subroutine RunOneLightLevel(light_frac_val, ilight)
    !
    ! DESCRIPTION:
    ! Run one light level's independent year/day/sub-daily trajectory, from a
    ! freshly built recruitment-size cohort, writing daily output and the annual
    ! light profile

    ! ARGUMENTS:
    real(r8), intent(in) :: light_frac_val ! this light level's incident light fraction [0-1]
    integer,  intent(in) :: ilight         ! light-level index

    ! LOCALS:
    type(fates_cohort_type), pointer :: cohort                                    ! cohort for this light level
    type(light_env_type)             :: light_env                                 ! prescribed light environment for this light level
    type(cohort_phys_type)           :: phys                                      ! per-leaf-layer photosynthetic capacity/leaf physics
    real(r8)                         :: diagnostic_gross_assim(n_ppfd_diagnostic) ! one year's gross-assimilation output for the light response curve [kgC/indiv/s]
    real(r8)                         :: diagnostic_leaf_resp(n_ppfd_diagnostic)   ! one year's leaf dark respiration output for the light response curve[kgC/indiv/s]
    real(r8)                         :: diagnostic_total_resp(n_ppfd_diagnostic)  ! one year's total-respiration output for the light response curve[kgC/indiv/s]
    real(r8)                         :: diagnostic_leaf_anet(n_ppfd_diagnostic)   ! one year's leaf-level net-photosynthesis output for the light response curve [umolC/m2 leaf/s]
    real(r8)                         :: hour_of_day                               ! hour of day at the midpoint of the current substep [0-24]
    real(r8)                         :: frac_store                                ! ratio of storage carbon to target_leaf_c [-]
    real(r8)                         :: npp_acc_to_prt                            ! carbon_balance passed to PARTEH via cohort%npp_acc, net of growth respiration
    real(r8)                         :: cmort                                     ! carbon starvation mortality rate [indiv/year]
    real(r8)                         :: gpp_tstep                                 ! one substep's GPP [kgC/indiv/s]
    real(r8)                         :: rdark_tstep                               ! one substep's leaf dark respiration [kgC/indiv/s]
    real(r8)                         :: nonleaf_mr_tstep                          ! one substep's non-leaf maintenance respiration [kgC/indiv/s]
    real(r8)                         :: daily_net_c                               ! GPP - leaf dark resp - non-leaf maintenance resp, integrated over the day [kgC/indiv/day]
    real(r8)                         :: daily_gpp                                 ! GPP alone, integrated over the day [kgC/indiv/day]
    real(r8)                         :: daily_rdark                               ! leaf dark respiration, integrated over the day [kgC/indiv/day]
    real(r8)                         :: daily_livestem_mr                         ! live stem maintenance resp, integrated over the day [kgC/indiv/day]
    real(r8)                         :: daily_livecroot_mr                        ! live coarse root maintenance resp, integrated over the day [kgC/indiv/day]
    real(r8)                         :: daily_froot_mr                            ! fine root maintenance resp, integrated over the day [kgC/indiv/day] 
    real(r8)                         :: growth_resp                               ! one day's growth respiration [kgC/indiv/day]
    real(r8)                         :: leaf_turnover                             ! one day's leaf turnover [kgC/indiv/day]
    real(r8)                         :: fnrt_turnover                             ! one day's fineroot turnover [kgC/indiv/day]
    real(r8)                         :: sapw_turnover                             ! one day's sapwood turnover [kgC/indiv/day]
    real(r8)                         :: struct_turnover                           ! one day's structural turnover [kgC/indiv/day]
    real(r8)                         :: storage_c                                 ! current storage carbon [kgC/indiv]
    real(r8)                         :: par_toc                                   ! current substep's incident PAR at the top of the crown [W/m2]
    real(r8)                         :: par_beam                                  ! current substep's direct PAR at the top of the crown [W/m2]
    real(r8)                         :: par_diff                                  ! current substep's diffuse PAR at the top of the crown [W/m2]
    real(r8)                         :: daily_absorbed_par                        ! whole-plant absorbed PAR (parsun_z+parsha_z summed over layers), per unit leaf area, integrated over the day [J/m2 leaf]
    real(r8)                         :: daily_incident_par                        ! incident PAR at the top of the crown, integrated over the day [J/m2 crown footprint] 
    real(r8)                         :: daily_absorbed_par_indiv                  ! whole-plant absorbed PAR per individual, integrated over the day [J/indiv/day] 
    real(r8)                         :: light_intercept_eff                       ! one day's light interception efficiency
    real(r8)                         :: maintresp_reduction_factor                ! one days's storage-based maintenance-respiration factor [0-1]
    real(r8)                         :: mean_solve_iter                           ! this light level's mean Ci-solver iteration count (sum_solve_iter/n_photo_calls), computed once at the end of this trajectory
    integer                          :: iyear, iday                               ! year/day looping indices
    integer                          :: iday_all                                  ! day index within this light level's trajectory (1..n_days_total)
    integer                          :: isubday                                   ! sub-daily looping index
    integer                          :: n_photo_calls                             ! Ci solver diagnostic
    integer                          :: n_bisection_calls                         ! Ci solver diagnostic
    integer                          :: max_solve_iter                            ! Ci solver diagnostic 
    integer                          :: sum_solve_iter                            ! Ci solver diagnostic
    integer                          :: term_reason                               ! reason a cohort was terminated
    integer                          :: term_year                                 ! year a cohort was terminated
    
    
    n_photo_calls = 0
    n_bisection_calls = 0
    max_solve_iter = 0
    sum_solve_iter = 0

    ! create a new cohort at recruitment size
    allocate(cohort)
    call CohortFactory(cohort, pft, can_tlai, dbh=dbh_recruit, number=n_indiv,      &
      patch_area=patch_area, age=coh_age, site_spread=site_spread)
    cohort%nv = GetNVegLayers(cohort%treelai + cohort%treesai)

    ! initialize the light environment
    call light_env%Init(cohort%treelai, cohort%treesai, cohort%height, pft)

    ! prescribed atmospheric/soil boundary conditions
    call env%Init()

    ! allocate the history arrays once, from the first light level
    if (ilight == 1) then
      call hist%Init(n_days_total, nlevleaf, n_light_levels, nyears, n_ppfd_diagnostic, &
        reduced_output)
    end if

    term_year = 0
    
    ! loop on years and days
    year_loop: do iyear = 1, nyears
      day_loop: do iday = 1, days_per_year

        iday_all = (iyear - 1)*days_per_year + iday

        daily_net_c = 0.0_r8
        daily_gpp = 0.0_r8
        daily_rdark = 0.0_r8
        daily_livestem_mr = 0.0_r8
        daily_livecroot_mr = 0.0_r8
        daily_froot_mr = 0.0_r8
        daily_absorbed_par = 0.0_r8
        daily_incident_par = 0.0_r8
        daily_absorbed_par_indiv = 0.0_r8

        ! once-per-day setup: maintenance respiration factor, sapwood/fine-root N, and the per-layer
        ! nitrogen-scaling factor 
        call phys%DailySetup(cohort, pft, frac_store)

        ! today's storage-based maintenance-respiration factor
        ! recalculated for output
        ! TODO: can we avoid recalculating this somehow??
        call LowstorageMainRespReduction(frac_store, pft, maintresp_reduction_factor)

        do isubday = 1, n_substeps_per_day

          ! sub-daily environmental variables (temperature, VPD, light, etc.)
          hour_of_day = (real(isubday, r8) - 0.5_r8) * hrs_per_substep
          call env%SetHour(iday, hour_of_day)
          
          ! get par
          call env%CalculatePAR(par_beam, par_diff, par_toc, light_frac=light_frac_val)
          
          ! attenuate light through the canopy
          call light_env%AttenuateCanopy(par_beam, par_diff, env%coszen, &
            light_env%parsun_z, light_env%parsha_z, light_env%laisun_z, light_env%laisha_z)

          ! per-year light-profile snapshot: first day of each year, substep
          ! nearest solar noon
          if (iday == 1 .and. isubday == noon_substep) then
            call hist%RecordLightProfile(iyear, ilight, cohort, light_env)
          end if

          ! sub-daily carbon uptake 
          call phys%SubdailyStep(cohort, pft, env, light_env, lnc_top, step_size,  &
            n_photo_calls, n_bisection_calls, max_solve_iter, sum_solve_iter,      &
            gpp_tstep, rdark_tstep, nonleaf_mr_tstep)
          
          ! accumulate metrics for light interception efficiency
          if (light_env%treelai > nearzero) then
            daily_absorbed_par = daily_absorbed_par +                              &
              sum(light_env%parsun_z(:) + light_env%parsha_z(:)) / light_env%treelai * step_size
            daily_incident_par = daily_incident_par + par_toc * step_size
            daily_absorbed_par_indiv = daily_absorbed_par_indiv +                  &
              sum(light_env%parsun_z(:) + light_env%parsha_z(:)) * cohort%c_area / cohort%n * step_size
          end if

          ! instantaneous whole-plant light-response
          if (iday == 1 .and. isubday == noon_substep) then
            
            ! calculate gross_assim and total_resp across the PPFD values
            call LightResponseSweep(cohort, pft, light_env, phys, env,           &
              diagnostic_ppfd, light_frac_val, iday, hour_of_day,                &
              diagnostic_gross_assim, diagnostic_leaf_resp, diagnostic_total_resp)
              
            ! write to output
            call hist%RecordLightResponse(iyear, ilight, diagnostic_gross_assim, &
              diagnostic_leaf_resp, diagnostic_total_resp)

            ! leaf-level (LCPleaf) companion diagnostic
            call phys%LeafNetAssimSweep(pft, env, diagnostic_ppfd,               &
              leaf_lcp_layer, diagnostic_leaf_anet)
            call hist%RecordLeafNetAssim(iyear, ilight, diagnostic_leaf_anet)
          end if

          ! daily net carbon = GPP - leaf dark resp - non-leaf MR, integrated over
          ! the sub-daily steps [kgC/indiv/day]
          daily_net_c = daily_net_c + (gpp_tstep - rdark_tstep - nonleaf_mr_tstep) * step_size
          daily_gpp = daily_gpp + gpp_tstep * step_size
          daily_rdark = daily_rdark + rdark_tstep * step_size
          daily_livestem_mr = daily_livestem_mr + cohort%livestem_mr * step_size
          daily_livecroot_mr = daily_livecroot_mr + cohort%livecroot_mr * step_size
          daily_froot_mr = daily_froot_mr + cohort%froot_mr * step_size

        end do

        ! calculate light interception efficiency
        if (daily_incident_par > nearzero) then
          light_intercept_eff = daily_absorbed_par / (daily_incident_par * par_absorptance)
        else
          light_intercept_eff = 0.0_r8
        end if

        ! calculate daily means of environmental variables
        call env%UpdateDailyMeans()

        ! prescribed leaf-on/leaf-off phenology (cold-deciduous PFTs only)
        call DailyPhenology(cohort, pft, iday)

        ! the daily growth sequence and carbon starvation mortality
        call DailyGrowthAndMortality(cohort, daily_net_c, frac_store,            &
          npp_acc_to_prt, cmort, growth_resp, leaf_turnover, fnrt_turnover,      &
          sapw_turnover, struct_turnover)

        ! refresh crown area, treelai/treesai/nv, and the light environment's
        ! canopy structure to reflect today's growth
        call carea_allom(cohort%dbh, cohort%n, site_spread, pft,         &
          cohort%crowndamage, cohort%c_area)
        call UpdateCohortLAI(cohort, can_tlai, patch_area)
        call light_env%Refresh(cohort%treelai, cohort%treesai, cohort%height)

        ! capture today's daily time series row
        call hist%RecordDay(iday_all, ilight, cohort, daily_net_c, daily_gpp,      &
          daily_rdark, daily_livestem_mr, daily_livecroot_mr, daily_froot_mr,      &
          growth_resp, leaf_turnover, fnrt_turnover, sapw_turnover,                &
          struct_turnover, npp_acc_to_prt, frac_store, cmort,                      &
          light_intercept_eff, maintresp_reduction_factor, daily_incident_par,     &
          daily_absorbed_par_indiv, env)
          
        ! stop running this cohort if it meets criteria for being terminated
        term_reason = Terminate(cohort)
        if (term_reason /= term_none) then 
          term_year = iyear 
          exit year_loop
        end if 

      end do day_loop

      storage_c = cohort%prt%GetState(store_organ, carbon12_element)

    end do year_loop

    n_photo_calls_total = n_photo_calls_total + n_photo_calls
    n_bisection_calls_total = n_bisection_calls_total + n_bisection_calls
    max_solve_iter_total = max(max_solve_iter_total, max_solve_iter)

    ! this light level's whole-trajectory Ci-solver summary - mean iteration
    ! count and bisection-fallback count, over every LeafLayerPhotosynthesis
    ! call across all n_days_total days
    mean_solve_iter = real(sum_solve_iter, r8) / real(n_photo_calls, r8)
    call hist%RecordLightLevelSummary(ilight, mean_solve_iter, n_bisection_calls)
    
    call hist%RecordTermination(ilight, term_year, term_reason)

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
    ! Prescribed cold-deciduous (ihard_season_decid) leaf-on/leaf-off phenology

    ! ARGUMENTS:
    type(fates_cohort_type), intent(inout) :: cohort      ! cohort to update
    integer,                  intent(in)   :: pft         ! plant functional type index
    integer,                  intent(in)   :: day_of_year ! day of year [1-365]

    ! LOCALS:
    real(r8) :: elong_factor_today       ! today's leaf elongation factor [0 or 1]
    real(r8) :: leaf_c                   ! current leaf carbon [kgC]
    real(r8) :: fnrt_c                   ! current fineroot carbon [kgC]
    real(r8) :: sapw_c                   ! current sapwood carbon [kgC]
    real(r8) :: struct_c                 ! current structural carbon [kgC]
    real(r8) :: store_c                  ! current storage carbon [kgC]
    real(r8) :: leaf_deficit_c           ! leaf flush-time carbon deficit relative to target [kgC]
    real(r8) :: fnrt_deficit_c           ! fineroot flush-time carbon deficit relative to target [kgC]
    real(r8) :: sapw_deficit_c           ! sapwood flush-time carbon deficit relative to target [kgC]
    real(r8) :: struct_deficit_c         ! structural flush-time carbon deficit relative to target [kgC]
    real(r8) :: total_deficit_c          ! total flush-time carbon deficit relative to target [kgC]
    real(r8) :: eff_leaf_drop_fraction   ! leaf shed-time drop fractions [0-1]
    real(r8) :: eff_fnrt_drop_fraction   ! fineroot shed-time drop fractions [0-1]
    real(r8) :: eff_sapw_drop_fraction   ! sapwood shed-time drop fraction [0-1]
    real(r8) :: eff_struct_drop_fraction ! structural shed-time drop fraction [0-1]
    real(r8) :: target_leaf_c            ! target leaf carbon at today's elongation factor [kgC]
    real(r8) :: target_fnrt_c            ! target fineroot carbon at today's elongation factor [kgC]
    real(r8) :: target_sapw_c            ! target sapwood carbon at today's elongation factor [kgC]
    real(r8) :: target_agw_c             ! target aboveground carbon at today's elongation factor [kgC]
    real(r8) ::  target_bgw_c            ! target belowground carbon at today's elongation factor [kgC]
    real(r8) :: target_struct_c          ! target structural carbon at today's elongation factor [kgC]
    real(r8) :: sapw_area                ! unused diagnostic output of bsap_allom [m2]
    real(r8) :: store_c_transfer_frac    ! fraction of storage carbon transferred to flush tissues [-]
    logical  :: is_leaf_on_season        ! true on days within [leaf_on_doy, leaf_off_doy)
    logical  :: is_flushing_time         ! true only on the day leaves transition off -> on
    logical  :: is_shedding_time         ! true only on the day leaves transition on -> off
  
    ! CONSTANTS:
    real(r8), parameter :: carbon_store_buffer = 0.10_r8 ! matches phenology_leafonoff's identical local constant

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
    call bleaf(cohort%dbh, pft, cohort%crowndamage, cohort%canopy_trim,        &
      cohort%efleaf_coh, target_leaf_c)
    call bfineroot(cohort%dbh, pft, cohort%canopy_trim, cohort%l2fr,            &
      cohort%effnrt_coh, target_fnrt_c)
    call bsap_allom(cohort%dbh, pft, cohort%crowndamage, cohort%canopy_trim,    &
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

      leaf_c = cohort%prt%GetState(leaf_organ,   carbon12_element)
      fnrt_c = cohort%prt%GetState(fnrt_organ,   carbon12_element)
      sapw_c = cohort%prt%GetState(sapw_organ,   carbon12_element)
      struct_c = cohort%prt%GetState(struct_organ, carbon12_element)

      eff_leaf_drop_fraction = max(0.0_r8, min(1.0_r8, 1.0_r8 - target_leaf_c/max(leaf_c, nearzero)))
      eff_fnrt_drop_fraction = max(0.0_r8, min(1.0_r8, 1.0_r8 - target_fnrt_c/max(fnrt_c, nearzero)))
      eff_sapw_drop_fraction = max(0.0_r8, min(1.0_r8, 1.0_r8 - target_sapw_c/max(sapw_c, nearzero)))
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
    light_frac_val, day_of_year, hour_of_day, gross_assim, leaf_resp, total_resp, &
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
    ! advancing any state. This is a separate calculation from the daily/sub-daily loop
    !
    ! Illumination defaults to pure diffuse (all incident PPFD as diffuse, no
    ! direct beam) - matching Sterck et al. (2013)'s definition of light
    ! interception efficiency, the comparison target for this diagnostic.
    ! use_beam_illumination switches to the opposite extreme (all incident
    ! PPFD as direct beam at zenith) so the two can be compared. Sun-angle
    ! geometry is fixed at coszen=1 (sun directly overhead)

    ! ARGUMENTS:
    type(fates_cohort_type), intent(in)    :: cohort         ! cohort to diagnose (frozen - read-only)
    type(light_env_type),    intent(inout) :: light_env      ! this light level's light environment (mutated only transiently; restored before return)
    type(cohort_phys_type),  intent(in)    :: phys           ! today's already-computed leaf physiology (frozen - read-only)
    type(environment_type),  intent(in)    :: env            ! prescribed atmospheric/soil boundary conditions (frozen - read-only)
    real(r8),                intent(in)    :: ppfd_values(:) ! incident PPFD values to sweep, at the top of the crown [umol/m2/s]
    real(r8),                intent(in)    :: light_frac_val ! this light level's real incident light fraction [0-1] - used only to restore light_env afterward
    real(r8),                intent(in)    :: hour_of_day    ! real hour of day [0-24] - used only to restore light_env afterward
    integer,                 intent(in)    :: pft            ! plant functional type index
    integer,                 intent(in)    :: day_of_year    ! real day of year [1-365] - used only to restore light_env afterward
    real(r8),                intent(out)   :: gross_assim(:) ! whole-plant gross assimilation at each swept PPFD [kgC/indiv/s], size(ppfd_values)
    real(r8),                intent(out)   :: total_resp(:)  ! whole-plant total respiration at each swept PPFD [kgC/indiv/s], size(ppfd_values)
    real(r8),                intent(out)   :: leaf_resp(:)   ! leaf dark respiration at each swept PPFD [kgC/indiv/s], size(ppfd_values)
    logical,                 intent(in), optional :: use_beam_illumination ! if true, sweep with all incident PPFD as direct beam at zenith instead of the default pure diffuse (Sterck et al. 2013)

    ! LOCALS:
    real(r8), allocatable :: parsun_z(:) ! this sweep step's attenuated PAR profile (sun)
    real(r8), allocatable :: parsha_z(:) ! this sweep step's attenuated PAR profile (shade) 
    real(r8), allocatable :: laisun_z(:) ! this sweep step's lai profile (sun) 
    real(r8), allocatable :: laisha_z(:) ! this sweep step's lai profile (shade) 
    real(r8)              :: par_toc     ! this sweep step's total incident PAR at the top of the crown [W/m2]
    real(r8)              :: par_beam    ! this sweep step's direct-beam incident PAR at the top of the crown [W/m2]
    real(r8)              :: par_diff    ! this sweep step's diffuse incident PAR at the top of the crown [W/m2]
    integer               :: ippfd       ! PPFD-sweep looping index
    logical               :: use_beam    ! resolved value of use_beam_illumination (defaults to false: pure diffuse)

    ! CONSTANTS:
    real(r8), parameter :: sweep_coszen = 1.0_r8 ! fixed sun-angle geometry for the sweep (sun directly overhead)
  
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
        
      ! maintance respiration factor is set to 1.0 here; we want these outputs to be 
      ! properties of the plant's size and physiology rather than of how storage-depleted 
      ! it happened to be that day.
      call phys%GrossAssimAndResp(cohort, pft, env, parsun_z, parsha_z,              &
        laisun_z, laisha_z, 1.0_r8, step_size, gross_assim(ippfd), leaf_resp(ippfd), &
        total_resp(ippfd))
    end do

    deallocate(parsun_z, parsha_z, laisun_z, laisha_z)

    ! restore light_env's real per-substep state (see header comment)
    call env%CalculatePAR(par_beam, par_diff, par_toc, light_frac=light_frac_val)
    call light_env%AttenuateCanopy(par_beam, par_diff, env%coszen, &
      light_env%parsun_z, light_env%parsha_z, light_env%laisun_z, light_env%laisha_z)

  end subroutine LightResponseSweep

  ! ==========================================================================

  subroutine DailyGrowthAndMortality(cohort, daily_net_c, frac_store,          &
    npp_acc_to_prt, cmort, growth_resp, leaf_turnover, fnrt_turnover,          &
    sapw_turnover, struct_turnover)
    !
    ! DESCRIPTION:
    ! The daily growth sequence followed by carbon starvation mortality

    ! ARGUMENTS:
    type(fates_cohort_type), intent(inout) :: cohort          ! cohort to grow and apply mortality to
    real(r8),                 intent(in)   :: daily_net_c     ! today's net carbon (GPP - leaf Rd - nonleaf MR) [kgC/indiv/day]
    real(r8),                 intent(in)   :: frac_store      ! storage carbon / target leaf carbon, computed pre-growth at the top of today's loop [-]
    real(r8),                 intent(out)  :: npp_acc_to_prt  ! carbon balance handed to PARTEH (net of growth respiration) [kgC/indiv/day]
    real(r8),                 intent(out)  :: cmort           ! carbon starvation mortality rate [indiv/year]
    real(r8),                 intent(out)  :: growth_resp     ! today's growth respiration [kgC/indiv/day]
    real(r8),                 intent(out)  :: leaf_turnover   ! today's leaf turnover loss, from PRTMaintTurnover [kgC/indiv/day]
    real(r8),                 intent(out)  :: fnrt_turnover   ! today's fine root turnover loss, from PRTMaintTurnover [kgC/indiv/day]
    real(r8),                 intent(out)  :: sapw_turnover   ! today's sapwood turnover loss, from PRTMaintTurnover [kgC/indiv/day]
    real(r8),                 intent(out)  :: struct_turnover ! today's structural turnover loss, from PRTMaintTurnover [kgC/indiv/day]

    ! LOCALS:
    real(r8) :: resp_g_acc      ! growth respiration: a tax on the positive part of daily_net_c [kgC/indiv/day]
    real(r8) :: delta_dbh       ! unused output of EvaluateAndCorrectDBH
    real(r8) :: delta_height    ! unused output of EvaluateAndCorrectDBH
    real(r8) :: leaf_c_before   ! leaf carbon immediately before PRTMaintTurnover [kgC/indiv]
    real(r8) :: fnrt_c_before   ! fineroot carbon immediately before PRTMaintTurnover [kgC/indiv]
    real(r8) :: sapw_c_before   ! sapwood carbon immediately before PRTMaintTurnover [kgC/indiv]
    real(r8) :: struct_c_before ! structural carbon immediately before PRTMaintTurnover [kgC/indiv]

    ! ---------------------------------------------------------------------
    ! NSC update / PRT allocation - once per day, after the sub-daily loop
    ! ---------------------------------------------------------------------

    ! growth respiration: a tax on the positive part of today's net carbon
    resp_g_acc = prt_params%grperc(pft) * max(0.0_r8, daily_net_c)
    growth_resp = resp_g_acc

    ! carbon_balance boundary condition
    cohort%npp_acc = daily_net_c - resp_g_acc
    npp_acc_to_prt = cohort%npp_acc

    ! maintenance turnover moves senesced tissue mass to storage/litter
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
    ! carbon starvation mortality. Uses frac_store as computed at the top of today's loop 
    ! (pre-growth), matching production's Mortality_Derivative, which runs before
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
    ! to n 
    cohort%n = max(0.0_r8, cohort%n * (1.0_r8 - cmort/real(days_per_year, r8)))

  end subroutine DailyGrowthAndMortality
  
  ! ==========================================================================

  integer function Terminate(cohort) result(reason)
    !
    ! DESCRIPTION:
    ! 
    ! Check whether this cohort should be terminated for some reason

    ! ARGUMENTS:
    type(fates_cohort_type), intent(in) :: cohort ! cohort

    ! LOCALS:
    real(r8) :: leaf_c   ! leaf carbon [kgC/indiv]
    real(r8) :: store_c  ! storage carbon [kgC/indiv]
    real(r8) :: sapw_c   ! sapwood carbon [kgC/indiv]
    real(r8) :: fnrt_c   ! fine root carbon [kgC/indiv]
    real(r8) :: struct_c ! structural carbon [kgC/indiv]

    ! CONSTANTS:
    real(r8), parameter :: pool_min = 1.0e-10_r8  ! depletion threshold [kgC/indiv]

    reason = term_none

    leaf_c = cohort%prt%GetState(leaf_organ, carbon12_element)
    store_c = cohort%prt%GetState(store_organ, carbon12_element)
    sapw_c = cohort%prt%GetState(sapw_organ, carbon12_element)
    fnrt_c = cohort%prt%GetState(fnrt_organ, carbon12_element)
    struct_c = cohort%prt%GetState(struct_organ, carbon12_element)

    ! ordered so the earliest-firing criterion is reported: in a carbon-starved
    ! cohort storage crosses zero years before the live pools run out
    if (store_c < pool_min) then
      reason = term_storage
    else if ((sapw_c + leaf_c + fnrt_c) < pool_min) then
      reason = term_live_pools
    else if ((struct_c + sapw_c + leaf_c + fnrt_c + store_c) < 0.0_r8) then
      reason = term_neg_biomass
    else if (cohort%n < min_n_safemath) then
      reason = term_num_dens
    end if

  end function Terminate

end program FatesSingleCohort
