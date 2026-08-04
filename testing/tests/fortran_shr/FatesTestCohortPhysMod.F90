module FatesTestCohortPhysMod
  !
  ! DESCRIPTION:
  ! Per-leaf-layer photosynthetic capacity/dark-respiration working arrays for a
  ! single cohort, plus the two pieces of leaf physics a standalone, patch-less/
  ! site-less cohort test driver needs each simulated day: a once-per-day setup
  ! (DailySetup) that caches the per-layer nitrogen-scaling factor (nscaler_z -
  ! depends on cumulative LAI above each layer, which only changes daily), and a
  ! per-substep carbon uptake step (SubdailyStep) that both refreshes the
  ! temperature-dependent capacity/dark-respiration rates from nscaler_z and
  ! this substep's env%tempk (matching production's per-timestep
  ! LeafLayerBiophysicalRates pattern, needed now that env%tempk has a real
  ! diurnal cycle - see FatesTestEnvironmentMod) and consumes them for
  ! photosynthesis.
  !
  ! Leaf photosynthesis (SubdailyStep) is wired in below: one LeafLayerPhotosynthesis
  ! call per leaf layer per sunlit/shaded category, then scaled up to a
  ! per-individual carbon flux (cohort%gpp_tstep) mirroring
  ! FatesPlantRespPhotosynthMod's ScaleLeafLayerFluxToCohort (which is private to
  ! that module, so this reimplements its scaling formula rather than reusing it).
  ! Leaf dark respiration (cohort%rdark) falls out of this as an unavoidable side
  ! effect, since it is a required input to LeafLayerPhotosynthesis itself. Non-leaf
  ! (stem/root) maintenance respiration is wired in via NonleafMaintenanceRespiration.
  !

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesConstantsMod,           only : nearzero
  use FatesConstantsMod,           only : umolC_to_kgC
  use FatesConstantsMod,           only : wm2_to_umolm2s
  use FatesCohortMod,              only : fates_cohort_type
  use FatesAllometryMod,           only : VegAreaLayer
  use FatesAllometryMod,           only : bleaf
  use FatesAllometryMod,           only : storage_fraction_of_target
  use PRTParametersMod,            only : prt_params
  use PRTGenericMod,                only : store_organ, sapw_organ, fnrt_organ, carbon12_element
  use FatesTestEnvironmentMod,     only : environment_type
  use FatesTestLightEnvMod,        only : light_env_type
  use FatesPlantRespPhotosynthMod, only : NonleafMaintenanceRespiration
  use LeafBiophysicsMod,           only : LeafLayerBiophysicalRates
  use LeafBiophysicsMod,           only : LeafLayerMaintenanceRespiration_Ryan_1991
  use LeafBiophysicsMod,           only : LeafLayerPhotosynthesis
  use LeafBiophysicsMod,           only : DecayCoeffVcmax
  use LeafBiophysicsMod,           only : LowstorageMainRespReduction
  use LeafBiophysicsMod,           only : GetCanopyGasParameters

  implicit none
  private

  real(r8), parameter :: ci_tol = 0.5_r8               ! convergence tolerance for intracellular CO2 [Pa], matches FatesPlantRespPhotosynthMod
  real(r8), parameter :: min_la_to_solve = 1.0e-10_r8  ! minimum leaf area to bother solving photosynthesis for [m2], matches ConvertPar in FatesPlantRespPhotosynthMod
  integer,  parameter :: newton_max_iters = 10          ! matches LeafLayerPhotosynthesis's private max_iters - solve_iter beyond this means the call fell back to CiBisection

  type, public :: cohort_phys_type

     private

     real(r8), allocatable :: vcmax_z(:), jmax_z(:), kp_z(:) ! per-leaf-layer photosynthetic capacity
     real(r8), allocatable :: gs0_z(:), gs1_z(:), gs2_z(:)   ! per-leaf-layer stomatal conductance slope/intercept
     real(r8), allocatable :: lmr_z(:)                        ! per-leaf-layer leaf dark respiration [umolCO2/m2 leaf/s]
     real(r8), allocatable :: nscaler_z(:)                    ! per-leaf-layer nitrogen-scaling factor [0-1], cached daily by DailySetup, consumed every substep by SubdailyStep

     real(r8) :: live_stem_n, live_croot_n   ! aboveground/belowground sapwood N [kgN/plant] - feed NonleafMaintenanceRespiration
     real(r8) :: fnrt_n                      ! fine root N [kgN/plant] - feeds NonleafMaintenanceRespiration
     real(r8) :: maintresp_reduction_factor  ! storage-based throttle on maintenance respiration [0-1]

   contains

     procedure, public :: DailySetup
     procedure, public :: SubdailyStep
     procedure, public :: GrossAssimAndResp

  end type cohort_phys_type

contains

  ! ==========================================================================

  subroutine EnsureAllocated(this, nv)
    !
    ! DESCRIPTION:
    ! (Re)allocate the per-leaf-layer working arrays if the cohort's occupied
    ! leaf-layer count has changed since the last call.

    ! ARGUMENTS:
    class(cohort_phys_type), intent(inout) :: this ! cohort physiology object
    integer,                  intent(in)    :: nv   ! cohort's current number of occupied leaf layers

    if (allocated(this%vcmax_z)) then
      if (size(this%vcmax_z) /= nv) then
        deallocate(this%vcmax_z, this%jmax_z, this%kp_z)
        deallocate(this%gs0_z, this%gs1_z, this%gs2_z, this%lmr_z)
        deallocate(this%nscaler_z)
        allocate(this%vcmax_z(nv), this%jmax_z(nv), this%kp_z(nv))
        allocate(this%gs0_z(nv), this%gs1_z(nv), this%gs2_z(nv))
        allocate(this%lmr_z(nv))
        allocate(this%nscaler_z(nv))
      end if
    else
      allocate(this%vcmax_z(nv), this%jmax_z(nv), this%kp_z(nv))
      allocate(this%gs0_z(nv), this%gs1_z(nv), this%gs2_z(nv))
      allocate(this%lmr_z(nv))
      allocate(this%nscaler_z(nv))
    end if

  end subroutine EnsureAllocated

  ! ==========================================================================

  subroutine DailySetup(this, cohort, pft, frac_store)
    !
    ! DESCRIPTION:
    ! Once-per-day setup: the storage-based maintenance-respiration throttle,
    ! sapwood/fine-root N (feed NonleafMaintenanceRespiration in SubdailyStep), and
    ! the per-leaf-layer nitrogen-scaling factor (nscaler_z - depends on
    ! cumulative LAI above each layer, which only changes daily). The
    ! temperature-dependent capacity/dark-respiration rates that consume
    ! nscaler_z are computed per-substep in SubdailyStep instead, from this
    ! day's nscaler_z and that substep's env%tempk.

    ! ARGUMENTS:
    class(cohort_phys_type), intent(inout) :: this       ! cohort physiology object
    type(fates_cohort_type), intent(in)    :: cohort      ! cohort to set up for today
    integer,                  intent(in)    :: pft         ! plant functional type index
    real(r8),                 intent(out)   :: frac_store  ! ratio of storage carbon to target_leaf_c [-]

    ! LOCALS:
    real(r8) :: target_leaf_c  ! reference leaf biomass when fully flushed [kgC] (see EDMortalityFunctionsMod.F90's comment on this same calculation - the target intentionally ignores the plant's current elongation state)
    real(r8) :: sapw_c, fnrt_c ! sapwood/fine root carbon [kgC/plant]
    integer  :: iv             ! leaf-layer looping index
    real(r8) :: vai_top, vai_bot       ! vegetation area index bounds of the current leaf layer
    real(r8) :: elai_layer, esai_layer ! exposed leaf/stem area index of the current leaf layer
    real(r8) :: cumulative_lai         ! LAI above the middle of the current leaf layer
    real(r8) :: lai_above              ! running LAI above the current leaf layer
    real(r8) :: kn                     ! nitrogen vertical-scaling decay coefficient

    ! storage-based throttle on maintenance respiration. The elongation factor
    ! is intentionally 1.0 here (fully flushed), not the cohort's current
    ! elongf_leaf - matches EDMortalityFunctionsMod.F90's identical calculation,
    ! whose comment explains the reference target must stay at full-flush so the
    ! ratio remains meaningful through leaf-off periods for deciduous PFTs.
    ! cohort%dbh/treelai/nv/c_area are already current here - either from
    ! CohortFactory (day 1) or from the refresh at the end of the previous day
    call bleaf(cohort%dbh, pft, cohort%crowndamage, cohort%canopy_trim, 1.0_r8, &
      target_leaf_c)
    call storage_fraction_of_target(target_leaf_c,                             &
      cohort%prt%GetState(store_organ, carbon12_element), frac_store)
    call LowstorageMainRespReduction(frac_store, pft, this%maintresp_reduction_factor)

    ! aboveground/belowground sapwood N and fine root N - feed
    ! NonleafMaintenanceRespiration in SubdailyStep. This driver never damages the
    ! crown (crowndamage=1, undamaged), so the aboveground/belowground sapwood
    ! split is just allom_agb_frac with no damage correction, unlike the full
    ! crown-damage-aware split in FatesPlantRespPhotosynthMod.F90
    sapw_c = cohort%prt%GetState(sapw_organ, carbon12_element)
    fnrt_c = cohort%prt%GetState(fnrt_organ, carbon12_element)
    this%live_stem_n  = sapw_c * prt_params%allom_agb_frac(pft) *              &
      prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(sapw_organ))
    this%live_croot_n = sapw_c * (1.0_r8 - prt_params%allom_agb_frac(pft)) *   &
      prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(sapw_organ))
    this%fnrt_n = fnrt_c * prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(fnrt_organ))

    call EnsureAllocated(this, cohort%nv)

    lai_above = 0.0_r8
    do iv = 1, cohort%nv
      call VegAreaLayer(cohort%treelai, cohort%treesai, cohort%height, iv,     &
        cohort%nv, pft, 0.0_r8, vai_top, vai_bot, elai_layer, esai_layer)
      cumulative_lai = lai_above + 0.5_r8*elai_layer
      lai_above = lai_above + elai_layer

      kn = DecayCoeffVcmax(cohort%vcmax25top,                                  &
        prt_params%leafn_vert_scaler_coeff1(pft),                             &
        prt_params%leafn_vert_scaler_coeff2(pft))
      this%nscaler_z(iv) = exp(-kn*cumulative_lai)
    end do

  end subroutine DailySetup

  ! ==========================================================================

  subroutine SubdailyStep(this, cohort, pft, env, light_env, lnc_top, step_size, &
    n_photo_calls, n_bisection_calls, max_solve_iter, sum_solve_iter,        &
    gpp_tstep, rdark_tstep, nonleaf_mr_tstep)
    !
    ! DESCRIPTION:
    ! A single sub-daily carbon uptake step: refresh the temperature-dependent
    ! leaf photosynthetic capacity/dark-respiration rates from today's cached
    ! nscaler_z (DailySetup) and this substep's env%tempk (matching production's
    ! per-timestep LeafLayerBiophysicalRates pattern), then prescribed absorbed
    ! PAR (light_env) -> leaf photosynthesis -> a per-individual GPP and leaf
    ! dark respiration, plus non-leaf (stem/root) maintenance respiration for
    ! the same substep.

    ! ARGUMENTS:
    class(cohort_phys_type), intent(inout) :: this        ! cohort physiology object
    type(fates_cohort_type), intent(inout) :: cohort       ! cohort to step
    integer,                  intent(in)    :: pft          ! plant functional type index
    type(environment_type),  intent(in)    :: env          ! prescribed atmospheric/soil boundary conditions
    type(light_env_type),    intent(in)    :: light_env    ! this substep's absorbed-PAR/LAI profile
    real(r8),                 intent(in)    :: lnc_top      ! leaf N content at the canopy top [gN/m2 leaf]
    real(r8),                 intent(in)    :: step_size    ! model time step [s]
    integer,                  intent(inout) :: n_photo_calls     ! running count of LeafLayerPhotosynthesis calls
    integer,                  intent(inout) :: n_bisection_calls ! running count of calls that fell back to CiBisection
    integer,                  intent(inout) :: max_solve_iter    ! running max Ci-solver iteration count
    integer,                  intent(inout) :: sum_solve_iter    ! running sum of Ci-solver iteration counts, for a mean once divided by n_photo_calls
    real(r8),                 intent(out)   :: gpp_tstep        ! this substep's GPP [kgC/indiv/s]
    real(r8),                 intent(out)   :: rdark_tstep      ! this substep's leaf dark respiration [kgC/indiv/s]
    real(r8),                 intent(out)   :: nonleaf_mr_tstep ! this substep's non-leaf (stem/root) maintenance respiration [kgC/indiv/s]

    ! LOCALS:
    real(r8) :: mm_kco2      ! Michaelis-Menten constant for CO2 at this substep's env%tempk [Pa]
    real(r8) :: mm_ko2       ! Michaelis-Menten constant for O2 at this substep's env%tempk [Pa]
    real(r8) :: co2_cpoint   ! CO2 compensation point at this substep's env%tempk [Pa]
    integer  :: iv                         ! leaf-layer looping index
    real(r8) :: par_abs                    ! absorbed PAR per unit leaf area [umol photons/m2 leaf/s]
    real(r8) :: agross_sun, agross_sha     ! gross photosynthesis, sunlit/shaded [umolC/m2 leaf/s]
    real(r8) :: anet_sun, anet_sha         ! net photosynthesis, sunlit/shaded (unused diagnostics)
    real(r8) :: gs_sun, gs_sha             ! stomatal conductance, sunlit/shaded (unused diagnostics)
    real(r8) :: ci_sun, ci_sha             ! intracellular CO2, sunlit/shaded (unused diagnostics)
    real(r8) :: c13disc_sun, c13disc_sha   ! carbon-13 discrimination, sunlit/shaded (unused diagnostics)
    integer  :: solve_iter                 ! Ci-solver iteration count - >newton_max_iters means it fell back to CiBisection
    real(r8) :: fsun                       ! sunlit fraction of leaf area in the current layer
    real(r8) :: psn_layer                  ! area-weighted mean gross photosynthesis for the layer [umolC/m2 leaf/s]
    real(r8) :: cohort_layer_eleaf_area    ! cohort's effective leaf area in the current layer [m2]
    real(r8) :: gpp_sum, rdark_sum         ! running per-substep GPP/dark-respiration accumulators [umolC/s]

    ! CO2/O2 kinetics are temperature-dependent (Arrhenius-type - see
    ! GetCanopyGasParameters), so with env%tempk now varying diurnally/annually
    ! (FatesTestEnvironmentMod) these must be re-evaluated every substep rather
    ! than once at a fixed reference temperature
    call GetCanopyGasParameters(env%can_press, env%can_o2_ppress, env%tempk,    &
      mm_kco2, mm_ko2, co2_cpoint)

    ! prescribed absorbed PAR -> leaf photosynthesis: one LeafLayerPhotosynthesis
    ! call per leaf layer per sunlit/shaded category (converting parsun_z/parsha_z
    ! from W/m2 ground to umol photons/m2 leaf/s the same way
    ! FatesPlantRespPhotosynthMod's ConvertPar does), then combined into an
    ! area-weighted per-layer rate and scaled up to a per-individual carbon flux -
    ! reimplements ScaleLeafLayerFluxToCohort's formula, which is private to
    ! FatesPlantRespPhotosynthMod and so not directly reusable here
    gpp_sum = 0.0_r8
    rdark_sum = 0.0_r8
    do iv = 1, cohort%nv

      ! refresh this layer's temperature-dependent capacity/dark-respiration
      ! rates from today's cached nscaler_z (DailySetup) and this substep's
      ! env%tempk - see this subroutine's header comment
      call LeafLayerBiophysicalRates(pft, cohort%vcmax25top, cohort%jmax25top, &
        cohort%kp25top, this%nscaler_z(iv), env%tempk, env%dayl_factor,       &
        env%t_growth, env%t_home, env%btran, this%vcmax_z(iv),               &
        this%jmax_z(iv), this%kp_z(iv), this%gs0_z(iv), this%gs1_z(iv),      &
        this%gs2_z(iv))

      call LeafLayerMaintenanceRespiration_Ryan_1991(lnc_top,                  &
        this%nscaler_z(iv), pft, env%tempk, this%lmr_z(iv))

      if (light_env%laisun_z(iv) > min_la_to_solve .and.                       &
        light_env%parsun_z(iv) > nearzero) then
        par_abs = light_env%parsun_z(iv)/light_env%laisun_z(iv)*wm2_to_umolm2s
      else
        par_abs = 0.0_r8
      end if
      call LeafLayerPhotosynthesis(par_abs, pft, this%vcmax_z(iv),             &
        this%jmax_z(iv), this%kp_z(iv), this%gs0_z(iv), this%gs1_z(iv),        &
        this%gs2_z(iv), env%tempk, env%can_press, env%can_co2_ppress,          &
        env%can_o2_ppress, env%veg_esat, env%gb, env%can_vpress, mm_kco2,      &
        mm_ko2, co2_cpoint, this%lmr_z(iv), ci_tol, agross_sun, gs_sun,        &
        anet_sun, c13disc_sun, ci_sun, solve_iter)
      n_photo_calls = n_photo_calls + 1
      if (solve_iter > newton_max_iters) n_bisection_calls = n_bisection_calls + 1
      max_solve_iter = max(max_solve_iter, solve_iter)
      sum_solve_iter = sum_solve_iter + solve_iter

      if (light_env%laisha_z(iv) > min_la_to_solve .and.                       &
        light_env%parsha_z(iv) > nearzero) then
        par_abs = light_env%parsha_z(iv)/light_env%laisha_z(iv)*wm2_to_umolm2s
      else
        par_abs = 0.0_r8
      end if
      call LeafLayerPhotosynthesis(par_abs, pft, this%vcmax_z(iv),             &
        this%jmax_z(iv), this%kp_z(iv), this%gs0_z(iv), this%gs1_z(iv),        &
        this%gs2_z(iv), env%tempk, env%can_press, env%can_co2_ppress,          &
        env%can_o2_ppress, env%veg_esat, env%gb, env%can_vpress, mm_kco2,      &
        mm_ko2, co2_cpoint, this%lmr_z(iv), ci_tol, agross_sha, gs_sha,        &
        anet_sha, c13disc_sha, ci_sha, solve_iter)
      n_photo_calls = n_photo_calls + 1
      if (solve_iter > newton_max_iters) n_bisection_calls = n_bisection_calls + 1
      max_solve_iter = max(max_solve_iter, solve_iter)
      sum_solve_iter = sum_solve_iter + solve_iter

      ! area-weighted mean gross photosynthesis for the layer, sunlit + shaded
      if (light_env%laisun_z(iv) + light_env%laisha_z(iv) > nearzero) then
        fsun = light_env%laisun_z(iv) /                                        &
          (light_env%laisun_z(iv) + light_env%laisha_z(iv))
      else
        fsun = 0.0_r8
      end if
      psn_layer = fsun*agross_sun + (1.0_r8 - fsun)*agross_sha

      ! [umolC/m2 leaf/s] * [m2 leaf] -> [umolC/s], leaf area per this cohort
      ! in this layer = layer's exposed LAI * the cohort's own crown area
      cohort_layer_eleaf_area = (light_env%laisun_z(iv) +                      &
        light_env%laisha_z(iv)) * cohort%c_area
      gpp_sum   = gpp_sum   + psn_layer * cohort_layer_eleaf_area
      rdark_sum = rdark_sum + this%lmr_z(iv) * cohort_layer_eleaf_area

    end do

    ! [umolC/s] -> [kgC/indiv/s]
    gpp_tstep = gpp_sum * umolC_to_kgC / cohort%n
    rdark_tstep = rdark_sum * umolC_to_kgC * this%maintresp_reduction_factor / cohort%n
    cohort%gpp_tstep = gpp_tstep
    cohort%rdark = rdark_tstep

    ! non-leaf (stem/root) maintenance respiration - Q10-scaled by env%tempk
    ! (stem) and env%t_soil (roots), throttled by the same
    ! maintresp_reduction_factor applied to rdark above, matching how
    ! FatesPlantRespPhotosynthMod combines the two (resp_m_tstep = livestem_mr +
    ! livecroot_mr + froot_mr + rdark)
    call NonleafMaintenanceRespiration(pft, env%tempk, env%nlevsoil,           &
      [env%t_soil], [env%rootfr_ft], this%live_stem_n, this%live_croot_n,      &
      this%fnrt_n, this%maintresp_reduction_factor, step_size,                 &
      cohort%livestem_mr, cohort%livecroot_mr, cohort%froot_mr,                &
      cohort%sym_nfix_tstep)
    nonleaf_mr_tstep = cohort%livestem_mr + cohort%livecroot_mr + cohort%froot_mr

  end subroutine SubdailyStep

  ! ==========================================================================

  subroutine GrossAssimAndResp(this, cohort, pft, env, parsun_z, parsha_z,       &
    laisun_z, laisha_z, step_size, gross_assim, total_resp)
    !
    ! DESCRIPTION:
    ! Whole-plant gross assimilation and total respiration (leaf dark + non-leaf
    ! maintenance) at an arbitrary already-attenuated PAR profile
    ! (parsun_z/parsha_z/laisun_z/laisha_z - typically from
    ! FatesTestLightEnvMod's AttenuateCanopy at a diagnostic incident PPFD, not
    ! this substep's real light_env%parsun_z etc.), using this cohort's current
    ! ("frozen") per-layer capacity (vcmax_z/jmax_z/kp_z/gs0_z/gs1_z/gs2_z/lmr_z)
    ! exactly as SubdailyStep left them - NOT recomputed here, since a diagnostic
    ! sweep across incident PPFD is not supposed to advance any state. Unlike
    ! SubdailyStep, nothing here is written to cohort (gpp_tstep/rdark/
    ! livestem_mr/livecroot_mr/froot_mr/sym_nfix_tstep all stay untouched) or to
    ! this (vcmax_z etc.), so there is nothing to restore afterward - this is a
    ! pure read of already-existing state, safe to call repeatedly without
    ! perturbing the driver's actual daily/sub-daily trajectory.

    ! ARGUMENTS:
    class(cohort_phys_type), intent(in) :: this          ! cohort physiology object (read-only: today's already-computed capacity)
    type(fates_cohort_type), intent(in) :: cohort         ! cohort to diagnose (read-only: c_area, n, nv)
    integer,                  intent(in) :: pft            ! plant functional type index
    type(environment_type),  intent(in) :: env            ! prescribed atmospheric/soil boundary conditions
    real(r8),                 intent(in) :: parsun_z(:)    ! absorbed PAR, sunlit, per leaf layer [W/m2 ground] - from AttenuateCanopy, not light_env%parsun_z
    real(r8),                 intent(in) :: parsha_z(:)    ! absorbed PAR, shaded, per leaf layer [W/m2 ground] - from AttenuateCanopy, not light_env%parsha_z
    real(r8),                 intent(in) :: laisun_z(:)    ! sunlit LAI per leaf layer [m2/m2] - from AttenuateCanopy, not light_env%laisun_z
    real(r8),                 intent(in) :: laisha_z(:)    ! shaded LAI per leaf layer [m2/m2] - from AttenuateCanopy, not light_env%laisha_z
    real(r8),                 intent(in) :: step_size      ! model time step [s] - passed through to NonleafMaintenanceRespiration; does not affect the carbon-flux rates read here (only its own sym_nfix_tstep output, which this diagnostic discards)
    real(r8),                 intent(out) :: gross_assim   ! whole-plant gross assimilation at this PAR profile [kgC/indiv/s]
    real(r8),                 intent(out) :: total_resp     ! whole-plant total respiration (leaf dark + non-leaf maintenance) at this PAR profile [kgC/indiv/s]

    ! LOCALS:
    real(r8) :: mm_kco2, mm_ko2, co2_cpoint ! Michaelis-Menten constants for CO2/O2, CO2 compensation point at env%tempk [Pa]
    integer  :: iv                         ! leaf-layer looping index
    real(r8) :: par_abs                    ! absorbed PAR per unit leaf area [umol photons/m2 leaf/s]
    real(r8) :: agross_sun, agross_sha     ! gross photosynthesis, sunlit/shaded [umolC/m2 leaf/s]
    real(r8) :: anet_sun, anet_sha         ! net photosynthesis, sunlit/shaded (unused diagnostics)
    real(r8) :: gs_sun, gs_sha             ! stomatal conductance, sunlit/shaded (unused diagnostics)
    real(r8) :: ci_sun, ci_sha             ! intracellular CO2, sunlit/shaded (unused diagnostics)
    real(r8) :: c13disc_sun, c13disc_sha   ! carbon-13 discrimination, sunlit/shaded (unused diagnostics)
    integer  :: solve_iter                 ! Ci-solver iteration count (unused diagnostic here - see SubdailyStep for the tracked version)
    real(r8) :: fsun                       ! sunlit fraction of leaf area in the current layer
    real(r8) :: psn_layer                  ! area-weighted mean gross photosynthesis for the layer [umolC/m2 leaf/s]
    real(r8) :: cohort_layer_eleaf_area    ! cohort's effective leaf area in the current layer [m2]
    real(r8) :: gross_assim_sum, leaf_resp_sum ! running whole-plant accumulators [umolC/s]
    real(r8) :: livestem_mr, livecroot_mr, froot_mr, sym_nfix_tstep ! NonleafMaintenanceRespiration outputs (local - not cohort%livestem_mr etc.)

    call GetCanopyGasParameters(env%can_press, env%can_o2_ppress, env%tempk,    &
      mm_kco2, mm_ko2, co2_cpoint)

    gross_assim_sum = 0.0_r8
    leaf_resp_sum = 0.0_r8
    do iv = 1, cohort%nv

      if (laisun_z(iv) > min_la_to_solve .and. parsun_z(iv) > nearzero) then
        par_abs = parsun_z(iv)/laisun_z(iv)*wm2_to_umolm2s
      else
        par_abs = 0.0_r8
      end if
      call LeafLayerPhotosynthesis(par_abs, pft, this%vcmax_z(iv),              &
        this%jmax_z(iv), this%kp_z(iv), this%gs0_z(iv), this%gs1_z(iv),        &
        this%gs2_z(iv), env%tempk, env%can_press, env%can_co2_ppress,          &
        env%can_o2_ppress, env%veg_esat, env%gb, env%can_vpress, mm_kco2,      &
        mm_ko2, co2_cpoint, this%lmr_z(iv), ci_tol, agross_sun, gs_sun,        &
        anet_sun, c13disc_sun, ci_sun, solve_iter)

      if (laisha_z(iv) > min_la_to_solve .and. parsha_z(iv) > nearzero) then
        par_abs = parsha_z(iv)/laisha_z(iv)*wm2_to_umolm2s
      else
        par_abs = 0.0_r8
      end if
      call LeafLayerPhotosynthesis(par_abs, pft, this%vcmax_z(iv),              &
        this%jmax_z(iv), this%kp_z(iv), this%gs0_z(iv), this%gs1_z(iv),        &
        this%gs2_z(iv), env%tempk, env%can_press, env%can_co2_ppress,          &
        env%can_o2_ppress, env%veg_esat, env%gb, env%can_vpress, mm_kco2,      &
        mm_ko2, co2_cpoint, this%lmr_z(iv), ci_tol, agross_sha, gs_sha,        &
        anet_sha, c13disc_sha, ci_sha, solve_iter)

      if (laisun_z(iv) + laisha_z(iv) > nearzero) then
        fsun = laisun_z(iv) / (laisun_z(iv) + laisha_z(iv))
      else
        fsun = 0.0_r8
      end if
      psn_layer = fsun*agross_sun + (1.0_r8 - fsun)*agross_sha

      cohort_layer_eleaf_area = (laisun_z(iv) + laisha_z(iv)) * cohort%c_area
      gross_assim_sum = gross_assim_sum + psn_layer * cohort_layer_eleaf_area
      leaf_resp_sum   = leaf_resp_sum   + this%lmr_z(iv) * cohort_layer_eleaf_area

    end do

    ! [umolC/s] -> [kgC/indiv/s]
    gross_assim = gross_assim_sum * umolC_to_kgC / cohort%n

    call NonleafMaintenanceRespiration(pft, env%tempk, env%nlevsoil,           &
      [env%t_soil], [env%rootfr_ft], this%live_stem_n, this%live_croot_n,      &
      this%fnrt_n, this%maintresp_reduction_factor, step_size,                 &
      livestem_mr, livecroot_mr, froot_mr, sym_nfix_tstep)

    total_resp = leaf_resp_sum * umolC_to_kgC * this%maintresp_reduction_factor / cohort%n + &
      livestem_mr + livecroot_mr + froot_mr

  end subroutine GrossAssimAndResp

end module FatesTestCohortPhysMod
