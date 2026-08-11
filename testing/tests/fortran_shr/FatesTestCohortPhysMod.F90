module FatesTestCohortPhysMod
  !
  ! DESCRIPTION:
  ! Per-leaf-layer photosynthetic capacity/dark-respiration working arrays for a
  ! single cohort, plus the two pieces of leaf physics a standalone, patch-less/
  ! site-less cohort test driver needs each simulated day: a once-per-day setup
  ! (DailySetup) that caches the per-layer nitrogen-scaling factor (nscaler_z -
  ! depends on cumulative LAI above each layer, which only changes daily, and
  ! is computed by FatesTestLeafPhotoMod's shared LeafLayerNitrogenScaling so a
  ! prescribed-canopy driver with no cohort can derive it identically), and a
  ! per-substep carbon uptake step (SubdailyStep) that both refreshes the
  ! temperature-dependent capacity/dark-respiration rates from nscaler_z and
  ! this substep's env%tempk (matching production's per-timestep
  ! LeafLayerBiophysicalRates pattern, needed now that env%tempk has a real
  ! diurnal cycle - see FatesTestEnvironmentMod) and consumes them for
  ! photosynthesis.
  !
  ! Leaf photosynthesis (SubdailyStep) is wired in below: one
  ! FatesTestLeafPhotoMod::LeafLayerCapacity + LeafLayerSunShade pair per leaf
  ! layer (the latter evaluating both the sunlit and shaded categories at that
  ! capacity and area-weighting them), then scaled up to a per-individual carbon flux
  ! (cohort%gpp_tstep) mirroring FatesPlantRespPhotosynthMod's
  ! ScaleLeafLayerFluxToCohort (which is private to that module, so this
  ! reimplements its scaling formula rather than reusing it). Leaf dark
  ! respiration (cohort%rdark) falls out of this as an unavoidable side
  ! effect, since it is a required input to LeafLayerPhotosynthesis itself.
  ! Non-leaf (stem/root) maintenance respiration is wired in via
  ! NonleafMaintenanceRespiration.
  !

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesConstantsMod,           only : nearzero
  use FatesConstantsMod,           only : umolC_to_kgC
  use FatesConstantsMod,           only : wm2_to_umolm2s
  use FatesCohortMod,              only : fates_cohort_type
  use FatesAllometryMod,           only : bleaf
  use FatesAllometryMod,           only : storage_fraction_of_target
  use PRTParametersMod,            only : prt_params
  use PRTGenericMod,               only : store_organ, sapw_organ, fnrt_organ, carbon12_element
  use FatesTestEnvironmentMod,     only : environment_type
  use FatesTestLightEnvMod,        only : light_env_type
  use FatesPlantRespPhotosynthMod, only : NonleafMaintenanceRespiration
  use LeafBiophysicsMod,           only : LeafLayerPhotosynthesis
  use LeafBiophysicsMod,           only : LowstorageMainRespReduction
  use LeafBiophysicsMod,           only : GetCanopyGasParameters
  use FatesTestLeafPhotoMod,       only : leaf_capacity_type
  use FatesTestLeafPhotoMod,       only : LeafLayerCapacity
  use FatesTestLeafPhotoMod,       only : LeafLayerSunShade
  use FatesTestLeafPhotoMod,       only : LeafLayerNitrogenScaling
  use FatesTestLeafPhotoMod,       only : ci_tol

  implicit none
  private

  integer,  parameter :: newton_max_iters = 10

  type, public :: cohort_phys_type

    private

    type(leaf_capacity_type), allocatable :: cap_z(:) ! per-leaf-layer photosynthetic capacity, stomatal conductance terms and dark respiration
    real(r8), allocatable :: nscaler_z(:)             ! per-leaf-layer nitrogen-scaling factor [0-1]

    real(r8) :: live_stem_n, live_croot_n   ! aboveground/belowground sapwood N [kgN/plant]
    real(r8) :: fnrt_n                      ! fine root N [kgN/plant]
    real(r8) :: maintresp_reduction_factor  ! storage-based throttle on maintenance respiration [0-1]

   contains

    procedure, public :: DailySetup
    procedure, public :: SubdailyStep
    procedure, public :: GrossAssimAndResp
    procedure, public :: LeafNetAssimSweep

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
    integer,                  intent(in)    :: nv  ! cohort's current number of occupied leaf layers

    if (allocated(this%cap_z)) then
      if (size(this%cap_z) /= nv) then
        deallocate(this%cap_z, this%nscaler_z)
        allocate(this%cap_z(nv), this%nscaler_z(nv))
      end if
    else
      allocate(this%cap_z(nv), this%nscaler_z(nv))
    end if

  end subroutine EnsureAllocated

  ! ==========================================================================

  subroutine DailySetup(this, cohort, pft, frac_store)
    !
    ! DESCRIPTION:
    ! Calculate a daily setup: storage-based maintenance-respiration throttle,
    ! sapwood/fine-root N, and per-leaf-layer nitrogen-scaling factor (nscaler_z depends on
    ! cumulative LAI above each layer, which only changes daily)
    !

    ! ARGUMENTS:
    class(cohort_phys_type), intent(inout) :: this       ! cohort physiology object
    type(fates_cohort_type), intent(in)    :: cohort     ! cohort to set up for today
    integer,                  intent(in)   :: pft        ! plant functional type index
    real(r8),                 intent(out)  :: frac_store ! ratio of storage carbon to target_leaf_c [-]

    ! LOCALS:
    real(r8) :: target_leaf_c  ! reference leaf biomass when fully flushed [kgC]
    real(r8) :: sapw_c, fnrt_c ! sapwood/fine root carbon [kgC/plant]

    ! storage-based maintenance respiration. The elongation factor
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

    ! aboveground/belowground sapwood N and fine root N
    sapw_c = cohort%prt%GetState(sapw_organ, carbon12_element)
    fnrt_c = cohort%prt%GetState(fnrt_organ, carbon12_element)
    this%live_stem_n  = sapw_c * prt_params%allom_agb_frac(pft) *              &
      prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(sapw_organ))
    this%live_croot_n = sapw_c * (1.0_r8 - prt_params%allom_agb_frac(pft)) *   &
      prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(sapw_organ))
    this%fnrt_n = fnrt_c * prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(fnrt_organ))

    call EnsureAllocated(this, cohort%nv)

    call LeafLayerNitrogenScaling(cohort%treelai, cohort%treesai,              &
      cohort%height, cohort%nv, pft, cohort%vcmax25top, this%nscaler_z)

  end subroutine DailySetup

  ! ==========================================================================

  subroutine SubdailyStep(this, cohort, pft, env, light_env, lnc_top, step_size, &
    n_photo_calls, n_bisection_calls, max_solve_iter, sum_solve_iter,        &
    gpp_tstep, rdark_tstep, nonleaf_mr_tstep)
    !
    ! DESCRIPTION:
    ! A single sub-daily carbon uptake step: refresh the temperature-dependent
    ! leaf photosynthetic capacity/dark-respiration rates from today's
    ! nscaler_z (from DailySetup) and this substep's env%tempk (matching production's
    ! per-timestep LeafLayerBiophysicalRates pattern), then prescribed absorbed
    ! PAR (light_env) -> leaf photosynthesis -> a per-individual GPP and leaf
    ! dark respiration, plus non-leaf (stem/root) maintenance respiration for
    ! the same substep.

    ! ARGUMENTS:
    class(cohort_phys_type), intent(inout) :: this              ! cohort physiology object
    type(fates_cohort_type), intent(inout) :: cohort            ! cohort to step
    integer,                 intent(in)    :: pft               ! plant functional type index
    type(environment_type),  intent(in)    :: env               ! prescribed atmospheric/soil boundary conditions
    type(light_env_type),    intent(in)    :: light_env         ! this substep's absorbed-PAR/LAI profile
    real(r8),                intent(in)    :: lnc_top           ! leaf N content at the canopy top [gN/m2 leaf]
    real(r8),                intent(in)    :: step_size         ! model time step [s]
    integer,                 intent(inout) :: n_photo_calls     ! running count of LeafLayerPhotosynthesis calls
    integer,                 intent(inout) :: n_bisection_calls ! running count of calls that fell back to CiBisection
    integer,                 intent(inout) :: max_solve_iter    ! running max Ci-solver iteration count
    integer,                 intent(inout) :: sum_solve_iter    ! running sum of Ci-solver iteration counts, for a mean once divided by n_photo_calls
    real(r8),                intent(out)   :: gpp_tstep         ! this substep's GPP [kgC/indiv/s]
    real(r8),                intent(out)   :: rdark_tstep       ! this substep's leaf dark respiration [kgC/indiv/s]
    real(r8),                intent(out)   :: nonleaf_mr_tstep  ! this substep's non-leaf (stem/root) maintenance respiration [kgC/indiv/s]

    ! LOCALS:
    integer  :: iv                          ! leaf-layer looping index
    real(r8) :: mm_kco2, mm_ko2, co2_cpoint ! Michaelis-Menten constants for CO2/O2, CO2 compensation point at env%tempk [Pa]
    integer  :: solve_iter_sun              ! Ci-solver iteration count, sunlit call - >newton_max_iters means it fell back to CiBisection
    integer  :: solve_iter_sha              ! Ci-solver iteration count, shaded call - >newton_max_iters means it fell back to CiBisection
    real(r8) :: psn_layer                   ! area-weighted mean gross photosynthesis for the layer [umolC/m2 leaf/s]
    real(r8) :: anet_layer                  ! area-weighted mean net photosynthesis for the layer [umolC/m2 leaf/s] - unused here, since a whole-plant carbon balance integrates GROSS assimilation and accounts for respiration separately
    real(r8) :: lai_layer                   ! this layer's total leaf area index [m2 leaf/m2 ground]
    real(r8) :: cohort_layer_eleaf_area     ! cohort's effective leaf area in the current layer [m2]
    real(r8) :: gpp_sum, rdark_sum          ! running per-substep GPP/dark-respiration accumulators [umolC/s]

    ! prescribed absorbed PAR -> leaf photosynthesis, one leaf layer at a time:
    ! this layer's temperature-dependent capacity/dark-respiration rates are
    ! refreshed (LeafLayerCapacity) from env%tempk/dayl_factor/t_growth/t_home/
    ! btran and today's nscaler_z (DailySetup), matching production's
    ! per-timestep LeafLayerBiophysicalRates pattern, then the layer's sunlit
    ! and shaded halves are evaluated at that capacity and area-weighted into a
    ! single per-layer rate (LeafLayerSunShade, which also converts
    ! parsun_z/parsha_z from W/m2 ground to umol photons/m2 leaf/s the same way
    ! FatesPlantRespPhotosynthMod's ConvertPar does). That per-layer rate is
    ! scaled up to a per-individual carbon flux below - reimplements
    ! ScaleLeafLayerFluxToCohort's formula, which is private to
    ! FatesPlantRespPhotosynthMod and so not directly reusable here.
    !
    
    ! The Michaelis-Menten constants/CO2
    call GetCanopyGasParameters(env%can_press, env%can_o2_ppress, env%tempk,   &
      mm_kco2, mm_ko2, co2_cpoint)

    gpp_sum = 0.0_r8
    rdark_sum = 0.0_r8
    do iv = 1, cohort%nv

      ! capacity is independent of sunlit vs. shaded, so it is computed once
      ! per layer and saved to cap_z(iv)
      call LeafLayerCapacity(pft, env%tempk, env%t_growth, env%t_home,         &
        this%nscaler_z(iv), env%dayl_factor, env%btran, cohort%vcmax25top,     &
        cohort%jmax25top, cohort%kp25top, lnc_top, this%cap_z(iv))

      call LeafLayerSunShade(pft, this%cap_z(iv), light_env%parsun_z(iv),      &
        light_env%parsha_z(iv), light_env%laisun_z(iv),                        &
        light_env%laisha_z(iv), env%tempk, env%can_press,                      &
        env%can_co2_ppress, env%can_o2_ppress, env%veg_esat, env%can_vpress,   &
        env%gb, mm_kco2, mm_ko2, co2_cpoint, psn_layer, anet_layer,            &
        lai_layer, solve_iter_sun, solve_iter_sha)

      n_photo_calls = n_photo_calls + 2
      if (solve_iter_sun > newton_max_iters) n_bisection_calls = n_bisection_calls + 1
      if (solve_iter_sha > newton_max_iters) n_bisection_calls = n_bisection_calls + 1
      max_solve_iter = max(max_solve_iter, solve_iter_sun, solve_iter_sha)
      sum_solve_iter = sum_solve_iter + solve_iter_sun + solve_iter_sha

      ! [umolC/m2 leaf/s] * [m2 leaf] -> [umolC/s], leaf area per this cohort
      ! in this layer = layer's exposed LAI * the cohort's own crown area
      cohort_layer_eleaf_area = lai_layer * cohort%c_area
      gpp_sum   = gpp_sum   + psn_layer * cohort_layer_eleaf_area
      rdark_sum = rdark_sum + this%cap_z(iv)%lmr * cohort_layer_eleaf_area

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
    laisun_z, laisha_z, maintresp_reduction_factor, step_size, gross_assim, total_resp)
    !
    ! DESCRIPTION:
    ! Whole-plant gross assimilation and total respiration (leaf dark + non-leaf
    ! maintenance) at an arbitrary already-attenuated PAR profile
    ! (parsun_z/parsha_z/laisun_z/laisha_z - typically from
    ! FatesTestLightEnvMod's AttenuateCanopy at a diagnostic incident PPFD, not
    ! this substep's real light_env%parsun_z etc.), using this cohort's current
    ! ("frozen") per-layer capacity (cap_z)
    ! exactly as SubdailyStep left them - NOT recomputed here, since a diagnostic
    ! sweep across incident PPFD is not supposed to advance any state. Unlike
    ! SubdailyStep, nothing here is written to cohort (gpp_tstep/rdark/
    ! livestem_mr/livecroot_mr/froot_mr/sym_nfix_tstep all stay untouched) or to
    ! this (cap_z), so there is nothing to restore afterward - this is a
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
    real(r8),                 intent(in) :: maintresp_reduction_factor  ! storage-based throttle on maintenance respiration [0-1]
    real(r8),                 intent(in) :: step_size      ! model time step [s] - passed through to NonleafMaintenanceRespiration; does not affect the carbon-flux rates read here (only its own sym_nfix_tstep output, which this diagnostic discards)
    real(r8),                 intent(out) :: gross_assim   ! whole-plant gross assimilation at this PAR profile [kgC/indiv/s]
    real(r8),                 intent(out) :: total_resp     ! whole-plant total respiration (leaf dark + non-leaf maintenance) at this PAR profile [kgC/indiv/s]

    ! LOCALS:
    real(r8) :: mm_kco2, mm_ko2, co2_cpoint ! Michaelis-Menten constants for CO2/O2, CO2 compensation point at env%tempk [Pa]
    integer  :: iv                         ! leaf-layer looping index
    real(r8) :: psn_layer                  ! area-weighted mean gross photosynthesis for the layer [umolC/m2 leaf/s]
    real(r8) :: anet_layer                 ! area-weighted mean net photosynthesis for the layer [umolC/m2 leaf/s] - unused here, since this diagnostic integrates GROSS assimilation and accounts for respiration separately below
    real(r8) :: lai_layer                  ! this layer's total leaf area index [m2 leaf/m2 ground]
    real(r8) :: cohort_layer_eleaf_area    ! cohort's effective leaf area in the current layer [m2]
    real(r8) :: gross_assim_sum, leaf_resp_sum ! running whole-plant accumulators [umolC/s]
    real(r8) :: livestem_mr, livecroot_mr, froot_mr, sym_nfix_tstep ! NonleafMaintenanceRespiration outputs (local - not cohort%livestem_mr etc.)

    call GetCanopyGasParameters(env%can_press, env%can_o2_ppress, env%tempk,    &
      mm_kco2, mm_ko2, co2_cpoint)

    gross_assim_sum = 0.0_r8
    leaf_resp_sum = 0.0_r8
    do iv = 1, cohort%nv

      call LeafLayerSunShade(pft, this%cap_z(iv), parsun_z(iv), parsha_z(iv),   &
        laisun_z(iv), laisha_z(iv), env%tempk, env%can_press,                  &
        env%can_co2_ppress, env%can_o2_ppress, env%veg_esat, env%can_vpress,   &
        env%gb, mm_kco2, mm_ko2, co2_cpoint, psn_layer, anet_layer, lai_layer)

      cohort_layer_eleaf_area = lai_layer * cohort%c_area
      gross_assim_sum = gross_assim_sum + psn_layer * cohort_layer_eleaf_area
      leaf_resp_sum   = leaf_resp_sum   + this%cap_z(iv)%lmr * cohort_layer_eleaf_area

    end do

    ! [umolC/s] -> [kgC/indiv/s]
    gross_assim = gross_assim_sum * umolC_to_kgC / cohort%n

    call NonleafMaintenanceRespiration(pft, env%tempk, env%nlevsoil,           &
      [env%t_soil], [env%rootfr_ft], this%live_stem_n, this%live_croot_n,      &
      this%fnrt_n, maintresp_reduction_factor, step_size,                 &
      livestem_mr, livecroot_mr, froot_mr, sym_nfix_tstep)

    total_resp = leaf_resp_sum * umolC_to_kgC * maintresp_reduction_factor / cohort%n + &
      livestem_mr + livecroot_mr + froot_mr

  end subroutine GrossAssimAndResp

  ! ==========================================================================

  subroutine LeafNetAssimSweep(this, pft, env, ppfd_values, iv, anet)
    !
    ! DESCRIPTION:
    ! Leaf-level light-response sweep at a single, fixed canopy position (iv),
    ! evaluated at that layer's already-"frozen" capacity (cap_z(iv), exactly
    ! as SubdailyStep left it). Unlike
    ! GrossAssimAndResp (whole-plant assimilation integrated over a caller-
    ! supplied, already-self-shaded canopy PAR profile), this sweeps a caller-
    ! supplied incident PPFD directly onto one leaf layer's photosynthesis, with
    ! no canopy attenuation or self-shading at all - this is Sterck et al.
    ! (2013)'s LCPleaf protocol: an IRGA light-response curve on a single leaf,
    ! giving Aarea (area-based net photosynthesis) as a function of PPFD
    ! incident directly on that leaf, as opposed to LCPplant (GrossAssimAndResp/
    ! LightResponseSweep), where the swept PPFD is whole-canopy incident and
    ! self-shading is resolved by the two-stream solver before reaching each
    ! layer.
    !
    ! iv=1 (top of crown) is the natural choice of canopy position for the
    ! caller to pass, since nscaler_z(1) is calibrated to (almost) zero
    ! cumulative LAI above it - this driver has not been checked against Sterck
    ! et al. (2013)'s actual IRGA measurement protocol (which leaf age/canopy
    ! position they sampled), so that choice is a documented assumption, not a
    ! confirmed match.

    ! ARGUMENTS:
    class(cohort_phys_type), intent(in)  :: this          ! cohort physiology object (read-only: today's already-computed capacity)
    integer,                  intent(in)  :: pft            ! plant functional type index
    type(environment_type),  intent(in)  :: env            ! prescribed atmospheric/soil boundary conditions
    real(r8),                 intent(in)  :: ppfd_values(:) ! incident PPFD values swept directly onto the leaf, no canopy attenuation [umol photons/m2 leaf/s]
    integer,                  intent(in)  :: iv             ! canopy layer to evaluate at (1 = top of crown - see header comment)
    real(r8),                 intent(out) :: anet(:)        ! leaf-level net photosynthesis (Aarea) at each swept PPFD [umolC/m2 leaf/s], size(ppfd_values)

    ! LOCALS:
    real(r8) :: mm_kco2, mm_ko2, co2_cpoint ! Michaelis-Menten constants for CO2/O2, CO2 compensation point at env%tempk [Pa]
    real(r8) :: agross     ! gross photosynthesis (unused diagnostic here - see anet) [umolC/m2 leaf/s]
    real(r8) :: gs         ! stomatal conductance (unused diagnostic)
    real(r8) :: ci         ! intracellular CO2 (unused diagnostic) [Pa]
    real(r8) :: c13disc    ! carbon-13 discrimination (unused diagnostic)
    integer  :: solve_iter ! Ci-solver iteration count (unused diagnostic here - see SubdailyStep for the tracked version)
    integer  :: ippfd      ! PPFD-sweep looping index

    call GetCanopyGasParameters(env%can_press, env%can_o2_ppress, env%tempk,    &
      mm_kco2, mm_ko2, co2_cpoint)

    do ippfd = 1, size(ppfd_values)
      call LeafLayerPhotosynthesis(ppfd_values(ippfd), pft,                     &
        this%cap_z(iv)%vcmax, this%cap_z(iv)%jmax, this%cap_z(iv)%kp,          &
        this%cap_z(iv)%gs0, this%cap_z(iv)%gs1, this%cap_z(iv)%gs2,            &
        env%tempk, env%can_press, env%can_co2_ppress, env%can_o2_ppress,       &
        env%veg_esat, env%gb, env%can_vpress, mm_kco2, mm_ko2, co2_cpoint,     &
        this%cap_z(iv)%lmr, ci_tol, agross, gs, anet(ippfd), c13disc, ci,      &
        solve_iter)
    end do

  end subroutine LeafNetAssimSweep

end module FatesTestCohortPhysMod
