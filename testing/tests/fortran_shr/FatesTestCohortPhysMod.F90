module FatesTestCohortPhysMod
  !
  ! DESCRIPTION:
  ! Methods for calculating physics for a single cohort
  !
  
  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesConstantsMod,           only : nearzero
  use FatesConstantsMod,           only : umolC_to_kgC
  use FatesConstantsMod,           only : wm2_to_umolm2s
  use FatesConstantsMod,           only : g_per_kg
  use FatesConstantsMod,           only : ihard_season_decid
  use FatesConstantsMod,           only : ihard_stress_decid, isemi_stress_decid
  use FatesConstantsMod,           only : ievergreen
  use EDPftvarcon,                 only : EDPftvarcon_inst
  use EDParamsMod,                 only : nclmax, nlevleaf
  use EDParamsMod,                 only : GetNVegLayers
  use EDParamsMod,                 only : dlower_vai, dinc_vai
  use LeafBiophysicsMod,           only : DecayCoeffVcmax
  use FatesCohortMod,              only : fates_cohort_type
  use FatesAllometryMod,           only : bleaf, carea_allom, tree_lai_sai
  use FatesAllometryMod,           only : bfineroot
  use FatesAllometryMod,           only : storage_fraction_of_target
  use PRTParametersMod,            only : prt_params
  use PRTGenericMod,               only : store_organ, sapw_organ, fnrt_organ, carbon12_element, leaf_organ
  use FatesTestEnvironmentMod,     only : environment_type
  use FatesTestLightEnvMod,        only : light_env_type
  use FatesPlantRespPhotosynthMod, only : NonleafMaintenanceRespiration
  use LeafBiophysicsMod,           only : LeafLayerPhotosynthesis
  use LeafBiophysicsMod,           only : LowstorageMainRespReduction
  use LeafBiophysicsMod,           only : GetCanopyGasParameters
  use FatesTestLeafPhotoMod,       only : leaf_capacity_type
  use FatesTestLeafPhotoMod,       only : LeafLayerCapacity
  use FatesTestLeafPhotoMod,       only : LeafLayerSunShade
  use FatesTestLeafPhotoMod,       only : LeafLayerVerticalScaling
  use FatesTestLeafPhotoMod,       only : ci_tol
  
  implicit none
  private

  integer,  parameter :: newton_max_iters = 10
  real(r8), parameter :: decid_leaf_long_max = 1.0_r8

  type, public :: cohort_phys_type

    private

    type(leaf_capacity_type), allocatable :: cap_z(:)     ! per-leaf-layer photosynthetic capacity, stomatal conductance terms and dark respiration
    real(r8),                 allocatable :: nscaler_z(:) ! per-leaf-layer nitrogen-scaling factor [0-1]
    real(r8),                 allocatable :: rdark_scaler_z(:) ! per-leaf-layer respiration-scaling factor [0-1]
    real(r8)                              :: live_stem_n  ! aboveground sawpood N [kgN/plant]
    real(r8)                              :: live_croot_n ! belowground sapwood N [kgN/plant]
    real(r8)                              :: fnrt_n       ! fine root N [kgN/plant]
    real(r8)                              :: maintresp_reduction_factor  ! storage-based factor on maintenance respiration [0-1]

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
        deallocate(this%cap_z, this%nscaler_z, this%rdark_scaler_z)
        allocate(this%cap_z(nv), this%nscaler_z(nv), this%rdark_scaler_z(nv))
      end if
    else
      allocate(this%cap_z(nv), this%nscaler_z(nv), this%rdark_scaler_z(nv))
    end if

  end subroutine EnsureAllocated

  ! ==========================================================================

  subroutine DailySetup(this, cohort, pft, lai_above, frac_store)
    !
    ! DESCRIPTION:
    ! Calculate a daily setup: storage-based maintenance-respiration factor,
    ! sapwood/fine-root N, and the per-leaf-layer vertical-scaling factors
    ! (nscaler_z/rdark_scaler_z depend on cumulative LAI above each layer,
    ! which only changes daily)
    !

    ! ARGUMENTS:
    class(cohort_phys_type), intent(inout) :: this       ! cohort physiology object
    type(fates_cohort_type), intent(in)    :: cohort     ! cohort to set up for today
    integer,                 intent(in)    :: pft        ! plant functional type index
    real(r8),                intent(in)    :: lai_above  ! lai above the cohort [m2/m2]
    real(r8),                intent(out)   :: frac_store ! ratio of storage carbon to target_leaf_c [-]

    ! LOCALS:
    real(r8) :: target_leaf_c  ! reference leaf biomass when fully flushed [kgC]
    real(r8) :: sapw_c         ! sawpood carbon [kgC/plant
    real(r8) :: fnrt_c         ! fine root carbon [kgC/plant]
    
    ! CONSTANTS:
    real(r8) :: elongation_factor = 1.0_r8 ! elongation factor, set to 1.0 because to get target_bleaf, we assume fully flushed leaves
    
    ! calculate storage-based maintenance respiration
    call bleaf(cohort%dbh, pft, cohort%crowndamage, cohort%canopy_trim, &
      elongation_factor, target_leaf_c)
    call storage_fraction_of_target(target_leaf_c,                    &
      cohort%prt%GetState(store_organ, carbon12_element), frac_store)
    call LowstorageMainRespReduction(frac_store, pft, this%maintresp_reduction_factor)

    ! calculate aboveground/belowground sapwood N and fine root N
    sapw_c = cohort%prt%GetState(sapw_organ, carbon12_element)
    fnrt_c = cohort%prt%GetState(fnrt_organ, carbon12_element)
    this%live_stem_n  = sapw_c * prt_params%allom_agb_frac(pft) *              &
      prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(sapw_organ))
    this%live_croot_n = sapw_c * (1.0_r8 - prt_params%allom_agb_frac(pft)) *   &
      prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(sapw_organ))
    this%fnrt_n = fnrt_c * prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(fnrt_organ))

    call EnsureAllocated(this, cohort%nv)

    call LeafLayerVerticalScaling(cohort%treelai, cohort%treesai,              &
      cohort%height, cohort%nv, pft, cohort%vcmax25top, lai_above,             &
      this%nscaler_z, this%rdark_scaler_z)

  end subroutine DailySetup

  ! ==========================================================================

  subroutine RefreshCapacity(this, cohort, pft, env, lnc_top)
    !
    ! DESCRIPTION:
    ! Refresh every leaf layer's temperature-dependent capacity and dark
    ! respiration from today's vertical-scaling factors and the current env%tempk

    ! ARGUMENTS:
    class(cohort_phys_type), intent(inout) :: this    ! cohort physiology object
    type(fates_cohort_type), intent(in)    :: cohort  ! cohort to refresh
    integer,                 intent(in)    :: pft     ! plant functional type index
    type(environment_type),  intent(in)    :: env     ! prescribed atmospheric/soil boundary conditions
    real(r8),                intent(in)    :: lnc_top ! leaf N content at the canopy top [gN/m2 leaf]

    ! LOCALS:
    integer :: iv ! leaf-layer looping index

    ! capacity is independent of sunlit vs. shaded, so it is computed once
    ! per layer and saved to cap_z(iv)
    do iv = 1, cohort%nv
      call LeafLayerCapacity(pft, env%tempk, env%t_growth, env%t_home,         &
        this%nscaler_z(iv), this%rdark_scaler_z(iv), env%dayl_factor,          &
        env%btran, cohort%vcmax25top, cohort%jmax25top, cohort%kp25top,        &
        lnc_top, this%cap_z(iv))
    end do

  end subroutine RefreshCapacity

  ! ==========================================================================

  subroutine IntegrateLeafLayers(cap_z, cohort, pft, env, parsun_z, parsha_z,  &
    laisun_z, laisha_z, gross_assim_sum, leaf_resp_sum, n_photo_calls,         &
    n_bisection_calls, max_solve_iter, sum_solve_iter, anet_z)
    !
    ! DESCRIPTION:
    ! Integrate leaf photosynthesis down this cohort's leaf layers at an
    ! already-attenuated PAR profile, at the capacity held in cap_z. Reads
    ! cap_z rather than refreshing it, so a diagnostic sweep advances no state
    ! (see RefreshCapacity).
    !
    ! The solver-diagnostic counters are optional and accumulated only when
    ! present.

    ! ARGUMENTS:
    type(leaf_capacity_type), intent(in)  :: cap_z(:)        ! per-leaf-layer capacity and dark respiration
    type(fates_cohort_type),  intent(in)  :: cohort          ! cohort to integrate (read-only: c_area, nv)
    integer,                  intent(in)  :: pft             ! plant functional type index
    type(environment_type),   intent(in)  :: env             ! prescribed atmospheric/soil boundary conditions
    real(r8),                 intent(in)  :: parsun_z(:)     ! absorbed PAR, sunlit, per leaf layer [W/m2 crown footprint]
    real(r8),                 intent(in)  :: parsha_z(:)     ! absorbed PAR, shaded, per leaf layer [W/m2 crown footprint]
    real(r8),                 intent(in)  :: laisun_z(:)     ! sunlit LAI per leaf layer [m2 leaf/m2 crown footprint]
    real(r8),                 intent(in)  :: laisha_z(:)     ! shaded LAI per leaf layer [m2 leaf/m2 crown footprint]
    real(r8),                 intent(out) :: gross_assim_sum ! whole-cohort gross assimilation [umolC/s]
    real(r8),                 intent(out) :: leaf_resp_sum   ! whole-cohort leaf dark respiration, before maintresp_reduction_factor [umolC/s]
    integer,  intent(inout), optional :: n_photo_calls     ! running count of LeafLayerPhotosynthesis calls
    integer,  intent(inout), optional :: n_bisection_calls ! running count of calls that fell back to CiBisection
    integer,  intent(inout), optional :: max_solve_iter    ! running max Ci-solver iteration count
    integer,  intent(inout), optional :: sum_solve_iter    ! running sum of Ci-solver iteration counts
    real(r8), intent(out),   optional :: anet_z(:)         ! net assimilation by leaf layer

    ! LOCALS:
    integer  :: iv                          ! leaf-layer looping index
    real(r8) :: mm_kco2, mm_ko2, co2_cpoint ! Michaelis-Menten constants for CO2/O2, CO2 compensation point at env%tempk [Pa]
    integer  :: solve_iter_sun              ! Ci-solver iteration count, sunlit call
    integer  :: solve_iter_sha              ! Ci-solver iteration count, shaded call
    real(r8) :: psn_layer                   ! area-weighted mean gross photosynthesis for a layer [umolC/m2 leaf/s]
    real(r8) :: anet_layer                  ! area-weighted mean net photosynthesis for a layer [umolC/m2 leaf/s]
    real(r8) :: lai_layer                   ! individual layer total leaf area index [m2 leaf/m2 crown footprint]
    real(r8) :: cohort_layer_eleaf_area     ! cohort's effective leaf area in the current layer [m2]

    call GetCanopyGasParameters(env%can_press, env%can_o2_ppress, env%tempk,   &
      mm_kco2, mm_ko2, co2_cpoint)
    
    if (present(anet_z)) anet_z(:) = 0.0_r8
    gross_assim_sum = 0.0_r8
    leaf_resp_sum = 0.0_r8
    do iv = 1, cohort%nv

      call LeafLayerSunShade(pft, cap_z(iv), parsun_z(iv), parsha_z(iv),       &
        laisun_z(iv), laisha_z(iv), env%tempk, env%can_press,                  &
        env%can_co2_ppress, env%can_o2_ppress, env%veg_esat, env%can_vpress,   &
        env%gb, mm_kco2, mm_ko2, co2_cpoint, psn_layer, anet_layer,            &
        lai_layer, solve_iter_sun, solve_iter_sha)

      if (present(n_photo_calls)) n_photo_calls = n_photo_calls + 2
      if (present(n_bisection_calls)) then
        if (solve_iter_sun > newton_max_iters) n_bisection_calls = n_bisection_calls + 1
        if (solve_iter_sha > newton_max_iters) n_bisection_calls = n_bisection_calls + 1
      end if
      if (present(max_solve_iter)) max_solve_iter = max(max_solve_iter, solve_iter_sun, solve_iter_sha)
      if (present(sum_solve_iter)) sum_solve_iter = sum_solve_iter + solve_iter_sun + solve_iter_sha

      ! [umolC/m2 leaf/s] * [m2 leaf] -> [umolC/s], leaf area per this cohort
      ! in this layer = layer's exposed LAI * the cohort's own crown area
      cohort_layer_eleaf_area = lai_layer * cohort%c_area
      gross_assim_sum = gross_assim_sum + psn_layer * cohort_layer_eleaf_area
      leaf_resp_sum = leaf_resp_sum + cap_z(iv)%lmr * cohort_layer_eleaf_area
      if (present(anet_z)) anet_z(iv) = anet_layer

    end do

  end subroutine IntegrateLeafLayers

  ! ==========================================================================

  subroutine SubdailyStep(this, cohort, pft, env, light_env, lnc_top, step_size, &
    n_photo_calls, n_bisection_calls, max_solve_iter, sum_solve_iter,        &
    gpp_tstep, rdark_tstep, nonleaf_mr_tstep)
    !
    ! DESCRIPTION:
    ! A single sub-daily carbon uptake step: refresh the temperature-dependent
    ! leaf photosynthetic capacity/dark-respiration rates, then prescribed absorbed
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
    real(r8) :: anet_z(nlevleaf)   ! per-leaf-layer net assimilation [umolC/m2 leaf/s]
    real(r8) :: gpp_sum, rdark_sum ! running per-substep GPP/dark-respiration accumulators [umolC/s]

    call RefreshCapacity(this, cohort, pft, env, lnc_top)

    call IntegrateLeafLayers(this%cap_z, cohort, pft, env, light_env%parsun_z, &
      light_env%parsha_z, light_env%laisun_z, light_env%laisha_z, gpp_sum,     &
      rdark_sum, n_photo_calls, n_bisection_calls, max_solve_iter,             &
      sum_solve_iter, anet_z)
      
    ! per-leaf-layer net uptake, for annual canopy trim 
    cohort%ts_net_uptake(:) = anet_z(:) * umolC_to_kgC * step_size

    ! [umolC/s] -> [kgC/indiv/s]
    gpp_tstep = gpp_sum * umolC_to_kgC / cohort%n
    rdark_tstep = rdark_sum * umolC_to_kgC * this%maintresp_reduction_factor / cohort%n
    cohort%gpp_tstep = gpp_tstep
    cohort%rdark = rdark_tstep

    ! non-leaf (stem/root) maintenance respiration
    call NonleafMaintenanceRespiration(pft, env%tempk, env%nlevsoil,           &
      [env%t_soil], [env%rootfr_ft], this%live_stem_n, this%live_croot_n,      &
      this%fnrt_n, this%maintresp_reduction_factor, step_size,                 &
      cohort%livestem_mr, cohort%livecroot_mr, cohort%froot_mr,                &
      cohort%sym_nfix_tstep)
    nonleaf_mr_tstep = cohort%livestem_mr + cohort%livecroot_mr + cohort%froot_mr

  end subroutine SubdailyStep

  ! ==========================================================================

  subroutine GrossAssimAndResp(this, cohort, pft, env, parsun_z, parsha_z,   &
    laisun_z, laisha_z, maintresp_reduction_factor, step_size, gross_assim,  &
    leaf_resp, total_resp)
    !
    ! DESCRIPTION:
    ! Whole-plant gross assimilation and total respiration (leaf dark + non-leaf
    ! maintenance) at an arbitrary already-attenuated PAR profile, using this cohort's 
    ! current per-layer capacity (cap_z), not recomputed here.
    ! Nothing is saved to any objects here.

    ! ARGUMENTS:
    class(cohort_phys_type),  intent(in)  :: this        ! cohort physiology object 
    type(environment_type),   intent(in)  :: env         ! prescribed atmospheric/soil boundary conditions
    type(fates_cohort_type),  intent(in)  :: cohort      ! cohort to diagnose (read-only: c_area, n, nv)    
    real(r8),                 intent(in)  :: parsun_z(:) ! absorbed PAR, sunlit, per leaf layer [W/m2 crown footprint]
    real(r8),                 intent(in)  :: parsha_z(:) ! absorbed PAR, shaded, per leaf layer [W/m2 crown footprint]
    real(r8),                 intent(in)  :: laisun_z(:) ! sunlit LAI per leaf layer [m2 leaf/m2 crown footprint]
    real(r8),                 intent(in)  :: laisha_z(:) ! shaded LAI per leaf layer [m2 leaf/m2 crown footprint]
    real(r8),                 intent(in)  :: step_size   ! model time step [s]
    real(r8),                 intent(out) :: gross_assim ! whole-plant gross assimilation at this PAR profile [kgC/indiv/s]
    real(r8),                 intent(out) :: leaf_resp   ! leaf dark respiration at this PAR profile [kgC/indiv/s]
    real(r8),                 intent(out) :: total_resp  ! whole-plant total respiration (leaf dark + non-leaf maintenance) at this PAR profile [kgC/indiv/s]
    real(r8),                 intent(in)  :: maintresp_reduction_factor  ! storage-based factor on maintenance respiration [0-1]
    integer,                  intent(in)  :: pft         ! plant functional type index

    ! LOCALS:
    real(r8) :: gross_assim_sum, leaf_resp_sum ! running whole-plant accumulators [umolC/s]
    real(r8) :: livestem_mr, livecroot_mr, froot_mr, sym_nfix_tstep ! NonleafMaintenanceRespiration outputs

    call IntegrateLeafLayers(this%cap_z, cohort, pft, env, parsun_z, parsha_z, &
      laisun_z, laisha_z, gross_assim_sum, leaf_resp_sum)

    ! [umolC/s] -> [kgC/indiv/s]
    gross_assim = gross_assim_sum * umolC_to_kgC / cohort%n

    call NonleafMaintenanceRespiration(pft, env%tempk, env%nlevsoil,           &
      [env%t_soil], [env%rootfr_ft], this%live_stem_n, this%live_croot_n,      &
      this%fnrt_n, maintresp_reduction_factor, step_size,                 &
      livestem_mr, livecroot_mr, froot_mr, sym_nfix_tstep)

    leaf_resp  = leaf_resp_sum * umolC_to_kgC * maintresp_reduction_factor / cohort%n
    total_resp = leaf_resp + livestem_mr + livecroot_mr + froot_mr

  end subroutine GrossAssimAndResp

  ! ==========================================================================

  subroutine LeafNetAssimSweep(this, pft, env, ppfd_values, iv, anet)
    !
    ! DESCRIPTION:
    ! Leaf-level light-response sweep at a single, fixed canopy position (iv),
    ! evaluated at that layer's current leaf capacity
    !

    ! ARGUMENTS:
    class(cohort_phys_type), intent(in)  :: this           ! cohort physiology object 
    integer,                 intent(in)  :: pft            ! plant functional type index
    type(environment_type),  intent(in)  :: env            ! prescribed atmospheric/soil boundary conditions
    real(r8),                intent(in)  :: ppfd_values(:) ! incident PPFD values swept directly onto the leaf, no canopy attenuation [umol photons/m2 leaf/s]
    integer,                 intent(in)  :: iv             ! canopy layer to evaluate at
    real(r8),                intent(out) :: anet(:)        ! leaf-level net photosynthesis (Aarea) at each swept PPFD [umolC/m2 leaf/s]

    ! LOCALS:
    real(r8) :: mm_kco2    ! Michaelis-Menten constants for CO2 at env%tempk [Pa]
    real(r8) :: mm_ko2     ! Michaelis-Menten constants for O2 at env%tempk [Pa]
    real(r8) :: co2_cpoint ! CO2 compensation point at env%tempk [Pa]
    real(r8) :: agross     ! gross photosynthesis [umolC/m2 leaf/s]
    real(r8) :: gs         ! stomatal conductance 
    real(r8) :: ci         ! intracellular CO2 [Pa]
    real(r8) :: c13disc    ! carbon-13 discrimination
    integer  :: solve_iter ! Ci-solver iteration count
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
  
  ! ==========================================================================
  
end module FatesTestCohortPhysMod
