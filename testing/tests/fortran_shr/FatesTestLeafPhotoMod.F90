module FatesTestLeafPhotoMod
  !
  ! DESCRIPTION:
  ! Stateless leaf-level photosynthesis evaluator for standalone test drivers -
  ! wraps the current production call sequence (GetCanopyGasParameters ->
  ! LeafLayerBiophysicalRates -> LeafLayerMaintenanceRespiration_Ryan_1991 ->
  ! LeafLayerPhotosynthesis) so a test can evaluate leaf gas exchange at any
  ! prescribed, arbitrary driver conditions. Every argument to
  ! EvaluateLeafPhotosynthesis below is recomputed from scratch on every call,
  ! since a leaf-level sensitivity sweep needs to vary things (temperature,
  ! nscaler, btran) that cohort_phys_type's GrossAssimAndResp/LeafNetAssimSweep
  ! hold fixed at whatever a cohort's daily setup already left them.
  !
  ! FatesTestCohortPhysMod's SubdailyStep also uses this call sequence, once
  ! per leaf layer per sunlit/shaded category - it calls
  ! EvaluateLeafPhotosynthesis too rather than reimplementing the sequence,
  ! passing the cohort's own (possibly acclimated) vcmax25top/jmax25top/
  ! kp25top instead of the flat PFT defaults test_LeafLevelPhoto.F90 passes,
  ! and reading back the optional vcmax/jmax/kp/gs0/gs1/gs2/lmr outputs so it
  ! can cache them for GrossAssimAndResp/LeafNetAssimSweep to reuse "frozen"
  ! later - those two evaluate LeafLayerPhotosynthesis directly at that cached
  ! capacity, with no BiophysicalRates/MaintenanceRespiration recompute, so
  ! they are a narrower operation this module does not wrap.
  !
  ! That capacity half of the sequence (LeafLayerBiophysicalRates ->
  ! LeafLayerMaintenanceRespiration_Ryan_1991) is available on its own as
  ! LeafLayerCapacity, returning a leaf_capacity_type - EvaluateLeafPhotosynthesis
  ! is itself written as LeafLayerCapacity followed by LeafLayerPhotosynthesis, so
  ! there is one implementation of that chain rather than one per caller. A caller
  ! evaluating the same leaf layer more than once at a single capacity (sunlit and
  ! shaded, or a diagnostic sweep across incident PPFD) can hold the
  ! leaf_capacity_type and drive LeafLayerPhotosynthesis directly, which is what
  ! makes the "frozen capacity" distinction visible at the call site instead of
  ! implied by which routine was called.
  !

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : wm2_to_umolm2s
  use LeafBiophysicsMod, only : GetCanopyGasParameters
  use LeafBiophysicsMod, only : LeafLayerBiophysicalRates
  use LeafBiophysicsMod, only : LeafLayerMaintenanceRespiration_Ryan_1991
  use LeafBiophysicsMod, only : LeafLayerPhotosynthesis
  use LeafBiophysicsMod, only : DecayCoeffVcmax
  use FatesAllometryMod, only : VegAreaLayer
  use PRTParametersMod,  only : prt_params
  use PRTGenericMod,     only : leaf_organ

  implicit none
  private

  ! convergence tolerance for intracellular CO2 [Pa] - matches production's
  ! own use. Public so FatesTestCohortPhysMod's GrossAssimAndResp/
  ! LeafNetAssimSweep (which call LeafLayerPhotosynthesis directly, at
  ! already-cached capacity, rather than through EvaluateLeafPhotosynthesis)
  ! can share this same value instead of redeclaring it.
  real(r8), public, parameter :: ci_tol = 0.5_r8

  ! minimum leaf area to bother solving photosynthesis for [m2] - matches
  ! ConvertPar in FatesPlantRespPhotosynthMod. Public so every caller
  ! converting a per-ground-area absorbed PAR to a per-leaf-area PPFD shares
  ! one threshold rather than redeclaring or guessing it.
  real(r8), public, parameter :: min_la_to_solve = 1.0e-10_r8

  ! ------------------------------------------------------------------------------------
  ! One leaf layer's photosynthetic capacity and dark respiration, evaluated at a
  ! specific set of leaf conditions - everything LeafLayerPhotosynthesis needs beyond
  ! the instantaneous driver state (absorbed PAR, temperature, humidity, and the gas
  ! partial pressures/Michaelis-Menten constants).
  !
  ! Bundled into a type because callers split into two groups that must stay
  ! distinguishable at the call site: those recomputing capacity from current
  ! conditions on every evaluation (via LeafLayerCapacity), and those deliberately
  ! reusing a previously computed ("frozen") capacity so that a diagnostic sweep
  ! advances no state (FatesTestCohortPhysMod's GrossAssimAndResp/LeafNetAssimSweep).
  !
  ! The gs0/gs1/gs2 comments follow LeafBiophysicsMod's own wording for the
  ! corresponding LeafLayerBiophysicalRates outputs. Production states units for gs0
  ! and gs2 but leaves gs1 dimensionally undocumented, so none is invented for it here
  ! ------------------------------------------------------------------------------------
  type, public :: leaf_capacity_type
     real(r8) :: vcmax ! maximum rate of carboxylation [umolCO2/m2 leaf/s]
     real(r8) :: jmax  ! maximum electron transport rate [umol electrons/m2 leaf/s]
     real(r8) :: kp    ! initial slope of the CO2 response curve, C4 plants [umol/m2 leaf/s]
     real(r8) :: gs0   ! effective stomatal conductance intercept [umol H2O/m2 leaf/s]
     real(r8) :: gs1   ! effective stomatal conductance slope, possibly multiplied by btran (units not stated in LeafBiophysicsMod)
     real(r8) :: gs2   ! alternative btran term applied to the whole non-intercept side of the Medlyn conductance equation - either 1 or btran [-]
     real(r8) :: lmr   ! leaf maintenance (dark) respiration [umolCO2/m2 leaf/s]
  end type leaf_capacity_type

  public :: LeafLayerCapacity
  public :: LeafLayerSunShade
  public :: EvaluateLeafPhotosynthesis
  public :: LeafNitrogenContent
  public :: LeafLayerNitrogenScaling
  public :: LayerParPerLeafArea
  public :: SunlitFraction

contains

  ! ==========================================================================

  function LayerParPerLeafArea(par_z, lai_z) result(par_abs)
    !
    ! DESCRIPTION:
    ! Converts one leaf layer's absorbed PAR from per unit GROUND area
    ! [W/m2 ground] to the per unit LEAF area PPFD [umol photons/m2 leaf/s]
    ! that LeafLayerPhotosynthesis expects - the same conversion
    ! FatesPlantRespPhotosynthMod's ConvertPar performs in production.
    !
    ! Call this once per sunlit/shaded category per layer, passing that
    ! category's own par/lai pair (parsun_z/laisun_z, or parsha_z/laisha_z).
    ! Returns zero rather than dividing when the layer holds effectively no
    ! leaf area of that category, or receives effectively no light - the
    ! guard that keeps an empty sunlit fraction at the bottom of a deep
    ! canopy from producing a division by (near) zero.

    ! ARGUMENTS:
    real(r8), intent(in) :: par_z   ! this layer/category's absorbed PAR [W/m2 ground]
    real(r8), intent(in) :: lai_z   ! this layer/category's leaf area index [m2 leaf/m2 ground]
    real(r8)             :: par_abs ! absorbed PAR per unit leaf area [umol photons/m2 leaf/s]

    if (lai_z > min_la_to_solve .and. par_z > nearzero) then
      par_abs = par_z/lai_z*wm2_to_umolm2s
    else
      par_abs = 0.0_r8
    end if

  end function LayerParPerLeafArea

  ! ==========================================================================

  function SunlitFraction(laisun_z, laisha_z) result(fsun)
    !
    ! DESCRIPTION:
    ! Sunlit fraction of one leaf layer's leaf area, for area-weighting that
    ! layer's sunlit and shaded photosynthesis into a single per-layer rate.
    ! Returns zero for a layer holding effectively no leaf area at all (rather
    ! than dividing by it), in which case the weighting is moot - the layer
    ! contributes nothing either way.

    ! ARGUMENTS:
    real(r8), intent(in) :: laisun_z ! this layer's sunlit leaf area index [m2 leaf/m2 ground]
    real(r8), intent(in) :: laisha_z ! this layer's shaded leaf area index [m2 leaf/m2 ground]
    real(r8)             :: fsun     ! sunlit fraction of this layer's leaf area [0-1]

    if (laisun_z + laisha_z > nearzero) then
      fsun = laisun_z / (laisun_z + laisha_z)
    else
      fsun = 0.0_r8
    end if

  end function SunlitFraction

  ! ==========================================================================

  subroutine LeafLayerNitrogenScaling(treelai, treesai, height, nv, pft,        &
    vcmax25top, nscaler_z)
    !
    ! DESCRIPTION:
    ! Per-leaf-layer nitrogen-scaling factor (nscaler), the vertical decay of
    ! photosynthetic capacity with cumulative leaf area above a layer - the
    ! nscaler argument EvaluateLeafPhotosynthesis takes. Each layer is
    ! evaluated at the cumulative LAI above its own MIDPOINT (lai_above +
    ! 0.5*elai_layer), not its top or bottom edge.
    !
    ! Depends only on canopy structure (treelai/treesai/height/nv) and the
    ! reference canopy-top capacity, so a caller with a prescribed canopy and
    ! no cohort can use it identically to one driving a real cohort - this is
    ! shared by FatesTestCohortPhysMod's DailySetup (a cohort's own
    ! allometry-derived canopy, refreshed daily) and test_CanopyLevelPhoto.F90
    ! (a prescribed, fixed LAI), rather than either reimplementing it.
    !
    ! Snow depth is fixed at zero here, matching every current standalone test
    ! driver - none of them model snow, and VegAreaLayer's only use of it is
    ! the above-snow exposed fraction.

    ! ARGUMENTS:
    real(r8), intent(in)  :: treelai      ! in-crown leaf area index [m2 leaf/m2 crown footprint]
    real(r8), intent(in)  :: treesai      ! in-crown stem area index [m2 stem/m2 crown footprint]
    real(r8), intent(in)  :: height       ! plant/canopy height [m]
    integer,  intent(in)  :: nv           ! number of occupied leaf layers
    integer,  intent(in)  :: pft          ! plant functional type index
    real(r8), intent(in)  :: vcmax25top   ! reference (25C, canopy-top) maximum carboxylation rate [umol/m2/s]
    real(r8), intent(out) :: nscaler_z(:) ! per-leaf-layer nitrogen-scaling factor [0-1], first nv entries filled

    ! LOCALS:
    integer  :: iv                     ! leaf-layer looping index
    real(r8) :: vai_top, vai_bot       ! vegetation area index bounds of the current leaf layer
    real(r8) :: elai_layer, esai_layer ! exposed leaf/stem area index of the current leaf layer
    real(r8) :: cumulative_lai         ! LAI above the middle of the current leaf layer
    real(r8) :: lai_above              ! running LAI above the current leaf layer
    real(r8) :: kn                     ! nitrogen vertical-scaling decay coefficient

    real(r8), parameter :: snow_depth = 0.0_r8 ! no snow modeled in any standalone test driver [m]

    ! loop-invariant: depends only on pft and the reference canopy-top
    ! capacity, neither of which varies by layer
    kn = DecayCoeffVcmax(vcmax25top, prt_params%leafn_vert_scaler_coeff1(pft), &
      prt_params%leafn_vert_scaler_coeff2(pft))

    lai_above = 0.0_r8
    do iv = 1, nv
      call VegAreaLayer(treelai, treesai, height, iv, nv, pft, snow_depth,     &
        vai_top, vai_bot, elai_layer, esai_layer)
      cumulative_lai = lai_above + 0.5_r8*elai_layer
      lai_above = lai_above + elai_layer

      nscaler_z(iv) = exp(-kn*cumulative_lai)
    end do

  end subroutine LeafLayerNitrogenScaling

  ! ==========================================================================

  function LeafNitrogenContent(pft) result(lnc_top)
    !
    ! DESCRIPTION:
    ! Leaf nitrogen content at the canopy top [gN/m2 leaf], needed by
    ! LeafLayerMaintenanceRespiration_Ryan_1991. Shared by
    ! test_LeafLevelPhoto.F90 and test_SingleCohort.F90 for their own
    ! lnc_top. Requires PRTDerivedParams() to have already populated
    ! prt_params%organ_param_id (see test_LeafLevelPhoto.F90's setup
    ! sequence)

    ! ARGUMENTS:
    integer, intent(in) :: pft ! plant functional type index

    ! LOCALS:
    real(r8) :: lnc_top ! leaf N content at the canopy top [gN/m2 leaf]

    lnc_top = prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(leaf_organ)) / &
      prt_params%slatop(pft)

  end function LeafNitrogenContent

  ! ==========================================================================

  subroutine LeafLayerCapacity(pft, veg_tempk, t_growth, t_home, nscaler,         &
    dayl_factor, btran, vcmax25top, jmax25top, kp25top, lnc_top, cap)
    !
    ! DESCRIPTION:
    ! One leaf layer's photosynthetic capacity and dark respiration at the given
    ! leaf conditions - the LeafLayerBiophysicalRates ->
    ! LeafLayerMaintenanceRespiration_Ryan_1991 half of the production call
    ! sequence EvaluateLeafPhotosynthesis wraps. That routine is written in terms
    ! of this one (it calls this, then LeafLayerPhotosynthesis), so this is the
    ! single implementation of that chain rather than a parallel copy of it.
    !
    ! Separated out so a caller evaluating the same leaf layer more than once at
    ! one capacity - sunlit and shaded, or a diagnostic sweep across incident PPFD
    ! - can compute capacity once and hold it. Both routines called here are
    ! closed-form Arrhenius-type evaluations rather than solves, so this exists to
    ! let the frozen-capacity callers say so structurally (see
    ! leaf_capacity_type's comment), not to save their cost.
    !
    ! Requires ReadParameters, InitializeGlobals and PRTDerivedParams to have
    ! already been called, and lb_params' model switches to have been set by the
    ! calling driver - see EvaluateLeafPhotosynthesis's header comment for both
    ! dependencies.

    ! ARGUMENTS:
    integer,  intent(in)  :: pft         ! plant functional type index
    real(r8), intent(in)  :: veg_tempk   ! instantaneous leaf temperature [K]
    real(r8), intent(in)  :: t_growth    ! 10-day running-mean growth temperature [K]
    real(r8), intent(in)  :: t_home      ! long-term running-mean home temperature [K]
    real(r8), intent(in)  :: nscaler     ! leaf nitrogen vertical-scaling factor [0-1]
    real(r8), intent(in)  :: dayl_factor ! day-length photosynthetic-capacity acclimation factor [0-1]
    real(r8), intent(in)  :: btran       ! soil moisture stress factor [0-1]
    real(r8), intent(in)  :: vcmax25top  ! reference (25C, canopy-top) maximum carboxylation rate [umol/m2/s]
    real(r8), intent(in)  :: jmax25top   ! reference (25C, canopy-top) maximum electron transport rate [umol/m2/s]
    real(r8), intent(in)  :: kp25top     ! reference (25C, canopy-top) initial slope of C4 CO2 response [umol/m2/s]
    real(r8), intent(in)  :: lnc_top     ! leaf N content at the canopy top [gN/m2 leaf] (see LeafNitrogenContent)
    type(leaf_capacity_type), intent(out) :: cap ! this layer's capacity/dark respiration at these conditions

    call LeafLayerBiophysicalRates(pft, vcmax25top, jmax25top, kp25top, nscaler,   &
      veg_tempk, dayl_factor, t_growth, t_home, btran, cap%vcmax, cap%jmax,       &
      cap%kp, cap%gs0, cap%gs1, cap%gs2)

    call LeafLayerMaintenanceRespiration_Ryan_1991(lnc_top, nscaler, pft,          &
      veg_tempk, cap%lmr)

  end subroutine LeafLayerCapacity

  ! ==========================================================================

  subroutine LeafLayerSunShade(pft, cap, par_sun_z, par_sha_z, lai_sun_z,        &
    lai_sha_z, veg_tempk, can_press, can_co2_ppress, can_o2_ppress, veg_esat,    &
    can_vpress, gb, mm_kco2, mm_ko2, co2_cpoint, agross_layer, anet_layer,       &
    lai_layer, solve_iter_sun, solve_iter_sha)
    !
    ! DESCRIPTION:
    ! One leaf layer's sunlit and shaded photosynthesis at a given capacity,
    ! area-weighted by that layer's sunlit fraction into a single per-layer rate
    ! per unit LEAF area - the inner step every canopy-integrating caller shares,
    ! whatever it then does with the result.
    !
    ! Absorbed PAR arrives per unit GROUND area (parsun_z/parsha_z, straight off
    ! FatesTestLightEnvMod's AttenuateCanopy) and is converted per sunlit/shaded
    ! category by LayerParPerLeafArea, exactly as production's ConvertPar does.
    !
    ! Both agross_layer and anet_layer are returned because they fall out of the
    ! same two LeafLayerPhotosynthesis calls at no extra cost, and callers differ
    ! on which they want: a canopy-flux protocol integrates NET assimilation
    ! (leaf dark respiration included, non-leaf respiration excluded), while a
    ! whole-plant carbon-balance caller integrates GROSS assimilation and accounts
    ! for respiration separately. Returning both is what keeps that a call-site
    ! choice rather than a reason to write a second copy of this loop body.
    !
    ! lai_layer (this layer's total leaf area index) is returned as the caller's
    ! area weight, computed here as lai_sun_z + lai_sha_z so that every caller
    ! weights by the same quantity. Scaling it further - by crown area, to reach
    ! a per-individual flux - is left to the caller, since that step is what
    ! distinguishes a per-ground-area canopy flux from a per-individual one.
    !
    ! Capacity is passed in rather than derived: a caller refreshing capacity
    ! every evaluation calls LeafLayerCapacity first, while one deliberately
    ! holding capacity fixed (a diagnostic sweep that must advance no state)
    ! simply passes what it already had. That distinction stays visible at the
    ! call site - see leaf_capacity_type's comment.
    !
    ! The Michaelis-Menten constants/CO2 compensation point are likewise passed
    ! in rather than derived per layer: they depend only on can_press,
    ! can_o2_ppress and veg_tempk, none of which vary down a canopy within one
    ! evaluation, so the caller hoists GetCanopyGasParameters out of its layer
    ! loop.

    ! ARGUMENTS:
    integer,  intent(in)  :: pft            ! plant functional type index
    type(leaf_capacity_type), intent(in) :: cap ! this layer's capacity/dark respiration (see LeafLayerCapacity)
    real(r8), intent(in)  :: par_sun_z      ! absorbed PAR, sunlit leaves, this layer [W/m2 ground]
    real(r8), intent(in)  :: par_sha_z      ! absorbed PAR, shaded leaves, this layer [W/m2 ground]
    real(r8), intent(in)  :: lai_sun_z      ! sunlit leaf area index, this layer [m2 leaf/m2 ground]
    real(r8), intent(in)  :: lai_sha_z      ! shaded leaf area index, this layer [m2 leaf/m2 ground]
    real(r8), intent(in)  :: veg_tempk      ! leaf temperature [K]
    real(r8), intent(in)  :: can_press      ! air pressure at the leaf surface [Pa]
    real(r8), intent(in)  :: can_co2_ppress ! CO2 partial pressure at the leaf surface [Pa]
    real(r8), intent(in)  :: can_o2_ppress  ! O2 partial pressure at the leaf surface [Pa]
    real(r8), intent(in)  :: veg_esat       ! saturation vapor pressure at veg_tempk [Pa]
    real(r8), intent(in)  :: can_vpress     ! vapor pressure of the canopy air [Pa]
    real(r8), intent(in)  :: gb             ! leaf boundary layer conductance [umol/m2/s]
    real(r8), intent(in)  :: mm_kco2        ! Michaelis-Menten constant for CO2 at veg_tempk [Pa]
    real(r8), intent(in)  :: mm_ko2         ! Michaelis-Menten constant for O2 at veg_tempk [Pa]
    real(r8), intent(in)  :: co2_cpoint     ! CO2 compensation point at veg_tempk [Pa]
    real(r8), intent(out) :: agross_layer   ! area-weighted gross photosynthesis for this layer [umolC/m2 leaf/s]
    real(r8), intent(out) :: anet_layer     ! area-weighted net photosynthesis for this layer [umolC/m2 leaf/s]
    real(r8), intent(out) :: lai_layer      ! this layer's total leaf area index [m2 leaf/m2 ground]
    integer,  intent(out), optional :: solve_iter_sun ! Ci-solver iteration count, sunlit call, for callers tracking solver diagnostics
    integer,  intent(out), optional :: solve_iter_sha ! Ci-solver iteration count, shaded call, for callers tracking solver diagnostics

    ! LOCALS:
    real(r8) :: par_abs                  ! absorbed PAR per unit leaf area [umol photons/m2 leaf/s]
    real(r8) :: agross_sun, agross_sha   ! gross photosynthesis, sunlit/shaded [umolC/m2 leaf/s]
    real(r8) :: anet_sun, anet_sha       ! net photosynthesis, sunlit/shaded [umolC/m2 leaf/s]
    real(r8) :: gs_sun, gs_sha           ! stomatal conductance, sunlit/shaded (unused diagnostic here)
    real(r8) :: ci_sun, ci_sha           ! intracellular CO2, sunlit/shaded (unused diagnostic here)
    real(r8) :: c13disc_sun, c13disc_sha ! carbon-13 discrimination, sunlit/shaded (unused diagnostic here)
    integer  :: solve_iter_sun_l         ! Ci-solver iteration count, sunlit call
    integer  :: solve_iter_sha_l         ! Ci-solver iteration count, shaded call
    real(r8) :: fsun                     ! sunlit fraction of this layer's leaf area [0-1]

    par_abs = LayerParPerLeafArea(par_sun_z, lai_sun_z)
    call LeafLayerPhotosynthesis(par_abs, pft, cap%vcmax, cap%jmax, cap%kp,       &
      cap%gs0, cap%gs1, cap%gs2, veg_tempk, can_press, can_co2_ppress,           &
      can_o2_ppress, veg_esat, gb, can_vpress, mm_kco2, mm_ko2, co2_cpoint,      &
      cap%lmr, ci_tol, agross_sun, gs_sun, anet_sun, c13disc_sun, ci_sun,        &
      solve_iter_sun_l)

    par_abs = LayerParPerLeafArea(par_sha_z, lai_sha_z)
    call LeafLayerPhotosynthesis(par_abs, pft, cap%vcmax, cap%jmax, cap%kp,       &
      cap%gs0, cap%gs1, cap%gs2, veg_tempk, can_press, can_co2_ppress,           &
      can_o2_ppress, veg_esat, gb, can_vpress, mm_kco2, mm_ko2, co2_cpoint,      &
      cap%lmr, ci_tol, agross_sha, gs_sha, anet_sha, c13disc_sha, ci_sha,        &
      solve_iter_sha_l)

    fsun = SunlitFraction(lai_sun_z, lai_sha_z)
    agross_layer = fsun*agross_sun + (1.0_r8 - fsun)*agross_sha
    anet_layer   = fsun*anet_sun + (1.0_r8 - fsun)*anet_sha
    lai_layer    = lai_sun_z + lai_sha_z

    if (present(solve_iter_sun)) solve_iter_sun = solve_iter_sun_l
    if (present(solve_iter_sha)) solve_iter_sha = solve_iter_sha_l

  end subroutine LeafLayerSunShade

  ! ==========================================================================

  subroutine EvaluateLeafPhotosynthesis(pft, par_abs, veg_tempk, t_growth, t_home, &
    can_press, can_co2_ppress, can_o2_ppress, veg_esat, can_vpress, gb, nscaler,   &
    dayl_factor, btran, vcmax25top, jmax25top, kp25top, lnc_top, agross, anet, gs, ci, &
    vcmax, jmax, kp, gs0, gs1, gs2, lmr, solve_iter_out)
    !
    ! DESCRIPTION:
    ! Evaluates leaf-level photosynthesis at arbitrary prescribed driver
    ! conditions - the full current production call sequence
    ! (GetCanopyGasParameters -> LeafLayerBiophysicalRates ->
    ! LeafLayerMaintenanceRespiration_Ryan_1991 -> LeafLayerPhotosynthesis, the
    ! middle two reached via LeafLayerCapacity above),
    ! with no self-shading and no canopy structure: every argument is a
    ! single instantaneous scalar, and every biophysical rate is recomputed
    ! from these arguments alone on every call, with nothing cached between
    ! calls (see FatesTestCohortPhysMod.F90's SubdailyStep, which calls this
    ! twice per leaf layer - once sunlit, once shaded - and so recomputes
    ! GetCanopyGasParameters redundantly within a single substep; that is
    ! deliberate, not an oversight - it is a closed-form Arrhenius-type
    ! evaluation, not a solve, so the redundant calls are negligible next to
    ! LeafLayerPhotosynthesis's own Ci iteration).
    !
    ! vcmax25top/jmax25top/kp25top are the caller's choice of reference
    ! (25C, canopy-top) photosynthetic capacity - test_LeafLevelPhoto.F90
    ! passes the flat PFT default (EDPftvarcon_inst%vcmax25top(pft,1),
    ! param_derived%jmax25top/kp25top(pft,1)), while FatesTestCohortPhysMod's
    ! SubdailyStep passes a cohort's own, possibly N/light-acclimated,
    ! cohort%vcmax25top/jmax25top/kp25top.
    !
    ! vcmax/jmax/kp/gs0/gs1/gs2/lmr are optional outputs of the intermediate
    ! LeafLayerBiophysicalRates/LeafLayerMaintenanceRespiration_Ryan_1991
    ! results - test_LeafLevelPhoto.F90 has no use for them and leaves them
    ! absent, while SubdailyStep reads them back to cache as this leaf
    ! layer's "frozen" capacity for GrossAssimAndResp/LeafNetAssimSweep to
    ! reuse later. solve_iter is similarly optional, for SubdailyStep's own
    ! Ci-solver diagnostic counters (n_bisection_calls/max_solve_iter/
    ! sum_solve_iter) - test_LeafLevelPhoto.F90 has no use for it either.
    !
    ! Note this means LeafLayerBiophysicalRates and
    ! LeafLayerMaintenanceRespiration_Ryan_1991 (both closed-form, not
    ! iterative solves) run once per call here, so SubdailyStep calling this
    ! twice per leaf layer (sunlit, then shaded) recomputes them twice per
    ! layer where the original inlined version computed them once and reused
    ! the result for both - deliberate, not an oversight, since the actual
    ! Ci-solve inside LeafLayerPhotosynthesis (which still only runs once per
    ! sunlit/shaded category either way) dominates the per-call cost.
    !
    ! Requires ReadParameters, InitializeGlobals, and PRTDerivedParams to have
    ! already been called (see test_LeafLevelPhoto.F90's setup sequence, and
    ! LeafNitrogenContent's header comment for the specific dependency on
    ! PRTDerivedParams) - this evaluates production physics directly, with no
    ! fallback defaults of its own for anything those calls populate.
    !
    ! lb_params' model switches (electron_transport_model, stomatal_model,
    ! stomatal_assim_model, photo_tempsens_model) are HLM-namelist-controlled
    ! in production and so are not set by any call in this module - the
    ! calling driver must set them explicitly first (see
    ! test_LeafLevelPhoto.F90's setup, which matches test_SingleCohort.F90's
    ! choices for consistency across this test suite), or
    ! LeafLayerBiophysicalRates' temperature-acclimation branch will hit its
    ! unhandled case-default and abort the run.

    ! ARGUMENTS:
    integer,  intent(in)  :: pft             ! plant functional type index
    real(r8), intent(in)  :: par_abs         ! absorbed PAR per unit leaf area [umol photons/m2 leaf/s]
    real(r8), intent(in)  :: veg_tempk       ! instantaneous leaf temperature [K]
    real(r8), intent(in)  :: t_growth        ! 10-day running-mean growth temperature [K]
    real(r8), intent(in)  :: t_home          ! long-term running-mean home temperature [K]
    real(r8), intent(in)  :: can_press       ! air pressure at the leaf surface [Pa]
    real(r8), intent(in)  :: can_co2_ppress  ! CO2 partial pressure at the leaf surface [Pa]
    real(r8), intent(in)  :: can_o2_ppress   ! O2 partial pressure at the leaf surface [Pa]
    real(r8), intent(in)  :: veg_esat        ! saturation vapor pressure at veg_tempk [Pa]
    real(r8), intent(in)  :: can_vpress      ! vapor pressure of the canopy air [Pa]
    real(r8), intent(in)  :: gb              ! leaf boundary layer conductance [umol/m2/s]
    real(r8), intent(in)  :: nscaler         ! leaf nitrogen vertical-scaling factor [0-1]
    real(r8), intent(in)  :: dayl_factor     ! day-length photosynthetic-capacity acclimation factor [0-1]
    real(r8), intent(in)  :: btran           ! soil moisture stress factor [0-1]
    real(r8), intent(in)  :: vcmax25top      ! reference (25C, canopy-top) maximum carboxylation rate [umol/m2/s]
    real(r8), intent(in)  :: jmax25top       ! reference (25C, canopy-top) maximum electron transport rate [umol/m2/s]
    real(r8), intent(in)  :: kp25top         ! reference (25C, canopy-top) initial slope of C4 CO2 response [umol/m2/s]
    real(r8), intent(in)  :: lnc_top         ! leaf N content at the canopy top [gN/m2 leaf] (see LeafNitrogenContent)
    real(r8), intent(out) :: agross          ! gross photosynthesis [umolC/m2/s]
    real(r8), intent(out) :: anet            ! net photosynthesis [umolC/m2/s]
    real(r8), intent(out) :: gs              ! leaf stomatal conductance [umol H2O/m2/s]
    real(r8), intent(out) :: ci              ! intracellular CO2 [Pa]
    real(r8), intent(out), optional :: vcmax ! leaf biophysical capacity at these conditions, for callers that cache it (e.g. SubdailyStep)
    real(r8), intent(out), optional :: jmax  ! leaf biophysical capacity at these conditions, for callers that cache it
    real(r8), intent(out), optional :: kp    ! leaf biophysical capacity at these conditions, for callers that cache it
    real(r8), intent(out), optional :: gs0   ! stomatal conductance slope/intercept term, for callers that cache it
    real(r8), intent(out), optional :: gs1   ! stomatal conductance slope/intercept term, for callers that cache it
    real(r8), intent(out), optional :: gs2   ! stomatal conductance slope/intercept term, for callers that cache it
    real(r8), intent(out), optional :: lmr   ! leaf maintenance (dark) respiration [umolCO2/m2/s], for callers that cache it
    integer,  intent(out), optional :: solve_iter_out ! Ci-solver iteration count, for callers tracking solver diagnostics (e.g. SubdailyStep)

    ! LOCALS:
    real(r8) :: mm_kco2, mm_ko2, co2_cpoint ! Michaelis-Menten constants for CO2/O2, CO2 compensation point [Pa]
    type(leaf_capacity_type) :: cap         ! leaf biophysical capacity/dark respiration at these conditions
    real(r8) :: c13disc                     ! carbon-13 discrimination (unused diagnostic here)
    integer  :: solve_iter                  ! Ci-solver iteration count (unused diagnostic here)

    call GetCanopyGasParameters(can_press, can_o2_ppress, veg_tempk, mm_kco2, mm_ko2, co2_cpoint)

    call LeafLayerCapacity(pft, veg_tempk, t_growth, t_home, nscaler, dayl_factor, &
      btran, vcmax25top, jmax25top, kp25top, lnc_top, cap)

    call LeafLayerPhotosynthesis(par_abs, pft, cap%vcmax, cap%jmax, cap%kp,        &
      cap%gs0, cap%gs1, cap%gs2, veg_tempk, can_press, can_co2_ppress,            &
      can_o2_ppress, veg_esat, gb, can_vpress, mm_kco2, mm_ko2, co2_cpoint,       &
      cap%lmr, ci_tol, agross, gs, anet, c13disc, ci, solve_iter)

    if (present(vcmax)) vcmax = cap%vcmax
    if (present(jmax))  jmax  = cap%jmax
    if (present(kp))    kp    = cap%kp
    if (present(gs0))   gs0   = cap%gs0
    if (present(gs1))   gs1   = cap%gs1
    if (present(gs2))   gs2   = cap%gs2
    if (present(lmr))   lmr   = cap%lmr
    if (present(solve_iter_out)) solve_iter_out = solve_iter

  end subroutine EvaluateLeafPhotosynthesis

end module FatesTestLeafPhotoMod
