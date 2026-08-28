module FatesTestLeafPhotoMod
  !
  ! DESCRIPTION:
  ! Helper methods for running leaf-level photosynthesis
  !

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesConstantsMod,           only : nearzero
  use FatesConstantsMod,           only : wm2_to_umolm2s
  use FatesInterfaceTypesMod,      only : hlm_maintresp_leaf_model
  use FatesConstantsMod,           only : lmrmodel_ryan_1991
  use FatesConstantsMod,           only : lmrmodel_atkin_etal_2017
  use LeafBiophysicsMod,           only : GetCanopyGasParameters
  use LeafBiophysicsMod,           only : LeafLayerBiophysicalRates
  use LeafBiophysicsMod,           only : LeafLayerMaintenanceRespiration_Ryan_1991
  use LeafBiophysicsMod,           only : LeafLayerMaintenanceRespiration_Atkin_etal_2017
  use LeafBiophysicsMod,           only : LeafLayerPhotosynthesis
  use LeafBiophysicsMod,           only : DecayCoeffVcmax
  use FatesPlantRespPhotosynthMod, only : ConvertPar
  use FatesAllometryMod,           only : VegAreaLayer
  use EDPftvarcon,                 only : EDPftvarcon_inst
  use PRTParametersMod,            only : prt_params
  use PRTGenericMod,               only : leaf_organ

  implicit none
  private

  ! convergence tolerance for intracellular CO2 [Pa]
  real(r8), public, parameter :: ci_tol = 0.5_r8

  ! minimum leaf area to bother solving photosynthesis for [m2]
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
  ! ------------------------------------------------------------------------------------
  type, public :: leaf_capacity_type
     real(r8) :: vcmax ! maximum rate of carboxylation [umolCO2/m2 leaf/s]
     real(r8) :: jmax  ! maximum electron transport rate [umol electrons/m2 leaf/s]
     real(r8) :: kp    ! initial slope of the CO2 response curve, C4 plants [umol/m2 leaf/s]
     real(r8) :: gs0   ! effective stomatal conductance intercept [umol H2O/m2 leaf/s]
     real(r8) :: gs1   ! effective stomatal conductance slope
     real(r8) :: gs2   ! alternative btran term applied to the whole non-intercept side of the Medlyn conductance equation
     real(r8) :: lmr   ! leaf maintenance (dark) respiration [umolCO2/m2 leaf/s]
  end type leaf_capacity_type

  public :: LeafLayerCapacity
  public :: LeafLayerSunShade
  public :: EvaluateLeafPhotosynthesis
  public :: LeafNitrogenContent
  public :: LeafLayerVerticalScaling
  public :: SunlitFraction

contains

  ! ==========================================================================
  
  function LeafNitrogenContent(pft) result(lnc_top)
    !
    ! DESCRIPTION:
    ! Leaf nitrogen content at the canopy top [gN/m2 leaf], needed by both leaf
    ! maintenance respiration models. Requires PRTDerivedParams() to have
    ! already populated prt_params%organ_param_id

    ! ARGUMENTS:
    integer, intent(in) :: pft ! plant functional type index

    ! LOCALS:
    real(r8) :: lnc_top ! leaf N content at the canopy top [gN/m2 leaf]

    lnc_top = prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(leaf_organ)) / &
      prt_params%slatop(pft)

  end function LeafNitrogenContent

  ! ==========================================================================
  
  subroutine EvaluateLeafPhotosynthesis(pft, par_abs, veg_tempk, t_growth, t_home, &
    can_press, can_co2_ppress, can_o2_ppress, veg_esat, can_vpress, gb, nscaler,   &
    rdark_scaler, dayl_factor, btran, vcmax25top, jmax25top, kp25top, lnc_top,     &
    agross, anet, gs, ci)
    !
    ! DESCRIPTION:
    ! Evaluates leaf-level photosynthesis at arbitrary prescribed driver
    ! conditions. This reproduces the full current production call sequence:
    !
    !   GetCanopyGasParameters -> LeafLayerBiophysicalRates ->
    !   LeafLayerMaintenanceRespiration_* -> LeafLayerPhotosynthesis
    !
    ! lb_params' model switches (electron_transport_model, stomatal_model,
    ! stomatal_assim_model, photo_tempsens_model) and hlm_maintresp_leaf_model
    ! are HLM-namelist-controlled in production and so are not set by any call
    ! in this module. The calling driver must set them explicitly
  
    ! ARGUMENTS:
    integer,  intent(in)  :: pft            ! plant functional type index
    real(r8), intent(in)  :: par_abs        ! absorbed PAR per unit leaf area [umol photons/m2 leaf/s]
    real(r8), intent(in)  :: veg_tempk      ! instantaneous leaf temperature [K]
    real(r8), intent(in)  :: t_growth       ! 10-day running-mean growth temperature [K]
    real(r8), intent(in)  :: t_home         ! long-term running-mean home temperature [K]
    real(r8), intent(in)  :: can_press      ! air pressure at the leaf surface [Pa]
    real(r8), intent(in)  :: can_co2_ppress ! CO2 partial pressure at the leaf surface [Pa]
    real(r8), intent(in)  :: can_o2_ppress  ! O2 partial pressure at the leaf surface [Pa]
    real(r8), intent(in)  :: veg_esat       ! saturation vapor pressure at veg_tempk [Pa]
    real(r8), intent(in)  :: can_vpress     ! vapor pressure of the canopy air [Pa]
    real(r8), intent(in)  :: gb             ! leaf boundary layer conductance [umol/m2/s]
    real(r8), intent(in)  :: nscaler        ! leaf nitrogen vertical-scaling factor [0-1]
    real(r8), intent(in)  :: rdark_scaler   ! leaf respiration vertical-scaling factor [0-1], Atkin only
    real(r8), intent(in)  :: dayl_factor    ! day-length photosynthetic-capacity acclimation factor [0-1]
    real(r8), intent(in)  :: btran          ! soil moisture stress factor [0-1]
    real(r8), intent(in)  :: vcmax25top     ! reference (25C, canopy-top) maximum carboxylation rate [umol/m2/s]
    real(r8), intent(in)  :: jmax25top      ! reference (25C, canopy-top) maximum electron transport rate [umol/m2/s]
    real(r8), intent(in)  :: kp25top        ! reference (25C, canopy-top) initial slope of C4 CO2 response [umol/m2/s]
    real(r8), intent(in)  :: lnc_top        ! leaf N content at the canopy top [gN/m2 leaf] (see LeafNitrogenContent)
    real(r8), intent(out) :: agross         ! gross photosynthesis [umolC/m2/s]
    real(r8), intent(out) :: anet           ! net photosynthesis [umolC/m2/s]
    real(r8), intent(out) :: gs             ! leaf stomatal conductance [umol H2O/m2/s]
    real(r8), intent(out) :: ci             ! intracellular CO2 [Pa]

    ! LOCALS:
    real(r8) :: mm_kco2, mm_ko2, co2_cpoint ! Michaelis-Menten constants for CO2/O2, CO2 compensation point [Pa]
    type(leaf_capacity_type) :: cap         ! leaf biophysical capacity/dark respiration at these conditions
    real(r8) :: c13disc                     ! carbon-13 discrimination (unused diagnostic here)
    integer  :: solve_iter                  ! Ci-solver iteration count (unused diagnostic here)

    call GetCanopyGasParameters(can_press, can_o2_ppress, veg_tempk, mm_kco2, mm_ko2, co2_cpoint)

    call LeafLayerCapacity(pft, veg_tempk, t_growth, t_home, nscaler,             &
      rdark_scaler, dayl_factor, btran, vcmax25top, jmax25top, kp25top, lnc_top,  &
      cap)

    call LeafLayerPhotosynthesis(par_abs, pft, cap%vcmax, cap%jmax, cap%kp,       &
      cap%gs0, cap%gs1, cap%gs2, veg_tempk, can_press, can_co2_ppress,            &
      can_o2_ppress, veg_esat, gb, can_vpress, mm_kco2, mm_ko2, co2_cpoint,       &
      cap%lmr, ci_tol, agross, gs, anet, c13disc, ci, solve_iter)
      
  end subroutine EvaluateLeafPhotosynthesis
  
  ! ==========================================================================

  subroutine LeafLayerCapacity(pft, veg_tempk, t_growth, t_home, nscaler,         &
    rdark_scaler, dayl_factor, btran, vcmax25top, jmax25top, kp25top, lnc_top, cap)
    !
    ! DESCRIPTION:
    ! One leaf layer's photosynthetic capacity and dark respiration at the given
    ! leaf conditions
    !
    ! The leaf maintenance respiration model is selected by
    ! hlm_maintresp_leaf_model, which is HLM-namelist-controlled in production
    ! and so must be set explicitly by the calling driver
    !
    ! t_growth doubles as the Atkin et al. (2017) acclimation temperature, the
    ! same quantity production passes there (currentPatch%tveg_lpa%GetMean())
    !

    ! ARGUMENTS:
    integer,  intent(in)  :: pft                 ! plant functional type index
    real(r8), intent(in)  :: veg_tempk           ! instantaneous leaf temperature [K]
    real(r8), intent(in)  :: t_growth            ! 10-day running-mean growth temperature [K]
    real(r8), intent(in)  :: t_home              ! long-term running-mean home temperature [K]
    real(r8), intent(in)  :: nscaler             ! leaf nitrogen vertical-scaling factor [0-1]
    real(r8), intent(in)  :: rdark_scaler        ! leaf respiration vertical-scaling factor [0-1], Atkin only
    real(r8), intent(in)  :: dayl_factor         ! day-length photosynthetic-capacity acclimation factor [0-1]
    real(r8), intent(in)  :: btran               ! soil moisture stress factor [0-1]
    real(r8), intent(in)  :: vcmax25top          ! reference (25C, canopy-top) maximum carboxylation rate [umol/m2/s]
    real(r8), intent(in)  :: jmax25top           ! reference (25C, canopy-top) maximum electron transport rate [umol/m2/s]
    real(r8), intent(in)  :: kp25top             ! reference (25C, canopy-top) initial slope of C4 CO2 response [umol/m2/s]
    real(r8), intent(in)  :: lnc_top             ! leaf N content at the canopy top [gN/m2 leaf] (see LeafNitrogenContent)
    type(leaf_capacity_type), intent(out) :: cap ! this layer's capacity/dark respiration at these conditions

    call LeafLayerBiophysicalRates(pft, vcmax25top, jmax25top, kp25top, nscaler,   &
      veg_tempk, dayl_factor, t_growth, t_home, btran, cap%vcmax, cap%jmax,       &
      cap%kp, cap%gs0, cap%gs1, cap%gs2)

    select case (hlm_maintresp_leaf_model)

    case (lmrmodel_ryan_1991)

      call LeafLayerMaintenanceRespiration_Ryan_1991(lnc_top, nscaler, pft,      &
        veg_tempk, cap%lmr)

    case (lmrmodel_atkin_etal_2017)

      call LeafLayerMaintenanceRespiration_Atkin_etal_2017(lnc_top,              &
        rdark_scaler, pft, veg_tempk, t_growth, cap%lmr)

    case default

      write(*,*) 'LeafLayerCapacity: unrecognized leaf respiration model: ',     &
        hlm_maintresp_leaf_model
      error stop

    end select

  end subroutine LeafLayerCapacity
  
  ! ==========================================================================
  
  subroutine LeafLayerVerticalScaling(treelai, treesai, height, nv, pft,        &
    vcmax25top, lai_above_in, nscaler_z, rdark_scaler_z, vert_scaler_min)
    !
    ! DESCRIPTION:
    ! Per-leaf-layer vertical-scaling factors, the decay of photosynthetic
    ! capacity (nscaler) and of leaf maintenance respiration (rdark_scaler)
    ! with cumulative leaf area above a layer
    !
    ! Both profiles share the same cumulative LAI and differ only in their
    ! decay coefficient: nscaler follows leafn_vert_scaler_coeff1/2 and
    ! rdark_scaler follows maintresp_leaf_vert_scaler_coeff1/2. rdark_scaler
    ! is used by the Atkin et al. (2017) respiration model only, and is
    ! computed regardless of which model is selected
    !
    ! Depends only on canopy structure (treelai/treesai/height/nv) and the
    ! reference canopy-top capacity
    !
    ! The prescribed exponentials have no lower bound, so at large cumulative
    ! LAI they extrapolate below any viable leaf nitrogen content.
    ! vert_scaler_min, when supplied, floors both profiles at a minimum viable
    ! fraction of canopy-top capacity. The two decay coefficients differ, so the
    ! floor binds at different cumulative LAI in each
    !
    ! Snow depth is fixed at zero
    !

    ! ARGUMENTS:
    real(r8), intent(in)  :: treelai           ! in-crown leaf area index [m2 leaf/m2 crown footprint]
    real(r8), intent(in)  :: treesai           ! in-crown stem area index [m2 stem/m2 crown footprint]
    real(r8), intent(in)  :: height            ! plant/canopy height [m]
    real(r8), intent(in)  :: lai_above_in      ! lai above the cohort
    integer,  intent(in)  :: nv                ! number of occupied leaf layers
    integer,  intent(in)  :: pft               ! plant functional type index
    real(r8), intent(in)  :: vcmax25top        ! reference (25C, canopy-top) maximum carboxylation rate [umol/m2/s]
    real(r8), intent(out) :: nscaler_z(:)      ! per-leaf-layer nitrogen-scaling factor [0-1], first nv entries filled
    real(r8), intent(out) :: rdark_scaler_z(:) ! per-leaf-layer respiration-scaling factor [0-1], first nv entries filled
    real(r8), intent(in), optional :: vert_scaler_min ! floor on both scaling factors, unbounded if absent [0-1]

    ! LOCALS:
    integer  :: iv                     ! leaf-layer looping index
    real(r8) :: vai_top, vai_bot       ! vegetation area index bounds of the current leaf layer
    real(r8) :: elai_layer, esai_layer ! exposed leaf/stem area index of the current leaf layer
    real(r8) :: cumulative_lai         ! LAI above the middle of the current leaf layer
    real(r8) :: lai_above              ! running LAI above the current leaf layer
    real(r8) :: kn                     ! nitrogen vertical-scaling decay coefficient
    real(r8) :: kn_rdark               ! respiration vertical-scaling decay coefficient

    real(r8), parameter :: snow_depth = 0.0_r8 ! no snow modeled in any standalone test driver [m]

    ! loop-invariant: depends only on pft and the reference canopy-top
    ! capacity, neither of which varies by layer
    kn = DecayCoeffVcmax(vcmax25top, prt_params%leafn_vert_scaler_coeff1(pft), &
      prt_params%leafn_vert_scaler_coeff2(pft))
    kn_rdark = DecayCoeffVcmax(vcmax25top,                                    &
      EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff1(pft),                &
      EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff2(pft))

    lai_above = lai_above_in
    do iv = 1, nv
      call VegAreaLayer(treelai, treesai, height, iv, nv, pft, snow_depth,     &
        vai_top, vai_bot, elai_layer, esai_layer)
      cumulative_lai = lai_above + 0.5_r8*elai_layer
      lai_above = lai_above + elai_layer
      nscaler_z(iv) = exp(-kn*cumulative_lai)
      rdark_scaler_z(iv) = exp(-kn_rdark*cumulative_lai)
      if (present(vert_scaler_min)) then
        nscaler_z(iv) = max(nscaler_z(iv), vert_scaler_min)
        rdark_scaler_z(iv) = max(rdark_scaler_z(iv), vert_scaler_min)
      end if
    end do

  end subroutine LeafLayerVerticalScaling
  
  ! ==========================================================================

  function SunlitFraction(laisun_z, laisha_z) result(fsun)
    !
    ! DESCRIPTION:
    ! Sunlit fraction of one leaf layer's leaf area, for area-weighting that
    ! layer's sunlit and shaded photosynthesis into a single per-layer rate.
    ! Returns zero for a layer holding effectively no leaf area at all

    ! ARGUMENTS:
    real(r8), intent(in) :: laisun_z ! this layer's sunlit leaf area index [m2 leaf/m2 crown footprint]
    real(r8), intent(in) :: laisha_z ! this layer's shaded leaf area index [m2 leaf/m2 crown footprint]
    real(r8)             :: fsun     ! sunlit fraction of this layer's leaf area [0-1]

    if (laisun_z + laisha_z > nearzero) then
      fsun = laisun_z / (laisun_z + laisha_z)
    else
      fsun = 0.0_r8
    end if

  end function SunlitFraction

  ! ==========================================================================

  subroutine LeafLayerSunShade(pft, cap, par_sun_z, par_sha_z, lai_sun_z,        &
    lai_sha_z, veg_tempk, can_press, can_co2_ppress, can_o2_ppress, veg_esat,    &
    can_vpress, gb, mm_kco2, mm_ko2, co2_cpoint, agross_layer, anet_layer,       &
    lai_layer, solve_iter_sun, solve_iter_sha)
    !
    ! DESCRIPTION:
    ! One leaf layer's sunlit and shaded photosynthesis at a given capacity,
    ! area-weighted by that layer's sunlit fraction into a single per-layer rate
    ! per unit leaf area

    ! ARGUMENTS:
    integer,  intent(in)                 :: pft            ! plant functional type index
    type(leaf_capacity_type), intent(in) :: cap            ! this layer's capacity/dark respiration (see LeafLayerCapacity)
    real(r8), intent(in)                 :: par_sun_z      ! absorbed PAR, sunlit leaves, this layer [W/m2 crown footprint]
    real(r8), intent(in)                 :: par_sha_z      ! absorbed PAR, shaded leaves, this layer [W/m2 crown footprint]
    real(r8), intent(in)                 :: lai_sun_z      ! sunlit leaf area index, this layer [m2 leaf/m2 crown footprint]
    real(r8), intent(in)                 :: lai_sha_z      ! shaded leaf area index, this layer [m2 leaf/m2 crown footprint]
    real(r8), intent(in)                 :: veg_tempk      ! leaf temperature [K]
    real(r8), intent(in)                 :: can_press      ! air pressure at the leaf surface [Pa]
    real(r8), intent(in)                 :: can_co2_ppress ! CO2 partial pressure at the leaf surface [Pa]
    real(r8), intent(in)                 :: can_o2_ppress  ! O2 partial pressure at the leaf surface [Pa]
    real(r8), intent(in)                 :: veg_esat       ! saturation vapor pressure at veg_tempk [Pa]
    real(r8), intent(in)                 :: can_vpress     ! vapor pressure of the canopy air [Pa]
    real(r8), intent(in)                 :: gb             ! leaf boundary layer conductance [umol/m2/s]
    real(r8), intent(in)                 :: mm_kco2        ! Michaelis-Menten constant for CO2 at veg_tempk [Pa]
    real(r8), intent(in)                 :: mm_ko2         ! Michaelis-Menten constant for O2 at veg_tempk [Pa]
    real(r8), intent(in)                 :: co2_cpoint     ! CO2 compensation point at veg_tempk [Pa]
    real(r8), intent(out)                :: agross_layer   ! area-weighted gross photosynthesis for this layer [umolC/m2 leaf/s]
    real(r8), intent(out)                :: anet_layer     ! area-weighted net photosynthesis for this layer [umolC/m2 leaf/s]
    real(r8), intent(out)                :: lai_layer      ! this layer's total leaf area index [m2 leaf/m2 crown footprint]
    integer,  intent(out), optional      :: solve_iter_sun ! Ci-solver iteration count, sunlit call, for callers tracking solver diagnostics
    integer,  intent(out), optional      :: solve_iter_sha ! Ci-solver iteration count, shaded call, for callers tracking solver diagnostics

    ! LOCALS:
    real(r8) :: par_abs                  ! absorbed PAR per unit leaf area [umol photons/m2 leaf/s]
    real(r8) :: agross_sun, agross_sha   ! gross photosynthesis, sunlit/shaded [umolC/m2 leaf/s]
    real(r8) :: anet_sun, anet_sha       ! net photosynthesis, sunlit/shaded [umolC/m2 leaf/s]
    real(r8) :: gs_sun, gs_sha           ! stomatal conductance, sunlit/shaded (unused diagnostic here)
    real(r8) :: ci_sun, ci_sha           ! intracellular CO2, sunlit/shaded (unused diagnostic here)
    real(r8) :: c13disc_sun, c13disc_sha ! carbon-13 discrimination, sunlit/shaded (unused diagnostic here)
    real(r8) :: fsun                     ! sunlit fraction of this layer's leaf area [0-1]
    integer  :: solve_iter_sun_l         ! Ci-solver iteration count, sunlit call
    integer  :: solve_iter_sha_l         ! Ci-solver iteration count, shaded call

    par_abs = ConvertPar(lai_sun_z, par_sun_z)
    call LeafLayerPhotosynthesis(par_abs, pft, cap%vcmax, cap%jmax, cap%kp,       &
      cap%gs0, cap%gs1, cap%gs2, veg_tempk, can_press, can_co2_ppress,           &
      can_o2_ppress, veg_esat, gb, can_vpress, mm_kco2, mm_ko2, co2_cpoint,      &
      cap%lmr, ci_tol, agross_sun, gs_sun, anet_sun, c13disc_sun, ci_sun,        &
      solve_iter_sun_l)

    par_abs = ConvertPar(lai_sha_z, par_sha_z)
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

end module FatesTestLeafPhotoMod
