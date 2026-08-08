module FatesTestLeafPhotoMod
  !
  ! DESCRIPTION:
  ! Stateless leaf-level photosynthesis evaluator for standalone test drivers -
  ! wraps the current production call sequence (GetCanopyGasParameters ->
  ! LeafLayerBiophysicalRates -> LeafLayerMaintenanceRespiration_Ryan_1991 ->
  ! LeafLayerPhotosynthesis) so a test can evaluate leaf gas exchange at any
  ! prescribed, arbitrary driver conditions - no cohort, no canopy layers, no
  ! day/year loop, unlike FatesTestCohortPhysMod's cohort_phys_type (which
  ! caches per-layer capacity across a simulated day/year and is tightly
  ! coupled to a cohort's dbh/height/crown geometry). Every argument to
  ! EvaluateLeafPhotosynthesis below is recomputed from scratch on every call,
  ! since a leaf-level sensitivity sweep needs to vary things (temperature,
  ! nscaler, btran) that cohort_phys_type's GrossAssimAndResp/LeafNetAssimSweep
  ! hold fixed at whatever a cohort's daily setup already left them.
  !
  ! Replaces this repo's original leaf-level photosynthesis test support
  ! module (FatesTestPhotosynthesisMod.F90, not present in this repo - it
  ! predates the reorganization described below and no longer compiles). Its
  ! one piece of real logic, a bespoke saturation-vapor-pressure polynomial
  ! (CalcVaporPressure/sat_vapor_press, Bonan 2019 textbook constants), is
  ! superseded here by LeafBiophysicsMod's own QSat (an exact clone of CTSM's
  ! QSatMod.F90) - reused directly rather than reimplemented, and already the
  ! convention FatesTestEnvironmentMod uses for the same purpose.
  !
  ! What else has changed since the original test was last updated (see
  ! EvaluateLeafPhotosynthesis's header comment for the call-by-call detail):
  ! GetCanopyGasParameters/LeafLayerBiophysicalRates/LeafLayerPhotosynthesis
  ! moved from FatesPlantRespPhotosynthMod to LeafBiophysicsMod; stomatal
  ! conductance is now a proper (gs0,gs1,gs2) triple computed inside
  ! LeafLayerBiophysicalRates rather than a hand-rolled btran-scaled
  ! intercept; leaf dark respiration is now the real, PFT-parameterized
  ! LeafLayerMaintenanceRespiration_Ryan_1991 rather than a hardcoded
  ! 0.015/0.025*vcmax25top coefficient; leaf boundary-layer conductance (gb)
  ! is a direct input rather than derived from a resistance inside
  ! GetCanopyGasParameters; and the default temperature-acclimation model
  ! (Kumarathunge et al. 2019) needs two independent reference temperatures
  ! (t_growth, t_home) rather than the one canopy temperature the old
  ! interface reused for both.
  !

  use FatesConstantsMod, only : r8 => fates_r8
  use LeafBiophysicsMod, only : GetCanopyGasParameters
  use LeafBiophysicsMod, only : LeafLayerBiophysicalRates
  use LeafBiophysicsMod, only : LeafLayerMaintenanceRespiration_Ryan_1991
  use LeafBiophysicsMod, only : LeafLayerPhotosynthesis
  use EDPftvarcon,       only : EDPftvarcon_inst
  use PRTParametersMod,  only : prt_params
  use PRTGenericMod,     only : leaf_organ
  use FatesParameterDerivedMod, only : param_derived

  implicit none
  private

  ! convergence tolerance for intracellular CO2 [Pa] - matches
  ! FatesTestCohortPhysMod.F90's identical parameter and production's own use
  real(r8), parameter :: ci_tol = 0.5_r8

  public :: EvaluateLeafPhotosynthesis
  public :: LeafNitrogenContent

contains

  ! ==========================================================================

  function LeafNitrogenContent(pft) result(lnc_top)
    !
    ! DESCRIPTION:
    ! Leaf nitrogen content at the canopy top [gN/m2 leaf] - the same
    ! calculation test_SingleCohort.F90 uses for its own lnc_top, needed by
    ! LeafLayerMaintenanceRespiration_Ryan_1991. Requires PRTDerivedParams()
    ! to have already populated prt_params%organ_param_id (see
    ! test_LeafLevelPhoto.F90's setup sequence)

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
    dayl_factor, btran, lnc_top, agross, anet, gs, ci)
    !
    ! DESCRIPTION:
    ! Evaluates leaf-level photosynthesis at arbitrary prescribed driver
    ! conditions - the full current production call sequence
    ! (GetCanopyGasParameters -> LeafLayerBiophysicalRates ->
    ! LeafLayerMaintenanceRespiration_Ryan_1991 -> LeafLayerPhotosynthesis),
    ! with no self-shading and no canopy structure: every argument is a
    ! single instantaneous scalar, and every biophysical rate is recomputed
    ! from these arguments alone (nothing is cached between calls).
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
    real(r8), intent(in)  :: lnc_top         ! leaf N content at the canopy top [gN/m2 leaf] (see LeafNitrogenContent)
    real(r8), intent(out) :: agross          ! gross photosynthesis [umolC/m2/s]
    real(r8), intent(out) :: anet            ! net photosynthesis [umolC/m2/s]
    real(r8), intent(out) :: gs              ! leaf stomatal conductance [umol H2O/m2/s]
    real(r8), intent(out) :: ci              ! intracellular CO2 [Pa]

    ! LOCALS:
    real(r8) :: mm_kco2, mm_ko2, co2_cpoint ! Michaelis-Menten constants for CO2/O2, CO2 compensation point [Pa]
    real(r8) :: vcmax, jmax, kp             ! leaf biophysical capacity at these conditions
    real(r8) :: gs0, gs1, gs2               ! stomatal conductance slope/intercept terms
    real(r8) :: lmr                         ! leaf maintenance (dark) respiration [umolCO2/m2/s]
    real(r8) :: c13disc                     ! carbon-13 discrimination (unused diagnostic here)
    integer  :: solve_iter                  ! Ci-solver iteration count (unused diagnostic here)

    call GetCanopyGasParameters(can_press, can_o2_ppress, veg_tempk, mm_kco2, mm_ko2, co2_cpoint)

    call LeafLayerBiophysicalRates(pft, EDPftvarcon_inst%vcmax25top(pft,1),        &
      param_derived%jmax25top(pft,1), param_derived%kp25top(pft,1), nscaler,      &
      veg_tempk, dayl_factor, t_growth, t_home, btran, vcmax, jmax, kp, gs0, gs1, gs2)

    call LeafLayerMaintenanceRespiration_Ryan_1991(lnc_top, nscaler, pft, veg_tempk, lmr)

    call LeafLayerPhotosynthesis(par_abs, pft, vcmax, jmax, kp, gs0, gs1, gs2,     &
      veg_tempk, can_press, can_co2_ppress, can_o2_ppress, veg_esat, gb,          &
      can_vpress, mm_kco2, mm_ko2, co2_cpoint, lmr, ci_tol, agross, gs, anet,     &
      c13disc, ci, solve_iter)

  end subroutine EvaluateLeafPhotosynthesis

end module FatesTestLeafPhotoMod
