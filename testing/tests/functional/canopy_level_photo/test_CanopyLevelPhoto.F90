program FatesCanopyLevelPhoto
  !
  ! DESCRIPTION:
  ! Canopy-level photosynthesis sensitivity sweeps, replicating the canopy
  ! panels of Rogers et al. (2017, New Phytologist 213:22-42, "A roadmap for
  ! improving the representation of photosynthesis in Earth system models") -
  ! the canopy-scale companion to test_LeafLevelPhoto.F90, which already
  ! replicates that paper's leaf-level panels.
  !
  ! The canopy here is PRESCRIBED, not grown: leaf area index is an
  ! experimental treatment (Rogers' Fig. 1 spans LAI = 1, 3 and 7), not an
  ! allometric consequence of a plant's size. There is therefore no cohort,
  ! no allometry and no patch/site anywhere in this driver - canopy structure
  ! enters as plain treelai/treesai/height scalars, handed straight to
  ! FatesTestLightEnvMod's light_env_type (whose Init/Refresh take exactly
  ! those scalars for this reason) and, once photosynthesis is wired in, to
  ! FatesTestLeafPhotoMod's LeafLayerNitrogenScaling. Both are shared with
  ! test_SingleCohort.F90, which drives them from a real cohort's allometry
  ! instead - identical two-stream attenuation and nitrogen-decay physics,
  ! differing only in where the canopy structure came from.
  !
  ! *** THIS IS AN INTERMEDIATE STAGE. *** It resolves the within-canopy light
  ! profile and integrates canopy net photosynthesis at Rogers' STANDARD
  ! CONDITIONS only (25 C, 1 kPa VPD, 380 umol/mol CO2, Q = 1500 umol/m2/s).
  ! None of Rogers' five sweeps are run yet - that is the next step, and the
  ! per-layer machinery exercised here is exactly what they will drive.
  !
  ! CANOPY INTEGRATION: per leaf layer, absorbed PAR is converted from a per-
  ! ground-area flux to the per-leaf-area PPFD leaf photosynthesis expects
  ! (LayerParPerLeafArea), and EvaluateLeafPhotosynthesis is called twice -
  ! once sunlit, once shaded - at that layer's nitrogen-scaled capacity
  ! (LeafLayerNitrogenScaling). The two are combined by sunlit fraction
  ! (SunlitFraction) into one per-layer rate [umolC/m2 leaf/s], multiplied by
  ! that layer's leaf area index [m2 leaf/m2 ground] and summed down the
  ! canopy, giving canopy net photosynthesis per unit GROUND area
  ! [umolC/m2 ground/s] - the y-axis of Rogers' canopy panels, as distinct
  ! from the per-leaf-area y-axis of the leaf panels test_LeafLevelPhoto.F90
  ! reproduces.
  !
  ! Rogers' canopy A is NET leaf-level assimilation integrated over the
  ! canopy: leaf dark respiration is included (it is internal to anet), while
  ! non-leaf (stem/coarse root/fine root) maintenance respiration is NOT -
  ! this is a canopy flux, not a whole-plant carbon balance. That is why
  ! FatesTestCohortPhysMod's GrossAssimAndResp is deliberately not reused
  ! here: it returns gross assimilation and whole-plant respiration scaled
  ! per individual, at frozen daily capacity, none of which matches this
  ! protocol (which must recompute capacity per swept point).
  !
  ! PFT: target_pft below is 6, broadleaf_colddecid_extratrop_tree, matching
  ! Rogers' stated "generic temperate broad leaved deciduous tree".
  ! test_LeafLevelPhoto.F90 uses the same PFT so the leaf and canopy panels
  ! describe the same plant and can be shown as a pair. The choice matters
  ! far more at canopy scale than at leaf scale: PFT 1
  ! (broadleaf_evergreen_tropical_tree) carries much steeper Lloyd et al.
  ! nitrogen-decay coefficients (0.0079/1.357, kn = 0.41) than PFT 6
  ! (0.00963/2.43, kn = 0.15), which suppresses deep-canopy capacity by ~5x
  ! at LAI = 7 while leaving leaf-level results untouched (leaf level never
  ! applies nscaler at all).
  !
  ! Vcmax25top is this PFT's own native value (52 umol/m2/s), NOT Rogers'
  ! standardized 60. This driver therefore answers "what does FATES do"
  ! rather than "what does FATES do under Rogers' standardized protocol",
  ! and its curves are expected to sit somewhat BELOW the published model
  ! spread for that reason alone - a deliberate choice, not a discrepancy to
  ! chase.
  !
  ! ILLUMINATION GEOMETRY: Rogers' protocol specifies the incident quantum
  ! flux density (Q) but states neither a solar zenith angle nor a direct/
  ! diffuse split, so both are prescribed here as documented assumptions:
  !   - coszen is fixed at 1 (sun directly overhead) for every sweep. This
  !     matches test_SingleCohort.F90's LightResponseSweep, and deliberately
  !     decouples the response-curve shape from solar geometry so curves
  !     differ only via canopy structure - the point of a reference sweep.
  !   - the direct/diffuse split reuses FatesTestLightEnvMod's own clear-sky
  !     direct_frac (85% direct/15% diffuse), so this driver's light
  !     environment is consistent with what test_SingleCohort.F90 already
  !     assumes rather than introducing a second, conflicting convention,
  !     and because a mostly-direct beam is the physically realistic case
  !     for the full-sun end of Rogers' Q sweep.
  ! Neither is a value taken from the paper; both are choices this driver
  ! makes on the paper's behalf. The split is a live sensitivity, not a
  ! detail: canopy Anet at LAI = 7 spans roughly 21 (pure beam) to 39 (pure
  ! diffuse) umolC/m2 ground/s, so it should be revisited before any
  ! quantitative claim about agreement with a specific published model.
  !
  ! CANOPY HEIGHT is likewise not part of Rogers' protocol, but VegAreaLayer
  ! (reached via light_env_type) requires one. With no snow modeled it enters
  ! only through the crown-depth calculation, so it is close to inert here;
  ! canopy_height below is a plausible value for the paper's "generic
  ! temperate broad leaved deciduous tree", not a fitted or paper-sourced
  ! number.
  !
  ! OUTPUT:
  !   canopy_anet - canopy net photosynthesis [umolC/m2 ground/s], per LAI
  !   the per-leaf-layer profile at each prescribed LAI, dimensioned
  !   (nlevleaf, lai): parsun_z/parsha_z (absorbed PAR, sunlit/shaded
  !   [W/m2 ground]), laisun_z/laisha_z (sunlit/shaded leaf area index
  !   [m2/m2]), nscaler_z (nitrogen-scaling factor [0-1]) and anet_z (that
  !   layer's area-weighted net photosynthesis [umolC/m2 leaf/s])
  ! nlevleaf (EDParamsMod) is a compile-time maximum layer count; layers
  ! above a given LAI's actual nv are left at fates_unset_r8 (registered as
  ! each variable's _FillValue), so the array survives nv differing across
  ! the swept LAI values.
  !
  ! COMMAND LINE: ./CanopyLevelPhoto_exe <parameter_file.json>
  !

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesConstantsMod,           only : fates_unset_r8
  use FatesConstantsMod,           only : wm2_to_umolm2s
  use FatesConstantsMod,           only : t_water_freeze_k_1atm
  use EDParamsMod,                 only : nlevleaf
  use EDParamsMod,                 only : GetNVegLayers
  use EDPftvarcon,                 only : EDPftvarcon_inst
  use PRTParametersMod,            only : prt_params
  use PRTInitParamsFatesMod,       only : PRTDerivedParams
  use FatesParameterDerivedMod,    only : param_derived
  use FatesUnitTestParamReaderMod, only : ReadParameters
  use FatesArgumentUtils,          only : command_line_arg
  use FatesFactoryMod,             only : InitializeGlobals
  use FatesGlobals,                only : FatesGlobalsInit
  use FatesInterfaceTypesMod,      only : numpft
  use FatesTwoStreamUtilsMod,      only : TransferRadParams
  use LeafBiophysicsMod,           only : lb_params
  use LeafBiophysicsMod,           only : FvCB1980, medlyn_model, net_assim_model
  use LeafBiophysicsMod,           only : photosynth_acclim_model_kumarathunge_etal_2019
  use LeafBiophysicsMod,           only : QSat
  use FatesTestEnvironmentMod,     only : environment_type
  use FatesTestLightEnvMod,        only : light_env_type
  use FatesTestLightEnvMod,        only : direct_frac, diffuse_frac
  use FatesTestLeafPhotoMod,       only : EvaluateLeafPhotosynthesis
  use FatesTestLeafPhotoMod,       only : LeafNitrogenContent
  use FatesTestLeafPhotoMod,       only : LeafLayerNitrogenScaling
  use FatesTestLeafPhotoMod,       only : LayerParPerLeafArea
  use FatesTestLeafPhotoMod,       only : SunlitFraction
  use FatesUnitTestIOMod,          only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod,          only : WriteVar, RegisterVar, EndNCDef
  use FatesUnitTestIOMod,          only : RegisterFillValue
  use FatesUnitTestIOMod,          only : type_double, type_int

  implicit none

  ! LOCALS:
  character(len=:), allocatable :: param_file ! input parameter file
  type(environment_type)        :: env        ! prescribed atmospheric boundary conditions
  type(light_env_type)          :: light_env  ! prescribed light environment for the current LAI
  real(r8) :: lnc_top        ! leaf N content at the canopy top [gN/m2 leaf]
  real(r8) :: vcmax25top_pft ! reference (25C, canopy-top) maximum carboxylation rate [umol/m2/s]
  real(r8) :: jmax25top_pft  ! reference (25C, canopy-top) maximum electron transport rate [umol/m2/s]
  real(r8) :: kp25top_pft    ! reference (25C, canopy-top) initial slope of C4 CO2 response [umol/m2/s]

  ! per-leaf-layer light profile at each prescribed LAI, (nlevleaf, n_lai)
  real(r8), allocatable :: parsun_z_out(:,:) ! absorbed PAR, sunlit [W/m2 ground]
  real(r8), allocatable :: parsha_z_out(:,:) ! absorbed PAR, shaded [W/m2 ground]
  real(r8), allocatable :: laisun_z_out(:,:) ! sunlit leaf area index [m2/m2]
  real(r8), allocatable :: laisha_z_out(:,:) ! shaded leaf area index [m2/m2]
  real(r8), allocatable :: nscaler_z_out(:,:) ! per-layer nitrogen-scaling factor [0-1]
  real(r8), allocatable :: anet_z_out(:,:)    ! per-layer area-weighted net photosynthesis [umolC/m2 leaf/s]
  integer,  allocatable :: nv_out(:)         ! number of occupied leaf layers at each prescribed LAI
  integer,  allocatable :: layer_index(:)    ! leaf-layer index coordinate [1..nlevleaf]

  ! canopy-integrated output, per prescribed LAI
  real(r8), allocatable :: canopy_anet(:) ! canopy net photosynthesis [umolC/m2 ground/s]

  ! per-layer working arrays, sized to the current LAI's nv
  real(r8), allocatable :: nscaler_z(:) ! per-leaf-layer nitrogen-scaling factor [0-1]

  ! the shared default vapor-pressure state at the standard reference
  ! condition - leaf temperature and VPD both fixed, so these are constant
  ! across the whole run at this stage
  real(r8) :: default_veg_esat   ! saturation vapor pressure at default_veg_tempk [Pa]
  real(r8) :: default_can_vpress ! canopy air vapor pressure at the default VPD [Pa]
  real(r8) :: qs_dummy           ! saturation specific humidity output from QSat (unused here)

  ! per-layer photosynthesis results
  real(r8) :: par_abs                  ! absorbed PAR per unit leaf area [umol photons/m2 leaf/s]
  real(r8) :: agross_sun, agross_sha   ! gross photosynthesis, sunlit/shaded [umolC/m2 leaf/s] (unused diagnostics)
  real(r8) :: anet_sun, anet_sha       ! net photosynthesis, sunlit/shaded [umolC/m2 leaf/s]
  real(r8) :: gs_sun, gs_sha           ! stomatal conductance, sunlit/shaded (unused diagnostics)
  real(r8) :: ci_sun, ci_sha           ! intracellular CO2, sunlit/shaded (unused diagnostics)
  real(r8) :: fsun                     ! sunlit fraction of this layer's leaf area [0-1]
  real(r8) :: anet_layer               ! area-weighted net photosynthesis for this layer [umolC/m2 leaf/s]

  integer :: ilai ! prescribed-LAI looping index
  integer :: iv   ! leaf-layer looping index

  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'canopy_level_photo_out.nc' ! output file
  integer,  parameter :: target_pft = 6 ! PFT index to evaluate (1-based)

  ! prescribed canopy - Rogers' Fig. 1 canopy panels (b, c, d)
  integer,  parameter :: n_lai = 3
  real(r8), parameter :: lai_vals(n_lai) = [1.0_r8, 3.0_r8, 7.0_r8] ! prescribed leaf area index [m2 leaf/m2 ground]
  real(r8), parameter :: canopy_sai    = 0.0_r8  ! prescribed stem area index [m2 stem/m2 ground] - Rogers' protocol is a leaf-area experiment with no stem-area treatment, so stems are excluded rather than invented
  real(r8), parameter :: canopy_height = 20.0_r8 ! prescribed canopy height [m] - see the header comment; near-inert with no snow

  ! illumination geometry - see the header comment for why this is prescribed
  ! here rather than taken from Rogers
  real(r8), parameter :: diagnostic_coszen = 1.0_r8 ! cosine of solar zenith angle (sun directly overhead)

  ! incident PAR this intermediate stage resolves the profile at - the
  ! standard-condition saturating Q of Rogers' non-light sweeps
  ! (1500 umol/m2/s), converted to the W/m2 light_env works in
  real(r8), parameter :: reference_par = 1500.0_r8/wm2_to_umolm2s ! [W/m2]

  ! Rogers' standard reference conditions - identical to
  ! test_LeafLevelPhoto.F90's, so the canopy panels sit on exactly the same
  ! reference state as the leaf panels. CO2/O2/dayl_factor/btran come from
  ! FatesTestEnvironmentMod's shared reference-atmosphere defaults (env%Init)
  ! rather than being redeclared here
  real(r8), parameter :: default_veg_tempk = 25.0_r8 + t_water_freeze_k_1atm ! [K]
  real(r8), parameter :: default_vpd       = 1000.0_r8 ! [Pa] leaf-to-air VPD, esat(Tleaf) - eair


  ! read in parameter file name from command line
  param_file = command_line_arg(1)

  ! read in parameter file
  call ReadParameters(param_file)

  ! step_size is unused by anything this test exercises (no time-stepping
  ! occurs), but InitializeGlobals requires a value
  call InitializeGlobals(86400.0_r8)

  ! initialize FATES logging (the two-stream module's debug/error paths write
  ! to this unit and it is otherwise never set in a standalone driver) and the
  ! two-stream radiation parameters (leaf/stem optical properties, per pft)
  numpft = size(prt_params%wood_density, dim=1)
  call FatesGlobalsInit(6, .false.)
  call TransferRadParams()

  ! derive the organ_id -> parameter-file-index reverse lookup map
  ! (prt_params%organ_param_id) - normally done by the host model's own
  ! interface setup (FatesInterfaceMod.F90), which this standalone driver
  ! bypasses entirely
  call PRTDerivedParams()

  ! host-model-namelist-controlled leaf biophysics switches - matches
  ! test_LeafLevelPhoto.F90/test_SingleCohort.F90's choices for consistency
  ! across this test suite. Not consumed at this light-profile-only stage, but
  ! set here so the setup sequence is complete and correct as photosynthesis
  ! is added on top
  lb_params%electron_transport_model = FvCB1980
  lb_params%stomatal_model           = medlyn_model
  lb_params%stomatal_assim_model     = net_assim_model
  lb_params%photo_tempsens_model     = photosynth_acclim_model_kumarathunge_etal_2019

  ! leaf N content and reference (25C, canopy-top) photosynthetic capacity for
  ! target_pft - the flat PFT defaults (no cohort/acclimation state exists in
  ! this driver). Constant for the whole run; consumed once photosynthesis is
  ! wired in
  lnc_top        = LeafNitrogenContent(target_pft)
  vcmax25top_pft = EDPftvarcon_inst%vcmax25top(target_pft,1)
  jmax25top_pft  = param_derived%jmax25top(target_pft,1)
  kp25top_pft    = param_derived%kp25top(target_pft,1)

  ! prescribed atmospheric boundary conditions (canopy pressure, background
  ! CO2/O2 partial pressure, leaf boundary-layer conductance, dayl_factor,
  ! btran - see FatesTestEnvironmentMod's shared reference-atmosphere defaults)
  call env%Init()

  ! ---------------------------------------------------------------------------------------
  ! resolve the within-canopy light profile at each prescribed LAI
  ! ---------------------------------------------------------------------------------------

  allocate(parsun_z_out(nlevleaf, n_lai), parsha_z_out(nlevleaf, n_lai))
  allocate(laisun_z_out(nlevleaf, n_lai), laisha_z_out(nlevleaf, n_lai))
  allocate(nscaler_z_out(nlevleaf, n_lai), anet_z_out(nlevleaf, n_lai))
  allocate(nv_out(n_lai), layer_index(nlevleaf), canopy_anet(n_lai))

  ! layers above a given LAI's actual nv are never written, and keep this fill
  ! value (registered as each variable's _FillValue below)
  parsun_z_out(:,:)  = fates_unset_r8
  parsha_z_out(:,:)  = fates_unset_r8
  laisun_z_out(:,:)  = fates_unset_r8
  laisha_z_out(:,:)  = fates_unset_r8
  nscaler_z_out(:,:) = fates_unset_r8
  anet_z_out(:,:)    = fates_unset_r8

  do iv = 1, nlevleaf
    layer_index(iv) = iv
  end do

  ! the standard reference condition's vapor-pressure state - leaf
  ! temperature and VPD are both fixed at this stage, so this is computed
  ! once for the whole run
  call QSat(default_veg_tempk, env%can_press, qs_dummy, default_veg_esat)
  default_can_vpress = default_veg_esat - default_vpd

  do ilai = 1, n_lai

    ! build the prescribed canopy's two-stream light environment - no cohort:
    ! Init takes canopy structure as scalars precisely so a prescribed-LAI
    ! driver like this one can use it
    call light_env%Init(lai_vals(ilai), canopy_sai, canopy_height, target_pft)

    nv_out(ilai) = GetNVegLayers(lai_vals(ilai) + canopy_sai)

    ! this canopy's per-layer nitrogen scaling - the vertical decay of
    ! photosynthetic capacity with cumulative leaf area above each layer,
    ! from the same shared routine test_SingleCohort.F90 drives off a real
    ! cohort's allometry
    if (allocated(nscaler_z)) deallocate(nscaler_z)
    allocate(nscaler_z(nv_out(ilai)))
    call LeafLayerNitrogenScaling(lai_vals(ilai), canopy_sai, canopy_height,   &
      nv_out(ilai), target_pft, vcmax25top_pft, nscaler_z)

    ! the reference case: this driver's prescribed clear-sky beam fraction
    call CanopyNetAssim(nv_out(ilai), nscaler_z, reference_par, direct_frac,   &
      default_veg_tempk, default_veg_esat, default_can_vpress,                 &
      env%can_co2_ppress, env%btran, vcmax25top_pft, canopy_anet(ilai),        &
      store_profile=.true., ilai_store=ilai)

    write(*, '(a, f6.2, a, i3, a, f8.3, a)') ' LAI = ', lai_vals(ilai),        &
      ' (', nv_out(ilai), ' leaf layers): canopy Anet = ',                     &
      canopy_anet(ilai), ' umolC/m2 ground/s'

    call light_env%Free()

  end do

  call WriteOutput()

contains

  ! ==========================================================================

  subroutine CanopyNetAssim(nv, nscaler_z, par_toc, beam_frac, veg_tempk,      &
    veg_esat, can_vpress, can_co2_ppress, btran, vcmax25top, canopy_anet_out,  &
    store_profile, ilai_store)
    !
    ! DESCRIPTION:
    ! Integrates leaf photosynthesis down an already-Init'd canopy
    ! (light_env, holding this LAI's structure) to a canopy net assimilation
    ! per unit ground area - see the program header for the method.
    !
    ! Every argument that Rogers varies in one of his five sweeps is an
    ! argument here (incident PAR, leaf temperature, canopy air vapor
    ! pressure, CO2 partial pressure, btran, vcmax25top), so each sweep is a
    ! loop over calls to this routine with one argument moving.
    !
    ! Writes nothing except canopy_anet_out unless store_profile is set, in
    ! which case the per-layer diagnostics for column ilai_store of the
    ! output arrays are filled as well.

    ! ARGUMENTS:
    integer,  intent(in)  :: nv              ! number of occupied leaf layers
    real(r8), intent(in)  :: nscaler_z(:)    ! per-leaf-layer nitrogen-scaling factor [0-1]
    real(r8), intent(in)  :: par_toc         ! incident PAR at the top of the canopy [W/m2]
    real(r8), intent(in)  :: beam_frac       ! fraction of par_toc arriving as direct beam [0-1]
    real(r8), intent(in)  :: veg_tempk       ! leaf temperature [K]
    real(r8), intent(in)  :: veg_esat        ! saturation vapor pressure at veg_tempk [Pa]
    real(r8), intent(in)  :: can_vpress      ! canopy air vapor pressure [Pa]
    real(r8), intent(in)  :: can_co2_ppress  ! CO2 partial pressure at the leaf surface [Pa]
    real(r8), intent(in)  :: btran           ! soil moisture stress factor [0-1]
    real(r8), intent(in)  :: vcmax25top      ! reference (25C, canopy-top) maximum carboxylation rate [umol/m2/s]
    real(r8), intent(out) :: canopy_anet_out ! canopy net photosynthesis [umolC/m2 ground/s]
    logical,  intent(in), optional :: store_profile ! also fill the per-layer output arrays (default .false.)
    integer,  intent(in), optional :: ilai_store    ! output-array column to fill when store_profile is set

    ! LOCALS:
    logical :: do_store ! resolved store_profile

    do_store = .false.
    if (present(store_profile)) do_store = store_profile

    ! attenuate the incident PAR through this canopy at the prescribed
    ! overhead-sun geometry. AttenuateCanopy (rather than Profile) is called
    ! directly because Profile derives its incident PAR from a light fraction
    ! and a day/hour solar cycle, which this static, prescribed-Q protocol
    ! has no equivalent of
    call light_env%AttenuateCanopy(beam_frac*par_toc,                          &
      (1.0_r8 - beam_frac)*par_toc, diagnostic_coszen,                         &
      light_env%parsun_z, light_env%parsha_z, light_env%laisun_z,              &
      light_env%laisha_z)

    canopy_anet_out = 0.0_r8
    do iv = 1, nv

      ! sunlit leaves in this layer
      par_abs = LayerParPerLeafArea(light_env%parsun_z(iv), light_env%laisun_z(iv))
      call EvaluateLeafPhotosynthesis(target_pft, par_abs, veg_tempk,          &
        veg_tempk, veg_tempk, env%can_press, can_co2_ppress,                   &
        env%can_o2_ppress, veg_esat, can_vpress, env%gb, nscaler_z(iv),        &
        env%dayl_factor, btran, vcmax25top, jmax25top_pft, kp25top_pft,        &
        lnc_top, agross_sun, anet_sun, gs_sun, ci_sun)

      ! shaded leaves in this layer
      par_abs = LayerParPerLeafArea(light_env%parsha_z(iv), light_env%laisha_z(iv))
      call EvaluateLeafPhotosynthesis(target_pft, par_abs, veg_tempk,          &
        veg_tempk, veg_tempk, env%can_press, can_co2_ppress,                   &
        env%can_o2_ppress, veg_esat, can_vpress, env%gb, nscaler_z(iv),        &
        env%dayl_factor, btran, vcmax25top, jmax25top_pft, kp25top_pft,        &
        lnc_top, agross_sha, anet_sha, gs_sha, ci_sha)

      ! area-weighted mean net photosynthesis for the layer, then scaled by
      ! this layer's leaf area index: [umolC/m2 leaf/s] * [m2 leaf/m2 ground]
      ! -> [umolC/m2 ground/s], accumulated down the canopy
      fsun = SunlitFraction(light_env%laisun_z(iv), light_env%laisha_z(iv))
      anet_layer = fsun*anet_sun + (1.0_r8 - fsun)*anet_sha
      canopy_anet_out = canopy_anet_out + anet_layer *                         &
        (light_env%laisun_z(iv) + light_env%laisha_z(iv))

      if (do_store) then
        parsun_z_out(iv, ilai_store)  = light_env%parsun_z(iv)
        parsha_z_out(iv, ilai_store)  = light_env%parsha_z(iv)
        laisun_z_out(iv, ilai_store)  = light_env%laisun_z(iv)
        laisha_z_out(iv, ilai_store)  = light_env%laisha_z(iv)
        nscaler_z_out(iv, ilai_store) = nscaler_z(iv)
        anet_z_out(iv, ilai_store)    = anet_layer
      end if

    end do

  end subroutine CanopyNetAssim

  ! ==========================================================================

  subroutine WriteOutput()
    !
    ! DESCRIPTION:
    ! Writes the prescribed LAI values and the per-leaf-layer light profile
    ! resolved at each of them to netCDF.

    ! LOCALS:
    integer           :: ncid          ! netcdf file id
    character(len=20) :: dim_names(2)  ! dimension names
    integer           :: dimIDs(2)     ! dimension IDs
    integer           :: laiID, layerID, nvID, canopyanetID
    integer           :: parsunID, parshaID, laisunID, laishaID
    integer           :: nscalerID, anetzID

    dim_names = [character(len=20) :: 'layer', 'lai']

    call OpenNCFile(trim(out_file), ncid, 'readwrite')
    call RegisterNCDims(ncid, dim_names, (/nlevleaf, n_lai/), 2, dimIDs)

    call RegisterVar(ncid, 'layer', dimIDs(1:1), type_int,                     &
      [character(len=20) :: 'units', 'long_name'],                             &
      [character(len=150) :: '-', 'leaf layer index, 1 = top of canopy'], 2, layerID)
    call RegisterVar(ncid, 'lai', dimIDs(2:2), type_double,                    &
      [character(len=20) :: 'units', 'long_name'],                             &
      [character(len=150) :: 'm2 m-2', 'prescribed canopy leaf area index'], 2, laiID)
    call RegisterVar(ncid, 'nv', dimIDs(2:2), type_int,                        &
      [character(len=20) :: 'units', 'long_name'],                             &
      [character(len=150) :: '-', 'number of occupied leaf layers at each prescribed LAI'], 2, nvID)
    call RegisterVar(ncid, 'canopy_anet', dimIDs(2:2), type_double,            &
      [character(len=20) :: 'units', 'long_name'],                             &
      [character(len=150) :: 'umolC m-2 s-1', 'canopy net photosynthesis per unit ground area'], 2, canopyanetID)

    call RegisterVar(ncid, 'parsun_z', (/dimIDs(1), dimIDs(2)/), type_double,  &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],              &
      [character(len=150) :: 'layer lai', 'W m-2', 'absorbed PAR per unit ground area, sunlit leaves'], 3, parsunID)
    call RegisterFillValue(ncid, parsunID, fates_unset_r8)
    call RegisterVar(ncid, 'parsha_z', (/dimIDs(1), dimIDs(2)/), type_double,  &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],              &
      [character(len=150) :: 'layer lai', 'W m-2', 'absorbed PAR per unit ground area, shaded leaves'], 3, parshaID)
    call RegisterFillValue(ncid, parshaID, fates_unset_r8)
    call RegisterVar(ncid, 'laisun_z', (/dimIDs(1), dimIDs(2)/), type_double,  &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],              &
      [character(len=150) :: 'layer lai', 'm2 m-2', 'sunlit leaf area index per layer'], 3, laisunID)
    call RegisterFillValue(ncid, laisunID, fates_unset_r8)
    call RegisterVar(ncid, 'laisha_z', (/dimIDs(1), dimIDs(2)/), type_double,  &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],              &
      [character(len=150) :: 'layer lai', 'm2 m-2', 'shaded leaf area index per layer'], 3, laishaID)
    call RegisterFillValue(ncid, laishaID, fates_unset_r8)
    call RegisterVar(ncid, 'nscaler_z', (/dimIDs(1), dimIDs(2)/), type_double, &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],              &
      [character(len=150) :: 'layer lai', '-', 'nitrogen-scaling factor per layer'], 3, nscalerID)
    call RegisterFillValue(ncid, nscalerID, fates_unset_r8)
    call RegisterVar(ncid, 'anet_z', (/dimIDs(1), dimIDs(2)/), type_double,    &
      [character(len=20) :: 'coordinates', 'units', 'long_name'],              &
      [character(len=150) :: 'layer lai', 'umolC m-2 s-1', 'area-weighted net photosynthesis per unit leaf area, per layer'], 3, anetzID)
    call RegisterFillValue(ncid, anetzID, fates_unset_r8)

    call EndNCDef(ncid)

    call WriteVar(ncid, layerID, layer_index(:))
    call WriteVar(ncid, laiID, lai_vals(:))
    call WriteVar(ncid, nvID, nv_out(:))
    call WriteVar(ncid, canopyanetID, canopy_anet(:))
    call WriteVar(ncid, parsunID, parsun_z_out(:,:))
    call WriteVar(ncid, parshaID, parsha_z_out(:,:))
    call WriteVar(ncid, laisunID, laisun_z_out(:,:))
    call WriteVar(ncid, laishaID, laisha_z_out(:,:))
    call WriteVar(ncid, nscalerID, nscaler_z_out(:,:))
    call WriteVar(ncid, anetzID, anet_z_out(:,:))

    call CloseNCFile(ncid)

  end subroutine WriteOutput

end program FatesCanopyLevelPhoto
