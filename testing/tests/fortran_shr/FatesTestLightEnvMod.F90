module FatesTestLightEnvMod
  !
  ! DESCRIPTION:
  ! Prescribed light environment for standalone, patch-less/site-less test
  ! drivers. Attenuates a prescribed incident PAR through a canopy's own leaf
  ! layers using FATES's two-stream radiation solver (TwoStreamMLPEMod), with a
  ! single scattering element (one canopy layer, one column, occupying 100% of its
  ! own footprint) standing in for that canopy - no fates_patch_type/ed_site_type
  ! is built or required.
  !
  ! Canopy structure enters as plain treelai/treesai/height scalars (Init/
  ! Refresh), not a fates_cohort_type, so this serves both a cohort-driven
  ! caller whose leaf area follows from allometry (test_SingleCohort.F90) and
  ! one whose LAI is prescribed outright as an experimental treatment
  ! (test_CanopyLevelPhoto.F90) - see Init's header comment.
  !
  ! Reference full-sun PAR (at cosz=1) and the direct/diffuse split are
  ! assumptions with no existing precedent elsewhere in the repo to draw from;
  ! ground albedo is a typical soil/litter PAR value. The diurnal/annual cycle
  ! itself - solar declination from day of year (Cooper 1969 single-term
  ! sinusoidal approximation of Earth's obliquity, no eccentricity/perihelion
  ! correction), and coszen(hour) from latitude and declination via the
  ! standard hour-angle formula - follows real solar geometry, driven by the
  ! site latitude prescribed in FatesTestSiteMod (shared with
  ! FatesTestEnvironmentMod's temperature cycle, so both respond consistently
  ! to a "different site" experiment).
  !

  use FatesConstantsMod,   only : r8 => fates_r8
  use FatesConstantsMod,   only : pi_const
  use FatesConstantsMod,   only : rad_per_deg
  use FatesConstantsMod,   only : wm2_to_umolm2s
  use EDParamsMod,         only : GetNVegLayers
  use FatesAllometryMod,   only : VegAreaLayer
  use FatesRadiationMemMod, only : ivis
  use TwoStreamMLPEMod,    only : twostream_type
  use TwoStreamMLPEMod,    only : normalized_upper_boundary
  use FatesTestSiteMod,    only : latitude_deg

  implicit none
  private

  ! ------------------------------------------------------------------------------------
  ! PRESCRIBED LIGHT ENVIRONMENT ASSUMPTIONS
  ! ------------------------------------------------------------------------------------
  real(r8), public, parameter :: ref_par_full_sun = 2000.0_r8/wm2_to_umolm2s ! reference full-sun incident PAR at cosz=1 [W/m2] (~2000 umol/m2/s)
  ! public so a driver prescribing its own incident PAR outright (rather than
  ! going through Profile's light-fraction/solar-cycle path) can reuse this
  ! same clear-sky split instead of introducing a second, conflicting one -
  ! see test_CanopyLevelPhoto.F90
  real(r8), public, parameter :: direct_frac  = 0.85_r8  ! fraction of incident PAR that is direct beam (typical clear sky)
  real(r8), public, parameter :: diffuse_frac = 1.0_r8 - direct_frac ! fraction of incident PAR that is diffuse
  real(r8), parameter :: max_declin_deg    = 23.45_r8 ! Earth's obliquity, used as the declination amplitude (Cooper 1969) [deg]
  real(r8), parameter :: ground_albedo_par = 0.10_r8  ! soil/litter PAR albedo (diffuse and beam)
  real(r8), parameter :: frac_snow         = 0.0_r8   ! canopy snow-covered fraction (no snow)
  real(r8), parameter :: snow_depth        = 0.0_r8   ! physical snow depth [m] (no snow)

  type, public :: light_env_type

     private

     integer  :: pft            ! plant functional type index
     integer  :: nv              ! number of occupied leaf layers
     real(r8), public :: treelai         ! cached cohort total leaf area index [m2/m2]
     real(r8) :: treesai         ! cached cohort total stem area index [m2/m2]
     real(r8) :: height          ! cached cohort height [m]
     type(twostream_type) :: twostr ! two-stream radiation object (one element)

     real(r8), public, allocatable :: parsun_z(:) ! absorbed PAR, sunlit leaves, per leaf layer [W/m2 ground]
     real(r8), public, allocatable :: parsha_z(:) ! absorbed PAR, shaded leaves, per leaf layer [W/m2 ground]
     real(r8), public, allocatable :: laisun_z(:) ! sunlit LAI per leaf layer [m2/m2]
     real(r8), public, allocatable :: laisha_z(:) ! shaded LAI per leaf layer [m2/m2]

   contains

     procedure, public :: Init
     procedure, public :: Refresh
     procedure, public :: Profile
     procedure, public :: AttenuateCanopy
     procedure, public :: Free

  end type light_env_type

contains

  ! ==========================================================================

  subroutine Init(this, treelai, treesai, height, pft)
    !
    ! DESCRIPTION:
    ! Allocate the two-stream object for a single canopy, set up its single
    ! scattering element and ground albedo, and allocate the per-leaf-layer arrays.
    !
    ! Canopy structure is taken as plain scalars rather than a
    ! fates_cohort_type: nothing in this module needs anything from a cohort
    ! beyond treelai/treesai/height, and taking them directly is what lets a
    ! driver with a PRESCRIBED canopy and no cohort at all
    ! (test_CanopyLevelPhoto.F90, whose LAI is an experimental treatment
    ! rather than an allometric consequence of a dbh) reuse this identical
    ! two-stream attenuation rather than reimplementing it. A cohort-driven
    ! caller simply passes cohort%treelai/treesai/height (see
    ! test_SingleCohort.F90).

    ! ARGUMENTS:
    class(light_env_type),   intent(inout) :: this    ! light environment object
    real(r8),                intent(in)    :: treelai ! in-crown leaf area index [m2 leaf/m2 crown footprint]
    real(r8),                intent(in)    :: treesai ! in-crown stem area index [m2 stem/m2 crown footprint]
    real(r8),                intent(in)    :: height  ! plant/canopy height [m]
    integer,                 intent(in)    :: pft     ! plant functional type index

    this%pft     = pft
    this%treelai = treelai
    this%treesai = treesai
    this%height  = height
    this%nv      = GetNVegLayers(this%treelai + this%treesai)

    call this%twostr%AllocInitTwoStream((/ivis/), 1, 1)
    this%twostr%scelg(1,1)%pft  = this%pft
    this%twostr%scelg(1,1)%area = 1.0_r8
    this%twostr%scelg(1,1)%lai  = this%treelai
    this%twostr%scelg(1,1)%sai  = this%treesai
    this%twostr%n_col(1) = 1
    call this%twostr%GetNSCel()
    this%twostr%band(ivis)%albedo_grnd_diff = ground_albedo_par
    this%twostr%band(ivis)%albedo_grnd_beam = ground_albedo_par
    this%twostr%force_prep = .true.
    call this%twostr%CanopyPrep(frac_snow)

    allocate(this%parsun_z(this%nv), this%parsha_z(this%nv))
    allocate(this%laisun_z(this%nv), this%laisha_z(this%nv))

  end subroutine Init

  ! ==========================================================================

  subroutine Refresh(this, treelai, treesai, height)
    !
    ! DESCRIPTION:
    ! Re-sync the scattering element's canopy structure to the caller's current
    ! treelai/treesai/height. The two-stream element was otherwise only ever built
    ! once, at recruitment - this is the fix for crown structure silently going
    ! stale once PRT allocation starts changing leaf area daily. Takes scalars
    ! rather than a cohort, for the same reason Init does (see its header).

    ! ARGUMENTS:
    class(light_env_type),   intent(inout) :: this    ! light environment object
    real(r8),                intent(in)    :: treelai ! in-crown leaf area index [m2 leaf/m2 crown footprint]
    real(r8),                intent(in)    :: treesai ! in-crown stem area index [m2 stem/m2 crown footprint]
    real(r8),                intent(in)    :: height  ! plant/canopy height [m]

    ! LOCALS:
    integer :: nv_new ! number of occupied leaf layers, recomputed from current lai/sai

    this%treelai = treelai
    this%treesai = treesai
    this%height  = height

    this%twostr%scelg(1,1)%lai = this%treelai
    this%twostr%scelg(1,1)%sai = this%treesai
    this%twostr%force_prep = .true.
    call this%twostr%CanopyPrep(frac_snow)

    nv_new = GetNVegLayers(this%treelai + this%treesai)
    if (nv_new /= this%nv) then
      deallocate(this%parsun_z, this%parsha_z, this%laisun_z, this%laisha_z)
      this%nv = nv_new
      allocate(this%parsun_z(this%nv), this%parsha_z(this%nv))
      allocate(this%laisun_z(this%nv), this%laisha_z(this%nv))
    end if

  end subroutine Refresh

  ! ==========================================================================

  subroutine Profile(this, light_frac, day_of_year, hour_of_day, par_toc_out)
    !
    ! DESCRIPTION:
    ! Prescribe incident PAR at the top of the crown for the given light fraction,
    ! day of year, and hour of day, and attenuate it through the cohort's own leaf
    ! layers via AttenuateCanopy, filling this%parsun_z/parsha_z/laisun_z/laisha_z
    ! (this type's own persistent per-substep state). par_toc_out optionally
    ! returns the incident PAR this substep actually used, for callers (e.g. the
    ! driver's light-interception-efficiency diagnostic) that need it alongside
    ! the attenuated profile rather than re-deriving it themselves.

    ! ARGUMENTS:
    class(light_env_type), intent(inout) :: this        ! light environment object
    real(r8),               intent(in)   :: light_frac  ! incident light fraction [0-1]
    integer,                intent(in)   :: day_of_year ! day of year [1-365], drives the solar declination
    real(r8),               intent(in)   :: hour_of_day ! hour of day [0-24]
    real(r8),               intent(out), optional :: par_toc_out ! total incident PAR at the top of the crown [W/m2]

    ! LOCALS:
    real(r8) :: declin   ! solar declination angle, today [radians]
    real(r8) :: coszen   ! cosine of the solar zenith angle [-1 to 1]; <=0 is night
    real(r8) :: par_toc  ! total incident PAR at the top of the crown [W/m2]
    real(r8) :: par_beam ! direct-beam incident PAR at the top of the crown [W/m2]
    real(r8) :: par_diff ! diffuse incident PAR at the top of the crown [W/m2]

    declin = SolarDeclination(day_of_year)
    coszen = CosSolarZenith(declin, hour_of_day)
    ! clamped at zero: coszen<0 (sun below horizon) means no incident PAR, not
    ! negative incident PAR - AttenuateCanopy's own par_beam+par_diff<=0 check
    ! already handles this correctly internally regardless, but par_toc_out is
    ! returned as-is to callers, so the physical quantity must be clamped here
    ! rather than left to go negative
    par_toc = max(0.0_r8, light_frac * ref_par_full_sun * coszen)
    par_beam = direct_frac * par_toc
    par_diff = diffuse_frac * par_toc
    if (present(par_toc_out)) par_toc_out = par_toc

    call this%AttenuateCanopy(par_beam, par_diff, coszen, this%parsun_z,        &
      this%parsha_z, this%laisun_z, this%laisha_z)

  end subroutine Profile

  ! ==========================================================================

  subroutine AttenuateCanopy(this, par_beam, par_diff, coszen, parsun_z, parsha_z, laisun_z, laisha_z)
    !
    ! DESCRIPTION:
    ! Attenuate given incident beam/diffuse PAR (par_beam/par_diff, at the given
    ! coszen) through this cohort's own leaf layers via the two-stream solver,
    ! filling the caller-supplied parsun_z/parsha_z/laisun_z/laisha_z arrays
    ! (size this%nv). Factored out of Profile - which derives par_beam/par_diff
    ! from (light_frac, day_of_year, hour_of_day) via the fixed direct_frac/
    ! diffuse_frac split - so a diagnostic instantaneous light-response sweep
    ! (an arbitrary incident PPFD, at an arbitrary beam/diffuse split and
    ! coszen, independent of light_frac/day_of_year/hour_of_day) can reuse the
    ! identical attenuation physics without touching this%parsun_z/parsha_z/
    ! laisun_z/laisha_z at all - nothing to restore there. this%twostr's
    ! internal solved state is still touched (ZenithPrep/Solve are not
    ! side-effect-free), but Profile unconditionally re-solves it fresh on every
    ! call regardless of what it was left at, so a caller that needs the real
    ! per-substep state restored afterward should simply call Profile again with
    ! the real light_frac/day_of_year/hour_of_day once done.

    ! ARGUMENTS:
    class(light_env_type), intent(inout) :: this        ! light environment object
    real(r8),               intent(in)   :: par_beam     ! direct-beam incident PAR at the top of the crown [W/m2]
    real(r8),               intent(in)   :: par_diff     ! diffuse incident PAR at the top of the crown [W/m2]
    real(r8),               intent(in)   :: coszen       ! cosine of the solar zenith angle to use for this solve [-1 to 1]; <=0 is treated as no incident PAR
    real(r8),               intent(out)  :: parsun_z(:)  ! absorbed PAR, sunlit leaves, per leaf layer [W/m2 ground]
    real(r8),               intent(out)  :: parsha_z(:)  ! absorbed PAR, shaded leaves, per leaf layer [W/m2 ground]
    real(r8),               intent(out)  :: laisun_z(:)  ! sunlit LAI per leaf layer [m2/m2]
    real(r8),               intent(out)  :: laisha_z(:)  ! shaded LAI per leaf layer [m2/m2]

    ! LOCALS:
    real(r8) :: vai_top, vai_bot           ! vegetation area index bounds of the current leaf layer
    real(r8) :: elai_layer, esai_layer     ! exposed leaf/stem area index of the current leaf layer
    real(r8) :: Rb_abs, Rd_abs             ! total absorbed beam/diffuse radiation, current layer [W/m2 ground]
    real(r8) :: Rb_abs_leaf, Rd_abs_leaf   ! absorbed beam/diffuse radiation from leaves, current layer [W/m2 ground]
    real(r8) :: R_abs_stem, R_abs_snow     ! absorbed radiation from stems/snow, current layer [W/m2 ground]
    real(r8) :: leaf_sun_frac              ! sunlit fraction of leaf area in the current layer
    logical  :: call_fail                  ! GetAbsRad failure flag
    real(r8) :: taulamb(2)     ! two-stream solve scratch space (size = 2*n_scel = 2)
    real(r8) :: omega(2,2)     ! two-stream solve scratch space
    integer  :: ipiv(2)        ! two-stream solve scratch space (LAPACK pivots)
    real(r8) :: albedo_beam, albedo_diff, consv_err                      ! two-stream solve outputs (unused diagnostics)
    real(r8) :: frac_abs_beam, frac_abs_diff                             ! two-stream solve outputs (unused diagnostics)
    real(r8) :: frac_beam_grnd, frac_diff_grnd_beam, frac_diff_grnd_diff ! two-stream solve outputs (unused diagnostics)
    integer  :: iv ! leaf-layer looping index

    ! No incident PAR: skip the solve entirely rather than calling ZenithPrep with a
    ! non-positive coszen, and treat all leaf area as shaded (there is no sun to
    ! define a sunlit fraction against).
    if (par_beam + par_diff <= 0.0_r8 .or. coszen <= 0.0_r8) then
      parsun_z(:) = 0.0_r8
      parsha_z(:) = 0.0_r8
      laisun_z(:) = 0.0_r8
      do iv = 1, this%nv
        call VegAreaLayer(this%treelai, this%treesai, this%height, iv, this%nv,     &
          this%pft, snow_depth, vai_top, vai_bot, elai_layer, esai_layer)
        laisha_z(iv) = elai_layer
      end do
      return
    end if

    ! ZenithPrep + Solve are always run with the literal normalized boundary
    ! (Rbeam_atm=Rdiff_atm=1) - Solve's internal energy-conservation check is
    ! hardwired to that convention (see production usage at
    ! FatesRadiationDriveMod.F90:170-183). This produces normalized scattering
    ! coefficients (per unit incident beam/diffuse) that do not themselves depend on
    ! the real light magnitude. The real incident PAR is applied afterward by
    ! assigning it directly to twostr%band(ivis)%Rbeam_atm/Rdiff_atm (bypassing
    ! Solve entirely for that step, mirroring FatesRadiationDriveMod.F90:363-364),
    ! which is what GetAbsRad actually reads to scale the normalized profile to
    ! real absorbed radiation.
    call this%twostr%ZenithPrep(coszen)
    call this%twostr%Solve(ivis, normalized_upper_boundary, 1.0_r8, 1.0_r8,         &
      taulamb, omega, ipiv, albedo_beam, albedo_diff, consv_err,                    &
      frac_abs_beam, frac_abs_diff, frac_beam_grnd, frac_diff_grnd_beam,            &
      frac_diff_grnd_diff)

    this%twostr%band(ivis)%Rbeam_atm = par_beam
    this%twostr%band(ivis)%Rdiff_atm = par_diff

    do iv = 1, this%nv
      call VegAreaLayer(this%treelai, this%treesai, this%height, iv, this%nv,       &
        this%pft, snow_depth, vai_top, vai_bot, elai_layer, esai_layer)
      call this%twostr%GetAbsRad(1, 1, ivis, vai_top, vai_bot, Rb_abs, Rd_abs,       &
        Rd_abs_leaf, Rb_abs_leaf, R_abs_stem, R_abs_snow, leaf_sun_frac, call_fail)
      ! per unit GROUND area, matching GetAbsRad's own convention and how
      ! FatesPlantRespPhotosynthMod.F90 (par_per_sunla/par_per_shala, ~line 668)
      ! consumes it - NOT per unit leaf area. The per-leaf-area, W/m2->umol/m2/s
      ! conversion for LeafLayerPhotosynthesis's par_abs happens later, at the
      ! photosynthesis call site (see FatesPlantRespPhotosynthMod.F90's ConvertPar),
      ! not here.
      parsun_z(iv) = Rb_abs_leaf + Rd_abs_leaf*leaf_sun_frac
      parsha_z(iv) = Rd_abs_leaf*(1.0_r8 - leaf_sun_frac)
      laisun_z(iv) = elai_layer*leaf_sun_frac
      laisha_z(iv) = elai_layer*(1.0_r8 - leaf_sun_frac)
    end do

  end subroutine AttenuateCanopy

  ! ==========================================================================

  function SolarDeclination(day_of_year) result(declin)
    !
    ! DESCRIPTION:
    ! Solar declination angle for the given day of year, via the Cooper (1969)
    ! single-term sinusoidal approximation: a fixed present-day obliquity
    ! (max_declin_deg) modulated by a sine wave phased so the northern-hemisphere
    ! winter solstice falls near day 355 (day 284 + 10 days to solstice at the
    ! 1/365 day-angle scale). No orbital eccentricity/perihelion correction.

    ! ARGUMENTS:
    integer, intent(in) :: day_of_year ! day of year [1-365]

    ! RESULT:
    real(r8) :: declin ! solar declination angle [radians]

    declin = max_declin_deg * rad_per_deg *                                      &
      sin(rad_per_deg * (360.0_r8/365.0_r8) * (284.0_r8 + real(day_of_year, r8)))

  end function SolarDeclination

  ! ==========================================================================

  function CosSolarZenith(declin, hour_of_day) result(coszen)
    !
    ! DESCRIPTION:
    ! Cosine of the solar zenith angle from the prescribed site latitude, today's
    ! solar declination, and hour of day, via the standard hour-angle formula
    ! (equivalent to CIME's shr_orb_cosz, share/src/shr_orb_mod.F90, re-expressed
    ! in a local-solar-time hour angle measured from solar noon rather than a
    ! Julian-day fraction + longitude). May be negative (sun below the horizon);
    ! callers must not pass a non-positive result to ZenithPrep.

    ! ARGUMENTS:
    real(r8), intent(in) :: declin      ! solar declination angle [radians]
    real(r8), intent(in) :: hour_of_day ! hour of day [0-24]

    ! RESULT:
    real(r8) :: coszen ! cosine of the solar zenith angle [-1 to 1]

    ! LOCALS:
    real(r8) :: lat_rad     ! prescribed site latitude [radians]
    real(r8) :: hour_angle  ! solar hour angle, zero at solar noon [radians]

    lat_rad = latitude_deg * rad_per_deg
    hour_angle = (pi_const/12.0_r8) * (hour_of_day - 12.0_r8)

    coszen = sin(lat_rad)*sin(declin) + cos(lat_rad)*cos(declin)*cos(hour_angle)

  end function CosSolarZenith

  ! ==========================================================================

  subroutine Free(this)
    !
    ! DESCRIPTION:
    ! Tear down the two-stream object and per-leaf-layer arrays.

    ! ARGUMENTS:
    class(light_env_type), intent(inout) :: this ! light environment object

    call this%twostr%DeallocTwoStream()
    deallocate(this%parsun_z, this%parsha_z, this%laisun_z, this%laisha_z)

  end subroutine Free

end module FatesTestLightEnvMod
