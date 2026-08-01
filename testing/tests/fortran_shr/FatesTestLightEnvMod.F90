module FatesTestLightEnvMod
  !
  ! DESCRIPTION:
  ! Prescribed light environment for standalone, patch-less/site-less cohort test
  ! drivers. Attenuates a prescribed incident PAR through a lone cohort's own leaf
  ! layers using FATES's two-stream radiation solver (TwoStreamMLPEMod), with a
  ! single scattering element (one canopy layer, one column, occupying 100% of its
  ! own footprint) standing in for the cohort - no fates_patch_type/ed_site_type is
  ! built or required.
  !
  ! Reference full-sun PAR, the direct/diffuse split, and the diurnal shape
  ! (including the matching idealized coszen arc) are assumptions with no existing
  ! precedent elsewhere in the repo to draw from; ground albedo is a typical
  ! soil/litter PAR value.
  !

  use FatesConstantsMod,   only : r8 => fates_r8
  use FatesConstantsMod,   only : pi_const
  use FatesConstantsMod,   only : wm2_to_umolm2s
  use EDParamsMod,         only : GetNVegLayers
  use FatesAllometryMod,   only : VegAreaLayer
  use FatesCohortMod,      only : fates_cohort_type
  use FatesRadiationMemMod, only : ivis
  use TwoStreamMLPEMod,    only : twostream_type
  use TwoStreamMLPEMod,    only : normalized_upper_boundary

  implicit none
  private

  ! ------------------------------------------------------------------------------------
  ! PRESCRIBED LIGHT ENVIRONMENT ASSUMPTIONS
  ! ------------------------------------------------------------------------------------
  real(r8), public, parameter :: ref_par_full_sun = 2000.0_r8/wm2_to_umolm2s ! reference full-sun incident PAR at solar noon [W/m2] (~2000 umol/m2/s)
  real(r8), parameter :: direct_frac       = 0.85_r8  ! fraction of incident PAR that is direct beam (typical clear sky)
  real(r8), parameter :: diffuse_frac      = 1.0_r8 - direct_frac ! fraction of incident PAR that is diffuse
  real(r8), parameter :: daylength_hours   = 12.0_r8  ! fixed photoperiod - no seasonal cycle
  real(r8), parameter :: sunrise_hour      = 12.0_r8 - 0.5_r8*daylength_hours ! local hour of sunrise
  real(r8), parameter :: sunset_hour       = 12.0_r8 + 0.5_r8*daylength_hours ! local hour of sunset
  real(r8), parameter :: coszen_noon       = 0.87_r8  ! prescribed cosine of the zenith angle at solar noon
  real(r8), parameter :: ground_albedo_par = 0.10_r8  ! soil/litter PAR albedo (diffuse and beam)
  real(r8), parameter :: frac_snow         = 0.0_r8   ! canopy snow-covered fraction (no snow)
  real(r8), parameter :: snow_depth        = 0.0_r8   ! physical snow depth [m] (no snow)

  type, public :: light_env_type

     private

     integer  :: pft            ! plant functional type index
     integer  :: nv              ! number of occupied leaf layers
     real(r8) :: treelai         ! cached cohort total leaf area index [m2/m2]
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
     procedure, public :: Free

  end type light_env_type

contains

  ! ==========================================================================

  subroutine Init(this, cohort, pft)
    !
    ! DESCRIPTION:
    ! Allocate the two-stream object for a lone cohort, set up its single
    ! scattering element and ground albedo, and allocate the per-leaf-layer arrays.

    ! ARGUMENTS:
    class(light_env_type),   intent(inout) :: this   ! light environment object
    type(fates_cohort_type), intent(in)    :: cohort ! the cohort this light environment tracks
    integer,                 intent(in)    :: pft    ! plant functional type index

    this%pft     = pft
    this%treelai = cohort%treelai
    this%treesai = cohort%treesai
    this%height  = cohort%height
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

  subroutine Refresh(this, cohort)
    !
    ! DESCRIPTION:
    ! Re-sync the scattering element's canopy structure from the cohort's current
    ! treelai/treesai/height. The two-stream element was otherwise only ever built
    ! once, at recruitment - this is the fix for crown structure silently going
    ! stale once PRT allocation starts changing leaf area daily.

    ! ARGUMENTS:
    class(light_env_type),   intent(inout) :: this   ! light environment object
    type(fates_cohort_type), intent(in)    :: cohort ! the cohort this light environment tracks

    ! LOCALS:
    integer :: nv_new ! number of occupied leaf layers, recomputed from current lai/sai

    this%treelai = cohort%treelai
    this%treesai = cohort%treesai
    this%height  = cohort%height

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

  subroutine Profile(this, light_frac, hour_of_day)
    !
    ! DESCRIPTION:
    ! Prescribe incident PAR at the top of the crown for the given light fraction
    ! and hour of day, attenuate it through the cohort's own leaf layers via the
    ! two-stream solver, and fill parsun_z/parsha_z/laisun_z/laisha_z. Skips the
    ! solve and zeros the arrays outright when there is no incident PAR, so coszen
    ! never reaches ZenithPrep at zero.

    ! ARGUMENTS:
    class(light_env_type), intent(inout) :: this        ! light environment object
    real(r8),               intent(in)   :: light_frac  ! incident light fraction [0-1]
    real(r8),               intent(in)   :: hour_of_day ! hour of day [0-24]

    ! LOCALS:
    real(r8) :: diurnal_shape ! prescribed diurnal light shape [0-1]
    real(r8) :: coszen        ! prescribed cosine of the solar zenith angle [0-1]
    real(r8) :: par_toc       ! total incident PAR at the top of the crown [W/m2]
    real(r8) :: par_beam      ! direct-beam incident PAR at the top of the crown [W/m2]
    real(r8) :: par_diff      ! diffuse incident PAR at the top of the crown [W/m2]
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

    if (hour_of_day > sunrise_hour .and. hour_of_day < sunset_hour) then
      diurnal_shape = sin(pi_const * (hour_of_day - sunrise_hour)/daylength_hours)
    else
      diurnal_shape = 0.0_r8
    end if

    par_toc = light_frac * ref_par_full_sun * diurnal_shape

    ! No incident PAR: skip the solve entirely rather than calling ZenithPrep with a
    ! coszen of zero, and treat all leaf area as shaded (there is no sun to define a
    ! sunlit fraction against).
    if (par_toc <= 0.0_r8) then
      this%parsun_z(:) = 0.0_r8
      this%parsha_z(:) = 0.0_r8
      this%laisun_z(:) = 0.0_r8
      do iv = 1, this%nv
        call VegAreaLayer(this%treelai, this%treesai, this%height, iv, this%nv,     &
          this%pft, snow_depth, vai_top, vai_bot, elai_layer, esai_layer)
        this%laisha_z(iv) = elai_layer
      end do
      return
    end if

    par_beam = direct_frac * par_toc
    par_diff = diffuse_frac * par_toc
    coszen   = coszen_noon * diurnal_shape

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
      this%parsun_z(iv) = Rb_abs_leaf + Rd_abs_leaf*leaf_sun_frac
      this%parsha_z(iv) = Rd_abs_leaf*(1.0_r8 - leaf_sun_frac)
      this%laisun_z(iv) = elai_layer*leaf_sun_frac
      this%laisha_z(iv) = elai_layer*(1.0_r8 - leaf_sun_frac)
    end do

  end subroutine Profile

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
