module FatesTestLightEnvMod
  !
  ! DESCRIPTION:
  ! Prescribed light environment for functional tests. Attenuates a prescribed incident 
  ! PAR through a single canopy's own leaf layers using FATES's two-stream radiation solver (TwoStreamMLPEMod)
  !
  ! Reference full-sun PAR (at cosz=1) and the direct/diffuse split are
  ! assumptions. Ground albedo is a typical soil/litter PAR value.
  !  
  ! The diurnal/annual cycle is derived using solar declination from day of year 
  ! (Cooper 1969 single-term sinusoidal approximation of Earth's obliquity, 
  ! no eccentricity/perihelion correction), and coszen(hour) from latitude and declination 
  ! via the standard hour-angle formula
  !

  use FatesConstantsMod,    only : r8 => fates_r8
  use EDParamsMod,          only : GetNVegLayers
  use FatesAllometryMod,    only : VegAreaLayer
  use FatesRadiationMemMod, only : ivis
  use TwoStreamMLPEMod,     only : twostream_type
  use TwoStreamMLPEMod,     only : normalized_upper_boundary

  implicit none
  private

  ! ------------------------------------------------------------------------------------
  ! PRESCRIBED LIGHT ENVIRONMENT ASSUMPTIONS
  ! ------------------------------------------------------------------------------------
  real(r8), parameter :: ground_albedo_par = 0.10_r8 ! soil/litter PAR albedo (diffuse and beam)
  real(r8), parameter :: frac_snow         = 0.0_r8  ! canopy snow-covered fraction (no snow)
  real(r8), parameter :: snow_depth        = 0.0_r8  ! physical snow depth [m] (no snow)

  type, public :: light_env_type

     private

     integer              :: pft     ! plant functional type index
     integer              :: nv      ! number of occupied leaf layers
     real(r8), public     :: treelai ! saved in-crown total leaf area index [m2 leaf/m2 crown footprint]
     real(r8)             :: treesai ! saved in-crown total stem area index [m2 stem/m2 crown footprint]
     real(r8)             :: height  ! saved cohort height [m]
     type(twostream_type) :: twostr  ! two-stream radiation object (one element)

     real(r8), public, allocatable :: parsun_z(:) ! absorbed PAR, sunlit leaves, per leaf layer [W/m2 crown footprint]
     real(r8), public, allocatable :: parsha_z(:) ! absorbed PAR, shaded leaves, per leaf layer [W/m2 crown footprint]
     real(r8), public, allocatable :: laisun_z(:) ! sunlit LAI per leaf layer [m2 leaf/m2 crown footprint]
     real(r8), public, allocatable :: laisha_z(:) ! shaded LAI per leaf layer [m2 leaf/m2 crown footprint]

   contains

     procedure, public :: Init
     procedure, public :: Refresh
     procedure, public :: AttenuateCanopy
     procedure, public :: Free

  end type light_env_type

contains

  ! ==========================================================================

  subroutine Init(this, treelai, treesai, height, pft)
    !
    ! DESCRIPTION:
    ! Allocate the two-stream object for a canopy, set up its single
    ! scattering element and ground albedo, and allocate the per-leaf-layer arrays.
    !

    ! ARGUMENTS:
    class(light_env_type),   intent(inout) :: this    ! light environment object
    real(r8),                intent(in)    :: treelai ! in-crown leaf area index [m2 leaf/m2 crown footprint]
    real(r8),                intent(in)    :: treesai ! in-crown stem area index [m2 stem/m2 crown footprint]
    real(r8),                intent(in)    :: height  ! plant/canopy height [m]
    integer,                 intent(in)    :: pft     ! plant functional type index

    this%pft = pft
    this%treelai = treelai
    this%treesai = treesai
    this%height = height
    this%nv = GetNVegLayers(this%treelai + this%treesai)

    call this%twostr%AllocInitTwoStream((/ivis/), 1, 1)
    this%twostr%scelg(1,1)%pft = this%pft
    this%twostr%scelg(1,1)%area = 1.0_r8
    this%twostr%scelg(1,1)%lai = this%treelai
    this%twostr%scelg(1,1)%sai = this%treesai
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
    ! treelai/treesai/height.

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
  
  subroutine AttenuateCanopy(this, par_beam, par_diff, coszen, parsun_z, parsha_z, laisun_z, laisha_z)
    !
    ! DESCRIPTION:
    ! Attenuate given incident beam/diffuse PAR (par_beam/par_diff, at the given
    ! coszen) through this canopy's leaf layers via the two-stream solver,
    ! filling the caller-supplied parsun_z/parsha_z/laisun_z/laisha_z arrays

    ! ARGUMENTS:
    class(light_env_type),  intent(inout) :: this        ! light environment object
    real(r8),               intent(in)    :: par_beam    ! direct-beam incident PAR at the top of the crown [W/m2]
    real(r8),               intent(in)    :: par_diff    ! diffuse incident PAR at the top of the crown [W/m2]
    real(r8),               intent(in)    :: coszen      ! cosine of the solar zenith angle to use for this solve [-1 to 1]; <=0 is treated as no incident PAR
    real(r8),               intent(out)   :: parsun_z(:) ! absorbed PAR, sunlit leaves, per leaf layer [W/m2 crown footprint]
    real(r8),               intent(out)   :: parsha_z(:) ! absorbed PAR, shaded leaves, per leaf layer [W/m2 crown footprint]
    real(r8),               intent(out)   :: laisun_z(:) ! sunlit LAI per leaf layer [m2 leaf/m2 crown footprint]
    real(r8),               intent(out)   :: laisha_z(:) ! shaded LAI per leaf layer [m2 leaf/m2 crown footprint]

    ! LOCALS:
    real(r8) :: vai_top, vai_bot         ! vegetation area index bounds of the current leaf layer
    real(r8) :: elai_layer, esai_layer   ! exposed leaf/stem area index of the current leaf layer
    real(r8) :: Rb_abs, Rd_abs           ! total absorbed beam/diffuse radiation, current layer [W/m2 crown footprint]
    real(r8) :: Rb_abs_leaf, Rd_abs_leaf ! absorbed beam/diffuse radiation from leaves, current layer [W/m2 crown footprint]
    real(r8) :: R_abs_stem, R_abs_snow   ! absorbed radiation from stems/snow, current layer [W/m2 crown footprint]
    real(r8) :: leaf_sun_frac            ! sunlit fraction of leaf area in the current layer
    real(r8) :: taulamb(2)               ! two-stream solve scratch space (size = 2*n_scel = 2)
    real(r8) :: omega(2,2)               ! two-stream solve scratch space
    real(r8) :: albedo_beam              ! two-stream solve outputs (unused diagnostics)
    real(r8) :: albedo_diff, consv_err   ! two-stream solve outputs (unused diagnostics)
    real(r8) :: frac_abs_beam            ! two-stream solve outputs (unused diagnostics)
    real(r8) :: frac_abs_diff            ! two-stream solve outputs (unused diagnostics)
    real(r8) :: frac_beam_grnd           ! two-stream solve outputs (unused diagnostics)
    real(r8) :: frac_diff_grnd_beam      ! two-stream solve outputs (unused diagnostics)
    real(r8) :: frac_diff_grnd_diff      ! two-stream solve outputs (unused diagnostics)
    integer  :: ipiv(2)                  ! two-stream solve scratch space (LAPACK pivots)
    integer  :: iv                       ! leaf-layer looping index
    logical  :: call_fail                ! GetAbsRad failure flag

    ! no incident PAR: skip the solve and treat all leaf area as shaded
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
      parsun_z(iv) = Rb_abs_leaf + Rd_abs_leaf*leaf_sun_frac
      parsha_z(iv) = Rd_abs_leaf*(1.0_r8 - leaf_sun_frac)
      laisun_z(iv) = elai_layer*leaf_sun_frac
      laisha_z(iv) = elai_layer*(1.0_r8 - leaf_sun_frac)
    end do

  end subroutine AttenuateCanopy

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
