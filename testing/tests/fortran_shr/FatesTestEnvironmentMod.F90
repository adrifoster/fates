module FatesTestEnvironmentMod
  !
  ! DESCRIPTION:
  ! Prescribed atmospheric and soil boundary conditions for standalone, patch-less/
  ! site-less cohort test drivers. Held fixed for the entire run by Init below - the
  ! only thing that varies in these drivers is light (across levels, and diurnally
  ! within a day; see FatesTestLightEnvMod). Photosynthesis (LeafLayerPhotosynthesis/
  ! LeafLayerBiophysicalRates) and non-leaf maintenance respiration read these
  ! directly. Bundling them into a type means a future seasonal cycle (e.g. a
  ! seasonally varying temperature or soil moisture) only touches this module.
  !

  use FatesConstantsMod, only : r8 => fates_r8

  implicit none
  private

  type, public :: environment_type

    real(r8) :: tempk         ! vegetation/leaf temperature [K]
    real(r8) :: can_press     ! air pressure at the leaf surface [Pa]
    real(r8) :: can_co2_ppress ! CO2 partial pressure at the leaf surface [Pa]
    real(r8) :: can_o2_ppress  ! O2 partial pressure at the leaf surface [Pa]
    real(r8) :: veg_esat      ! saturation vapor pressure at tempk [Pa]
    real(r8) :: can_vpress    ! vapor pressure of the canopy air [Pa]
    real(r8) :: gb            ! leaf boundary layer conductance [umol/m2/s]
    real(r8) :: btran         ! soil moisture stress factor [0-1]
    real(r8) :: dayl_factor   ! day-length photosynthetic-capacity acclimation factor [0-1]
    real(r8) :: t_growth      ! 10-day running-mean growth temperature [K]
    real(r8) :: t_home        ! long-term running-mean home temperature [K]
    real(r8) :: rootfr_ft     ! this pft's root fraction in the (single, implicit) soil layer [0-1]
    real(r8) :: t_soil        ! soil temperature [K]
    integer  :: nlevsoil      ! number of (implicit) soil layers

  contains

    procedure, public :: Init

  end type environment_type

contains

  ! ==========================================================================

  subroutine Init(this)
    !
    ! DESCRIPTION:
    ! Set the prescribed atmospheric and soil boundary conditions.

    ! ARGUMENTS:
    class(environment_type), intent(out) :: this ! environment object

    this%tempk         = 298.15_r8   ! [K] (25 C)
    this%can_press     = 101325.0_r8 ! [Pa] (sea level)
    this%can_co2_ppress = 42.6_r8    ! [Pa] (~420 ppm at can_press)
    this%can_o2_ppress  = 21177.0_r8 ! [Pa] (20.9% at can_press)
    this%veg_esat      = 3169.0_r8   ! [Pa] (saturation vapor pressure at tempk)
    this%can_vpress    = 2218.3_r8   ! [Pa] (assumed 70% RH at tempk)
    this%gb            = 2.0e6_r8    ! [umol/m2/s] (well-ventilated leaf, ~20 s/m equivalent)
    this%dayl_factor   = 1.0_r8      ! [0-1] (no seasonal daylength change assumed - distinct from the diurnal light cycle itself)
    this%t_growth      = this%tempk  ! (steady climate assumed)
    this%t_home        = this%tempk  ! (steady climate assumed)
    this%btran         = 1.0_r8      ! [0-1] (non-limiting water assumed)
    this%rootfr_ft     = 1.0_r8      ! [0-1]
    this%t_soil        = this%tempk  ! (same steady-climate assumption as tempk)
    this%nlevsoil      = 1           ! matches rootfr_ft holding 100% of roots in one layer

  end subroutine Init

end module FatesTestEnvironmentMod
