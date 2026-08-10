module FatesTestEnvironmentMod
  !
  ! DESCRIPTION:
  ! Prescribed atmospheric and soil boundary conditions for standalone, patch-less/
  ! site-less cohort test drivers. Everything except temperature is held fixed for
  ! the entire run by Init below. Temperature follows a real annual + diurnal cycle
  ! (SetHour, called once per substep alongside FatesTestLightEnvMod's Profile),
  ! parameterized by FatesTestSiteMod's t_annual_mean/t_annual_amp/
  ! t_annual_peak_doy/t_diurnal_amp/t_diurnal_peak_hour (fit to 2003-2016 hourly
  ! TBOT from the BCI DATM forcing set, ~/Documents/02_Projects/06_FATES/BCI_datm
  ! - see that module's header for the fitting/phase-alignment method): a
  ! single-term (Cooper-style) sinusoid for the annual cycle, and a separate
  ! single-term sinusoid for the diurnal anomaly. The diurnal cycle at BCI is
  ! actually asymmetric (fast morning rise, slower overnight decline, typical of
  ! a convective tropical climate) and a single sinusoid smooths over that shape,
  ! but matches the same single-harmonic simplification already used for solar
  ! declination. Sharing these site descriptors with FatesTestLightEnvMod (which
  ! uses FatesTestSiteMod's latitude_deg) means a "different site" experiment
  ! only requires swapping FatesTestSiteMod's values, not this module's.
  !
  ! T_growth (10-day running mean) and T_home (long-term running mean) - the
  ! acclimation temperatures consumed by
  ! photosynth_acclim_model_kumarathunge_etal_2019 - are genuine running means
  ! accumulated day by day as the simulation progresses (UpdateDailyMeans, called
  ! once per day after the sub-daily loop), not copies of the instantaneous
  ! tempk: production FATES only ever consumes these as boundary conditions (see
  ! LeafBiophysicsMod.F90), with no in-repo host-model code computing them, so
  ! there is no existing convention to reuse here.
  !
  ! can_press/can_co2_ppress/can_o2_ppress/gb/dayl_factor/btran are this
  ! module's consolidated reference-atmosphere defaults - the single source
  ! every standalone driver (single_cohort, leaf_level_photo) shares, rather
  ! than each independently hardcoding its own copy (which had drifted: this
  ! module previously baked in ~420 ppm CO2/20.9% O2 while leaf_level_photo
  ! independently computed 380 ppm/21.0%). See the default_* parameters below
  ! Init - leaf_level_photo's original numbers were kept as the shared
  ! default. Init's optional co2_molfrac argument is the override point for a
  ! driver that deliberately needs to diverge (none currently do).
  !

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : pi_const
  use LeafBiophysicsMod, only : QSat
  use FatesTestSiteMod,  only : t_annual_mean
  use FatesTestSiteMod,  only : t_annual_amp
  use FatesTestSiteMod,  only : t_annual_peak_doy
  use FatesTestSiteMod,  only : t_diurnal_amp
  use FatesTestSiteMod,  only : t_diurnal_peak_hour
  use FatesTestSiteMod,  only : relative_humidity

  implicit none
  private

  integer, parameter :: growth_window_days = 10 ! width of the T_growth running-mean window [days]

  ! ------------------------------------------------------------------------------------
  ! Shared reference-atmosphere defaults - the single source of truth for every
  ! standalone test driver (single_cohort, leaf_level_photo), consolidated here
  ! rather than each driver independently hardcoding its own copy. Values are
  ! leaf_level_photo's original defaults (380 ppm CO2, 21.0% O2 - a literature-
  ! standard ambient reference), not this module's own previous ~420 ppm/20.9%
  ! guess. can_co2_ppress is exposed as a mole fraction, converted to partial
  ! pressure at Init() time from the actual can_press, rather than a
  ! hand-rounded partial-pressure constant that could silently drift out of
  ! sync with can_press if it were ever changed. Init()'s optional co2_molfrac
  ! argument is the only override point any driver currently needs; add
  ! further optional arguments here (following the same pattern) if a driver
  ! ever needs to deliberately diverge on one of the others
  ! ------------------------------------------------------------------------------------
  real(r8), public, parameter :: default_co2_molfrac = 380.0e-6_r8 ! [mol/mol] (380 umol/mol)
  real(r8), public, parameter :: default_o2_molfrac  = 210.0e-3_r8 ! [mol/mol] (210 mmol/mol, 21.0%)
  real(r8), public, parameter :: default_dayl_factor = 1.0_r8      ! [0-1] (no seasonal daylength change assumed)
  real(r8), public, parameter :: default_btran       = 1.0_r8      ! [0-1] (non-limiting water assumed)

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

    ! T_growth/T_home running-mean bookkeeping - see the module header
    real(r8) :: daily_temp_buffer(growth_window_days) ! circular buffer of the last (up to) growth_window_days daily-mean temperatures [K]
    integer  :: buffer_next     ! next slot in daily_temp_buffer to overwrite
    integer  :: buffer_count    ! number of days recorded so far, capped at growth_window_days
    real(r8) :: today_temp_sum  ! running sum of today's sub-daily tempk samples [K]
    integer  :: today_n_samples ! number of sub-daily samples accumulated into today_temp_sum so far
    real(r8) :: home_temp_sum   ! cumulative sum of every daily-mean temperature since this trajectory began [K]
    integer  :: home_n_days     ! number of days accumulated into home_temp_sum

  contains

    procedure, public :: Init
    procedure, public :: SetHour
    procedure, public :: UpdateDailyMeans

  end type environment_type

contains

  ! ==========================================================================

  subroutine Init(this, co2_molfrac)
    !
    ! DESCRIPTION:
    ! Set the prescribed atmospheric and soil boundary conditions, and reset the
    ! T_growth/T_home running-mean bookkeeping for a fresh trajectory. tempk itself
    ! is seeded at the annual-cycle mean as a placeholder; SetHour overwrites it
    ! (and everything derived from it) before the first substep is used.

    ! ARGUMENTS:
    class(environment_type), intent(out) :: this ! environment object
    real(r8), optional, intent(in) :: co2_molfrac ! override for can_co2_ppress's mole fraction [mol/mol]; defaults to default_co2_molfrac if omitted

    ! LOCALS:
    real(r8) :: qs_dummy       ! saturation specific humidity output from QSat (unused here)
    real(r8) :: co2_molfrac_local ! this call's actual CO2 mole fraction [mol/mol]

    co2_molfrac_local = default_co2_molfrac
    if (present(co2_molfrac)) co2_molfrac_local = co2_molfrac

    this%can_press      = 101325.0_r8 ! [Pa] (sea level)
    this%can_co2_ppress = co2_molfrac_local * this%can_press ! [Pa]
    this%can_o2_ppress  = default_o2_molfrac * this%can_press ! [Pa]
    this%gb            = 2.0e6_r8         ! [umol/m2/s] (well-ventilated leaf, ~20 s/m equivalent)
    this%dayl_factor   = default_dayl_factor ! [0-1] (distinct from the diurnal temperature/light cycles themselves)
    this%btran         = default_btran       ! [0-1]
    this%rootfr_ft     = 1.0_r8      ! [0-1]
    this%nlevsoil      = 1           ! matches rootfr_ft holding 100% of roots in one layer

    ! placeholder tempk (and everything derived from it) - overwritten by the
    ! first call to SetHour
    this%tempk = t_annual_mean
    call QSat(this%tempk, this%can_press, qs_dummy, this%veg_esat)
    this%can_vpress = this%veg_esat * relative_humidity
    this%t_soil     = this%tempk

    ! T_growth/T_home cold-start seed: the annual-mean temperature is the best
    ! available prior for a running mean that hasn't accumulated any days yet,
    ! since this synthetic climate is a stationary, repeating annual cycle with
    ! no long-term trend to be biased away from
    this%t_growth = t_annual_mean
    this%t_home   = t_annual_mean

    this%daily_temp_buffer = 0.0_r8
    this%buffer_next     = 0
    this%buffer_count    = 0
    this%today_temp_sum   = 0.0_r8
    this%today_n_samples  = 0
    this%home_temp_sum    = 0.0_r8
    this%home_n_days      = 0

  end subroutine Init

  ! ==========================================================================

  subroutine SetHour(this, day_of_year, hour_of_day)
    !
    ! DESCRIPTION:
    ! Prescribe tempk (and everything derived from it: veg_esat via QSat,
    ! can_vpress at the fixed relative_humidity, and t_soil) for the given day of
    ! year and local solar hour of day, via the annual + diurnal sinusoid fit to
    ! BCI TBOT (see the module header). Accumulates this sample into today's
    ! running sum for UpdateDailyMeans to consume at day's end.

    ! ARGUMENTS:
    class(environment_type), intent(inout) :: this        ! environment object
    integer,                  intent(in)    :: day_of_year ! day of year [1-365]
    real(r8),                 intent(in)    :: hour_of_day ! local solar hour of day [0-24]

    ! LOCALS:
    real(r8) :: annual_term  ! this day's offset from t_annual_mean [K]
    real(r8) :: diurnal_term ! this hour's offset from today's mean [K]
    real(r8) :: qs_dummy     ! saturation specific humidity output from QSat (unused here)

    annual_term = t_annual_amp * cos((2.0_r8*pi_const/365.0_r8) *                &
      (real(day_of_year, r8) - t_annual_peak_doy))
    diurnal_term = t_diurnal_amp * cos((pi_const/12.0_r8) *                      &
      (hour_of_day - t_diurnal_peak_hour))

    this%tempk = t_annual_mean + annual_term + diurnal_term
    call QSat(this%tempk, this%can_press, qs_dummy, this%veg_esat)
    this%can_vpress = this%veg_esat * relative_humidity
    this%t_soil     = this%tempk

    this%today_temp_sum  = this%today_temp_sum + this%tempk
    this%today_n_samples = this%today_n_samples + 1

  end subroutine SetHour

  ! ==========================================================================

  subroutine UpdateDailyMeans(this)
    !
    ! DESCRIPTION:
    ! Once-per-day update of the T_growth/T_home running means, from the mean of
    ! today's sub-daily tempk samples (accumulated by SetHour). T_growth is a
    ! growth_window_days-day boxcar (an expanding window until that many days have
    ! been recorded); T_home is an expanding-window mean over the whole trajectory
    ! so far. Resets today's accumulator for tomorrow.

    ! ARGUMENTS:
    class(environment_type), intent(inout) :: this ! environment object

    ! LOCALS:
    real(r8) :: today_mean ! mean of today's sub-daily tempk samples [K]

    today_mean = this%today_temp_sum / real(this%today_n_samples, r8)

    this%buffer_next = mod(this%buffer_next, growth_window_days) + 1
    this%daily_temp_buffer(this%buffer_next) = today_mean
    this%buffer_count = min(this%buffer_count + 1, growth_window_days)
    this%t_growth = sum(this%daily_temp_buffer(1:this%buffer_count)) /          &
      real(this%buffer_count, r8)

    this%home_temp_sum = this%home_temp_sum + today_mean
    this%home_n_days   = this%home_n_days + 1
    this%t_home = this%home_temp_sum / real(this%home_n_days, r8)

    this%today_temp_sum  = 0.0_r8
    this%today_n_samples = 0

  end subroutine UpdateDailyMeans

end module FatesTestEnvironmentMod
