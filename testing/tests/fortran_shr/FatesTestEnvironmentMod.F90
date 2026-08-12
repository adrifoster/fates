module FatesTestEnvironmentMod
  !
  ! DESCRIPTION:
  ! Prescribed atmospheric and soil boundary conditions for functional tests drivers.
  ! 
  ! Non-default temperature can additionally follow a real annual + diurnal cycle, 
  ! parameterized by FatesTestSiteMod namelist (or default) variables
  !
  ! T_growth (10-day running mean) and T_home (long-term running mean) (the
  ! acclimation temperatures used by photosynth_acclim_model_kumarathunge_etal_2019) are 
  ! running means accumulated day by day as the simulation progresses
  !

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : pi_const
  use FatesConstantsMod, only : t_water_freeze_k_1atm
  use FatesConstantsMod, only : rad_per_deg
  use FatesConstantsMod, only : wm2_to_umolm2s
  use FatesGlobals,      only : fates_log
  use LeafBiophysicsMod, only : QSat
  use LeafBiophysicsMod, only : GetConstrainedVPress
  use FatesTestSiteMod,  only : t_annual_mean
  use FatesTestSiteMod,  only : t_annual_amp
  use FatesTestSiteMod,  only : t_annual_peak_doy
  use FatesTestSiteMod,  only : t_diurnal_amp
  use FatesTestSiteMod,  only : t_diurnal_peak_hour
  use FatesTestSiteMod,  only : t_diurnal2_amp
  use FatesTestSiteMod,  only : t_diurnal2_peak_hour
  use FatesTestSiteMod,  only : latitude_deg
  
  implicit none
  private

  ! CONSTANTS:
  integer,  parameter :: growth_window_days = 10                       ! width of the T_growth running-mean window [days]
  real(r8), parameter :: sea_level_press    = 101325.0_r8              ! [Pa]
  real(r8), parameter :: gb_well_ventilated = 2.0e6_r8                 ! [umol/m2/s] (~20 s/m equivalent)
  real(r8), parameter :: max_declin_deg     = 23.45_r8                 ! Earth's obliquity, used as the declination amplitude (Cooper 1969) [deg]
  real(r8), parameter :: ref_par_full_sun   = 2000.0_r8/wm2_to_umolm2s ! reference full-sun incident PAR at cosz=1 [W/m2] (~2000 umol/m2/s)

  ! shared reference-atmosphere defaults
  real(r8), public, parameter :: default_vpd         = 1000.0_r8   ! leaf-to-air VPD [Pa]
  real(r8), public, parameter :: default_co2_molfrac = 380.0e-6_r8 ! [mol/mol] (380 umol/mol)
  real(r8), public, parameter :: default_o2_molfrac  = 210.0e-3_r8 ! [mol/mol] (210 mmol/mol)
  real(r8), public, parameter :: default_dayl_factor = 1.0_r8      ! [0-1] (no seasonal daylength change assumed)
  real(r8), public, parameter :: default_btran       = 1.0_r8      ! [0-1] (non-limiting water assumed)
  real(r8), public, parameter :: default_veg_tempk   = 25.0_r8 + t_water_freeze_k_1atm ! [K]
  real(r8), public, parameter :: default_par         = 1500.0_r8  ! [umol/m2/s]
  real(r8), public, parameter :: default_nscaler     = 1.0_r8     ! [0-1]
  real(r8), public, parameter :: default_rootfr      = 1.0_r8     ! [0-1]
  real(r8), public, parameter :: direct_frac         = 0.85_r8    ! fraction of incident PAR that is direct beam (typical clear sky)
  real(r8), public, parameter :: diffuse_frac        = 1.0_r8 - direct_frac ! fraction of incident PAR that is diffuse
  integer,  public, parameter :: default_nlevsoil    = 1

  public :: BtranFromSMP
  public :: SoilMatricPotential
  public :: CanopyVaporPressure

  type, public :: environment_type

    real(r8) :: tempk          ! vegetation/leaf temperature [K]
    real(r8) :: can_press      ! air pressure at the leaf surface [Pa]
    real(r8) :: can_co2_ppress ! CO2 partial pressure at the leaf surface [Pa]
    real(r8) :: can_o2_ppress  ! O2 partial pressure at the leaf surface [Pa]
    real(r8) :: veg_esat       ! saturation vapor pressure at tempk [Pa]
    real(r8) :: can_vpress     ! vapor pressure of the canopy air [Pa]
    real(r8) :: gb             ! leaf boundary layer conductance [umol/m2/s]
    real(r8) :: btran          ! soil moisture stress factor [0-1]
    real(r8) :: dayl_factor    ! day-length photosynthetic-capacity acclimation factor [0-1]
    real(r8) :: t_growth       ! 10-day running-mean growth temperature [K]
    real(r8) :: t_home         ! long-term running-mean home temperature [K]
    real(r8) :: rootfr_ft      ! this pft's root fraction in the (single, implicit) soil layer [0-1]
    real(r8) :: t_soil         ! soil temperature [K]
    integer  :: nlevsoil       ! number of (implicit) soil layers

    ! T_growth/T_home running-mean bookkeeping
    real(r8) :: daily_temp_buffer(growth_window_days) ! circular buffer of the last (up to) growth_window_days daily-mean temperatures [K]
    integer  :: buffer_next     ! next slot in daily_temp_buffer to overwrite
    integer  :: buffer_count    ! number of days recorded so far, capped at growth_window_days
    real(r8) :: today_temp_sum  ! running sum of today's sub-daily tempk samples [K]
    integer  :: today_n_samples ! number of sub-daily samples accumulated into today_temp_sum so far
    real(r8) :: home_temp_sum   ! cumulative sum of every daily-mean temperature since this trajectory began [K]
    integer  :: home_n_days     ! number of days accumulated into home_temp_sum

    ! Daily forcing diagnostics - what the run was actually driven with, so that a
    real(r8) :: daily_temp        ! mean of today's sub-daily tempk samples [K]
    real(r8) :: daily_veg_esat    ! mean of today's sub-daily veg_esat samples [Pa]
    real(r8) :: daily_can_vpress  ! mean of today's sub-daily can_vpress samples [Pa]
    real(r8) :: midday_temp       ! tempk at today's substep nearest solar noon [K]
    real(r8) :: midday_veg_esat   ! veg_esat at that substep [Pa]
    real(r8) :: midday_can_vpress ! can_vpress at that substep [Pa]
    real(r8) :: solar_decl        ! solar declination angle [radians]
    real(r8) :: coszen            ! solar zenith angle [radians]
    integer  :: n_vpress_constrained ! substeps today at which GetConstrainedVPress altered the requested vapor pressure [-]

    real(r8) :: today_esat_sum    ! running sum of today's sub-daily veg_esat samples [Pa]
    real(r8) :: today_vpress_sum  ! running sum of today's sub-daily can_vpress samples [Pa]
    real(r8) :: midday_offset     ! |hour_of_day - 12| of the closest-to-noon sample seen so far today [h]
    integer  :: today_n_constrained ! running count of today's constrained substeps [-]

  contains

    procedure, public :: Init
    procedure, public :: SetHour
    procedure, public :: UpdateDailyMeans
    procedure, public :: CalculatePar

  end type environment_type

contains

  ! ==========================================================================
  
    subroutine Init(this, co2_molfrac, tempk)
    !
    ! DESCRIPTION:
    ! Set the prescribed atmospheric and soil boundary conditions, and reset the
    ! T_growth/T_home running-mean bookkeeping for a fresh trajectory.

    ! ARGUMENTS:
    class(environment_type), intent(out) :: this        ! environment object
    real(r8), optional,      intent(in)  :: co2_molfrac ! override for can_co2_ppress's mole fraction [mol/mol]
    real(r8), optional,      intent(in)  :: tempk       ! override for t_annual_mean [K]

    ! LOCALS:
    real(r8) :: qs_dummy          ! saturation specific humidity output from QSat (unused here)
    real(r8) :: co2_molfrac_local ! this call's actual CO2 mole fraction [mol/mol]
    real(r8) :: tempk_local       ! this call's actual temperature [K]

    ! optinally set the co2_molfrac if supplied
    co2_molfrac_local = default_co2_molfrac
    if (present(co2_molfrac)) co2_molfrac_local = co2_molfrac

    ! initialize defaults
    this%can_press = sea_level_press ! [Pa] 
    this%can_co2_ppress = co2_molfrac_local * this%can_press ! [Pa]
    this%can_o2_ppress = default_o2_molfrac * this%can_press ! [Pa]
    this%gb = gb_well_ventilated ! [umol/m2/s]
    this%dayl_factor = default_dayl_factor ! [0-1]
    this%btran = default_btran ! [0-1]
    this%rootfr_ft = default_rootfr ! [0-1]
    this%nlevsoil = default_nlevsoil
    this%solar_decl = SolarDeclination(1)
    this%coszen = CosSolarZenith(this%solar_decl, 0.0_r8, latitude_deg)

    ! placeholder tempk (and everything derived from it)
    tempk_local = t_annual_mean
    if (present(tempk)) tempk_local = tempk
    this%tempk = tempk_local
    call QSat(this%tempk, this%can_press, qs_dummy, this%veg_esat)
    this%can_vpress = CanopyVaporPressure(this%veg_esat)
    this%t_soil = this%tempk

    ! T_growth/T_home cold-start initialization
    this%t_growth = tempk_local
    this%t_home = tempk_local

    this%daily_temp_buffer = 0.0_r8
    this%buffer_next = 0
    this%buffer_count = 0
    this%today_temp_sum = 0.0_r8
    this%today_n_samples = 0
    this%home_temp_sum = 0.0_r8
    this%home_n_days = 0

    ! daily forcing diagnostics
    ! these are overwritten by the first UpdateDailyMeans
    this%daily_temp = this%tempk
    this%daily_veg_esat = this%veg_esat
    this%daily_can_vpress = this%can_vpress
    this%midday_temp = this%tempk
    this%midday_veg_esat = this%veg_esat
    this%midday_can_vpress = this%can_vpress
    this%n_vpress_constrained = 0

    this%today_esat_sum = 0.0_r8
    this%today_vpress_sum = 0.0_r8
    this%today_n_constrained = 0
    this%midday_offset = huge(1.0_r8)

  end subroutine Init

  ! ==========================================================================
  
  function CanopyVaporPressure(veg_esat, vpress_unconstrained, vpd) result(can_vpress)
    !
    ! DESCRIPTION:
    ! Prescribed canopy-air vapor pressure
    !
    ! This value is constrained through GetConstrainedVPress
    !

    ! ARGUMENTS:
    real(r8), intent(in)            :: veg_esat             ! saturation vapor pressure at the current tempk [Pa]
    real(r8), intent(out), optional :: vpress_unconstrained ! uncontrained vapor pressure
    real(r8), intent(in),  optional :: vpd                  ! leaf-to-air VPD to impose [Pa]; default_vpd if absent

    ! RESULT:
    real(r8) :: can_vpress ! canopy air vapor pressure [Pa]

    ! LOCALS:
    real(r8) :: vpd_local ! this call's imposed VPD [Pa]

    vpd_local = default_vpd
    if (present(vpd)) vpd_local = vpd

    can_vpress = veg_esat - vpd_local
    if (present(vpress_unconstrained)) vpress_unconstrained = can_vpress
    can_vpress = GetConstrainedVPress(can_vpress, veg_esat)

  end function CanopyVaporPressure
  
  ! ==========================================================================

  pure function BtranFromSMP(smp, smpsc, smpso) result(btran)
    !
    ! DESCRIPTION:
    ! Soil moisture stress factor (btran) from soil matric potential, for a
    ! driver with no soil column.
    !
    ! Production subroutine EDBtranMod.F90::btran_ed builds btran as a rootfrac-weighted
    ! sum of a per-layer term over the whole soil column:
    !
    !   smp_node = max(smpsc, smp_sl(j))
    !   rresis   = min( (eff_porosity_sl(j)/watsat_sl(j)) *
    !                   (smp_node - smpsc)/(smpso - smpsc), 1 )
    !   btran    = SUM_j rootfrac(j)*rresis(j)
    !
    ! evaluated only over layers holding liquid water.
    !
    ! With a single unfrozen layer at root fraction 1, eff_porosity equals watsat 
    ! (effective porosity is saturated water content minus the ice volume, and we assume
    ! no ice here), so that factor = 1.0, so we can remove it here
    !

    ! ARGUMENTS:
    real(r8), intent(in) :: smp   ! soil matric potential [mm], negative
    real(r8), intent(in) :: smpsc ! soil matric potential at full stomatal closure [mm], negative
    real(r8), intent(in) :: smpso ! soil matric potential at full stomatal opening [mm], negative

    ! RESULT:
    real(r8) :: btran ! soil moisture stress factor [0-1]

    ! LOCALS:
    real(r8) :: smp_node ! smp clamped at the full-closure threshold [mm]

    smp_node = max(smpsc, smp)
    btran = min((smp_node - smpsc)/(smpso - smpsc), 1.0_r8)

  end function BtranFromSMP

  ! ==========================================================================

  pure function SoilMatricPotential(soilfrac, smpsc) result(smp)
    !
    ! DESCRIPTION:
    ! Soil matric potential at a given soil water content, expressed as a
    ! fraction of saturation.
    !
    ! This function interpolates linearly between smpsc at zero water
    ! content and that at saturation

    ! ARGUMENTS:
    real(r8), intent(in) :: soilfrac ! soil water content, fraction of saturation [0-1]
    real(r8), intent(in) :: smpsc    ! soil matric potential at full stomatal closure [mm], negative

    ! RESULT:
    real(r8) :: smp ! soil matric potential [mm], negative

    smp = smpsc*(1.0_r8 - soilfrac)

  end function SoilMatricPotential

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
  
  function CosSolarZenith(declin, hour_of_day, latitude_deg) result(coszen)
    !
    ! DESCRIPTION:
    ! Cosine of the solar zenith angle from the prescribed site latitude, today's
    ! solar declination, and hour of day, via the standard hour-angle formula
    ! (equivalent to CIME's shr_orb_cosz, share/src/shr_orb_mod.F90, re-expressed
    ! in a local-solar-time hour angle measured from solar noon rather than a
    ! Julian-day fraction + longitude). May be negative (sun below the horizon);
    ! callers must not pass a non-positive result to ZenithPrep.

    ! ARGUMENTS:
    real(r8), intent(in) :: declin       ! solar declination angle [radians]
    real(r8), intent(in) :: hour_of_day  ! hour of day [0-24]
    real(r8), intent(in) :: latitude_deg ! latitude [degrees]

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
  
  subroutine CalculatePAR(this, par_beam, par_diff, par_toc, light_frac)
    !
    ! DESCRIPTION:
    ! Calculate PAR (beam and diffuse, based on input )
    
    ! ARGUMENTS:
    class(environment_type), intent(in)            :: this       ! environment object
    real(r8),                intent(out)           :: par_beam   ! direct PAR
    real(r8),                intent(out)           :: par_diff   ! diffuse PAR
    real(r8),                intent(out)           :: par_toc    ! output total PAR
    real(r8),                intent(in),  optional :: light_frac ! optional light fraction
    
    ! LOCALS:
    real(r8) :: light_frac_local ! local override for light_frac
    
    ! coszen < 0.0 means sun below horizon so no PAR
    light_frac_local = 1.0_r8 
    if (present(light_frac)) light_frac_local = light_frac
    par_toc = max(0.0_r8, light_frac_local * ref_par_full_sun * this%coszen)
    par_beam = direct_frac * par_toc
    par_diff = diffuse_frac * par_toc
  
  end subroutine CalculatePAR
  
  ! ==========================================================================

  subroutine SetHour(this, day_of_year, hour_of_day)
    !
    ! DESCRIPTION:
    ! Prescribe tempk (and everything derived from it: veg_esat via QSat,
    ! can_vpress via CanopyVaporPressure, and t_soil) for the given day of
    ! year and local solar hour of day, via the annual + diurnal harmonic fit to
    ! BCI TBOT (see FatesTestSiteMod). Accumulates this sample into today's
    ! running sum for UpdateDailyMeans to consume at day's end.

    ! ARGUMENTS:
    class(environment_type), intent(inout) :: this        ! environment object
    integer,                 intent(in)    :: day_of_year ! day of year [1-365]
    real(r8),                intent(in)    :: hour_of_day ! local solar hour of day [0-24]

    ! LOCALS:
    real(r8) :: annual_term   ! this day's offset from t_annual_mean [K]
    real(r8) :: diurnal_term  ! this hour's offset from today's mean, fundamental [K]
    real(r8) :: diurnal2_term ! this hour's offset from today's mean, overtone [K]
    real(r8) :: qs_dummy      ! saturation specific humidity output from QSat (unused here)
    real(r8) :: noon_offset   ! how far this substep sits from solar noon [h]
    real(r8) :: vpress_uncon  ! canopy air vapor pressure before constraint [Pa]

    annual_term = t_annual_amp * cos((2.0_r8*pi_const/365.0_r8) *                &
      (real(day_of_year, r8) - t_annual_peak_doy))
    diurnal_term = t_diurnal_amp * cos((pi_const/12.0_r8) *                      &
      (hour_of_day - t_diurnal_peak_hour))
    ! the overtone has a 12 h period, hence pi/6
    diurnal2_term = t_diurnal2_amp * cos((pi_const/6.0_r8) *                     &
      (hour_of_day - t_diurnal2_peak_hour))

    this%tempk = t_annual_mean + annual_term + diurnal_term + diurnal2_term
    call QSat(this%tempk, this%can_press, qs_dummy, this%veg_esat)

    ! VPD is held at default_vpd as tempk varies, so canopy air vapor pressure
    ! follows veg_esat
    this%can_vpress = CanopyVaporPressure(this%veg_esat, vpress_uncon)
    if (vpress_uncon /= this%can_vpress) then
      this%today_n_constrained = this%today_n_constrained + 1
    end if

    this%t_soil = this%tempk
    this%today_temp_sum  = this%today_temp_sum + this%tempk
    this%today_n_samples = this%today_n_samples + 1
    
    ! light environment variables
    this%solar_decl = SolarDeclination(day_of_year)
    this%coszen = CosSolarZenith(this%solar_decl, hour_of_day, latitude_deg)

    ! daily forcing diagnostics (see the type definition). The vapor pressure
    ! sums are accumulated over the same samples as the temperature sum, so
    ! every daily mean this module reports is on one consistent sampling.
    this%today_esat_sum   = this%today_esat_sum + this%veg_esat
    this%today_vpress_sum = this%today_vpress_sum + this%can_vpress

    ! keep the sample closest to solar noon. Ties go to the later substep (<=
    ! rather than <), which for an even number of substeps per day picks the
    ! first sample after noon - matching the noon_substep convention the drivers
    ! use for their own once-per-year snapshots
    noon_offset = abs(hour_of_day - 12.0_r8)
    if (noon_offset <= this%midday_offset) then
      this%midday_offset     = noon_offset
      this%midday_temp       = this%tempk
      this%midday_veg_esat   = this%veg_esat
      this%midday_can_vpress = this%can_vpress
    end if

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

    ! today's forcing diagnostics, finalized so that a driver can record it
    this%daily_temp = today_mean
    this%daily_veg_esat = this%today_esat_sum / real(this%today_n_samples, r8)
    this%daily_can_vpress = this%today_vpress_sum / real(this%today_n_samples, r8)
    this%n_vpress_constrained = this%today_n_constrained
    
    ! now reset
    this%today_temp_sum  = 0.0_r8
    this%today_n_samples = 0
    this%today_esat_sum = 0.0_r8
    this%today_vpress_sum = 0.0_r8
    this%today_n_constrained = 0
    this%midday_offset = huge(1.0_r8)

  end subroutine UpdateDailyMeans
  
  ! ==========================================================================

end module FatesTestEnvironmentMod
