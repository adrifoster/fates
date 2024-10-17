module SFCanadianFireWeatherMod

  use FatesConstantsMod,      only : r8 => fates_r8
  use FatesConstantsMod,      only : nearzero
  use FatesInterfaceTypesMod, only : hlm_current_month
  use FatesGlobals,           only : endrun => fates_endrun
  use FatesGlobals,           only : fates_log
  use SFFireWeatherMod,       only : fire_weather
  use shr_log_mod,            only : errMsg => shr_log_errMsg
 
  implicit none
  private

  type, public, extends(fire_weather) :: canadian_fire_weather
  
    real(r8) :: dmc ! duff moisture code

    contains 

      procedure, public :: Init => init_canadian_fire_weather
      procedure, public :: UpdateIndex => update_cfw_index
      procedure, public :: update_dmc
            
  end type canadian_fire_weather

  contains 

    subroutine init_canadian_fire_weather(this)
      !
      !  DESCRIPTION:
      !  Initializes class attributes
      
      ! ARGUMENTS
      class(canadian_fire_weather), intent(inout) :: this ! canadian fire weather index extended class

      ! initialize values to 0.0 (or other initialization)
      this%fire_weather_index   = 0.0_r8
      this%effective_windspeed  = 0.0_r8
      this%dmc                  = 6.0_r8

    end subroutine init_canadian_fire_weather

    !-------------------------------------------------------------------------------------

    subroutine update_cfw_index(this, temp_C, precip, rh, wind, latitude)
      !
      !  DESCRIPTION:
      !  Updates Canadian Fire Weather Indices
      
      ! ARGUMENTS
      class(canadian_fire_weather), intent(inout) :: this     ! canadian index extended class
      real(r8),                     intent(in)    :: temp_C   ! daily averaged temperature [degrees C]
      real(r8),                     intent(in)    :: precip   ! daily precipitation [mm]
      real(r8),                     intent(in)    :: rh       ! daily relative humidity [%]
      real(r8),                     intent(in)    :: wind     ! daily wind speed [m/min]
      real(r8),                     intent(in)    :: latitude ! latitude of site [degress]
      
      call this%update_dmc(temp_C, rh, precip, hlm_current_month, latitude)

    end subroutine update_cfw_index

    !-------------------------------------------------------------------------------------
    
    subroutine update_dmc(this, temp_C, rh, precip, month, latitude)
        !
        !  DESCRIPTION:
        !  Calculates the Canadian Forest Fire Rating System (CFRS) duff moisture code (dmc)
        !
        !  Adapted from Wang 2015 NRC Information Report NOR-X-424
        !
        ! ARGUMENTS:
        class(canadian_fire_weather), intent(inout) :: this      ! canadian index extended class
        real(r8),                     intent(in)    :: temp_C    ! daily average air temperature [degC]
        real(r8),                     intent(in)    :: rh        ! daily average relative humidity (%)
        real(r8),                     intent(in)    :: precip    ! precipitation over last 24 hrs [mm]
        real(r8),                     intent(in)    :: latitude  ! site latitude [degrees]
        integer,                      intent(in)    :: month     ! month of simulation
  
        ! LOCALS
        real(r8), dimension(12) :: day_length_adj    ! day length adjustment
        real(r8)                :: rw                ! net rainfall (mm)
        real(r8)                :: b                 ! temporary variable to calculate moisture content after rain
        real(r8)                :: temp_constrained  ! corrected temperature (degC)
        real(r8)                :: wmi               ! previous day's duff moisture (%)
        real(r8)                :: wmr               ! moisture content after rain (%)
        real(r8)                :: pr                ! corrected rainfall (mm)
        real(r8)                :: rk                ! log drying rate
        real(r8)                :: previous_dmc      ! previous day's dmc

        ! CONSTANTS:

        ! Day length adjustments
        ! For latitude near equator (-10, 10 degrees), use a factor of 9.0 for all months
        
        ! latititude >= 30 N
        real(r8), dimension(12), parameter :: ELL01 = [6.5_r8, 7.5_r8, 9.0_r8, 12.8_r8,  &
          13.9_r8, 13.9_r8, 12.4_r8, 10.9_r8, 9.4_r8, 8.0_r8, 7.0_r8, 6.0_r8]
          
        ! 30 > latitude >= 10
        real(r8), dimension(12), parameter :: ELL02 = [7.9_r8, 8.4_r8, 8.9_r8, 9.5_r8,   &
          9.9_r8, 10.2_r8, 10.1_r8, 9.7_r8, 9.1_r8, 8.6_r8, 8.1_r8, 7.8_r8] 
        
        ! -10 > latitude >= -30
        real(r8), dimension(12), parameter :: ELL03 = [10.1_r8, 9.6_r8, 9.1_r8, 8.5_r8,  &
          8.1_r8, 7.8_r8, 7.9_r8, 8.3_r8, 8.9_r8, 9.4_r8, 9.9_r8, 10.2_r8]
        
        ! latitude < -30
        real(r8), dimension(12), parameter :: ELL04 = [11.5_r8, 10.5_r8, 9.2_r8, 7.9_r8, 6.8_r8,  &
          6.2_r8, 6.5_r8, 7.4_r8, 8.7_r8, 10.0_r8, 11.2_r8, 11.8_r8]
          
        ! save previous day's value of dmc
        previous_dmc = this%dmc

        ! constrain low end of temperature
        if (temp_C < -1.1_r8) then
            temp_constrained = -1.1_r8
        else
            temp_constrained = temp_C
        end if

        ! determine day length adjustment based on latitude
        if (latitude > 30.0_r8) then
            day_length_adj = ELL01
        else if (latitude <= 30.0_r8 .and. latitude > 10.0_r8) then
            day_length_adj = ELL02
        else if (latitude <= 10.0_r8 .and. latitude > -10.0_r8) then
            day_length_adj = 9.0_r8
        else if (latitude <= -10.0_r8 .and. latitude > -30.0_r8) then
            day_length_adj = ELL03
        else if (latitude <= -30.0_r8) then
            day_length_adj = ELL04
        end if

        ! Log drying rate
        rk = 1.894_r8*(temp_constrained + 1.1_r8)*(100.0_r8 - rh)*(day_length_adj(month)*0.0001_r8)

        ! net rainfall (mm)
        rw = 0.92_r8*precip - 1.27_r8

        ! Alteration to Eq. 12 to calculate more accurately
        wmi = 20.0_r8 + 280.0_r8/exp(0.023_r8*previous_dmc)

        ! Eqs 13a-c
        if (previous_dmc <= 33.0_r8) then
            b = 100.0_r8/(0.5_r8 + 0.3_r8*previous_dmc)
        else if (previous_dmc > 33.0_r8 .and. previous_dmc <= 65.0_r8) then
            b = 14.0_r8 - 1.3_r8*log(previous_dmc)
        else
            b = 6.2_r8*log(previous_dmc) - 17.2_r8
        end if

        ! Moisture content after rain
        wmr = wmi + 1000.0_r8*rw/(48.77_r8 + b*rw)

        ! Constrain p
        if (precip <= 1.5_r8) then
            pr = previous_dmc
        else
            pr = 43.43_r8*(5.6348_r8 - log(wmr - 20.0_r8))
        end if

        if (pr < 0.0_r8) pr = 0.0_r8

        ! Calculate dmc
        this%dmc = pr + rk
        if (this%dmc < nearzero) this%dmc = 0.0_r8

    end subroutine update_dmc
    
    !-------------------------------------------------------------------------------------

end module SFCanadianFireWeatherMod