module FatesTestSiteMod
  !
  ! DESCRIPTION:
  ! Prescribed site characteristics for running synthetic environmental conditions
  !
  ! Defaults below are for BCI, Panama (9.153N, 280.154E):
  ! temperature terms are fit to 2003-2016 hourly TBOT from that BCI DATM set (see
  ! FatesTestEnvironmentMod's header for the fitting method), and
  ! relative_humidity is a simple fixed assumption, not fit to data.
  !
  ! ReadSiteNamelist reads a &fates_test_site namelist group and overwrites whichever of 
  ! these variables it lists, leaving any variable the namelist file omits at its BCI default 
  !

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesGlobals,      only : fates_log

  implicit none
  private

  ! ------------------------------------------------------------------------------------
  ! SITE DATA - defaults are for BCI, Panama; override via ReadSiteNamelist
  ! ------------------------------------------------------------------------------------
  real(r8), public :: latitude_deg        = 9.15_r8     ! site latitude, drives the annual solar cycle [deg N]
  real(r8), public :: t_annual_mean       = 299.1172_r8 ! annual-mean temperature [K]
  real(r8), public :: t_annual_amp        = 0.6527_r8   ! amplitude of the annual temperature cycle [K]
  real(r8), public :: t_annual_peak_doy   = 120.555_r8  ! day of year of the annual temperature peak [1-365]
  real(r8), public :: t_diurnal_amp       = 2.2547_r8   ! amplitude of the diurnal temperature cycle [K]
  real(r8), public :: t_diurnal_peak_hour = 13.450_r8   ! local solar hour of the diurnal temperature peak [0-24]
  real(r8), public :: relative_humidity   = 0.80_r8     ! canopy air RH, held fixed [0-1]

  public :: ReadSiteNamelist

contains

  ! ==========================================================================

  subroutine ReadSiteNamelist(nl_filename)
    !
    ! DESCRIPTION:
    ! Reads a &fates_test_site namelist group from nl_filename, overwriting
    ! whichever of the site variables above it lists. Variables the namelist
    ! file omits keep their BCI default (set at declaration above), so a
    ! namelist only needs to specify what differs from BCI.

    ! ARGUMENTS:
    character(len=*), intent(in) :: nl_filename ! path to a &fates_test_site namelist file

    ! LOCALS:
    integer :: unit_nl ! namelist file unit
    integer :: ios     ! iostat return value

    namelist /fates_test_site/ latitude_deg, t_annual_mean, t_annual_amp,      &
      t_annual_peak_doy, t_diurnal_amp, t_diurnal_peak_hour, relative_humidity

    open(newunit=unit_nl, file=trim(nl_filename), status='old', action='read', &
      iostat=ios)
    if (ios /= 0) then
      write(fates_log(),*) 'ReadSiteNamelist: could not open namelist file: ', &
        trim(nl_filename)
      call abort()
    end if

    read(unit_nl, nml=fates_test_site, iostat=ios)
    if (ios /= 0) then
      write(fates_log(),*) 'ReadSiteNamelist: error reading &fates_test_site ' // &
        'namelist group from: ', trim(nl_filename)
      call abort()
    end if

    close(unit_nl)

  end subroutine ReadSiteNamelist

end module FatesTestSiteMod
