module FatesTestSiteMod
  !
  ! DESCRIPTION:
  ! Prescribed site descriptors, shared between FatesTestLightEnvMod (latitude,
  ! for solar declination/coszen) and FatesTestEnvironmentMod (everything else,
  ! for the annual/diurnal temperature cycle and canopy-air humidity). Kept in
  ! one place, as an explicit bundle, so a "different site" experiment (e.g. a
  ! false-latitude sensitivity test) means swapping in a different self-consistent
  ! set of these values together, rather than deriving some of them from others
  ! via an invented latitude-sensitivity formula.
  !
  ! Values below are for BCI, Panama (9.153N, 280.154E - the coordinates baked
  ! into ~/Documents/02_Projects/06_FATES/BCI_datm's forcing files): the
  ! temperature terms are fit to 2003-2016 hourly TBOT from that DATM set (see
  ! FatesTestEnvironmentMod's header for the fitting method), and
  ! relative_humidity is a simple fixed assumption, not fit to data.
  !

  use FatesConstantsMod, only : r8 => fates_r8

  implicit none
  private

  ! ------------------------------------------------------------------------------------
  ! SITE DATA - BCI, Panama
  ! ------------------------------------------------------------------------------------
  real(r8), public, parameter :: latitude_deg        = 9.15_r8    ! site latitude, drives the annual solar cycle [deg N]
  real(r8), public, parameter :: t_annual_mean       = 299.1172_r8 ! annual-mean temperature [K]
  real(r8), public, parameter :: t_annual_amp        = 0.6527_r8   ! amplitude of the annual temperature cycle [K]
  real(r8), public, parameter :: t_annual_peak_doy   = 120.555_r8  ! day of year of the annual temperature peak [1-365]
  real(r8), public, parameter :: t_diurnal_amp       = 2.2547_r8   ! amplitude of the diurnal temperature cycle [K]
  real(r8), public, parameter :: t_diurnal_peak_hour = 13.450_r8   ! local solar hour of the diurnal temperature peak [0-24]
  real(r8), public, parameter :: relative_humidity   = 0.70_r8     ! canopy air RH, held fixed [0-1]

end module FatesTestSiteMod
