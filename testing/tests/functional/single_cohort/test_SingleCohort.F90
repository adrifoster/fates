program FatesSingleCohort  
  !
  ! DESCRIPTION:
  ! A test of a single cohort across fractional light levels
  ! For each light level, a single cohort is created at recruitment size and simulated
  ! for a specified number of years, including radiation, photosynthesis, and daily
  ! allocation, phenology, and growth. Each light level is an independent trajectory. 
  !
  ! Only cold-deciduous phenology is available currently for non-evergreen trees (i.e.,
  ! no drought deciduous)
  !
  ! Only cabon starvation mortality is calculated and decreases cohort%n from the initial 
  ! 1.0 value

  use FatesConstantsMod, only : r8 => fates_r8

  implicit none

  ! LOCALS:
  character(len=:), allocatable :: param_file   ! input parameter file

  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'my_output_file.nc' ! output file
  
  print *, "Hello SingleCohort"
  
end program FatesSingleCohort

! ----------------------------------------------------------------------------------------
