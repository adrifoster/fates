program FatesTestROS
  
  use FatesConstantsMod,           only : r8 => fates_r8 
  use FatesTestFireMod,            only : SetUpFuel, WriteROSData
  use FatesArgumentUtils,          only : command_line_arg
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use SyntheticFuelModels,         only : fuel_models_array_class
  use FatesFuelMod,                only : fuel_type
  use FatesFuelClassesMod,         only : num_fuel_classes
  use SFParamsMod,                 only : SF_val_SAV
  use SFParamsMod,                 only : SF_val_FBD, SF_val_part_dens
  use SFParamsMod,                 only : SF_val_miner_total
  use SFMainMod,                   only : HeatOfPreignition, EffectiveHeatingNumber
  use SFMainMod,                   only : PhiWind, PropagatingFlux
  use SFMainMod,                   only : RateOfSpread, ReactionIntensity
  
  implicit none
  
  ! LOCALS:
  type(fates_unit_test_param_reader)             :: param_reader         ! param reader instance
  type(fuel_models_array_class)                  :: fuel_models_array    ! array of fuel models
  type(fuel_type),                   allocatable :: fuel(:)              ! fuel objects                          
  character(len=:),                  allocatable :: param_file           ! input parameter file
  real(r8),                          allocatable :: wind_speed(:)        ! midflame wind speed [m/min]
  real(r8)                                       :: wind_adj             ! adjusted windspeed [m/min]
  real(r8),                          allocatable :: fuel_moisture(:,:)   ! fuel moisture [m3/m3]
  real(r8),                          allocatable :: beta(:)              ! packing ratio [dimensionless]
  real(r8),                          allocatable :: beta_op(:)           ! optimum packing ratio [dimensionless]
  real(r8),                          allocatable :: eps(:)               ! effective heating number [dimensionless]
  real(r8),                          allocatable :: prop_flux(:)         ! propagating flux ratio [dimensionless]
  real(r8),                          allocatable :: q_ig(:)              ! heat of preignition [kJ/kg]
  real(r8),                          allocatable :: i_r(:,:)             ! reaction intensity [kJ/m2/min]
  real(r8),                          allocatable :: phi_wind(:,:)        ! wind factor [dimensionless]
  real(r8),                          allocatable :: ros(:,:,:)           ! rate of spread [m/min]
  real(r8)                                       :: beta_ratio           ! ratio of packing ratio to optimum packing ratio [dimensionless]
  real(r8)                                       :: non_mineral_loading  ! non mineral loading [kgC/m2]
  integer                                        :: num_fuel_models      ! number of fuel models to simulate
  character(len=100),                allocatable :: fuel_names(:)        ! names of fuel models
  character(len=2),                  allocatable :: carriers(:)          ! carriers of fuel models
  integer                                        :: num_wind             ! array size
  integer                                        :: num_moist            ! number of moisture classes
  integer                                        :: i, f, m              ! looping indices
  
  ! CONSTANTS:
  character(len=*), parameter               :: out_file = 'ros_out.nc' ! output file 
  integer,          parameter, dimension(5) :: fuel_models = (/101, 102, 103, 104, 105/)   ! fuel models to test
  
  ! fuel moisture values to test - from Scott & Bergen 2005
  real(r8),         parameter, dimension(4) :: fuel_moisture_1hr =                       &
                                              (/3.0_r8, 6.0_r8, 9.0_r8, 12.0_r8/) ! dead fuel moisture values for 1-hr fuels [%]
  real(r8),         parameter, dimension(4) :: fuel_moisture_10hr =                      &
                                              (/4.0_r8, 7.0_r8, 10.0_r8, 13.0_r8/) ! dead fuel moisture values for 10-hr fuels [%]
  real(r8),         parameter, dimension(4) :: fuel_moisture_100hr =                     &
                                              (/5.0_r8, 8.0_r8, 11.0_r8, 14.0_r8/) ! dead fuel moisture values for 100-hr fuels [%]
  real(r8),         parameter, dimension(4) :: fuel_moisture_live_herb =                 &
                                              (/30.0_r8, 60.0_r8, 90.0_r8, 120.0_r8/) ! live herbaceous fuel moisture values [%]
  real(r8),         parameter, dimension(4) :: fuel_moisture_live_woody =                &
                                              (/60.0_r8, 90.0_r8, 120.0_r8, 150.0_r8/) ! live woody fuel moisture values [%]
  real(r8), parameter :: min_wind = 0.0_r8  ! minimum wind speed to test [m/s]
  real(r8), parameter :: max_wind = 540.0_r8 ! maximum wind speed to test [m/s]
  real(r8), parameter :: wind_inc = 1.0_r8  ! wind increment to use [m/s]
                                           
  ! read in parameter file name and DATM file from command line
  param_file = command_line_arg(1)
  
  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()
  
  num_wind = int((max_wind - min_wind)/wind_inc + 1)
  num_moist = size(fuel_moisture_1hr)
  num_fuel_models = size(fuel_models)
  
  allocate(wind_speed(num_wind))
  allocate(q_ig(num_moist))
  allocate(fuel(num_fuel_models))
  allocate(fuel_names(num_fuel_models))
  allocate(carriers(num_fuel_models))
  allocate(beta(num_fuel_models))
  allocate(beta_op(num_fuel_models))
  allocate(eps(num_fuel_models))
  allocate(prop_flux(num_fuel_models))
  allocate(fuel_moisture(num_fuel_models, num_moist))
  allocate(i_r(num_fuel_models, num_moist))
  allocate(phi_wind(num_fuel_models, num_wind))
  allocate(ros(num_fuel_models, num_moist, num_wind))
  
  ! set up fuel objects and calculate loading
  call fuel_models_array%GetFuelModels()
    
  ! calculate ROS -------
  do f = 1, num_fuel_models
    
    ! uses data from fuel_models to initialize fuel
    call SetUpFuel(fuel(f), fuel_models_array, fuel_models(f), fuel_names(f), carriers(f))
    
    ! sum up fuel and calculate loading
    call fuel(f)%SumLoading()
    call fuel(f)%CalculateFractionalLoading()
     
    beta(f) = fuel(f)%bulk_density_notrunks/SF_val_part_dens
    print *, SF_val_part_dens
    beta_op(f) = 0.200395_r8*(fuel(f)%SAV_notrunks**(-0.8189_r8))
    beta_ratio = beta(f)/beta_op(f)
    
    eps(f) = EffectiveHeatingNumber(fuel(f)%SAV_notrunks)
    prop_flux(f) = PropagatingFlux(fuel(f)%SAV_notrunks, beta(f))

    non_mineral_loading = fuel(f)%non_trunk_loading*(1.0_r8 - SF_val_miner_total)
    
    do m = 1, num_moist

      call fuel_models_array%fuel_models(f)%CalculateMoisture(fuel_moisture_1hr(m),      &
        fuel_moisture_10hr(m), fuel_moisture_100hr(m), fuel_moisture_live_herb(m),       &
        fuel_moisture_live_woody(m), fuel_moisture(f,m)) 
      
      q_ig(m) = HeatOfPreignition(fuel_moisture(f,m))
      
      i_r(f,m) = ReactionIntensity(non_mineral_loading, fuel(f)%SAV_notrunks,            &
        beta_ratio, fuel_moisture(f,m), fuel(f)%MEF_notrunks)
        
      do i = 1, num_wind
        
        wind_speed(i) = min_wind + wind_inc*(i-1)
        wind_adj = wind_speed(i)*(fuel_models_array%fuel_models(f)%wind_adj_factor)
        
        phi_wind(f,i) = PhiWind(wind_adj, beta_ratio, fuel(f)%SAV_notrunks)
        
        ros(f,m,i) = RateOfSpread(fuel(f)%bulk_density_notrunks, eps(f), q_ig(m),        &
          i_r(f,m), prop_flux(f), phi_wind(f,i))
      
      end do
    end do
  end do
  
  call WriteROSData(out_file, num_wind, num_moist, num_fuel_models, wind_speed, beta,    &
      beta_op, eps, prop_flux, q_ig, fuel_moisture, i_r, phi_wind, ros, fuel_models)
  
  if (allocated(wind_speed)) deallocate(wind_speed)
  if (allocated(fuel)) deallocate(fuel)
  if (allocated(fuel_names)) deallocate(fuel_names)
  if (allocated(carriers)) deallocate(carriers)
  if (allocated(beta)) deallocate(beta)
  if (allocated(beta_op)) deallocate(beta_op)
  if (allocated(eps)) deallocate(eps)
  if (allocated(prop_flux)) deallocate(prop_flux)
  if (allocated(q_ig)) deallocate(q_ig)
  if (allocated(fuel_moisture)) deallocate(fuel_moisture)
  if (allocated(i_r)) deallocate(i_r)
  if (allocated(phi_wind)) deallocate(phi_wind)
  if (allocated(ros)) deallocate(ros)
  
end program FatesTestROS
