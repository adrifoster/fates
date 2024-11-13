module SyntheticFuelModels
  
  use FatesConstantsMod, only : r8 => fates_r8
  use SFParamsMod,       only : SF_val_part_dens

  implicit none
  private
  
  ! Fuel model numbers come from Scott and Burgen (2005) RMRS-GTR-153
  integer, parameter, public, dimension(52) :: all_fuel_models = (/1, 2, 101, 102, 104,  &
                                                      107, 121, 122, 3, 103, 105, 106,   &
                                                      108, 109, 123, 124, 4, 5, 6, 141,  &
                                                      142, 145, 147, 161, 164, 10, 7,    &
                                                      143, 144, 146, 148, 149, 162,      &
                                                      163, 8, 9, 181, 182, 183, 184,     &
                                                      185, 186, 187, 188, 189, 11, 12,   &
                                                      13, 201, 202, 203, 204/)
  
  integer,  parameter :: chunk_size = 10
  real(r8), parameter :: ustons_to_kg = 907.185_r8
  real(r8), parameter :: acres_to_m2 = 4046.86_r8
  real(r8), parameter :: ustons_acre_to_kgC_m2 = ustons_to_kg/acres_to_m2 !*0.45_r8
  real(r8), parameter :: ft_to_m = 0.3048_r8
  
  ! holds data for fake fuel models that can be used for functional
  ! testing of the FATES fire model
  ! these are taken from the fire behavior fuel models in Scott & Burgan 2005
  type, public :: synthetic_fuel_model
    
    integer            :: fuel_model_index       ! fuel model index
    character(len=2)   :: carrier                ! carrier ('GR', 'GS', etc.)
    character(len=5)   :: fuel_model_code        ! carrier plus fuel model
    character(len=100) :: fuel_model_name        ! long name of fuel model
    real(r8)           :: wind_adj_factor        ! wind adjustment factor
    real(r8)           :: dead_fuel_loading(3)   ! loading for 1 hr, 10 hr, and 100 hr fuels [kg/m2]
    real(r8)           :: live_fuel_loading(2)   ! loading for live woody and herbaceous fuels [kg/m2]
    real(r8)           :: fuel_depth             ! fuel bed depth [m]
    real(r8)           :: dead_fuel_sav(3)       ! surface area to volume ratio of dead fuels [/cm]
    real(r8)           :: live_fuel_sav(2)       ! surface area to volume ratio of live fuels [/cm]
    real(r8)           :: dead_fuel_weighting(3) ! weighting factor for dead fuels [/cm]
    real(r8)           :: live_fuel_weighting(2) ! weighting factor for dead fuels [/cm]
    real(r8)           :: moist_extinct          ! dead fuel extinction moisture [m3/m3]
    real(r8)           :: bulk_density           ! bulk density of fuel [kg/m3]
    real(r8)           :: sav                    ! surface area to volume ratio of fuel [/cm]
    real(r8)           :: mean_live_surf_area    ! average surface area of live fuels [m2]
    real(r8)           :: mean_dead_surf_area    ! average surface area of dead fuels [m2]
    real(r8)           :: net_loading_live       ! net fuel loading of live fuels [kg/m2]
    real(r8)           :: net_loading_dead       ! net fuel loading of dead fuels [kg/m2]
    real(r8)           :: q_ig_dead(3)
    real(r8)           :: q_ig_live(2)
    real(r8)           :: mean_weighting_live
    real(r8)           :: mean_weighting_dead

    contains 
    
      procedure :: InitFuelModel
      procedure :: CalculateMoisture
      procedure :: CalculateWeightingFactors
      procedure :: CalculateSAV
      procedure :: CalculateNetFuelLoad
      procedure :: CalculateDeadLiveRatio
      procedure :: CalculateHeatSink
  
  end type synthetic_fuel_model
  
  ! --------------------------------------------------------------------------------------
  
  ! a class to just hold an array of these fuel models
  type, public :: fuel_models_array_class
  
    type(synthetic_fuel_model), allocatable :: fuel_models(:)  ! array of fuel models
    integer                                 :: num_fuel_models ! number of total fuel models
    
    contains 
      
      procedure :: AddFuelModel
      procedure :: GetFuelModels
      procedure :: FuelModelPosition
  
  end type fuel_models_array_class
  
  ! --------------------------------------------------------------------------------------
  
  contains 
    
  subroutine InitFuelModel(this, fuel_model_index, carrier, fuel_model_name,             &
    wind_adj_factor, hr1_loading, hr10_loading, hr100_loading, live_herb_loading,        &
    live_woody_loading, fuel_depth, hr1_sav, live_herb_sav, live_woody_sav, moist_extinct)
    !
    ! DESCRIPTION:
    ! Initializes the fuel model with input characteristics
    ! Also converts units as needed
    !
    ! NOTE THE UNITS ON INPUTS
    !
    
    ! ARGUMENTS:
    class(synthetic_fuel_model), intent(inout) :: this
    integer,                     intent(in)    :: fuel_model_index   ! fuel model index
    character(len=2),            intent(in)    :: carrier            ! main carrier
    character(len=*),            intent(in)    :: fuel_model_name    ! fuel model long name
    real(r8),                    intent(in)    :: wind_adj_factor    ! wind adjustment factor
    real(r8),                    intent(in)    :: hr1_loading        ! loading for 1-hr fuels [US tons/acre]
    real(r8),                    intent(in)    :: hr10_loading       ! loading for 10-hr fuels [US tons/acre]
    real(r8),                    intent(in)    :: hr100_loading      ! loading for 100-hr fuels [US tons/acre]
    real(r8),                    intent(in)    :: live_herb_loading  ! loading for live herbacious fuels [US tons/acre]
    real(r8),                    intent(in)    :: live_woody_loading ! loading for live woody fuels [US tons/acre]
    real(r8),                    intent(in)    :: fuel_depth         ! fuel bed depth [ft]
    real(r8),                    intent(in)    :: hr1_sav            ! surface area to volume ratio of 1 hour fuels [/ft]
    real(r8),                    intent(in)    :: live_herb_sav      ! surface area to volume ratio of live herbacious fuels [/ft]
    real(r8),                    intent(in)    :: live_woody_sav     ! surface area to volume ratio of live woody fuels [/ft]
    real(r8),                    intent(in)    :: moist_extinct      ! dead fuel extinction moisture [%]
    
    ! LOCALS:
    real(r8) :: fine_fuel_loading ! fine fuel loading [kg/m2]
        
    this%fuel_model_index = fuel_model_index
    this%carrier = carrier 
    this%fuel_model_name = fuel_model_name
    this%wind_adj_factor = wind_adj_factor
    this%dead_fuel_loading(1) = hr1_loading*ustons_acre_to_kgC_m2 ! convert to kgC/m2
    this%dead_fuel_loading(2) = hr10_loading*ustons_acre_to_kgC_m2 ! convert to kgC/m2
    this%dead_fuel_loading(3) = hr100_loading*ustons_acre_to_kgC_m2 ! convert to kgC/m2
    this%live_fuel_loading(1) = live_herb_loading*ustons_acre_to_kgC_m2  ! convert to kgC/m2
    this%live_fuel_loading(2) = live_woody_loading*ustons_acre_to_kgC_m2  ! convert to kgC/m2
    this%fuel_depth = fuel_depth*ft_to_m ! convert to m
    this%dead_fuel_sav(1) = hr1_sav/ft_to_m/100.0_r8 ! convert to cm
    this%dead_fuel_sav(2) = 109.0_r8/ft_to_m/100.0_r8 ! convert to cm
    this%dead_fuel_sav(3) = 30.0_r8/ft_to_m/100.0_r8 ! convert to cm
    this%live_fuel_sav(1) = live_herb_sav/ft_to_m/100.0_r8 ! convert to cm
    this%live_fuel_sav(2) = live_woody_sav/ft_to_m/100.0_r8 ! convert to cm
    this%moist_extinct = moist_extinct/100.0_r8 ! convert to [m3/m3]
    
    this%bulk_density = (sum(this%dead_fuel_loading(:)) + sum(this%live_fuel_loading(:)))/this%fuel_depth
    
    call this%CalculateWeightingFactors()
    call this%CalculateSAV()
    call this%CalculateNetFuelLoad(0.0555_r8)
    
  end subroutine InitFuelModel
  
  ! --------------------------------------------------------------------------------------
  
  subroutine CalculateHeatSink(this, bulk_density, heat_sink)
    !
    ! DESCRIPTION:
    ! Calculates the heat sink for the ROS equation as in Rothermel
    !
    
    ! ARGUMENTS:
    class(synthetic_fuel_model), intent(in)  :: this
    real(r8),                    intent(in)  :: bulk_density
    real(r8),                    intent(out) :: heat_sink
    
    ! LOCALS:
    integer :: i
    real(r8) :: heat_sink_live
    real(r8) :: heat_sink_dead 
    
    heat_sink_dead = 0.0_r8 
    do i = 1, 3
      if (this%dead_fuel_loading(i) > 0.0_r8) then 
        heat_sink_dead = heat_sink_dead + this%dead_fuel_weighting(i)*(exp(-4.528_r8/this%dead_fuel_sav(i)))*this%q_ig_dead(i)
      end if 
    end do
  
    heat_sink_live = 0.0_r8 
    do i = 1, 2
      if (this%live_fuel_weighting(i) > 0.0_r8) then 
        heat_sink_live = heat_sink_live + this%live_fuel_weighting(i)*(exp(-4.528_r8/this%live_fuel_sav(i)))*this%q_ig_live(i)
      end if
    end do
    
    heat_sink = bulk_density*(this%mean_weighting_dead*heat_sink_dead + this%mean_weighting_live*heat_sink_live)
    
  end subroutine CalculateHeatSink
  
  ! --------------------------------------------------------------------------------------

  real(r8) function FuelSurfaceArea(load, sav, particle_density)
    !
    ! DESCRIPTION:
    ! Calculates oven-dry fuel load, including noncombustable mineral fraction
    !
    
    ! ARGUMENTS:
    real(r8), intent(in) :: load             ! mean load of size class [kg/m2]
    real(r8), intent(in) :: sav              ! surface are to volume ratio [/cm]
    real(r8), intent(in) :: particle_density ! particle density [kg/m3]   
    
    FuelSurfaceArea = (sav*load)/particle_density
  
  end function FuelSurfaceArea
  
  ! --------------------------------------------------------------------------------------
  
  subroutine CalculateWeightingFactors(this)
    !
    ! DESCRIPTION:
    ! Calculates fuel surface area and weighting factors
    !
    
    ! ARGUMENTS:
    class(synthetic_fuel_model), intent(inout) :: this
    
    ! LOCALS:
    integer  :: i                    ! looping index
    real(r8) :: surface_area_live(2) ! surface area of live fuels [m2]
    real(r8) :: surface_area_dead(3) ! surface area of dead fuels [m2]
    real(r8) :: mean_live_sa         ! mean surface area of live fuels [m2]
    real(r8) :: mean_dead_sa         ! mean surface area of dead fuels [m2]
    
    mean_live_sa = 0.0_r8
    mean_dead_sa = 0.0_r8
    do i = 1, 3
      if (this%dead_fuel_loading(i) > 0.0_r8) then 
        surface_area_dead(i) = FuelSurfaceArea(this%dead_fuel_loading(i),                  &
          this%dead_fuel_sav(i), SF_val_part_dens)
      else 
        surface_area_dead(i) = 0.0_r8
      end if
      mean_dead_sa = mean_dead_sa + surface_area_dead(i)
    end do 
    do i = 1, 2
      if (this%live_fuel_loading(i) > 0.0_r8) then 
        surface_area_live(i) = FuelSurfaceArea(this%live_fuel_loading(i),                  &
          this%live_fuel_sav(i), SF_val_part_dens)
      else
        surface_area_live(i) = 0.0_r8
      end if
      mean_live_sa = mean_live_sa + surface_area_live(i)
    end do
    
    if (mean_dead_sa > 0.0_r8) then 
      this%dead_fuel_weighting(:) = surface_area_dead(:)/mean_dead_sa
    else
      this%dead_fuel_weighting(:) = 0.0_r8 
    end if 
    if (mean_live_sa > 0.0_r8) then 
      this%live_fuel_weighting(:) = surface_area_live(:)/mean_live_sa
    else 
      this%live_fuel_weighting(:) = 0.0_r8 
    end if
    
    this%mean_live_surf_area = mean_live_sa
    this%mean_dead_surf_area = mean_dead_sa
    
    ! calculate weighting factors for live and dead
    this%mean_weighting_live = this%mean_live_surf_area/(this%mean_live_surf_area + this%mean_dead_surf_area)
    this%mean_weighting_dead = this%mean_dead_surf_area/(this%mean_live_surf_area + this%mean_dead_surf_area)

  end subroutine CalculateWeightingFactors
  
  ! --------------------------------------------------------------------------------------
  
  subroutine CalculateSAV(this)
    !
    ! DESCRIPTION:
    ! Calculates fuel surface area and weighting factors
    !
    
    ! ARGUMENTS:
    class(synthetic_fuel_model), intent(inout) :: this
    
    ! LOCALS:
    integer  :: i             ! looping index
    real(r8) :: mean_live_sav ! mean surface area to volume ratio of live fuels [/cm]
    real(r8) :: mean_dead_sav ! mean surface area to volume ratio of dead fuels [/cm]
    real(r8) :: live_wf       ! live weighting factor
    real(r8) :: dead_wf       ! dead weighting factor
    

    
    ! aclculate surface-area-to-volume ratio of live and dead fuels
    mean_live_sav = 0.0_r8
    mean_dead_sav = 0.0_r8
    do i = 1, 3
      if (this%dead_fuel_loading(i) > 0.0_r8) then 
        mean_dead_sav = mean_dead_sav + this%dead_fuel_weighting(i)*this%dead_fuel_sav(i)
      end if
    end do 
    do i = 1, 2
      if (this%live_fuel_loading(i) > 0.0_r8) then 
        mean_live_sav = mean_live_sav + this%live_fuel_weighting(i)*this%live_fuel_sav(i)
      end if
    end do
    
    ! calculate average surface area to volume ratio
    this%sav = mean_dead_sav*this%mean_weighting_dead + mean_live_sav*this%mean_weighting_live
    
  end subroutine CalculateSAV
  
  ! --------------------------------------------------------------------------------------
  
  subroutine CalculateNetFuelLoad(this, mineral_fraction)
    !
    ! DESCRIPTION:
    ! Calculates fuel surface area and weighting factors
    !
    
    ! ARGUMENTS:
    class(synthetic_fuel_model), intent(inout) :: this
    real(r8),                    intent(in)    :: mineral_fraction ! total mineral content [fraction]
    
    ! LOCALS:
    integer  :: i                  ! looping index
    real(r8) :: net_fuel_live(2)   ! net fuel loading of live fuels [kg/m2]
    real(r8) :: net_fuel_dead(3)   ! net fuel loading of dead fuels [kg/m2]
    real(r8) :: total_live_loading ! total loading for live fuels [kg/m2]
    real(r8) :: total_dead_loading ! total loading for dead fuels [kg/m2]
    real(r8) :: weighting_factor
    
    total_dead_loading = 0.0_r8
    total_live_loading = 0.0_r8
    
    do i = 1, 3
      if (this%dead_fuel_loading(i) > 0.0_r8) then 
        net_fuel_dead(i) = this%dead_fuel_loading(i)*(1.0_r8 - mineral_fraction)
        if (this%dead_fuel_sav(i) < 0.525_r8) then 
          weighting_factor = 0.0_r8 
        else 
          weighting_factor = this%dead_fuel_weighting(i)
        end if
        total_dead_loading = total_dead_loading + net_fuel_dead(i)*weighting_factor
      end if
    end do
    
    do i = 1, 2
      if (this%live_fuel_loading(i) > 0.0_r8) then 
        net_fuel_live(i) = this%live_fuel_loading(i)*(1.0_r8 - mineral_fraction)
        if (this%live_fuel_sav(i) < 0.525_r8) then 
          weighting_factor = 0.0_r8 
        else 
          weighting_factor = this%live_fuel_weighting(i)
        end if
        total_live_loading = total_live_loading + net_fuel_live(i)*weighting_factor
      end if
    end do
    
    this%net_loading_live = total_live_loading
    this%net_loading_dead = total_dead_loading
    
  end subroutine CalculateNetFuelLoad
  
  ! --------------------------------------------------------------------------------------

  subroutine CalculateMoisture(this, hr1_moisture, hr10_moisture, hr100_moisture,        & 
    live_herb_moisture, live_woody_moisture, live_moisture, dead_moisture, mef_live)
    !
    ! DESCRIPTION:
    ! Calculates average fuel moisture and dead fuel moisture of extinction for the 
    ! synthetic fuel based on input values
    !
    
    ! ARGUMENTS:
    class(synthetic_fuel_model), intent(inout)  :: this
    real(r8),                    intent(in)  :: hr1_moisture        ! 1-hr fuel moisture [%]
    real(r8),                    intent(in)  :: hr10_moisture       ! 10-hr fuel moisture [%]
    real(r8),                    intent(in)  :: hr100_moisture      ! 100-hr fuel moisture [%]
    real(r8),                    intent(in)  :: live_herb_moisture  ! live herb fuel moisture [%]
    real(r8),                    intent(in)  :: live_woody_moisture ! live woody fuel moisture [%]
    real(r8),                    intent(out) :: live_moisture       ! average live fuel moisture [m3/m3]
    real(r8),                    intent(out) :: dead_moisture       ! average dead fuel moisture [m3/m3]
    real(r8),                    intent(out) :: mef_live            ! average live fuel moisture of extinction [m3/m3]
    
    ! LOCALS:
    integer  :: i 
    real(r8) :: dead_live_ratio
    real(r8) :: fine_dead_moisture
    real(r8) :: fine_dead_fuel_moisture_top
    real(r8) :: fine_dead_fuel_moisture_bottom
    real(r8) :: dead_moisture_all(3)
    real(r8) :: live_moisture_all(2)
    
    call this%CalculateDeadLiveRatio(dead_live_ratio)

    dead_moisture_all(1) = hr1_moisture/100.0_r8
    dead_moisture_all(2) = hr10_moisture/100.0_r8
    dead_moisture_all(3) = hr100_moisture/100.0_r8
    live_moisture_all(1) = live_herb_moisture/100.0_r8
    live_moisture_all(2) = live_woody_moisture/100.0_r8
    
    ! calculate average live and dead moisture
    dead_moisture = 0.0_r8 
    do i = 1, 3
      if (this%dead_fuel_loading(i) > 0.0_r8) then 
        dead_moisture = dead_moisture + this%dead_fuel_weighting(i)*dead_moisture_all(i)
      end if
    end do
    live_moisture = 0.0_r8 
    do i = 1, 2
      if (this%live_fuel_loading(i) > 0.0_r8) then 
        live_moisture = live_moisture + this%live_fuel_weighting(i)*live_moisture_all(i)
      end if
    end do
  
    ! calculate live fuel moisture of extinction
    fine_dead_fuel_moisture_top = 0.0_r8 
    fine_dead_fuel_moisture_bottom = 0.0_r8
    do i = 1, 3
      if (this%dead_fuel_loading(i) > 0.0_r8) then 
        fine_dead_fuel_moisture_top = fine_dead_fuel_moisture_top +                        &
          dead_moisture_all(i)*this%dead_fuel_loading(i)*exp(-4.528/this%dead_fuel_sav(i))
        fine_dead_fuel_moisture_bottom = fine_dead_fuel_moisture_bottom +                  &
          this%dead_fuel_loading(i)*exp(-4.528/this%dead_fuel_sav(i))
      end if
    end do
    
    fine_dead_moisture = fine_dead_fuel_moisture_top/fine_dead_fuel_moisture_bottom
    mef_live = 2.9_r8*dead_live_ratio*(1.0_r8 - fine_dead_moisture/this%moist_extinct) - 0.226_r8
    if (mef_live < this%moist_extinct) mef_live = this%moist_extinct
    
    ! calculate qig
    do i = 1, 3
      if (this%dead_fuel_loading(i) > 0.0_r8) then 
        this%q_ig_dead(i) =  581.0_r8 + 2594.0_r8*dead_moisture_all(i)
      else 
        this%q_ig_dead(i) = 0.0_r8
      end if
    end do
    do i = 1, 2
      if (this%live_fuel_loading(i) > 0.0_r8) then 
        this%q_ig_live(i) =  581.0_r8 + 2594.0_r8*live_moisture_all(i)
      else 
        this%q_ig_live(i) = 0.0_r8 
      end if
    end do
    
  end subroutine CalculateMoisture
  
  ! --------------------------------------------------------------------------------------
  
  subroutine CalculateDeadLiveRatio(this, dead_live_ratio)
    !
    ! DESCRIPTION:
    ! Calculates dead:live ratio of fuel
    !
    
    ! ARGUMENTS:
    class(synthetic_fuel_model), intent(in)  :: this
    real(r8),                    intent(out) :: dead_live_ratio
    
    ! LOCALS:
    integer  :: i                  ! looping index
    real(r8) :: live_loading       ! live fuel loading, weighted by SAV 
    real(r8) :: dead_loading       ! dead fuel loading, weighted by SAV
    
    ! calculate dead:live ratio
    dead_loading = 0.0_r8 
    live_loading = 0.0_r8
    do i = 1, 3
      if (this%dead_fuel_loading(i) > 0.0_r8) then
        dead_loading = dead_loading + this%dead_fuel_weighting(i)*exp(-4.528*this%dead_fuel_sav(i))
      end if
    end do 
    
    do i = 1, 2
      if (this%live_fuel_loading(i) > 0.0_r8) then
        live_loading = live_loading + this%live_fuel_weighting(i)*exp(-16.4_r8*this%live_fuel_sav(i))
      end if
    end do
    
    if (live_loading > 0.0_r8) then 
      dead_live_ratio = dead_loading/live_loading
    else
      dead_live_ratio = 0.0_r8 
    end if
    
  end subroutine CalculateDeadLiveRatio
  
  ! --------------------------------------------------------------------------------------
    
  subroutine AddFuelModel(this, fuel_model_index, carrier, fuel_model_name,              &
    wind_adj_factor, hr1_loading, hr10_loading, hr100_loading, live_herb_loading,        &
    live_woody_loading, fuel_depth, hr1_sav, live_herb_sav, live_woody_sav,              &
    moist_extinct)
    !
    ! DESCRIPTION:
    ! Adds a fuel model to the dynamic array
    !
    ! NOTE THE UNITS ON INPUTS
    !
    
    ! ARGUMENTS:
    class(fuel_models_array_class), intent(inout) :: this               ! array of fuel models
    integer,                        intent(in)    :: fuel_model_index   ! fuel model index
    character(len=2),               intent(in)    :: carrier            ! main carrier
    character(len=*),               intent(in)    :: fuel_model_name    ! fuel model long name
    real(r8),                       intent(in)    :: wind_adj_factor    ! wind adjustment factor
    real(r8),                       intent(in)    :: hr1_loading        ! loading for 1-hr fuels [US tons/acre]
    real(r8),                       intent(in)    :: hr10_loading       ! loading for 10-hr fuels [US tons/acre]
    real(r8),                       intent(in)    :: hr100_loading      ! loading for 100-hr fuels [US tons/acre]
    real(r8),                       intent(in)    :: live_herb_loading  ! loading for live herbacious fuels [us tons/acre]
    real(r8),                       intent(in)    :: live_woody_loading ! loading for live woody fuels [US tons/acre]
    real(r8),                       intent(in)    :: fuel_depth         ! fuel bed depth [ft]
    real(r8),                       intent(in)    :: hr1_sav            ! surface area to volume ratio of 1 hour fuels [/ft]
    real(r8),                       intent(in)    :: live_herb_sav      ! surface area to volume ratio of live herbacious fuels [/ft]
    real(r8),                       intent(in)    :: live_woody_sav     ! surface area to volume ratio of live woody fuels [/ft]
    real(r8),                       intent(in)    :: moist_extinct      ! dead fuel extinction moisture [%]
    
    ! LOCALS:
    type(synthetic_fuel_model)              :: fuel_model         ! fuel model
    type(synthetic_fuel_model), allocatable :: temporary_array(:) ! temporary array to hold data while re-allocating
    
    ! first make sure we have enough space in the array
    if (allocated(this%fuel_models)) then
      ! already allocated to some size
      if (this%num_fuel_models == size(this%fuel_models)) then 
        ! need to add more space
        allocate(temporary_array(size(this%fuel_models) + chunk_size))
        temporary_array(1:size(this%fuel_models)) = this%fuel_models
        call move_alloc(temporary_array, this%fuel_models)
      end if 
      
      this%num_fuel_models = this%num_fuel_models + 1
  
    else 
      ! first element in array 
      allocate(this%fuel_models(chunk_size))
      this%num_fuel_models = 1
    end if 
    
    call fuel_model%InitFuelModel(fuel_model_index, carrier, fuel_model_name,            &
      wind_adj_factor, hr1_loading, hr10_loading, hr100_loading, live_herb_loading,      &
      live_woody_loading, fuel_depth, hr1_sav, live_herb_sav, live_woody_sav, moist_extinct)
    
    this%fuel_models(this%num_fuel_models) = fuel_model
      
  end subroutine AddFuelModel
  
  ! --------------------------------------------------------------------------------------
  
  integer function FuelModelPosition(this, fuel_model_index)
    !
    ! DESCRIPTION:
    ! Returns the index of a desired fuel model
    !
    
    ! ARGUMENTS:
    class(fuel_models_array_class), intent(in)  :: this             ! array of fuel models
    integer,                        intent(in)  :: fuel_model_index ! desired fuel model index
        
    ! LOCALS:
    integer :: i ! looping index 
    
    do i = 1, this%num_fuel_models
      if (this%fuel_models(i)%fuel_model_index == fuel_model_index) then
        FuelModelPosition = i
        return
      end if
    end do
    write(*, '(a, i2, a)') "Cannot find the fuel model index ", fuel_model_index, "."
    stop
  
  end function FuelModelPosition
  
  ! --------------------------------------------------------------------------------------
  
  subroutine GetFuelModels(this)
    !
    ! DESCRIPTION:
    ! Returns an array of hard-coded fuel models
    ! these are taken from the fire behavior fuel models in Scott & Burgan 2005
    !
    
    ! ARGUMENTS:
    class(fuel_models_array_class), intent(inout) :: this ! array of fuel models
    
    ! governed by fine herbaceous fuels that are cured or nearly cured
    ! fires move rapidly through cured grass and associated material
    ! very little shrub or timber present
    ! grasslands/savanna are represented along with stubble, grass tundra, grass-shrub combinations
    ! annual and perennial grasses are included fuels
    call this%AddFuelModel(fuel_model_index=1, carrier='GR', fuel_model_name='short grass',     &
      wind_adj_factor=0.36_r8, hr1_loading=0.74_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,   &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=1.0_r8,                   &
      hr1_sav=3500.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=30.0_r8)
      
    ! fire spread primarily through herbaceous fuels, either curing or dead
    ! surface fires
    call this%AddFuelModel(fuel_model_index=2, carrier='GR', fuel_model_name='timber and grass understory', &
      wind_adj_factor=0.36_r8, hr1_loading=2.0_r8, hr10_loading=1.0_r8, hr100_loading=0.5_r8,               &
      live_herb_loading=0.5_r8, live_woody_loading=0.0_r8, fuel_depth=1.0_r8,                               &
      hr1_sav=3000.0_r8, live_herb_sav=1500.0_r8, live_woody_sav=0.0_r8, moist_extinct=15.0_r8)
    
    ! most intense of grass group - high rates of spread under influence of wind
    call this%AddFuelModel(fuel_model_index=3, carrier='GR', fuel_model_name='tall grass',    &
      wind_adj_factor=0.44_r8, hr1_loading=3.0_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8, &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=2.5_r8,                 &
      hr1_sav=1500.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=25.0_r8)
    
    ! fire intensity and fast-spreading fires involve foliage and live and dead fine woody material in shrub layer
    ! dead woody material contributes significantly to fire intensity
    ! deep litter layer may confound suppression efforts
    call this%AddFuelModel(fuel_model_index=4, carrier='SH', fuel_model_name='chapparal',     &
      wind_adj_factor=0.55_r8, hr1_loading=5.0_r8, hr10_loading=4.0_r8, hr100_loading=2.0_r8, &
      live_herb_loading=0.0_r8, live_woody_loading=5.0_r8, fuel_depth=6.0_r8,                 &
      hr1_sav=2000.0_r8, live_herb_sav=0.0_r8, live_woody_sav=1500.0_r8, moist_extinct=20.0_r8)
    
    ! primary carrier is litter from shrubs, and grasses/forbs in understory
    ! shrubs generally not tall, but have nearly total coverage of area
    ! young, green stands with no deadwood
    call this%AddFuelModel(fuel_model_index=5, carrier='SH', fuel_model_name='brush',         &
      wind_adj_factor=0.42_r8, hr1_loading=1.0_r8, hr10_loading=0.5_r8, hr100_loading=0.0_r8, &
      live_herb_loading=0.0_r8, live_woody_loading=2.0_r8, fuel_depth=2.0_r8,                 &
      hr1_sav=2000.0_r8, live_herb_sav=0.0_r8, live_woody_sav=1500.0_r8, moist_extinct=20.0_r8)
    
    ! fire carries through shrub layer, requiring at least moderate winds
    ! fire will drop to ground at low wind speeds or openings in stand
    ! shrubs are older
    call this%AddFuelModel(fuel_model_index=6, carrier='SH', fuel_model_name='dormant brush', &
      wind_adj_factor=0.44_r8, hr1_loading=1.5_r8, hr10_loading=2.5_r8, hr100_loading=2.0_r8, &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=2.5_r8,                 &
      hr1_sav=1750.0_r8, live_herb_sav=0.0_r8, live_woody_sav=1500.0_r8, moist_extinct=25.0_r8)
    
    ! fires burn through surface and shrub strata with equal ease and can occur at higher dead
    ! fuel moisture contents because of flammable nature of live foliage and other live material
    ! stands of shrubs generally between 2-6 ft high
    call this%AddFuelModel(fuel_model_index=7, carrier='SH', fuel_model_name='southern rough', &
      wind_adj_factor=0.44_r8, hr1_loading=1.1_r8, hr10_loading=1.9_r8, hr100_loading=1.0_r8,  &
      live_herb_loading=0.0_r8, live_woody_loading=0.4_r8, fuel_depth=2.5_r8,                  &
      hr1_sav=1750.0_r8, live_herb_sav=0.0_r8, live_woody_sav=1500.0_r8, moist_extinct=40.0_r8)
    
    ! slow-burning ground fires with low flame heights
    ! layer mainly needles, leaves, and some twigs
    call this%AddFuelModel(fuel_model_index=8, carrier='TL', fuel_model_name='compact timber litter', &
      wind_adj_factor=0.28_r8, hr1_loading=1.5_r8, hr10_loading=1.0_r8, hr100_loading=2.5_r8,         &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.2_r8,                         &
      hr1_sav=2000.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=30.0_r8)
    
    ! long-needle conifer and hardwood stands typical
    call this%AddFuelModel(fuel_model_index=9, carrier='TL', fuel_model_name='hardwood litter', &
      wind_adj_factor=0.28_r8, hr1_loading=2.9_r8, hr10_loading=0.4_r8, hr100_loading=0.2_r8,   &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.2_r8,                   &
      hr1_sav=2500.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=25.0_r8)
    
    ! dead down fuels include greater quantities of 3-inch or larger limbwood resulting from over 
    ! maturity or natural events that create a large load of dead material
    ! crown fire and spotting is more frequent
    call this%AddFuelModel(fuel_model_index=10, carrier='TU', fuel_model_name='timber and litter understorey', &
      wind_adj_factor=0.46_r8, hr1_loading=3.0_r8, hr10_loading=2.0_r8, hr100_loading=5.0_r8,                  &
      live_herb_loading=0.0_r8, live_woody_loading=2.0_r8, fuel_depth=1.0_r8,                                  &
      hr1_sav=2000.0_r8, live_herb_sav=0.0_r8, live_woody_sav=1500.0_r8, moist_extinct=25.0_r8)
    
    ! fires are fairly active in the slash and intermixed herbaceous material
    ! spacing of the rather light fuel load, shading from overstory, or the aging of the 
    ! fine fuels can contribute to limiting the fire potential 
    call this%AddFuelModel(fuel_model_index=11, carrier='SB', fuel_model_name='light slash',  &
      wind_adj_factor=0.36_r8, hr1_loading=1.5_r8, hr10_loading=4.5_r8, hr100_loading=5.5_r8, &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=1.0_r8,                 &
      hr1_sav=1500.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=15.0_r8)
    
    ! rapidly spreading fires with high intensities capable of generating firebrands can occur
    ! when fire starts, it is generally sustained until a fuel break or change in fuels is encountered
    ! visual impression is dominated by slash, most of it less than 3 inches in diameter
    call this%AddFuelModel(fuel_model_index=12, carrier='SB', fuel_model_name='medium slash',   &
      wind_adj_factor=0.43_r8, hr1_loading=4.0_r8, hr10_loading=14.0_r8, hr100_loading=16.5_r8, &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=2.3_r8,                 &
      hr1_sav=1500.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=20.0_r8)
    
    ! fire is generally carried across the area by a continuous layer of slash 
    ! large quantities of greater-than-3-inch material are present
    ! active flaming is sustained for long periods and firebrands of various sizes may be generated
    ! these contribute to spotting problems
    call this%AddFuelModel(fuel_model_index=13, carrier='SB', fuel_model_name='heavy slash',    &
      wind_adj_factor=0.46_r8, hr1_loading=7.0_r8, hr10_loading=23.0_r8, hr100_loading=28.1_r8, &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=3.0_r8,                 &
      hr1_sav=1500.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=25.0_r8)
      
      
    ! Scott and Burgen -----------
    
    ! primary carrier of fire is sparse grass - small amounts of dead fuel may be present
    ! grass is short (naturally or through grazing)
    call this%AddFuelModel(fuel_model_index=101, carrier='GR', fuel_model_name='short, sparse dry climate grass', &
      wind_adj_factor=0.31_r8, hr1_loading=0.1_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                     &
      live_herb_loading=0.3_r8, live_woody_loading=0.0_r8, fuel_depth=0.4_r8,                                     &
      hr1_sav=2200.0_r8, live_herb_sav=2000.0_r8, live_woody_sav=0.0_r8, moist_extinct=15.0_r8)
    
    ! primary carrier of fire is grass, small amounts of dead fuel may be present
    ! fuelbed more continuous than GR101
    ! shrubs do not affect fire behavior
    call this%AddFuelModel(fuel_model_index=102, carrier='GR', fuel_model_name='low load dry climate grass',  &
      wind_adj_factor=0.36_r8, hr1_loading=0.1_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                 &
      live_herb_loading=1.0_r8, live_woody_loading=0.0_r8, fuel_depth=1.0_r8,                                 &
      hr1_sav=2000.0_r8, live_herb_sav=1800.0_r8, live_woody_sav=0.0_r8, moist_extinct=15.0_r8)
    
    ! primary carrier of fire is continuous, coarse, humid-climate grass
    ! grass/herb fuel load light
    ! shrubs do not affect fire behavior
    call this%AddFuelModel(fuel_model_index=103, carrier='GR', fuel_model_name='low load very coarse humid climate grass',  &
      wind_adj_factor=0.42_r8, hr1_loading=0.1_r8, hr10_loading=0.4_r8, hr100_loading=0.0_r8,                               &
      live_herb_loading=1.5_r8, live_woody_loading=0.0_r8, fuel_depth=2.0_r8,                                               &
      hr1_sav=1500.0_r8, live_herb_sav=1300.0_r8, live_woody_sav=0.0_r8, moist_extinct=30.0_r8)
    
    ! primary carrier is continuous, dry-climate grass; load/depth greater than GR102
    call this%AddFuelModel(fuel_model_index=104, carrier='GR', fuel_model_name='moderate load dry climate grass',  &
      wind_adj_factor=0.42_r8, hr1_loading=0.25_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                      &
      live_herb_loading=1.9_r8, live_woody_loading=0.0_r8, fuel_depth=2.0_r8,                                      &
      hr1_sav=2000.0_r8, live_herb_sav=1800.0_r8, live_woody_sav=0.0_r8, moist_extinct=15.0_r8)
    
    ! primary carrier is humid-climate grass
    ! load greater than GR103, depth is smaller
    call this%AddFuelModel(fuel_model_index=105, carrier='GR', fuel_model_name='low load humid climate grass',  &
      wind_adj_factor=0.39_r8, hr1_loading=0.4_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                   &
      live_herb_loading=2.5_r8, live_woody_loading=0.0_r8, fuel_depth=1.5_r8,                                   &
      hr1_sav=1800.0_r8, live_herb_sav=1600.0_r8, live_woody_sav=0.0_r8, moist_extinct=40.0_r8)
    
    ! primary carrier is continuous humid-climate grass; load greater than GR105, depth about the same
    ! grass is less coarse
    call this%AddFuelModel(fuel_model_index=106, carrier='GR', fuel_model_name='moderate load humid climate grass',  &
      wind_adj_factor=0.39_r8, hr1_loading=0.1_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                        &
      live_herb_loading=3.4_r8, live_woody_loading=0.0_r8, fuel_depth=1.5_r8,                                        &
      hr1_sav=2200.0_r8, live_herb_sav=2000.0_r8, live_woody_sav=0.0_r8, moist_extinct=40.0_r8)
    
    ! primary carrier is continuous dry-climate grass, grass about 3 ft tall
    call this%AddFuelModel(fuel_model_index=107, carrier='GR', fuel_model_name='high load dry climate grass',  &
      wind_adj_factor=0.46_r8, hr1_loading=1.0_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                  &
      live_herb_loading=5.4_r8, live_woody_loading=0.0_r8, fuel_depth=3.0_r8,                                  &
      hr1_sav=2000.0_r8, live_herb_sav=1800.0_r8, live_woody_sav=0.0_r8, moist_extinct=15.0_r8)
    
    ! primary carrier is continuous, very coarse, humid-climate grass
    ! spread rate and flame length can be extreme if grass is fully cured
    call this%AddFuelModel(fuel_model_index=108, carrier='GR', fuel_model_name='high load humid climate grass', &
      wind_adj_factor=0.49_r8, hr1_loading=0.5_r8, hr10_loading=1.0_r8, hr100_loading=0.0_r8,                   &
      live_herb_loading=7.3_r8, live_woody_loading=0.0_r8, fuel_depth=4.0_r8,                                   &
      hr1_sav=1500.0_r8, live_herb_sav=1300.0_r8, live_woody_sav=0.0_r8, moist_extinct=30.0_r8)
    
    ! primary carrier is dense, tall, humid-climate grass
    ! spread rate and flame length can be extreme if grass is fully or mostly cured
    call this%AddFuelModel(fuel_model_index=109, carrier='GR', fuel_model_name='very high load humid climate grass-shrub', &
      wind_adj_factor=0.52_r8, hr1_loading=1.0_r8, hr10_loading=1.0_r8, hr100_loading=0.0_r8,                              &
      live_herb_loading=9.0_r8, live_woody_loading=0.0_r8, fuel_depth=5.0_r8,                                              &
      hr1_sav=1800.0_r8, live_herb_sav=1600.0_r8, live_woody_sav=0.0_r8, moist_extinct=40.0_r8)
      
    ! primary carrier is grass and shrubs combined, shrubs about 1 ft, grass load is low
    ! spread rate is moderate, flame length low; low moisture of extinction
    call this%AddFuelModel(fuel_model_index=121, carrier='GS', fuel_model_name='low load dry climate grass-shrub',  &
      wind_adj_factor=0.35_r8, hr1_loading=0.2_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                       &
      live_herb_loading=0.5_r8, live_woody_loading=0.65_r8, fuel_depth=0.9_r8,                                      &
      hr1_sav=2000.0_r8, live_herb_sav=1800.0_r8, live_woody_sav=1800.0_r8, moist_extinct=15.0_r8)
    
    ! primary carrier is grass and shrubs combined, shrubs 1-3 ft; grass load moderate
    ! spread rate is high, flame length moderate; low moisture of extinction
    call this%AddFuelModel(fuel_model_index=122, carrier='GS', fuel_model_name='moderate load dry climate grass-shrub',  &
      wind_adj_factor=0.39_r8, hr1_loading=0.5_r8, hr10_loading=0.5_r8, hr100_loading=0.0_r8,                            &
      live_herb_loading=0.6_r8, live_woody_loading=1.0_r8, fuel_depth=1.5_r8,                                            &
      hr1_sav=2000.0_r8, live_herb_sav=1800.0_r8, live_woody_sav=1800.0_r8, moist_extinct=15.0_r8)
    
    ! primary carrier is grass and shrubs combined, moderate grass/shrub load - average depth less than 2 ft
    ! spread rate high; flame length moderate; high moisture of extinction
    call this%AddFuelModel(fuel_model_index=123, carrier='GS', fuel_model_name='moderate load humid climate grass-shrub', &
      wind_adj_factor=0.41_r8, hr1_loading=0.3_r8, hr10_loading=0.25_r8, hr100_loading=0.0_r8,                            &
      live_herb_loading=1.45_r8, live_woody_loading=1.25_r8, fuel_depth=1.8_r8,                                           &
      hr1_sav=1800.0_r8, live_herb_sav=1600.0_r8, live_woody_sav=1600.0_r8, moist_extinct=40.0_r8)
    
    ! primary carrier is grass and shrubs combined; heavy grass/shrub load - greater than 2 ft depth
    ! spread rate high; flame length very high; high moisture of extinction
    call this%AddFuelModel(fuel_model_index=124, carrier='GS', fuel_model_name='high load humid climate grass-shrub',  &
      wind_adj_factor=0.42_r8, hr1_loading=1.9_r8, hr10_loading=0.3_r8, hr100_loading=0.1_r8,                          &
      live_herb_loading=3.4_r8, live_woody_loading=7.1_r8, fuel_depth=2.1_r8,                                          &
      hr1_sav=1800.0_r8, live_herb_sav=1600.0_r8, live_woody_sav=1600.0_r8, moist_extinct=40.0_r8)
    
    ! primary carrier is woody shrubs and shrub litter
    ! low shrub fuel load, some grass may be present
    ! spread rate is very low; flame length very low
    call this%AddFuelModel(fuel_model_index=141, carrier='SH', fuel_model_name='low load dry climate shrub', &
      wind_adj_factor=0.36_r8, hr1_loading=0.25_r8, hr10_loading=0.25_r8, hr100_loading=0.0_r8,              &
      live_herb_loading=0.15_r8, live_woody_loading=1.3_r8, fuel_depth=1.0_r8,                               &
      hr1_sav=2000.0_r8, live_herb_sav=1800.0_r8, live_woody_sav=1600.0_r8, moist_extinct=15.0_r8)
    
    ! primary carrier is woody shrubs and shrub litter
    ! moderate fuel load, no grass fuel
    ! spread rate is low, flame length is low
    call this%AddFuelModel(fuel_model_index=142, carrier='SH', fuel_model_name='moderate load dry climate shrub', &
      wind_adj_factor=0.36_r8, hr1_loading=1.35_r8, hr10_loading=2.4_r8, hr100_loading=0.75_r8,                   &
      live_herb_loading=0.0_r8, live_woody_loading=3.85_r8, fuel_depth=1.0_r8,                                    &
      hr1_sav=2000.0_r8, live_herb_sav=0.0_r8, live_woody_sav=1600.0_r8, moist_extinct=15.0_r8)
    
    ! primary carrier is woody shrubs and shrub litter
    ! moderate shrub load, possibly pine overstory and herbaceous fuel
    ! spread rate is low; flame length low
    call this%AddFuelModel(fuel_model_index=143, carrier='SH', fuel_model_name='moderate load humid climate shrub', &
      wind_adj_factor=0.44_r8, hr1_loading=0.45_r8, hr10_loading=3.0_r8, hr100_loading=0.0_r8,                       &
      live_herb_loading=0.0_r8, live_woody_loading=6.2_r8, fuel_depth=2.4_r8,                                        &
      hr1_sav=1600.0_r8, live_herb_sav=0.0_r8, live_woody_sav=1400.0_r8, moist_extinct=40.0_r8)
    
    ! primary carrier is woody shrubs and shrub litter
    ! low to moderate shrub and litter load, possibly with pine overstory
    ! spread rate high; flame length moderate
    call this%AddFuelModel(fuel_model_index=144, carrier='SH', fuel_model_name='low load humid climate timber-shrub', &
      wind_adj_factor=0.46_r8, hr1_loading=0.85_r8, hr10_loading=1.15_r8, hr100_loading=0.2_r8,                       &
      live_herb_loading=0.0_r8, live_woody_loading=2.55_r8, fuel_depth=3.0_r8,                                        &
      hr1_sav=2000.0_r8, live_herb_sav=1800.0_r8, live_woody_sav=1600.0_r8, moist_extinct=30.0_r8)
    
    ! primary carrier is woody shrubs and shrub litter
    ! heavy shrub load
    ! spread rate very high; flame length very high
    call this%AddFuelModel(fuel_model_index=145, carrier='SH', fuel_model_name='high load dry climate shrub', &
      wind_adj_factor=0.55_r8, hr1_loading=3.6_r8, hr10_loading=2.1_r8, hr100_loading=0.0_r8,                 &
      live_herb_loading=0.0_r8, live_woody_loading=2.9_r8, fuel_depth=6.0_r8,                                 &
      hr1_sav=750.0_r8, live_herb_sav=0.0_r8, live_woody_sav=1600.0_r8, moist_extinct=15.0_r8)  
    
    ! primary carrier is woody shrubs and shrub litter
    ! dense shrubs; little or no herbaceous fuel
    ! spread rate high; flame length high
    call this%AddFuelModel(fuel_model_index=146, carrier='SH', fuel_model_name='low load humid climate shrub', &
      wind_adj_factor=0.42_r8, hr1_loading=2.9_r8, hr10_loading=1.45_r8, hr100_loading=0.0_r8,                 &
      live_herb_loading=0.0_r8, live_woody_loading=1.4_r8, fuel_depth=2.0_r8,                                  &
      hr1_sav=750.0_r8, live_herb_sav=0.0_r8, live_woody_sav=1600.0_r8, moist_extinct=30.0_r8)   
    
    ! primary carrier is woody shrubs and shrub litter
    ! very heavy shrub load
    ! spread rate high, flame length very high
    call this%AddFuelModel(fuel_model_index=147, carrier='SH', fuel_model_name='very high load dry climate shrub', &
      wind_adj_factor=0.55_r8, hr1_loading=3.5_r8, hr10_loading=5.3_r8, hr100_loading=2.2_r8,                      &
      live_herb_loading=0.0_r8, live_woody_loading=3.4_r8, fuel_depth=6.0_r8,                                      &
      hr1_sav=750.0_r8, live_herb_sav=0.0_r8, live_woody_sav=1600.0_r8, moist_extinct=15.0_r8)
    
    ! primary carrier is woody shrubs and shrub litter
    ! dense shrubs, little or no herbaceous fuel
    ! spread rate high; flame length high
    call this%AddFuelModel(fuel_model_index=148, carrier='SH', fuel_model_name='high load humid climate shrub', &
      wind_adj_factor=0.46_r8, hr1_loading=2.05_r8, hr10_loading=3.4_r8, hr100_loading=0.85_r8,                 &
      live_herb_loading=0.0_r8, live_woody_loading=4.35_r8, fuel_depth=3.0_r8,                                  &
      hr1_sav=750.0_r8, live_herb_sav=0.0_r8, live_woody_sav=1600.0_r8, moist_extinct=40.0_r8)
    
    ! primary carrier woody shrubs and shrub litter
    ! dense, finely branched shrubs with significant fine dead fuel; some herbaceous fuel may be present
    ! spread rate high; flame length very high
    call this%AddFuelModel(fuel_model_index=149, carrier='SH', fuel_model_name='very high load humid climate shrub', &
      wind_adj_factor=0.5_r8, hr1_loading=4.5_r8, hr10_loading=2.45_r8, hr100_loading=0.0_r8,                        &
      live_herb_loading=1.55_r8, live_woody_loading=7.0_r8, fuel_depth=4.4_r8,                                       &
      hr1_sav=750.0_r8, live_herb_sav=1800.0_r8, live_woody_sav=1500.0_r8, moist_extinct=40.0_r8) 
    
    ! primary carrier is low load grass and/or shrub with litter
    ! spread rate low; flame length low
    call this%AddFuelModel(fuel_model_index=161, carrier='TU', fuel_model_name='light load dry climate timber-grass-shrub', &
      wind_adj_factor=0.33_r8, hr1_loading=0.2_r8, hr10_loading=0.9_r8, hr100_loading=1.5_r8,                               &
      live_herb_loading=0.2_r8, live_woody_loading=0.9_r8, fuel_depth=0.6_r8,                                               &
      hr1_sav=2000.0_r8, live_herb_sav=1800.0_r8, live_woody_sav=1600.0_r8, moist_extinct=20.0_r8)
    
    ! primary carrier is moderate litter load with shrub component
    ! spread rate moderate; flame length low
    call this%AddFuelModel(fuel_model_index=162, carrier='TU', fuel_model_name='moderate load humid climate timber-shrub', &
      wind_adj_factor=0.36_r8, hr1_loading=0.95_r8, hr10_loading=1.8_r8, hr100_loading=1.25_r8,                            &
      live_herb_loading=0.0_r8, live_woody_loading=0.2_r8, fuel_depth=1.0_r8,                                              &
      hr1_sav=2000.0_r8, live_herb_sav=0.0_r8, live_woody_sav=1600.0_r8, moist_extinct=20.0_r8)
    
    ! primary carrier is moderate forest litter with grass and shrub components
    ! spread rate high; flame length moderate
    call this%AddFuelModel(fuel_model_index=163, carrier='TU', fuel_model_name='moderate load humid climate timber-grass-shrub', &
      wind_adj_factor=0.38_r8, hr1_loading=1.1_r8, hr10_loading=0.15_r8, hr100_loading=0.25_r8,                                  &
      live_herb_loading=0.65_r8, live_woody_loading=1.1_r8, fuel_depth=1.3_r8,                                                   &
      hr1_sav=1800.0_r8, live_herb_sav=1600.0_r8, live_woody_sav=1400.0_r8, moist_extinct=30.0_r8)
    
    ! primary carrier is short conifer trees with grass or moss understory
    ! spread rate moderate; flame length moderate
    call this%AddFuelModel(fuel_model_index=164, carrier='TU', fuel_model_name='dwarf conifer with understory', &
      wind_adj_factor=0.32_r8, hr1_loading=4.5_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                   &
      live_herb_loading=0.0_r8, live_woody_loading=2.0_r8, fuel_depth=0.5_r8,                                   &
      hr1_sav=2300.0_r8, live_herb_sav=0.0_r8, live_woody_sav=2000.0_r8, moist_extinct=12.0_r8)
 
    ! primary carrier is heavy forest litter with a shrub or small tree understory
    ! spread rate moderate; flame length moderate
    call this%AddFuelModel(fuel_model_index=165, carrier='TU', fuel_model_name='very high load dry climate timber-shrub', &
      wind_adj_factor=0.33_r8, hr1_loading=4.0_r8, hr10_loading=4.0_r8, hr100_loading=3.0_r8,                             &
      live_herb_loading=0.0_r8, live_woody_loading=3.0_r8, fuel_depth=1.0_r8,                                             &
      hr1_sav=1500.0_r8, live_herb_sav=0.0_r8, live_woody_sav=750.0_r8, moist_extinct=25.0_r8)
    
    ! primary carrier is compact forest litter; light to moderate load
    ! may be used to represent a recently burned forest
    ! spread rate is very low; flame length very low
    call this%AddFuelModel(fuel_model_index=181, carrier='TL', fuel_model_name='low load compact conifer litter', &
      wind_adj_factor=0.28_r8, hr1_loading=1.0_r8, hr10_loading=2.2_r8, hr100_loading=3.6_r8,                     &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.2_r8,                                     &
      hr1_sav=2000.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=30.0_r8)
    
    ! primary carrier is broadleaf (hardwood) litter; low load
    ! compact broadleaf litter
    ! spread rate is very low; flame length is very low
    call this%AddFuelModel(fuel_model_index=182, carrier='TL', fuel_model_name='low load broadleaf litter', &
      wind_adj_factor=0.28_r8, hr1_loading=1.4_r8, hr10_loading=2.3_r8, hr100_loading=2.2_r8,               &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.2_r8,                               &
      hr1_sav=2000.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=25.0_r8)
    
    ! primary carrier is moderate load conifer litter, light load of coarse fuels
    ! spread rate is very low; flame length low
    call this%AddFuelModel(fuel_model_index=183, carrier='TL', fuel_model_name='moderate load conifer litter', &
      wind_adj_factor=0.29_r8, hr1_loading=0.5_r8, hr10_loading=2.2_r8, hr100_loading=2.8_r8,                  &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.3_r8,                                  &
      hr1_sav=2000.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=20.0_r8)
    
    ! primary carrier is moderate load of fine litter and coarse fuels
    ! includes small diameter downed logs
    ! spread rate is low; flame length low
    call this%AddFuelModel(fuel_model_index=184, carrier='TL', fuel_model_name='small downed logs', &
      wind_adj_factor=0.31_r8, hr1_loading=0.5_r8, hr10_loading=1.5_r8, hr100_loading=4.2_r8,       &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.4_r8,                       &
      hr1_sav=2000.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=25.0_r8)
    
    ! primary carrier is high load of conifer litter; light slash or mortality fuel
    ! spread rate is low; flame length low
    call this%AddFuelModel(fuel_model_index=185, carrier='TL', fuel_model_name='high load conifer litter', &
      wind_adj_factor=0.33_r8, hr1_loading=1.15_r8, hr10_loading=2.5_r8, hr100_loading=4.4_r8,             &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.6_r8,                              &
      hr1_sav=2000.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=25.0_r8)
    
    ! primary carrier is moderate load broadleaf litter, less compact than TL182
    ! spread rate is moderate; flame length low
    call this%AddFuelModel(fuel_model_index=186, carrier='TL', fuel_model_name='moderate load broadleaf litter', &
      wind_adj_factor=0.29_r8, hr1_loading=2.4_r8, hr10_loading=1.2_r8, hr100_loading=1.2_r8,                    &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.3_r8,                                    &
      hr1_sav=2000.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=25.0_r8)
      
    ! primary carrier is heavy load forest litter; includes large diameter downed logs
    ! spread rate low; flame length low
    call this%AddFuelModel(fuel_model_index=187, carrier='TL', fuel_model_name='large downed logs', &
      wind_adj_factor=0.31_r8, hr1_loading=0.3_r8, hr10_loading=1.4_r8, hr100_loading=8.1_r8,       &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.4_r8,                       &
      hr1_sav=2000.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=25.0_r8)
    
    ! primary carrier is moderate load long-needle pine litter, may include a small amount of herbaceous load
    ! spread rate is moderate; flame length low
    call this%AddFuelModel(fuel_model_index=188, carrier='TL', fuel_model_name='long-needle litter', &
      wind_adj_factor=0.29_r8, hr1_loading=5.8_r8, hr10_loading=1.4_r8, hr100_loading=1.1_r8,        &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.3_r8,                        &
      hr1_sav=1800.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=35.0_r8)
    
    ! primary carrier is very high load, fluffy broadleaf litter
    ! can also be used to represent heavy needle-drape
    ! spread rate is moderate; flame length moderate
    call this%AddFuelModel(fuel_model_index=189, carrier='TL', fuel_model_name='very high load broadleaf litter', &
      wind_adj_factor=0.33_r8, hr1_loading=6.65_r8, hr10_loading=3.3_r8, hr100_loading=4.15_r8,                   &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.6_r8,                                     &
      hr1_sav=1800.0_r8, live_herb_sav=0.0_r8, live_woody_sav=1600.0_r8, moist_extinct=35.0_r8)
    
    ! primary carrier is light dead and down activity fuel
    ! fine fuel load is 10-20 us t/acre, weighted towards 1-3 inch. diameter class
    ! spread rate is moderate, flame length low
    call this%AddFuelModel(fuel_model_index=201, carrier='SB', fuel_model_name='low load activity fuel', &
      wind_adj_factor=0.36_r8, hr1_loading=1.5_r8, hr10_loading=3.0_r8, hr100_loading=11.1_r8,           &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=1.0_r8,                            &
      hr1_sav=2000.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=25.0_r8)
    
    ! primary carrier is moderate dead and down activity fuel or light blowdown
    ! fine fuel load is evenly distributed across diameter classes
    ! blowdown is scattered with many trees still standing
    ! spread rate is moderate; flame length moderate
    call this%AddFuelModel(fuel_model_index=202, carrier='SB', fuel_model_name='moderate load activity fuel or low load blowdown', &
      wind_adj_factor=0.36_r8, hr1_loading=4.5_r8, hr10_loading=4.25_r8, hr100_loading=4.0_r8,                                     &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=1.0_r8,                                                      &
      hr1_sav=2000.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=25.0_r8)
    
    ! primary carrier is heavy dead and down activity fuel or moderate blowdown
    ! fine fuel load is weighted towards 0-0.25 inch diameter class
    ! blowdown is moderate, trees compacted to near the ground
    ! spread rate is high; flame length high
    call this%AddFuelModel(fuel_model_index=203, carrier='SB', fuel_model_name='high load activity fuel or moderate load blowdown', &
      wind_adj_factor=0.38_r8, hr1_loading=5.5_r8, hr10_loading=2.75_r8, hr100_loading=3.0_r8,                                      &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=1.2_r8,                                                       &
      hr1_sav=2000.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=25.0_r8)
    
    ! primary carrier is heavy blowdown fuel
    ! blowdown is total, fuelbed not compacted; most foliage and fine fuel still attached to blowdown
    ! spread rate very high; flame length very high
    call this%AddFuelModel(fuel_model_index=204, carrier='SB', fuel_model_name='high load blowdown', &
      wind_adj_factor=0.45_r8, hr1_loading=5.25_r8, hr10_loading=3.5_r8, hr100_loading=5.25_r8,      &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=2.7_r8,                        &
      hr1_sav=2000.0_r8, live_herb_sav=0.0_r8, live_woody_sav=0.0_r8, moist_extinct=25.0_r8)
      
    end subroutine GetFuelModels
  
  ! --------------------------------------------------------------------------------------
    
end module SyntheticFuelModels
