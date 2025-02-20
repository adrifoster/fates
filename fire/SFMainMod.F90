module SFMainMod

  ! ============================================================================
  ! All subroutines related to the SPITFIRE fire routine. 
  ! Code originally developed by Allan Spessa & Rosie Fisher as part of the NERC-QUEST project.  
  ! ============================================================================

  use FatesConstantsMod,      only : r8 => fates_r8
  use FatesConstantsMod,      only : itrue, ifalse
  use FatesConstantsMod,      only : pi_const
  use FatesConstantsMod,      only : nocomp_bareground, nearzero
  use FatesGlobals,           only : fates_log
  use FatesInterfaceTypesMod, only : hlm_masterproc 
  use FatesInterfaceTypesMod, only : hlm_spitfire_mode
  use FatesInterfaceTypesMod, only : hlm_sf_nofire_def
  use FatesInterfaceTypesMod, only : hlm_sf_scalar_lightning_def
  use FatesInterfaceTypesMod, only : hlm_sf_successful_ignitions_def
  use FatesInterfaceTypesMod, only : hlm_sf_anthro_ignitions_def
  use FatesInterfaceTypesMod, only : bc_in_type
  use EDPftvarcon,            only : EDPftvarcon_inst
  use PRTParametersMod,       only : prt_params
  use PRTGenericMod,          only : element_pos
  use EDtypesMod,             only : ed_site_type
  use FatesPatchMod,          only : fates_patch_type
  use FatesCohortMod,         only : fates_cohort_type
  use EDtypesMod,             only : AREA
  use FatesLitterMod,         only : litter_type
  use FatesFuelClassesMod,    only : num_fuel_classes
  use PRTGenericMod,          only : leaf_organ
  use PRTGenericMod,          only : carbon12_element
  use PRTGenericMod,          only : sapw_organ
  use PRTGenericMod,          only : struct_organ
  use FatesInterfaceTypesMod, only : numpft
  use FatesAllometryMod,      only : CrownDepth
  use FatesFuelClassesMod,    only : fuel_classes
  
  implicit none
  private

  public :: DailyFireModel
  public :: UpdateFuelCharacteristics

  ! ======================================================================================

contains

  subroutine DailyFireModel(currentSite, bc_in)
    !
    !  DESCRIPTION:
    !  Runs the daily fire model

    ! ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite ! site object
    type(bc_in_type),   intent(in)            :: bc_in       ! BC in object

    ! LOCALS:  
    type (fates_patch_type), pointer :: currentPatch ! patch object
        
    if (hlm_spitfire_mode > hlm_sf_nofire_def) then
      call UpdateFireWeather(currentSite, bc_in)
      call UpdateFuelCharacteristics(currentSite)
      call CalculateIgnitionsandFDI(currentSite, bc_in)
      call CalculateSurfaceRateOfSpread(currentSite)
      call CalculateSurfaceFireIntensity(currentSite)
      call CalculateAreaBurnt(currentSite)
      
      
      ! call crown_scorching(currentSite) ! calculates SH
      ! call crown_damage(currentSite)    ! calculates cambial damage rate
      ! call cambial_damage_kill(currentSite)  ! cambial damage probi
      ! call post_fire_mortality(currentSite)
    end if

  end subroutine DailyFireModel

  !---------------------------------------------------------------------------------------
 
  subroutine UpdateFireWeather(currentSite, bc_in)
    !
    !  DESCRIPTION:
    !  Updates the site's fire weather index and calculates effective windspeed based on 
    !   vegetation characteristics
    !
    !  Currently we use tree and grass fraction averaged over whole grid (site) to 
    !  prevent extreme divergence

    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    use FatesConstantsMod, only : sec_per_day, sec_per_min
    use EDTypesMod,        only : CalculateTreeGrassAreaSite

    ! ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite
    type(bc_in_type),   intent(in)            :: bc_in

    ! LOCALS:  
    type(fates_patch_type), pointer :: currentPatch   ! patch object
    real(r8)                        :: temp_C         ! daily averaged temperature [deg C]
    real(r8)                        :: precip         ! daily precip [mm/day]
    real(r8)                        :: rh             ! daily relative humidity [%]
    real(r8)                        :: wind           ! wind speed [m/s]
    real(r8)                        :: tree_fraction  ! site-level tree fraction [0-1]
    real(r8)                        :: grass_fraction ! site-level grass fraction [0-1]
    real(r8)                        :: bare_fraction  ! site-level bare ground fraction [0-1]
    integer                         :: iofp           ! index of oldest the fates patch

    ! NOTE that the boundary conditions of temperature, precipitation and relative humidity
    ! are available at the patch level. We are currently using a simplification where the whole site
    ! is simply using the values associated with the first patch.
    ! which probably won't have much impact, unless we decide to ever calculated fire weather for each patch.  

    currentPatch => currentSite%oldest_patch

    ! If the oldest patch is a bareground patch (i.e. nocomp mode is on) use the first vegetated patch
    ! for the iofp index (i.e. the next younger patch)
    if (currentPatch%nocomp_pft_label == nocomp_bareground) then
      currentPatch => currentPatch%younger
    endif

    iofp = currentPatch%patchno
    temp_C = currentPatch%tveg24%GetMean() - tfrz
    precip = bc_in%precip24_pa(iofp)*sec_per_day
    rh = bc_in%relhumid24_pa(iofp)
    wind = bc_in%wind24_pa(iofp)

    ! convert to m/min 
    currentSite%wind = wind*sec_per_min

    ! update fire weather index
    call currentSite%fireWeather%UpdateIndex(temp_C, precip, rh, wind)

    ! calculate site-level tree, grass, and bare fraction
    call CalculateTreeGrassAreaSite(currentSite, tree_fraction, grass_fraction, bare_fraction)

    ! update effective wind speed
    call currentSite%fireWeather%UpdateEffectiveWindSpeed(wind*sec_per_min, tree_fraction, &
      grass_fraction, bare_fraction)

  end subroutine UpdateFireWeather

  !---------------------------------------------------------------------------------------
  
  subroutine UpdateFuelCharacteristics(currentSite)
    !
    !  DESCRIPTION:
    !  Updates fuel characteristics on each patch of the site
    !

    use SFParamsMod, only : SF_val_drying_ratio, SF_val_SAV, SF_val_FBD

    ! ARGUMENTS:
    type(ed_site_type), intent(in), target :: currentSite  ! site object

    ! LOCALS:
    type(fates_patch_type), pointer :: currentPatch ! FATES patch 
    type(litter_type),      pointer :: litter       ! pointer to patch litter class
    real(r8) :: MEF_trunks, fuel_moisture_trunks
    
    currentPatch => currentSite%oldest_patch 
    do while(associated(currentPatch))  

      if (currentPatch%nocomp_pft_label /= nocomp_bareground) then

        ! calculate live grass [kgC/m2]
        call currentPatch%UpdateLiveGrass()

        ! update fuel loading [kgC/m2]
        litter => currentPatch%litter(element_pos(carbon12_element))
        call currentPatch%fuel%UpdateLoading(sum(litter%leaf_fines(:)),                  &
          litter%ag_cwd(1), litter%ag_cwd(2), litter%ag_cwd(3), litter%ag_cwd(4),        &
          currentPatch%livegrass)
            
        ! sum up fuel classes and calculate fractional loading for each
        call currentPatch%fuel%SumLoading()
        call currentPatch%fuel%CalculateFractionalLoading()
          
        ! calculate fuel moisture [m3/m3]
        call currentPatch%fuel%UpdateFuelMoisture(SF_val_SAV, SF_val_drying_ratio,       &
          currentSite%fireWeather)
        
        ! calculate geometric properties
        call currentPatch%fuel%AverageBulkDensity_NoTrunks(SF_val_FBD)
        call currentPatch%fuel%AverageSAV_NoTrunks(SF_val_SAV)
            
      end if 
      currentPatch => currentPatch%younger
    end do 

  end subroutine UpdateFuelCharacteristics

  !---------------------------------------------------------------------------------------
  
  !****************************************************************
  subroutine  characteristics_of_crown ( currentSite )
  !****************************************************************.  

    !returns the live crown fuel characteristics within each patch.
    ! passive_crown_FI is minimum fire intensity to ignite canopy crown fuel

    use SFParamsMod,    only : SF_VAL_CWD_FRAC

    type(ed_site_type), intent(in), target :: currentSite

    type(ed_patch_type) , pointer :: currentPatch
    type(ed_cohort_type), pointer :: currentCohort

    ! ARGUMENTS

    ! LOCAL
    real(r8) ::  crown_depth          ! depth of crown (m)
    real(r8) ::  height_cbb           ! clear branch bole height or crown base height (m)
    real(r8) ::  max_height           ! max cohort on patch (m)
    real(r8) ::  crown_ignite_energy  ! heat yield for crown (kJ/kg)
    real(r8) ::  tree_sapw_struct_c   ! above-ground tree struct and sap biomass in cohort (kgC)
    real(r8) ::  leaf_c                  ! leaf carbon (kgC)
    real(r8) ::  sapw_c                  ! sapwood carbon (kgC)
    real(r8) ::  struct_c                ! structure carbon (kgC)
    real(r8) ::  twig_sapw_struct_c      ! above-ground twig sap and struct in cohort (kgC)
    real(r8) ::  crown_fuel_c            ! biomass of 1 hr fuels (leaves,twigs) in cohort (kg C)
    real(r8) ::  crown_fuel_biomass      ! biomass of crown fuel in cohort (kg biomass)
    real(r8) ::  crown_fuel_per_m        ! crown fuel per 1m section in cohort
    real(r8) ::  height_base_canopy      ! lowest height of fuels in patch to carry fire in crown

    integer  ::  ih                      ! counter

    real, dimension(70):: biom_matrix   ! matrix to track biomass from bottom to 70m
    real(r8),parameter :: min_density_canopy_fuel = 0.011_r8 !min canopy fuel density (kg/m3) sufficient to
                                                             !propogate fire vertically through canopy
                                                             !Scott and Reinhardt 2001 RMRS-RP-29
    real(r8),parameter :: foliar_moist_content = 1.0_r8      !foliar moisture content default 100% Scott & Reinhardt 2001


    !returns the live crown fuel characteristics within each patch.
    ! passive_crown_FI is the required minimum fire intensity to ignite canopy crown fuel

    currentPatch => currentSite%oldest_patch

    do while(associated(currentPatch))
       !zero Patch level variables
       height_base_canopy                   = 0.0_r8
       currentPatch%canopy_fuel_load        = 0.0_r8
       currentPatch%passive_crown_FI        = 0.0_r8
       currentPatch%canopy_bulk_density     = 0.0_r8
       max_height = 0._r8
       biom_matrix(:) = 0._r8

          currentCohort=>currentPatch%tallest
          do while(associated(currentCohort))

             !zero cohort level variables
             tree_sapw_struct_c                   = 0.0_r8
             leaf_c                               = 0.0_r8
             sapw_c                               = 0.0_r8
             struct_c                             = 0.0_r8
             twig_sapw_struct_c                   = 0.0_r8
             crown_fuel_c                         = 0.0_r8
             crown_fuel_biomass                   = 0.0_r8
             crown_fuel_per_m                     = 0.0_r8

             ! Calculate crown 1hr fuel biomass (leaf, twig sapwood, twig structural biomass)
             if ( int(prt_params%woody(currentCohort%pft)) == itrue) then !trees

                call CrownDepth(currentCohort%hite,currentCohort%pft,crown_depth)
                height_cbb   = currentCohort%hite - crown_depth

                !find patch max height for stand canopy fuel
                if (currentCohort%hite > max_height) then
                   max_height = currentCohort%hite
                endif

                leaf_c   = currentCohort%prt%GetState(leaf_organ, all_carbon_elements)
                sapw_c   = currentCohort%prt%GetState(sapw_organ, all_carbon_elements)
                struct_c = currentCohort%prt%GetState(struct_organ, all_carbon_elements)

                tree_sapw_struct_c =  currentCohort%n * &
                        (prt_params%allom_agb_frac(currentCohort%pft)*(sapw_c + struct_c))

                twig_sapw_struct_c =  tree_sapw_struct_c * SF_VAL_CWD_frac(1)   !only 1hr fuel

                crown_fuel_c = (currentCohort%n * leaf_c) + twig_sapw_struct_c  !crown fuel (kgC)

                crown_fuel_biomass = crown_fuel_c / 0.45_r8            ! crown fuel (kg biomass)

                crown_fuel_per_m = crown_fuel_biomass / crown_depth    ! kg biomass per m

                !sort crown fuel into bins from bottom to top of crown
                !accumulate across cohorts to find density within canopy 1m sections
                do ih = max(1, nint(height_cbb)), max(1, nint(currentCohort%hite))
                   biom_matrix(ih) = biom_matrix(ih) + crown_fuel_per_m
                end do

                !accumulate available canopy fuel for patch (kg biomass)
                ! use this in CFB (crown fraction burn) calculation and FI final
                currentPatch%canopy_fuel_load = currentPatch%canopy_fuel_load + crown_fuel_biomass  !canopy fuel in patch

             endif !trees only

             currentCohort => currentCohort%shorter;

          enddo !end cohort loop

          biom_matrix(:) = biom_matrix(:) / currentPatch%area    !kg biomass/m3

          !loop from 1m to 70m to find bin with total density = 0.011 kg/m3
          !min canopy fuel density to propogate fire vertically in canopy across patch
          do ih=1,70
             if (biom_matrix(ih) > min_density_canopy_fuel) then
                height_base_canopy = float(ih)
                exit
             end if
          end do

          !canopy_bulk_denisty (kg/m3) for Patch
          if (max_height - height_base_canopy > 0._r8) then
             currentPatch%canopy_bulk_density = sum(biom_matrix) / (max_height - height_base_canopy)
          else
             currentPatch%canopy_bulk_density = 0._r8
          end if

          ! Note: crown_ignition_energy to be calculated based on PFT foliar moisture content from FATES-Hydro
          ! or create foliar moisture % based on BTRAN
          ! Use foliar_moisture(currentCohort%pft) and compute weighted PFT average with Eq 3 Van Wagner 1977
          ! in place of foliar_moist_content parameter

          ! Eq 3 Van Wagner 1977, Eq 11 Scott & Reinhardt 2001
          ! h = 460.0 + 25.9*m
          ! h = crown_ignite_energy (kJ/kg), m = foliar moisture content based on dry fuel (%)
          crown_ignite_energy = 460.0 + 25.9 * foliar_moist_content

          ! Crown fuel ignition potential (kW/m), Eq 4 Van Wagner 1977, Eq 11 Scott & Reinhardt 2001
          ! FI = (Czh)**3/2 where z=canopy base height,h=heat of crown ignite energy, FI=fire intensity
          ! 0.01 = C, empirical constant Van Wagner 1977 Eq 4 for 6m canopy base height, 100% FMC, FI 2500kW/m
          ! passive_crown_FI = min fire intensity to ignite canopy fuel (kW/m or kJ/m/s)
          currentPatch%passive_crown_FI = (0.01_r8 * height_base_canopy * crown_ignite_energy)**1.5_r8
      
      currentPatch => currentPatch%younger;

    enddo !end patch loop

  end subroutine characteristics_of_crown
  
  !---------------------------------------------------------------------------------------
  
  !*****************************************************************
  subroutine  active_crown_fire ( currentSite)
  !*****************************************************************

    !evaluates if there will be an active crown fire based on canopy fuel and rate of spread
    !returns final rate of spread and fire intensity in patch with added fuel from active crown fire.
    !currentCohort%fraction_crown_burned is the proportion of crown affected by fire

    use EDParamsMod, only : active_crown_fire_switch
    use SFParamsMod, only  : SF_val_miner_total, SF_val_part_dens, SF_val_miner_damp, &
                             SF_val_fuel_energy, SF_val_drying_ratio

 
    type(ed_site_type), intent(in), target :: currentSite

    type(ed_patch_type) , pointer :: currentPatch
    type(ed_cohort_type), pointer :: currentCohort

    ! ARGUMENTS

    ! LOCAL VARIABLES
    ! Active crown Rothermel fire spread model parameters using FM 10
    real(r8) beta,beta_op         ! weighted average of packing ratio (unitless)
    real(r8) ir                   ! reaction intensity (kJ/m2/min)
    real(r8) xi,eps,phi_wind      ! all are unitless
    real(r8) q_ig                 ! heat of pre-ignition (kJ/kg)
    real(r8) reaction_v_opt,reaction_v_max !reaction velocity (per min)!optimum and maximum
    real(r8) moist_damp,mw_weight ! moisture dampening coefficient and ratio fuel moisture to extinction
    real(r8) beta_ratio           ! ratio of beta/beta_op
    real(r8) a_beta               ! dummy variable for product of a* beta_ratio for react_v_opt equation
    real(r8) a,b,c,e              ! function of fuel sav
    real(r8) total_fuel           ! total fuel (kg biomass/m2)
    real(r8) net_fuel             ! net fuel (kg biomass/m2) without minerals
    real(r8) fuel_depth           ! fuel depth (m)
    real(r8) fuel_bd              ! fuel bulk density (kg biomass/m3)
    real(r8) fuel_sav             ! fuels average sav 
    real(r8) fuel_eff_moist       ! fuels effective moisture
    real(r8) fuel_moist1hr        ! moisture 1 hour fuels
    real(r8) fuel_moist10hr       ! moisture 10 hour fuels
    real(r8) fuel_moist100hr      ! moisture 100 hour fuels
    real(r8) fuel_moistlive       ! moisture live fuels
    real(r8) fuel_1hr             ! FM 10 1-hr fuel loading (kg/m2)
    real(r8) fuel_10hr            ! FM 10 10-hr fuel loading (kg/m2)
    real(r8) fuel_100hr           ! FM 10 100-hr fuel loading (kg/m2)
    real(r8) fuel_live            ! FM 10 live fuel loading (kg/m2)
    real(r8) SAV_1hr              ! surface area to volume 1 hour fuels (twigs)
    real(r8) SAV_10hr             ! surface area to volume 10 hour fuels (small branches)
    real(r8) SAV_100hr            ! surface area to volume 100 hour fuels (large branches)
    real(r8) SAV_live             ! surface area to volume live fuels
    real(r8) midflame_wind        ! 40% of open wind speed, Scott & Reinhardt 2001 
    real(r8) db                   ! distance fire has traveld backward (m)
    real(r8) df                   ! distance fire has travelled forward (m)
    real(r8) AB                   ! daily area burnt (m2 per km2)  
    real(r8) size_of_fire         ! in m2
    real(r8) ROS_active           ! actual rate of spread (m/min) using FM 10 fuels
    real(r8) ROS_active_min       ! minimum rate of spread to ignite active crown fire
    real(r8) CI_temp              ! temporary variable to calculate wind_active_min
    real(r8) wind_active_min      ! open windspeed to sustain active crown fire where ROS_SA = ROS_active_min
    real(r8) ROS_SA               ! rate of spread for surface fire with wind_active_min
    real(r8) canopy_frac_burnt    ! fraction of canopy fuels consumed (0, surface fire to 1,active crown fire) 
    real(r8) ROS_final            ! final rate of spread for combined surface and canopy spread (m/min)
    real(r8) FI_final             ! final fireline intensity (kW/m or kJ/m/sec) with canopy consumption 

    real(r8),parameter :: q_dry = 581.0_r8                 !heat of pre-ignition of dry fuels (kJ/kg)
    ! fuel loading, MEF, and depth from Anderson 1982 Aids to determining fuel models for fire behavior
    ! SAV values from BEHAVE model Burgan & Rothermel (1984)
    real(r8),parameter :: fuel_1hr_ton   = 3.01_r8           ! FM 10 1-hr fuel loading (US tons/acre)
    real(r8),parameter :: fuel_10hr_ton  = 2.0_r8            ! FM 10 10-hr fuel loading (US tons/acre)
    real(r8),parameter :: fuel_100hr_ton = 5.01_r8           ! FM 10 100-hr fuel loading (US tons/acre)
    real(r8),parameter :: fuel_live_ton  = 2.0_r8            ! FM 10 live fuel loading (US tons/acre)
    real(r8),parameter :: fuel_mef     = 0.25_r8             ! FM 10 moisture of extinction (volumetric)
    real(r8),parameter :: fuel_depth_ft= 1.0_r8              ! FM 10 fuel depth (ft)
    real(r8),parameter :: sav_1hr_ft   = 2000.0_r8           ! FM 10 1-hr SAV (ft2/ft3)
    real(r8),parameter :: sav_10hr_ft  = 109.0_r8            ! FM 10 10-hr SAV (ft2/ft3)             
    real(r8),parameter :: sav_100hr_ft = 30.0_r8             ! FM 10 100-hr SAV (ft2/ft3)
    real(r8),parameter :: sav_live_ft  = 1650.0_r8           ! FM 10 live SAV (ft2/ft3)
    real(r8),parameter :: tonnes_acre_to_kg_m2 = 0.2241701   ! convert tons/acre to kg/m2
    real(r8),parameter :: sqft_cubicft_to_sqm_cubicm = 0.03280844 !convert ft2/ft3 to m2/m3
    real(r8),parameter :: canopy_ignite_energy = 18000.0_r8  ! heat yield for canopy fuels (kJ/kg)
    real(r8),parameter :: critical_mass_flow_rate = 0.05_r8  ! critical mass flow rate (kg/m2/sec)for crown fire
    real(r8),parameter :: km2_to_m2 = 1000000.0_r8           ! area conversion for square km to square m

    integer  :: passive_canopy_fuel_flg                    ! flag if canopy fuel true for vertical spread


    currentPatch => currentSite%oldest_patch

    do while(associated(currentPatch))

       if (currentPatch%fire == 1) then
          passive_canopy_fuel_flg = 0         !does patch have canopy fuels for vertical spread?
          ROS_active = 0.0_r8

          ! check initiation of passive crown fire
          if (currentPatch%FI >= currentPatch%passive_crown_FI) then
             passive_canopy_fuel_flg = 1      !enough passive canopy fuels for vertical spread

             ! Calculate rate of spread using FM 10 as in Rothermel 1977 
             ! fuel characteristics 
             fuel_1hr   = fuel_1hr_ton * tonnes_acre_to_kg_m2
             fuel_10hr  = fuel_10hr_ton * tonnes_acre_to_kg_m2
             fuel_100hr = fuel_100hr_ton * tonnes_acre_to_kg_m2
             fuel_live  = fuel_live_ton * tonnes_acre_to_kg_m2

             total_fuel = fuel_1hr + fuel_10hr + fuel_100hr + fuel_live  !total fuel (kg/m2)

             SAV_1hr   = sav_1hr_ft * sqft_cubicft_to_sqm_cubicm
             SAV_10hr  = sav_10hr_ft * sqft_cubicft_to_sqm_cubicm
             SAV_100hr = sav_100hr_ft * sqft_cubicft_to_sqm_cubicm
             SAV_live  = sav_live_ft * sqft_cubicft_to_sqm_cubicm

             fuel_moist1hr    = exp(-1.0_r8 * ((SAV_1hr/SF_val_drying_ratio) * currentSite%acc_NI))
             fuel_moist10hr   = exp(-1.0_r8 * ((SAV_10hr/SF_val_drying_ratio) * currentSite%acc_NI))
             fuel_moist100hr  = exp(-1.0_r8 * ((SAV_100hr/SF_val_drying_ratio) * currentSite%acc_NI))
             fuel_moistlive   = exp(-1.0_r8 * ((SAV_live/SF_val_drying_ratio) * currentSite%acc_NI))

             fuel_depth       = fuel_depth_ft *0.3048           !convert to meters
             fuel_bd          = total_fuel / fuel_depth         !fuel bulk density (kg/m3)

             fuel_sav         = SAV_1hr *(fuel_1hr/total_fuel) + SAV_10hr*(fuel_10hr/total_fuel) + & 
                                 SAV_100hr*(fuel_100hr/total_fuel) + SAV_live*(fuel_live/total_fuel)

             fuel_eff_moist = fuel_moist1hr *(fuel_1hr/total_fuel) + fuel_moist10hr*(fuel_10hr/total_fuel) + & 
                               fuel_moist100hr*(fuel_100hr/total_fuel) + fuel_moistlive*(fuel_live/total_fuel)

             ! remove mineral content from net fuel load
             net_fuel = total_fuel * (1.0_r8 - SF_val_miner_total) !net of minerals

             ! ---start spreading---
             !beta = packing ratio (unitless)
             beta = fuel_bd / SF_val_part_dens
             beta_op = 0.200395_r8 *(fuel_sav**(-0.8189_r8))
             beta_ratio = beta/beta_op  

             ! -- heat of pre-ignition --
             q_ig = q_dry + 2594.0_r8 * fuel_eff_moist

             ! ---effective heating number---
             ! Eq A3 in Thonicke et al. 2010.  
             eps = exp(-4.528_r8 / fuel_sav)     
             ! Eq A7 in Thonicke et al. 2010 per Eq 49, Rothermel 1972
             b = 0.15988_r8 * (fuel_sav**0.54_r8)
             ! Eq A8 in Thonicke et al. 2010 per Eq 48, Rothermel 1972 
             c = 7.47_r8 * (exp(-0.8711_r8 * (fuel_sav**0.55_r8))) 
             ! Eq A9 in Thonicke et al. 2010. (typo in Eq A9, using coefficient Eq 50, Rothermel 1972)
             e = 0.715_r8 * (exp(-0.01094_r8 * fuel_sav))

             midflame_wind = currentSite%wind *0.40_r8  !Scott & Reinhardt 2001 40% open wind speed

             ! Eq A5 in Thonicke et al. 2010
             ! include convert wind from m/min to ft/min for Rothermel ROS eqn
             phi_wind = c * ((3.281_r8*midflame_wind)**b)*(beta_ratio**(-e)) !unitless

             ! ---propagating flux = xi (unitless) 
             ! Eq A2 in Thonicke et al.2010 and Eq 42 Rothermel 1972
             xi = (exp((0.792_r8 + 3.7597_r8 * (fuel_sav**0.5_r8)) * (beta+0.1_r8))) / &
                (192_r8+7.9095_r8 * fuel_sav) 

             ! ---reaction intensity----
             ! Eq in table A1 Thonicke et al. 2010. 
             a = 8.9033_r8 * (fuel_sav**(-0.7913_r8))
             a_beta = exp(a*(1.0_r8-beta_ratio))  !dummy variable for reaction_v_opt equation
  
             ! Eq in table A1 Thonicke et al. 2010.
             ! reaction_v_max and reaction_v_opt = reaction velocity in units of per min
             ! reaction_v_max = Eq 36 in Rothermel 1972 and Fig 12 
             reaction_v_max  = 1.0_r8 / (0.0591_r8 + 2.926_r8* (fuel_sav**(-1.5_r8)))
             ! reaction_v_opt =  Eq 38 in Rothermel 1972 and Fig 11
             reaction_v_opt = reaction_v_max*(beta_ratio**a)*a_beta

             ! mw_weight = relative fuel moisture/fuel moisture of extinction
             mw_weight = fuel_eff_moist/fuel_mef
       
             ! Eq in table A1 Thonicke et al. 2010. (unitless)
             moist_damp = max(0.0_r8,(1.0_r8 - (2.59_r8 * mw_weight) + (5.11_r8 * (mw_weight**2.0_r8)) - &
               (3.52_r8*(mw_weight**3.0_r8))))

             ! ir = reaction intenisty in kJ/m2/min
             ! sum_fuel as kgBiomass/m2 for ir calculation
             ir = reaction_v_opt*(net_fuel)*SF_val_fuel_energy*moist_damp*SF_val_miner_damp  
 
             ! actual ROS (m/min) for FM 10 fuels for open windspeed, Eq 8 Scott & Reinhardt 2001
             ROS_active = 3.34_r8 * ((ir*xi*(1.0_r8+phi_wind)) / (fuel_bd * eps * q_ig))

             ! critical min rate of spread (m/min) for active crowning
             ROS_active_min = (critical_mass_flow_rate / fuel_bd) * 60.0_r8

             ! check threshold intensity and rate of spread
             if (active_crown_fire_switch .and. &
                 ROS_active >= ROS_active_min) then
                currentPatch%active_crown_fire_flg = 1  ! active crown fire ignited
                !ROS_final = ROS_surface+CFB(ROS_active - ROS_surface), Eq 21 Scott & Reinhardt 2001
                !with active crown fire CFB (canopy fraction burned) = 100%
                canopy_frac_burnt = 1.0_r8

             else 
                currentPatch%active_crown_fire_flg = 0  ! only passive crown fire with partial crown burnt

                ! phi_slope is not used yet. consider adding with later
                ! development
                ! calculate open wind speed critical to sustain active crown
                ! fire Eq 20 Scott & Reinhardt
                if (ir > 0._r8 .and. currentPatch%canopy_bulk_density > 0._r8) then
                   CI_temp = ((164.8_r8 * eps * q_ig)/(ir * currentPatch%canopy_bulk_density)) - 1.0_r8
                else
                   CI_temp = 0._r8
                end if

                wind_active_min = 0.0457_r8 * (CI_temp / 0.001612_r8)**0.7_r8

                ! use open wind speed "wind_active_min" for ROS surface fire
                ! where ROS_SA=ROS_active_min
                ROS_SA =  (ir * xi * (1.0_r8 + wind_active_min)) / (fuel_bd * eps * q_ig)

                ! canopy fraction burnt, Eq 28 Scott & Reinhardt Appendix A
                canopy_frac_burnt = max(0._r8, min(1.0_r8, &
                   (currentPatch%ROS_front - ROS_active_min) / (ROS_SA - ROS_active_min)))
                
             endif !check intensity & ROS for active crown fire thresholds

             !ROS_final = ROS_surface+CFB(ROS_active - ROS_surface), Eq 21 Scott & Reinhardt 2001
             ROS_final = currentPatch%ROS_front + &
                canopy_frac_burnt * (ROS_active - currentPatch%ROS_front)

             ! recalculate area burned with new ROS_front value from ROS_final
             ! ---- re-calculate length of major axis for df using new ROS_front value from ROS final---
             db = currentPatch%ROS_back  * currentPatch%FD !(m) 
             df = ROS_final * currentPatch%FD              !(m) update with ROS final 
             
             ! update ROS_front with ROS_final for output variable 
             ! if changing, expect only an increase in ROS_front
             currentPatch%ROS_front = max(ROS_final, currentPatch%ROS_front)

             ! --- calculate updated area burnt using df from ROS final---
             if (currentPatch%lb > 0.0_r8) then

                ! Eq 1 in Thonicke et al. 2010
                ! To Do: Connect here with the Li & Levis GDP fire suppression algorithm. 
                ! Eq 16 in arora and boer model JGR 2005
                ! AB = AB *3.0_r8

                !size of fire = Eq 14 Arora and Boer JGR 2005 (area of an ellipse)
                size_of_fire = (df + db) * (df + db) * pi_const / (4.0_r8 * currentPatch%lb)

                ! AB = daily area burnt = size fires in m2 * num ignitions per day per km2 * prob ignition starts fire
                ! AB = m2 per km2 per day
                ! the denominator in the units of currentSite%NF is total gridcell area, but since we assume that ignitions 
                ! are equally probable across patches, currentSite%NF is equivalently per area of a given patch
                ! thus AB has units of m2 burned area per km2 patch area per day
                AB = size_of_fire * currentSite%NF * currentSite%FDI

                ! frac_burnt 
                ! just a unit conversion from AB, to become area burned per area patch per day, 
                ! or just the fraction of the patch burned on that day
                currentPatch%frac_burnt = (min(0.99_r8, AB / km2_to_m2))
             
                if(write_SF == itrue)then
                   if ( hlm_masterproc == itrue ) write(fates_log(),*) 'frac_burnt',currentPatch%frac_burnt
                endif

             else  
                currentPatch%frac_burnt = 0.0_r8
             endif ! lb

             !final fireline intensity (kJ/m/sec or kW/m), Eq 22 Scott & Reinhardt 2001
             FI_final = (currentPatch%heat_per_area + &
  currentPatch%canopy_fuel_load * canopy_ignite_energy * canopy_frac_burnt) * &
  currentPatch%ROS_front / 60._r8
             ! update patch FI to adjust according to potential canopy fuel consumed (passive and active)
             ! if changing, expect only an increase in FI
             currentPatch%FI = max(FI_final, currentPatch%FI)

           endif !check if passive crown fire?
       endif !fire?

       currentPatch => currentPatch%younger;

    enddo !end patch loop

  end subroutine active_crown_fire
  
  !---------------------------------------------------------------------------------------
  
  subroutine CalculateIgnitionsandFDI(currentSite, bc_in)
    !
    !  DESCRIPTION:
    !  Calculates ignitions and fire danger index (FDI) for a site
    !
    
    use FatesInterfaceTypesMod, only : hlm_spitfire_mode
    use EDParamsMod,            only : cg_strikes
    use EDParamsMod,            only : ED_val_nignitions
    use SFParamsMod,            only : SF_val_fdi_alpha
    use FatesConstantsMod,      only : years_per_day

    ! ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite ! site object
    type(bc_in_type),   intent(in)            :: bc_in       ! BC in object
    
    ! LOCALS:
    type(fates_patch_type), pointer :: currentPatch            ! patch object
    real(r8)                        :: cloud_to_ground_strikes ! fraction of cloud-to-ground strikes [0-1]
    real(r8)                        :: anthro_ignitions        ! anthropogenic ignitions [count/km2/day]
    integer                         :: iofp                    ! patch index
    
    ! CONSTANTS:
    real(r8), parameter :: alpha = 0.0035_r8  ! potential human ignition counts (alpha in Li et al. 2012) (#/person/month)

    ! initialize site parameters to zero
    currentSite%NF_successful = 0.0_r8

    ! Equation 7 from Venevsky et al GCB 2002 (modification of equation 8 in Thonicke et al. 2010) 
    ! FDI 0.1 = low, 0.3 moderate, 0.75 high, and 1 = extreme ignition potential for alpha 0.000337
    if (hlm_spitfire_mode == hlm_sf_successful_ignitions_def) then
      ! READING "SUCCESSFUL IGNITION" DATA
      ! force ignition potential to be extreme
      ! cloud_to_ground_strikes = 1 means using 100% of incoming observed ignitions
      currentSite%FDI = 1.0_r8  
      cloud_to_ground_strikes = 1.0_r8   
    else  
      ! USING LIGHTNING STRIKE DATA
      currentSite%FDI  = 1.0_r8 - exp(-SF_val_fdi_alpha*currentSite%fireWeather%fire_weather_index)
      cloud_to_ground_strikes = cg_strikes
    end if

    ! if the oldest patch is a bareground patch (i.e. nocomp mode is on) use the first vegetated patch
    ! for the iofp index (i.e. the next younger patch)
    currentPatch => currentSite%oldest_patch
    if(currentPatch%nocomp_pft_label .eq. nocomp_bareground)then
      currentPatch => currentPatch%younger
    endif
    iofp = currentPatch%patchno

    ! NF = number of lighting strikes per day per km2 scaled by cloud to ground strikes
    if (hlm_spitfire_mode == hlm_sf_scalar_lightning_def) then
      currentSite%NF = ED_val_nignitions*years_per_day*cloud_to_ground_strikes
    else    
      ! use external daily lightning ignition data
      currentSite%NF = bc_in%lightning24(iofp)*cloud_to_ground_strikes
    end if

    ! calculate anthropogenic ignitions according to Li et al. (2012)
    ! add to ignitions by lightning
    if (hlm_spitfire_mode == hlm_sf_anthro_ignitions_def) then
      ! anthropogenic ignitions (count/km2/day)
      !           =  (ignitions/person/month)*6.8*population_density**0.43/approximate days per month
      anthro_ignitions = alpha*6.8_r8*bc_in%pop_density(iofp)**0.43_r8/30.0_r8
      currentSite%NF = currentSite%NF + anthro_ignitions
    end if

  end subroutine CalculateIgnitionsandFDI
  
  !---------------------------------------------------------------------------------------
  
  subroutine CalculateSurfaceRateOfSpread(currentSite) 
    !
    !  DESCRIPTION:
    !  Calculates potential rate of spread based on fuel characteristics for 
    !  each patch of a site
    !

    use SFParamsMod,    only : SF_val_miner_total, SF_val_part_dens
    use SFEquationsMod, only : OptimumPackingRatio, ReactionIntensity
    use SFEquationsMod, only : HeatofPreignition, EffectiveHeatingNumber
    use SFEquationsMod, only : WindFactor, PropagatingFlux
    use SFEquationsMod, only : ForwardRateOfSpread, BackwardRateOfSpread

    ! ARGUMENTS:
    type(ed_site_type), intent(in), target :: currentSite ! site object

    ! LOCALS:
    type(fates_patch_type), pointer :: currentPatch ! patch object 
    real(r8)                        :: beta         ! packing ratio [unitless]
    real(r8)                        :: beta_op      ! optimum packing ratio [unitless]
    real(r8)                        :: beta_ratio   ! relative packing ratio [unitless]
    real(r8)                        :: i_r          ! reaction intensity [kJ/m2/min]
    real(r8)                        :: xi           ! propagating flux ratio [unitless]
    real(r8)                        :: eps          ! effective heating number [unitless]
    real(r8)                        :: phi_wind     ! wind factor [unitless]
    real(r8)                        :: q_ig         ! heat of pre-ignition [kJ/kg]

    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))
      if (currentPatch%nocomp_pft_label /= nocomp_bareground .and.                       &
        currentPatch%fuel%non_trunk_loading > nearzero) then
        
        ! fraction of fuel array volume occupied by fuel, i.e. compactness of fuel bed [unitless]
        ! Rothermel 1972 Eq. 31
        beta = currentPatch%fuel%bulk_density_notrunks/SF_val_part_dens
        
        ! optimum packing ratio [unitless]
        beta_op = OptimumPackingRatio(currentPatch%fuel%SAV_notrunks)
        
        ! relative packing ratio [unitless]
        if (beta_op < nearzero) then 
          beta_ratio = 0.0_r8
        else
          beta_ratio = beta/beta_op 
        end if
        
        ! remove mineral content from fuel load per Thonicke 2010 
        currentPatch%fuel%non_trunk_loading = currentPatch%fuel%non_trunk_loading*(1.0_r8 - SF_val_miner_total) 
        
        ! reaction intensity [kJ/m2/min]
        i_r = ReactionIntensity(currentPatch%fuel%non_trunk_loading/0.45_r8,             &
          currentPatch%fuel%SAV_notrunks, beta_ratio,                                    &
          currentPatch%fuel%average_moisture_notrunks, currentPatch%fuel%MEF_notrunks)
   
        ! heat of preignition [kJ/kg] 
        q_ig = HeatofPreignition(currentPatch%fuel%average_moisture_notrunks)

        ! effective heating number [unitless]
        eps = EffectiveHeatingNumber(currentPatch%fuel%SAV_notrunks)
        
        ! wind factor [unitless]
        phi_wind = WindFactor(currentSite%fireWeather%effective_windspeed, beta_ratio,      &
         currentPatch%fuel%SAV_notrunks)

        ! propagating flux [unitless]       
        xi = PropagatingFlux(beta, currentPatch%fuel%SAV_notrunks)
        
        ! forward rate of spread [m/min]
        currentPatch%ROS_front = ForwardRateOfSpread(currentPatch%fuel%bulk_density_notrunks, &
         eps, q_ig, i_r, xi, phi_wind)

        ! backwards rate of spread [m/min]
        !  backward ROS wind not changed by vegetation - so use wind, not effective_windspeed
        currentPatch%ROS_back = BackwardRateOfSpread(currentPatch%ROS_front,             &
         currentSite%wind)

      end if 
      currentPatch => currentPatch%younger
    end do

  end subroutine CalculateSurfaceRateOfSpread
  
  !---------------------------------------------------------------------------------------
  
  subroutine CalculateSurfaceFireIntensity(currentSite)
    !
    !  DESCRIPTION:
    !  Calculates surface fireline intensity for each patch of a site
    !  Right now also calculates the area burnt...
    !
    use FatesConstantsMod, only : m2_per_km2
    use SFEquationsMod,    only : FireDuration, LengthToBreadth
    use SFEquationsMod,    only : AreaBurnt, FireSize, FireIntensity
    use SFParamsMod,       only : SF_val_fire_threshold

    ! ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite
    
    ! LOCALS:
    type(fates_patch_type), pointer :: currentPatch                    ! patch object
    real(r8)                        :: fuel_consumed(num_fuel_classes) ! fuel consumed [kgC/m2]
    real(r8)                        :: tree_fraction_patch             ! treed fraction on patch [0-1]
    real(r8)                        :: length_to_breadth               ! length to breadth ratio of fire ellipse (unitless)
    real(r8)                        :: fire_size                       ! size of fire [m2]
    real(r8)                        :: area_burnt                      ! area burnt [m2/km2]
    
    currentPatch => currentSite%oldest_patch 
    do while (associated(currentPatch))
      
      currentPatch%fuel%frac_burnt(:) = 0.0_r8

      if (currentPatch%nocomp_pft_label /= nocomp_bareground) then

        call currentPatch%fuel%CalculateFuelBurnt(fuel_consumed)
        call currentPatch%fuel%CalculateResidenceTime(currentPatch%tau_l)

        ! calculate overall fuel consumed by spreading fire
        ! ignore 1000-hr fuels (i.e. trunks)
        currentPatch%TFC_ROS = sum(fuel_consumed) - fuel_consumed(fuel_classes%trunks())  

        ! initialize patch parameters to zero
        currentPatch%FI = 0.0_r8 
        currentPatch%fire = 0
        
        if (currentSite%NF > 0.0_r8) then
          
          ! fire intensity [kW/m]
          currentPatch%FI = FireIntensity(currentPatch%TFC_ROS/0.45_r8, currentPatch%ROS_front/60.0_r8)

          ! track fires greater than kW/m energy threshold
          if (currentPatch%FI > SF_val_fire_threshold) then 
            currentPatch%fire = 1 
            currentSite%NF_successful = currentSite%NF_successful + &
              currentSite%NF*currentSite%FDI*currentPatch%area/area
          end if
          
        end if
      end if
      currentPatch => currentPatch%younger
    end do    

  end subroutine CalculateSurfaceFireIntensity
   
  !---------------------------------------------------------------------------------------
  
  subroutine CalculateAreaBurnt(currentSite)
    !
    !  DESCRIPTION:
    !  Calculates area burnt for each patch of a site
    !
    use FatesConstantsMod, only : m2_per_km2
    use SFEquationsMod,    only : FireDuration, LengthToBreadth
    use SFEquationsMod,    only : AreaBurnt, FireSize
    use SFParamsMod,       only : SF_val_fire_threshold

    ! ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite
    
    ! LOCALS:
    type(fates_patch_type), pointer :: currentPatch                    ! patch object
    real(r8)                        :: tree_fraction_patch             ! treed fraction on patch [0-1]
    real(r8)                        :: length_to_breadth               ! length to breadth ratio of fire ellipse (unitless)
    real(r8)                        :: fire_size                       ! size of fire [m2]
    real(r8)                        :: area_burnt                      ! area burnt [m2/km2]
    
    currentPatch => currentSite%oldest_patch 
    do while (associated(currentPatch))

      if (currentPatch%nocomp_pft_label /= nocomp_bareground) then

        ! initialize patch parameters to zero
        currentPatch%FD = 0.0_r8
        currentPatch%frac_burnt = 0.0_r8

        if (currentSite%NF > 0.0_r8 .and. currentPatch%FI > SF_val_fire_threshold) then

          ! fire duration [min]
          currentPatch%FD = FireDuration(currentSite%FDI)
          
          ! length-to-breadth ratio of fire ellipse [unitless]
          tree_fraction_patch  = currentPatch%total_tree_area/currentPatch%area
          length_to_breadth = LengthToBreadth(currentSite%fireWeather%effective_windspeed, tree_fraction_patch)

          ! fire size [m2]
          fire_size = FireSize(length_to_breadth, currentPatch%ROS_back,              &
              currentPatch%ROS_front, currentPatch%FD)

          ! area burnt [m2/km2]
          area_burnt = AreaBurnt(fire_size, currentSite%NF, currentSite%FDI)
          
          ! convert to area burned per area patch per day
          ! i.e., fraction of the patch burned on that day
          currentPatch%frac_burnt = min(0.99_r8, area_burnt/m2_per_km2)
          
        end if
      end if
      currentPatch => currentPatch%younger
    end do    

  end subroutine CalculateAreaBurnt
   
  !---------------------------------------------------------------------------------------
  
  subroutine CalculatePostFireMortality(currentSite)
    !
    !  DESCRIPTION:
    !  Calculates mortality (for woody PFTs) due to fire from crown scorching and cambial damage
    !
    use SFEquationsMod, only : ScorchHeight, CrownFireMortality, CrownFractionBurnt
    use SFEquationsMod, only : CambialMortality, TotalFireMortality
    
    ! ARGUMENTS:
    type(ed_site_type), intent(in), target :: currentSite ! site object
    
    ! LOCALS:
    type(fates_patch_type),  pointer :: currentPatch    ! patch object
    type(fates_cohort_type), pointer :: currentCohort   ! cohort object
    real(r8)                         :: tree_ag_biomass ! total amount of above-ground tree biomass in patch [kgC/m2]
    
    currentPatch => currentSite%oldest_patch
    do while (associated(currentPatch))
      if (currentPatch%nocomp_pft_label /= nocomp_bareground .and. currentPatch%fire == 1) then
        
        ! sum up woody agoveground biomass on patch
        tree_ag_biomass = 0.0_r8
        currentCohort => currentPatch%tallest
        do while (associated(currentCohort))
          if (prt_params%woody(currentCohort%pft) == itrue) then
            leaf_c = currentCohort%prt%GetState(leaf_organ, carbon12_element)
            sapw_c = currentCohort%prt%GetState(sapw_organ, carbon12_element)
            struct_c = currentCohort%prt%GetState(struct_organ, carbon12_element)
            tree_ag_biomass = tree_ag_biomass + currentCohort%n*(leaf_c + prt_params%allom_agb_frac(currentCohort%pft)*(sapw_c + struct_c))
          end if 
          currentCohort => currentCohort%shorter
        end do
        
        ! calculate scorch height [m]
        currentPatch%Scorch_ht(:) = 0.0_r8
        if (tree_ag_biomass > 0.0_r8) then 
          do i_pft = 1, numpft
            if (prt_params%woody(i_pft) == itrue) then 
              currentPatch%Scorch_ht(i_pft) = ScorchHeight(EDPftvarcon_inst%fire_alpha_SH(i_pft), currentPatch%FI)
            end if
          end do
        end if 
        
        ! now calculate actual mortality
        do while (associated(currentCohort))  
          
          currentCohort%fraction_crown_burned = 0.0_r8
          currentCohort%fire_mort = 0.0_r8
          currentCohort%crownfire_mort = 0.0_r8
          currentCohort%cambial_damage_kill = 0.0_r8
          
          if (prt_params%woody(currentCohort%pft) == itrue) then 
            
            call CrownDepth(currentCohort%height, currentCohort%pft, crown_depth)
            
            ! calculate crown fraction burned
            currentCohort%fraction_crown_burned = CrownFractionBurnt(currentPatch%Scorch_ht(currentCohort%pft),  &
              currentCohort%height, crown_depth)
              
            ! calculate cambial mortality            
            currentCohort%cambial_damage_kill = CambialMortality(EDPftvarcon_inst%bark_scaler(currentCohort%pft), & 
              currentCohort%dbh, currentPatch%tau_l)
            
            ! calculate crown fire mortality
            currentCohort%crownfire_mort = CrownFireMortality(EDPftvarcon_inst%crown_kill(currentCohort%pft), &
              currentCohort%fraction_crown_burned)
              
            ! total fire mortality
            currentCohort%fire_mort = TotalFireMortality(currentCohort%crownfire_mort,   &
              currentCohort%cambial_mort)
              
          end if
          currentCohort => currentCohort%shorter
        end do 
      end if
      currentPatch => currentPatch%younger
    end do 
  
  end subroutine CalculatePostFireMortality
  
  !---------------------------------------------------------------------------------------
  
  subroutine crown_scorching(currentSite) 

  type(ed_site_type), intent(in), target :: currentSite

  type(fates_patch_type),  pointer :: currentPatch
  type(fates_cohort_type), pointer :: currentCohort

  real(r8) ::  tree_ag_biomass ! total amount of above-ground tree biomass in patch. kgC/m2
  real(r8) ::  leaf_c          ! leaf carbon      [kg]
  real(r8) ::  sapw_c          ! sapwood carbon   [kg]
  real(r8) ::  struct_c        ! structure carbon [kg]
  integer  ::  i_pft

  currentPatch => currentSite%oldest_patch
  do while (associated(currentPatch)) 

    if (currentPatch%nocomp_pft_label /= nocomp_bareground) then

      tree_ag_biomass = 0.0_r8
      if (currentPatch%fire == 1) then
      
        currentCohort => currentPatch%tallest
        do while (associated(currentCohort))  
          if (prt_params%woody(currentCohort%pft) == itrue) then
            leaf_c = currentCohort%prt%GetState(leaf_organ, carbon12_element)
            sapw_c = currentCohort%prt%GetState(sapw_organ, carbon12_element)
            struct_c = currentCohort%prt%GetState(struct_organ, carbon12_element)
            tree_ag_biomass = tree_ag_biomass + currentCohort%n*(leaf_c + prt_params%allom_agb_frac(currentCohort%pft)*(sapw_c + struct_c))
          end if 
          currentCohort => currentCohort%shorter
        end do 

        do i_pft = 1, numpft
          if (tree_ag_biomass > 0.0_r8 .and. prt_params%woody(i_pft) == itrue) then 
            ! Equation 16 in Thonicke et al. 2010 !Van Wagner 1973 EQ8 !2/3 Byram (1959)
            currentPatch%Scorch_ht(i_pft) = EDPftvarcon_inst%fire_alpha_SH(i_pft)*(currentPatch%FI**0.667_r8)
          else
            currentPatch%Scorch_ht(i_pft) = 0.0_r8
          end if
        end do
      end if 
    end if 
    currentPatch => currentPatch%younger
  end do

  end subroutine crown_scorching
  
  !---------------------------------------------------------------------------------------
  
  subroutine crown_damage(currentSite)

    type(ed_site_type), intent(in), target :: currentSite
    type(fates_patch_type) , pointer :: currentPatch
    type(fates_cohort_type), pointer :: currentCohort
    real(r8)                         :: crown_depth    ! Depth of crown in meters

    currentPatch => currentSite%oldest_patch
    do while (associated(currentPatch)) 

      if (currentPatch%nocomp_pft_label /= nocomp_bareground) then
        if (currentPatch%fire == 1) then
          currentCohort => currentPatch%tallest
          do while (associated(currentCohort))  
            currentCohort%fraction_crown_burned = 0.0_r8
            if (prt_params%woody(currentCohort%pft) == itrue) then 
              ! Flames lower than bottom of canopy. 
              ! c%height is height of cohort
              call CrownDepth(currentCohort%height, currentCohort%pft, crown_depth)

              if (currentPatch%Scorch_ht(currentCohort%pft) < (currentCohort%height-crown_depth)) then 
                currentCohort%fraction_crown_burned = 0.0_r8
              else
                ! Flames part of way up canopy. 
                ! Equation 17 in Thonicke et al. 2010. 
                ! flames over bottom of canopy but not over top.
                if((currentCohort%height > 0.0_r8) .and. (currentPatch%Scorch_ht(currentCohort%pft) >=  &
                  (currentCohort%height-crown_depth))) then 
                  currentCohort%fraction_crown_burned = (currentPatch%Scorch_ht(currentCohort%pft) - &
                    (currentCohort%height - crown_depth))/crown_depth
                else 
                  ! Flames over top of canopy. 
                  currentCohort%fraction_crown_burned =  1.0_r8 
                end if
              endif
              ! Check for strange values. 
              currentCohort%fraction_crown_burned = min(1.0_r8, max(0.0_r8,currentCohort%fraction_crown_burned))              
            end if !trees only
            !shrink canopy to account for burnt section.     
            !currentCohort%canopy_trim = min(currentCohort%canopy_trim,(1.0_r8-currentCohort%fraction_crown_burned))
            currentCohort => currentCohort%shorter
          end do 
        end if 
      end if
      currentPatch => currentPatch%younger
    end do

  end subroutine crown_damage
  
  !---------------------------------------------------------------------------------------

  subroutine cambial_damage_kill(currentSite)
    ! routine description.
    ! returns the probability that trees dies due to cambial char
    ! currentPatch%tau_l = duration of lethal stem heating (min). Calculated at patch level.

    type(ed_site_type), intent(in), target :: currentSite

    type(fates_patch_type) , pointer :: currentPatch
    type(fates_cohort_type), pointer :: currentCohort

    real(r8) :: tau_c !critical time taken to kill cambium (minutes) 
    real(r8) :: bt    !bark thickness in cm.

    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))
      if (currentPatch%nocomp_pft_label /= nocomp_bareground) then
        if (currentPatch%fire == 1) then
          currentCohort => currentPatch%tallest
          do while(associated(currentCohort))  
            if ( prt_params%woody(currentCohort%pft) == itrue) then 
              
              ! Equation 21 in Thonicke et al 2010
              bt = EDPftvarcon_inst%bark_scaler(currentCohort%pft)*currentCohort%dbh 
              ! Equation 20 in Thonicke et al. 2010. 
              tau_c = 2.9_r8*bt**2.0_r8
              ! Equation 19 in Thonicke et al. 2010
              
              if ((currentPatch%tau_l/tau_c) >= 2.0_r8) then
                currentCohort%cambial_mort = 1.0_r8
              else
                if ((currentPatch%tau_l/tau_c) > 0.22_r8) then
                  currentCohort%cambial_mort = (0.563_r8*(currentPatch%tau_l/tau_c)) - 0.125_r8
                else
                  currentCohort%cambial_mort = 0.0_r8
                endif
              end if
            end if 
            currentCohort => currentCohort%shorter
          end do 
        end if 
      end if
      currentPatch => currentPatch%younger
    end do
    
  end subroutine cambial_damage_kill
  
  !---------------------------------------------------------------------------------------

  subroutine post_fire_mortality(currentSite)

    type(ed_site_type), intent(in), target :: currentSite
    type(fates_patch_type),  pointer :: currentPatch
    type(fates_cohort_type), pointer :: currentCohort

    currentPatch => currentSite%oldest_patch
    do while (associated(currentPatch)) 

      if (currentPatch%nocomp_pft_label /= nocomp_bareground) then
        if (currentPatch%fire == 1) then 
          currentCohort => currentPatch%tallest
          do while(associated(currentCohort))  
            
            currentCohort%fire_mort = 0.0_r8
            currentCohort%crownfire_mort = 0.0_r8
            
            if (prt_params%woody(currentCohort%pft) == itrue) then
              ! Equation 22 in Thonicke et al. 2010. 
              currentCohort%crownfire_mort = EDPftvarcon_inst%crown_kill(currentCohort%pft)*currentCohort%fraction_crown_burned**3.0_r8
              ! Equation 18 in Thonicke et al. 2010. 
              currentCohort%fire_mort = max(0._r8,min(1.0_r8,currentCohort%crownfire_mort+currentCohort%cambial_mort- &
                (currentCohort%crownfire_mort*currentCohort%cambial_mort)))  !joint prob.   
              else
                currentCohort%fire_mort = 0.0_r8
            end if 
            currentCohort => currentCohort%shorter
          end do
        end if 
      end if 
      currentPatch => currentPatch%younger
    end do

  end subroutine post_fire_mortality

  ! ======================================================================================
  
end module SFMainMod
