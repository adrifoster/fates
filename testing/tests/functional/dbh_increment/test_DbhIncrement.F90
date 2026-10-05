program FatesDbhIncrement

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesConstantsMod,           only : g_per_kg, sec_per_day, days_per_year
  use FatesConstantsMod,           only : t_water_freeze_k_1atm, umolC_to_kgC
  use FatesConstantsMod,           only : lmrmodel_ryan_1991, fates_check_param_set
  use FatesArgumentUtils,          only : command_line_arg
  use FatesUnitTestParamReaderMod, only : ReadParameters
  use FatesUnitTestParamReaderMod, only : CheckLeafRespParams
  use FatesFactoryMod,             only : InitializeGlobals
  use FatesGlobals,                only : FatesGlobalsInit
  use FatesInterfaceTypesMod,      only : numpft
  use FatesInterfaceTypesMod,      only : hlm_maintresp_leaf_model
  use EDPftvarcon,                 only : EDPftvarcon_inst
  use FatesParameterDerivedMod,    only : param_derived
  use PRTParametersMod,            only : prt_params
  use PRTInitParamsFatesMod,       only : PRTDerivedParams
  use PRTGenericMod,               only : sapw_organ, fnrt_organ
  use FatesAllometryMod,           only : h_allom, bagw_allom, bleaf
  use FatesAllometryMod,           only : carea_allom, bsap_allom, bbgw_allom
  use FatesAllometryMod,           only : bfineroot, bstore_allom, bdead_allom
  use FatesPlantRespPhotosynthMod, only : NonleafMaintenanceRespiration
  use FatesTestLeafPhotoMod,       only : EvaluateLeafPhotosynthesis, LeafNitrogenContent
  use LeafBiophysicsMod,           only : lb_params
  use LeafBiophysicsMod,           only : FvCB1980, medlyn_model, net_assim_model
  use LeafBiophysicsMod,           only : photosynth_acclim_model_kumarathunge_etal_2019
  use LeafBiophysicsMod,           only : QSat
  use FatesTestEnvironmentMod,     only : environment_type
  use FatesTestEnvironmentMod,     only : default_nscaler, default_rdark_scaler
  use FatesTestEnvironmentMod,     only : CanopyVaporPressure

  implicit none

  ! LOCALS:
  character(len=:), allocatable :: param_file    ! input parameter file
  character(len=:), allocatable :: out_file      ! output file
  logical                       :: a_net_prescribed ! a_net read from fates_dev_arbitrary_pft
  type(environment_type)        :: env           ! prescribed atmospheric boundary conditions
  real(r8),         allocatable :: dbh(:)        ! dbh [cm]
  real(r8),         allocatable :: dbh_incr(:)   ! annual dbh increment [cm yr-1]
  real(r8),         allocatable :: c_leaf(:)     ! leaf carbon [kg]
  real(r8),         allocatable :: c_fnrt(:)     ! fineroot carbon [kg]
  real(r8),         allocatable :: c_sapw(:)     ! sapwood carbon [kg]
  real(r8),         allocatable :: c_agw(:)      ! aboveground carbon [kg]
  real(r8),         allocatable :: c_bgw(:)      ! belowground carbon [kg]
  real(r8),         allocatable :: c_struct(:)   ! structural carbon [kg]
  real(r8),         allocatable :: c_store(:)    ! storage carbon [kg]
  real(r8),         allocatable :: dleafdd_yr(:) ! leaf carbon derivative wrt dbh [kgC cm-1]
  real(r8),         allocatable :: dtotaldd_yr(:) ! total target carbon derivative wrt dbh [kgC cm-1]
  real(r8),         allocatable :: dbh_grid(:)      ! fixed dbh grid [cm]
  real(r8),         allocatable :: c_leaf_grid(:)   ! leaf carbon on dbh grid [kgC]
  real(r8),         allocatable :: dleafdd_grid(:)  ! leaf carbon derivative wrt dbh on dbh grid [kgC cm-1]
  real(r8),         allocatable :: dtotaldd_grid(:) ! total target carbon derivative wrt dbh on dbh grid [kgC cm-1]
  real(r8),         allocatable :: c_fnrt_grid(:)   ! fineroot carbon on dbh grid [kgC]
  real(r8),         allocatable :: c_sapw_grid(:)   ! sapwood carbon on dbh grid [kgC]
  real(r8),         allocatable :: c_agw_grid(:)    ! aboveground woody carbon on dbh grid [kgC]
  real(r8),         allocatable :: c_bgw_grid(:)    ! belowground woody carbon on dbh grid [kgC]
  real(r8),         allocatable :: c_store_grid(:)  ! storage carbon on dbh grid [kgC]
  real(r8),         allocatable :: c_struct_grid(:) ! structural carbon on dbh grid [kgC]
  real(r8),         allocatable :: height_grid(:)   ! plant height on dbh grid [m]
  real(r8),         allocatable :: a_net(:)      ! net assimilation net of leaf dark respiration [kgC m-2 leaf day-1]
  real(r8)                      :: dbh_now       ! working dbh [cm]
  real(r8)                      :: dbh_start     ! dbh at start of year [cm]
  real(r8)                      :: c_leaf_now    ! working leaf carbon [kg]
  real(r8)                      :: c_fnrt_now    ! working fineroot carbon [kg]
  real(r8)                      :: c_sapw_now    ! working sapwood carbon [kg]
  real(r8)                      :: c_agw_now     ! working aboveground carbon [kg]
  real(r8)                      :: c_bgw_now     ! working belowground carbon [kg]
  real(r8)                      :: c_struct_now  ! working structural carbon [kg]
  real(r8)                      :: c_store_now   ! working storage carbon [kg]
  real(r8)                      :: sapw_area     ! sapwood area [m2]
  real(r8)                      :: height_now    ! working height [m]
  real(r8)                      :: l2fr          ! leaf to fineroot ratio
  real(r8)                      :: dleafdd       ! leaf derivative
  real(r8)                      :: dfnrtdd       ! fineroot derivative
  real(r8)                      :: dsapwdd       ! sapwood derivative
  real(r8)                      :: dagwdd        ! aboveground derivative
  real(r8)                      :: dbgwdd        ! belowground derivative
  real(r8)                      :: dstructdd     ! structural derivative
  real(r8)                      :: dstoredd      ! storage derivative
  real(r8)                      :: dtotaldd      ! total derivative
  real(r8)                      :: lnc_top       ! leaf N content at the canopy top [gN m-2 leaf]
  real(r8)                      :: vcmax25top    ! canopy-top vcmax at 25degC [umol m-2 s-1]
  real(r8)                      :: jmax25top     ! canopy-top jmax at 25degC [umol m-2 s-1]
  real(r8)                      :: kp25top       ! canopy-top C4 initial slope at 25degC [umol m-2 s-1]
  real(r8)                      :: veg_esat      ! saturation vapor pressure at t_ref [Pa]
  real(r8)                      :: can_vpress    ! canopy air vapor pressure [Pa]
  real(r8)                      :: qs_dummy      ! saturation specific humidity (unused)
  real(r8)                      :: agross_light  ! gross photosynthesis, light [umolC m-2 s-1]
  real(r8)                      :: anet_light    ! net photosynthesis, light [umolC m-2 s-1]
  real(r8)                      :: agross_dark   ! gross photosynthesis, dark [umolC m-2 s-1]
  real(r8)                      :: anet_dark     ! net photosynthesis, dark [umolC m-2 s-1]
  real(r8)                      :: gs            ! stomatal conductance (unused) [umol H2O m-2 s-1]
  real(r8)                      :: ci            ! intracellular CO2 (unused) [Pa]
  real(r8)                      :: a_net_now     ! daily net assimilation net of leaf dark respiration [kgC m-2 leaf day-1]
  real(r8)                      :: leaf_area     ! leaf area [m2]
  real(r8)                      :: c_gross       ! leaf net assimilation [kgC day-1]
  real(r8)                      :: live_stem_n   ! aboveground sapwood N [kgN]
  real(r8)                      :: live_croot_n  ! belowground sapwood N [kgN]
  real(r8)                      :: fnrt_n        ! fineroot N [kgN]
  real(r8)                      :: livestem_mr   ! live stem MR [kgC s-1]
  real(r8)                      :: livecroot_mr  ! live coarse root MR [kgC s-1]
  real(r8)                      :: froot_mr      ! fineroot MR [kgC s-1]
  real(r8)                      :: sym_nfix      ! symbiotic N fixation [kgN day-1]
  real(r8)                      :: resp_m        ! non-leaf maintenance respiration [kgC day-1]
  real(r8)                      :: resp_g        ! growth respiration [kgC day-1]
  real(r8)                      :: npp           ! net primary production [kgC day-1]
  real(r8)                      :: turnover      ! turnover replacement [kgC day-1]
  real(r8)                      :: c_growth      ! carbon available for growth [kgC day-1]
  real(r8)                      :: c_growth_ann  ! annual carbon available for growth [kgC yr-1]
  real(r8)                      :: repro_fraction ! fraction of growth carbon to reproduction [0-1]
  integer                       :: iyr           ! year index
  integer                       :: iday          ! day index
  integer                       :: igrid         ! dbh grid index

  ! CONSTANTS:
  real(r8),         parameter :: init_dbh = 2.5_r8            ! initial DBH [cm]
  real(r8),         parameter :: canopy_trim = 1.0_r8         ! no trimming
  real(r8),         parameter :: elongf_leaf = 1.0_r8         ! fully elongated leaves
  real(r8),         parameter :: elongf_fnrt = 1.0_r8         ! fully elongated roots
  real(r8),         parameter :: elongf_stem = 1.0_r8         ! fully elongated stem
  real(r8),         parameter :: par_abs = 2500.0_r8          ! absorbed PAR, light period [umol photons m-2 leaf s-1]
  real(r8),         parameter :: light_seconds = 43200.0_r8   ! length of light period [s day-1]
  real(r8),         parameter :: t_ref = 25.0_r8 + t_water_freeze_k_1atm ! vegetation and soil temperature [K]
  real(r8),         parameter :: maintresp_reduction = 1.0_r8 ! no MR reduction (storage at target)
  integer,          parameter :: crown_damage = 1             ! undamaged
  integer,          parameter :: pft = 6                      ! plant functional type to simulate
  integer,          parameter :: num_years = 400              ! number of years to simulate
  integer,          parameter :: num_days = 365               ! days per year
  integer,          parameter :: num_grid = 200               ! number of dbh grid points
  real(r8),         parameter :: dbh_grid_min = 2.5_r8        ! smallest grid dbh [cm]
  real(r8),         parameter :: dbh_grid_max = 150.0_r8      ! largest grid dbh [cm]

  interface

    subroutine WriteDbhIncrementData(out_file, num_years, dbh, dbh_incr, c_leaf, c_fnrt, &
        c_sapw, c_agw, c_bgw, c_struct, c_store, a_net, dleafdd, dtotaldd, num_grid,     &
        dbh_grid, c_leaf_grid, dleafdd_grid, dtotaldd_grid, c_fnrt_grid, c_sapw_grid,     &
        c_agw_grid, c_bgw_grid, c_store_grid, c_struct_grid, height_grid)

      use FatesConstantsMod, only : r8 => fates_r8
      implicit none

      character(len=*), intent(in) :: out_file
      integer,          intent(in) :: num_years
      real(r8),         intent(in) :: dbh(:)
      real(r8),         intent(in) :: dbh_incr(:)
      real(r8),         intent(in) :: c_leaf(:)
      real(r8),         intent(in) :: c_fnrt(:)
      real(r8),         intent(in) :: c_sapw(:)
      real(r8),         intent(in) :: c_agw(:)
      real(r8),         intent(in) :: c_bgw(:)
      real(r8),         intent(in) :: c_struct(:)
      real(r8),         intent(in) :: c_store(:)
      real(r8),         intent(in) :: a_net(:)
      real(r8),         intent(in) :: dleafdd(:)
      real(r8),         intent(in) :: dtotaldd(:)
      integer,          intent(in) :: num_grid
      real(r8),         intent(in) :: dbh_grid(:)
      real(r8),         intent(in) :: c_leaf_grid(:)
      real(r8),         intent(in) :: dleafdd_grid(:)
      real(r8),         intent(in) :: dtotaldd_grid(:)
      real(r8),         intent(in) :: c_fnrt_grid(:)
      real(r8),         intent(in) :: c_sapw_grid(:)
      real(r8),         intent(in) :: c_agw_grid(:)
      real(r8),         intent(in) :: c_bgw_grid(:)
      real(r8),         intent(in) :: c_store_grid(:)
      real(r8),         intent(in) :: c_struct_grid(:)
      real(r8),         intent(in) :: height_grid(:)
    end subroutine WriteDbhIncrementData

  end interface

  ! read in parameter file name from command line
  param_file = command_line_arg(1)

  ! optional 2nd command-line argument: output file
  if (command_argument_count() >= 2) then
    out_file = trim(command_line_arg(2))
  else
    out_file = 'dbh_incr_out.nc'
  end if

  ! read in parameter file
  call ReadParameters(param_file)

  call InitializeGlobals(sec_per_day)
  numpft = size(prt_params%wood_density, dim=1)
  call FatesGlobalsInit(6, .false.)
  call PRTDerivedParams()

  ! leaf biophysics switches
  hlm_maintresp_leaf_model = lmrmodel_ryan_1991
  lb_params%electron_transport_model = FvCB1980
  lb_params%stomatal_model = medlyn_model
  lb_params%stomatal_assim_model = net_assim_model
  lb_params%photo_tempsens_model = photosynth_acclim_model_kumarathunge_etal_2019
  call CheckLeafRespParams()

  ! canopy-top photosynthetic capacity
  lnc_top = LeafNitrogenContent(pft)
  vcmax25top = EDPftvarcon_inst%vcmax25top(pft, 1)
  jmax25top = param_derived%jmax25top(pft, 1)
  kp25top = param_derived%kp25top(pft, 1)

  call env%Init(tempk=t_ref)

  ! allocate arrays
  allocate(dbh(num_years))
  allocate(dbh_incr(num_years))
  allocate(c_leaf(num_years))
  allocate(c_fnrt(num_years))
  allocate(c_sapw(num_years))
  allocate(c_agw(num_years))
  allocate(c_bgw(num_years))
  allocate(c_struct(num_years))
  allocate(c_store(num_years))
  allocate(a_net(num_years))
  allocate(dleafdd_yr(num_years))
  allocate(dtotaldd_yr(num_years))
  allocate(dbh_grid(num_grid))
  allocate(c_leaf_grid(num_grid))
  allocate(dleafdd_grid(num_grid))
  allocate(dtotaldd_grid(num_grid))
  allocate(c_fnrt_grid(num_grid))
  allocate(c_sapw_grid(num_grid))
  allocate(c_agw_grid(num_grid))
  allocate(c_bgw_grid(num_grid))
  allocate(c_store_grid(num_grid))
  allocate(c_struct_grid(num_grid))
  allocate(height_grid(num_grid))

  l2fr = prt_params%allom_l2fr(pft)
  
  ! daily net assimilation per unit leaf area, prescribed if fates_dev_arbitrary_pft is set [kgC m-2 leaf day-1]
  a_net_prescribed = (EDPftvarcon_inst%dev_arbitrary_pft(pft) < fates_check_param_set)
  if (a_net_prescribed) then
    a_net_now = EDPftvarcon_inst%dev_arbitrary_pft(pft)
  else
    call QSat(env%tempk, env%can_press, qs_dummy, veg_esat)
    can_vpress = CanopyVaporPressure(veg_esat)

    call EvaluateLeafPhotosynthesis(pft, par_abs, env%tempk, env%tempk, env%tempk,   &
      env%can_press, env%can_co2_ppress, env%can_o2_ppress, veg_esat, can_vpress,    &
      env%gb, default_nscaler, default_rdark_scaler, env%dayl_factor, env%btran,     &
      vcmax25top, jmax25top, kp25top, lnc_top, agross_light, anet_light, gs, ci)

    call EvaluateLeafPhotosynthesis(pft, 0.0_r8, env%tempk, env%tempk, env%tempk,    &
      env%can_press, env%can_co2_ppress, env%can_o2_ppress, veg_esat, can_vpress,    &
      env%gb, default_nscaler, default_rdark_scaler, env%dayl_factor, env%btran,     &
      vcmax25top, jmax25top, kp25top, lnc_top, agross_dark, anet_dark, gs, ci)

    a_net_now = (anet_light*light_seconds + anet_dark*(sec_per_day - light_seconds))*umolC_to_kgC
  end if

  dbh_now = init_dbh
  do iyr = 1, num_years
    dbh_start = dbh_now
    c_growth_ann = 0.0_r8
    do iday = 1, num_days
      call bleaf(dbh_now, pft, crown_damage, canopy_trim, elongf_leaf, c_leaf_now, dbldd=dleafdd)
      call bfineroot(dbh_now, pft, canopy_trim, l2fr, elongf_fnrt, c_fnrt_now, dfnrtdd)
      call bsap_allom(dbh_now, pft, crown_damage, canopy_trim, elongf_stem, sapw_area, c_sapw_now, dsapwdd)
      call bagw_allom(dbh_now, pft, crown_damage, elongf_stem, c_agw_now, dagwdd)
      call bbgw_allom(dbh_now, pft, elongf_stem, c_bgw_now, dbgwdd)
      call bdead_allom(c_agw_now, c_bgw_now, c_sapw_now, pft, c_struct_now, dagwdd, dbgwdd, dsapwdd, dstructdd)
      call bstore_allom(dbh_now, pft, crown_damage, canopy_trim, c_store_now, dstoredd)

      dtotaldd = dleafdd + dfnrtdd + dsapwdd + dstructdd + dstoredd

      ! leaf net assimilation
      leaf_area = c_leaf_now*prt_params%slatop(pft)*g_per_kg
      c_gross = a_net_now*leaf_area

      ! non-leaf maintenance respiration
      live_stem_n = prt_params%allom_agb_frac(pft)*c_sapw_now*                             &
        prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(sapw_organ))
      live_croot_n = (1.0_r8 - prt_params%allom_agb_frac(pft))*c_sapw_now*                 &
        prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(sapw_organ))
      fnrt_n = c_fnrt_now*prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(fnrt_organ))

      call NonleafMaintenanceRespiration(pft, t_ref, 1, [t_ref], [1.0_r8], live_stem_n,    &
        live_croot_n, fnrt_n, maintresp_reduction, sec_per_day, livestem_mr, livecroot_mr, &
        froot_mr, sym_nfix)

      resp_m = (livestem_mr + livecroot_mr + froot_mr)*sec_per_day

      ! growth respiration
      resp_g = prt_params%grperc(pft)*max(0.0_r8, c_gross - resp_m)

      npp = c_gross - resp_m - resp_g

      ! turnover replacement
      turnover = c_leaf_now/(prt_params%leaf_long(pft, size(prt_params%leaf_long, dim=2))*days_per_year) + &
        c_fnrt_now/(prt_params%root_long(pft)*days_per_year) +                             &
        (c_sapw_now + c_struct_now + c_store_now)/(prt_params%branch_long(pft)*days_per_year)

      c_growth = npp - turnover

      ! reproductive allocation
      if (dbh_now <= prt_params%dbh_repro_threshold(pft)) then
        repro_fraction = prt_params%seed_alloc(pft)
      else
        repro_fraction = prt_params%seed_alloc(pft) + prt_params%seed_alloc_mature(pft)
      end if

      if (c_growth > 0.0_r8) then
        dbh_now = dbh_now + (1.0_r8 - repro_fraction)*c_growth/dtotaldd
      end if
      c_growth_ann = c_growth_ann + (1.0_r8 - repro_fraction)*c_growth
    end do

    ! end-of-year allometry
    call bleaf(dbh_now, pft, crown_damage, canopy_trim, elongf_leaf, c_leaf_now, dbldd=dleafdd)
    call bfineroot(dbh_now, pft, canopy_trim, l2fr, elongf_fnrt, c_fnrt_now, dfnrtdd)
    call bsap_allom(dbh_now, pft, crown_damage, canopy_trim, elongf_stem, sapw_area, c_sapw_now, dsapwdd)
    call bagw_allom(dbh_now, pft, crown_damage, elongf_stem, c_agw_now, dagwdd)
    call bbgw_allom(dbh_now, pft, elongf_stem, c_bgw_now, dbgwdd)
    call bdead_allom(c_agw_now, c_bgw_now, c_sapw_now, pft, c_struct_now, dagwdd, dbgwdd, dsapwdd, dstructdd)
    call bstore_allom(dbh_now, pft, crown_damage, canopy_trim, c_store_now, dstoredd)

    dtotaldd = dleafdd + dfnrtdd + dsapwdd + dstructdd + dstoredd

    dbh(iyr) = dbh_now
    dbh_incr(iyr) = dbh_now - dbh_start
    c_leaf(iyr) = c_leaf_now
    c_fnrt(iyr) = c_fnrt_now
    c_sapw(iyr) = c_sapw_now
    c_agw(iyr) = c_agw_now
    c_bgw(iyr) = c_bgw_now
    c_struct(iyr) = c_struct_now
    c_store(iyr) = c_store_now
    a_net(iyr) = c_growth_ann
    dleafdd_yr(iyr) = dleafdd
    dtotaldd_yr(iyr) = dtotaldd
  end do

  ! allometry on a fixed log-spaced dbh grid
  do igrid = 1, num_grid
    dbh_grid(igrid) = exp(log(dbh_grid_min) + real(igrid - 1, r8)*                     &
      (log(dbh_grid_max) - log(dbh_grid_min))/real(num_grid - 1, r8))

    call bleaf(dbh_grid(igrid), pft, crown_damage, canopy_trim, elongf_leaf, c_leaf_now, dbldd=dleafdd)
    call bfineroot(dbh_grid(igrid), pft, canopy_trim, l2fr, elongf_fnrt, c_fnrt_now, dfnrtdd)
    call bsap_allom(dbh_grid(igrid), pft, crown_damage, canopy_trim, elongf_stem, sapw_area, c_sapw_now, dsapwdd)
    call bagw_allom(dbh_grid(igrid), pft, crown_damage, elongf_stem, c_agw_now, dagwdd)
    call bbgw_allom(dbh_grid(igrid), pft, elongf_stem, c_bgw_now, dbgwdd)
    call bdead_allom(c_agw_now, c_bgw_now, c_sapw_now, pft, c_struct_now, dagwdd, dbgwdd, dsapwdd, dstructdd)
    call bstore_allom(dbh_grid(igrid), pft, crown_damage, canopy_trim, c_store_now, dstoredd)
    call h_allom(dbh_grid(igrid), pft, height_now)

    dtotaldd = dleafdd + dfnrtdd + dsapwdd + dstructdd + dstoredd

    c_leaf_grid(igrid) = c_leaf_now
    dleafdd_grid(igrid) = dleafdd
    dtotaldd_grid(igrid) = dtotaldd
    c_fnrt_grid(igrid) = c_fnrt_now
    c_sapw_grid(igrid) = c_sapw_now
    c_agw_grid(igrid) = c_agw_now
    c_bgw_grid(igrid) = c_bgw_now
    c_store_grid(igrid) = c_store_now
    c_struct_grid(igrid) = c_struct_now
    height_grid(igrid) = height_now
  end do

  ! write out data to netcdf file
  call WriteDbhIncrementData(out_file, num_years, dbh, dbh_incr, c_leaf, c_fnrt, c_sapw, &
    c_agw, c_bgw, c_struct, c_store, a_net, dleafdd_yr, dtotaldd_yr, num_grid, dbh_grid, &
    c_leaf_grid, dleafdd_grid, dtotaldd_grid, c_fnrt_grid, c_sapw_grid, c_agw_grid,       &
    c_bgw_grid, c_store_grid, c_struct_grid, height_grid)

  ! deallocate arrays
  if (allocated(dbh)) deallocate(dbh)
  if (allocated(dbh_incr)) deallocate(dbh_incr)
  if (allocated(c_leaf)) deallocate(c_leaf)
  if (allocated(c_fnrt)) deallocate(c_fnrt)
  if (allocated(c_sapw)) deallocate(c_sapw)
  if (allocated(c_agw)) deallocate(c_agw)
  if (allocated(c_bgw)) deallocate(c_bgw)
  if (allocated(c_struct)) deallocate(c_struct)
  if (allocated(c_store)) deallocate(c_store)
  if (allocated(a_net)) deallocate(a_net)
  if (allocated(dleafdd_yr)) deallocate(dleafdd_yr)
  if (allocated(dtotaldd_yr)) deallocate(dtotaldd_yr)
  if (allocated(dbh_grid)) deallocate(dbh_grid)
  if (allocated(c_leaf_grid)) deallocate(c_leaf_grid)
  if (allocated(dleafdd_grid)) deallocate(dleafdd_grid)
  if (allocated(dtotaldd_grid)) deallocate(dtotaldd_grid)
  if (allocated(c_fnrt_grid)) deallocate(c_fnrt_grid)
  if (allocated(c_sapw_grid)) deallocate(c_sapw_grid)
  if (allocated(c_agw_grid)) deallocate(c_agw_grid)
  if (allocated(c_bgw_grid)) deallocate(c_bgw_grid)
  if (allocated(c_store_grid)) deallocate(c_store_grid)
  if (allocated(c_struct_grid)) deallocate(c_struct_grid)
  if (allocated(height_grid)) deallocate(height_grid)

end program FatesDbhIncrement

! ----------------------------------------------------------------------------------------

subroutine WriteDbhIncrementData(out_file, num_years, dbh, dbh_incr, c_leaf, c_fnrt,     &
  c_sapw, c_agw, c_bgw, c_struct, c_store, a_net, dleafdd, dtotaldd, num_grid, dbh_grid,  &
  c_leaf_grid, dleafdd_grid, dtotaldd_grid, c_fnrt_grid, c_sapw_grid, c_agw_grid,         &
  c_bgw_grid, c_store_grid, c_struct_grid, height_grid)
  !
  ! DESCRIPTION:
  ! Writes out data from the dbh increment test
  !
  use FatesConstantsMod,  only : r8 => fates_r8
  use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod, only : WriteVar
  use FatesUnitTestIOMod, only : RegisterVarAtts
  use FatesUnitTestIOMod, only : EndNCDef
  use FatesUnitTestIOMod, only : type_double, type_int

  implicit none

  ! ARGUMENTS:
  character(len=*), intent(in) :: out_file    ! output file name
  integer,          intent(in) :: num_years   ! number of years
  real(r8),         intent(in) :: dbh(:)      ! dbh [cm]
  real(r8),         intent(in) :: dbh_incr(:) ! annual dbh increment [cm yr-1]
  real(r8),         intent(in) :: c_leaf(:)   ! leaf carbon [kgC]
  real(r8),         intent(in) :: c_fnrt(:)   ! fineroot carbon [kgC]
  real(r8),         intent(in) :: c_sapw(:)   ! sapwood carbon [kgC]
  real(r8),         intent(in) :: c_agw(:)    ! aboveground woody carbon [kgC]
  real(r8),         intent(in) :: c_bgw(:)    ! belowground woody carbon [kgC]
  real(r8),         intent(in) :: c_struct(:) ! structural carbon [kgC]
  real(r8),         intent(in) :: c_store(:)  ! storage carbon [kgC]
  real(r8),         intent(in) :: a_net(:)    ! net assimilation net of leaf dark respiration [kgC m-2 leaf day-1]
  real(r8),         intent(in) :: dleafdd(:)  ! leaf carbon derivative wrt dbh [kgC cm-1]
  real(r8),         intent(in) :: dtotaldd(:) ! total target carbon derivative wrt dbh [kgC cm-1]
  integer,          intent(in) :: num_grid         ! number of dbh grid points
  real(r8),         intent(in) :: dbh_grid(:)      ! fixed dbh grid [cm]
  real(r8),         intent(in) :: c_leaf_grid(:)   ! leaf carbon on dbh grid [kgC]
  real(r8),         intent(in) :: dleafdd_grid(:)  ! leaf carbon derivative wrt dbh on dbh grid [kgC cm-1]
  real(r8),         intent(in) :: dtotaldd_grid(:) ! total target carbon derivative wrt dbh on dbh grid [kgC cm-1]
  real(r8),         intent(in) :: c_fnrt_grid(:)   ! fineroot carbon on dbh grid [kgC]
  real(r8),         intent(in) :: c_sapw_grid(:)   ! sapwood carbon on dbh grid [kgC]
  real(r8),         intent(in) :: c_agw_grid(:)    ! aboveground woody carbon on dbh grid [kgC]
  real(r8),         intent(in) :: c_bgw_grid(:)    ! belowground woody carbon on dbh grid [kgC]
  real(r8),         intent(in) :: c_store_grid(:)  ! storage carbon on dbh grid [kgC]
  real(r8),         intent(in) :: c_struct_grid(:) ! structural carbon on dbh grid [kgC]
  real(r8),         intent(in) :: height_grid(:)   ! plant height on dbh grid [m]

  ! LOCALS:
  integer, allocatable :: years(:)     ! array of years to write out
  integer              :: i            ! looping index
  integer              :: ncid         ! netcdf file id
  character(len=8)     :: dim_names(2) ! dimension names
  integer              :: dimIDs(2)    ! dimension IDs
  integer              :: yearID
  integer              :: dbhID, dbhincrID
  integer              :: leafID, fnrtID
  integer              :: sapwID, agwID
  integer              :: bgwID, structID
  integer              :: storeID, anetID
  integer              :: dleafddID, dtotalddID
  integer              :: dbhgridID, leafgridID
  integer              :: dleafddgridID, dtotalddgridID
  integer              :: fnrtgridID, sapwgridID
  integer              :: agwgridID, bgwgridID
  integer              :: storegridID, structgridID
  integer              :: heightgridID

  ! create years
  allocate(years(num_years))
  do i = 1, num_years
    years(i) = i
  end do

  ! dimension names
  dim_names = [character(len=8) :: 'year', 'dbh_grid']

  ! open file
  call OpenNCFile(trim(out_file), ncid, 'readwrite')

  ! register dimensions
  call RegisterNCDims(ncid, dim_names, (/num_years, num_grid/), 2, dimIDs)

  ! register year
  call RegisterVarAtts(ncid, dim_names(1), dimIDs(1:1), type_int, 'yr',                 &
    'simulation year', yearID)

  ! register dbh
  call RegisterVarAtts(ncid, 'dbh', dimIDs(1:1), type_double, 'cm',                     &
    'diameter at breast height at end of year', dbhID)

  ! register dbh increment
  call RegisterVarAtts(ncid, 'dbh_incr', dimIDs(1:1), type_double, 'cm yr-1',           &
    'annual diameter increment', dbhincrID)

  ! register leaf carbon
  call RegisterVarAtts(ncid, 'c_leaf', dimIDs(1:1), type_double, 'kgC',                 &
    'leaf carbon', leafID)

  ! register fineroot carbon
  call RegisterVarAtts(ncid, 'c_fnrt', dimIDs(1:1), type_double, 'kgC',                 &
    'fineroot carbon', fnrtID)

  ! register sapwood carbon
  call RegisterVarAtts(ncid, 'c_sapw', dimIDs(1:1), type_double, 'kgC',                 &
    'sapwood carbon', sapwID)

  ! register aboveground woody carbon
  call RegisterVarAtts(ncid, 'c_agw', dimIDs(1:1), type_double, 'kgC',                  &
    'aboveground woody carbon', agwID)

  ! register belowground woody carbon
  call RegisterVarAtts(ncid, 'c_bgw', dimIDs(1:1), type_double, 'kgC',                  &
    'belowground woody carbon', bgwID)

  ! register structural carbon
  call RegisterVarAtts(ncid, 'c_struct', dimIDs(1:1), type_double, 'kgC',               &
    'structural carbon', structID)

  ! register storage carbon
  call RegisterVarAtts(ncid, 'c_store', dimIDs(1:1), type_double, 'kgC',                &
    'storage carbon', storeID)

  ! register net assimilation
  call RegisterVarAtts(ncid, 'a_net', dimIDs(1:1), type_double, 'kgC m-2 day-1',        &
    'daily net assimilation per unit leaf area, net of leaf dark respiration', anetID)

  ! register leaf carbon derivative
  call RegisterVarAtts(ncid, 'dleafdd', dimIDs(1:1), type_double, 'kgC cm-1',           &
    'leaf carbon derivative wrt dbh', dleafddID)

  ! register total target carbon derivative
  call RegisterVarAtts(ncid, 'dtotaldd', dimIDs(1:1), type_double, 'kgC cm-1',          &
    'total target carbon derivative wrt dbh', dtotalddID)

  ! register dbh grid
  call RegisterVarAtts(ncid, dim_names(2), dimIDs(2:2), type_double, 'cm',              &
    'fixed dbh grid', dbhgridID)

  ! register leaf carbon on dbh grid
  call RegisterVarAtts(ncid, 'c_leaf_grid', dimIDs(2:2), type_double, 'kgC',            &
    'leaf carbon on dbh grid', leafgridID)

  ! register leaf carbon derivative on dbh grid
  call RegisterVarAtts(ncid, 'dleafdd_grid', dimIDs(2:2), type_double, 'kgC cm-1',      &
    'leaf carbon derivative wrt dbh on dbh grid', dleafddgridID)

  ! register total target carbon derivative on dbh grid
  call RegisterVarAtts(ncid, 'dtotaldd_grid', dimIDs(2:2), type_double, 'kgC cm-1',     &
    'total target carbon derivative wrt dbh on dbh grid', dtotalddgridID)

  ! register fineroot carbon on dbh grid
  call RegisterVarAtts(ncid, 'c_fnrt_grid', dimIDs(2:2), type_double, 'kgC',            &
    'fineroot carbon on dbh grid', fnrtgridID)

  ! register sapwood carbon on dbh grid
  call RegisterVarAtts(ncid, 'c_sapw_grid', dimIDs(2:2), type_double, 'kgC',            &
    'sapwood carbon on dbh grid', sapwgridID)

  ! register aboveground woody carbon on dbh grid
  call RegisterVarAtts(ncid, 'c_agw_grid', dimIDs(2:2), type_double, 'kgC',             &
    'aboveground woody carbon on dbh grid', agwgridID)

  ! register belowground woody carbon on dbh grid
  call RegisterVarAtts(ncid, 'c_bgw_grid', dimIDs(2:2), type_double, 'kgC',             &
    'belowground woody carbon on dbh grid', bgwgridID)

  ! register storage carbon on dbh grid
  call RegisterVarAtts(ncid, 'c_store_grid', dimIDs(2:2), type_double, 'kgC',           &
    'storage carbon on dbh grid', storegridID)

  ! register structural carbon on dbh grid
  call RegisterVarAtts(ncid, 'c_struct_grid', dimIDs(2:2), type_double, 'kgC',          &
    'structural carbon on dbh grid', structgridID)

  ! register plant height on dbh grid
  call RegisterVarAtts(ncid, 'height_grid', dimIDs(2:2), type_double, 'm',              &
    'plant height on dbh grid', heightgridID)

  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, yearID, years(:))
  call WriteVar(ncid, dbhID, dbh(:))
  call WriteVar(ncid, dbhincrID, dbh_incr(:))
  call WriteVar(ncid, leafID, c_leaf(:))
  call WriteVar(ncid, fnrtID, c_fnrt(:))
  call WriteVar(ncid, sapwID, c_sapw(:))
  call WriteVar(ncid, agwID, c_agw(:))
  call WriteVar(ncid, bgwID, c_bgw(:))
  call WriteVar(ncid, structID, c_struct(:))
  call WriteVar(ncid, storeID, c_store(:))
  call WriteVar(ncid, anetID, a_net(:))
  call WriteVar(ncid, dleafddID, dleafdd(:))
  call WriteVar(ncid, dtotalddID, dtotaldd(:))
  call WriteVar(ncid, dbhgridID, dbh_grid(:))
  call WriteVar(ncid, leafgridID, c_leaf_grid(:))
  call WriteVar(ncid, dleafddgridID, dleafdd_grid(:))
  call WriteVar(ncid, dtotalddgridID, dtotaldd_grid(:))
  call WriteVar(ncid, fnrtgridID, c_fnrt_grid(:))
  call WriteVar(ncid, sapwgridID, c_sapw_grid(:))
  call WriteVar(ncid, agwgridID, c_agw_grid(:))
  call WriteVar(ncid, bgwgridID, c_bgw_grid(:))
  call WriteVar(ncid, storegridID, c_store_grid(:))
  call WriteVar(ncid, structgridID, c_struct_grid(:))
  call WriteVar(ncid, heightgridID, height_grid(:))

  ! close file
  call CloseNCFile(ncid)

end subroutine WriteDbhIncrementData
