module FatesTestHistoryMod
  !
  ! DESCRIPTION:
  ! Output accumulation and netCDF writing for standalone, patch-less/site-less
  ! cohort test drivers that sweep a set of independent trajectories (e.g. a light
  ! level sweep). Two groups of variables are accumulated in memory (via RecordDay/
  ! RecordLightProfile) and written once at the end (via Write):
  !   - Daily whole-cohort time series, dimensioned (time, light_level), where time
  !     is the day index within each light level's independent trajectory: dbh,
  !     height, treelai, crown area, the PARTEH carbon pools (leaf/fine root/
  !     sapwood/structure/storage), daily net carbon, daily GPP, the carbon balance
  !     handed to PARTEH (npp_acc, net of growth respiration), frac_store (storage
  !     as a fraction of target leaf carbon), the carbon starvation mortality rate
  !     (cmort, indiv/year), and cohort number density (n).
  !   - A per-leaf-layer light profile (parsun_z/parsha_z/laisun_z/laisha_z),
  !     dimensioned (nlevleaf, year, light_level), captured once per year. nlevleaf
  !     is a compile-time maximum leaf+stem layer count; layers above a given
  !     year's cohort%nv are left at the fates_unset_r8 fill value (registered as
  !     each variable's _FillValue attribute), so the array survives nv changing
  !     over time and differing across light levels.
  !

  use FatesConstantsMod,  only : r8 => fates_r8
  use FatesConstantsMod,  only : fates_unset_r8
  use FatesCohortMod,     only : fates_cohort_type
  use PRTGenericMod,      only : leaf_organ, fnrt_organ, sapw_organ, struct_organ, store_organ
  use PRTGenericMod,      only : carbon12_element
  use FatesTestLightEnvMod, only : light_env_type
  use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod, only : WriteVar, RegisterVar, RegisterFillValue, EndNCDef
  use FatesUnitTestIOMod, only : type_double, type_int

  implicit none
  private

  type, public :: history_type

     private

     integer :: n_time  ! days per light level's trajectory (n_days_total)
     integer :: n_leaf   ! compile-time max leaf+stem layers (nlevleaf)
     integer :: n_light  ! number of light levels
     integer :: n_year   ! number of simulated years

     ! daily whole-cohort time series, dimensioned (time, light_level)
     real(r8), allocatable :: dbh(:,:)         ! dbh [cm]
     real(r8), allocatable :: height(:,:)      ! height [m]
     real(r8), allocatable :: treelai(:,:)     ! total leaf area index [m2/m2]
     real(r8), allocatable :: crown_area(:,:)  ! crown area [m2]
     real(r8), allocatable :: leaf_c(:,:)      ! leaf carbon [kgC/indiv]
     real(r8), allocatable :: fnrt_c(:,:)      ! fine root carbon [kgC/indiv]
     real(r8), allocatable :: sapw_c(:,:)      ! sapwood carbon [kgC/indiv]
     real(r8), allocatable :: struct_c(:,:)    ! structural carbon [kgC/indiv]
     real(r8), allocatable :: storage_c(:,:)   ! storage carbon [kgC/indiv]
     real(r8), allocatable :: daily_net_c(:,:) ! daily net C (GPP - leaf Rd - nonleaf MR) [kgC/indiv/day]
     real(r8), allocatable :: daily_gpp(:,:)   ! daily GPP [kgC/indiv/day]
     real(r8), allocatable :: npp_acc(:,:)     ! carbon balance handed to PARTEH (npp_acc, net of growth respiration) [kgC/indiv/day]
     real(r8), allocatable :: frac_store(:,:)  ! storage carbon as a fraction of target leaf carbon [-]
     real(r8), allocatable :: cmort(:,:)       ! carbon starvation mortality rate [indiv/year]
     real(r8), allocatable :: n(:,:)           ! cohort number density, surviving fraction of the original recruitment cohort [indiv]

     ! per-leaf-layer light profile, dimensioned (nlevleaf, year, light_level)
     real(r8), allocatable :: parsun_z(:,:,:) ! absorbed PAR, sunlit leaves [W/m2 ground]
     real(r8), allocatable :: parsha_z(:,:,:) ! absorbed PAR, shaded leaves [W/m2 ground]
     real(r8), allocatable :: laisun_z(:,:,:) ! sunlit LAI [m2/m2]
     real(r8), allocatable :: laisha_z(:,:,:) ! shaded LAI [m2/m2]

   contains

     procedure, public :: Init
     procedure, public :: RecordDay
     procedure, public :: RecordLightProfile
     procedure, public :: Write

  end type history_type

contains

  ! ==========================================================================

  subroutine Init(this, n_time, n_leaf, n_light, n_year)
    !
    ! DESCRIPTION:
    ! Allocate the time-series arrays and pre-fill the light-profile arrays with
    ! the fill value - only entries at leaf layers 1:nv, for each (year, light
    ! level) actually recorded, ever get overwritten by RecordLightProfile.

    ! ARGUMENTS:
    class(history_type), intent(inout) :: this    ! history object
    integer,              intent(in)    :: n_time  ! days per light level's trajectory
    integer,              intent(in)    :: n_leaf   ! compile-time max leaf+stem layers (nlevleaf)
    integer,              intent(in)    :: n_light  ! number of light levels
    integer,              intent(in)    :: n_year   ! number of simulated years

    this%n_time  = n_time
    this%n_leaf  = n_leaf
    this%n_light = n_light
    this%n_year  = n_year

    allocate(this%dbh(n_time, n_light))
    allocate(this%height(n_time, n_light))
    allocate(this%treelai(n_time, n_light))
    allocate(this%crown_area(n_time, n_light))
    allocate(this%leaf_c(n_time, n_light))
    allocate(this%fnrt_c(n_time, n_light))
    allocate(this%sapw_c(n_time, n_light))
    allocate(this%struct_c(n_time, n_light))
    allocate(this%storage_c(n_time, n_light))
    allocate(this%daily_net_c(n_time, n_light))
    allocate(this%daily_gpp(n_time, n_light))
    allocate(this%npp_acc(n_time, n_light))
    allocate(this%frac_store(n_time, n_light))
    allocate(this%cmort(n_time, n_light))
    allocate(this%n(n_time, n_light))

    allocate(this%parsun_z(n_leaf, n_year, n_light))
    allocate(this%parsha_z(n_leaf, n_year, n_light))
    allocate(this%laisun_z(n_leaf, n_year, n_light))
    allocate(this%laisha_z(n_leaf, n_year, n_light))
    this%parsun_z(:,:,:) = fates_unset_r8
    this%parsha_z(:,:,:) = fates_unset_r8
    this%laisun_z(:,:,:) = fates_unset_r8
    this%laisha_z(:,:,:) = fates_unset_r8

  end subroutine Init

  ! ==========================================================================

  subroutine RecordDay(this, iday_all, ilight, cohort, daily_net_c, daily_gpp, &
    npp_acc, frac_store, cmort)
    !
    ! DESCRIPTION:
    ! Capture one day's daily time series row.

    ! ARGUMENTS:
    class(history_type),     intent(inout) :: this        ! history object
    integer,                  intent(in)    :: iday_all    ! day index within this light level's trajectory (1..n_time)
    integer,                  intent(in)    :: ilight      ! light-level index
    type(fates_cohort_type), intent(in)    :: cohort      ! cohort to record
    real(r8),                 intent(in)    :: daily_net_c ! daily net C (GPP - leaf Rd - nonleaf MR) [kgC/indiv/day]
    real(r8),                 intent(in)    :: daily_gpp   ! daily GPP [kgC/indiv/day]
    real(r8),                 intent(in)    :: npp_acc     ! carbon balance handed to PARTEH (net of growth respiration) [kgC/indiv/day]
    real(r8),                 intent(in)    :: frac_store  ! storage carbon as a fraction of target leaf carbon [-]
    real(r8),                 intent(in)    :: cmort       ! carbon starvation mortality rate [indiv/year]

    this%dbh(iday_all, ilight) = cohort%dbh
    this%height(iday_all, ilight) = cohort%height
    this%treelai(iday_all, ilight) = cohort%treelai
    this%crown_area(iday_all, ilight) = cohort%c_area
    this%leaf_c(iday_all, ilight) = cohort%prt%GetState(leaf_organ, carbon12_element)
    this%fnrt_c(iday_all, ilight) = cohort%prt%GetState(fnrt_organ, carbon12_element)
    this%sapw_c(iday_all, ilight) = cohort%prt%GetState(sapw_organ, carbon12_element)
    this%struct_c(iday_all, ilight) = cohort%prt%GetState(struct_organ, carbon12_element)
    this%storage_c(iday_all, ilight) = cohort%prt%GetState(store_organ, carbon12_element)
    this%daily_net_c(iday_all, ilight) = daily_net_c
    this%daily_gpp(iday_all, ilight) = daily_gpp
    this%npp_acc(iday_all, ilight) = npp_acc
    this%frac_store(iday_all, ilight) = frac_store
    this%cmort(iday_all, ilight) = cmort
    this%n(iday_all, ilight) = cohort%n

  end subroutine RecordDay

  ! ==========================================================================

  subroutine RecordLightProfile(this, iyear, ilight, cohort, light_env)
    !
    ! DESCRIPTION:
    ! Capture the per-leaf-layer light profile snapshot for the given year and
    ! light level. cohort%nv <= n_leaf always (n_leaf is a compile-time maximum),
    ! so layers nv+1:n_leaf are simply never touched here and keep the fill value
    ! set by Init.

    ! ARGUMENTS:
    class(history_type),     intent(inout) :: this      ! history object
    integer,                  intent(in)    :: iyear     ! simulated year index (1..n_year)
    integer,                  intent(in)    :: ilight    ! light-level index
    type(fates_cohort_type), intent(in)    :: cohort    ! cohort this light profile belongs to
    type(light_env_type),    intent(in)    :: light_env ! light environment holding this substep's profile

    this%parsun_z(1:cohort%nv, iyear, ilight) = light_env%parsun_z(:)
    this%parsha_z(1:cohort%nv, iyear, ilight) = light_env%parsha_z(:)
    this%laisun_z(1:cohort%nv, iyear, ilight) = light_env%laisun_z(:)
    this%laisha_z(1:cohort%nv, iyear, ilight) = light_env%laisha_z(:)

  end subroutine RecordLightProfile

  ! ==========================================================================

  subroutine Write(this, out_file, light_frac)
    !
    ! DESCRIPTION:
    ! Writes out the daily whole-cohort time series (dbh, height, treelai, crown
    ! area, PARTEH carbon pools, daily net C/GPP, the carbon balance handed to
    ! PARTEH, frac_store, carbon starvation mortality rate, and cohort number
    ! density), dimensioned (time, light_level), and the annual per-leaf-layer
    ! light profile snapshot, dimensioned (nlevleaf, year, light_level) with
    ! unoccupied layers filled with fates_unset_r8, both across the light sweep.

    ! ARGUMENTS:
    class(history_type), intent(in) :: this          ! history object
    character(len=*),     intent(in) :: out_file      ! output file name
    real(r8),              intent(in) :: light_frac(:) ! swept incident light fractions [0-1]

    ! LOCALS:
    integer, allocatable :: time_idx(:)   ! day index, 1 = first simulated day
    integer, allocatable :: leaf_layer(:) ! leaf-layer index, 1 = top of crown
    integer, allocatable :: year_idx(:)   ! year index, 1 = first simulated year
    integer               :: i            ! looping index
    integer               :: ncid         ! netcdf file id
    character(len=20)     :: dim_names(4) ! dimension names
    integer               :: dimIDs(4)    ! dimension IDs
    integer               :: timeID, leaflayerID, lightfracID, yearID
    integer               :: dbhID, heightID, treelaiID, crownareaID
    integer               :: leafcID, fnrtcID, sapwcID, structcID, storagecID
    integer               :: dailynetcID, dailygppID, nppaccID, fracstoreID
    integer               :: cmortID, nID
    integer               :: parsunID, parshaID, laisunID, laishaID

    ! create day, leaf layer, and year indices
    allocate(time_idx(this%n_time))
    do i = 1, this%n_time
      time_idx(i) = i
    end do
    allocate(leaf_layer(this%n_leaf))
    do i = 1, this%n_leaf
      leaf_layer(i) = i
    end do
    allocate(year_idx(this%n_year))
    do i = 1, this%n_year
      year_idx(i) = i
    end do

    ! dimension names
    dim_names = [character(len=20) :: 'time', 'nlevleaf', 'light_level', 'year']

    ! open file
    call OpenNCFile(trim(out_file), ncid, 'readwrite')

    ! register dimensions
    call RegisterNCDims(ncid, dim_names,                                       &
      (/this%n_time, this%n_leaf, this%n_light, this%n_year/), 4, dimIDs)

    ! register day index
    call RegisterVar(ncid, dim_names(1), dimIDs(1:1), type_int,                &
      [character(len=20)  :: 'units', 'long_name'],                           &
      [character(len=150) :: '', 'day index within this light level''s trajectory, 1 = first simulated day'], &
      2, timeID)

    ! register leaf layer
    call RegisterVar(ncid, dim_names(2), dimIDs(2:2), type_int,                &
      [character(len=20)  :: 'units', 'long_name'],                           &
      [character(len=150) :: '', 'leaf layer index, 1 = top of crown (compile-time max nlevleaf; unoccupied layers filled)'], &
      2, leaflayerID)

    ! register light level
    call RegisterVar(ncid, dim_names(3), dimIDs(3:3), type_double,             &
      [character(len=20)  :: 'units', 'long_name'],                           &
      [character(len=150) :: 'fraction', 'incident light fraction, relative to full sun'], &
      2, lightfracID)

    ! register year
    call RegisterVar(ncid, dim_names(4), dimIDs(4:4), type_int,                &
      [character(len=20)  :: 'units', 'long_name'],                           &
      [character(len=150) :: '', 'simulated year index, 1 = first year'], 2, yearID)

    ! register daily whole-cohort time series, dimensioned (time, light_level)
    call RegisterVar(ncid, 'dbh', (/dimIDs(1), dimIDs(3)/), type_double,       &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'cm', 'dbh'], 3, dbhID)

    call RegisterVar(ncid, 'height', (/dimIDs(1), dimIDs(3)/), type_double,    &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'm', 'height'], 3, heightID)

    call RegisterVar(ncid, 'treelai', (/dimIDs(1), dimIDs(3)/), type_double,   &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'm2 m-2', 'total leaf area index'], &
      3, treelaiID)

    call RegisterVar(ncid, 'crown_area', (/dimIDs(1), dimIDs(3)/), type_double, &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'm2', 'crown area'], 3, crownareaID)

    call RegisterVar(ncid, 'leaf_c', (/dimIDs(1), dimIDs(3)/), type_double,    &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'kgC indiv-1', 'leaf carbon'], &
      3, leafcID)

    call RegisterVar(ncid, 'fnrt_c', (/dimIDs(1), dimIDs(3)/), type_double,    &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'kgC indiv-1', 'fine root carbon'], &
      3, fnrtcID)

    call RegisterVar(ncid, 'sapw_c', (/dimIDs(1), dimIDs(3)/), type_double,    &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'kgC indiv-1', 'sapwood carbon'], &
      3, sapwcID)

    call RegisterVar(ncid, 'struct_c', (/dimIDs(1), dimIDs(3)/), type_double,  &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'kgC indiv-1', 'structural carbon'], &
      3, structcID)

    call RegisterVar(ncid, 'storage_c', (/dimIDs(1), dimIDs(3)/), type_double, &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'kgC indiv-1', 'storage carbon'], &
      3, storagecID)

    call RegisterVar(ncid, 'daily_net_c', (/dimIDs(1), dimIDs(3)/), type_double, &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'kgC indiv-1 day-1',          &
      'daily net carbon (GPP - leaf dark resp - nonleaf MR)'], 3, dailynetcID)

    call RegisterVar(ncid, 'daily_gpp', (/dimIDs(1), dimIDs(3)/), type_double, &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'kgC indiv-1 day-1', 'daily GPP'], &
      3, dailygppID)

    call RegisterVar(ncid, 'npp_acc', (/dimIDs(1), dimIDs(3)/), type_double,   &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'kgC indiv-1 day-1',          &
      'carbon balance handed to PARTEH (net of growth respiration)'], 3, nppaccID)

    call RegisterVar(ncid, 'frac_store', (/dimIDs(1), dimIDs(3)/), type_double, &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', '-',                        &
      'storage carbon as a fraction of target leaf carbon'], 3, fracstoreID)

    call RegisterVar(ncid, 'cmort', (/dimIDs(1), dimIDs(3)/), type_double,     &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'indiv yr-1',                &
      'carbon starvation mortality rate'], 3, cmortID)

    call RegisterVar(ncid, 'n', (/dimIDs(1), dimIDs(3)/), type_double,         &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'time light_level', 'indiv',                    &
      'cohort number density (surviving fraction of the original recruitment cohort)'], &
      3, nID)

    ! register the annual per-leaf-layer light profile, dimensioned
    ! (nlevleaf, year, light_level), with unoccupied layers above that year's nv
    ! filled with fates_unset_r8 (registered below as each variable's _FillValue)
    call RegisterVar(ncid, 'parsun_z', (/dimIDs(2), dimIDs(4), dimIDs(3)/), type_double, &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'nlevleaf year light_level', 'W m-2',            &
      'absorbed PAR, sunlit leaves, per unit ground area (first day of year, solar noon)'], &
      3, parsunID)
    call RegisterFillValue(ncid, parsunID, fates_unset_r8)

    call RegisterVar(ncid, 'parsha_z', (/dimIDs(2), dimIDs(4), dimIDs(3)/), type_double, &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'nlevleaf year light_level', 'W m-2',            &
      'absorbed PAR, shaded leaves, per unit ground area (first day of year, solar noon)'], &
      3, parshaID)
    call RegisterFillValue(ncid, parshaID, fates_unset_r8)

    call RegisterVar(ncid, 'laisun_z', (/dimIDs(2), dimIDs(4), dimIDs(3)/), type_double, &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'nlevleaf year light_level', 'm2 m-2',           &
      'sunlit LAI (first day of year, solar noon)'], 3, laisunID)
    call RegisterFillValue(ncid, laisunID, fates_unset_r8)

    call RegisterVar(ncid, 'laisha_z', (/dimIDs(2), dimIDs(4), dimIDs(3)/), type_double, &
      [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
      [character(len=150) :: 'nlevleaf year light_level', 'm2 m-2',           &
      'shaded LAI (first day of year, solar noon)'], 3, laishaID)
    call RegisterFillValue(ncid, laishaID, fates_unset_r8)

    ! finish defining variables
    call EndNCDef(ncid)

    ! write out data
    call WriteVar(ncid, timeID, time_idx(:))
    call WriteVar(ncid, leaflayerID, leaf_layer(:))
    call WriteVar(ncid, lightfracID, light_frac(:))
    call WriteVar(ncid, yearID, year_idx(:))
    call WriteVar(ncid, dbhID, this%dbh(:,:))
    call WriteVar(ncid, heightID, this%height(:,:))
    call WriteVar(ncid, treelaiID, this%treelai(:,:))
    call WriteVar(ncid, crownareaID, this%crown_area(:,:))
    call WriteVar(ncid, leafcID, this%leaf_c(:,:))
    call WriteVar(ncid, fnrtcID, this%fnrt_c(:,:))
    call WriteVar(ncid, sapwcID, this%sapw_c(:,:))
    call WriteVar(ncid, structcID, this%struct_c(:,:))
    call WriteVar(ncid, storagecID, this%storage_c(:,:))
    call WriteVar(ncid, dailynetcID, this%daily_net_c(:,:))
    call WriteVar(ncid, dailygppID, this%daily_gpp(:,:))
    call WriteVar(ncid, nppaccID, this%npp_acc(:,:))
    call WriteVar(ncid, fracstoreID, this%frac_store(:,:))
    call WriteVar(ncid, cmortID, this%cmort(:,:))
    call WriteVar(ncid, nID, this%n(:,:))
    call WriteVar(ncid, parsunID, this%parsun_z(:,:,:))
    call WriteVar(ncid, parshaID, this%parsha_z(:,:,:))
    call WriteVar(ncid, laisunID, this%laisun_z(:,:,:))
    call WriteVar(ncid, laishaID, this%laisha_z(:,:,:))

    ! close the file
    call CloseNCFile(ncid)

  end subroutine Write

end module FatesTestHistoryMod
