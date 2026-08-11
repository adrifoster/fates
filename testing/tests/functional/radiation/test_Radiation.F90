program FatesRadiation

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesInterfaceTypesMod,      only : numpft
  use PRTParametersMod,            only : prt_params
  use FatesRadiationMemMod,        only : num_swb, ivis, inir
  use TwoStreamMLPEMod,            only : twostream_type, air_ft, normalized_upper_boundary
  use FatesTwoStreamUtilsMod,      only : TransferRadParams
  use FatesUnitTestParamReaderMod, only : ReadParameters
  use FatesArgumentUtils,          only : command_line_arg

  implicit none

  ! LOCALS:
  character(len=:), allocatable :: param_file   ! input parameter file
  type(twostream_type)          :: twostream    ! two-stream radiation solver instance
  integer                       :: n_scr        ! two-stream solve scratch space size
  integer                       :: iv           ! vai looping index
  integer                       :: ib           ! band looping index (matches ivis, inir)
  real(r8), allocatable         :: taulamb(:)   ! two-stream solve scratch space
  real(r8), allocatable         :: omega(:,:)   ! two-stream solve scratch space
  integer,  allocatable         :: ipiv(:)      ! two-stream solve scratch space (LAPACK pivots)
  real(r8)                      :: vai_top      ! vai bin upper (top) bound [m2/m2]
  real(r8)                      :: vai_bot      ! vai bin lower (bottom) bound [m2/m2]
  real(r8)                      :: rb_abs       ! total absorbed beam radiation, this bin (dummy) [W/m2 ground]
  real(r8)                      :: rd_abs       ! total absorbed diffuse radiation, this bin (dummy) [W/m2 ground]
  real(r8)                      :: r_abs_snow   ! absorbed radiation by snow, this bin (dummy, no snow) [W/m2 ground]
  logical                       :: call_fail    ! GetAbsRad failure flag

  ! two-stream solve outputs not needed for the VAI profile (canopy-integrated diagnostics)
  real(r8) :: albedo_beam, albedo_diff, consv_err
  real(r8) :: frac_abs_beam, frac_abs_diff
  real(r8) :: frac_beam_grnd, frac_diff_grnd_beam, frac_diff_grnd_diff

  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'radiation_out.nc' ! output file
  integer,  parameter :: pft           = 1     ! plant functional type of the cohort
  real(r8), parameter :: lai           = 2.0_r8 ! cohort leaf area index [m2/m2]
  real(r8), parameter :: sai           = 0.5_r8 ! cohort stem area index [m2/m2]
  real(r8), parameter :: area_frac     = 0.9_r8 ! fraction of ground covered by the cohort [-]
  real(r8), parameter :: ground_albedo = 0.1_r8 ! ground albedo, diffuse and beam [-]
  real(r8), parameter :: frac_snow     = 0.0_r8 ! canopy snow-covered fraction [-]
  real(r8), parameter :: cosz          = 1.0_r8 ! cosine of solar zenith angle (overhead sun) [-]
  real(r8), parameter :: r_beam_atm    = 1.0_r8 ! normalized incident beam radiation at canopy top [-]
  real(r8), parameter :: r_diff_atm    = 1.0_r8 ! normalized incident diffuse radiation at canopy top [-]
  integer,  parameter :: n_vai         = 100    ! number of vai increments to resolve

  ! band-resolved profile output, dimensioned (vai, band)
  real(r8) :: vai(n_vai)                     ! cumulative vai at the bottom of each bin [m2/m2]
  real(r8) :: r_beam(n_vai, num_swb)         ! normalized downwelling beam intensity [-]
  real(r8) :: r_diff_dn(n_vai, num_swb)      ! normalized downwelling diffuse intensity [-]
  real(r8) :: r_diff_up(n_vai, num_swb)      ! normalized upwelling diffuse intensity [-]
  real(r8) :: rb_abs_leaf(n_vai, num_swb)    ! beam radiation absorbed by leaves, this bin [W/m2 ground]
  real(r8) :: rd_abs_leaf(n_vai, num_swb)    ! diffuse radiation absorbed by leaves, this bin [W/m2 ground]
  real(r8) :: r_abs_stem(n_vai, num_swb)     ! beam+diffuse radiation absorbed by stems, this bin [W/m2 ground]
  real(r8) :: leaf_sun_frac(n_vai, num_swb)  ! sunlit fraction of leaf area, this bin [-]

  interface

    subroutine WriteRadiationData(out_file, n_vai, num_swb, vai, r_beam, r_diff_dn, &
      r_diff_up, rb_abs_leaf, rd_abs_leaf, r_abs_stem, leaf_sun_frac)

      use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
      use FatesUnitTestIOMod, only : WriteVar, RegisterVarAtts
      use FatesUnitTestIOMod, only : type_double, type_int
      use FatesConstantsMod,  only : r8 => fates_r8
      implicit none

      character(len=*), intent(in) :: out_file
      integer,          intent(in) :: n_vai, num_swb
      real(r8),         intent(in) :: vai(:)
      real(r8),         intent(in) :: r_beam(:,:)
      real(r8),         intent(in) :: r_diff_dn(:,:)
      real(r8),         intent(in) :: r_diff_up(:,:)
      real(r8),         intent(in) :: rb_abs_leaf(:,:)
      real(r8),         intent(in) :: rd_abs_leaf(:,:)
      real(r8),         intent(in) :: r_abs_stem(:,:)
      real(r8),         intent(in) :: leaf_sun_frac(:,:)
    end subroutine WriteRadiationData

  end interface

  ! read in parameter file name from command line
  param_file = command_line_arg(1)

  ! read in parameter file, including per-pft optical properties
  ! (rhol/rhos/taul/taus/xl/clumping_index) via EDPftvarcon
  ! also sets numpft
  call ReadParameters(param_file)

  ! TransferRadParams (radiation/FatesTwoStreamUtilsMod.F90) sizes rad_params off
  ! of the global FatesInterfaceTypesMod::numpft, set in ReadParameters
  call TransferRadParams()

  ! Build a single-cohort canopy: one vegetated column plus an "air" column
  ! standing in for the fraction of ground the cohort does not cover.
  call twostream%AllocInitTwoStream((/ivis,inir/), 1, 2)
  twostream%scelg(1,1)%pft  = pft
  twostream%scelg(1,1)%area = area_frac
  twostream%scelg(1,1)%lai  = lai
  twostream%scelg(1,1)%sai  = sai
  twostream%scelg(1,2)%pft  = air_ft
  twostream%scelg(1,2)%area = 1._r8 - area_frac
  twostream%scelg(1,2)%lai  = 0._r8
  twostream%scelg(1,2)%sai  = 0._r8
  twostream%n_col(1) = 2
  call twostream%GetNSCel()

  twostream%band(ivis)%albedo_grnd_diff = ground_albedo
  twostream%band(ivis)%albedo_grnd_beam = ground_albedo
  twostream%band(inir)%albedo_grnd_diff = ground_albedo
  twostream%band(inir)%albedo_grnd_beam = ground_albedo

  twostream%force_prep = .true.
  call twostream%CanopyPrep(frac_snow)
  call twostream%ZenithPrep(cosz)

  ! Scratch space for the linear solve, sized as in production usage
  ! (radiation/FatesTwoStreamUtilsMod.F90:FatesConstructRadElements)
  n_scr = 2*twostream%n_scel + 8
  allocate(taulamb(n_scr), omega(n_scr,n_scr), ipiv(n_scr))

  ! vai bin edges from the top of the canopy (vai=0) to the bottom (vai=lai+sai)
  do iv = 1, n_vai
    vai(iv) = (lai + sai) * real(iv, r8) / real(n_vai, r8)
  end do

  do ib = 1, num_swb

    call twostream%Solve(ib, normalized_upper_boundary, r_beam_atm, r_diff_atm, &
      taulamb, omega, ipiv, albedo_beam, albedo_diff, consv_err,                &
      frac_abs_beam, frac_abs_diff, frac_beam_grnd, frac_diff_grnd_beam,        &
      frac_diff_grnd_diff)

    ! Solve() unsets Rbeam_atm/Rdiff_atm for a normalized solution (see
    ! TwoStreamMLPEMod.F90:Solve)
    ! reassign explicitly so GetRdDn/GetRdUp/GetRb/GetAbsRad, which read them directly, 
    ! return the normalized profile.
    twostream%band(ib)%Rbeam_atm = r_beam_atm
    twostream%band(ib)%Rdiff_atm = r_diff_atm

    do iv = 1, n_vai

      r_diff_dn(iv,ib) = twostream%GetRdDn(1, 1, ib, vai(iv))
      r_diff_up(iv,ib) = twostream%GetRdUp(1, 1, ib, vai(iv))
      r_beam(iv,ib)    = twostream%GetRb(1, 1, ib, vai(iv))

      vai_bot = vai(iv)
      if (iv == 1) then
        vai_top = 0._r8
      else
        vai_top = vai(iv-1)
      end if

      call twostream%GetAbsRad(1, 1, ib, vai_top, vai_bot, rb_abs, rd_abs,        &
        rd_abs_leaf(iv,ib), rb_abs_leaf(iv,ib), r_abs_stem(iv,ib), r_abs_snow,    &
        leaf_sun_frac(iv,ib), call_fail)

      if (call_fail) then
        write(*,*) 'GetAbsRad failed: band=', ib, ' vai_top=', vai_top, ' vai_bot=', vai_bot
        call abort()
      end if

    end do

  end do

  ! write out data to netcdf file
  call WriteRadiationData(out_file, n_vai, num_swb, vai, r_beam, r_diff_dn,      &
    r_diff_up, rb_abs_leaf, rd_abs_leaf, r_abs_stem, leaf_sun_frac)

end program FatesRadiation

! ----------------------------------------------------------------------------------------

subroutine WriteRadiationData(out_file, n_vai, num_swb, vai, r_beam, r_diff_dn, &
  r_diff_up, rb_abs_leaf, rd_abs_leaf, r_abs_stem, leaf_sun_frac)
  !
  ! DESCRIPTION:
  ! Writes out data from the radiation test
  !
  use FatesConstantsMod,  only : r8 => fates_r8
  use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod, only : WriteVar
  use FatesUnitTestIOMod, only : RegisterVarAtts
  use FatesUnitTestIOMod, only : EndNCDef
  use FatesUnitTestIOMod, only : type_double, type_int

  implicit none

  ! ARGUMENTS:
  character(len=*), intent(in) :: out_file             ! output file name
  integer,          intent(in) :: n_vai, num_swb        ! size of vai and band dimensions
  real(r8),         intent(in) :: vai(:)                ! cumulative vai at bottom of each bin [m2/m2]
  real(r8),         intent(in) :: r_beam(:,:)            ! normalized downwelling beam intensity [-]
  real(r8),         intent(in) :: r_diff_dn(:,:)         ! normalized downwelling diffuse intensity [-]
  real(r8),         intent(in) :: r_diff_up(:,:)         ! normalized upwelling diffuse intensity [-]
  real(r8),         intent(in) :: rb_abs_leaf(:,:)       ! beam radiation absorbed by leaves [W/m2 ground]
  real(r8),         intent(in) :: rd_abs_leaf(:,:)       ! diffuse radiation absorbed by leaves [W/m2 ground]
  real(r8),         intent(in) :: r_abs_stem(:,:)        ! beam+diffuse radiation absorbed by stems [W/m2 ground]
  real(r8),         intent(in) :: leaf_sun_frac(:,:)     ! sunlit fraction of leaf area [-]

  ! LOCALS:
  integer, allocatable :: band_indices(:) ! array of band indices to write out (1=vis, 2=nir)
  integer               :: i              ! looping index
  integer               :: ncid           ! netcdf file id
  character(len=8)      :: dim_names(2)   ! dimension names
  integer               :: dimIDs(2)      ! dimension IDs
  integer               :: vaiID, bandID  ! variable IDs for dimensions
  integer               :: rbeamID, rdiffdnID, rdiffupID
  integer               :: rbabsleafID, rdabsleafID
  integer               :: rabsstemID, sunfracID

  ! create band indices
  allocate(band_indices(num_swb))
  do i = 1, num_swb
    band_indices(i) = i
  end do

  ! dimension names
  dim_names = [character(len=8) :: 'vai', 'band']

  ! open file
  call OpenNCFile(trim(out_file), ncid, 'readwrite')

  ! register dimensions
  call RegisterNCDims(ncid, dim_names, (/n_vai, num_swb/), 2, dimIDs)

  ! register vai
  call RegisterVarAtts(ncid, dim_names(1), dimIDs(1:1), type_double, 'm2 m-2',          &
    'cumulative vegetation area index from canopy top (bin bottom)', vaiID)

  ! register band
  call RegisterVarAtts(ncid, dim_names(2), dimIDs(2:2), type_int, '',                   &
    'shortwave band (1=visible, 2=near-infrared)', bandID)

  ! register normalized downwelling beam intensity
  call RegisterVarAtts(ncid, 'r_beam', dimIDs(1:2), type_double, '-',                   &
    'normalized downwelling beam radiation intensity', rbeamID, coordinates='vai band')

  ! register normalized downwelling diffuse intensity
  call RegisterVarAtts(ncid, 'r_diff_dn', dimIDs(1:2), type_double, '-',                &
    'normalized downwelling diffuse radiation intensity', rdiffdnID,                    &
    coordinates='vai band')

  ! register normalized upwelling diffuse intensity
  call RegisterVarAtts(ncid, 'r_diff_up', dimIDs(1:2), type_double, '-',                &
    'normalized upwelling diffuse radiation intensity', rdiffupID,                      &
    coordinates='vai band')

  ! register beam radiation absorbed by leaves
  call RegisterVarAtts(ncid, 'rb_abs_leaf', dimIDs(1:2), type_double, 'W m-2',          &
    'beam radiation absorbed by leaves, per vai bin', rbabsleafID,                      &
    coordinates='vai band')

  ! register diffuse radiation absorbed by leaves
  call RegisterVarAtts(ncid, 'rd_abs_leaf', dimIDs(1:2), type_double, 'W m-2',          &
    'diffuse radiation absorbed by leaves, per vai bin', rdabsleafID,                   &
    coordinates='vai band')

  ! register radiation absorbed by stems
  call RegisterVarAtts(ncid, 'r_abs_stem', dimIDs(1:2), type_double, 'W m-2',           &
    'beam+diffuse radiation absorbed by stems, per vai bin', rabsstemID,                &
    coordinates='vai band')

  ! register sunlit fraction of leaves
  call RegisterVarAtts(ncid, 'leaf_sun_frac', dimIDs(1:2), type_double, '-',            &
    'sunlit fraction of leaf area, per vai bin', sunfracID, coordinates='vai band')

  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, vaiID, vai(:))
  call WriteVar(ncid, bandID, band_indices(:))
  call WriteVar(ncid, rbeamID, r_beam(:,:))
  call WriteVar(ncid, rdiffdnID, r_diff_dn(:,:))
  call WriteVar(ncid, rdiffupID, r_diff_up(:,:))
  call WriteVar(ncid, rbabsleafID, rb_abs_leaf(:,:))
  call WriteVar(ncid, rdabsleafID, rd_abs_leaf(:,:))
  call WriteVar(ncid, rabsstemID, r_abs_stem(:,:))
  call WriteVar(ncid, sunfracID, leaf_sun_frac(:,:))

  ! close the file
  call CloseNCFile(ncid)

end subroutine WriteRadiationData
