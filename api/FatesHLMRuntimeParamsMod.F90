module FatesHLMRuntimeParamsMod
  ! --------------------------------------------------------------------------------------
  ! This module stores a class with parameters that are dictated by the Host Land Model
  ! These are private to the HLMRuntimeParams class, and can only be accessed via 
  ! setters and getters.
  ! --------------------------------------------------------------------------------------

  use FatesConstantsMod,   only : r8 => fates_r8
  use FatesConstantsMod,   only : nearzero
  use FatesConstantsMod,   only : maxSWb
  use FatesConstantsMod,   only : fates_check_param_set
  use FatesConstantsMod,   only : ivis, inir
  use FatesGlobals,        only : fates_log
  use FatesGlobals,        only : endrun => fates_endrun
  use FatesFinalParamType, only : real_param, character_param
  use FatesFinalParamType, only : logical_param, integer_param
  use FatesFinalParamType, only : char_len
  use PRTGenericMod,       only : prt_cnp_flex_allom_hyp

  use shr_log_mod,         only : errMsg => shr_log_errMsg

  implicit none 
  private

  character(len=*), parameter :: sourcefile = __FILE__

  type :: hlm_runtime_params
    
    type(integer_param),   private :: num_swb                     ! number of broadbands in the shortwave radiation spectrum to track
                                                                     !   typically 2 as a default, VIS & NIR
    type(integer_param),   private :: ivis                        ! the HLM's array index for the VIS portion of the spectrum in SW radiation arrays
    type(integer_param),   private :: inir                        ! the HLM's array index for the NIR portion of the spectrum in SW radiation arrays
    type(integer_param),   private :: maxlevsoil                  ! maximum number of soil layers
    type(logical_param),   private :: is_restart                  ! is the HLM signaling that this is a restart? 
    type(character_param), private :: hlm_name                    ! character string passed by the HLM
                                                                     !   is used during the processing of IO data, 
                                                                     !   so that FATES knows which IO variables it 
                                                                     !   should prepare.  For instance
                                                                     !   ATS, ALM and CLM will only want variables 
                                                                     !   specficially packaged for them.
                                                                     !   This string sets which filter is enacted.
    type(character_param), private :: decomp_scheme               ! which soil decomposition scheme is active
                                                                     !    expected values are one of CENTURY,MIMICS,CTC
    type(character_param), private :: nutrient_scheme             ! which soil nutrient scheme is active
                                                                     !    current options with:
                                                                     !    E3SM: RD, ECA
                                                                     !    CESM: NONE
                                                                     !    ATS: ?
                                                                     !    NORESM: ?
    type(integer_param),   private :: nitrogen_spec               ! which nitrogen species active (if any):
                                                                     !    0: none
                                                                     !    1: NH4 only
                                                                     !    2: NH4 and NO3
    type(logical_param),   private :: phosphorus_spec             ! is phosphorous is turned on in the HLM 
    type(real_param),      private :: stepsize                    ! step-size of the host land model
                                                                     !    moreover, this is the shortest main-model timestep
                                                                     !    at which FATES will be called on the main model integration loop
    type(real_param),      private :: hio_ignore_val              ! value flushed to history diagnostics which the HLM will ignore in aggregation functions
    type(logical_param),   private :: masterproc                  ! is this the master processor? 
                                                                     !   typically useful for knowing if current processor should print log messages
    type(integer_param),   private :: ipedof                      ! HLM pedotransfer index
                                                                     !   only used by the plant hydraulics submodule to check and/or enable
                                                                     !   consistency between the pedotransfer functions of the HLM
                                                                     !   and how it moves and stores water in its rhizosphere shells
    type(integer_param),   private :: parteh_mode                 ! which Plant Allocation and Reactive Transport (exensible) Hypothesis (PARTEH) to use
    type(logical_param),   private :: use_ch4                     ! whether the methane model in ELM/CLM is active 
                                                                     !   if active, boundary conditions need to be prepped
    type(logical_param),   private :: use_vertsoilc               ! whether or not the HLM is using vertically discretized soil carbon 
    type(integer_param),   private :: spitfire_mode               ! SPITFIRE (fire model) mode
                                                                     !    0: off
                                                                     !    1: constant ignitions
                                                                     !   >1: ignitions from external data sources (lightning and/or anthropogenic)
    type(logical_param),   private :: use_lu_harvest              ! whether or not to use harvest data from the HLM 
                                                                     !    if 1, automatically sets use_logging to 1
    type(integer_param),   private :: num_lu_harvest_cats         ! number of HLM harvest categories (e.g. primary forest harvest, secondary young forest harvest, etc.)
                                                                     !    this is the first dimension of harvest_rates in dynHarvestMod and 
                                                                     !    bc_in%hlm_harvest_rates and bc_in%hlm_harvest_catnames
    type(integer_param),   private :: sf_nofire_def               ! definition of a no-fire case for spitfire_mode
    type(integer_param),   private :: sf_scalar_lightning_def     ! definition of a scalar-lightning case for spitfire_mode
    type(integer_param),   private :: sf_successful_ignitions_def ! definition of a successful-ignition dataset case for spitfire_mode
    type(integer_param),   private :: sf_anthro_ignitions_def     ! definition of an anthropogenic-ignition dataset case for spitfire_mode
    type(logical_param),   private :: use_logging                 ! whether or not to use the logging module
                                                                     !    if use_lu_harvest=0, logging is determined by 
                                                                     !    the FATES parameter file, otherwise this flag automatically
                                                                     !    set to 1 and logging determined by land use harvest input from HLM
    type(logical_param),   private :: use_planthydro              ! whether or not to use plant hydraulics (bchristo/xu methods)
    type(logical_param),   private :: use_cohort_age_tracking     ! whether or not to use cohort age tracking 
    type(logical_param),   private :: use_tree_damage             ! whether or not to turn on the tree damage model 
    type(logical_param),   private :: use_ed_st3                  ! whether or not to use ST(atic) ST(and) ST(ructure) mode (ST3) 
                                                                     !    essentially, this gives us the ability to turn off "dynamics," i.e.
                                                                     !    growth, disturbance, and mortality
                                                                     !    EXPERIMENTAL! default should be FALSE, cannot be true with perscribed_phys
    type(logical_param),   private :: use_ed_prescribed_phys      ! whether or not to use prescribed physiology (i.e. opposite of ST3) 
                                                                     !    turns off fast processes like photosynthesis and respiration, prescribe NPP
                                                                     !    deafult should be FALSE, cannot be true with ST3 on
    type(logical_param),   private :: use_inventory_init          ! whether or not to initialize simulation from inventory file 
                                                                     !    if on, an inventory control file must be specified
    type(character_param), private :: inventory_ctrl_file         ! full path to the inventory controle file that specifies available inventory
                                                                     !    datasets, their locations, and their formats
                                                                     !    only needs to be defined when use_inventory_init = 1
    type(logical_param),   private :: use_fixed_biogeog           !  flag to use FATES fixed biogeography mode 
    type(logical_param),   private :: use_nocomp                  !  flag to use FATES no competition mode 
    type(logical_param),   private :: use_sp                      !  flag to use FATES satellite phenology (LAI) mode 

    contains 

      procedure :: check_all_set
      procedure :: check_compatibility
      procedure :: set_num_swb, get_num_swb
      procedure :: set_ivis, get_ivis
      procedure :: set_inir, get_inir
      procedure :: set_maxlevsoil, get_maxlevsoil
      procedure :: set_is_restart, get_is_restart
      procedure :: set_hlm_name, get_hlm_name
      procedure :: set_decomp_scheme, get_decomp_scheme
      procedure :: set_nutrient_scheme, get_nutrient_scheme
      procedure :: set_nitrogen_spec, get_nitrogen_spec
      procedure :: set_phosphorus_spec, get_phosphorus_spec
      procedure :: set_stepsize, get_stepsize
      procedure :: set_hio_ignore_val, get_hio_ignore_val
      procedure :: set_masterproc, get_masterproc
      procedure :: set_ipedof, get_ipedof
      procedure :: set_parteh_mode, get_parteh_mode
      procedure :: set_use_ch4, get_use_ch4
      procedure :: set_use_vertsoilc, get_use_vertsoilc
      procedure :: set_spitfire_mode, get_spitfire_mode
      procedure :: set_use_lu_harvest, get_use_lu_harvest
      procedure :: set_num_lu_harvest_cats, get_num_lu_harvest_cats
      procedure :: set_sf_nofire_def, get_sf_nofire_def
      procedure :: set_sf_scalar_lightning_def, get_sf_scalar_lightning_def
      procedure :: set_sf_successful_ignitions_def, get_sf_successful_ignitions_def
      procedure :: set_sf_anthro_ignitions_def, get_sf_anthro_ignitions_def
      procedure :: set_use_logging, get_use_logging
      procedure :: set_use_planthydro, get_use_planthydro
      procedure :: set_use_cohort_age_tracking, get_use_cohort_age_tracking
      procedure :: set_use_tree_damage, get_use_tree_damage
      procedure :: set_use_ed_st3, get_use_ed_st3
      procedure :: set_use_ed_prescribed_phys, get_use_ed_prescribed_phys
      procedure :: set_use_inventory_init, get_use_inventory_init
      procedure :: set_inventory_ctrl_file, get_inventory_ctrl_file
      procedure :: set_use_fixed_biogeog, get_use_fixed_biogeog
      procedure :: set_use_nocomp, get_use_nocomp
      procedure :: set_use_sp, get_use_sp
    
  end type hlm_runtime_params

  ! allow this instance to be passed around...
  ! TODO - change so that only the top-level drivers can access this
  type(hlm_runtime_params), public :: hlm_runtime_params_inst

  ! ======================================================================================

  private :: throw_if_unset
  private :: throw_if_incompatible

  contains 

    ! ------------------------------------------------------------------------------------

    subroutine check_all_set(this)
      !
      !  DESCRIPTION:
      !  check that all runtime parameters have been set
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this ! runtime params object

      if (.not. this%num_swb%already_set()) then
        call throw_if_unset('num_swb')
      else if (.not. this%ivis%already_set()) then
        call throw_if_unset('ivis')
      else if (.not. this%inir%already_set()) then
        call throw_if_unset('inir')
      else if (.not. this%maxlevsoil%already_set()) then
        call throw_if_unset('maxlevsoil')
      else if (.not. this%is_restart%already_set()) then
        call throw_if_unset('is_restart')
      else if (.not. this%hlm_name%already_set()) then
        call throw_if_unset('hlm_name')
      else if (.not. this%decomp_scheme%already_set()) then
        call throw_if_unset('decomp_scheme')
      else if (.not. this%nutrient_scheme%already_set()) then
        call throw_if_unset('nutrient_scheme')
      else if (.not. this%nitrogen_spec%already_set()) then
        call throw_if_unset('nitrogen_spec')
      else if (.not. this%phosphorus_spec%already_set()) then
        call throw_if_unset('phosphorus_spec')
      else if (.not. this%stepsize%already_set()) then
        call throw_if_unset('stepsize')
      else if (.not. this%hio_ignore_val%already_set()) then
        call throw_if_unset('hio_ignore_val')
      else if (.not. this%masterproc%already_set()) then
        call throw_if_unset('masterproc')
      else if (.not. this%ipedof%already_set()) then
        call throw_if_unset('ipedof')
      else if (.not. this%parteh_mode%already_set()) then
        call throw_if_unset('parteh_mode')
      else if (.not. this%use_ch4%already_set()) then
        call throw_if_unset('use_ch4')
      else if (.not. this%use_vertsoilc%already_set()) then
        call throw_if_unset('use_vertsoilc')
      else if (.not. this%spitfire_mode%already_set()) then
        call throw_if_unset('spitfire_mode')
      else if (.not. this%use_lu_harvest%already_set()) then
        call throw_if_unset('use_lu_harvest')
      else if (.not. this%num_lu_harvest_cats%already_set()) then
        call throw_if_unset('num_lu_harvest_cats')
      else if (.not. this%sf_nofire_def%already_set()) then
        call throw_if_unset('sf_nofire_def')
      else if (.not. this%sf_scalar_lightning_def%already_set()) then
        call throw_if_unset('sf_scalar_lightning_def')
      else if (.not. this%sf_successful_ignitions_def%already_set()) then
        call throw_if_unset('sf_successful_ignitions_def')
      else if (.not. this%sf_anthro_ignitions_def%already_set()) then
        call throw_if_unset('sf_anthro_ignitions_def')
      else if (.not. this%use_logging%already_set()) then
        call throw_if_unset('use_logging')
      else if (.not. this%use_planthydro%already_set()) then
        call throw_if_unset('use_planthydro')
      else if (.not. this%use_cohort_age_tracking%already_set()) then
        call throw_if_unset('use_cohort_age_tracking')
      else if (.not. this%use_tree_damage%already_set()) then
        call throw_if_unset('use_tree_damage')
      else if (.not. this%use_ed_st3%already_set()) then
        call throw_if_unset('use_ed_st3')
      else if (.not. this%use_ed_prescribed_phys%already_set()) then
        call throw_if_unset('use_ed_prescribed_phys')
      else if (.not. this%use_inventory_init%already_set()) then
        call throw_if_unset('use_inventory_init')
      else if (.not. this%inventory_ctrl_file%already_set()) then
        call throw_if_unset('inventory_ctrl_file')
      else if (.not. this%use_fixed_biogeog%already_set()) then
        call throw_if_unset('use_fixed_biogeog')
      else if (.not. this%use_nocomp%already_set()) then
        call throw_if_unset('use_nocomp')
      else if (.not. this%use_sp%already_set()) then
        call throw_if_unset('use_sp')
      end if
  
    end subroutine check_all_set

    ! ------------------------------------------------------------------------------------

    subroutine check_compatibility(this, prescribed_puptake, prescribed_nuptake,         &
        mort_ip_age_senescence)
      !
      !  DESCRIPTION:
      !  check that parameter values are compatible with one another
      !

      ! ARGUMENTS:
      class(hlm_runtime_params),              intent(inout) :: this                      ! runtime params object
      real(r8),                  allocatable, intent(in)    :: prescribed_puptake(:)     ! prescribed uptake rate for phosphorous, this is the fraction of plant demand
      real(r8),                  allocatable, intent(in)    :: prescribed_nuptake(:)     ! prescribed uptake rate for nitrogen, this is the fraction of plant demand
      real(r8),                  allocatable, intent(in)    :: mort_ip_age_senescence(:) ! inflection point of age-dependent senescence

      ! first make sure all parameters are set
      call this%check_all_set()

      if (this%hlm_name%get_value() == 'CLM' .and. this%parteh_mode%get_value() == 2) then 
        if (sum(abs(prescribed_puptake(:))) < nearzero .and.                             &
          sum(abs(prescribed_nuptake(:))) < nearzero) then
          write(fates_log(),*) 'PARTEH hypothesis 2 is only viable'                                                     
          write(fates_log(),*) 'with forced boundary conditions for CLM (currently)'      
          write(fates_log(),*) 'prescribed_puptake or prescribed_nuptake must > 0'
          call throw_if_incompatible()
        end if
      end if 

      if (this%use_tree_damage%get_value() .and.                                         &
        this%parteh_mode%get_value() == prt_cnp_flex_allom_hyp) then 
        write(fates_log(),*) 'FATES tree damage (use_fates_tree_damage = .true.) is not' 
        write(fates_log(),*) '(yet) compatible with CNP allocation (fates_parteh_mode = 2)'                                      
        call throw_if_incompatible()
      end if 

      if (any(mort_ip_age_senescence < fates_check_param_set) .and.                      &
        .not. this%use_cohort_age_tracking%get_value()) then
        write(fates_log(),*) 'Age-dependent mortality cannot be on if'                 
        write(fates_log(),*) 'cohort age tracking is off.'                                   
        write(fates_log(),*) 'Set use_fates_cohort_age_tracking = .true.'                   
        write(fates_log(),*) 'in FATES namelist options.'                                   
        call throw_if_incompatible()
      end if
      
      if (this%use_inventory_init%get_value() .and.                                      &
        this%use_cohort_age_tracking%get_value()) then 
        write(fates_log(),*) 'Fates inventory init cannot be used with age dependent mortality' 
        write(fates_log(),*) 'Set use_fates_cohort_age_tracking to 0 or turn off inventory init'                                      
        call throw_if_incompatible()
      end if 

      if (this%use_ed_prescribed_phys%get_value() .and. this%use_ed_st3%get_value()) then
        write(fates_log(),*) 'FATES ST3 and prescribed physiology cannot both be turned on.' 
        write(fates_log(),*) 'Review the namelist entries'
        call throw_if_incompatible()
      end if 

      if (.not. this%use_fixed_biogeog%get_value() .and. this%use_sp%get_value()) then 
        write(fates_log(),*) 'SP cannot be on if fixed biogeog mode is off.'
        call throw_if_incompatible()
      end if 

      if (.not. this%use_nocomp%get_value() .and. this%use_sp%get_value()) then 
        write(fates_log(),*) 'SP cannot be on if nocomp mode is off.'
        call throw_if_incompatible()
      end if 

    end subroutine check_compatibility

    ! ------------------------------------------------------------------------------------

    subroutine throw_if_incompatible()
      !
      !  DESCRIPTION:
      !  throws an error (stops model) and if any runtime parameters are incompatible
      !

      write(fates_log(),*) 'Aborting.'
      call endrun(msg=errMsg(sourcefile, __LINE__))

    end subroutine throw_if_incompatible

    ! ------------------------------------------------------------------------------------

    subroutine throw_if_unset(param)
      !
      !  DESCRIPTION:
      !  throws an error and prints a message if a runtime parameter is unset
      !

      ! ARGUMENTS:
      character(len=*), intent(in) :: param ! parameter name

      write(fates_log(), *) 'FATES HLM runtime parameter unset: ' // param
      write(fates_log(),*) 'Aborting.'
      call endrun(msg=errMsg(sourcefile, __LINE__))

    end subroutine throw_if_unset

    ! ------------------------------------------------------------------------------------

    subroutine set_num_swb(this, num_swb)
      !
      !  DESCRIPTION:
      !  sets value for num_swb
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this    ! runtime params object
      integer,                   intent(in)    :: num_swb ! number of broadbands in the shortwave radiation spectrum to track
      
      if (num_swb > maxSWb) then
        write(fates_log(),*) 'FATES sets a maximum number of shortwave bands'
        write(fates_log(),*) 'for some scratch-space, maxSWb'
        write(fates_log(),*) 'it defaults to 2, but can be increased as needed'
        write(fates_log(),*) 'your driver or host model is intending to drive'
        write(fates_log(),*) 'FATES with:', num_swb, ' bands.'
        write(fates_log(),*) 'please increase maxSWb in EDTypes to match'
        write(fates_log(),*) 'or exceed this value'
        call endrun(msg=errMsg(sourcefile, __LINE__))
      else 
        call this%num_swb%set_value(num_swb)
      end if

    end subroutine set_num_swb

    ! ------------------------------------------------------------------------------------

    integer function get_num_swb(this)
      !
      !  DESCRIPTION:
      !  gets num_swb value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_num_swb = this%num_swb%get_value()

    end function get_num_swb

    ! ------------------------------------------------------------------------------------
    
    subroutine set_ivis(this, ivis_in)
      !
      !  DESCRIPTION:
      !  sets value for ivis
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this    ! runtime params object
      integer,                   intent(in)    :: ivis_in ! the HLM's array index for the VIS portion of the spectrum in SW radiation arrays
      
      if (ivis_in .ne. ivis) then 
        write(fates_log(),*) 'FATES assumption about the index of visible shortwave'
        write(fates_log(),*) 'radiation is different from the HLM, exiting'
      else
        call this%ivis%set_value(ivis_in)
      end if 

    end subroutine set_ivis

    ! ------------------------------------------------------------------------------------

    integer function get_ivis(this)
      !
      !  DESCRIPTION:
      !  gets ivis value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_ivis = this%ivis%get_value()

    end function get_ivis

    ! ------------------------------------------------------------------------------------
    
    subroutine set_inir(this, inir_in)
      !
      !  DESCRIPTION:
      !  sets value for inir
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this    ! runtime params object
      integer,                   intent(in)    :: inir_in ! the HLM's array index for the NIR portion of the spectrum in SW radiation arrays
      
      if (inir_in .ne. inir) then 
        write(fates_log(),*) 'FATES assumption about the index of NIR shortwave'
        write(fates_log(),*) 'radiation is different from the HLM, exiting'
      else 
        call this%inir%set_value(inir_in)
      end if 

    end subroutine set_inir

    ! ------------------------------------------------------------------------------------

    integer function get_inir(this)
      !
      !  DESCRIPTION:
      !  gets inir value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_inir = this%inir%get_value()

    end function get_inir

    ! ------------------------------------------------------------------------------------
    
    subroutine set_maxlevsoil(this, maxlevsoil)
      !
      !  DESCRIPTION:
      !  sets value for maxlevsoil
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this       ! runtime params object
      integer,                   intent(in)    :: maxlevsoil ! maximum number of soil layers
      
      call this%maxlevsoil%set_value(maxlevsoil)

    end subroutine set_maxlevsoil

    ! ------------------------------------------------------------------------------------

    integer function get_maxlevsoil(this)
      !
      !  DESCRIPTION:
      !  gets maxlevsoil value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_maxlevsoil = this%maxlevsoil%get_value()

    end function get_maxlevsoil

    ! ------------------------------------------------------------------------------------
    
    subroutine set_is_restart(this, is_restart)
      !
      !  DESCRIPTION:
      !  sets value for is_restart
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this       ! runtime params object
      logical,                   intent(in)    :: is_restart ! is the HLM signaling that this is a restart? [1=TRUE
      
      call this%is_restart%set_value(is_restart)

    end subroutine set_is_restart

    ! ------------------------------------------------------------------------------------

    logical function get_is_restart(this)
      !
      !  DESCRIPTION:
      !  gets is_restart value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_is_restart = this%is_restart%get_value()

    end function get_is_restart

    ! ------------------------------------------------------------------------------------
    
    subroutine set_hlm_name(this, hlm_name)
      !
      !  DESCRIPTION:
      !  sets value for hlm_name
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this     ! runtime params object
      character(char_len),       intent(in)    :: hlm_name ! character string passed by the HLM used during the processing of IO data
      
      call this%hlm_name%set_value(hlm_name)

    end subroutine set_hlm_name

    ! ------------------------------------------------------------------------------------

    character(char_len) function get_hlm_name(this)
      !
      !  DESCRIPTION:
      !  gets hlm_name value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_hlm_name = this%hlm_name%get_value()

    end function get_hlm_name

    ! ------------------------------------------------------------------------------------
    
    subroutine set_decomp_scheme(this, decomp_scheme)
      !
      !  DESCRIPTION:
      !  sets value for decomp_scheme
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this          ! runtime params object
      character(char_len),       intent(in)    :: decomp_scheme ! which soil decomposition scheme is active
      
      if (.not. ((trim(decomp_scheme) == 'MIMICS') .or.                                  &
        (trim(decomp_scheme) == 'CENTURY') .or.                                          &
        (trim(decomp_scheme) == 'CTC') .or.                                              &
        (trim(decomp_scheme) == 'NONE'))) then
        write(fates_log(),*) 'Invalid decomposition scheme.'
        write(fates_log(),*) 'valid: NONE, MIMICS, CENTURY, CTC, yours: ', trim(decomp_scheme)
      else 
        call this%decomp_scheme%set_value(decomp_scheme)
      end if 

    end subroutine set_decomp_scheme

    ! ------------------------------------------------------------------------------------

    character(char_len) function get_decomp_scheme(this)
      !
      !  DESCRIPTION:
      !  gets decomp_scheme value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_decomp_scheme = this%decomp_scheme%get_value()

    end function get_decomp_scheme

    ! ------------------------------------------------------------------------------------
    
    subroutine set_nutrient_scheme(this, nutrient_scheme)
      !
      !  DESCRIPTION:
      !  sets value for nutrient_scheme
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this            ! runtime params object
      character(char_len),       intent(in)    :: nutrient_scheme ! which soil nutrient scheme is active
      
      call this%nutrient_scheme%set_value(nutrient_scheme)

    end subroutine set_nutrient_scheme

    ! ------------------------------------------------------------------------------------

    character(char_len) function get_nutrient_scheme(this)
      !
      !  DESCRIPTION:
      !  gets nutrient_scheme value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_nutrient_scheme = this%nutrient_scheme%get_value()

    end function get_nutrient_scheme

    ! ------------------------------------------------------------------------------------
    
    subroutine set_nitrogen_spec(this, nitrogen_spec)
      !
      !  DESCRIPTION:
      !  sets value for nitrogen_spec
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this          ! runtime params object
      integer,                   intent(in)    :: nitrogen_spec ! which nitrogen species active (if any)
      
      call this%nitrogen_spec%set_value(nitrogen_spec)

    end subroutine set_nitrogen_spec

    ! ------------------------------------------------------------------------------------

    integer function get_nitrogen_spec(this)
      !
      !  DESCRIPTION:
      !  gets nitrogen_spec value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_nitrogen_spec = this%nitrogen_spec%get_value()

    end function get_nitrogen_spec

    ! ------------------------------------------------------------------------------------
    
    subroutine set_phosphorus_spec(this, phosphorus_spec)
      !
      !  DESCRIPTION:
      !  sets value for phosphorus_spec
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this            ! runtime params object
      logical,                   intent(in)    :: phosphorus_spec ! is phosphorous is turned on in the HLM
      
      call this%phosphorus_spec%set_value(phosphorus_spec)

    end subroutine set_phosphorus_spec

    ! ------------------------------------------------------------------------------------

    logical function get_phosphorus_spec(this)
      !
      !  DESCRIPTION:
      !  gets phosphorus_spec value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_phosphorus_spec = this%phosphorus_spec%get_value()

    end function get_phosphorus_spec

    ! ------------------------------------------------------------------------------------
    
    subroutine set_stepsize(this, stepsize)
      !
      !  DESCRIPTION:
      !  sets value for stepsize
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this     ! runtime params object
      real(r8),                  intent(in)    :: stepsize ! step-size of the host land model
      
      call this%stepsize%set_value(stepsize)

    end subroutine set_stepsize

    ! ------------------------------------------------------------------------------------

    real(r8) function get_stepsize(this)
      !
      !  DESCRIPTION:
      !  gets stepsize value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_stepsize = this%stepsize%get_value()

    end function get_stepsize

    ! ------------------------------------------------------------------------------------
    
    subroutine set_hio_ignore_val(this, hio_ignore_val)
      !
      !  DESCRIPTION:
      !  sets value for hio_ignore_val
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this           ! runtime params object
      real(r8),                  intent(in)    :: hio_ignore_val ! value flushed to history diagnostics which the HLM will ignore in aggregation functions
      
      call this%hio_ignore_val%set_value(hio_ignore_val)

    end subroutine set_hio_ignore_val

    ! ------------------------------------------------------------------------------------

    real(r8) function get_hio_ignore_val(this)
      !
      !  DESCRIPTION:
      !  gets hio_ignore_val value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_hio_ignore_val = this%hio_ignore_val%get_value()

    end function get_hio_ignore_val

    ! ------------------------------------------------------------------------------------
    
    subroutine set_masterproc(this, masterproc)
      !
      !  DESCRIPTION:
      !  sets value for masterproc
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this       ! runtime params object
      logical,                   intent(in)    :: masterproc ! is this the master processor? 
      
      call this%masterproc%set_value(masterproc)

    end subroutine set_masterproc

    ! ------------------------------------------------------------------------------------

    logical function get_masterproc(this)
      !
      !  DESCRIPTION:
      !  gets masterproc value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_masterproc = this%masterproc%get_value()

    end function get_masterproc

    ! ------------------------------------------------------------------------------------
    
    subroutine set_ipedof(this, ipedof)
      !
      !  DESCRIPTION:
      !  sets value for ipedof
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this   ! runtime params object
      integer,                   intent(in)    :: ipedof ! HLM pedotransfer index
      
      call this%ipedof%set_value(ipedof)

    end subroutine set_ipedof

    ! ------------------------------------------------------------------------------------

    integer function get_ipedof(this)
      !
      !  DESCRIPTION:
      !  gets ipedof value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_ipedof = this%ipedof%get_value()

    end function get_ipedof

    ! ------------------------------------------------------------------------------------
    
    subroutine set_parteh_mode(this, parteh_mode)
      !
      !  DESCRIPTION:
      !  sets value for parteh_mode
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this        ! runtime params object
      integer,                   intent(in)    :: parteh_mode ! which Plant Allocation and Reactive Transport (exensible) Hypothesis (PARTEH) to use
      
      call this%parteh_mode%set_value(parteh_mode)
    
    end subroutine set_parteh_mode

    ! ------------------------------------------------------------------------------------

    integer function get_parteh_mode(this)
      !
      !  DESCRIPTION:
      !  gets parteh_mode value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_parteh_mode = this%parteh_mode%get_value()

    end function get_parteh_mode

    ! ------------------------------------------------------------------------------------
    
    subroutine set_use_ch4(this, use_ch4)
      !
      !  DESCRIPTION:
      !  sets value for use_ch4
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this    ! runtime params object
      logical,                   intent(in)    :: use_ch4 ! whether the methane model in HLM is active 
      
      call this%use_ch4%set_value(use_ch4)

    end subroutine set_use_ch4

    ! ------------------------------------------------------------------------------------

    logical function get_use_ch4(this)
      !
      !  DESCRIPTION:
      !  gets use_ch4 value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_use_ch4 = this%use_ch4%get_value()

    end function get_use_ch4

    ! ------------------------------------------------------------------------------------
    
    subroutine set_use_vertsoilc(this, use_vertsoilc)
      !
      !  DESCRIPTION:
      !  sets value for use_vertsoilc
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this          ! runtime params object
      logical,                   intent(in)    :: use_vertsoilc ! whether or not the HLM is using vertically discretized soil carbon 
      
      call this%use_vertsoilc%set_value(use_vertsoilc)

    end subroutine set_use_vertsoilc

    ! ------------------------------------------------------------------------------------

    logical function get_use_vertsoilc(this)
      !
      !  DESCRIPTION:
      !  gets use_vertsoilc value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_use_vertsoilc = this%use_vertsoilc%get_value()

    end function get_use_vertsoilc

    ! ------------------------------------------------------------------------------------
    
    subroutine set_spitfire_mode(this, spitfire_mode)
      !
      !  DESCRIPTION:
      !  sets value for spitfire_mode
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this          ! runtime params object
      integer,                   intent(in)    :: spitfire_mode ! SPITFIRE (fire model) mode
      
      call this%spitfire_mode%set_value(spitfire_mode)

    end subroutine set_spitfire_mode

    ! ------------------------------------------------------------------------------------

    integer function get_spitfire_mode(this)
      !
      !  DESCRIPTION:
      !  gets spitfire_mode value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_spitfire_mode = this%spitfire_mode%get_value()

    end function get_spitfire_mode

    ! ------------------------------------------------------------------------------------
    
    subroutine set_use_lu_harvest(this, use_lu_harvest)
      !
      !  DESCRIPTION:
      !  sets value for use_lu_harvest
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this           ! runtime params object
      logical,                   intent(in)    :: use_lu_harvest ! whether or not to use harvest data from the HLM 
      
      call this%use_lu_harvest%set_value(use_lu_harvest)

    end subroutine set_use_lu_harvest

    ! ------------------------------------------------------------------------------------

    logical function get_use_lu_harvest(this)
      !
      !  DESCRIPTION:
      !  gets use_lu_harvest value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_use_lu_harvest = this%use_lu_harvest%get_value()

    end function get_use_lu_harvest

    ! ------------------------------------------------------------------------------------
    
    subroutine set_num_lu_harvest_cats(this, num_lu_harvest_cats)
      !
      !  DESCRIPTION:
      !  sets value for num_lu_harvest_cats
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this                ! runtime params object
      integer,                   intent(in)    :: num_lu_harvest_cats ! number of HLM harvest categories 
      
      if (num_lu_harvest_cats < 0) then 
        write(fates_log(), *) 'The FATES number of hlm harvest cats must be >= 0, exiting'
        call endrun(msg=errMsg(sourcefile, __LINE__))
      else 
        call this%num_lu_harvest_cats%set_value(num_lu_harvest_cats)
      end if 

    end subroutine set_num_lu_harvest_cats

    ! ------------------------------------------------------------------------------------

    integer function get_num_lu_harvest_cats(this)
      !
      !  DESCRIPTION:
      !  gets num_lu_harvest_cats value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_num_lu_harvest_cats = this%num_lu_harvest_cats%get_value()

    end function get_num_lu_harvest_cats

    ! ------------------------------------------------------------------------------------
    
    subroutine set_sf_nofire_def(this, sf_nofire_def)
      !
      !  DESCRIPTION:
      !  sets value for sf_nofire_def
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this          ! runtime params object
      integer,                   intent(in)    :: sf_nofire_def !  definition of a no-fire case for spitfire_mode
      
      call this%sf_nofire_def%set_value(sf_nofire_def)

    end subroutine set_sf_nofire_def

    ! ------------------------------------------------------------------------------------

    integer function get_sf_nofire_def(this)
      !
      !  DESCRIPTION:
      !  gets sf_nofire_def value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_sf_nofire_def = this%sf_nofire_def%get_value()

    end function get_sf_nofire_def

    ! ------------------------------------------------------------------------------------
    
    subroutine set_sf_scalar_lightning_def(this, sf_scalar_lightning_def)
      !
      !  DESCRIPTION:
      !  sets value for sf_scalar_lightning_def
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this                    ! runtime params object
      integer,                   intent(in)    :: sf_scalar_lightning_def ! definition of a scalar-lightning case for spitfire_mode
      
      call this%sf_scalar_lightning_def%set_value(sf_scalar_lightning_def)

    end subroutine set_sf_scalar_lightning_def

    ! ------------------------------------------------------------------------------------

    integer function get_sf_scalar_lightning_def(this)
      !
      !  DESCRIPTION:
      !  gets sf_scalar_lightning_def value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_sf_scalar_lightning_def = this%sf_scalar_lightning_def%get_value()

    end function get_sf_scalar_lightning_def

    ! ------------------------------------------------------------------------------------
    
    subroutine set_sf_successful_ignitions_def(this, sf_successful_ignitions_def)
      !
      !  DESCRIPTION:
      !  sets value for sf_successful_ignitions_def
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this                        ! runtime params object
      integer,                   intent(in)    :: sf_successful_ignitions_def ! definition of a successful-ignition dataset case for spitfire_mode
      
      call this%sf_successful_ignitions_def%set_value(sf_successful_ignitions_def)

    end subroutine set_sf_successful_ignitions_def

    ! ------------------------------------------------------------------------------------

    integer function get_sf_successful_ignitions_def(this)
      !
      !  DESCRIPTION:
      !  gets sf_successful_ignitions_def value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_sf_successful_ignitions_def = this%sf_successful_ignitions_def%get_value()

    end function get_sf_successful_ignitions_def

    ! ------------------------------------------------------------------------------------
    
    subroutine set_sf_anthro_ignitions_def(this, sf_anthro_ignitions_def)
      !
      !  DESCRIPTION:
      !  sets value for sf_anthro_ignitions_def
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this                    ! runtime params object
      integer,                   intent(in)    :: sf_anthro_ignitions_def ! definition of an anthropogenic-ignition dataset case for spitfire_mode
      
      call this%sf_anthro_ignitions_def%set_value(sf_anthro_ignitions_def)

    end subroutine set_sf_anthro_ignitions_def

    ! ------------------------------------------------------------------------------------

    integer function get_sf_anthro_ignitions_def(this)
      !
      !  DESCRIPTION:
      !  gets sf_anthro_ignitions_def value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_sf_anthro_ignitions_def = this%sf_anthro_ignitions_def%get_value()

    end function get_sf_anthro_ignitions_def

    ! ------------------------------------------------------------------------------------
    
    subroutine set_use_logging(this, use_logging)
      !
      !  DESCRIPTION:
      !  sets value for use_logging
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this        ! runtime params object
      logical,                   intent(in)    :: use_logging ! whether or not to use the logging module
      
      call this%use_logging%set_value(use_logging)

    end subroutine set_use_logging

    ! ------------------------------------------------------------------------------------

    logical function get_use_logging(this)
      !
      !  DESCRIPTION:
      !  gets use_logging value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_use_logging = this%use_logging%get_value()

    end function get_use_logging

    ! ------------------------------------------------------------------------------------
    
    subroutine set_use_planthydro(this, use_planthydro)
      !
      !  DESCRIPTION:
      !  sets value for use_planthydro
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this           ! runtime params object
      logical,                   intent(in)    :: use_planthydro ! whether or not to use plant hydraulics 
      
      call this%use_planthydro%set_value(use_planthydro)
      if (use_planthydro) then 
          write(fates_log(), *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          write(fates_log(), *) ''
          write(fates_log(), *) ' use_planthydro is an      EXPERIMENTAL FEATURE        '
          write(fates_log(), *) ' please see header of fates/biogeophys/FatesHydraulicsMod.F90'
          write(fates_log(), *) ' for more information.'
          write(fates_log(), *) ''
          write(fates_log(), *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      end if 

    end subroutine set_use_planthydro

    ! ------------------------------------------------------------------------------------

    logical function get_use_planthydro(this)
      !
      !  DESCRIPTION:
      !  gets use_planthydro value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_use_planthydro = this%use_planthydro%get_value()

    end function get_use_planthydro

    ! ------------------------------------------------------------------------------------
    
    subroutine set_use_cohort_age_tracking(this, use_cohort_age_tracking)
      !
      !  DESCRIPTION:
      !  sets value for use_cohort_age_tracking
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this                    ! runtime params object
      logical,                   intent(in)    :: use_cohort_age_tracking ! whether or not to use cohort age tracking 
      
      call this%use_cohort_age_tracking%set_value(use_cohort_age_tracking)

    end subroutine set_use_cohort_age_tracking

    ! ------------------------------------------------------------------------------------

    logical function get_use_cohort_age_tracking(this)
      !
      !  DESCRIPTION:
      !  gets use_cohort_age_tracking value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_use_cohort_age_tracking = this%use_cohort_age_tracking%get_value()

    end function get_use_cohort_age_tracking

    ! ------------------------------------------------------------------------------------
    
    subroutine set_use_tree_damage(this, use_tree_damage)
      !
      !  DESCRIPTION:
      !  sets value for use_tree_damage
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this            ! runtime params object
      logical,                   intent(in)    :: use_tree_damage ! whether or not to turn on the tree damage model 
      
      call this%use_tree_damage%set_value(use_tree_damage)

    end subroutine set_use_tree_damage

    ! ------------------------------------------------------------------------------------

    logical function get_use_tree_damage(this)
      !
      !  DESCRIPTION:
      !  gets use_tree_damage value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_use_tree_damage = this%use_tree_damage%get_value()

    end function get_use_tree_damage

    ! ------------------------------------------------------------------------------------
    
    subroutine set_use_ed_st3(this, use_ed_st3)
      !
      !  DESCRIPTION:
      !  sets value for use_ed_st3
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this       ! runtime params object
      logical,                   intent(in)    :: use_ed_st3 ! whether or not to use ST(atic) ST(and) ST(ructure) mode (ST3) 
      
      call this%use_ed_st3%set_value(use_ed_st3)

    end subroutine set_use_ed_st3

    ! ------------------------------------------------------------------------------------

    logical function get_use_ed_st3(this)
      !
      !  DESCRIPTION:
      !  gets use_ed_st3 value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_use_ed_st3 = this%use_ed_st3%get_value()

    end function get_use_ed_st3

    ! ------------------------------------------------------------------------------------
    
    subroutine set_use_ed_prescribed_phys(this, use_ed_prescribed_phys)
      !
      !  DESCRIPTION:
      !  sets value for use_ed_prescribed_phys
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this                   ! runtime params object
      logical,                   intent(in)    :: use_ed_prescribed_phys ! whether or not to use prescribed physiology (i.e. opposite of ST3) 
      
      call this%use_ed_prescribed_phys%set_value(use_ed_prescribed_phys)

    end subroutine set_use_ed_prescribed_phys

    ! ------------------------------------------------------------------------------------

    logical function get_use_ed_prescribed_phys(this)
      !
      !  DESCRIPTION:
      !  gets use_ed_prescribed_phys value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_use_ed_prescribed_phys = this%use_ed_prescribed_phys%get_value()

    end function get_use_ed_prescribed_phys

    ! ------------------------------------------------------------------------------------
    
    subroutine set_use_inventory_init(this, use_inventory_init)
      !
      !  DESCRIPTION:
      !  sets value for use_inventory_init
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this               ! runtime params object
      logical,                   intent(in)    :: use_inventory_init ! whether or not to initialize simulation from inventory file 
      
      call this%use_inventory_init%set_value(use_inventory_init)

    end subroutine set_use_inventory_init

    ! ------------------------------------------------------------------------------------

    logical function get_use_inventory_init(this)
      !
      !  DESCRIPTION:
      !  gets use_inventory_init value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_use_inventory_init = this%use_inventory_init%get_value()

    end function get_use_inventory_init

    ! ------------------------------------------------------------------------------------
    
    subroutine set_inventory_ctrl_file(this, inventory_ctrl_file)
      !
      !  DESCRIPTION:
      !  sets value for inventory_ctrl_file
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this                ! runtime params object
      character(char_len),       intent(in)    :: inventory_ctrl_file ! full path to the inventory controle file that specifies available inventory
      
      call this%inventory_ctrl_file%set_value(inventory_ctrl_file)

    end subroutine set_inventory_ctrl_file

    ! ------------------------------------------------------------------------------------

    character(char_len) function get_inventory_ctrl_file(this)
      !
      !  DESCRIPTION:
      !  gets inventory_ctrl_file value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_inventory_ctrl_file = this%inventory_ctrl_file%get_value()

    end function get_inventory_ctrl_file

    ! ------------------------------------------------------------------------------------
    
    subroutine set_use_fixed_biogeog(this, use_fixed_biogeog)
      !
      !  DESCRIPTION:
      !  sets value for use_fixed_biogeog
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this              ! runtime params object
      logical,                   intent(in)    :: use_fixed_biogeog ! flag to use FATES fixed biogeography mode 
      
      call this%use_fixed_biogeog%set_value(use_fixed_biogeog)

    end subroutine set_use_fixed_biogeog

    ! ------------------------------------------------------------------------------------

    logical function get_use_fixed_biogeog(this)
      !
      !  DESCRIPTION:
      !  gets use_fixed_biogeog value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_use_fixed_biogeog = this%use_fixed_biogeog%get_value()

    end function get_use_fixed_biogeog

    ! ------------------------------------------------------------------------------------
    
    subroutine set_use_nocomp(this, use_nocomp)
      !
      !  DESCRIPTION:
      !  sets value for use_nocomp
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this       ! runtime params object
      logical,                   intent(in)    :: use_nocomp ! flag to use FATES no competition mode 
      
      call this%use_nocomp%set_value(use_nocomp)

    end subroutine set_use_nocomp

    ! ------------------------------------------------------------------------------------

    logical function get_use_nocomp(this)
      !
      !  DESCRIPTION:
      !  gets use_nocomp value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_use_nocomp = this%use_nocomp%get_value()

    end function get_use_nocomp

    ! ------------------------------------------------------------------------------------
    
    subroutine set_use_sp(this, use_sp)
      !
      !  DESCRIPTION:
      !  sets value for use_sp
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(inout) :: this   ! runtime params object
      logical,                   intent(in)    :: use_sp ! flag to use FATES satellite phenology (LAI) mode
      
      call this%use_sp%set_value(use_sp)

    end subroutine set_use_sp

    ! ------------------------------------------------------------------------------------

    logical function get_use_sp(this)
      !
      !  DESCRIPTION:
      !  gets use_sp value
      !

      ! ARGUMENTS:
      class(hlm_runtime_params), intent(in) :: this ! runtime params object

      get_use_sp = this%use_sp%get_value()

    end function get_use_sp

    ! ------------------------------------------------------------------------------------

end module FatesHLMRuntimeParamsMod