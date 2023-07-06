module FatesFinalParamType
  ! --------------------------------------------------------------------------------------
  ! This defines a type for holding, getting, and setting runtime parameters that should 
  ! only be set once
  ! --------------------------------------------------------------------------------------

  use FatesConstantsMod, only : fates_unset_r8
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesGlobals,      only : endrun => fates_endrun
  use FatesGlobals,      only : fates_log

  use shr_log_mod,       only : errMsg => shr_log_errMsg

  implicit none 
  private

  character(len=*), parameter         :: sourcefile = __FILE__ ! for error printing

  ! ======================================================================================

  type, public :: param

    logical, private :: is_set = .false. ! once .true. cannot write to parameter
                                         ! must be .true. to read parameter
    contains

      procedure          :: already_set
      procedure, private :: throw_if_bad_set 
      procedure, private :: throw_if_bad_get

  end type param

  ! ======================================================================================

  type, public, extends(param) :: real_param

    real(r8), private :: value 

    contains 
    
      procedure, public :: set_value => set_real_value
      procedure, public :: get_value => get_real_value

  end type real_param

  ! ======================================================================================

  type, public, extends(param) :: integer_param

  integer, private :: value 

  contains 
  
    procedure, public :: set_value => set_int_value
    procedure, public :: get_value => get_int_value

  end type integer_param

! ========================================================================================
  type, public, extends(param) :: logical_param

  logical, private :: value 

  contains 
  
    procedure, public  :: set_value => set_logical_value
    procedure, public  :: get_value => get_logical_value

  end type logical_param

! ========================================================================================

  type, public, extends(param) :: character_param

  character(:), allocatable, private :: value 

  contains 
  
    procedure, public :: set_value => set_char_value
    procedure, public :: get_value => get_char_value

  end type character_param

! ========================================================================================

  ! restrict access to the actual procedure names
  private :: set_real_value, get_real_value
  private :: set_int_value, get_int_value
  private :: set_char_value, get_char_value
  private :: set_logical_value, get_logical_value

  contains

    ! ====================================================================================

    subroutine set_real_value(this, val)
      !
      !  DESCRIPTION:
      !  set the value of this parameter
      !

      ! ARGUMENTS:
      class(real_param), intent(inout) :: this ! real param object
      real(r8),          intent(in)    :: val  ! value to set it to

      if (.not. this%already_set()) then 
        this%value = val
        this%is_set = .true.
      else 
        call this%throw_if_bad_set()
      end if 

    end subroutine set_real_value

    ! ------------------------------------------------------------------------------------

    real function get_real_value(this)
      !
      !  DESCRIPTION:
      !  get the value of this parameter
      !

      ! ARGUMENTS:
      class(real_param), intent(in) :: this ! real param object
    
      if (this%already_set()) then 
        get_real_value = this%value 
      else 
        call this%throw_if_bad_get()
      end if 

    end function get_real_value

    ! ------------------------------------------------------------------------------------

    subroutine set_int_value(this, val)
      !
      !  DESCRIPTION:
      !  set the value of this parameter
      !

      ! ARGUMENTS:
      class(integer_param), intent(inout) :: this ! integer param object
      integer,              intent(in)    :: val  ! value to set it to

      if (.not. this%already_set()) then 
        this%value = val
        this%is_set = .true.
      else 
        call this%throw_if_bad_set()
      end if 

    end subroutine set_int_value

    ! ------------------------------------------------------------------------------------

    integer function get_int_value(this)
      !
      !  DESCRIPTION:
      !  get the value of this parameter
      !

      ! ARGUMENTS:
      class(integer_param), intent(in) :: this ! integer param object
    
      if (this%already_set()) then 
        get_int_value = this%value 
      else 
        call this%throw_if_bad_get()
      end if 

    end function get_int_value

    ! ------------------------------------------------------------------------------------

    subroutine set_logical_value(this, val)
      !
      !  DESCRIPTION:
      !  set the value of this parameter
      !

      ! ARGUMENTS:
      class(logical_param), intent(inout) :: this ! logical param object
      logical,              intent(in)    :: val  ! value to set it to

      if (.not. this%already_set()) then 
        this%value = val
        this%is_set = .true.
      else 
        call this%throw_if_bad_set()
      end if 

    end subroutine set_logical_value

    ! ------------------------------------------------------------------------------------

    logical function get_logical_value(this)
      !
      !  DESCRIPTION:
      !  get the value of this parameter
      !

      ! ARGUMENTS:
      class(logical_param), intent(in) :: this ! logical param object
    
      if (this%already_set()) then 
        get_logical_value = this%value 
      else 
        call this%throw_if_bad_get()
      end if 

    end function get_logical_value

    ! ------------------------------------------------------------------------------------

    subroutine set_char_value(this, val)
      !
      !  DESCRIPTION:
      !  set the value of this parameter
      !

      ! ARGUMENTS:
      class(character_param), intent(inout) :: this ! character param object
      character(len=*),       intent(in)    :: val  ! value to set it to

      if (.not. this%already_set()) then 
        this%value = trim(val)
        this%is_set = .true.
      else 
        call this%throw_if_bad_set()
      end if 

    end subroutine set_char_value

    ! ------------------------------------------------------------------------------------

    function get_char_value(this) result(res)
      !
      !  DESCRIPTION:
      !  get the value of this parameter
      !

      ! ARGUMENTS:
      character(:), allocatable          :: res  ! result
      class(character_param), intent(in) :: this ! character param object
    
      if (this%already_set()) then 
        res = this%value 
      else 
        call this%throw_if_bad_get()
      end if 

    end function get_char_value

    ! ------------------------------------------------------------------------------------
    
    subroutine throw_if_bad_set(this)
      !
      !  DESCRIPTION:
      !  kills model if parameter is asked to be set after already set
      !

      ! ARGUMENTS:
      class(param), intent(in) :: this ! param object

      write(fates_log(),*) 'Attempt to set runtime parameter that has already been set!'
      call endrun(msg=errMsg(sourcefile, __LINE__))

    end subroutine throw_if_bad_set

    ! ------------------------------------------------------------------------------------

    subroutine throw_if_bad_get(this)
      !
      !  DESCRIPTION:
      !  kills model if parameter has not yet been set
      !

      ! ARGUMENTS:
      class(param), intent(in) :: this ! param object

      write(fates_log(),*) 'Attempt to get runtime parameter that has not been set yet!'
      call endrun(msg=errMsg(sourcefile, __LINE__))

    end subroutine throw_if_bad_get

    ! ------------------------------------------------------------------------------------

    logical function already_set(this)
    !
    !  DESCRIPTION:
    !  checks if this variable has been set yet
    !

    ! ARGUMENTS:
    class(param), intent(in) :: this ! param object

    already_set = this%is_set

    end function already_set

    ! ------------------------------------------------------------------------------------

end module FatesFinalParamType