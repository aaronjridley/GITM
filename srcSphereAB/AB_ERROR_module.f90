!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!-*- F90 -*- so emacs thinks this is an f90 file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NAME: AB_ERROR_module
!
! PURPOSE: a module which implements the error reporting for AB.
!
! HISTORY:
!  1/24/01 Robert Oehmke: created
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module AB_ERROR_module
  implicit none

   integer, parameter, private :: max_msg_len=1024
   character (len=max_msg_len), private :: error_msg
   character (len=max_msg_len), private :: error_loc


   interface AB_ERROR_set
      module procedure AB_ERROR_set_s, AB_ERROR_set_sn
   end interface

contains

  subroutine AB_ERROR_set_s(func,msg)
    character (len=*), intent(in) :: func,msg

    error_msg=msg
    error_loc=func

  end subroutine AB_ERROR_set_s


  subroutine AB_ERROR_set_sn(func,msg,i)
    character (len=*), intent(in) :: func,msg
    integer, intent(in) :: i
    character (len=100) :: s_i
  
    ! convert integer to string
    write(s_i,*) i
    
    ! construct error message
    error_msg=msg//trim(s_i)

    ! record location
    error_loc=func
    
  end subroutine AB_ERROR_set_sn


  subroutine AB_ERROR_write()
    
    write(*,*) "Error: ",trim(error_msg)
    write(*,*) "In subroutine: ",trim(error_loc)
    
  end subroutine AB_ERROR_write

  

end module AB_ERROR_module


