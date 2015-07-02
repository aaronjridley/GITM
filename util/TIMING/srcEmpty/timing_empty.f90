!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================================
! Empty Timing module
!=============================================================================

subroutine timing_version(on,name,number)

  logical, intent(out) :: on
  character (len=40), intent(out) :: name
  real, intent(out) :: number

  on    =.false.
  name  ='TIMING EMPTY'
  number=0.0

end subroutine timing_version

!==============================================================================
subroutine timing_active(value)
  logical, intent(in) :: value

  if(value)write(*,*)'Warning: TIMING_EMPTY cannot be activated !!!' 
end subroutine timing_active

!==============================================================================
subroutine timing_step(value)
  integer, intent(in) :: value
end subroutine timing_step

!==============================================================================
subroutine timing_comp_proc(value1,value2)
  character (len=*), intent(in) :: value1
  integer, intent(in) :: value2
end subroutine timing_comp_proc

!==============================================================================
subroutine timing_iounit(value1)
  integer, intent(in) :: value1
end subroutine timing_iounit

!==============================================================================
subroutine timing_depth(value)
  integer, intent(in) :: value
end subroutine timing_depth

!==============================================================================
subroutine timing_report_style(value)
  character (LEN=*), intent(in) :: value
end subroutine timing_report_style

!==============================================================================
subroutine timing_param_put_i(name,value,error)
  character (LEN=*), intent(in) :: name
  integer, intent(in) :: value
  integer, intent(out):: error
end subroutine timing_param_put_i

!==============================================================================

function timing_func_d(func_name,iclock,name,parent_name)
  implicit none
  real(selected_real_kind(12)) timing_func_d
  character (LEN=*), intent(in):: func_name, name, parent_name
  integer, intent(in) :: iclock
  timing_func_d = -1
end function timing_func_d

!==============================================================================

subroutine timing_start(name)
  character (LEN=*), intent(in):: name
end subroutine timing_start

!==============================================================================
subroutine timing_stop(name)
  character (LEN=*), intent(in):: name
end subroutine timing_stop

!==============================================================================
subroutine timing_reset_all
  call timing_reset('#all',2)
end subroutine timing_reset_all

!==============================================================================
subroutine timing_reset(name,nclock)
  character (LEN=*), intent(in):: name
  integer, intent(in) :: nclock
end subroutine timing_reset

!==============================================================================
subroutine timing_show(name,iclock)
  character (LEN=*), intent(in):: name
  integer,           intent(in):: iclock
end subroutine timing_show

!==============================================================================
subroutine timing_report
end subroutine timing_report

!==============================================================================
subroutine timing_report_total
end subroutine timing_report_total

!==============================================================================
subroutine timing_tree(iclock,show_depth)
  integer, intent(in):: iclock, show_depth
end subroutine timing_tree

!==============================================================================
subroutine timing_sort(iclock,show_length,unique)
  integer, intent(in) :: iclock, show_length
  logical, intent(in) :: unique
end subroutine timing_sort

!==============================================================================


