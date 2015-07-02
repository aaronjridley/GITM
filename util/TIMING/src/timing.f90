!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
!BOP
!
!MODULE: TIMING LIBRARY - a library for timing and profiling F90 codes
!
!DESCRIPTION:
!
! This collection of subroutines facilitates simple and versatile
! timing and profiling of Fortran 90 codes.
! See the user manual for usage.
!
!REVISION HISTORY:
! 05/11/2001 G. Toth <gtoth@umich.edu> - initial version of TIMING
!            several bug fixes and improvements later on
!
!EOP
!
!BOP
!ROUTINE: timing_version - provide version of the TIMING
!INTERFACE:
subroutine timing_version(on,name,number)

  !OUTPUT ARGUMENTS:
  logical, intent(out)            :: on
  character (len=40), intent(out) :: name
  real, intent(out)               :: number
  !DESCRIPTION:
  ! Provide version information. 
  !EOP
  !BOC
  on    =.true.
  name  ='TIMING by G. Toth (2001)'
  number=1.2
  !EOC
end subroutine timing_version

!==============================================================================
!BOP
!ROUTINE: timing_active - activate timing on this processor
!INTERFACE:
subroutine timing_active(value)
  !USES:
  use ModTiming, ONLY: UseTiming
  implicit none
  !INPUT ARGUMENTS:
  logical, intent(in) :: value
  !DESCRIPTION:
  ! Only the active processors do timing. If a processor is not active,
  ! all other timing subroutines return without any action, while timing
  ! functions will return some impossible values. Usually only one processor
  ! does the timing, but it is possible to do parallel timing on multiple
  ! processors. The timing reports can be distinguished by setting the
  ! component name and processor number with the timing\_comp\_proc subroutine.
  !EOP
  !BOC ------------------------------------------------------------------------
  UseTiming=value
  !EOC
end subroutine timing_active

!==============================================================================
!BOP
!ROUTINE: timing_step - set the current step number 
!INTERFACE:
subroutine timing_step(value)
  !USES:
  use ModTiming, ONLY: step
  implicit none
  !INPUT ARGUMENTS:
  integer, intent(in) :: value
  !DESCRIPTION:
  ! Information about number of calls/step and CPU time/step is based on
  ! the (increasing) step number provided to the timing module.
  ! The step number should be initialized (usually to 0) 
  ! at the beginning of the execution and then
  ! set at the beginning of every time step/iteration.
  !EOP
  !BOC -----------------------------------------------------------------------
  step = value
  !EOC
end subroutine timing_step

!==============================================================================
!BOP
!ROUTINE: timing_comp_proc - set component name and processor rank
!INTERFACE:
subroutine timing_comp_proc(value1,value2)
  !USES:
  use ModTiming, ONLY: NameComp, iProc
  implicit none
  !INPUT ARGUMENTS:
  character (len=*), intent(in) :: value1
  integer, intent(in) :: value2
  !DESCRIPTION:
  ! Set the name of the component and the processor number to be used
  ! in the timing reports.
  !EOP
  !BOC ------------------------------------------------------------------------
  NameComp = value1
  iProc    = value2
  !EOC
end subroutine timing_comp_proc

!==============================================================================
!BOP
!ROUTINE: timing_iounit - set I/O unit for normal output
!INTERFACE:
subroutine timing_iounit(value1)
  !USES:
  use ModTiming, ONLY: iUnit
  implicit none
  !INPUT ARGUMENTS:
  integer, intent(in) :: value1
  !DESCRIPTION:
  ! Set the I/O unit for normal output for this processor.
  !EOP
  !BOC ------------------------------------------------------------------------
  iUnit = value1
  !EOC
end subroutine timing_iounit

!==============================================================================
!BOP
!ROUTINE: timing_depth - set the maximum depth in the timing tree
!INTERFACE:
subroutine timing_depth(value)
  !USES:
  use ModTiming, ONLY: max_depth
  implicit none
  !INPUT ARGUMENTS:
  integer, intent(in) :: value
  !DESCRIPTION:
  ! Timing is done in a nested tree. The maximum level of nesting
  ! can be set to elinimate low level timings. A negative value
  ! indicates that there is no limit in the timing depth.
  !EOP
  !BOC ------------------------------------------------------------------------
  max_depth = value
  !EOC
end subroutine timing_depth

!==============================================================================
!BOP
!ROUTINE: timing_report_style - set the style of the timing report
!INTERFACE:
subroutine timing_report_style(value)
  !USES:
  use ModTiming, ONLY: report_style
  implicit none
  !INPUT ARGUMENTS:
  character (LEN=*), intent(in) :: value
  !DESCRIPTION:
  ! Timing reports can be produced in three formats. The 'tree' format
  ! shows the nested timing tree as an indented list of timings,
  ! the 'list' format shows the timings as a linear list with the
  ! sorted by the CPU time (timings called from different places are
  ! distinguished by the name of their 'parent' timings).
  ! Finally 'cumu' results in a cumulative timing, where 
  ! all timings with the same name are added up and sorted.
  !EOP
  !BOC ------------------------------------------------------------------------
  report_style = value
  !EOC
end subroutine timing_report_style

!==============================================================================
!BOP
!ROUTINE: timing_param_put_i - set an integer parameter in ModTiming
!INTERFACE:
subroutine timing_param_put_i(name,value,error)
  !USES:
  use ModTiming
  implicit none
  !INPUT ARGUMENTS:
  character (LEN=*), intent(in) :: name
  integer, intent(in) :: value
  !OUTPUT ARGUMENTS:
  integer, intent(out):: error
  !DESCRIPTION:
  ! This is a general subroutine to set an integer parameter in ModTiming.
  ! The parameter is identifed by its 'name' and set to 'value'. A -1 value
  ! is returned in 'error' if the 'name' is not recognized.
  !EOP
  !BOC ------------------------------------------------------------------------
  error=0
  select case(name)
  case('step')
     step=value
  case('depth')
     max_depth=value
  case('verbose')
     lVerbose=value
  case default
     error=-1
  end select
  !EOC
end subroutine timing_param_put_i

!==============================================================================
!BOP
!ROUTINE: timing_func_d - get a double precision value from ModTiming
!INTERFACE:
function timing_func_d(func_name,iclock,name,parent_name)
  !USES:
  use ModTiming
  implicit none
  !RETURN VALUE:
  real(Real8_) :: timing_func_d
  !INPUT ARGUMENTS:
  character (LEN=*), intent(in):: func_name, name, parent_name
  integer, intent(in) :: iclock
  !DESCRIPTION:
  !  The first string argument 'func\_name' determines the function
  ! to be returned. The available values are 'sum', 'sum/iter' and 'sum/call'.
  ! The latter two functions only make sense for clocks 2 and 3.
  ! The second integer argument 'iclock' selects the clock.
  ! The third string argument 'name' selects the timing,
  ! which is further specified by the last string argument 'parent\_name'.
  ! The parent is the timing that was started last but not stopped when
  ! the timing for 'name' is started. The parent of the first timing
  ! is itself so 
  !\begin{verbatim}
  !  write(*,*)'Elapsed time=',timing_func_d('sum',1,'main','main')
  !\end{verbatim}
  ! prints out the total time spent by 'main' since the last reset.
  ! The parent is needed to distinguish between timings called from
  ! different places. If the names do not match, no output is produced.
  !EOP
  integer     :: i, qclock, qiter, qcall
  real(Real8_) :: qsum
  !----------------------------------------------------------------------------

  ! write(*,*)'func_name=',func_name,' name=',name,' parent_name=',parent_name

  if(.not.UseTiming)then
     timing_func_d=-1.0
     RETURN
  end if

  do i=1,ntiming
     if(sa_name(i)==name .and. sa_name(ia_parent(i))==parent_name)EXIT
  end do

  if(i>ntiming)then
     timing_func_d=-1.0
     RETURN
  end if

  qclock = min( max(iclock,1), maxclock)

  if(la_active(i))then
     qsum = da_sum(i,qclock) + timing_cpu() - da_start(i)
  else
     qsum = da_sum(i,qclock)
  end if

  select case(func_name)
  case('sum')
     timing_func_d=qsum
  case('sum/iter')
     qiter=ia_iter(i,qclock)
     if(qiter<1)qiter=-1
     timing_func_d=qsum/qiter
  case('sum/call')
     qcall=ia_call(i,qclock)
     if(qcall<1)qcall=-1
     timing_func_d=qsum/qcall
  case default
     timing_func_d=-1.0
  end select

end function timing_func_d

!==============================================================================
!BOP
!ROUTINE: timing_start - start timing for 'name'
!INTERFACE:
subroutine timing_start(name)
  !USES
  use ModTiming
  implicit none
  !INPUT ARGUMENTS:
  character (LEN=*), intent(in):: name
  !DESCRIPTION:
  ! Start timing for 'name', where 'name' is an arbitrary string.
  !EOP
  integer :: i
  !----------------------------------------------------------------------------
  if(.not.UseTiming)RETURN

  current_depth = current_depth + 1

  if(max_depth >= 0 .and. current_depth > max_depth) RETURN

  if(lVerbose>2)write(*,*)'timing_start for ',name

  ! Search for previous timings of the same name at the same depth
  do i = i_last+1, ntiming
     if(ia_depth(i)==current_depth .and. &
          sa_name(i)==name) goto 100
  end do
  ! New name
  if(ntiming==maxtiming-1)write(*,*) &
       'WARNING: number of timings has reached maxtiming in ModTiming'

  if(ntiming==maxtiming) RETURN ! Cannot add more timing

  ntiming=ntiming+1
  i=ntiming
  sa_name(i)    = name
  ia_step(i)    = -1

100  continue

  ia_call(i,2:maxclock) = ia_call(i,2:maxclock)+1
  if(ia_step(i) < step) &
       ia_iter(i,2:maxclock) = ia_iter(i,2:maxclock)+1
  ia_step(i)     = step
  ia_depth(i)    = current_depth
  ia_parent(i)   = i_last
  la_active(i)   = .true.
  da_start(i)    = timing_cpu()
  i_last         = i

  if(lVerbose>2)write(*,*)'index, start:',i,da_start(i)
  
end subroutine timing_start

!==============================================================================
!BOP
!ROUTINE: timing_stop - stop timing for name
!INTERFACE:
subroutine timing_stop(name)
  !USES:
  use ModTiming
  implicit none
  !INPUT ARGUMENTS:
  character (LEN=*), intent(in):: name
  !DESCRITPTION:
  ! Stop timing for 'name'. If 'name' does not match the 'name' used
  ! in the last timing\_start call, a warning message is produced.
  !EOP
  integer     :: i
  real(Real8_) :: qnow
  !----------------------------------------------------------------------------
  if(.not.UseTiming)RETURN

  current_depth = current_depth - 1

  if(max_depth >= 0 .and. current_depth >= max_depth) RETURN

  if(lVerbose>2)write(*,*)'timing_stop  for ',name

  ! Write a warning for non-matching name
  if(sa_name(i_last)/=name.and.ntiming<maxtiming)write(*,*) &
       'WARNING in timing: unexpected STOP requested for ',name

  i=i_last
  qnow   = timing_cpu()
  da_sum(i,1)          = qnow-da_start(i)
  da_sum(i,2:maxclock) = da_sum(i,2:maxclock) + da_sum(i,1)
  la_active(i)         = .false.
  i_last               = ia_parent(i)

  if(lVerbose>2)then
     write(*,*)'index, stop :',i,qnow
     write(*,*)'clocks:',da_sum(i,:)
  end if

end subroutine timing_stop

!==============================================================================
!BOP
!ROUTINE: timing_reset_all - reset all clocks
!INTERFACE:
subroutine timing_reset_all
  !DESCRIPTION:
  ! See timing\_reset.
  !EOP
  !BOC
  call timing_reset('#all',2)
  !EOC
end subroutine timing_reset_all

!==============================================================================
!BOP
!ROUTINE: timing_reset - reset some clocks for some names
!INTERFACE:
subroutine timing_reset(name,nclock)
  !USES:
  use ModTiming
  implicit none
  !INPUT ARGUMENT:
  character (LEN=*), intent(in):: name
  integer, intent(in) :: nclock
  !DESCRIPTION:
  ! Reset clocks 1 to 'nclock' for 'name'.
  ! If the 'name' argument is set to '\#all',
  ! then the clocks 1 to 'nclock' are reset for all names.
  !EOP
  real(Real8_) :: qnow
  integer     :: i, qclock
  !----------------------------------------------------------------------------
  if(.not.UseTiming)RETURN

  if(lVerbose>1)write(*,*)'timing_reset 1..nclock=',nclock

  qclock = min(max(nclock,1),maxclock)

  qnow = timing_cpu()

  do i=1,ntiming
     if(name/='#all'.and.sa_name(i)/=name)CYCLE

     if(la_active(i))then
        da_sum(i,nclock+1:maxclock) = da_sum(i,nclock+1:maxclock) &
             + qnow-da_start(i)
        da_start(i) = qnow
        ia_step(i)  = step
        ia_iter(i,2:nclock) = 1
        ia_call(i,2:nclock) = 1
        da_sum(i,1:nclock)  = 0.0
     else
        ia_step(i) = -1
        ia_iter(i,1:nclock) = 0
        ia_call(i,1:nclock) = 0
        da_sum(i,1:nclock)  = 0.0
     end if
  end do

  if(name=='#all')step_reset(2:nclock)=step

end subroutine timing_reset

!==============================================================================
!BOP
!ROUTINE: timing_show - show timing by a clock for a name
!INTERFACE:
subroutine timing_show(name,iclock)
  !USES:
  use ModTiming
  implicit none
  !INPUT ARGUMENTS:
  character (LEN=*), intent(in):: name
  integer,           intent(in):: iclock
  !DESCRIPTION:
  ! Report timing results by clock iclock for name. 
  ! For clock 1 the  name, the calling parent, and the very last timing 
  ! are shown.
  ! For clock 2 the cumulative timing since the last reset is given.
  ! All timings matching the name (but called from different parents) 
  ! are shown. The timing per iteration and per call and the percentage
  ! with respect to the parent are also shown.
  ! For clock 3 the total timing is reported with the same information
  ! as for clock 2.
  !EOP
  integer     :: i, i_parent, qclock, qiter, qcall
  real(Real8_) :: qnow, qsum, qsumparent
  character (LEN=40) :: s_parent
  !----------------------------------------------------------------------------
  if(.not.UseTiming)RETURN

  qclock = min(max(1,iclock),maxclock)

  do i=ntiming,1,-1
     if(sa_name(i)/=name)CYCLE

     i_parent=ia_parent(i)
     s_parent=sa_name(i_parent)
     qnow=timing_cpu()
     if(la_active(i))then
        qsum = da_sum(i,qclock) + qnow-da_start(i)
     else
        qsum = da_sum(i,qclock)
     end if

     if(qclock==1)then
        write(iUnit,'(5a,f8.2,a)')                            &
                'Timing for last ',name,                  &
                ' (',s_parent(1:len_trim(s_parent)),'):', &
                qsum,' sec'
        RETURN
     end if

     if(la_active(i_parent))then
        qsumparent = da_sum(i_parent,qclock) + qnow-da_start(i_parent)
     else
        qsumparent = da_sum(i_parent,qclock)
     end if
     if(qsumparent<=0.0)qsumparent=-1.0

     qiter=ia_iter(i,qclock); if(qiter<1)qiter=-1
     qcall=ia_call(i,qclock); if(qcall<1)qcall=-1

     write(iUnit,'(2a)',ADVANCE='NO')'Timing for ',name
     if(qclock==maxclock .or. step_reset(qclock)==-1)then
        write(iUnit,'(a,i8,a)') ' at step',step,' :'
     else
        write(iUnit,'(a,i8,a,i8,a)')&
             ' from step',step_reset(qclock),' to',step,' :'
     end if
     write(iUnit,'(f8.2,a,f8.3,a,f8.3,a,f8.2,2a)')           &
                qsum,' sec, ',                            &
                qsum/qiter,' s/iter',                     &
                qsum/qcall,' s/call',                     &
                100*qsum/qsumparent,' % of ',s_parent(1:len_trim(s_parent))
  end do

end subroutine timing_show

!==============================================================================
!BOP
!ROUTINE: timing_report - report timings since last reset
!INTERFACE:
subroutine timing_report
  !USES:
  use ModTiming, ONLY: report_style
  !DESCRIPTION:
  ! Report timings since the last reset using the style
  ! set by timing\_report\_style.
  !EOP
  !BOC -----------------------------------------------------------------------
  select case(report_style)
  case('tree')
     call timing_tree(2,-1)
  case('list')
     call timing_sort(2,-1,.false.)
  case default
     call timing_sort(2,-1,.true.)
  end select
  !EOC
end subroutine timing_report

!==============================================================================
!BOP
!ROUTINE: timing_report_total - report total timings
!INTERFACE:
subroutine timing_report_total
  !USES:
  use ModTiming, ONLY: report_style
  !DESCRIPTION:
  ! Report total timings using the style
  ! set by timing\_report\_style.
  !EOP
  !BOC -----------------------------------------------------------------------
  select case(report_style)
  case('tree')
     call timing_tree(3,-1)
  case('list')
     call timing_sort(3,-1,.false.)
  case default
     call timing_sort(3,-1,.true.)
  end select
  !EOC
end subroutine timing_report_total

!==============================================================================
!BOP
!ROUTINE: timing_tree - show timings as a nested tree
!INTERFACE:
subroutine timing_tree(iclock,show_depth)
  !USES:
  use ModTiming
  use ModIoUnit, ONLY: io_unit_new
  implicit none
  !INPUT ARGUMENTS:
  integer, intent(in):: iclock, show_depth
  !DESCRIPTION:
  ! Produce a timing report that indicates calling sequence with order
  ! and nesting of calls with indentation.
  ! Use clock iclock.
  ! If show\_depth is positive, show only the timings which are 
  ! at a nesting level not exceeding show\_depth.
  !EOP
  ! Indices
  integer :: i, i_parent

  ! Temporary variables
  logical :: DoSaveFile = .false.
  integer :: iUnitT,iUnitT0
  character*80 :: filename
  integer     :: qclock, qdepth, i_depth, qiter, qcall, indent
  real(Real8_) :: qsum, qsumparent, qnow

  !----------------------------------------------------------------------------
  if(.not.UseTiming)RETURN

  if(DoSaveFile)then
     if(iProc==0)then
        iUnitT0 = io_unit_new()
        write(filename,'(a,i8.8,a)')'STDOUT/timing_it',step,'.H'
        open(UNIT=iUnitT0,FILE=trim(filename),STATUS='unknown',FORM='formatted')
        write(iUnitT0,'(a,i8,a)') 'TITLE = "Timing at it=',step,'"'
        write(iUnitT0,'(a)') 'VARIABLES = '
        write(iUnitT0,'(a)') '  "PE"'
     end if
     iUnitT = io_unit_new()
     write(filename,'(a,i8.8,a,i6.6,a)')'STDOUT/timing_it',step,'_pe',iProc,'.tec'
     open(UNIT=iUnitT,FILE=trim(filename),STATUS='unknown',FORM='formatted')
     write(iUnitT,*) iProc
  end if

  qclock = min(max(iclock,2),maxclock)

  write(iUnit,'(a79)') sepline
  write(iUnit,'(a)',ADVANCE='NO')'TIMING TREE'
  if(show_depth>0)&
       write(iUnit,'(a,i2)',ADVANCE='NO')' of depth',show_depth
  if(step_reset(qclock)>=0)then
     write(iUnit,'(a,i8,a,i8)',ADVANCE='NO') &
          ' from step',step_reset(qclock),' to',step
  else
     write(iUnit,'(a,i8)',ADVANCE='NO')' at step',step
  end if
  write(iUnit,'(a,i4)')' '//NameComp//' on PE ',iProc

  write(iUnit,'(a20,a7,a8,a9,a9,a9,a9)') &
       'name'//spaces,'#iter','#calls','sec','s/iter','s/call','percent'
  write(iUnit,'(a79)')sepline
  qdepth = 0
  qnow   = timing_cpu()
  do i=1,ntiming
     
     if(show_depth>0 .and. ia_depth(i)>show_depth)CYCLE

     if(la_active(i))then
        da_sum(i,:) = da_sum(i,:) + qnow - da_start(i)
        da_start(i) = qnow
     endif
     qsum=da_sum(i,qclock)
     if(DoSaveFile)then
        if(ia_iter(i,qclock)<1)CYCLE
     else
        if(qsum < 0.0005)CYCLE
     end if

     qsumparent=da_sum(ia_parent(i),qclock)
     if(qsumparent < qsum)qsumparent=qsum
     qiter=ia_iter(i,qclock); if(qiter<1)qiter=-1
     qcall=ia_call(i,qclock); if(qcall<1)qcall=-1
     ! Negative indent results in error because repeat(' ',-2) fails
     ! and spaces(1:-2) (although correct F90) also fails due to a
     ! PGF90 compiler bug
     indent=max(0,ia_depth(i)*2-4)
     write(iUnit,'(a20,i7,i8,a,f9.2,f9.3,f9.3,f9.2)')                    &
          repeat(' ',indent)//sa_name(i),                            &
          qiter,                                                     &
          qcall,                                                     &
          repeat(' ',indent),                                        &
          qsum,                                                      &
          qsum/qiter,                                                &
          qsum/qcall,                                                &
          100*qsum/qsumparent
     if(DoSaveFile)then
        if(iProc==0)then
           write(iUnitT0,'(a)') '  ".'//repeat(' ',indent)//trim(sa_name(i))//'"'
        end if
        write(iUnitT,*) qsum
     end if

     if(ia_depth(i)==1) write(iUnit,'(a79)')sepline
     
     ! Add up times for this depth and report missing part

     i_depth=ia_depth(i)
     if(i_depth > qdepth)then
        da_sum_other(i_depth) = qsum
     else
        da_sum_other(i_depth) = da_sum_other(i_depth) + qsum
     end if

     i_parent=i
     do qdepth = i_depth, ia_depth(i+1) + 1, -1

        ! Find parent (of parent (of parent...))
        i_parent=ia_parent(i_parent)
        
        ! Unmeasured timing = parent timing - sum of same level timing:
        qsum       = da_sum(i_parent,qclock)-da_sum_other(qdepth)

        ! reset sum for this depth
        da_sum_other(qdepth) = 0.0

        if(qsum<0.001) CYCLE
        qsumparent = da_sum(i_parent,qclock)
        if(qsumparent < qsum)qsumparent=qsum
        qiter      = ia_iter(i_parent,qclock)
        if(qiter<1)qiter=-1

        indent = max(0,qdepth*2-4)
        write(iUnit,'(a,a35,f9.2,f9.3,f18.2)')    &
             repeat(' ',indent),              &
             '#others'//spaces,               &
             qsum,                            &
             qsum/qiter,                      &
             100*qsum/qsumparent
        if(DoSaveFile)then
           if(iProc==0)then
              write(iUnitT0,'(a)') '  ".'//repeat(' ',indent)//'#others'//'"'
           end if
           write(iUnitT,*) qsum
        end if

     end do
     qdepth = i_depth
  end do
  write(iUnit,'(a79)')sepline

  if(DoSaveFile)then
     write(iUnitT0,'(a)') 'ZONE T="timing"'
     write(iUnitT0,'(a)') ' I=NPEs, J=1, K=1, ZONETYPE=Ordered, DATAPACKING=POINT'
     close(iUnitT0)
     close(iUnitT)
  end if

end subroutine timing_tree

!==============================================================================
!BOP
!ROUTINE: timing_sort - show timings as sorted list
!INTERFACE:
subroutine timing_sort(iclock,show_length,unique)
  !USES:
  use ModTiming
  implicit none
  !INPUT ARGUMENTS:
  integer, intent(in) :: iclock, show_length
  logical, intent(in) :: unique
  !DESCRIPTION:
  ! Make a sorted report of timings made with clock 'iclock'.
  ! If 'show\_length' is positive, only that many items are listed.
  ! If 'unique' is true, all calls with the same name are added up.
  !EOP
  integer :: ia_sort(maxtiming) ! indirect index array for sorting
  integer :: i, j, k, qclock
  real(Real8_)  :: qsum, qsummax, qnow
  character (LEN=40) :: s_name, s_parent

  integer     :: qntiming, showntiming
  real(Real8_) :: da_qsum(maxtiming)
  integer     :: ia_qcall(maxtiming), ia_qiter(maxtiming)
  character (LEN=40) :: sa_qname(maxtiming)
  
  !----------------------------------------------------------------------------
  if(.not.UseTiming)RETURN

  qclock = min(max(iclock,1),maxclock)

  qnow=timing_cpu()
  if(unique)then
     ! Add up timing and number of iterations and calls for identical names
     da_qsum(1:ntiming) =0.0
     ia_qiter(1:ntiming)=0
     ia_qcall(1:ntiming)=0
     qntiming=0
     do i=1,ntiming
        ! Check if the same name occured before or not
        do j=1,qntiming
           if(sa_name(i)==sa_qname(j))EXIT
        end do
        sa_qname(j) = sa_name(i)
        if(la_active(i))then
           qsum = da_sum(i,qclock) + qnow - da_start(i)
        else
           qsum = da_sum(i,qclock)
        end if
        da_qsum(j)  = da_qsum(j)  + qsum
        ia_qiter(j) = ia_qiter(j) + ia_iter(i,qclock)
        ia_qcall(j) = ia_qcall(j) + ia_call(i,qclock)
        if (j > qntiming) qntiming = j
     end do
  else
     qntiming = ntiming
     sa_qname(1:ntiming) = sa_name(1:ntiming)
     where(la_active(1:ntiming))
        da_qsum(1:ntiming) = da_sum(1:ntiming,qclock) &
                           + qnow - da_start(1:ntiming)
     elsewhere
        da_qsum(1:ntiming) = da_sum(1:ntiming,qclock)
     end where
     ia_qiter(1:ntiming) = ia_iter(1:ntiming,qclock)
     ia_qcall(1:ntiming) = ia_call(1:ntiming,qclock)
  end if

  forall(i=1:qntiming)ia_sort(i)=i

  do i=1,qntiming-1
     do j=i+1,qntiming
        if(da_qsum(ia_sort(i)) < da_qsum(ia_sort(j)))then
           k=ia_sort(i)
           ia_sort(i)=ia_sort(j)
           ia_sort(j)=k
        end if
     end do
  end do

  qsummax=da_qsum(ia_sort(1))

  if(qsummax<=0.0)then
     write(*,*)'WARNING in timing_sort: Maximum timing is <= 0 !!!'
     RETURN
  end if

  write(iUnit,'(a79)')sepline

  write(iUnit,'(a,i8,a,i2)',ADVANCE='NO')'SORTED TIMING'
  if(show_length>0)&
       write(iUnit,'(a,i3)',ADVANCE='NO')' of length',show_length
  if(qclock>1 .and. step_reset(qclock)>=0)then
     write(iUnit,'(a,i8,a,i8)',ADVANCE='NO') &
          ' from step',step_reset(qclock),' to',step
  else
     write(iUnit,'(a,i8)',ADVANCE='NO')' at step',step
  end if
  write(iUnit,'(a,i4)')' '//NameComp//' on PE ',iProc

  write(iUnit,'(a20)',ADVANCE='NO')                'name'//spaces
  if(.not.unique)write(iUnit,'(a20)',ADVANCE='NO') '(parent)'//spaces
  write(iUnit,'(a10,a10)',ADVANCE='NO')              'sec','percent'
  if(qclock>1)write(iUnit,'(a10,a10)',ADVANCE='NO')  '#iter','#calls'
  write(iUnit,*)
  write(iUnit,'(a79)')sepline

  if(show_length>0)then
     showntiming=min(show_length,qntiming)
  else
     showntiming=qntiming
  endif

  do k=1,showntiming
     i=ia_sort(k)
     s_name  =sa_qname(i)
     qsum    =da_qsum(i)

     if(qsum < 0.001) CYCLE

     write(iUnit,'(a20)',ADVANCE='NO')       s_name

     if(.not.unique)then
        s_parent=sa_name(ia_parent(i))
        write(iUnit,'(a20)',ADVANCE='NO') &
          '('//s_parent(1:len_trim(s_parent))//')'//spaces
     end if

     write(iUnit,'(f10.2,f10.2)',ADVANCE='NO') qsum, 100*qsum/qsummax

     if(qclock>1)write(iUnit,'(i10,i10)',ADVANCE='NO') ia_qiter(i), ia_qcall(i)
     write(iUnit,*)
     if(ia_depth(i)==1) write(iUnit,'(a79)')sepline

  end do

  qsum=0
  do k=showntiming+1,qntiming
     i=ia_sort(k)
     qsum=qsum+da_qsum(i)
  end do
  if(qsum>0.0)then
     write(iUnit,'(a20)',ADVANCE='NO')'#others'//spaces
     if(.not.unique)write(iUnit,'(a20)',ADVANCE='NO')' '
     write(iUnit,'(f8.2,f8.2)')qsum, 100*qsum/qsummax
  end if

  write(iUnit,'(a79)')sepline

end subroutine timing_sort
