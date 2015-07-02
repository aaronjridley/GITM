!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
!QUOTE: \clearpage
!
!BOP
!
!MODULE: ModFreq - frequency related methods
!
!DESCRIPTION:
! The module provides 
! a data type TypeFreq and the is\_time\_to function for an easy handling
! of actions done with some frequency. The frequency can be defined in
! terms of time steps or simulation time. The initial step or time can
! be also given, which is adjusted to the current step and time with
! the adjust\_freq subroutine.

!INTERFACE:
module ModFreq

  !USES:

  implicit none

  save

  !PUBLIC TYPES:

  ! Derived type to store data about the frequency of some action
  type FreqType
     logical           :: DoThis ! Do the action or not
     integer           :: Dn     ! Frequency in terms of time steps
     real              :: Dt     ! Frequency in terms of seconds
     integer           :: nNext  ! The next time step the action should be done
     real              :: tNext  ! The next time the action should be done
  end type FreqType

  !PUBLIC MEMBER FUNCTIONS:
  public :: adjust_freq      ! Adjust initial step and time to current values
  public :: check_freq       ! Check frequency settings for correctness
  public :: is_time_to       ! Returns true it is time to act and updates act

  !PUBLIC DATA MEMBERS:

  !REVISION HISTORY:
  ! 01Aug03 Aaron Ridley and G. Toth - initial implementation
  ! 22Aug03 G. Toth - added TypeFreq and is_time_to function
  ! 25Aug03 G. Toth - added adjust_freq subroutine
  ! 23Mar04 G. Toth - splitting CON_time into CON_time, ModTime, ModTimeFreq
  ! 19May08 G. Toth - added check_freq subroutine
  !EOP

  character(len=*), parameter, private :: NameMod='ModFreq'

contains

  !BOP ========================================================================
  !IROUTINE: check_freq - check frequency settings
  !INTERFACE:
  subroutine check_freq(NameAct, Act, DoTimeAccurate)

    !INPUT ARGUMENTS:
    character(len=*), intent(in) :: NameAct        ! Name of action
    type(FreqType),   intent(in) :: Act            ! Frequency of some action
    logical,          intent(in) :: DoTimeAccurate ! Time accurate run?

    !DESCRIPTION:
    ! Check if Dt > 0 in time accurate mode and Dn > 0 in steady state mode
    !EOP
    character(len=*), parameter :: NameSub = NameMod//'::check_freq'
    !------------------------------------------------------------------------

    if(.not. Act % DoThis) RETURN

    if(Act % Dt <= 0.0 .and. DoTimeAccurate) then

       write(*,*)NameSub,' ERROR: Dt=',Act % Dt,' for action ',NameAct
       call CON_stop(NameSub//' Dt must be positive in time accurate mode!')

    else if(Act % Dn <= 0.0 .and. .not. DoTimeAccurate) then

       write(*,*)NameSub,' ERROR: Dn=',Act % Dn,' for action ',NameAct
       call CON_stop(NameSub//' Dn must be positive in steady state mode!')

    end if

  end subroutine check_freq

  !BOP ========================================================================
  !IROUTINE: adjust_freq - adjust initial values to current values
  !INTERFACE:
  subroutine adjust_freq(Act, nStep, tSim, DoTimeAccurate)

    !INPUT/OUTPUT ARGUMENTS:
    type(FreqType), intent(inout) :: Act    ! Frequency of some action

    !INPUT ARGUMENTS:
    integer,           intent(in) :: nStep  ! Current step number
    real,              intent(in) :: tSim   ! Current simulation time
    logical,           intent(in) :: DoTimeAccurate ! Time accurate run?

    !DESCRIPTION:
    ! Adjust the nNext and tNext fields of Act based on the
    ! current time step nStep and current simulation time tSim.
    !EOP
    character(len=*), parameter :: NameSub = NameMod//'::adjust_freq'
    !------------------------------------------------------------------------

    if(.not. Act % DoThis) RETURN

    if(Act % Dn > 0 .and. Act % nNext < nStep) &
         Act % nNext = nStep - mod(nStep - Act % nNext, Act % Dn) + Act % Dn

    if(Act % Dt > 0 .and. Act % tNext < tSim ) &
         Act % tNext = tSim  - mod(tSim  - Act % tNext, Act % Dt) + Act % Dt

  end subroutine adjust_freq

  !BOP ========================================================================
  !IROUTINE: is_time_to - is it time to do something
  !INTERFACE:
  function is_time_to(Act,nStep,tSimulation,DoTimeAccurate,DoCheckOnly) &
       result(IsTimeTo)

    !INPUT/OUTPUT ARGUMENTS:
    type(FreqType), intent(inout) :: Act

    !INPUT ARGUMENTS:
    integer,           intent(in) :: nStep          ! Step number
    real,              intent(in) :: tSimulation    ! Simulation time
    logical,           intent(in) :: DoTimeAccurate ! Time accurate run?
    logical, optional, intent(in) :: DoCheckOnly    ! Only checking?

    !RETURN VALUE:
    logical :: IsTimeTo

    !DESCRIPTION:
    ! Based on the next step/time info in Act and the tSimulation and 
    ! nStep values (and the DoTimeAccurate variable) decide 
    ! if Act should be done. If the answer is yes and the optional
    ! parameter DoCheckOnly is not present, then the next step/time 
    ! values in Act are increased with the step/time frequencies.
    ! If the increased values do not reach the tSimulation and
    ! nStep values then increase the next step/time relative to the
    ! values of tSimulation and nStep.
    !EOP
    !------------------------------------------------------------------------
    if(.not.Act % DoThis)then
       IsTimeTo = .false.
       RETURN
    end if

    if(Act % Dn == 0)then
       IsTimeTo = .true.
    else
       IsTimeTo = .false.

       if(DoTimeAccurate)then
          if(Act % Dt >=0) then
             IsTimeTo = tSimulation >= Act % tNext
          elseif(Act % Dn > 0)then
             IsTimeTo = nStep >= Act % nNext
          end if
       else
          if(Act % Dn > 0)then
             IsTimeTo = nStep >= Act % nNext
          end if
       end if
    end if

    if(IsTimeTo .and. .not.present(DoCheckOnly))then

       if(Act % Dt > 0)then
          Act % tNext = Act % tNext + Act % Dt
          if(Act % tNext < tSimulation) &
               Act % tNext = tSimulation + Act % Dt
       end if

       if(Act % Dn > 0) then
          Act % nNext = Act % nNext + Act % Dn
          if(Act % nNext <= nStep) &
               Act % nNext = nStep + Act % Dn
       end if
    end if

  end function is_time_to

  !BOP ========================================================================
  !IROUTINE: test_freq - test the methods in this module
  !INTERFACE:
  subroutine test_freq

    !LOCAL VARIABLES:
    logical, parameter :: IsVerbose = .false.
    type(FreqType)     :: Act
    integer, parameter :: nAct = 10
    logical            :: DoAct(nAct)
    logical            :: DoTimeAccurate = .true.
    logical, parameter :: F=.false., T=.true.
    !EOP
    !-------------------------------------------------------------------------
    !BOC
    write(*,*)'Testing is_time_to function'

    DoTimeAccurate = .false.

    Act = FreqType(.false., 1,  1.0, 0, 0.0); call check_act
    if(any(DoAct))stop 'Error: DoThis=.false. should yield all false'

    Act = FreqType(.true., -1, -1.0, 0, 0.0); call check_act
    if(any(DoAct))stop 'Error: frequency -1 -1.0 should yield all false'

    Act = FreqType(.true., 0, -1.0, 2, 0.0); call check_act
    if(.not.all(DoAct))stop 'Error: frequencies 0 -1.0 should yield all true'

    Act = FreqType(.true., 2, -1.0, 2, 0.0); call check_act
    if( any(DoAct .neqv. (/F,T,F,T,F,T,F,T,F,T/)) ) &
       stop 'Error with frequencies 2 -1.0 start at 2'

    Act = FreqType(.true., 4, -1.0, -100, 0.0); call check_act
    if( any(DoAct .neqv. (/T,F,F,F,T,F,F,F,T,F/)) ) &
       stop 'Error with frequency 4 -1.0 start at -100'

    DoTimeAccurate = .true.

    Act = FreqType(.true., -1, 0.0, 0, 0.0); call check_act
    if(.not.all(DoAct))stop 'Error: frequency -1 0.0 should yield all true'

    Act = FreqType(.true., -1, 1.0, 0, 0.0); call check_act
    if(.not.all(DoAct))stop 'Error: frequency -1 1.0 should yield all true'

    Act = FreqType(.true., -1, 2.0, 0, 1.9); call check_act
    if( any(DoAct .neqv. (/F,T,F,T,F,T,F,T,F,T/)) ) &
       stop 'Error with frequency 2.0 start at 1.9'

    Act = FreqType(.true., -1, 4.0, 0, -100.0); call check_act
    if( any(DoAct .neqv. (/T,F,F,F,T,F,F,F,T,F/)) ) &
       stop 'Error with frequency 4.0 start at -100.0'

    write(*,*)'Testing adjust_freq subroutine'

    Act = FreqType(.true.,10,0.0,0,0.0)
    call adjust_freq(Act,15,0.0,DoTimeAccurate)
    if( Act % nNext /= 20) &
         stop 'Error with freq=10,0.0,0,0.0 adjust_freq(15,0.0)'


    Act = FreqType(.true.,10,0.0,3,0.0)
    call adjust_freq(Act,15,0.0,DoTimeAccurate)
    if( Act % nNext /= 23) &
         stop 'Error with freq=10,0.0,3,0.0 adjust_freq(15,0.0)'

    Act = FreqType(.true.,10,0.0,3,0.0)
    call adjust_freq(Act,0,0.0,DoTimeAccurate)
    if( Act % nNext /= 3) &
         stop 'Error with freq=10,0.0,0,0.0 adjust_freq(0,0.0)'

    write(*,*)'Successful'

    !EOC

  contains

    !======================================================================

    subroutine check_act
      integer :: iAct
      real    :: t
      t=0.0
      do iAct=1,nAct
         t=t+1.0           
         DoAct(iAct)=is_time_to(Act,iAct,t,DoTimeAccurate)
      end do
      if(IsVerbose)then
         write(*,*)'DoTimeAccurate=',DoTimeAccurate
         write(*,*)'Act  =',Act
         write(*,*)'DoAct=',DoAct
      end if
    end subroutine check_act

  end subroutine test_freq

  !============================================================================

end module ModFreq
