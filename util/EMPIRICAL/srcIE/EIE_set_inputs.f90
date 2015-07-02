!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine EIE_set_inputs(StringInputLines)

  use ModEIE_Interface
  use ModEIEFiles
  use ModExtras

  implicit none

  character (len=100), dimension(*), intent(in) :: StringInputLines
  character (len=100) :: StringLine
  logical :: IsDone
  integer :: iLine, iDebugProc

  iLine = 1

  IsDone = .false.

  do while (.not.IsDone)

     StringLine = StringInputLines(iLine)

     if (StringLine(1:1) == "#") then

        if (index(StringLine,"#BACKGROUND") > 0) then

           call read_in_string(EIE_NameOfModelDir)
           call read_in_string(EIE_NameOfEFieldModel)
           call read_in_string(EIE_NameOfAuroralModel)
           call read_in_string(EIE_NameOfSolarModel)

           if (index(EIE_NameOfAuroralModel,'IHP') > 0) &
                EIE_NameOfAuroralModel = 'ihp'
           if (index(EIE_NameOfAuroralModel,'PEM') > 0) &
                EIE_NameOfAuroralModel = 'pem'

           if (index(EIE_NameOfEFieldModel,'AMIE') > 0) &
                EIE_NameOfEFieldModel = 'amie'

           if (index(EIE_NameOfEFieldModel,'weimer01') > 0) &
                EIE_NameOfEFieldModel = 'weimer01'
           if (index(EIE_NameOfEFieldModel,'Weimer01') > 0) &
                EIE_NameOfEFieldModel = 'weimer01'
           if (index(EIE_NameOfEFieldModel,'WEIMER01') > 0) &
                EIE_NameOfEFieldModel = 'weimer01'

           if (index(EIE_NameOfEFieldModel,'weimer05') > 0) &
                EIE_NameOfEFieldModel = 'weimer05'
           if (index(EIE_NameOfEFieldModel,'Weimer05') > 0) &
                EIE_NameOfEFieldModel = 'weimer05'
           if (index(EIE_NameOfEFieldModel,'WEIMER05') > 0) &
                EIE_NameOfEFieldModel = 'weimer05'

!           if (index(EIE_NameOfEFieldModel,'weimer') > 0 .and. &
!                index(EIE_NameOfEFieldModel,'01') == 0) &
!                EIE_NameOfEFieldModel = 'weimer96'
!           if (index(EIE_NameOfEFieldModel,'Weimer') > 0 .and. &
!                index(EIE_NameOfEFieldModel,'01') == 0) &
!                EIE_NameOfEFieldModel = 'weimer96'
!           if (index(EIE_NameOfEFieldModel,'WEIMER') > 0 .and. &
!                index(EIE_NameOfEFieldModel,'01') == 0) &
!                EIE_NameOfEFieldModel = 'weimer96'

           if (index(EIE_NameOfEFieldModel,'weimer96') > 0) &
                EIE_NameOfEFieldModel = 'weimer96'
           if (index(EIE_NameOfEFieldModel,'Weimer96') > 0) &
                EIE_NameOfEFieldModel = 'weimer96'
           if (index(EIE_NameOfEFieldModel,'WEIMER96') > 0) &
                EIE_NameOfEFieldModel = 'weimer96'

           if (index(EIE_NameOfEFieldModel,'SAMIE') > 0) &
                EIE_NameOfEFieldModel = 'samie'

           UseGridBasedEIE = .false.
           UAl_UseGridBasedEIE = .false.

        endif

        if (index(StringLine,"#AMIEFILES") > 0) then
           call read_in_string(AMIEFileNorth)
           call read_in_string(AMIEFileSouth)
           EIE_NameOfEFieldModel = "amie"
           UseGridBasedEIE = .true.
           UAl_UseGridBasedEIE = .true.
        endif

        if (index(StringLine,"#DEBUG") > 0) then
           call read_in_int(iDebugLevel)
           call read_in_int(iDebugProc)
           if (iDebugProc >= 0 .and. iProc /= iDebugProc) then
              iDebugLevel = -1
           endif
        endif

        if (index(StringLine,"#FIXTILT") > 0) then
           call read_in_logical(IsFixedTilt)
        endif

        if (index(StringLine,"#END") > 0) then
           IsDone = .true.
        endif

        if (iLine >= MaxInputLines) then
           IsDone = .true.
        endif

     else

        iLine = iLine + 1

     endif

  enddo

contains

  subroutine read_in_int(variable)
    integer, intent(out) :: variable
    iline = iline + 1
    read(StringInputLines(iline),*) variable
  end subroutine read_in_int

  subroutine read_in_logical(variable)
    logical, intent(out) :: variable
    iline = iline + 1
    read(StringInputLines(iline),*) variable
  end subroutine read_in_logical

  subroutine read_in_string(variable)
    character (len=100), intent(out) :: variable
    iline = iline + 1
    variable = StringInputLines(iline)
  end subroutine read_in_string

  subroutine read_in_real(variable)
    real :: variable
    iline = iline + 1
    read(StringInputLines(iline),*) variable
  end subroutine read_in_real

  subroutine read_in_time(variable)
    use ModTimeConvert, ONLY: time_int_to_real
    real*8 :: variable
    integer, dimension(7) :: itime
    iline = iline + 1
    read(StringInputLines(iline),*) itime(1:6)
    itime(7) = 0
    call time_int_to_real(itime, variable)
  end subroutine read_in_time

end subroutine EIE_set_inputs
