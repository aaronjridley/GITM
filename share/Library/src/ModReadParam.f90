!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!BOP
!
!QUOTE: \chapter{Fortran Libraries in share/ and util/}
!QUOTE: \section{share/Library}
!
!MODULE: ModReadParam - read, include, broadcast and distribute parameters
!INTERFACE:
module ModReadParam

  !DESCRIPTION:
  ! This is a library for reading parameters and distribute them between
  ! the components. It can also be used by a stand-alone mode.
  ! In the latter case the 'control component' corresponds to the main program
  ! of the stand-alone code.
  !
  ! {\bf Subroutine read\_file(NameFile)} reads the text from file NameFile
  ! which may include further files. The files can be included with the
  ! \begin{verbatim}
  ! #INCLUDE
  ! some/filename
  ! \end{verbatim}
  ! command. Include files can be nested up to MaxNestedFile=10 levels.
  !
  ! The text either ends at the end of the file, or at the 
  ! \begin{verbatim}
  ! #END
  ! \end{verbatim}
  ! command. After reading, the text is broadcast to all processors, or
  ! to the processors that belong to the MPI communicator iComm, which
  ! is an optional argument of subroutine read\_file.
  ! The text buffer contains at most MaxLine=1000 lines, which are at most
  ! lStringLine=200 character long. Normally only the control component
  ! calls {\bf read\_file}.
  !
  ! {\bf Subroutine read\_init} can select a part of the text buffer
  ! starting from line iLine+1 to the last line nLine.
  ! It also sets the session number, the name of the component which should
  ! read the selected text, and the output unit used by the compnent
  ! for echoing the parameters. Normally only the control component 
  ! calls {\bf read\_init}.
  !
  ! {\bf Function read\_line} reads the next line from the selected
  ! part of the text. If there are no more lines in the selected part,
  ! the logical function returns .false. The line itself and the line number
  ! are provided in optional output arguments.
  !
  ! \bigskip
  ! For the control component as well as for many physical components 
  ! the parameters are given in form of command lines followed by parameter 
  ! lines. The commands start with a \# character, which is usually 
  ! followed by upper case letters and numbers. Anything after a space
  ! or TAB (char(9)) is ignored. The number of parameter lines is determined 
  ! by the command.  Each parameter line contains either a logical, an
  ! integer, a real or a string variable. Comments placed after at least 
  ! 3 spaces or one TAB character are ignored. The lines which do not contain 
  ! a command or the corresponding parameters are ignored, and can be used 
  ! for remarks.  For example
  ! \begin{verbatim}
  ! #COMMAND1
  ! param1      Description1
  ! param2      Description2
  !
  ! some remark
  !
  ! #COMMAND2 some comment here
  ! param3      Description3
  ! \end{verbatim}
  ! {\bf Function read\_command} returns .true. if a command is found
  ! in the previously read line and it provides the name of the command
  ! in its output argument by stripping off anything behind a blank.
  !
  ! {\bf Subroutine read\_var} reads the parameter from the parameter line.
  ! The name of the parameter is provided as an input argument, the value
  ! is obtained from the output argument.
  ! 
  ! {\bf Subroutine read\_echo\_set} can be used to set the logical,
  ! which determines if the input parameters should be echoed back. 
  ! The default is no echo.
  !
  ! \bigskip
  ! The typical way to read parameters in form of commands and parameters
  ! in a component is the following
  ! \begin{verbatim}
  !   use ModReadParam
  !   implicit none
  !   character (len=lStringLine) :: NameCommand
  !   do
  !       if(.not.read_line() ) EXIT
  !       if(.not.read_command(NameCommand)) CYCLE
  !       select case(NameCommand)
  !       case("#IONODIR")
  !           call read_var("NameIonoDir", NameIonoDir)
  !       case("#IONOSPHERE")
  !           call read_var('iConductanceModel', iConductanceModel)
  !           call read_var('UseFullCurrent', UseFullCurrent)
  !           call read_var('F10.7 Flux', Flux107)
  !       ...
  !       end select
  !   end do
  ! \end{verbatim}
  ! If the component does not use the \#COMMAND parameters format,
  ! the lines can be read line by line like this
  ! \begin{verbatim}
  !   use ModReadParam
  !   implicit none
  !   character (len=lStringLine) :: StringLine
  !   do
  !       if(.not.read_line(StringLine) ) EXIT
  !       ! process StringLine
  !       read(StringLine,*) MyVariables
  !       ...
  !   end do
  ! \end{verbatim}
  ! If the component cannot process the parameters line by line,
  ! the following methods are available to read the text into a string array:
  !
  ! {\bf Function i\_line\_read()} returns the number of the line
  ! before the first line. 
  !
  ! {\bf Function n\_line\_read()} returns the line number of the last line.
  !
  ! {\bf Subroutine  read\_text} returns the selected part of the 
  ! text buffer in its output argument. 
  ! 
  ! \bigskip
  ! These methods can be used like this
  ! \begin{verbatim}
  !     use ModReadParam
  !     implicit none
  !     character(len=lStringLine), allocatable :: StringLine_I(:)
  !     allocate(StringLine_I(i_line_read()+1:n_line_read()))
  !     call read_text(StringLine_I)
  !     ! process text buffer
  !     ...
  ! \end{verbatim}

  !USES:
  use ModMpi
  use ModIoUnit, ONLY: io_unit_new, STDOUT_

  implicit none

  save

  private ! except

  !PUBLIC DATA MEMBERS:
  integer, parameter, public :: lStringLine=200 ! Max length of input lines

  !PUBLIC MEMBER FUNCTIONS:
  public :: read_file         ! Read text string from parameter file and bcast
  public :: read_init         ! Select the appropriate section of the text
  public :: read_line         ! Read next line, return false at the end
  public :: read_command      ! Read command, return false if not a command
  public :: read_var          ! Read scalar variable of any type with a name
  public :: read_in           ! Read scalar variable of any type without name
  public :: read_echo_set     ! Set if parameters should be echoed
  public :: i_line_read       ! Return the current line number
  public :: n_line_read       ! Return the last line number in the selected
  public :: i_session_read    ! Return the session number
  public :: i_line_command    ! Returns the line number for a command or -1
  public :: read_text         ! Provide the full text in the output argument

  !REVISION HISTORY:
  ! 01Sep03 G. Toth - initial implementation based on BATSRUS
  ! 31Oct04 G. Toth - added fractions 3/5 for reals, 
  !                   added read_in public method without a Name parameter,
  !                   replaced err=1 with iostat=iReadError.
  ! 27Nov06 G. Toth - added i_line_command function
  !EOP

  character(len=*), parameter :: NameMod='ModReadParam'

  ! Text buffer to hold all the input lines
  integer, parameter :: MaxLine=1000
  character (len=lStringLine) :: StringLine_I(MaxLine)

  character(len=lStringLine)  :: StringLine, StringParam

  character(len=2)  :: NameComp    = '  '  ! Name of the component
  character(len=3)  :: StringPrefix= '   ' ! Prefix for echo
  integer           :: iIoUnit     = -1    ! Unit number for echo
  integer           :: iLine=0             ! Current line number
  integer           :: nLine=0             ! Last line number
  integer           :: iSession=0          ! Session number
  logical           :: DoEcho = .false.    ! Do we echo parameters?

  ! Storage for all the commands for all sessions:
  ! the command index is increased from 1 to iCommand in every session
  ! and the name and line number is stored. The name contains the 
  ! component and command name together.
  integer, parameter :: MaxCommand = MaxLine/3, MaxSession = 10
  integer            :: iCommand = 0
  integer            :: iLineCommand_II(MaxCommand, MaxSession) = -1
  character(len=lStringLine) :: NameCommand_II(MaxCommand, MaxSession)

  interface read_var
     module procedure &
          read_var_c, read_var_l, read_var_r, read_var_i
  end interface

  interface read_in
     module procedure &
          read_integer, read_real, read_string, read_logical
  end interface

contains

  !BOP =======================================================================
  !IROUTINE: read_file - read parameter file
  !INTERFACE:
  subroutine read_file(NameFile, iCommIn, NameRestartFile)

    !INPUT ARGUMENTS:
    character (len=*), intent(in):: NameFile ! Name of the base param file
    integer, optional, intent(in):: iCommIn  ! MPI communicator for broadcast

    ! Name of the restart file to be read if a #RESTART command is found
    character (len=*), intent(in), optional :: NameRestartFile 

    !EOP
    integer, parameter :: MaxNestedFile = 10

    character(len=*), parameter :: NameSub=NameMod//'::read_file'

    character (len=lStringLine) :: NameCommand

    integer :: iUnit_I(MaxNestedFile)

    integer :: iFile, i, iError, iProc, iComm

    logical :: IsFound

    logical :: Done=.false., DoInclude
    !-----------------------------------------------------------------------
    if(Done)call CON_stop(NameSub//&
         ' ERROR: the parameter file should be read only once!')

    if(present(iCommIn))then
       iComm = iCommIn
    else
       iComm = MPI_COMM_WORLD
    end if

    ! Get processor rank
    call MPI_comm_rank(iComm,iProc,iError)

    !\
    ! Read all input file(s) into memory and broadcast
    !/
    if(iProc==0)then
       nLine=0
       inquire(file=NameFile,EXIST=IsFound)
       if(.not.IsFound)call CON_stop(NameSub//' SWMF_ERROR: '//&
            trim(NameFile)//" cannot be found")
       iFile=1
       iUnit_I(iFile)=io_unit_new()
       open(iUnit_I(iFile),file=NameFile,status="old")
       do
          read(iUnit_I(iFile),'(a)',ERR=100,END=100) StringLine
          NameCommand=StringLine
          i=index(NameCommand,' '); 
          if(i>0)NameCommand(i:len(NameCommand))=' '
          i=index(NameCommand,char(9)); 
          if(i>0)NameCommand(i:len(NameCommand))=' '
          if(NameCommand=='#INCLUDE')then
             ! Include text from file following the #INCLUDE command
             read(iUnit_I(iFile),'(a)')StringLine
             ! Remove anything after a space or TAB
             i=index(StringLine,' ')
             if(i>0)StringLine(i:len(StringLine))=' '
             i=index(StringLine,char(9))
             if(i>0)StringLine(i:len(StringLine))=' '
             DoInclude = .true.
          elseif(present(NameRestartFile) .and. NameCommand=='#RESTART')then
             ! If #RESTART command is followed by true then include
             ! the file named NameRestartFile
             read(iUnit_I(iFile),*,IOSTAT=iError)DoInclude
             if(iError>0)then
                write(*,*) NameSub,&
                     " ERROR: could not read logical after #RESTART command",&
                     " at line ",nLine+1
                call CON_stop("Correct "//trim(NameFile))
             end if
             if(DoInclude)then
                StringLine = NameRestartFile
             else
                StringLine = ' ' ! remove #RESTART command
             end if
          else
             DoInclude = .false.
          end if
          if(DoInclude)then
             iFile = iFile + 1
             if(iFile > MaxNestedFile)call CON_stop(NameSub// &
                  " SWMF_ERROR: more than MaxNestedFile nested files")
             inquire(file=StringLine,EXIST=IsFound)
             if(.not.IsFound)call CON_stop(NameSub// &
                  " SWMF_ERROR: include file cannot be found, name="//&
                  trim(StringLine))
             iUnit_I(iFile) = io_unit_new()
             open(iUnit_I(iFile),FILE=StringLine,STATUS="old")
             CYCLE
          else if(NameCommand/='#END')then
             ! Store line into buffer
             nLine=nLine+1
             if(nLine>maxline)call CON_stop(NameSub// &
                  " SWMF_ERROR: too many lines of input")
             StringLine_I(nLine)=StringLine
             CYCLE
          end if

100       continue
          close (iUnit_I(iFile))
          if(iFile > 1)then
             ! Continue reading the calling file
             iFile = iFile - 1
             CYCLE
          else
             ! The base file ended, stop reading
             EXIT
          end if
       end do
       if(nLine==0)call CON_stop(NameSub// &
            " SWMF_ERROR: no lines of input read")
    end if
    ! Broadcast the number of lines and the text itself to all processors
    call MPI_Bcast(nLine,1,MPI_INTEGER,0,iComm,iError)

    if(iError>0)call CON_stop(NameSub// &
         " MPI_ERROR: number of lines could not be broadcast")

    call MPI_Bcast(StringLine_I,lStringLine*nLine,MPI_CHARACTER,&
         0,iComm,iError)

    if(iError>0)call CON_stop(NameSub// &
         " MPI_ERROR: text could not be broadcast")

    if(iProc==0)write(*,'(a,i4,a)') NameSub// &
         ': read and broadcast nLine=',nLine,' lines of text'

    Done = .true.

  end subroutine read_file

  !===========================================================================

  subroutine read_init(NameCompIn, iSessionIn, iLineIn, nLineIn, iIoUnitIn)

    ! Initialize module variables

    character(len=2), optional, intent(in) :: NameCompIn
    integer,          optional, intent(in) :: iSessionIn
    integer,          optional, intent(in) :: iLineIn, nLineIn, iIoUnitIn

    integer:: iSessionNew

    character(len=*), parameter :: NameSub=NameMod//'::read_init'
    !------------------------------------------------------------------------
    if(present(iSessionIn))then
       iSessionNew = iSessionIn
    else
       iSessionNew = 1
    end if

    if(iSessionNew > MaxSession)call CON_stop(NameSub// &
         " ERROR: too many sessions in input")

    ! Set command counter to zero for a new session
    if(iSessionNew > iSession) iCommand = 0

    iSession     = iSessionNew

    if(present(NameCompIn))then
       NameComp = NameCompIn
    else
       NameComp = ''
    end if
    if(present(iLineIn)) iLine = iLineIn
    if(present(nLineIn))then
       nLine     = nLineIn
    else
       nLine     = size(StringLine_I)
    end if
    if(present(iIoUnitIn))then
       iIoUnit   = iIoUnitIn
    else
       iIoUnit   = STDOUT_
    end if
    if(iIoUnit==STDOUT_ .and. len_trim(NameComp)>0 )then
       StringPrefix = NameComp//': '
    else
       StringPrefix = ''
    end if
    
  end subroutine read_init

  !===========================================================================

  subroutine read_echo_set(DoEchoIn)

    logical, intent(in) :: DoEchoIn

    DoEcho = DoEchoIn

  end subroutine read_echo_set

  !BOP =======================================================================
  !IROUTINE: read_line - read the next line from the text buffer
  !INTERFACE:
  logical function read_line(StringLineOut, iLineOut)

    !OUTPUT ARGUMENTS:
    character (len=*), optional, intent(out) :: StringLineOut
    integer, optional, intent(out)           :: iLineOut

    !DESCRIPTION:
    ! Read the current line from StringLine\_I into StringLine,
    ! set the optional StringLineOut and iLineOut arguments.
    ! Return .true. if successful, otherwise (if there are
    ! no more lines in the selected part of the text buffer)
    ! return .false. and an empty string in StringLineOut if present.
    !EOP

    iLine=iLine+1
    if(present(iLineOut)) iLineOut = iLine
    if(iLine <= nLine)then
       StringLine = StringLine_I(iLine)
       if(present(StringLineOut)) StringLineOut = StringLine
       read_line  = .true.
    else
       if(present(StringLineOut)) StringLineOut = ''
       read_line = .false.
    endif

  end function read_line

  !BOP =======================================================================
  !IROUTINE: read_command - read the name of the command from the current line
  !INTERFACE:
  logical function read_command(NameCommand)

    !OUTPUT ARGUMENTS:
    character (len=*), intent(out) :: NameCommand

    !DESCRIPTION:
    ! If the current line contains a command name (starting with \#),
    ! return true, and put the name of the command into the 
    ! output argument. Otherwise return .false. and an empty string.
    !EOP

    integer :: i

    !------------------------------------------------------------------------
    if(StringLine(1:1)=="#")then

       if(DoEcho)then
          write(iIoUnit,'(a)')trim(StringPrefix)
          write(iIoUnit,'(a)')trim(StringPrefix)//trim(StringLine)
       end if

       ! Remove anything after a space or TAB
       i=index(StringLine,' ');     if(i>0)StringLine(i:len(StringLine))=' '
       i=index(StringLine,char(9)); if(i>0)StringLine(i:len(StringLine))=' '

       NameCommand = StringLine
       read_command = .true.

       ! Store command
       iCommand = iCommand + 1
       iLineCommand_II(iCommand, iSession) = iLine
       NameCommand_II(iCommand, iSession)  = NameComp//NameCommand

    else
       NameCommand  = ''
       read_command = .false.
    endif

  end function read_command

  !===========================================================================

  subroutine read_line_param(Type, Name, iError, DoReadWholeLine)

    ! read next line from text

    character (len=*), intent(in) :: Type, Name
    integer, optional, intent(out):: iError
    logical, optional, intent(in) :: DoReadWholeLine

    integer :: i, j
    !------------------------------------------------------------------------

    if(present(iError))iError=0
    iLine=iLine+1
    if(iLine>nLine)then
       if(present(iError))then
          iError=1
       else
          call CON_stop(&
               'Unexpected end of text after line='//StringLine)
       end if
    end if
    StringLine=StringLine_I(iLine)

    ! Get rid of leading spaces 
    StringParam = adjustl(StringLine)

    if(.not.present(DoReadWholeLine))then
       ! Get rid of trailing comments after a TAB character or 3 spaces
       i=index(StringParam,char(9))
       j=index(StringParam,'   ')
       if(i>1.and.j>1)then
          i = min(i,j)      ! Both TAB and 3 spaces found, take the closer one
       else
          i = max(i,j)      ! Take the one that was found (if any)
       end if
       if(i>0) StringParam(i:lStringLine)=' '
    end if

    if(len_trim(StringParam)==0)call read_error('missing',Name,iError)

  end subroutine read_line_param

  !===========================================================================

  subroutine read_echo(Name)

    character (len=*), intent(in)    :: Name
    !-------------------------------------------------------------------------

    if(index(StringLine,Name)<1) &
         StringLine=trim(StringLine)//char(9)//char(9)//Name
    write(iIoUnit,'(a)') trim(StringPrefix)//trim(StringLine)

  end subroutine read_echo

  !===========================================================================
  subroutine read_error(Type, Name, iError)

    ! Print error message for reading error of variable named Name of type Type

    character (len=*), intent(in) :: Type, Name
    integer, optional, intent(out):: iError

    !------------------------------------------------------------------------

    if(present(iError))then
       select case(Type)
       case('missing')
          iError = 2
       case('integer')
          iError = 3
       case('logical')
          iError = 4
       case('real')
          iError = 5
       case('character')
          iError = 6
       case default
          iError = -1
       end select
    else
       write(*,'(a,i3)')'Error in component '//NameComp//' in session',iSession
       call CON_stop('Error reading '//Type//' variable '//Name// &
            ' from line='//StringLine)
    end if

  end subroutine read_error

  !===========================================================================

  subroutine read_string(StringVar, iError, IsUpperCase, IsLowerCase, &
       DoReadWholeLine)

    ! Read a string variable
    ! Convert to upper or lower case if required.
    ! Arguments
    character (len=*), intent(out):: StringVar
    integer, optional, intent(out):: iError
    logical, optional, intent(in) :: IsUpperCase, IsLowerCase, DoReadWholeLine
    !-------------------------------------------------------------------------
    call read_var_c(' ', StringVar, iError, IsLowerCase, IsUpperCase, &
         DoReadWholeLine)

  end subroutine read_string

  !===========================================================================

  subroutine read_var_c(Name, StringVar, iError, &
       IsUpperCase, IsLowerCase, DoReadWholeLine)

    use ModUtilities, ONLY: lower_case, upper_case

    ! Read a string variable described by Name.
    ! Convert to upper or lower case if required.
    ! Arguments
    character (len=*), intent(in) :: Name
    character (len=*), intent(out):: StringVar
    integer, optional, intent(out):: iError
    logical, optional, intent(in) :: IsUpperCase, IsLowerCase, DoReadWholeLine
    !-------------------------------------------------------------------------

    call read_line_param('character', Name, iError, DoReadWholeLine)

    if(DoEcho)call read_echo(Name)

    StringVar=StringParam

    if(present(IsLowerCase))then
       if(IsLowerCase)call lower_case(StringVar)
    endif
    if(present(IsUpperCase))then
       if(IsUpperCase)call upper_case(StringVar)
    endif

  end subroutine read_var_c

  !===========================================================================
 
  subroutine read_integer(IntVar, iError)
    !OUTPUT ARGUMENTS:
    integer,           intent(out):: IntVar
    integer, optional, intent(out):: iError
    !-------------------------------------------------------------------------
    call read_var_i(' ', IntVar, iError)
  end subroutine read_integer

  !BOP ========================================================================
  !IROUTINE: read_var - read a variable following the command.
  !INTERFACE:
  subroutine read_var_i(Name, IntVar, iError)

    !INPUT ARGUMENTS:
    character (len=*), intent(in) :: Name
    !OUTPUT ARGUMENTS:
    integer,           intent(out):: IntVar
    integer, optional, intent(out):: iError

    !DESCRIPTION:
    ! Read a variable from the next line in the buffer.
    ! The variable name is given by the string Name, which is used in the
    ! echoing of the parameters as well as in error messages. 
    ! The value is returned in the non-optional output argument. 
    ! If the optional argument iError is present, it returns a non-zero
    ! value in case an error occurs. If iError is not present, all errors
    ! result in an error message and an abort of the run.
    ! There are four variants of this subroutine: for integer, real,
    ! character string and logical variable types.
    !EOP

    ! Local variable
    integer :: IntTmp, iReadError

    !-------------------------------------------------------------------------
    call read_line_param('integer', Name, iError)

    read(StringParam,*,iostat=iReadError) IntTmp
    if(iReadError/=0) call read_error('integer', Name, iError)
    if(DoEcho)        call read_echo(Name)
    IntVar=IntTmp

  end subroutine read_var_i

  !===========================================================================

  subroutine read_real(RealVar,iError)
    ! Read a real variable
    ! Arguments
    real, intent(out)             :: RealVar
    integer, optional, intent(out):: iError
    !-------------------------------------------------------------------------
    call read_var_r(' ', RealVar, iError)
  end subroutine read_real

  !===========================================================================

  subroutine read_var_r(Name, RealVar, iError)

    ! Read a real variable described by Name
    
    ! Arguments
    character (len=*), intent(in)   :: Name
    real, intent(out)               :: RealVar
    integer, optional, intent(out)  :: iError

    ! Local variable
    real :: Numerator, Denominator, RealTmp
    integer :: i, iReadError
    !-------------------------------------------------------------------------
    call read_line_param('real', Name, iError)

    ! Check for fraction
    i = index(StringParam,'/')
    if(i<1)then
       ! Simple real number
       read(StringParam,*,iostat=iReadError) RealTmp
       if(iReadError/=0)      call read_error('real',Name,iError)
    else
       ! Fraction: Numerator/Denominator
       read(StringParam(1:i-1),*,iostat=iReadError) Numerator
       if(iReadError/=0)      call read_error('numerator',Name,iError)
       read(StringParam(i+1:lStringLine),*,iostat=iReadError) Denominator
       if(iReadError/=0)      call read_error('denominator',Name,iError)
       if(Denominator == 0.0) call read_error('zero denominator',Name,iError)
       RealTmp = Numerator / Denominator
    end if
    if(DoEcho)call read_echo(Name)
    RealVar=RealTmp

  end subroutine read_var_r

  !===========================================================================

  subroutine read_logical(IsLogicVar, iError)
    ! Read a logical variable
    ! Arguments
    logical, intent(out)          :: IsLogicVar
    integer, optional, intent(out):: iError
    !-------------------------------------------------------------------------
    call read_var_l(' ',IsLogicVar, iError)
  end subroutine read_logical

  !===========================================================================

  subroutine read_var_l(Name, IsLogicVar, iError)

    ! Read a logical variable described by Name

    ! Arguments
    character (len=*), intent(in) :: Name
    logical, intent(out)          :: IsLogicVar
    integer, optional, intent(out):: iError

    ! Local variable
    logical :: IsLogicTmp
    integer :: iReadError

    !-------------------------------------------------------------------------
    call read_line_param('logical', Name, iError)

    read(StringParam,*,iostat=iReadError) IsLogicTmp
    if(iReadError/=0)call read_error('logical',Name,iError)
    if(DoEcho)call        read_echo(Name)
    IsLogicVar=IsLogicTmp

  end subroutine read_var_l

  !===========================================================================

  integer function i_line_read()
    i_line_read = iLine
  end function i_line_read

  !===========================================================================

  integer function n_line_read()
    n_line_read = nLine
  end function n_line_read

  !===========================================================================

  integer function i_session_read()
    i_session_read = iSession
  end function i_session_read

  !===========================================================================

  integer function i_line_command(NameCommandIn, iSessionIn)

    ! Search the command name for the current component.
    ! If iSessionIn is present, search for session iSessionIn, 
    ! otherwise search the currenst session. If the command is found
    ! return the line number in the input file, otherwise return -1.

    character(len=*), intent(in) :: NameCommandIn
    integer, optional, intent(in):: iSessionIn

    integer :: i, j
    character(len=lStringLine) :: NameCommand
    !------------------------------------------------------------------------
    NameCommand = NameComp // trim(NameCommandIn)

    j = iSession
    if(present(iSessionIn)) j = iSessionIn

    do i = 1, MaxCommand
       if(iLineCommand_II(i, j) < 0) EXIT
       if(NameCommand == NameCommand_II(i, j)) then
          i_line_command = iLineCommand_II(i, j)
          RETURN
       end if
    end do

    ! Could not find command, return -1
    i_line_command = -1
    
  end function i_line_command

  !BOP =======================================================================
  !IROUTINE: read_text - obtain selected text buffer 
  !INTERFACE:
  subroutine read_text(String_I)
    !OUTPUT ARGUMENTS:
    character(len=lStringLine), intent(out):: String_I(iLine+1:nLine)
    !EOP
    !BOC
    String_I = StringLine_I(iLine+1:nLine)
    !EOC
  end subroutine read_text

  !===========================================================================

end module ModReadParam
