!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
!BOP
!
!MODULE: ModUtilities - Simple Methods for CON and Components
!INTERFACE:
module ModUtilities

  !DESCRIPTION:
  ! Simple methods which are used by CON and can be used 
  ! by the science components too. 
  !
  ! This module is almost self contained. Only ModIoUnit, ModMpi, ModKind 
  ! and the external subroutine CON\_stop are used.
  !
  ! F77 and C++ codes need an F90 interface to access these utilities.
  !EOP

  implicit none

  private ! except

  public:: make_dir
  public:: check_dir
  public:: fix_dir_name
  public:: flush_unit
  public:: split_string
  public:: join_string
  public:: upper_case
  public:: lower_case
  public:: sleep
  public:: check_allocate
  public:: test_mod_utility

  logical, public :: DoFlush = .true.

  interface split_string
     module procedure split_string, split_string_simple
  end interface

  interface join_string
     module procedure join_string, join_string_simple
  end interface

contains


  !BOP ========================================================================
  !ROUTINE: make_dir - Create a directory
  !INTERFACE:
  subroutine make_dir(NameDir, iPermissionIn, iErrorOut, iErrorNumberOut)

    use iso_c_binding

    !INPUT ARGUMENTS:
    character(len=*), intent(in):: NameDir ! Directory name

    integer, intent(in), optional:: iPermissionIn ! Access permission

    !OUTPUT ARGUMENTS:
    integer, intent(out), optional:: iErrorOut    ! 0 on success, -1 otherwise
    integer, intent(out), optional:: iErrorNumberOut ! C error number

    !DESCRIPTION:
    ! Create the directory specified by NameDir. Trailing spaces are ignored.
    ! Nested directories are also allowed.
    ! If the directory already exists, this function does nothing.
    !
    ! The optional iPermissionIn parameter sets the permissions for the new
    ! directory. This should be specified in octal notation, which in
    ! Fortran is written with a capital O followed by the digits in quotes. 
    ! For instance, the permissions 0755 would be written as O'0775',
    ! which is the default value (drwxr-xr-x).
    ! The actual directory will have permissions modified by 
    ! the default mask (set by the Unix command umask). 
    !
    ! The optional iErrorOut is set to 0 if no error occured and -1 otherwise. 
    ! If iError is not present, the code stops with an error message.
    !
    ! The iErrorNumber contains the return value from the C mkdir() function.
    ! Values of iErrorNumber are found in errno.h and are not standardized,
    ! so exact values of these should not be relied upon.
    !EOP

    integer(c_int):: iPermission ! Octal permissions (value passed to C)
    integer:: iErrorNumber       ! Error number as returned from C
    integer:: iError             ! Return value as retrieved from C

    interface
       ! Interface for mkdir_wrapper implemented in ModUtility_c.c
       integer(kind=c_int) function make_dir_c(path, perm, errno) bind(C)
         use iso_c_binding
         character(kind=c_char), intent(in):: path
         integer(kind=c_int), intent(in), value:: perm
         integer(kind=c_int), intent(out):: errno
       end function make_dir_c

    end interface

    integer:: i, j, k, l

    character(len=*), parameter:: NameSub = 'make_dir'
    !------------------------------------------------------------------------
    if(.not. present(iPermissionIn)) then
       iPermission = O'0755'
    else
       iPermission = iPermissionIn
    endif

    ! Create the directory. 
    ! If NameDir contains one or more /, then create the top directory first
    ! and then the subdirectories.

    l = len_trim(NameDir)
    i = 1
    do
       ! Find the next '/' in NameDir
       j = l
       k = index(NameDir(i:l), '/')
       if(k > 1) j = i + k - 1

       ! Create (sub)directory
       iError = make_dir_c(NameDir(1:j)//C_NULL_CHAR, iPermission, iErrorNumber)

       ! Check for errors
       if(iError /= 0)then
          if(present(iErrorOut))then
             iErrorOut = iError
             if(present(iErrorNumberOut)) iErrorNumberOut = iErrorNumber
             RETURN
          else
             write(*,*) NameSub,' iError, iErrorNumber=', iError, iErrorNumber 
             call CON_stop(NameSub//' failed to open directory '//trim(NameDir))
          end if
       end if

       i = j + 1
       if(i > l) EXIT

    end do

    if(present(iErrorOut))       iErrorOut = iError
    if(present(iErrorNumberOut)) iErrorNumberOut = iErrorNumber

  end subroutine make_dir

  !BOP ========================================================================
  !ROUTINE: check_dir - check if a directory exists
  !INTERFACE:
  subroutine check_dir(NameDir)

    !USES:
    use ModIoUnit, ONLY: UNITTMP_

    !INPUT ARGUMENTS:
    character(len=*), intent(in) :: NameDir

    !DESCRIPTION:
    ! Check if a directory exists by trying to open a file in it.
    ! Die with an error message if the directory does not exist.
    ! Directory names are cached, so multiple calls with the same
    ! name do not repeat the check.
    !
    ! {\bf This subroutine should be called by the root PE of the component 
    ! only!}
    ! Calling the subroutine from multiple PE-s may result in a fatal error,
    ! namely one PE may delete the file written by the other PE, so the
    ! other PE thinks that the directory does not exist.
    !EOP

    integer, parameter :: MaxDir=100, lNameDir=100
    integer, save :: nDir=0

    character(len=lNameDir), save :: NameDir_I(MaxDir)
    integer :: iDir, iError

    character(len=*), parameter :: NameSub='check_dir'
    !--------------------------------------------------------------------------
    ! Only directory names shorter than lNameDir can be stored
    if(len_trim(NameDir) <= lNameDir) then
       ! Check if this directory has been checked already
       do iDir=1,nDir
          if(NameDir_I(iDir)==NameDir) RETURN
       end do

       ! Increase counter for different directory names
       nDir=nDir+1

       ! Store new name if possible. 
       if(nDir <= MaxDir) NameDir_I(nDir)=NameDir

       ! If not, warn once, and keep checking...
       if(nDir == MaxDir+1) write(*,'(a)')NameSub // &
            ' SWMF_WARNING: too many different directories!'
    end if

    ! Try to open a file in this directory
    open(UNITTMP_, file=trim(NameDir)//'.test', status='unknown', &
         iostat = iError)

    if (iError /= 0) then
       call CON_stop(NameSub//' ERROR: Cannot find/write into directory '&
            //trim(NameDir))
    else
       close(UNITTMP_, status = 'DELETE')
    endif

  end subroutine check_dir
  
  !BOP ========================================================================
  !ROUTINE: fix_dir_name - add a slash to the end of the directory name
  !INTERFACE:
  subroutine fix_dir_name(NameDir)

    !INPUT/OUTPUT ARGUMENTS:
    character(len=*), intent(inout) :: NameDir

    !DESCRIPTION:
    ! Append a '/' at the end of the directory name if it is not there
    ! and the directory name is not zero length (empty string).
    !
    ! {\bf This subroutine should be called by all PE-s of the component!}
    !EOP

    character(len=*), parameter :: NameSub='fix_dir_name'
    integer :: i
    !--------------------------------------------------------------------------
    i = len_trim(NameDir)
    if(i == 0) RETURN
    if(NameDir(i:i) == '/') RETURN

    if(i >= len(NameDir)) call CON_stop(NameSub// &
         "ERROR cannot append / to directory name "//NameDir)

    NameDir(i+1:i+1) = '/'

  end subroutine fix_dir_name

  !BOP ========================================================================
  !ROUTINE: flush_unit - flush output
  !INTERFACE:
  subroutine flush_unit(iUnit)

    !INPUT ARGUMENTS:
    integer, intent(in) :: iUnit

    !DESCRIPTION:
    ! Flush the I/O unit iUnit if DoFlush is true
    !EOP

    !-------------------------------------------------------------------------
    if(DoFlush) flush(iUnit)

  end subroutine flush_unit

  !BOP ========================================================================
  !ROUTINE: split_string_simple - split string into array of substrings
  !INTERFACE:
  subroutine split_string_simple(String, String_I, nStringOut, &
       StringSepIn, UseArraySyntaxIn, DoAddSeparator)

    !INPUT ARGUMENTS:
    character(len=*),    intent(in):: String    ! string to be split

    !OUTPUT ARGUMENTS:
    character (len=*), intent(out):: String_I(:) ! array of substrings

    !OPTIONAL ARGUMENTS
    integer,          optional, intent(out):: nStringOut ! number of substrings
    character(len=*), optional, intent(in):: StringSepIn ! separator string
    logical, optional, intent(in):: UseArraySyntaxIn     ! expand Var(10:20:2)
    logical, optional, intent(in):: DoAddSeparator !add separator to substrings

    call split_string(String, size(String_I), String_I, nStringOut, &
         StringSepIn, UseArraySyntaxIn, DoAddSeparator)

  end subroutine split_string_simple

  !BOP ========================================================================
  !ROUTINE: split_string - split string into array of substrings
  !INTERFACE:
  subroutine split_string(String, MaxString, String_I, nStringOut, &
       StringSepIn, UseArraySyntaxIn, DoAddSeparator)

    !INPUT ARGUMENTS:
    character(len=*),    intent(in):: String    ! string to be split
    integer,             intent(in):: MaxString ! maximum array size

    !OUTPUT ARGUMENTS:
    character (len=*), intent(out):: String_I(MaxString) ! array of substrings

    !OPTIONAL ARGUMENTS
    integer,          optional, intent(out):: nStringOut ! number of substrings
    character(len=*), optional, intent(in):: StringSepIn ! separator string
    logical, optional, intent(in):: UseArraySyntaxIn     ! expand Var(10:20:2)
    logical, optional, intent(in):: DoAddSeparator !add separator to substrings

    !DESCRIPTION:
    ! Cut the input string into an array of substrings. The separator
    ! string is either StringSepIn or a single space (default). 
    ! Multiple consecutive separator strings are treated as one.
    ! Leading and trailing spaces are ignored. For example
    !\begin{verbatim}
    ! ' IE  GM ' --> nString=2, String\_I=(/'IE','GM'/)
    !\end{verbatim}
    ! When UseArraySyntax is present, then expand strings containing
    ! parens into an array of substrings ending with numbers, e.g.
    !\begin{verbatim}
    ! 'Var(4)'      --> nString=4,  String\_I=(/'Var1','Var2','Var3','Var4'/)
    ! 'Var(11)'     --> nString=11, String\_I=(/'Var01','Var02',...,'Var11'/)
    ! 'Var(3:5)'    --> nString=3,  String\_I=(/'Var3','Var4','Var5'/)
    ! 'Var(7:11:2)' --> nString=3,  String\_I=(/'Var07','Var09','Var11'/)
    !\end{verbatim}
    !EOP

    integer          :: nString
    character(len=10):: StringSep
    logical          :: UseArraySyntax

    character(len=len(String)+10) :: StringTmp

    integer :: i, l, lSep, lKeep

    character(len=*), parameter :: NameSub = 'split_string'
    !--------------------------------------------------------------------------
    if(present(StringSepIn))then
       StringSep = StringSepIn
       lSep      = len(StringSepIn)
    else
       StringSep = ' '
       lSep      = 1
    end if

    UseArraySyntax = .false.
    if(present(UseArraySyntaxIn)) UseArraySyntax = UseArraySyntaxIn

    lKeep = 0
    if(present(DoAddSeparator))then
       if(DoAddSeparator) lKeep = lSep
    end if
    
    nString   = 0
    StringTmp = String
    l         = len_trim(StringTmp)
    StringTmp = trim(StringTmp) // StringSep(1:lSep) ! Add separator to the end
    do
       StringTmp = adjustl(StringTmp)          ! Remove leading spaces
       i = index(StringTmp, StringSep(1:lSep)) ! Find end of first part   
       if(i <= 1) EXIT                         ! Nothing before the separator
       nString = nString + 1                   ! Count parts

       if(lKeep>0)then                         ! Do not keep added separator
          if(i+lKeep >= len_trim(StringTmp)) lKeep = 0 
       end if

       String_I(nString) = StringTmp(1:i-1+lKeep) ! Put part into string array
       StringTmp = StringTmp(i+lSep:l+lSep) ! Delete part+separator from string

       if(UseArraySyntax) call expand_array(String_I(nString))

       if(nString == MaxString) EXIT        ! Check for maximum number of parts
    end do

    if(present(nStringOut)) nStringOut = nString

  contains
    !========================================================================
    subroutine expand_array(String1)

      ! Expand String1 if it contains array syntax, e.g.
      ! "Var(04)"     to   "Var01", "Var02", "Var03", "Var04"
      ! "Var(2:4)"    to   "Var2", "Var3", "Var4"
      ! "Var(8:12:2)" to   "Var08", "Var10", "Var12"

      character(len=*), intent(inout):: String1

      character(len=len(String1)) :: String2
      character(len=6):: StringFormat

      integer:: j, k, l, m, lNum, iFirst, iLast, Di, iNum, iError
      !---------------------------------------------------------------------
      ! Find the opening paren if any
      j = index(String1,'(')
      if(j < 1) RETURN
      k = index(String1,')')

      if(k < j) &
           call CON_stop(NameSub//' missing closing paren in String='//String)

      ! Check for colon
      l = index(String1,':')
      if(l > j)then
         ! read initial index value before the first colon        
         read(String1(j+1:l-1),*,IOSTAT=iError) iFirst
         if(iError /= 0 .or. iFirst < 1) call CON_stop(NameSub// &
              ' invalid initial index value in String='//String)
      else
         iFirst = 1
         l = j
      end if
      
      ! Check for a second colon
      m = index(String1,':',back=.true.)
      if(m > l)then
         ! read index stride value after the seecond colon        
         read(String1(m+1:k-1),*,IOSTAT=iError) Di
         if(iError /= 0 .or. Di < 1) call CON_stop(NameSub// &
              ' invalid index stride value in String='//String)
      else
         Di = 1
         m  = k
      end if

      ! read the last index value between the l and m positions
      read(String1(l+1:m-1),*,IOSTAT=iError) iLast
      if(iError /= 0 .or. iLast < iFirst) call CON_stop(NameSub// &
           ' invalid maximum index value in String='//String)

      ! Set length of numerical string and the format string
      lNum = m - l - 1
      write(StringFormat,'(a,i1,a,i1,a)') "(i",lNum,".",lNum,")"

      ! Set the beginning part of the string to the variable name
      String2 = ' '
      String2(1:j-1) = String1(1:j-1)

      ! Expand variable names by repating name and adding numerical value
      nString = nString - 1
      do iNum = iFirst, iLast, Di
 
         write(String2(j:j+lNum),StringFormat) iNum
         nString = nString + 1
         String_I(nString) = String2

         if(nString == MaxString) RETURN
         
      end do
    end subroutine expand_array

  end subroutine split_string

  !BOP ========================================================================
  !ROUTINE: join_string_simple - join string array into a single string
  !INTERFACE:

  subroutine join_string_simple(String_I, String, StringSepIn)

    !INPUT ARGUMENTS:
    character(len=*),    intent(in):: String_I(:)

    !OUTPUT ARGUMENTS:
    character(len=*), intent(out):: String

    !OPTIONAL ARGUMENTS:
    character(len=*), optional, intent(in):: StringSepIn

    call join_string(size(String_I), String_I, String, StringSepIn)

  end subroutine join_string_simple

  !BOP ========================================================================
  !ROUTINE: join_string - join string array into a single string
  !INTERFACE:

  subroutine join_string(nString, String_I, String, StringSepIn)

    !INPUT ARGUMENTS:
    integer,             intent(in):: nString
    character(len=*),    intent(in):: String_I(nString)

    !OUTPUT ARGUMENTS:
    character (len=*), intent(out):: String

    !OPTIONAL ARGUMENTS:
    character(len=*), optional, intent(in):: StringSepIn

    !DESCRIPTION:
    ! Join the input string array into one string. The parts are joined
    ! with spaces or the optional StringSepIn string (up to 10 characters).
    !EOP

    character(len=10):: StringSep

    integer :: i, l

    character(len=*), parameter :: NameSub = 'join_string'
    !--------------------------------------------------------------------------
    if(present(StringSepIn))then
       StringSep = StringSepIn
       l = len(StringSepIn)
    else
       StringSep = ' '
       l = 1
    endif
    String = String_I(1)

    do i = 2, nString
      String = trim(String) // StringSep(1:l) // String_I(i)
    end do

  end subroutine join_string

  !BOP ========================================================================
  !ROUTINE: upper_case - convert string to all upper case
  !INTERFACE:
  subroutine upper_case(String)

    !INPUT/OUTPUT ARGUMENTS:
    character (len=*), intent(inout) :: String

    !DESCRIPTION:
    ! Change characters to upper case in String
    !EOP

    integer, parameter :: iA=ichar('a'), iZ=ichar('z'), Di=ichar('A')-iA
    integer :: i, iC
    !--------------------------------------------------------------------------
    do i = 1, len_trim(String)
       iC = ichar(String(i:i))
       if(iC >= iA .and. iC <= iZ) String(i:i) = char(iC+Di)
    end do

  end subroutine upper_case

  !BOP ========================================================================
  !ROUTINE: lower_case - convert string to all lower case
  !INTERFACE:
  subroutine lower_case(String)

    !INPUT/OUTPUT ARGUMENTS:
    character (len=*), intent(inout) :: String

    !DESCRIPTION:
    ! Change characters to lower case in String
    !EOP

    integer, parameter :: iA=ichar('A'), iZ=ichar('Z'), Di=ichar('a')-iA
    integer :: i, iC
    !--------------------------------------------------------------------------
    do i = 1, len_trim(String)
       iC = ichar(String(i:i))
       if(iC >= iA .and. iC <= iZ) String(i:i) = char(iC+Di)
    end do

  end subroutine lower_case

  !BOP ========================================================================
  !ROUTINE: sleep - sleep a given number of seconds
  !INTERFACE:
  subroutine sleep(DtCpu)
    !USES
    use ModMpi, ONLY : MPI_wtime
    use ModKind

    !INPUT ARGUMENTS:
    real, intent(in) :: DtCpu  ! CPU time to sleep (in seconds)
    !LOCAL VARIABLES:
    real(Real8_) :: tCpu0
    !DESCRIPTION:
    ! This subroutine returns after the number of seconds
    ! given in its argument.
    !EOP
    !BOC
    tCpu0 = MPI_WTIME()
    do
       if(MPI_WTIME() > tCpu0 + DtCpu) RETURN
    end do
    !EOC
  end subroutine sleep

  !BOP ========================================================================
  !ROUTINE: check_allocate - check and stop for allocation errors
  !INTERFACE:
  subroutine check_allocate(iError,NameArray)

    !INPUT ARGUMENTS:
    integer,intent(in)::iError
    character(LEN=*),intent(in)::NameArray
    !EOP
    !BOC
    if (iError > 0) call CON_stop('check_allocate F90_ERROR '// &
         'Could not allocate array '//NameArray)
    !EOC
  end subroutine check_allocate

  !============================================================================
  subroutine test_mod_utility

    ! Test split_string, read a string containing separators 
    ! then print substrings
    ! Do this multiple times with various settings
   
    character(len=500):: String
    integer, parameter :: MaxString = 20
    integer :: nString
    character(len=30) :: String_I(MaxString)  
    integer :: iString
    integer:: iError

    character(len=*), parameter :: NameSub = 'test_mod_utility'
    !-----------------------------------------------------------------------
    write(*,'(a)') 'testing check_dir'
    write(*,'(a)') 'check directory "."'
    call check_dir('.')
    write(*,'(a)') 'check_dir returned successfully'
    write(*,'(a)') 'check directory "xxx/"'
    call check_dir('xxx/')

    write(*,'(a)') 'testing make_dir'
    write(*,'(a)') 'make directory "xxx/"'
    call make_dir('xxx')
    call check_dir('xxx/')
    write(*,'(a)') 'Making directory again (should produce "File exists" error)'
    call make_dir('xxx', iErrorOut=iError)
    write(*,*) iError

    write(*,'(/,a)') 'testing fix_dir_name'
    String = ' '
    call fix_dir_name(String)
    write(*,'(a)') 'fixed empty string=' // trim(String)
    String = 'GM/BATSRUS/data'
    write(*,'(a)') 'original    string=' // trim(String)
    call fix_dir_name(String)
    write(*,'(a)') 'fixed first string=' // trim(String)
    call fix_dir_name(String)
    write(*,'(a)') 'fixed again string=' // trim(String)

    write(*,'(/,a)') 'testing split_string'
    String = '  a(3)  bb(04:06) c,ddd ee,ff gg(8:12:2) '
    write(*,'(a)') 'String=' // trim(String)

    call split_string(String, MaxString, String_I, nString)
    write(*,'(a,i3,a)') 'with space separator split to', nString, ' parts:'
    do iString = 1, nString
       write(*,'(a)') trim(String_I(iString))
    end do

    call split_string(String, String_I, nString, ",")
    write(*,'(a,i3,a)') 'with "," separator split to', nString, ' parts:'
    do iString = 1, nString
       write(*,'(a)') trim(String_I(iString))
    end do

    call split_string(String, String_I, nString, ") ", DoAddSeparator=.true.)
    write(*,'(a,i3,a)') 'with ") " separator split to', nString, ' parts:'
    do iString = 1, nString
       write(*,'(a)') trim(String_I(iString))
    end do

    call split_string(String, String_I, nString, UseArraySyntaxIn=.true.)
    write(*,'(a,i3,a)') 'with UseArraySyntax split to', nString,' parts:'
    do iString = 1, nString
       write(*,'(a)') trim(String_I(iString))
    end do

    write(*,'(/,a)') 'testing join_string'
    call join_string(nString, String_I, String)
    write(*,'(a)') 'joined string='//trim(String)

    call join_string(String_I(1:4), String, '; ')
    write(*,'(a)') 'joined string='//trim(String)

    write(*,'(/,a)') 'testing upper_case and lower_case'
    String = 'abCD 123:'
    write(*,'(a)') 'mixed case string='//trim(String)
    call upper_case(String)
    write(*,'(a)') 'upper case string='//trim(String)
    call lower_case(String)
    write(*,'(a)') 'lower case string='//trim(String)


  end subroutine test_mod_utility

end module ModUtilities
