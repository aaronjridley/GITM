! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!==============================================================================

subroutine read_f107(iOutputError)

  use ModKind
  use ModIndices
  implicit none

  integer, intent(out)     :: iOutputError

  real, dimension(6,MaxIndicesEntries) :: tmp_all
  real :: tmp(7)
  integer :: iError, nPts, iPt
  logical :: done
  
  ! One line of input
  character (len=iCharLenIndices_) :: line
  
  !------------------------------------------------------------------------
  call init_mod_indices

  iOutputError = 0

  open(LunIndices_, file=NameOfIndexFile, status="old", iostat = ierror)

  if (ierror.ne.0) then
     iOutputError = 1
     return
  endif

  do iPt=1, 15
     read(LunIndices_, *) line
  enddo
  
  done = .false.

  iPt = 1
  do while (.not.done)
     read(LunIndices_,'(f4.0,5(1x,f2.0),1x,f3.0)', iostat = ierror ) tmp
     tmp_all(1:5, iPt) = tmp(1:5)
     tmp_all(6, iPt) = tmp(7)
     if (ierror /= 0) then 
        done = .true.
     endif
     iPt = iPt + 1
  enddo
  npts = iPt - 1

  call Insert_into_Indices_Array(tmp_all, f107_)
  close(LunIndices_)
  
end subroutine read_f107
