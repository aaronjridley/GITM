!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine EIE_End(iError)

  use ModEIE_Interface

  integer, intent(out) :: iError

  iError = 0
  if (allocated(EIEr3_HaveMLTs)) deallocate(EIEr3_HaveMLTs, stat = iError)
  if (allocated(EIEr3_HaveLats)) deallocate(EIEr3_HaveLats, stat = iError)
  if (allocated(EIEr3_HavePotential)) &
       deallocate(EIEr3_HavePotential, stat = iError)
  if (allocated(EIEr3_HaveEFlux)) deallocate(EIEr3_HaveEFlux, stat = iError)
  if (allocated(EIEr3_HaveAveE)) deallocate(EIEr3_HaveAveE, stat = iError)

  stop

end subroutine EIE_End
