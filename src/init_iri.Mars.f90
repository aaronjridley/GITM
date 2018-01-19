
subroutine init_iri

  implicit none

  call report("init_iri",1)

  write(*,*) "You can not use IRI with any planet except Earth!!!"
  write(*,*) "If you ARE running Earth, then make the code again, using"
  write(*,*) "configure Earth ; make"
  call stop_gitm("I can not continue...")

end subroutine init_iri
