!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!*** exp_tab.h ***
!*
!* IN : u , OUT: ey=exp(-u) , u=min(u,ex_max) ; 0 <= u 
!*
! ey=exp(-u)	! pour verification
! ex_y=ex_one-ey
!*		

	ex_u=min(ex_max, max(0.0, u))
  
	ex_i=int(ex_u*ex_sdu)

	ex_u = ex_u -ex_i*ex_du

	ey=ex_tab(ex_i)*(1.0 - ex_c*ex_u)
