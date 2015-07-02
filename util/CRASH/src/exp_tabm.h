!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
 ! *** exp_tabm.h ***
 ! *
 ! * IN : y , OUT: ey=exp(-Y), ex_y=1-ey , |Y|=min(|y|,ex_max) 
 ! *
	if(y>=0.0) then
	 ex_y = min(ex_max,y)
	 ex_u = ex_y*ex_sdu
	 ex_i = int(ex_u)
	 ex_u = ex_y -ex_i*ex_du
	 ey = ex_tab(ex_i)*(1.0 - ex_c*ex_u)
	else
	 ex_y = min(ex_max,-y)
	 ex_u = ex_y*ex_sdu
	 ex_i = int(ex_u)
	 ex_u = ex_y -ex_i*ex_du
	 ey = (1.0+ex_cm*ex_u)/ex_tab(ex_i)
	end if
	if(ex_i/=0) ex_y = 1.0 - ey
