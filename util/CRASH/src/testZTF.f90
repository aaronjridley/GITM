!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

 !------
	program testNLTE
	use M_TEST
	implicit none
	real(8) :: te,Etot,Ptot,zbar,cv

	numAto=79
	Atomass=179
	Efloor=1e-5
	Pcold=0
 
 1	print*,'ro,te,0  /or/  ro,0,Ee ?'
	read(*,*) ro,te,Etot
	if(ro.le.0) stop 'normal end'
	if(te.le.0 .and. Etot.le.0) stop 'Normal End'

 ! density dependant , will be passed thorugh MODULE:
	kBr_E= ro * ERGperEV
	kBr_P= ro * DYNEperEV

	if(Etot.eq.0) then
		print*,'call LTE_EOS_dir : ro,te=',ro,te
	 call LTE_EOS_dir(te,Etot,Ptot,Zbar,Cv)		! ro : in module m_TEST
	elseif(te.eq.0) then
		print*,'call LTE_EOS_inv : ro,ee=',ro,etot
	 call LTE_EOS_inv(te,Etot,Ptot,Zbar,Cv)		! ro : in module m_TEST
	endif
	print*,'ro=',ro,' te,ee,pe=',te,etot,ptot,' Z,Cv=',zbar,cv
	goto 1
	
	end
