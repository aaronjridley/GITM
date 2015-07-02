!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

	subroutine bracket_tz  ( ro ) ! ... )




	implicit none
	logical :: bracketed,converged
	real(8),intent(IN) :: ro
	real(8) :: facp,facm,teplus,teminus,teplus_last,teminus_last &
	   ,ftstp,ftstm ,ftstp_last,ftstm_last ,Tz,eff_Ee,e,tztst &
	   ,temin,temax ,cvloc,tzold,teold,fold,te_LTE ,dtzmin
	real(8) :: enflor=0,teflor=1
	integer :: fail
	real(8),parameter :: half=0.5d0,zero=0,one=1,two=2

	 bracketed = .false.
	 converged = .false.
	 facp = 2.0
	 facm = 0.5

	 teplus  = teold
	 teminus = teold
	 ftstm   = fold
	 teplus_last  = teold	! 090506
	 teminus_last = teold	! 090506
	 ftstp_last   = fold
	 ftstm_last   = fold


	 do while (.not. bracketed .and. &
            (teplus_last .lt.temax .or. teminus_last.gt.teflor))


	  teplus  = min(max(te_LTE,facp*teold),temax)	! 090506

	  if(teplus.gt.teplus_last)  then

	    eff_Ee = e - cvloc*(teplus-tztst)
	   do while (eff_Ee.le.enflor)
	    teplus= 0.75*teplus
	    tztst = 0.75*tztst
	    eff_Ee = e - cvloc*(teplus-tztst)
	   enddo 


            ftstp = tztst - Tz

            if(abs(ftstp).le.dtzmin*half*(tztst+Tz)) then


              converged = .true.
              bracketed = .true.
              teold     = teplus
              tzold     = half*(tztst+Tz)	! 090506
              Tz   = tzold			! 090506
              fail      = 0

            elseif((ftstp*ftstp_last).lt. zero)  then


              bracketed = .true.
              teminus   = teplus_last
              ftstm     = ftstp_last

            else


              ftstp_last  = ftstp
              teplus_last = teplus

            endif ! test on ftstp

          endif ! i have not run out the end on the plus side

          if(.not.bracketed) then


            teminus  = max(facm*teold,teflor)

            if(teminus.lt.teminus_last)  then

              eff_Ee = e - cvloc*(teminus-tztst)
               do while (eff_Ee.le.enflor)
                teminus= 0.75*teminus
                tztst = 0.75*tztst
                eff_Ee = e - cvloc*(teminus-tztst)
               enddo 


              ftstm = tztst - Tz

              if(abs(ftstm).le.dtzmin*half*(tztst+Tz)) then


                converged = .true.
                bracketed = .true.
                teold     = teminus
                tzold     = half*(tztst+Tz)	! 090506
                Tz   = tzold			! 090506
                fail      = 0

              elseif (ftstm_last*ftstm .lt. zero)  then


                bracketed = .true.
                teplus  = teminus_last
                ftstp   = ftstm_last

              else


                ftstm_last   = ftstm
                teminus_last = teminus
              endif  ! test on ftstm

            endif  ! i haven t run out the end from the minus side

          endif  ! bracketed from plus side


          facp = two*facp
          facm = half*facm

        enddo  ! while not bracketed

	return
	end subroutine bracket_tz
