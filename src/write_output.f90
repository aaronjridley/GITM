!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!-----------------------------------------------------------------------------
! $Id: write_output.f90,v 1.16 2013/10/12 04:01:00 kopmanis Exp $
!
! Author: Aaron Ridley, UMichigan
!
! Comments: Routines to write data output (both the desired model output files
!           and satellite output).  Also calls the logfile output routine.
!-----------------------------------------------------------------------------

subroutine write_output

  use ModTime
  use ModInputs
  use ModGITM

  implicit none

  real, external :: get_timing
  real :: ProjectedTime, CompletedTime, RealTime
  integer :: i, nMLTsTmp,nLatsTmp, iBlock
  logical :: IsDone

  if (floor((tSimulation-dt)/DtReport) /= &
       floor((tsimulation)/DtReport) .and. iDebugLevel >= 0) then
     if (IsFramework) then
        if(iProc==0)write(*,"(a,i6,a,3i2.2)") "UA:GITM2 iStep ", iStep, &
             ", Time : ",iTimeArray(4:6)
     else
        RealTime = get_timing("GITM")
        CompletedTime = (EndTime-CurrentTime)/(CurrentTime-RestartTime)
        ProjectedTime = RealTime * CompletedTime
        write(*,"(a,i6,a,3i2.2,a,f10.2,a,f10.2)") "iStep ", iStep, &
             ", Time : ",iTimeArray(4:6), &
             ", RealTime : ",RealTime/60.0," min, Projected : ",&
             ProjectedTime/60.0
     endif
  endif

  DtPlot = DtPlotSave
  if ( CurrentTime >= PlotTimeChangeStart .and. &
       CurrentTime <= PlotTimeChangeEnd) then 
     DtPlot = PlotTimeChangeDt
  endif

  IsDone = .false.
  do i = 1, nOutputTypes
     if (floor((tSimulation-dt)/DtPlot(i)) /= &
          floor((tsimulation)/DtPlot(i)) .or. tSimulation == 0.0) then
        do iBlock = 1, nBlocks
           call output("UA/data/",iBlock, i)
        enddo
     endif
  enddo

  call move_satellites

  if (floor((tSimulation-dt)/DtRestart) /= &
       floor((tsimulation)/DtRestart)) then
     call write_restart("UA/restartOUT/")
  endif

  if (floor((tSimulation-dt)/DtLogfile) /= &
       floor((tsimulation)/DtLogfile)) then
     call logfile("UA/data/")
  endif

end subroutine write_output

