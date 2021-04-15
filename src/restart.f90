!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine write_restart(dir)

  use ModGITM
  use ModInputs !alexey uses useDART and f107 variable from ModInputs as estimate updated by DART, f107a is assumed = f107
  use ModTime
  use ModSphereInterface, only: iStartBLK

  use ModIndicesInterfaces, only: get_f107            !alexey
  use ModMpi,               only: MPI_BCAST, MPI_REAL !alexey
  use ModUtilities,         only: sleep

  implicit none

  character (len=*), intent(in) :: dir
  character (len=5) :: cBlock

  integer :: iBlock, iSpecies, i

!alexey added from here
  integer :: iError !error for mpibcast statements
  integer :: ens_t, ens_c  !number of ensemble memebers total, ens. counter
  character :: ens_cs*5  !ensemble counter string (5 digits, so hopefully you don't do more than 10^5 ensembles :)
  real :: sigma_i2 !inflation value, as in \sigma_i^2 in eq (18) in GITM2/srcDoc/thermo.pdf and .tex, eg 49SFU^2, is a variance, not sd's
  real, allocatable :: f_tmp(:), f_tmp2(:) !f107 temporary storage variables (before inflation and after)
  real :: sum_f=0.0 !sum of f107 across the ensemble
  real :: sum_fsq=0.0 !sum of f107 squares
  real :: mean_f=0.0 !mean of f107 across ensemble
  real :: var_f=0.0 !variance of f107 across ens.
  logical :: IsThere=.false. !used to check if various files exist yet or not and suppressing writes
!to here to perform driver (f107) inflation

  call report("write_restart",1)

  if (iProc == 0) then

     open(unit=iRestartUnit_, file=dir//"header.rst", status="unknown")

     write(iRestartUnit_, *) ""
     write(iRestartUnit_, '(a)') "#ISTEP"
     write(iRestartUnit_, *) iStep

     write(iRestartUnit_, *) ""
     write(iRestartUnit_, '(a)') "#TSIMULATION"
     write(iRestartUnit_, *) tSimulation

     write(iRestartUnit_, *) ""
     write(iRestartUnit_, '(a)') "#TIMESTART"
     do i=1,6
        write(iRestartUnit_, *) iStartTime(i)
     enddo

     write(iRestartUnit_, *) ""
     write(iRestartUnit_, '(a)') "#SPHERE"
     write(iRestartUnit_, *) Is1D
     write(iRestartUnit_, *) IsFullSphere
     write(iRestartUnit_, *) LatStart
     write(iRestartUnit_, *) LatEnd
     write(iRestartUnit_, *) LonStart

     write(iRestartUnit_, *) ""
     write(iRestartUnit_, '(a)') "#END"
     write(iRestartUnit_, *) ""

     close(iRestartUnit_)

  endif

  do iBlock = 1, nBlocks

     write(cBlock,'(a1,i4.4)') "b",iBlock+iStartBLK
     open(unit=iRestartUnit_, file=dir//"/"//cBlock//".rst", &
          status="unknown", form="unformatted")

     write(iRestartUnit_) Longitude(:,iBlock)
     write(iRestartUnit_) Latitude(:,iBlock)
     if(UseTopography)then
        write(iRestartUnit_) Altitude_GB(:,:,:,iBlock)
     else
        write(iRestartUnit_) Altitude_GB(1,1,:,iBlock)
     end if

     do iSpecies=1,nSpeciesTotal
        write(iRestartUnit_) NDensityS(:,:,:,iSpecies,iBlock)
     enddo

     do iSpecies=1,nIons
        write(iRestartUnit_) IDensityS(:,:,:,iSpecies,iBlock)
     enddo

!alexey start
     MeanMajorMass(:,:,:)=0  
     do iSpecies = 1, nSpecies  
        MeanMajorMass(:,:,:) = MeanMajorMass(:,:,:) +   &  
             Mass(iSpecies) * NDensityS(:,:,:,iSpecies,iBlock)/   & 
             sum(NDensityS(:,:,:,1:3,iBlock),4) 
     enddo
     TempUnit(:,:,:) = &  
          MeanMajorMass(:,:,:)/ Boltzmanns_Constant 
     write(iRestartUnit_) Temperature(:,:,:,iBlock)*TempUnit(:,:,:) !alexey writes out temperature in Kelvin instead of scaled unit
!alexey end

     write(iRestartUnit_) ITemperature(:,:,:,iBlock)
     write(iRestartUnit_) eTemperature(:,:,:,iBlock)

     write(iRestartUnit_)  Velocity(:,:,:,:,iBlock)
     write(iRestartUnit_) IVelocity(:,:,:,:,iBlock)
     write(iRestartUnit_) VerticalVelocity(:,:,:,:,iBlock)

!alexey start of driver inflation code
     call get_f107(CurrentTime, f107, iError)

     !if DART is not used (useDART=0), then skip the following if statement
     if (useDART == 1) then !this is the DART master ensemble member - collect f107 from other ensemble members via text files (a.txt),
        !inflate f107 distribution, and distribute the new values via text files (b.txt)
        if (iProc==0) then !only process zero needs to do this, others will receive this via MPI_BCAST
           open(unit=103,file='a.txt')
           write(103,*) f107 
           close(103)

           open(unit=103,file='here.txt') !signify to the master that a.txt is written
           write(103,*) "I am here" !just so that file has something in it
           close(103)

           !read f107 from all of the ensemble members, inflate it, write it back out
           !must be in the advance_temp_e1 (the first ens. member run-directory)
           !other ens. members must have some kind of waiting-code

           open(unit=103,file='ens_size.txt',iostat=iError)  !find out how big the ensemble is 
           read(103,*) ens_t
           close(103)
           allocate(f_tmp(ens_t)) !allocate memory based on ens size
           allocate(f_tmp2(ens_t))
           if (iDebugLevel > 4) write(*,*) '=====> rst.f90: opened ens_size.txt with iostat, #_of_ens: ', iError, ens_t

           open(unit=103,file='sigma_i2.txt') !find out how much to inflate
           read(103,*) sigma_i2
           close(103)

           do ens_c=1,ens_t !cycle through all ensemble members to get each f107 estimate
              write(ens_cs,'(I5)') ens_c !convert ensemble counter integer into a string
              inquire(file='../advance_temp_e'//trim(adjustl(ens_cs))//'/here.txt',EXIST=IsThere) !see if this slave ready
              !adjustl moves the whole string left, trim trims the spaces on the right. 
              ! this is needed because string '1' has different length from string '11'

              do while (.not.IsThere) !wait for it to get ready
                 write(*,*) '=====> rst.f90: here.txt file NOT found. Waiting in advance_temp_e' & 
                      //trim(adjustl(ens_cs)) 
                 call sleep(1.0)
                 inquire(file='../advance_temp_e'//trim(adjustl(ens_cs))//'/here.txt',EXIST=IsThere)
              enddo

              if (iDebugLevel > 4) write(*,*) '=====> rst.f90: here.txt file found. Continuing  in advance_temp_e' &
                   //trim(adjustl(ens_cs))
              open(103, file = '../advance_temp_e'//trim(adjustl(ens_cs))//'/here.txt', status = 'OLD') !remove the readiness file
              close(103, status = 'DELETE')
              if (iDebugLevel > 4) write(*,*) '=====> rst.f90: here.txt removed in advance_temp_e'//trim(adjustl(ens_cs))

              open(unit=103,file='../advance_temp_e'//trim(adjustl(ens_cs))//'/a.txt', status = 'OLD') !read this ensemble member's estimate
              read(103,*) f_tmp(ens_c)
              close(103, status = 'DELETE')
              if (iDebugLevel > 4) write(*,*) '=====> rst.f90: f_tmp', f_tmp 
           enddo

  
           mean_f = samp_mean(f_tmp)
           var_f = samp_var(f_tmp)
           f_tmp2=sqrt(sigma_i2/var_f)*(f_tmp-mean_f)+mean_f !const variance, usually used for driver estimates
           !+ by constant variance I mean that the user picks a number he or she desires for the ensemble variance 
           !+ of the driver estimate. For my initial experiments with f107 the value of 7 SFU (this is standard deviation, 
           !+ so variance is 49 SFU^2) worked well when average F107 value was 140 SFU (so we are talking something like 5%)
           !+ GITM provides no dynamics for the drivers (just because it doesn't model them and hence doesn't update them).
           !+ Consequently, it is beneficial to make EAKF driver estimate equally uncertain at every step, and hence 
           !+ constantly keep searching for it and never let it get stuck at some fixed value.
           !+ Anyways, the formula above removes whatever spread the driver estimate had right now (hence division by var_f)
           !+ and replaces it with the desired spread (hence the multiplication by sigma_i2).
           if (iDebugLevel > 4) write(*,*) '=====> rst.f90: mean_f', mean_f
           if (iDebugLevel > 4) write(*,*) '=====> rst.f90: var_f', var_f
           if (iDebugLevel > 4) write(*,*) '=====> rst.f90: f_tmp2', f_tmp2


           do ens_c=1,ens_t !cycle through all ensemble members to get each f107 estimate
              write(ens_cs,'(I5)') ens_c !convert ensemble counter integer into a string
              open(unit=103,file='../advance_temp_e'//trim(adjustl(ens_cs))//'/b.txt',iostat=iError) 
              !write this ensemble member's inflated est.
              write(103,*) f_tmp2(ens_c)
              close(103)
              if (iDebugLevel > 4) write(*,*) '=====> rst.f90: wrote b.txt in ../advance_temp_e' & 
                   //trim(adjustl(ens_cs))//' with iostat of', iError
              
              open(unit=103,file='../advance_temp_e'//trim(adjustl(ens_cs))//'/here2.txt',iostat=iError) 
              !tell slave that inflated file is written
              write(103,*) 'Your inflated estimate is ready' !just so that file has something in it
              close(103)
              if (iDebugLevel > 4) write(*,*) '=====> rst.f90: wrote here2.txt in ../advance_temp_e' & 
                   //trim(adjustl(ens_cs))//' with iostat of', iError
           enddo

           f107=f_tmp2(1) !master should read his inflated f107
           !write(101,*) f107 !write f107 into the f107_out.txt file (for plotting) 
           if (iDebugLevel > 4) write(*,*) '=====> rst.f90: write_rst_f107', f107 
           
           deallocate(f_tmp) !clean up
           deallocate(f_tmp2)
        endif
        call MPI_BCAST(f107,1,MPI_Real,0,iCommGITM,iError) !+ this is b-casted so that all blocks have the same f107 
        if (iDebugLevel > 4) write(*,*) '=====> rst.f90: broadcasted to iproc, f107:',iProc, f107

     

     elseif (useDART == 2) then !this is the DART slave ensemble member - write out my estimate to a.txt, signal that it is written by
        !+ creating here.txt, wait for creation of a return handshake (here2.txt), signifying that the inflated file b.txt is ready
        if (iProc==0) then !only process zero needs to do this, others will receive this via MPI_BCAST
           open(unit=103,file='a.txt') !write my estimate
           write(103,*) f107 
           close(103)

           open(unit=103,file='here.txt') !signify to the master that my a.txt is written
           write(103,*) 'I am here' !just so that file has something in it
           close(103)

           inquire(file='here2.txt',EXIST=IsThere) !see if my inflater file is ready
           
           do while (.not.IsThere) !wait for it to get ready
              write(*,*) '=====> rst.f90: here2.txt file NOT found. Waiting.' 
              call sleep(1.0)
              inquire(file='here2.txt',EXIST=IsThere)
           enddo

           if (iDebugLevel > 4) write(*,*) '=====> rst.f90: here2.txt file found. Continuing.'
           open(103, file = 'here2.txt', status = 'OLD') !remove the readiness file
           close(103, status = 'DELETE')
           if (iDebugLevel > 4) write(*,*) '=====> rst.f90: here2.txt removed'
           
           open(unit=103,file='b.txt', status = 'OLD') !read my new (inflated) estimate
           read(103,*) f107
           close(103, status = 'DELETE')
           !write(101,*) f107 !write f107 into the f107_out.txt file (for plotting) 
           if (iDebugLevel > 4) write(*,*) '=====> rst.f90: write_rst_f107', f107 
           
        endif
        call MPI_BCAST(f107,1,MPI_Real,0,iCommGITM,iError) !+ this is b-casted so that all blocks have the same f107 
        if (iDebugLevel > 4) write(*,*) '=====> rst.f90: broadcasted to iproc, f107:',iProc, f107
     

     endif

     write(iRestartUnit_) f107 !alexey write f107 into the restart file (for DART)
     write(iRestartUnit_) Rho(:,:,:,iBlock) !alexey
!alexey end

     if (isMars) then
        write(iRestartUnit_) SurfaceTemp(:,:,iBlock)
        write(iRestartUnit_) SubSurfaceTemp(:,:,iBlock)
        write(iRestartUnit_) dSurfaceTemp(:,:,iBlock)
        write(iRestartUnit_) dSubSurfaceTemp(:,:,iBlock)
     endif
     close(iRestartUnit_)

  enddo

!routines alexey needed
contains 
  function samp_mean (array) !sample mean
    implicit none
    real, intent(in), dimension(:) :: array
    real :: samp_mean
    samp_mean = sum(array)/size(array)
  end function samp_mean

  function samp_var (array) !unbiased sample variance
    implicit none
    real, intent(in), dimension(:) :: array
    real :: samp_var, samp_mean
    samp_mean = sum(array)/size(array)
    samp_var = sum((array-samp_mean)**2)/(size(array)-1)
  end function samp_var
!end of alexey's routines

end subroutine write_restart

!=============================================================================

subroutine read_restart(dir)

  use ModGITM
  use ModInputs !alexey uses useDART and f107 variable as estimate updated by DART, f107a is set = f107 if useDART > 0
  use ModTime
  use ModSphereInterface, only: iStartBLK
 
  implicit none

  character (len=*), intent(in) :: dir
  character (len=5) :: cBlock
  integer :: iBlock, i, iSpecies, iAlt
  !---------------------------------------------------------------------------
  call report("read_restart",1)

  RestartTime = CurrentTime

  do iBlock = 1, nBlocks

     write(cBlock,'(a1,i4.4)') "b",iBlock+iStartBLK
     if (iDebugLevel > 2) write(*,*) "===> Reading block ",cBlock

     open(unit=iRestartUnit_, file=dir//"/"//cBlock//".rst", &
          status="old", form="unformatted")

     if (iDebugLevel > 4) write(*,*) "=====> Reading Longitude"
     read(iRestartUnit_) Longitude(:,iBlock)
     if (iDebugLevel > 4) write(*,*) "=====> Reading Latitude"
     read(iRestartUnit_) Latitude(:,iBlock)
     if (iDebugLevel > 4) write(*,*) "=====> Reading Altitude"

     if(UseTopography)then
        read(iRestartUnit_) Altitude_GB(:,:,:,iBlock)
     else
        read(iRestartUnit_) Altitude_GB(1,1,:,iBlock)
        do iAlt = -1, nAlts+2
           Altitude_GB(:,:,iAlt,iBlock) = Altitude_GB(1,1,iAlt,iBlock)
        end do
     end if
!     read(iRestartUnit_) Altitude_GB(:,:,:,iBlock)

     do iSpecies=1,nSpeciesTotal
        if (iDebugLevel > 3) &
             write(*,*) "====> Reading Species",iSpecies, Mass(iSpecies)
        read(iRestartUnit_) NDensityS(:,:,:,iSpecies,iBlock)
     enddo

     Rho(:,:,:,iBlock) = 0.0
     NDensity(:,:,:,iBlock) = 0.0
     do iSpecies=1,nSpecies
        Rho(:,:,:,iBlock) = Rho(:,:,:,iBlock) + &
          Mass(iSpecies)*NDensityS(:,:,:,iSpecies,iBlock)
        NDensity(:,:,:,iBlock) = NDensity(:,:,:,iBlock) + &
          NDensityS(:,:,:,iSpecies,iBlock)
     enddo

     do iSpecies=1,nIons
        if (iDebugLevel > 4) &
             write(*,*) "=====> Reading Ion Species",iSpecies, Mass(iSpecies)
        read(iRestartUnit_) IDensityS(:,:,:,iSpecies,iBlock)
     enddo

     if (iDebugLevel > 4) write(*,*) "=====> Reading Temperature"
     read(iRestartUnit_)  Temperature(:,:,:,iBlock)
     
!alexey start
     MeanMajorMass(:,:,:)=0  
     do iSpecies = 1, nSpecies  
        MeanMajorMass(:,:,:) = MeanMajorMass(:,:,:) +   &  
             Mass(iSpecies) * NDensityS(:,:,:,iSpecies,iBlock)/   &  
             sum(NDensityS(:,:,:,1:3,iBlock),4)  
     enddo  
     TempUnit(:,:,:) = &  
          MeanMajorMass(:,:,:)/ Boltzmanns_Constant  
     Temperature(:,:,:,iBlock) =  Temperature(:,:,:,iBlock)/TempUnit(:,:,:) 
!alexey end

     if (iDebugLevel > 4) write(*,*) "=====> Reading ITemperature"
     read(iRestartUnit_) ITemperature(:,:,:,iBlock)
     if (iDebugLevel > 4) write(*,*) "=====> Reading eTemperature"
     read(iRestartUnit_) eTemperature(:,:,:,iBlock)

     if (iDebugLevel > 4) write(*,*) "=====> Reading Velocity"
     read(iRestartUnit_)  Velocity(:,:,:,:,iBlock)
     if (iDebugLevel > 4) write(*,*) "=====> Reading IVelocity"
     read(iRestartUnit_) IVelocity(:,:,:,:,iBlock)
     if (iDebugLevel > 4) write(*,*) "=====> Reading VerticalVelocity"
     read(iRestartUnit_) VerticalVelocity(:,:,:,:,iBlock)

!alexey start
     read(iRestartUnit_) f107
     if (iDebugLevel > 4) write(*,*) "=====> Just read f107 = ", f107 

     if (useDART .ne. 0) then !this is a DART ensemble member - put f107estimate into f107 and f107a
        call IO_set_f107_single(f107) !replace f107 with the estimate
        call IO_set_f107a_single(f107) !replace f107 81-day average with the SAME estimate
     endif

     read(iRestartUnit_) Rho(:,:,:,iBlock) !read Rho, but hope DART will affect GITM via other vars instead
!alexey end

     if (isMars) then
        if (iDebugLevel > 4) write(*,*) "=====> Reading SurfaceTemp"
        read(iRestartUnit_) SurfaceTemp(:,:,iBlock)
        read(iRestartUnit_) SubSurfaceTemp(:,:,iBlock)
        read(iRestartUnit_) dSurfaceTemp(:,:,iBlock)
        read(iRestartUnit_) dSubSurfaceTemp(:,:,iBlock)
     endif
     close(iRestartUnit_)

  enddo

  if (.not. Is1D) call exchange_messages_sphere

end subroutine read_restart
