
subroutine check_for_nans_ions(cMarker)

  use ModSizeGitm
  use ModGITM
  use ModPlanet

  implicit none

  character (LEN=*), intent(in) :: cMarker
  integer :: iLon, iLat, iAlt, iIon

  do iLon=-1,nLons+2
     do iLat=-1,nLats+2
        do iAlt=-1,nAlts+2
           do iIon=1,nIons
              if (isnan(iDensityS(iLon,iLat,iAlt,iIon,1))) then
                 write(*,*) 'Nan found in iDensityS : '
                 write(*,*) cMarker
                 write(*,*) iLon,iLat,iAlt,iProc,iIon
              endif
           enddo
        enddo
     enddo
  enddo


end subroutine check_for_nans_ions

subroutine check_for_nans_neutrals(cMarker)

  use ModSizeGitm
  use ModGITM
  use ModPlanet

  implicit none

  character (LEN=*), intent(in) :: cMarker
  integer :: iLon, iLat, iAlt, iNeu

  do iLon=-1,nLons+2
     do iLat=-1,nLats+2
        do iAlt=-1,nAlts+2
           do iNeu=1,nSpecies
              if (isnan(nDensityS(iLon,iLat,iAlt,iNeu,1))) then
                 write(*,*) 'Nan found in nDensityS : '
                 write(*,*) cMarker
                 write(*,*) iLon,iLat,iAlt,iProc,iNeu
              endif
           enddo
        enddo
     enddo
  enddo

end subroutine check_for_nans_neutrals

subroutine check_for_nans_temps(cMarker)

  use ModSizeGitm
  use ModGITM
  use ModPlanet

  implicit none

  character (LEN=*), intent(in) :: cMarker
  integer :: iLon, iLat, iAlt, iNeu

  do iLon=-1,nLons+2
     do iLat=-1,nLats+2
        do iAlt=-1,nAlts+2
           if (isnan(Temperature(iLon,iLat,iAlt,1))) then
              write(*,*) 'Nan found in Temperature : '
              write(*,*) cMarker
              write(*,*) iLon,iLat,iAlt,iProc
           endif
        enddo
     enddo
  enddo

end subroutine check_for_nans_temps
