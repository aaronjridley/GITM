!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine divergence(InArray, OutArray, iBlock)

  ! This routine calculates the divergence of a scalar quantity.
  ! It assumes that it is working with a single block, but needs the
  ! variable "iBlock" because it has to figure out where it is in the
  ! grid.  This can be generalized by feeding in the specific grid
  ! for the given variable.  It is also assumed that the "InArray"
  ! and "OutArray" have the ghostcells defined.

  use ModSizeGitm
  use ModGITM, only: Latitude, CosLatitude, Longitude, &
       RadialDistance_GB, InvRadialDistance_GB
  use ModConstants

  implicit none

  real, dimension(-1:nLons+2,           &
                  -1:nLats+2,            &
                  nAlts), intent(in)  :: InArray

  real, dimension(-1:nLons+2,           &
                  -1:nLats+2,            &
                  nAlts), intent(out) :: OutArray

  real, dimension(-1:nLons+2,           &
                  -1:nLats+2,            &
                  nAlts) :: drm, drp, drmodrp, drr2, dr, &
                            dtm, dtp, dtmodtp, dtr2, dt

  integer, intent(in) :: iBlock

  real :: maxi, maxt, cosmaxt, tanmaxt

  integer :: i,j,k
  !----------------------------------------------------------------------------

  maxt = 80.0 * pi / 180.0
  cosmaxt = cos(maxt)
  tanmaxt = tan(maxt)

  do i = 0, nLons+1
     do j = 0, nLats+1

        k = 1
        dr(i,j,k) = RadialDistance_GB(i,j,k+1,iBlock) &
             -      RadialDistance_GB(i,j,k  ,iBlock)

        do k=2,nAlts-1
           drp(i,j,k) = RadialDistance_GB(i,j,k+1,iBlock) &
                -       RadialDistance_GB(i,j,k  ,iBlock)
           drm(i,j,k) = RadialDistance_GB(i,j,k  ,iBlock) &
                -       RadialDistance_GB(i,j,k-1,iBlock)
           drmodrp(i,j,k) = drm(i,j,k)/drp(i,j,k)
           drr2(i,j,k) = drmodrp(i,j,k) * drmodrp(i,j,k)
        enddo

        k = nAlts
        dr(i,j,k) = RadialDistance_GB(i,j,k  ,iBlock) &
             -      RadialDistance_GB(i,j,k-1,iBlock)
        
        do k = 1, nAlts
           dtm(i,j,k) = Latitude(j,iBlock) - Latitude(j-1,iBlock)
           dtp(i,j,k) = Latitude(j+1,iBlock) - Latitude(j,iBlock)
           dtmodtp(i,j,k) = dtm(i,j,k) / dtp(i,j,k)
           dtr2(i,j,k) = dtmodtp(i,j,k) * dtmodtp(i,j,k)
        enddo

     enddo
  enddo

  !\
  ! East First  (same as the gradient)
  !/

  do i = 0, nLons+1
     do j = 0, nLats+1
        do k = 1, nAlts
           maxi = max(CosLatitude(j,iBlock),cosmaxt) ! this is 80 degrees...
           OutArray(i,j,k) = ((InArray(i+1,j,k) - InArray(i-1,j,k)) / &
                (Longitude(i+1,iBlock) - Longitude(i-1,iBlock)))/ &
                (maxi*RadialDistance_GB(i,j,k,iBlock))
        enddo
     enddo
  enddo

  !\
  ! North Second
  !/

  do i = 0, nLons+1
     do j = 0, nLats+1
        do k = 1, nAlts

           ! North has 2 terms - gradient term and tan term

           OutArray(i,j,k) = OutArray(i,j,k) + &
                ((dtr2(i,j,k)*InArray(i,j+1,k) - InArray(i,j-1,k) - &
                  (dtr2(i,j,k)-1)*InArray(i,j,k)) / &
                 (dtr2(i,j,k)*dtp(i,j,k) + dtm(i,j,k))) &
                 *InvRadialDistance_GB(i,j,k,iBlock)

           maxi = max(tan(Latitude(j,iBlock)),tanmaxt) ! this is 80 degrees...
           OutArray(i,j,k) = OutArray(i,j,k) + &
                maxi*InArray(i,j,k)*InvRadialDistance_GB(i,j,k,iBlock)

        enddo
     enddo
  enddo

  !\
  ! Up Third
  !/

  do i = 0, nLons+1
     do j = 0, nLats+1

        ! Up has 2 terms also - gradient and 2/r term

        k = 1
        OutArray(i,j,k) = OutArray(i,j,k) + &
             (InArray(i,j,k+1) - InArray(i,j,k)) / &
             ( RadialDistance_GB(i,j,k+1,iBlock) &
             - RadialDistance_GB(i,j,k,iBlock))

        OutArray(i,j,k) = OutArray(i,j,k) + &
             2.0 * InArray(i,j,k) * InvRadialDistance_GB(i,j,k,iBlock)

        do k = 2, nAlts-1

           OutArray(i,j,k) = OutArray(i,j,k) + &
                (drr2(i,j,k)*InArray(i,j,k+1) - InArray(i,j,k-1) - &
                 (drr2(i,j,k)-1)*InArray(i,j,k)) / &
                (drr2(i,j,k)*drp(i,j,k) + drm(i,j,k))

           OutArray(i,j,k) = OutArray(i,j,k) + &
                2.0 * InArray(i,j,k) * InvRadialDistance_GB(i,j,k,iBlock)

        enddo

        k = nAlts
        OutArray(i,j,k) = OutArray(i,j,k) + &
             (InArray(i,j,k) - InArray(i,j,k-1)) / &
             ( RadialDistance_GB(i,j,k  ,iBlock) &
             - RadialDistance_GB(i,j,k-1,iBlock))

        OutArray(i,j,k) = OutArray(i,j,k) + &
             2.0 * InArray(i,j,k) * InvRadialDistance_GB(i,j,k,iBlock)

     enddo
  enddo

end subroutine divergence
