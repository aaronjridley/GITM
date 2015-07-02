      real*8 function GGFF (eavg,temp,zz,ze)
c
c     Function calculates the free-free Gaunt factor
c     a la Karzas and Latter
c
c__Calling arguments:
c     eavg    (i)  photon energy in Kev
c     temp    (i)  material temperature in Kev
c     zz      (i)  ion charge
c     ze      (i)  number of free electrons
c
c__Functions/subroutines required:
c     LUF
c
      implicit none
      real*8 eavg,temp,zz,ze
c
      integer*4 lev,ne,iu,iq
      parameter (lev = 7, ne = 16)
      real*8 uv(ne),qv(lev),gfv(ne,lev)
      real*8 u,q,r
      integer*4 LUF
c
      data uv /
     1  1.0d-4, 1.0d-3, 1.0d-2, 1.0d-1, 3.0d-1, 1.0d+0,
     2  1.5d+0, 2.0d+0, 3.5d+0, 5.0d+0, 7.0d+0, 1.0d+1,
     3  1.5d+1, 2.0d+1, 2.5d+1, 3.0d+1/
      data qv / 1.d-3, 1.d-2, 1.d-1, 1.d+0, 1.d+1, 1.d+2, 1.d+3/
      data gfv / 5.5d+0,  4.2d+0,  3.0d+0,  1.9d+0, 1.34d+0,
     1  .88d+0,  .78d+0,  .67d+0,  .52d+0,  .45d+0,  .37d+0,
     2  .33d+0,  .27d+0,  .24d+0,  .22d+0,  .20d+0,  5.5d+0,
     3  4.2d+0,  3.0d+0,  1.9d+0, 1.42d+0,  .97d+0,  .88d+0,
     4  .75d+0,  .60d+0,  .53d+0,  .46d+0,  .39d+0,  .33d+0,
     5  .28d+0,  .26d+0,  .24d+0,  5.4d+0,  4.1d+0,  3.0d+0,
     6  1.9d+0, 1.54d+0, 1.16d+0, 1.05d+0,  .98d+0,  .82d+0,
     7  .75d+0,  .66d+0,  .59d+0,  .51d+0,  .46d+0,  .42d+0,
     8  .38d+0,  5.0d+0,  3.8d+0,  2.7d+0,  1.8d+0,  1.5d+0,
     9  1.3d+0, 1.25d+0,  1.2d+0, 1.15d+0, 1.08d+0, 1.04d+0,
     1  .98d+0,  .90d+0,  .84d+0,  .80d+0,  .76d+0,  4.4d+0,
     2 3.25d+0,  2.2d+0,  1.5d+0, 1.32d+0,  1.2d+0, 1.18d+0,
     3 1.17d+0, 1.15d+0, 1.14d+0,1.135d+0, 1.13d+0, 1.12d+0,
     4 1.11d+0,  1.1d+0, 1.09d+0,  3.8d+0,  2.7d+0,  1.8d+0, 
     5  1.3d+0, 1.16d+0,  1.1d+0,1.095d+0, 1.09d+0,1.085d+0,
     6 1.08d+0, 1.08d+0, 1.08d+0, 1.09d+0, 1.09d+0, 1.09d+0,
     7 1.09d+0,  3.5d+0,  2.3d+0,  1.5d+0, 1.06d+0,1.055d+0,
     8 1.05d+0,1.047d+0,1.045d+0,1.042d+0, 1.04d+0, 1.04d+0,
     9 1.04d+0, 1.05d+0, 1.05d+0, 1.05d+0, 1.05d+0/
c
      u = max(1.d-4,min(4.d+1,eavg/temp))
      q = max(1.d-4,min(2.d+3,1.36059d-2*zz*zz/temp))
      r = max(1.d+0,1./(1.02575d+0-1.25d-3*min(1.15d+2,temp)-
     1   9.d-4*min(3.d+2,eavg)))
      r = r+temp*ze/(5.11d+2*zz*zz)*(8.5d-1+1.5d+0*sqrt(u)+3.2d-1*u)
      iu = LUF(u,uv(2),14)
      iq = LUF(q,qv(2),5)
      GGFF = r/((uv(iu+1)-uv(iu))*(qv(iq+1)-qv(iq)))
     1      *(gfv(iu+1,iq+1)*(u-uv(iu))*(q-qv(iq))
     2      +gfv(iu,iq+1)*(uv(iu+1)-u)*(q-qv(iq))
     3      +gfv(iu+1,iq)*(u-uv(iu))*(qv(iq+1)-q)
     4      +gfv(iu,iq)*(uv(iu+1)-u)*(qv(iq+1)-q))
      return
      end
