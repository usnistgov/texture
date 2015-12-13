c * pole figure projection
      subroutine agr2pol(agrain,miller,pole_sa)
      implicit None
      real*8 agrain(4), miller(3), pole_sa(3)
      integer i,j
      real*8 norm, ph,th,tm,r(3,3),ra(3,3)
Cf2py intent(agrain) agrain, miller
Cf2py intent(out) pole_sa

      ph = agrain(1)
      th = agrain(2)
      tm = agrain(3)
      call euler(1, ph,th,tm,r) ! ca<-sa
      do 10 i=1,3
      do j=1,e
         ra(i,j) = r(j,i)
 10   continue
      norm = 0.
      do i=1,3
         norm = norm + miller(i)**2
      enddo
      norm = sqrt(norm)
      do i=1,3
         miller(i) = miller(i) / norm
      enddo
      pole_sa(:) = 0.
      do 10 i=1,3
      do 20 j=1,3
         pole_sa(i) = pole_sa(i) + ra(i,j) * miller(j)
 20   continue
      return
      end subroutine


c *****************************************************************************
      subroutine euler (iopt,ph,th,tm,a)
c
c     CALCULATE THE EULER ANGLES ASSOCIATED WITH THE TRANSFORMATION
c     MATRIX A(I,J) IF IOPT=1 AND VICEVERSA IF IOPT=2
c     A(i,j) TRANSFORMS FROM SYSTEM sa TO SYSTEM ca.
c     ph,th,om ARE THE EULER ANGLES (in degrees) OF ca REFERRED TO sa.
c *****************************************************************************

      dimension a(3,3)
      pi=4.*atan(1.d0)

      if(iopt.eq.1) then
        th=acos(a(3,3))
        if(abs(a(3,3)).ge.0.9999) then
          tm=0.
          ph=atan2(a(1,2),a(1,1))
        else
          sth=sin(th)
          tm=atan2(a(1,3)/sth,a(2,3)/sth)
          ph=atan2(a(3,1)/sth,-a(3,2)/sth)
        endif
        th=th*180./pi
        ph=ph*180./pi
        tm=tm*180./pi
      else if(iopt.eq.2) then
        sph=sin(ph*pi/180.)
        cph=cos(ph*pi/180.)
        sth=sin(th*pi/180.)
        cth=cos(th*pi/180.)
        stm=sin(tm*pi/180.)
        ctm=cos(tm*pi/180.)
        a(1,1)=ctm*cph-sph*stm*cth
        a(2,1)=-stm*cph-sph*ctm*cth
        a(3,1)=sph*sth
        a(1,2)=ctm*sph+cph*stm*cth
        a(2,2)=-sph*stm+cph*ctm*cth
        a(3,2)=-sth*cph
        a(1,3)=sth*stm
        a(2,3)=ctm*sth
        a(3,3)=cth
      endif
      return
      end
