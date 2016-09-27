c -------------------------------------------------------
      subroutine projection(pole,x,y)
      real*8 pole(3), x, y, norm
      integer i
cf2py intent(in) pole
Cf2py intent(out) x,y
      norm = 0.
      do 10 i=1,3
         norm = norm + pole(i)**2
 10   continue
      norm = dsqrt(norm)
      do 20 i=1,3
         pole(i) = pole(i) / norm
 20   continue
      if (pole(3).eq.1) then
         X=0.d0
         Y=0.d0
      else
         X = pole(1)/(pole(3)-1.)
         y = pole(2)/(pole(3)-1.)
      endif
      return
      end subroutine

c -------------------------------------------------------
c     Based on grain's orientation <agrain> and <pole_ca>
c -------------------------------------------------------
c     * pole figure projection
      subroutine agr2pol(agrain,pole_ca,pole_sa)
      implicit None
      real*8 agrain(3), pole_ca(3), pole_sa(3)
      integer i,j
      real*8 norm, ph,th,tm,r(3,3),ra(3,3)
Cf2py intent(in) agrain, pole_ca
Cf2py intent(out) pole_sa
      ph = agrain(1)
      th = agrain(2)
      tm = agrain(3)
      call euler(2, ph, th, tm, r) ! ca<-sa
      do 10 i=1,3
      do 10 j=1,3
         ra(i,j) = r(j,i)
 10   continue
      norm = 0.
      do 12 i=1,3
         norm = norm + pole_ca(i)**2
 12   continue
      norm = sqrt(norm)
      do 15 i=1,3
         pole_ca(i) = pole_ca(i) / norm
 15   continue
      pole_sa(:) = 0.
      do 20 i=1,3
      do 20 j=1,3
         pole_sa(i) = pole_sa(i) + ra(i,j) * pole_ca(j)
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
      integer iopt
      real*8 ph,th,tm,a(3,3),pi,sth,sph,cph,cth,stm,ctm
Cf2py intent(in) iopt
Cf2py intent(out,in) ph,th,tm,a
      pi=4.*datan(1.d0)

      if(iopt.eq.1) then
        th=dacos(a(3,3))
        if(dabs(a(3,3)).ge.0.9999) then
          tm=0.
          ph=datan2(a(1,2),a(1,1))
        else
          sth=dsin(th)
          tm=datan2(a(1,3)/sth,a(2,3)/sth)
          ph=datan2(a(3,1)/sth,-a(3,2)/sth)
        endif
        th=th*180./pi
        ph=ph*180./pi
        tm=tm*180./pi
      else if(iopt.eq.2) then
        sph=dsin(ph*pi/180.)
        cph=dcos(ph*pi/180.)
        sth=dsin(th*pi/180.)
        cth=dcos(th*pi/180.)
        stm=dsin(tm*pi/180.)
        ctm=dcos(tm*pi/180.)
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

c *****************************************************************************
      subroutine grain2pole_sa(
     $     ngr,grains,npol,poles_ca,
     $     poles_sa,poles_wgt)
      integer ngr,npole,i,j,igr,n
      real*8 grains(ngr,4),poles_ca(npol,3),
     $     poles_sa(ngr,npol,3),poles_wgt(ngr,npol)
      real*8 phi1,phi,phi2,ag(3,3),aux3(3),agt(3,3)
cf2py intent(in) ngr, grains, npole, poles_ca
cf2py intent(out) poles_sa,poles_wgt
      do 10 igr=1,ngr
         phi1 = grains(igr,1)
         phi  = grains(igr,2)
         phi2 = grains(igr,3)
         call euler(2,phi1,phi,phi2,ag) ! ca<-sa R_ij p_i
         do 5 i=1,3
         do 5 j=1,3
            agt(i,j) = ag(j,i)
 5       continue
      do 10 n=1,npol
         do 20 i=1,3
            aux3(i) = 0.d0
         do 20 j=1,3
            aux3(i) = aux3(i) + agt(i,j) *
     $           poles_ca(n,j)
 20      continue
         poles_sa(igr,n,:) = aux3(:)
         poles_wgt(igr,n)  = grains(igr,4)
 10   continue
      return
      end subroutine
