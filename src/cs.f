c
c  put "! -*- f90 -*-" as the first liner  (e.g. when using f2py)
c----------------------------------------------------------------------c
c     Crystal Symmetry Library
c     written by Youngung Jeong (June 28, 2012)
c----------------------------------------------------------------------c
c     The below program 'du' was made for debugging purpose.
c     it can be used as, hopefully, a tutorial.
c$$$      program du
c$$$      implicit none
c$$$      real*8 n(3), b(3),  uniqset(3,2,48)
c$$$      integer isym, inb
c$$$      logical ipr
c$$$      ipr = .true.
c$$$      isym = 1
c$$$c     plane normal directions
c$$$      n(1) = 1.d0
c$$$      n(2) = 0.d0
c$$$      n(3) = 0.d0
c$$$c     directions on the planes.
c$$$      b(1) = 1.d0
c$$$      b(2) = 0.d0
c$$$      b(3) = 0.d0
c$$$c                   (n, b, isym, uniqset, inb, ipr)
c$$$      call uniqpnset(n, b, isym, uniqset, inb, .true.)
c$$$      end program du
c----------------------------------------------------------------------c
c     pn_dirc takes plane normal and directions and returns only uniq
c     combinations of the two.

      subroutine uniqpnset(n, b, isym, uniqset, inb, verbose)
      implicit none
      real*8 n(3), b(3), nuniq(3,48), buniq(3,48), uniqset(3,2,48),
     $     dum01(3), dum02(3), dum03(3,48), dum04(3,48), dum
      integer in, ib, i, j, isym, nnuniq, nbuniq, inb, inbc
      logical verbose, isperpen
cf2py intent(in) n, b, isym, verbose
cf2py intent(out) uniqset, inb
c     Given the plane normal direction and the direction on the plane,
c     Find the all possible crystallographically equivalent sets.

c     Find all equivalent, but unique planes.
c     Find all equivalent, but unique directions.

      if (verbose) write(*,*) 'Plane normal and direction'
      if (verbose) write(*,'(3f7.3, 2x, 3f7.3)') n, b
      if (isym.eq.1) then       ! cubic
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c        Given the vector v, finds crystallographically equivalent
c        vectors considering 'icen', the origin point symmetry flag.
c                        v, vuniq, n     , icen   , prt
         if (verbose) write(*,*) 'calling cubr. cubic_eqvect'
         call cubic_eqvect(n, nuniq, nnuniq, .false., .false.) ! plane
         if (verbose) then
            write(*,*) 'nnuniq', nnuniq
            do 8 in=1,nnuniq
            do 8 i=1,3
               write(*,'(3f7.3)') (nuniq(i,j),j=1,3)
 8          continue
         endif
         call cubic_eqvect(b, buniq, nbuniq, .true. , .false.) ! direct
         if (verbose) then
            write(*,*) 'nbuniq', nbuniq
            do 9 in=1, nbuniq
            do 9 i=1,3
               write(*,'(3f7.3)') (nuniq(i,j),j=1,3)
 9          continue
c            read(*,*)
         endif
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      else if(isym.eq.2) then   ! hexagonal
         write(*,*) "Nothing's made other than cubic yet... "
         stop
      else if(isym.eq.3) then   ! trigonal
         write(*,*) "Nothing's made other than cubic yet... "
         stop
      else if(isym.eq.4) then   ! tetragonal
         write(*,*) "Nothing's made other than cubic yet... "
         stop
      else if(isym.eq.5) then   ! orthorhombic
         write(*,*) "Nothing's made other than cubic yet... "
         stop
      else if(isym.eq.6) then   ! monoclinic
         write(*,*) "Nothing's made other than cubic yet... "
         stop
      else if(isym.eq.7) then   ! triclinic
         write(*,*) "Nothing's made other than cubic yet... "
         stop
      endif
c     read(*,*)
c     write(*,*)
c' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' c
c     print screen block
      if (verbose) then
         write(*,*) 'unique normals'
         write(*,*) 'nnuniq:', nnuniq, 'nbuniq', nbuniq
c         read(*,*)
         do i=1, nnuniq
            write(*,'(3f7.3)') (nuniq(j,i), j=1,3)
         enddo
         write(*,*) 'unique directions'
         do i=1, nbuniq
            write(*,'(3f7.3)') (buniq(j,i), j=1,3)
         enddo
c         read(*,*)
      endif
c' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' c
c     End of finding crystallographically equivalent vectors based on
c     the indiviudal structure's rotation symmtry groups.
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     Finds possible sets of plane-normal and its characteristics
c     direction. Note that the normality of suitable sets of them is
c     based by logical function isnorm.
      inb = 0
      do 1000 in=1, nnuniq    ! plane normal ends at the line # 1000
         do 10 i=1,3
            dum01(i) = nuniq(i,in)
 10      continue

      inbc = 0
      do 500 ib=1, nbuniq    ! direction
         do 20 i=1,3
            dum02(i) = buniq(i,ib) ! current direction
 20      continue
         if (verbose) write(*, '(3f7.3, 3x, 3f7.3)') (dum01(i), i=1,3),
     $        (dum02(i), i=1,3)
         dum = 0.d0
         do i=1,3
            dum = dum + dum01(i) * dum02(i)
         enddo
         write(*,*) 'dum = ', dum
         write(*,*) 'inbc= ', inbc
         if (isperpen(dum01, dum02, 3)) then
            inbc = inbc + 1     ! number of direction for this particular
c                                 plane (in) increases.
            do 30 i=1,3
c              Save those turned out to be on the plane dum01 to dum03.
c              For latter use.
               dum03(i,inbc) = dum02(i)
 30         continue
         endif
 500  continue ! over ib. It is not the end of loop over plane in.

c     if (verbose) write(*,*) 'inbc after putting things into dum03: ',
c$     inbc
c     if (verbose) read(*,*)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     Among the directions that turned out to be on the plane,
c     we want to remove the duplicated lines that are parallel (a=-b)
c         dum01(i) ! plane
c         dum03(i, ib)  ! the matched plane for
c                         the given current plane dum01
      if (inbc.eq.0) then
         write(*,*) 'inbc is wrong! at 153'
         stop
      endif
c
c     directions that found to be on the plane are in dum03 array.
c                         a -> , b    , n0  , n1, icen   , verbose
      call finduniquevect(dum03, dum04, inbc, ib, .false., .false.)
c     only unique directions are returned as dum04. and the number of
c     unique arrays are ib.
c' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' c
c     print screen block
      if (verbose) then
          write(*,*) 'before the removal of parallel vectors'
          do i=1,inbc
             write(*,'(3f7.3)') (dum03(j,i),j=1,3)
          enddo
          write(*,*) 'after'
          do i=1,ib
             write(*,'(3f7.3)') (dum04(j,i),j=1,3)
          enddo
          write(*,*) 'inbc: ', inbc
c          read(*,*)
       endif
c' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' c

      do 600 j=1,ib ! direction
         inb = inb + 1
      do 600 i=1,3
         uniqset(i,1,inb) = dum01(i)
         uniqset(i,2,inb) = dum04(i,j)
 600  continue                  ! over ib again.
 1000 continue                  ! loop over n & b set.
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      if (verbose) then
         write(*,*) 'inb:', inb
         do in=1,inb
            write(*,'(3f8.3, 2x, 3f8.3)') (uniqset(i,1,in),i=1,3),
     $           (uniqset(i,2,in),i=1,3)
         enddo
      endif
      return
      end subroutine uniqpnset
c----------------------------------------------------------------------c
c     crystal rotation symmetry operator calculation programs
c     the calculation of each element of operators are rigorously
c     implemented for facilitating general applications.
c----------------------------------------------------------------------c
      subroutine orthorhombic_rot_sym(h, n, icen)
c     Returns the rotation symmetry operators: h.
c     h_[ij] transforms a vector to another albeit crystallographically
c     equivalent one in ca for orthorhombic structure.
c     n returns the kinds of the transformation matrix h.
c     Symmetry operations:
c     1. 2-fold symmetry about [001]
c     2. Mirror by yz and xz planes, respectively.
c     3. Centro-symmetry by (0,0,0)
      implicit none
      real*8 h0(3,3), h(3,3,48), dum(3,3), hmiryz(3,3), hmirxz(3,3),
     $     hz180(3,3), hc(3,3), z(3), rz(3,3)
      integer i,j,n,in,iin
      logical icen
cf2py intent(in) icen
cf2py intent(out) h, n
      write(*,*) 'Not tested yet for orthorhombic_rot_sym'
      stop
      call imat(h0)
      do 10 i=1,3
      do 10 j=1,3
         h(i,j,1) = h0(i,j)
 10   continue
      call zmat(h0)
      call zvec(z)
      z(3) = 1.d0
      call vector_ang(z, 180.d0, hz180)
      call reflect(hmirxz, 0.d0) ! xz-plane reflection.
      call vector_ang(z, 90.d0, rz) ! rotation to convert xzR -> yzR
      call matrot(rz, hmirxz, hmiryz) ! rotate hmirxz into hmiryz
      n = 1
      iin = 0
      do 20 in = 1, n
         call take33mat(h0, 48, in, h)
         call matply(h0, hz180, dum)
         iin = iin + 1
         call add33mat(dum, 48, n + iin, h)
 20   continue
      n = n + iin

      iin = 0
      do 30 in = 1, n
         call take33mat(h0, 48, in, h)
         call matply(h0, hmirxz, dum)
         iin = iin + 1
         call add33mat(dum, 48, n + iin, h)
 30   continue
      n = n + iin

      iin = 0
      do 40 in = 1, n
         call take33mat(h0, 48, in, h)
         call matply(h0, hmiryz, dum)
         iin = iin + 1
         call add33mat(dum, 48, n+iin, h)
 40   continue
      n = n + iin
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     centro-symmetry
      if (icen) then
         call central(hc)       ! centro-sym op.
         iin = 0
         do 50 in = 1, n
            call take33mat(h0, 48, in, h)
            call matply(h0, hc, dum)
            iin = iin + 1
            call add33mat(dum, 48, n+iin, h)
 50      continue
         n = n + iin
      endif
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      return
      end subroutine orthorhombic_rot_sym
c----------------------------------------------------------------------c
      subroutine monoclinic_rot_sym(h,n,icen)
c     Returns the rotation symmetry operators: h.
c     h_[ij] transforms a vector to another albeit crystallographically
c     equivalent one in ca for monoclinic structure.
c     n returns the kinds of the transformation matrix h.
c     Symmetry operations:
c     1. 2-fold symmetry about [001]
c     2. Centro-symmetry (optional)
      implicit none
      real*8 h0(3,3), h(3,3,48), dum(3,3), hz180(3,3), hc(3,3), z(3)
      integer i,j,in,iin,n
      logical icen
cf2py intent(in) icen
cf2py intent(out) h,n
      write(*,*) 'Not tested yet for monoclinic_rot_sym'
      stop
      call imat(h0)
      do 10 i=1,3
      do 10 j=1,3
         h(i,j,1) = h0(i,j)
 10   continue
      call zmat(h0)

c      write(*,*) 'not completed!'
c      stop

      call zvec(z)
      z(3) = 1.d0
      call vector_ang(z, 180.d0, hz180)

      n = 1
      iin = 0
      do 20 in=1,n
         call take33mat(h0, 48, in, h)
         call matply(h0, hz180, dum)
         iin = iin + 1
         call add33mat(dum, 48, n+iin, h)
 20   continue
      n = n + iin

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      if (icen) then
         call central(hc)       ! centro-sym op.
         iin = 0
         do 30 in = 1, n
            call take33mat(h0, 48, in, h)
            call matply(h0, hc, dum)
            iin = iin + 1
            call add33mat(dum, 48, n+iin, h)
 30      continue
         n = n + iin
      endif
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      return
      end subroutine monoclinic_rot_sym
c----------------------------------------------------------------------c
      subroutine triclinic_rot_sym(h,n,icen)
c     Returns the rotation symmetry operators: h.
c     h_[ij] transforms a vector to another albeit crystallographically
c     equivalent one in ca for triclinic structure.
c     n returns the kinds of the transformation matrix h.
      implicit none
      real*8 h(3,3,48), hc(3,3), dum(3,3), h0(3,3)
      integer i,j,in,iin,n
      logical icen
cf2py intent(in) icen
cf2py intent(out) h,n
      write(*,*) 'Not tested yet for triclinic_rot_sym'
      stop
      call imat(h0)
      do 10 i=1,3
      do 10 j=1,3
         h(i,j,1) = h0(i,j)
 10   continue
      call zmat(h0)
      n = 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      if (icen) then
         call central(hc)       ! centro-sym op.
         iin = 0
         do in = 1, n
            call take33mat(h0, 48, in, h)
            call matply(h0, hc, dum)
            iin = iin + 1
            call add33mat(dum, 48, n+iin, h)
         enddo
         n = n + iin
      endif
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      return
      end subroutine triclinic_rot_sym
c----------------------------------------------------------------------c
      subroutine hexagonal_rot_sym(h,n,icen)
c     Returns the rotation symmetry operators: h.
c     h_[ij] transforms a vector to another albeit crystallographically
c     equivalent one in ca for orthorhombic structure.
c     n returns the kinds of the transformation matrix h.
c     Symmetry operations:
c     1. Mirror plane at 30 degree with respect to x1
c     2. 6-fold symmetry about [001]
      implicit none
      real*8 h(3,3,48), z(3), h0(3,3), r(3,3,6), r0(3,3), dum(3,3),
     $     hc(3,3), x(3), rxz(3,3), rz(3,3), hmirr(3,3)
      integer i,j,in,iin,n
      logical icen
cf2py intent(in) icen
cf2py intent(out) h, n
      write(*,*) 'Not tested yet for hexagonal_rot_sym'
      stop
      call imat(h0)
      do 10 i=1,3
      do 10 j=1,3
         h(i,j,1) = h0(i,j)
 10   continue
      call zmat(h0)

c     x and z unit vectors
      call zvec(x)
      call zvec(z)
      x(1) = 1.d0
      z(3) = 1.d0

c      write(*,*) 'not completed!'
c      stop

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     hmirr reflects by the plane at 30 deg with respect to x1
      call reflect(rxz, 0.d0)
c     rotate rxz by 180./6. degree about z
      call vector_ang(z, -30.d0, rz)
      call matrot(rz, rxz, hmirr)
      n = 1
      iin = 0
      do 20 in=1,n
         call take33mat(h0, 48, in, h)
         call matply(h0, hmirr, dum)
         iin = iin + 1
         call add33mat(dum, 48, n + iin, h)
 20   continue
      n = n + iin
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     Add 6-fold symmetry
      do 30 i=1,6
         call vector_ang(z, 360.d0/6. * i, dum)
         call add33mat(dum, 6, i, r)
 30   continue
      iin = 0
      do 40 in=1,n
         call take33mat(h0, 48, in, h)
      do 40 i=1,6
         call take33mat(r0, 6, i, r)
         call matply(h0, r0, dum)
         iin = iin + 1
         call add33mat(dum, 48, n+iin, h)
 40   continue
      n = n + iin
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     centro-symmetry
      if (icen) then
         call central(hc)       ! centro-sym op.
         iin = 0
         do 50 in = 1, n
            call take33mat(h0, 48, in, h)
            call matply(h0, hc, dum)
            iin = iin + 1
            call add33mat(dum, 48, n+iin, h)
 50      continue
         n = n + iin
      endif
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      return
      end subroutine hexagonal_rot_sym
c----------------------------------------------------------------------c
      subroutine tetragonal_rot_sym(h,n,icen)
      implicit none
      real*8 h0(3,3), h(3,3,48), dum(3,3), x(3), z(3), hc(3,3),
     $     rxz(3,3), hmirr(3,3), rz(3,3), r(3,3,4), r0(3,3)
      integer i,j,in,iin,n
      logical icen
cf2py intent(in) icen
cf2py intent(out) h, n
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     identity tensor
      write(*,*) 'Not tested yet for tetragonal_rot_sym'
      stop
      call imat(h0)
      do 10 i=1,3
      do 10 j=1,3
            h(i,j,1) = h0(i,j)
 10   continue
      call zmat(h0)

c     x and z unit vectors
      call zvec(x)
      call zvec(z)
      x(1) = 1.d0
      z(3) = 1.d0

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     hmirr reflects by the plane at 45 deg with respect to x1
      call reflect(rxz, 0.d0)
c     rotate rxz by 180./6. degree about z
      call vector_ang(z, -45.d0, rz)
      call matrot(rz, rxz, hmirr)

      n = 1
      iin = 0
      do 20 in = 1,n
         call take33mat(h0, 48, in, h)
         call matply(h0, hmirr, dum)
         iin = iin + 1
         call add33mat(dum, 48, n + iin, h)
 20   continue
      n = n + iin

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     Add 4-fold symmetry
      do 30 i=1,4
         call vector_ang(z, 360.d0/4. * i, dum)
         call add33mat(dum, 4, i, r)
 30   continue

      iin = 0
      do 40 in=1,n
         call take33mat(h0, 48, in, h)
      do 40 i=1,4
         call take33mat(r0, 4, i, r)
         call matply(h0, r0, dum)
         iin = iin + 1
         call add33mat(dum, 48, n+iin, h)
 40   continue
      n = n + iin
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     centro-symmetry
      if (icen) then
         call central(hc)       ! centro-sym op.
         iin = 0
         do 50 in = 1, n
            call take33mat(h0, 48, in, h)
            call matply(h0, hc, dum)
            iin = iin + 1
            call add33mat(dum, 48, n+iin, h)
 50      continue
         n = n + iin
      endif
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      return
      end subroutine tetragonal_rot_sym
c----------------------------------------------------------------------c
      subroutine trigonal_rot_sym(h,n,icen)
      implicit none
      real*8 h0(3,3), h(3,3,48), dum(3,3), x(3), z(3), hc(3,3),
     $     rxz(3,3), hmirr(3,3), rz(3,3), r(3,3,3), r0(3,3)
      integer i,j,in,iin,n
      logical icen
cf2py intent(in) icen
cf2py intent(out) h, n
      write(*,*) 'Not tested yet for trigonal_rot_sym'
      stop
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     identity tensor
      call imat(h0)
      do 10 i=1,3
      do 10 j=1,3
            h(i,j,1) = h0(i,j)
 10   continue
      call zmat(h0)

c     x and z unit vectors
      call zvec(x)
      call zvec(z)
      x(1) = 1.d0
      z(3) = 1.d0

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     hmirr reflects by the plane at 45 deg with respect to x1
      call reflect(rxz, 0.d0)
c     rotate rxz by 180./6. degree about z
      call vector_ang(z, -60.d0, rz)
      call matrot(rz, rxz, hmirr)
      n = 1
      iin = 0
      do 20 in = 1,n
         call take33mat(h0, 48, in, h)
         call matply(h0, hmirr, dum)
         iin = iin + 1
         call add33mat(dum, 48, n + iin, h)
 20   continue
      n = n + iin

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     Add 3-fold symmetry
      do 30 i=1,3
         call vector_ang(z, 360.d0/4. * i, dum)
         call add33mat(dum, 3, i, r)
 30   continue

      iin = 0
      do 40 in=1,n
         call take33mat(h0, 48, in, h)
      do 40 i=1,3
         call take33mat(r0, 3, i, r)
         call matply(h0, r0, dum)
         iin = iin + 1
         call add33mat(dum, 48, n+iin, h)
 40   continue
      n = n + iin
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     centro-symmetry
      if (icen) then
         call central(hc)       ! centro-sym op.
         iin = 0
         do 50 in = 1, n
            call take33mat(h0, 48, in, h)
            call matply(h0, hc, dum)
            iin = iin + 1
            call add33mat(dum, 48, n+iin, h)
 50      continue
         n = n + iin
      endif
      return
      end subroutine trigonal_rot_sym
c----------------------------------------------------------------------c
      subroutine cubic_orthrotmat(h, nvar)
      implicit none
      real*8 h0(3,3,48), h(3,3,48), aux1(3), aux2(3), aux3(3), aux4(3),
     $     aux5(3), aux6(3), aux33(3,3)
      integer i, j, ivar, n0, n, nvar
      logical isrotmat, tag(48)
c     Returns only orthogonal rotation matrices
      call cubic_rot_sym(h0, n0, .true.) ! Full 48 rotation operators.
      call finduniquemat(h0, h, n0, n, tag, .false., 0)

      call zvec(aux1)
      call zvec(aux2)
      call zvec(aux3)
      call zvec(aux4)
      call zvec(aux5)

      aux1(1) = 1.d0
      aux2(2) = 1.d0
      aux3(3) = 1.d0

      nvar = 0
      do ivar = 1, n
         call zvec(aux4)
         call zvec(aux5)
         call zvec(aux6)

         do 100 i=1,3
         do 100 j=1,3
            aux4(i) = aux4(i) + h(i,j,ivar) * aux1(j)
            aux5(i) = aux5(i) + h(i,j,ivar) * aux2(j)
            aux6(i) = aux6(i) + h(i,j,ivar) * aux3(j)
 100     continue

         call vnorm(aux4, 3)
         call vnorm(aux5, 3)
         call vnorm(aux6, 3)

         do i=1,3
            aux33(i,1) = aux4(i)
            aux33(i,2) = aux5(i)
            aux33(i,3) = aux6(i)
         enddo

         write(*,*) ivar
         if (isrotmat(aux33,3)) then
            nvar = nvar + 1
            do 200 i=1,3
            do 200j=1,3
               h(i,j,nvar) = aux33(i,j)
 200        continue
         endif
      enddo
      return
      end subroutine cubic_orthrotmat
c----------------------------------------------------------------------c
      subroutine hexagon_orthrotmat(h,nvar)
      implicit none
      real*8 h0(3,3,48), h(3,3,48), aux1(3), aux2(3), aux3(3), aux4(3),
     $     aux5(3), aux6(3), aux33(3,3)
      integer i, j, ivar, n0, n, nvar
      logical isrotmat, tag(48)
      call hexagon_rot_sym(h0,n0)
      call finduniquemat(h0,h,n0,n,tag,.false.,0)
      call zvec(aux1)
      call zvec(aux2)
      call zvec(aux3)
      call zvec(aux4)
      call zvec(aux5)

      aux1(1) = 1.d0
      aux1(2) = 1.d0
      aux1(3) = 1.d0

      nvar = 0
      do ivar = 1, n
         call zvec(aux4)
         call zvec(aux5)
         call zvec(aux6)

         do 100 i=1,3
         do 100 j=1,3
            aux4(i) = aux4(i) + h(i,j,ivar) * aux1(j)
            aux5(i) = aux5(i) + h(i,j,ivar) * aux2(j)
            aux6(i) = aux6(i) + h(i,j,ivar) * aux3(j)
 100     continue

         call vnorm(aux4, 3)
         call vnorm(aux5, 3)
         call vnorm(aux6, 3)

         do i=1,3
            aux33(i,1) = aux4(i)
            aux33(i,2) = aux5(i)
            aux33(i,3) = aux6(i)
         enddo

         write(*,*) ivar
         if (isrotmat(aux33,3)) then
            nvar = nvar + 1
            do 200 i=1,3
            do 200j=1,3
               h(i,j,nvar) = aux33(i,j)
 200        continue
         endif
      enddo

      return
      end subroutine hexagon_orthrotmat
c----------------------------------------------------------------------c
      subroutine cubic_rot_sym(h, n, icen)
c     Returns the rotation symmetry operators: h.
c     h_[ij] transforms a vector to another albeit crystallographically
c     equivalent one in ca for cubic structure.
c     n returns the kinds of the transformation matrix h.
c     Cubic symmetry has 48 including the centro-symmetry operation.
c     Symmetry operations:
c     1. 3-fold symmetry about [111] direction  - subr. rot111_120_240
c     2. 2-fold symmetry by (110) plane         - subr. mirror110
c     3. 4-fold symmetry about [100] direction  - subr. rot100_90_180_270
c     4. Centro-symmetry by (0,0,0) point       - subr. central

c     Only half of them can turn into 'orthogonal' rotation matrix.
c     This is embodied in subroutine cubic_orthrotmat.

c     the rotation operator h is not a rotation tensor... Mirroring
c     by (110) plane gives certain ones whose determinant is -1.

      implicit none
      real*8 h0(3,3), dum(3,3), h(3,3,48), h120(3,3), h240(3,3),
     $     hmirror(3,3), hx90(3,3), hx180(3,3), hx270(3,3), hc(3,3)
      integer i, j, n, in, iin
      logical icen
cf2py intent(out) h, n
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     identity tensor for the first symmetry element for h.
      call imat(h0)
      do 10 i=1,3
      do 10 j=1,3
         h(i,j,1) = h0(i,j)
 10   continue
      call zmat(h0)
      n = 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     Rotates 120 and 240 about [111] direction. 3 fold symmetry
c      call rot111_60_120(h120, h240)
      call rot111_120_240(h120, h240)
      iin = 0
      do 100 in=1, n
         call take33mat(h0, 48, in, h)
         call matply(h0, h120, dum)
         iin = iin + 1
         call add33mat(dum, 48, n+iin, h)
         call matply(h0, h240, dum)
         iin = iin + 1
         call add33mat(dum, 48, n+iin, h)
 100  continue
      n = n + iin
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      call mirror110(hmirror)
      iin = 0
      do 200 in=1, n
         call take33mat(h0, 48, in, h)
         iin = iin + 1
         call matply(h0, hmirror, dum)
         call add33mat(dum, 48, n+iin, h)
 200  continue
      n = n + iin
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     rotations of 90, 180, 270 around x3
c      call rot100_90_180_270(h90, h180, h270)
      call rot100_90_180_270(hx90, hx180, hx270)
c     $     hy90, hy180, hy270, hz90, hz180, hz270)
      iin = 0
      do 300 in=1, n
         call take33mat(h0, 48, in, h)
         iin = iin + 1
         call matply(h0, hx90, dum)
         call add33mat(dum, 48, n+iin, h)
         iin = iin + 1
         call matply(h0, hx180, dum)
         call add33mat(dum, 48, n+iin, h)
         iin = iin + 1
         call matply(h0, hx270, dum)
         call add33mat(dum, 48, n+iin, h)
 300  continue
      n = n + iin
c----------------------------------------------------------------------c
c     Central symmetry by the origin (0,0,0)
      if (icen) then
         call central(hc)
         iin = 0
         do 400 in=1,n
            call take33mat(h0, 48, in, h)
            iin = iin + 1
            call matply(h0, hc, dum)
            call add33mat(dum, 48, n+iin, h)
 400     continue
         n = n + iin
      endif
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c      write(*,'(a, i4)') 'n:', n
      return
      end subroutine cubic_rot_sym
c----------------------------------------------------------------------c
      subroutine hexagon_rot_sym(h,n)
c     Returns the rotation symmetry operators: h
c     h_[ij] transforms a vector to another albeit rystallographically
c     equivalent one in ca for hexagonal structure.
c     n returns the kinds of the transformation matrix h.

c     - mirror plane at 30 deg with respect to x1
c     - rotation of 2*pi/6 around axis <001>
      implicit none
      real*8 h0(3,3), h(3,3,24), ang, pi, hrot(3,3,5)
      integer i,j,k,n,nr,mn,ns
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      pi = dacos(0.d0) * 2.d0
c     identity tensor for the first symmetry element for h.
      call imat(h0)
      h(:,:,1) = h0(:,:)
      n = 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     mirror plane at 30 deg with respect to x1
      ang = pi / 6.0
      h(1,1,2) = dcos(ang)**2 - dsin(ang)**2
      h(2,2,2) = - h(1,1,2)
      h(3,3,2) = 1.d0
      h(1,2,2) = 2.d0 * dcos(ang) * dsin(ang)
      h(2,1,2) = h(1,2,2)
      n = n + 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     rotation of 2*pi/6 around axis <001>
      do nr=1,5
         ang = nr * 2. * pi / 6.
         hrot(1,1,nr) = dcos(ang)
         hrot(2,2,nr) = dcos(ang)
         hrot(3,3,nr) = 1.d0
         hrot(1,2,nr) =-dsin(ang)
         hrot(2,1,nr) = dsin(ang)
      enddo

      do 20 nr=1,5
      do 20 ns=1,n
         mn = nr * n + ns
      do 20 i=1,3
      do 20 j=1,3
      do 20 k=1,3
         h(i,j,mn) = h(i,j,mn) + hrot(i,k,nr) * h(k,j,ns)
 20   continue
      n = mn

      return
      end subroutine
c----------------------------------------------------------------------c
      subroutine sxfreader(fsx, icrysym, cdim, angs)
c     a replacement for crystal symmetry on ioption 1 in which crystal
c     symmetry is read from the given file and returns an integer label.
      implicit none
      integer icrysym, i
      character*80 fsx
      character*5 crysym

      real*8 cdim(3), angs(3)

      open(334, file=fsx, status='old')
      read(334,'(a)')
      read(334,'(a)') crysym
      read(334,  *) (cdim(i), i=1,3), (angs(i), i=1,3)
      close(334)

      icrysym=0
      if(crysym.eq.'cubic' .or. crysym.eq.'CUBIC') icrysym=1
      if(crysym.eq.'hexag' .or. crysym.eq.'HEXAG') icrysym=2
      if(crysym.eq.'trigo' .or. crysym.eq.'TRIGO') icrysym=3
      if(crysym.eq.'tetra' .or. crysym.eq.'TETRA') icrysym=4
      if(crysym.eq.'ortho' .or. crysym.eq.'ORTHO') icrysym=5
      if(crysym.eq.'monoc' .or. crysym.eq.'MONOC') icrysym=6
      if(crysym.eq.'tricl' .or. crysym.eq.'TRICL') icrysym=7
      if(icrysym.eq.0) then
         write(*,*) ' *** cannot recognize the crystal symmetry'
         stop
      endif
      return
      end subroutine sxfreader
c----------------------------------------------------------------------c
c     Calculates vectors in unit cell.
      subroutine unitcv(miller, vector, fsx, verbose)
c     miller : Miller indices
c     vector : to-be-returned vector in cartesian coordinate system
c     icrysym: crystal symmetry following the convention of CNT
c     angs   : angles of the unit cells.
c     cdim   : unit cell's dimension
c     verbose: logical flag to decide if more information to print out.
      implicit none
      integer icrysym, i, j
      logical verbose
      character*80 fsx
      real*8 miller(3), vector(3),angs(3), cdim(3), r(3,3), pi,
     $     a, b, c, al, be, ga, cvol
cf2py intent(in) miller, icrysym, angs, cdim, verbose
cf2py intent(out) vector
      pi = dacos(0.d0) * 2.d0

      call sxfreader(fsx, icrysym, cdim, angs)
      if (icrysym.ne.1) then
         write(*,*) 'Subr. unitcv is developed for cubic crystal.'
         stop
      endif

c     Follows the subroutine crystal_symmetry by CNT in VPSC.
      do i=1,3
         angs(i) = angs(i) * pi / 180.d0
      enddo
      al = angs(1)
      be = angs(2)
      ga = angs(3)

      a = cdim(1)
      b = cdim(2)
      c = cdim(3)

c     The relationship between the cartensian co-ordinate x,y,z and
c     the unit cell axes (a,b,c) is necessary.

c     The crystallographic unit cell coordinate system is a fractional
c     coordinate system which is characterized by the lengths of
c     each axes, a,b, and c, and the angles between them, alpha, beta
c     and gamma.


c     a case of it is as below when
c     from http://wwww.angelfire.com/linux/myp/FracCor/fraccor.html

c     a is taken along x-axis and b lines in x-y plane.


c     The point p in the cartensian cooridnate
c          |x|
c      p = |y|
c          |z|
c     and p` in the fractional coordinate
c          |x`|
c      p`= |y`|
c          |z`|
c
c     Fractional coordinates are weights, linear combination of which
c     with vectors a, b, c reproduces point p.
c     x` (a) + y` (b) + z` (c) = p
c     -> such transformation relationship:
c         p = Mp`  and M = [(a) (b) (c)]

c     To carry out such transformation we have to multiply p` by M
c     It is clear that       |a|            |b cos(gamma)|
c                      (a) = |0|  and (b) = |b sin(gamma)|
c                            |0|            |      0     |

c     And vector (c) is missing. Let's find each component of c
c     that is referred in the cartesian coordinate system.

c     First two coordinates of vector c are easily found using the dot
c     product:
c       (a) . (c) = |a||c|cos(beta) =  a1c1              (eq.1)
c       (b) . (c) = |b||c|cos(alpha) = b1c1 + b2c2       (eq.2)

c     from eq.1) c1 = |c| cos(beta)
c     from eq.2) c2 = {|b||c| cos(alpha) - b1c1} / b2

c     The last component of c is found using the cross product:
c     bxc=|b||c|sin(alpha)*(n), where n is normal to the plane bc (eq.3)
c     |bxc| = sqrt{b2^2 c3^2 + b1^2 c3^2 + (b1c2-b2c1)^2}

c     (bxc)_1 = b2c3 - b3c2 = b2c3
c     (bxc)_2 = b3c1 - b1c3 = - b1c3
c     (bxc)_3 = b1c2 - b2c1

c  |bxc| = sqrt(b2^2c3^2 + b1^2c3^2 + b1^2c2^2 + b2^2c1^2 - 2*b1b2c1c2)
c        = sqrt{b2^2c3^2 + b1^2c3^2 + (b1c2 - b2c1)^2)}
c        = sqrt{(b2^2+b1^2) c3^2 + (... )^2}

c     -> Rearrange it explicitly on c3
c     c3^2 = {|bxc|^2 - (...)^2} / (b2^2+b1^2)
c     c3   = (+-)sqrt[{|bxc|^2 - (...)^2} / (b2^2+b1^2)]
c     from (eq.3)
c      = (+-)sqrt[{(|b||c|sin(alpha))^2 -( b1c2-b2c1)^2}/(b1^2+b2^2)]

c     solving the above gives:
c     c3 = c/sin(gamma) {1-cos(alpha)^2 - cos(beta)^2 +
c                        2bc cos(alpha) cos(beta) cos(gamma)}^0.5
c        = V/{ab sin(gamma)}

c     alpha: al, beta: be, gamma: ga
c     where v is the volume of the unit parallelepiped
c     v = sqrt{1-cos(al)^2-cos(be)^2-cos(ga)+2cos(al)cos(be)cos(ga)}

c |a| |a bcos(gamma)                 ccos(beta)                  | |x|
c |b|=|0 bsin(gamma) c{cos(alpha)-cos(beta)cos(gamma)}/sin(gamma)| |y|
c |c| |0      0                    c{v/sin(gamma)}               | |z|

c     volume of the parallelepiped. - - - - - - - - - - - - - - - - - -c

      cvol = sqrt(1. - dcos(al)**2 - dcos(be)**2 - dcos(ga) +
     $     2 * dcos(al) * dcos(be) * dcos(ga))
      if (verbose) write(*,'(a, f7.3)') 'c vol: ', cvol

      r(1,1) = a
      r(1,2) = b * dcos(ga)
      r(1,3) = c * dcos(be)
      r(2,1) = 0.d0
      r(2,2) = b * dsin(ga)
      r(2,3) = c * (dcos(al) - dcos(be) * dcos(ga)) / dsin(ga)
      r(3,1) = 0.d0
      r(3,2) = 0.d0
      r(3,3) = c * (cvol / dsin(ga))

      if (verbose) then
         write(*,*) 'The transformation matrix:'
         do i=1,3
            write(*, '(3f7.3)') (r(i,j),j=1,3)
         enddo
         write(*,*)
      endif

c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -c
      if (icrysym.eq.1) then
c     (a,b,c)_unitcel -> (x,y,c)_cartesian
c     Remember a is taken along x
c     and b lies in the xy plane.
         do 1000 i=1,3
            vector(i) = 0.d0
         do 1000 j=1,3
            vector(i) = vector(i) +  r(i,j) * miller(i)
 1000    continue
      else
         write(*,*) 'not implemented yet!'
         stop
      endif

      call vnorm(vector, 3)

      if (verbose) then
         write(*,*) 'Transformed vector in the unit cell'
         write(*,'(3f7.3)') vector
      endif

      end subroutine unitcv
c----------------------------------------------------------------------c
      subroutine cubic_eqvect(v, vuniq, n, icen, ipr)
c     Given one vector v, finds only unique vectors considering central
c     symmetry flag (icen). Argument 'n' is supposed to be updated based
c     the number of unique sets of vectors.

cf2py intent(in) v
cf2py inetnt(out) vuniq, n
      implicit none
      real*8 h(3,3,48), v(3), vall(3,48), vuniq(3,48)
      integer i, j, n, ii, nuniq
      logical icen, ipr
      if (ipr) write(*,*) 'icen:', icen
      call cubic_rot_sym(h, n, .false.) ! n is now 24
c      open(660, file='cubic_h.out', status='unknown')
      do ii=1,n
         do i=1,3
            vall(i,ii) = 0.
c            write(660,'(3f7.3)') (h(i,j,ii), j=1,3)
            do j=1,3
               vall(i,ii) = vall(i,ii) + h(i,j,ii) * v(j)
            enddo
         enddo
         if (ipr) then
            write(*,'(i3,3f6.2)') ii, (vall(i,ii), i=1,3)
c            write(660,*)
         endif
      enddo

c      read(*,*)
      if (n.eq.0) then
         write(*,*) 'n is wrong! at 419'
         stop
      endif
      call finduniquevect(vall, vuniq, n, nuniq, icen, ipr)
      if (ipr) then
         write(*,*) 'n: ', n, ' in subroutine cubic_eqvect'
         write(*,*) 'nuniq: ', nuniq, ' in subroutine cubic_eqvect'
      endif
      n = nuniq
      if (ipr) then
         write(*,*) 'Unique vector sets are found:'
         do ii=1, nuniq
            write(*, '(3f7.3)') (vuniq(i, ii), i=1,3)
         enddo
c         read(*,*)
      endif
      return
      end subroutine cubic_eqvect
c----------------------------------------------------------------------c
      subroutine cubic_uniq_mat
      implicit none
      real*8 hh(3,3,48), b(3,3,48)
      integer n,  n1
      logical tag(48)
c      logical isin
      call cubic_rot_sym(hh,n,.true.)
      call finduniquemat(hh, b, 48, n1, tag, .false., 1)
      end subroutine cubic_uniq_mat

c     Crystal symmetry operators
c----------------------------------------------------------------------c
c     Returns the 3-fold symmetry operators about [111] direction
      subroutine rot111_120_240(h120, h240)
      implicit none
      real*8 h240(3,3), h120(3,3), u(3)
cf2py intent(out) h120, h240
c     Rotation of (pi/3) & (2*pi/3) about <111> direction
      u(1) = 1.d0
      u(2) = 1.d0
      u(3) = 1.d0
      call vector_ang(u, 120.d0, h120)
      call vector_ang(u, 240.d0, h240)
      return
      end subroutine rot111_120_240
c----------------------------------------------------------------------c
c     returns mirror rotation about (110) plane
      subroutine mirror110(hmirror)
      implicit none
      real*8 hmirror(3,3), z(3), rxz(3,3), rz(3,3)
      integer i, j, k, l
cf2py intent(out) hmirror
      z(1) = 0.
      z(2) = 0.
      z(3) = 1.
      call reflect(rxz, 0.d0)        ! xz-plane reflection. v` <- [rxz] v
      call vector_ang(z, -45.d0, rz) ! rotating about z by -45 deg.
c     rotates reflection angle to (110)
c     reflecting the vector by the rotated reflection matrix.
c     new matrix= rz^t . rxz . rz = rz_ki rxz_ik rz_kj
      do 100 i=1,3
      do 100 j=1,3
         hmirror(i,j) = 0.
      do 100 k=1,3
      do 100 l=1,3
c        rz is the rotation mat.
         hmirror(i,j) = hmirror(i,j) + rz(k,i) * rxz(k,l) * rz(l,j)
 100  continue
      return
      end subroutine mirror110
c----------------------------------------------------------------------c
c     returns rotation matrix that rotates 90, 180, 270 degrees around
c     [100] directions: 4-fold symmetry.
      subroutine rot100_90_180_270(rx90, rx180, rx270)
      implicit none
      real*8 x(3), rx90(3,3), rx180(3,3), rx270(3,3)
cf2py intent(out) rx90, rx180, rx270
c     rotations of 90, 180, 270 around x3
      x(1) = 0.d0
      x(2) = 0.d0
      x(3) = 1.d0
      call vector_ang(x,  90.d0, rx90 )
      call vector_ang(x, 180.d0, rx180)
      call vector_ang(x, 270.d0, rx270)
      return
      end subroutine rot100_90_180_270
c----------------------------------------------------------------------c
c     returns the matrix that rotate a vector by 180 degree about the
c     origin. in other words, a rotation matrix that flips a vector
c     to the opporiste direction. that is ... as below:
c     | -1  0  0 |
c     |  0 -1  0 |
c     |  0  0 -1 |
      subroutine central(r)
      implicit none
      real*8 r(3,3)
      integer i
      call zmat(r)
      do 100 i=1,3
         r(i,i) = -1.d0
 100  continue
      end subroutine central
c----------------------------------------------------------------------c
c     calculate the improper mirror rotation matrix
      subroutine reflect(r, th)
      implicit none
      real*8 r(3,3), th, pi, c2t, s2t
cf2py intent(in) th
cf2py intent(out) r
      pi = dacos(0.d0) * 2.d0
      c2t = dcos(2.d0 * th * pi / 180.d0)
      s2t = dsin(2.d0 * th * pi / 180.d0)
c     reflection matrix: (reflection by xz plane with t=0 or pi.)
c                        (reflection by yz plane with t=2/pi)
c     | cos(2t)  sin(2t)  0 |
c     | sin(2t) -cos(2t)  0 |
c     |    0        0     1 |
      call zmat(r)
      r(1,1) =  c2t
      r(1,2) =  s2t
      r(2,1) =  s2t
      r(2,2) =  c2t * (-1.d0)
      r(3,3) =  1.d0
      return
      end subroutine reflect
c----------------------------------------------------------------------c
c     returns the rotation matrix, r, that rotates the coordinate system
c     of others by th (in degree) about the vector u.
c
c     Realized that I have to indicated in which sense it rotates about
c     the given axis. So this note is added:
c     if the associated 3D space is 'right-handed', this rotation will be
c     counter-clockwise for an observer placed so that the axis u goes
c     in her/his direction.
      subroutine vector_ang(u, th, r)
      implicit none
      real*8 u(3), du(3), th, r(3,3), ct, st, cm(3,3), idx(3,3)
      real*8 pi
      integer i, j
cf2py intent(in) u, th
cf2py intent(out) r
      pi = dacos(0.d0) * 2.d0
      do i=1,3
         du(i) = u(i)
      enddo
      call vnorm(du, 3)
      call imat(idx)
      ct = dcos(th *  pi / 180.d0)
      st = dsin(th *  pi / 180.d0)
c     R = I cos(th) + sin(th) [u]_x + (1 - cos(th)) [f; u](u; f*)
c     r_ij = I_ij ct + st cm_ij + (1-ct) u_i u_j
      call crossop(du, cm)
      do 100 i=1,3
      do 100 j=1,3
         r(i,j) = idx(i,j) * ct + st * cm(i,j) +
     $        (1 - ct) * du(i) * du(j)
 100  continue
      return
      end subroutine vector_ang
c----------------------------------------------------------------------c
      logical function iseqrot(a,b)
c     check if the given rotation matrices are equivalent              c
      implicit none
      real*8 a(3,3), b(3,3), x0(3), x1(3), x0m(3)
      integer i, j, k, l, ncount, foundi(3)
      logical isvecteq, f(3)
c              |             |
c     matrix a | (x) (y) (z) |
c              |             |
c     matrix a rotates [100] into [x]
c                      [010] into [y]
c                      [001] into [z]
c     Another mat. consisting of any permutation of +-[x], +-[y], +-[z]
c     will be the same (?) -> NO
      do i=1,3
         f(i) = .false.
      enddo
      iseqrot = .false.
      do 90 i=1,3  ! for each of columns of mat a.
c        x0 is the current column
         do j=1,3
            x0(j) = a(i,j)
            x0m(j) = a(i,j) * (-1)
         enddo
c        write(*,*) 'x0, x0m'      'x(1,2,3)'
c        write(*,*)x0, x0m     or '-x(1,2,3)'
c        x1 is the current column of matrix b
         ncount = 0
         do k=1,3
            do l=1,3
               x1(l) = b(k,l)
            enddo
c            write(*,*) 'x1'
c            write(*,*) x1
            if (isvecteq(x0, x1).or.isvecteq(x0m, x1)) then
c               write(*,*) 'Equal!'
               foundi(i) = k
               ncount = ncount + 1
            endif
         enddo
         if (ncount.ne.1) then
            iseqrot = .false.
            goto 100
         else
            f(i) = .true.
         endif
 90   continue

      if (f(1) .and. f(2) .and. f(3)) then
         iseqrot = .true.
c         do i=1,3
c            write(*,'(3f5.2, 2x, 3f5.2)') (a(i,j), j=1,3),
c     $           (b(i,j), j=1,3)
c         enddo
c         write(*,'(3i2)') foundi
      endif
 100  return
      end function iseqrot
c----------------------------------------------------------------------c
      include 'mat_lib.f'
