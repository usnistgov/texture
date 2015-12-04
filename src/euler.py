import numpy as np
import math

def in_plane_rot(th):
    """
    Return an in-plane rotation matrix
    | cos(th)   -sin(th)   0 |
    | sin(th)    cos(th)   0 |
    | 0          0         1 |

    Argument
    --------
    th  (rotation angle in degree)
    """
    th = th * math.pi / 180.

    sth = math.sin(th)
    cth = math.cos(th)
    rot = np.zeros((3,3))
    rot[0,0] =  cth
    rot[0,1] = -sth
    rot[1,0] =  sth
    rot[1,1] =  cth
    rot[2,2] =  1.
    return rot

def euler(ph=None, th=None, tm=None, a=None, echo=True):
    """
    note:
          This is a pythonized fortran subroutine embedded in the VPSC7.sub
          done by youngung 2010-12-23

          Matrix A(i,j) is returned as a numpy array
          all Euler angles are in degree not in radian

          # A[i,j] transforms from SA to CA

          Thanks to python's non-hardwirable-arugment feature,
          if a matrix is not given, it automatically calculates
          the matrix with given angles.
          Vice versa, if matrix is given, given angle aurgments
          are ignored and new euler angles are returned.

    Nomenclature of Euler angle follows Bunge's.
          ph = phi1,
          th = phi
          tm = phi2
    """
    if type(a).__name__=='NoneType':  a=np.resize(np.array(()),(3,3));iopt=2
    else:
        if type(a).__name__=='ndarray':
            iopt = 1
            pass
        else:
            print 'Error: Unexpected Matrix a type'
            print 'It should be numpy.ndarry!'
            raise IOError

    if iopt==1 :
        th = math.acos(a[2,2])  #Radian
        if abs(a[2,2] > 0.99999):
            tm = 0.
            ph = math.atan2(a[0,1],a[0,0]) #Radian
        else:
            sth = math.sin(th)
            tm = math.atan2(a[0,2]/sth,a[1,2]/sth)
            ph = math.atan2(a[2,0]/sth,-a[2,1]/sth)
        th = th * 180./math.pi
        ph = ph * 180./math.pi
        tm = tm * 180./math.pi
        return [ph,th,tm] #phi1, phi, phi2

    elif (iopt == 2):
        angles = [ph,th,tm]
        if any(angles[i] == None for i in range(len(angles))):
            print 'Angles must be give if iopt==2'
            raise IOError

        """ Convert the angle into Radian"""
        ph = ph * math.pi / 180.
        th = th * math.pi / 180.
        tm = tm * math.pi / 180.

        sph = math.sin(ph)
        cph = math.cos(ph)
        sth = math.sin(th)
        cth = math.cos(th)
        stm = math.sin(tm)
        ctm = math.cos(tm)

        a[0,0] =  ctm * cph - sph * stm * cth
        a[1,0] = -stm * cph - sph * ctm * cth
        a[2,0] =  sph * sth
        a[0,1] =  ctm * sph + cph * stm * cth
        a[1,1] = -sph * stm + cph * ctm * cth
        a[2,1] = -sth * cph
        a[0,2] =  sth * stm
        a[1,2] =  ctm * sth
        a[2,2] =  cth

        if echo==True:
            print 'Matrix a is '
            print a
        return a

    else: print 'Error: Improper iopt input'; raise IOError
