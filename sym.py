"""
General crystallographic symmetry operators
Cubic and hexagonal operators are available.


-- List of symmery operators
def cubic()
def hexag()
def tetra()
def triclinic()
"""

import numpy as np
import time
## symmetry operators
def __60_120_rot111__(h):
    """
    For the given h operation,
    rotations of (pi/3) & (2*pi/3) around <111>
    are performed and returned

    *cubic
    """
    hx = h.copy()
    h60 = np.zeros((3,3)); h120 = np.zeros((3,3))
    h60[0,2] = 1.
    h60[1,0] = 1.
    h60[2,1] = 1.
    
    h120[0,1] = 1.
    h120[1,2] = 1.
    h120[2,0] = 1.
    return np.dot(h60,hx), np.dot(h120,hx)

def __mirror_110__(h):
    """
    Given the operation h, mirrored across the (110) plane returned

    *cubic
    """
    hx = h.copy()
    hm = np.zeros((3,3))
    hm[0,1] = 1.
    hm[1,0] = 1.
    hm[2,2] = 1.
    return np.dot(hm, hx)

def __rot_90_180_270__(h):
    """
    Given the operation h,
    the three rotated operations are returned

    *cubic
    """
    cos = np.cos; sin = np.sin; pi = np.pi
    hx = np.zeros((3,3,3))
    h_ = h.copy(); htemp = []
        
    for m in range(3):
        ang = pi/2. * float(m+1)
        hx[m,0,0] = cos(ang)
        hx[m,1,1] = cos(ang)
        hx[m,2,2] = 1.0
        hx[m,0,1] = -sin(ang)
        hx[m,1,0] = sin(ang)
        hx[m,0,2] = 0.
        hx[m,2,0] = 0.
        hx[m,1,2] = 0.
        hx[m,2,1] = 0.
        pass
    for m in range(3):
        htemp.append( np.dot(hx[m], h_) )
        pass
    return np.array(htemp)

def __rot_nrot_x1__(h,nrot):
    """
    Mirror plane at 30 or 60 or 45 deg with respect to x1
    
    *hexagonal, trigonal, tetragonal

    hexa: nrot = 6
    trig: nrot = 3
    tetr: nrot = 4
    """
    cos = np.cos; sin = np.sin; pi=np.pi
    hx = np.zeros((3,3))
    ang = pi/float(nrot)
    hx[0,0] = cos(ang)**2 - sin(ang)**2
    hx[1,1] = -h[0,0]
    hx[2,2] = 1.0
    hx[0,1] = 2.*cos(ang)*sin(ang)
    hx[1,0] = h[0,1]
    return np.dot(hx,h)

def __rot_nrot_001__(h, csym=None):
    """
    Rotations of 2*pi/nrot around axis <001>
    
    *hexagoanl, trigonal, tetragonal

    ---------
    Arguments
    h: symmetry operators
    csym: 'hexa'
    """
    if   csym=='hexag': nrot=6
    elif csym=='trigo': nrot=3
    elif csym=='tetra': nrot=4
    else: print 'Unexpected Error'; raise IOError

    cos = np.cos; sin = np.sin; pi = np.pi
    hx = np.zeros((nrot-1,3,3))
    h_ = h.copy(); htemp = []

    for nr in range(nrot-1):
        ang = (nr+1)*2.*pi/nrot
        hx[nr,0,0] = cos(ang)
        hx[nr,1,1] = cos(ang)
        hx[nr,2,2] = 1.0
        hx[nr,0,1] =-sin(ang)
        hx[nr,1,0] = sin(ang)
        
    for nr in range(nrot-1):
        htemp.append(np.dot(hx[nr], h_))

    return np.array(htemp)

def __trim0__(h):
    """
    if a value in the matrix h is fairly close to +-0.
    then returns zero. In that way, trimming is performed
    on the every component of given h matrix.
    """
    hx = h.copy()
    for i in range(len(hx)):
        for j in range(len(hx[i])):
            if abs(hx[i,j]) < 0.1**6:
                hx[i,j] = 0.
    return hx


""" -- mmm sample symmetry is found in COD_conv.py """
def __mmm__():
    m0 = [[ 1, 0, 0], [0, 1, 0], [0, 0, 1]]
    m1 = [[ 1, 0, 0], [0,-1, 0], [0, 0,-1]]
    m2 = [[-1, 0, 0], [0, 1, 0], [0, 0,-1]]
    m3 = [[-1, 0, 0], [0,-1, 0], [0, 0, 1]]
    h = np.array([m0,m1,m2,m3])
    return h

### deprecated ###
# def __ortho__(v): 
#     """
#     Orthogonal sample symmetry to a vector (v) in 3D
#     """
#     v1 = np.aray([a[0], a[1], a[2]])
#     v2 = np.array([-a[0], a[1], a[2]])
#     v3 = np.array([a[0], -a[1], a[2]])
#     v4 = np.array([-a[0], -a[1], a[2]])
    
#     v5 = v1.copy()* -1
#     v6 = v2.copy()* -1
#     v7 = v3.copy()* -1
#     v8 = v4.copy()* -1

#     return v1, v2, v3, v4


## cubic symmetry
def cubic():
    H = []     # H is the master list containing all the numpy arrays of operations
    H.append(np.identity(3))    # identity operation
    
    # rotations of (pi/3) & (2*pi/3) around <111>
    niter = len(H)
    for i in range(niter):
        h60, h120 = __60_120_rot111__(h=H[i].copy())
        h0 = h60.copy(); h1 = h120.copy()
        H.append(h0)
        H.append(h1)
        
    # mirror across the plane (110)
    niter = len(H)
    for i in range(niter):
        h = __mirror_110__(h=H[i].copy())
        H.append(h)

    # rotations of 90, 180, 270 around x3
    niter = len(H)
    for i in range(niter):
        h1, h2, h3 = __rot_90_180_270__(h=H[i].copy())
        h90 = h1.copy(); h180 = h2.copy(); h270 = h3.copy()
        H.append(h90)
        H.append(h180)
        H.append(h270)

    H = np.array(H) # Make the H as numpy array

    # trim the values.
    for i in range(len(H)):
        H[i] = __trim0__(h=H[i])
    return H    

def cubic_centro():
    h_old = cubic()
    h_new = []
    h_n = [[-1,0,0],[0,-1,0],[0,0,-1]]
    for i in range(len(h_old)):
        h_new.append(np.dot(h_old[i],h_n))
    return h_new
        

def triclinic():
    H = []
    H.append(np.identity(3))
    return H

## hexagonal
def hexag():
    H = []
    H.append(np.identity(3))

    #mirror plane at 30 degree with respect to x1
    nrot = 6
    niter = len(H)
    for i in range(niter):
        h = __rot_nrot_x1__(h=H[i].copy(),nrot=nrot)
        H.append(h)

    #rotations of 2*pi/6 around axis <001> for hexagonals.
    niter = len(H)
    for i in range(niter):
        h = __rot_nrot_001__(h=H[i],csym='hexag')
        for ix in range(len(h)):
            H.append(h[ix])

    for i in range(len(H)):
        H[i] = __trim0__(h=H[i])
    return H


## trigonal
def trigo():
    H = []
    H.append(np.identity(3))
    #mirror plane 60 degree with respect to x1
    nrot = 3
    niter = len(H)
    for i in range(niter):
        h = __rot_nrot_x1__(h=H[i].copy(), nrot=3)
        H.append(h)
    #rotations of 2*pi/3 around axis <001> for trigonals
    niter = len(H)
    for i in range(niter):
        h = __rot_nrot_001__(h=H[i], csym='trigo')
        H.append(h)
    
    for i in range(len(H)):
        H[i] = __trim0__(h=H[i])
    return H

def tetra():
    H = []
    H.append(np.identity(3))

    #mirror plane at 45 degree with respect to x1
    nrot = 4
    niter = len(H)
    for i in range(niter):
        h = __rot_nrot_x1__(h=H[i].copy(), nrot=nrot)
        H.append(h)

    #rotations of 2*pi/4 around axis <001> for hexagonals.
    niter = len(H)
    for i in range(niter):
        h = __rot_nrot_001__(h=H[i], csym='tetra')
        for ix in range(len(h)):
            H.append(h[ix])

    for i in range(len(H)):
        H[i] = __trim0__(h=H[i])
    return H

## 


def cvec(cdim=None, cang=None):
    """
    Generates and returns 'cvec[i,n]' of the unit cell that
    is characterized by unit cell dimension together with
    axes' angles

    ---------
    Arguments
      cdim=[1.,1.,1.]
      cang=[90.,90.,90.] : Should be in angle [90.,90.,90.] not radian
    """

    cdim = np.array(cdim)
    cang = np.array(cang)
    # angle to radian    
    cang = cang * np.pi/180.
    # cvec
    cvec = np.zeros((3,3))
    
    cvec[0,0] = np.sin(cang[1])
    cvec[1,0] = 0.
    cvec[2,0] = np.cos(cang[1])
    
    cvec[0,1] = (np.cos(cang[2])-np.cos(cang[0])\
                     *np.cos(cang[1]))/np.sin(cang[1])
    cvec[2,1] = np.cos(cang[0])
    cvec[1,1] = np.sqrt(1.-cvec[0,1]**2-cvec[2,1]**2)
    
    cvec[0,2] = 0.
    cvec[1,2] = 0.
    cvec[2,2] = 1.

    for i in range(3):
        for j in range(3):
            cvec[i,j] = cdim[j] * cvec[i,j]
    return cvec


def cv(pole, cdim=None, cang=None, csym=None):
    """
    Creats vector of the pole taking care of its unit cell's
    dimension and axes' angles.
    """
    sqrt = np.sqrt
    cvect = cvec(cdim=cdim, cang=cang)

    s = np.zeros((3,))
    s[2] = ( pole[2]                                     ) / cvect[2,2]
    s[0] = ( pole[0] - cvect[2,0]*s[2]                   ) / cvect[0,0]
    s[1] = ( pole[1] - cvect[0,1]*s[0] - cvect[2,1]*s[2] ) / cvect[1,1]

    norm = sqrt(s[0]**2 + s[1]**2 + s[2]**2)
    for i in range(3):
        s[i] = s[i]/norm
        if abs(s[i])<0.1**5: s[i]=0.
    return s
