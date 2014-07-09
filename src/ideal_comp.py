"""
rewriting the module text.
"""
from euler import euler as eul
import random, os, upf
import numpy as np
gauss = random.gauss

def crossop(u):
    """
    Cross product operator
    ** implementation of subroutine crossop in mat_lib.f
    """
    m = np.zeros((3,3))
    m[0,1] = -u[2]
    m[0,2] =  u[1]
    m[1,0] =  u[2]
    m[1,2] = -u[0]
    m[2,0] = -u[1]
    m[2,1] =  u[0]
    return m

def vector_ang(u, th):
    """
    ** Implementation of vector_ang in cs.f to python

    arguemnts:
      u : vector
      th: angle in degree

   returns:
      r : Rotation matrix that rotates the coordinate system
        of others by th (in degree) about the vector u.
        --
        If the associated 3D space is 'right-handed', this rotation will be
        counter-clockwise for an observer placed so that the axis u goes
        in her/his direction.
    """

    idx = [[1,0,0],[0,1,0],[0,0,1]]
    r = np.zeros((3,3))

    pi = np.pi
    ct = np.cos(th * pi / 180.)
    st = np.sin(th * pi / 180.)

    cm = crossop(u)
    for i in range(3):
        for j in range(3):
         r[i][j] = idx[i][j] * ct + st * cm[i][j] +\
             (1 - ct) * u[i] * u[j]
    return r

def polar2vect(delta, phi):
    """ Given delta and phi returns a vector"""
    import numpy as np
    x = np.cos(phi)
    y = np.sin(phi)
    z = np.sin(delta)
    vector = np.array([x,y,z])
    vector = vector / np.linalg.norm(vector)
    return vector

def rot_vectang(th,r):
    """ 
    Rotate the given rotation matrix r [ca<-sa] by a 
    random axis with th degree.
    """
    delta, phi = rot_axis() # random delt and phi
    #print 'delta,phi', delta,phi
    v = polar2vect(delta, phi)
    rot = vector_ang(u=v, th=th)
    newr = np.dot(rot, r)
    return newr

    
def sample_mmm():
    """
    orthonormal sample symmetry operators
    """
    m0 = [[ 1, 0, 0], [0, 1, 0], [0, 0, 1]] # itself
    m1 = [[ 1, 0, 0], [0,-1, 0], [0, 0,-1]]
    m2 = [[-1, 0, 0], [0, 1, 0], [0, 0,-1]]
    m3 = [[-1, 0, 0], [0,-1, 0], [0, 0, 1]]
    mmm = [m1,m2,m3]
    return mmm

def rot_axis():
    """
    Random rotation axis generator
    """
    import numpy as np
    pi = np.pi
    delta = random.uniform(-1.,1.)
    delta = np.arccos(delta)                 # arc
    phi = random.uniform(-pi, pi) # rotation
    return delta, phi   #radians...

def miller2euler(hkl, uvw):
    """
    RD // uvw  // 1
    ND // hkl  // 3

    3 x 1 = 2
    """
    # mat [ca<-sa]
    mat = miller2mat(hkl, uvw)
    phi1, phi, phi2 = eul(a=mat, echo=False)
    return phi1, phi, phi2

def miller2mat(hkl, uvw):
    """
    RD // uvw  // 1
    ND // hkl  // 3

    3 x 1 = 2
    """
    import numpy as np
    from euler import euler
    uvw = map(int, uvw)
    hkl = map(int, hkl)

    uvw = np.array(uvw)
    hkl = np.array(hkl)

    # x, y, z bases vectors.
    cx = uvw/np.linalg.norm(uvw)
    cz = hkl/np.linalg.norm(hkl)
    cy = np.cross(cz, cx)

    mat = np.zeros((3,3))
    for i in range(3):
        mat[i][0] = cx[i]
        mat[i][1] = cy[i]
        mat[i][2] = cz[i]

    # mat [ca<-sa]
    return mat

def rolling_fcc():
    ## copper
    ## brass
    ## S
    ## cube
    ## Goss
    pass

def main(hkl,uvw,w0,ngr=1,ifig=10):
    import upf
    import matplotlib.pyplot as plt
    phi1,phi,phi2 = miller2euler(hkl,uvw)
    mat0 = eul(phi1,phi,phi2,echo=False) # ca<-sa
    # print 'mat0:', mat0
    mmm = sample_mmm()
    mats = [mat0]
    # print mats
    for i in range(len(mmm)):
        mats.append(np.dot(mat0,mmm[i]))

    gauss = random.gauss
    r = random.random

    # print 'ind_mat:', ind_mat
    # print 'ngr:', ngr

    grs = []

    for i in range(ngr):
        ind_mat = int(r() * len(mats))
        if w0==0: w=0
        else: w = gauss(mu=0., sigma=w0)
        # print 'w:',w
        # print 'mats:', mats[ind_mat]
        C = rot_vectang(th=w, r=mats[ind_mat])
        # print C
        phi1, phi, phi2 = eul(a=C, echo=False)
        grs.append([phi1,phi,phi2,1])
        #print phi1,phi,phi2

    #print grs

    mypf = upf.polefigure(grains=grs,csym='cubic')
    mypf.pf(mode='contourf',cmode='gray_r',ifig=ifig)
    #mypf.pf(mode='dot',cmode='gray_r',ifig=ifig)
    plt.tight_layout()
    
    hkl = '%i%i%i'%(hkl[0],hkl[1],hkl[2])
    uvw = '%i%i%i'%(uvw[0],uvw[1],uvw[2])
    fn = 'hkl_%s_uvw_%s_th_%i_ngr_%i.txt'%(hkl,uvw,w0,ngr)
    f = open(fn,'w')
    f.writelines(' dum\n dum\n dum\n B  %i \n'%ngr)
    for i in range(ngr):
        f.writelines('%8.3f    %8.3f    %8.3f    %12.7f\n'%(
                grs[i][0],grs[i][1],grs[i][2],1./ngr))

    f.close()
    print '%s has created'%fn
