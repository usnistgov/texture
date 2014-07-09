## variable gamma fiber...

## model creation of gamma fiber (111) // ND
# (hkl) // ND, random vector // RD

import numpy as np
import random
from euler import euler
gauss = random.gauss
expov = random.expovariate
logno = random.lognormvariate
norma = random.normalvariate
from sym import __mmm__ as mmm

def main(ngrains=100,sigma=5.,iopt=1,ifig=1):
    """
    arguments
    =========
    ngrains = 100
    sigma   = 5.
    iopt    = 1 (1: gauss (recommended); 2: expov; 3: logno; 4: norma)
    ifig    = 1
    """
    import upf
    import matplotlib.pyplot as plt
    h = mmm()
    gr = []
    for i in range(ngrains):
        dth = random.uniform(-180., 180.)
        g = gen_gamma_gr(dth, sigma, iopt=iopt)
        for j in range(len(h)):
            temp = np.dot(g,h[j].T)
            phi1,phi,phi2 = euler(a=temp, echo=False)
            gr.append([phi1,phi,phi2,1./ngrains])

    fn = 'gam_fib_ngr%s_sigma%s.cmb'%(str(len(gr)).zfill(5),str(sigma).zfill(3))
    f = open(fn,'w')
    f.writelines('Artificial gamma fibered polycrystal aggregate\n')
    f.writelines('var_gam_fiber.py python script\n')
    f.writelines('distribution: %s')
    if iopt==1: f.writelines(' gauss')
    if iopt==2: f.writelines(' expov')
    if iopt==3: f.writelines(' logno')
    if iopt==4: f.writelines(' norma')
    f.writelines('\n')
    f.writelines('B %i\n'%ngrains)
    for i in range(len(gr)):
        f.writelines('%7.3f %7.3f %7.3f %13.4e\n'%(
            gr[i][0], gr[i][1], gr[i][2], 1./len(gr)))

    upf.cubgr(gr=gr,ifig=ifig)
    plt.figure(ifig).savefig('%s.pdf'%(fn.split('.cmb')[0]))

    return np.array(gr)

def gen_gamma_gr(th=0.,sigma=5.,iopt=1):
    """
    grain generator for gamma fiber

    =========
    Arguments

    th       = 0.0
    sigma    = 5.0
    iopt     = 1
    """
    from text import miller2mat

    ## gamma fiber: axisymmetry about ND
    hkl, uvw = hkl_gamma()
    g_casa = miller2mat(hkl, uvw)
    g_sasa = nd_rot(th)
    ##

    g = np.dot(g_casa,g_sasa)
    if iopt==1: dth=gauss(mu=0., sigma=sigma)
    if iopt==2: dth=expov(sigma)
    if iopt==3: dth=logno(mu=0., sigma=sigma)
    if iopt==4: dth=norma(mu=0., sigma=sigma)

    return rot_vectang(th=dth, r=g)

def rd_rot(thet):
    return vector_ang(u=[1,0,0],th=thet)
def td_rot(thet):
    return vector_ang(u=[0,1,0],th=thet)
def nd_rot(thet):
    return vector_ang(u=[0,0,1],th=thet)


def rot_vectang(th,r):
    """ Rotate the given rotation matrix r
    [ca<-sa] by a random axis with th degree.
    """
    #import cs
    delta, phi = rot_rand_axis()
    v = polar2vect(delta,phi)
    #rot = cs.vector_ang(u=v, th=th)
    rot = vector_ang(u=v, th=th)
    newr = np.dot(rot,r)
    return newr

def rot_rand_axis():
    delta = random.uniform(-1.,1.)
    delta = np.arccos(delta)
    phi   = random.uniform(-np.pi, np.pi)
    return delta, phi

def polar2vect(delta, phi):
    """ Given delta and phi returns a vector"""
    x = np.cos(phi)
    y = np.sin(phi)
    z = np.sin(delta)
    vector = np.array([x,y,z])
    vector = vector / np.linalg.norm(vector)
    return vector

def vector_ang(u,th):
    """
    implementation of subroutine vector_ang in cs.f

    arguments
    =========
    u[3] = vector axis about whic the rotation occurs
    th   = radian angle (degree of the rotation)
    """
    u = np.array(u)
    u = u/np.sqrt((u**2).sum()) # normalize

    idx = np.identity(3)

    ct = np.cos(th*np.pi/180.)
    st = np.sin(th*np.pi/180.)
    cm = crossop(u) # cross product operator

    # print th
    # print ct, st
    # print cm
    # raw_input()

    r = np.zeros((3,3))

    for i in range(3):
        for j in range(3):
            r[i,j] = idx[i,j] * ct + st * cm[i,j] + \
                     (1. - ct ) * u[i] * u[j]

    return r

def crossop(u):
    """
    Cross operator can be used to obtain the cross
    product as a matrix-vector product

    a x b  = [a_x]_ij [b]_j
    """
    m=np.zeros((3,3))
    m[0][1] = -u[2]
    m[0][2] =  u[1]
    m[1][0] =  u[2]
    m[1][2] = -u[0]
    m[2][0] = -u[1]
    m[2][1] =  u[0]

    return m

# hkl // ND
# uvw // RD
# xyz // TD

""" 
Only major diection is relevant
Any auxilary vector that is perpendicular to the
major will work it out.
 """
def hkl_gamma():
    hkl = [ 1, 1, 1] # major // ND
    uvw = [-1, 1, 0] # minor // RD
    return hkl, uvw

def hkl_alpha():
    uvw = [ 1, 1, 0] # major // RD
    hkl = [ 0, 0, 1] # minor // ND

def hkl_eta():
    uvw = [0,0,1] # major // RD
    hkl = [1,0,0] # minor // ND

def hkl_epsilon():
    xyz = [1,1,0] # major // TD
    uvw = [0,0,1] # minor // RD
