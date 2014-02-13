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
    grain generator

    =========
    Arguments

    th       = 0.0
    sigma    = 5.0
    iopt     = 1
    """
    import text
    hkl = [ 1, 1, 1]
    uvw = [-1, 1, 0] # auxilary vector that is perpedicular to hkl

    g_casa = text.miller2mat(hkl, uvw)
    g_sasa = in_plane_rot(th)
    g = np.dot(g_casa,g_sasa)

    if iopt==1: dth=gauss(mu=0., sigma=sigma)
    if iopt==2: dth=expov(sigma)
    if iopt==3: dth=logno(mu=0., sigma=sigma)
    if iopt==4: dth=norma(mu=0., sigma=sigma)

    return rot_vectang(th=dth, r=g)

def in_plane_rot(thet=0.):

    thet * np.pi / 180.
    amat = [[np.cos(thet), - np.sin(thet), 0.],
            [np.sin(thet),   np.cos(thet), 0.],
            [          0.,             0., 1.]]
    return amat

def rot_vectang(th,r):
    """ Rotate the given rotation matrix r
    [ca<-sa] by a random axis with th degree.
    """
    import cs
    delta, phi = rot_rand_axis()
    v = polar2vect(delta,phi)
    rot = cs.vector_ang(u=v, th=th)
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
