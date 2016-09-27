"""
Forked from var_gam_fiber.py to extend to more general BCC fibers.
- gamma, alpha, eta, and epsilon fibers.

## variable typical BCC rolling components

# (hkl) // ND, random vector // RD
# (uvw) // RD, random vector // ND
# (xyz) // TD, random vector // RD
"""

import numpy as np
import random
from euler import euler
gauss = random.gauss
expov = random.expovariate
logno = random.lognormvariate
norma = random.normalvariate
from sym import __mmm__ as mmm

def rolling(fibers=['gamma','alpha','eta','epsilon','sigma'],
            wgts=[0.7,0.1,0.1,0.1],ifig=3,sigma=15,
            ngr=100):
    """
    """
    fn='rolling_texture_%s.cmb'%str(ngr).zfill(5)
    f = open(fn,'w')
    f.writelines('Artifical rolling texture for bcc\n')
    f.writelines('Combinations of:')
    for i in xrange(len(fibers)):
        f.writelines('%3i)  %s; '%(i+1,fibers[i]))
    f.writelines('\n')
    f.writelines('weigts of comp :')
    for i in xrange(len(fibers)):
        f.writelines('%3.1f;  '%wgts[i])
    f.writelines('\n')

    import upf
    import matplotlib.pyplot as plt
    wgts = np.array(wgts)
    wgts = wgts/wgts.sum()

    total_pop = []

    for i in xrange(len(fibers)):
        gr = main(ngrains=ngr/len(fibers),sigma=sigma,
                  iopt=1,iexit=True,fiber=fibers[i])
        grt = gr.T
        grt[-1] = grt[-1] * wgts[i]
        gr = grt.T
        for j in xrange(len(gr)):
            total_pop.append(gr[j])
    total_pop = np.array(total_pop)

    f.writelines('B %i\n'%len(total_pop))
    for i in xrange(len(total_pop)):
        f.writelines('%+7.3f %+7.3f %+7.3f %+13.4e\n'%(
            total_pop[i][0], total_pop[i][1], total_pop[i][2], 1./len(total_pop)))

    upf.cubgr(gr=total_pop,ifig=ifig)
    return total_pop

def main(ngrains=100,sigma=5.,iopt=1,ifig=1,fiber='gamma',
         iexit=False,ipfplot=False):
    """
    Arguments
    =========
    ngrains = 100
    sigma   = 5.
    iopt    = 1 (1: gauss (recommended); 2: expov; 3: logno; 4: norma)
    ifig    = 1
    fiber   = 'gamma', 'alpha', 'eta', 'epsilon', 'sigma', 'random'
    ipfplot = False

    Returns
    -------
    It returns the poly-crystal aggregate in the form of numpy's ndarray.
    """
    import upf
    import matplotlib.pyplot as plt
    import cmb
    # h = mmm() ## m-m-m sample symmetry is applied.
    h = [np.identity(3)]
    gr = []
    for i in xrange(ngrains):
        dth = random.uniform(-180., 180.)
        if fiber in ['gamma','alpha', 'eta', 'epsilon', 'sigma']:
            g = gen_gr_fiber(dth,sigma,iopt,fiber)
        elif fiber=='random':
            g = cmb.randomGrain(360,360,360)

        for j in xrange(len(h)):
            temp = np.dot(g,h[j].T)
            phi1,phi,phi2 = euler(a=temp, echo=False)
            gr.append([phi1,phi,phi2,1./ngrains])

    if iexit: return np.array(gr)

    fn = '%s_fib_ngr%s_sigma%s.cmb'%(fiber,str(len(gr)).zfill(5),str(sigma).zfill(3))
    f = open(fn,'w')
    f.writelines('Artificial %s fibered polycrystal aggregate\n'%fiber)
    f.writelines('bcc_rolling_fiber.py python script\n')
    f.writelines('distribution: ')
    if iopt==1: f.writelines(' gauss')
    if iopt==2: f.writelines(' expov')
    if iopt==3: f.writelines(' logno')
    if iopt==4: f.writelines(' norma')
    f.writelines('\n')
    f.writelines('B %i\n'%ngrains)
    for i in xrange(len(gr)):
        f.writelines('%+7.3f %+7.3f %+7.3f %+13.4e\n'%(
            gr[i][0], gr[i][1], gr[i][2], 1./len(gr)))

    if ipfplot:
        mypf1 = upf.polefigure(grains=gr,csym='cubic')
        mypf1.pf_new(poles=[[1,0,0],[1,1,0],[1,1,1]],ix='RD',iy='TD')
        fig=plt.gcf()
        fig.tight_layout()
        print 'aggregate saved to %s'%fn
        fig.savefig(
            '%s_contf.pdf'%(fn.split('.cmb')[0]),
            bbox_inches='tight')
        fig.savefig(
            '%s_contf.png'%(fn.split('.cmb')[0]),
            bbox_inches='tight')
        fig.clf()
        plt.close(fig)

    return np.array(gr),fn

def gen_gr_fiber(th,sigma,iopt,fiber='gamma'):
    """
    Generate rotation matrix and
    'perturb' the orientation following the
    given 'distribution' function <iopt>
    with the given standard deviation <sigma>
    """
    from text import miller2mat, miller2mat_RT
    if fiber=='gamma':
        hkl, uvw = hkl_gamma()
        g_casa = miller2mat(hkl, uvw)
        g_sasa = nd_rot(th)
    elif fiber=='alpha':
        hkl, uvw = hkl_alpha()
        g_casa = miller2mat(hkl, uvw)
        g_sasa = rd_rot(th)
    elif fiber=='eta':
        hkl,uvw = hkl_eta()
        g_casa = miller2mat(hkl, uvw)
        g_sasa = rd_rot(th)
    elif fiber=='epsilon':
        xyz, uvw = hkl_epsilon()
        g_casa = miller2mat_RT(uvw=uvw,xyz=xyz)
        g_sasa = td_rot(th)
    elif fiber=='sigma':
        hkl, uvw = hkl_sigma()
        g_casa = miller2mat(hkl,uvw)
        g_sasa = nd_rot(th)

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
    u[3] = vector axis about which the rotation occurs
    th   = radian angle (degree of the rotation)
    """
    u = np.array(u)
    u = u/np.sqrt((u**2).sum()) # normalize
    idx = np.identity(3)
    ct = np.cos(th*np.pi/180.)
    st = np.sin(th*np.pi/180.)
    cm = crossop(u) # cross product operator
    r = np.zeros((3,3))
    for i in xrange(3):
        for j in xrange(3):
            r[i,j] = idx[i,j] * ct + st * cm[i,j] + \
                     (1. - ct ) * u[i] * u[j]
    return r

def crossop(u):
    """
    Cross operator can be used to obtain the cross
    product as a matrix-vector product

    a x b  = [a_x]_ij [b]_j

    Argument
    ========
    u
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
    return hkl, uvw

def hkl_eta():
    uvw = [0,0,1] # major // RD
    hkl = [1,0,0] # minor // ND
    return hkl, uvw

def hkl_epsilon():
    uvw = [0,0,1] # major // RD
    xyz = [1,1,0] # minor // TD
    return xyz, uvw

def hkl_sigma():
    hkl = [1, 0, 0] # major // ND
    uvw = [0, 0, 1] # minor // RD
    return hkl, uvw
