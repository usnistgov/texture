"""
Convert Sectioned OD to LABOTEXT format

-- version 2013 June 20
  * mmm sample symmetry fixed

- example
>>> import main(fn='A_steel.cod ', odfn='labo.txt')
"""

print __doc__

def reader(fn=None,nomen='b'):
    """
    Reads a file consisting of sections of 'phi2=constant'
    and return total section blocks and maximum phi2

    Argument
    ========
    fn = None
    """
    datlines = open(fn, 'r').readlines()
    blocks   = []
    temp     = []
    for i in range(len(datlines)):
        temp.append(datlines[i])
        if len(datlines[i]) < 3:
            blocks.append(temp)
            temp = []

    ## Hardwired.
    maxphi2=90.
    dphi2=5.

    if nomen=='b':
        pass
    elif nomen=='k':
        pass

    print 'maximum phi2', maxphi2
    print 'phi2 increment:', dphi2
    return blocks, maxphi2

def readblock(block,nomen='b'):
    """
    Use this method for the results obtained through 'Harmonics' method.

    Read a block and returns
    1) section's intensity mesh
    2) Phi  max
    3) phi1 max

    Argument
    ========
    block
    """
    import numpy as np
    header = block[0].split()
    dphi    = float(block[1][5:10])
    phimx   = float(block[1][10:15])
    dphi1   = float(block[1][15:20])
    phi1mx  = float(block[1][20:25])
    # print 'dphi, phimx, dphi1, phi1mx',
    # print dphi, phimx, dphi1, phi1mx

    if phi1mx==180: pass
    elif phi1mx==90: pass
    else: raise IOError, 'Not expected maximum phi value...'

    phi  = np.zeros(np.arange(0., phimx  + 0.001, dphi ).shape)
    phi1 = np.zeros(np.arange(0., phi1mx + 0.001, dphi1).shape)
    section = np.zeros((len(phi), len(phi1)))
    block = block[2:]

    if phi1mx==180:
        for i in range((len(block)-1)/2):
            arim = block[i*2][:18*4+1][1:] + block[i*2+1][:19*4+1][1:]
            for j in range(len(arim[::4])):
                section[i,j] = float(arim[4*j:4*j+4])
    elif phi1mx==90:
        for i in range(len(block)-1):
            arim = block[i][:19*4+1][1:]
            for j in range(len(arim[::4])):
                section[i,j] = float(arim[4*j:4*j+4])

    # # block = block[::-1][0:]
    # # block = block[::-1]

    # for i in range(len(block)-1):
    #     dum = block[i].split()
    #     section[i] = map(float, dum)

    if nomen=='b':
        section = section.T # (phi, phi1) -> (phi1, phi)
    elif nomen=='k':
        seciton ## (Theta, PSI)
    return section, phimx, phi1mx

def __cod2labo__(fn=None):
    """
    Argument
    =========
    fn = None
    nomen='b'
    """
    import numpy as np


    blocks, phi2mx = reader(fn=fn)
    sct, phimx, phi1mx = readblock(blocks[0])

    nphi2 = len(blocks)-1
    nphi1 = len(sct)
    nphi  = len(sct[0])

    cod = np.zeros((nphi1, nphi, nphi2))

    # len(blocks)
    for ib in range(len(blocks)-1): # nphi2
        section, pmx, pmx1 = readblock(blocks[ib])
        for i in range(len(section)): # nphi1
            for j in range(len(section[i])): #nphi
                cod[i,j,ib] = section[i,j]

    return cod, phi1mx, phimx, phi2mx

def main(fn=None, odfn='labo.txt',fnout='dum', rot=0,
         nomen='b'):
    """
    Writes COD to Labotex convention file

    arguments
    =========
    fn    = None
    odfn  = 'labo.txt'
    fnout = 'dum'
    rot   = 0 (in plane rotation about the normal direction (z)
    """
    import numpy as np
    cod, phi1mx, phimx, phi2mx = __cod2labo__(fn=fn)
    print phi1mx, phimx, phi2mx

    f = open(odfn, 'w')
    phi1 = np.linspace(0, phi1mx+0.001, len(cod))
    phi  = np.linspace(0, phimx+0.001, len(cod[0]))
    phi2 = np.linspace(0, phi2mx+0.001, len(cod[0][0]))

    cod = np.array(cod).swapaxes(2,0) # phi2, phi, phi1
    cod = np.array(cod).swapaxes(0,1) # phi, phi2, phi1

    mmm = False
    if phi1mx==90 and phimx==90 and phi2mx==90:
        print 'mmm sample symmetry must have been applied for cubic crystal'
        mmm = True

    print 'cod.shape', cod.shape

    f.write('phi1   phi2  phi  COD\n')
    for i in range(len(cod)): # phi
        if nomen=='b': p = phi[i]
        if nomen=='k': p = phi[i]-90.
        for j in range(len(cod[i])): #phi2
            if nomen=='b': p2 = phi2[j]
            if nomen=='k': p2 = -phi2[j]
            for k in range(len(cod[i][j])): #phi1
                if nomen=='b': p1 = phi1[k]
                if nomen=='k': p1 = -phi1[k]-90
                f.write('%4.1f  %4.1f  %4.1f  %6.3f \n'%(
                        p1, p2, p, cod[i][j][k]))

    f.close()

    ## write into weighted discrete grains files
    odf2dg(odfn=odfn, mmm=mmm, rot=rot, fnout=fnout)

def odf2dg(odfn=None, mmm=False, rot=0, fnout=None,ng=None):
    """
    write into weighted discrete grain files

    Arguments
    =========
    odfn = None
    mmm  = False
    rot  = 0
    """
    import os
    import matplotlib.pyplot as plt
    import numpy as np
    os.sys.path.append('/Users/yj/Dropbox/devel/EVPSC_f2py/pyscripts')
    os.sys.path.append('/Users/yj/Dropbox/devel/EVPSC_f2py/pyscripts/pf')

    import cmb, upf
    if type(ng)!=type(None):
        ng = [20000]

    if mmm: ng = np.array(ng) / 4
    for igr in range(len(ng)):
        cmb.main(odf=odfn, ngrain=ng[igr], outputfile='dum.tex', iplot=False)
        # sample symmetry recovery when mmm has been applied (reduced 90x90x90)
        if mmm: grains = sample_mmm(fn='dum.tex')
        else: grains = np.loadtxt('dum.tex', skiprows=4)

        os.remove('dum.tex')
        grains = inplanerot(rot=rot,grains=grains)

        fout = open('%s_%i.cmb'%(fnout, len(grains)), 'w')
        fout.writelines('dummy\ndummy\ndummy\n B  %i\n '%len(grains))

        for i in range(len(grains)):
            fout.writelines('%.5f %.5f %.5f  %.7e\n'%(
                grains[i][0],grains[i][1],grains[i][2],grains[i][3]))
        upf.cubgr(gr=grains,ifig=igr,poles=[[1,1,0],[2,0,0],[2,1,1]])
        plt.gcf().savefig('pf_%i.pdf'%ng[igr])

def sample_mmm(fn=None):
    """
    Apply the orthorhombic sample symmetry 2/m 2/m 2/m (shortened as mmm)
    360/2 = 180 rotation

    mmm contains: 4 operators (itself, m, m, m)
    """
    import numpy as np
    from euler import euler

    m0 = [[ 1, 0, 0], [0, 1, 0], [0, 0, 1]]
    m1 = [[ 1, 0, 0], [0,-1, 0], [0, 0,-1]]
    m2 = [[-1, 0, 0], [0, 1, 0], [0, 0,-1]]
    m3 = [[-1, 0, 0], [0,-1, 0], [0, 0, 1]]

    mmm = [m1,m2,m3]

    grains   = np.loadtxt(fn,skiprows=4)
    grains_t = grains.T
    wgt      = grains_t[-1][::]
    grains_t[-1] = grains_t[-1]/sum(wgt)
    grains   = grains_t.T
    grains_new = []
    for i in range(len(grains)):
        amat = euler(ph=grains[i][0], th=grains[i][1],
                     tm=grains[i][2], echo=None) # ca<-sa
        amat0 = [amat]
        #amat0.append(np.dot(amat, m1))
        for j in range(len(mmm)):
            #n = len(amat0)
            #for k in range(n):
            #amat0.append(np.dot(amat0[k],mmm[j]))
            amat0.append(np.dot(amat0[0],mmm[j]))

        for j in range(len(amat0)):
            phi1, phi, phi2 = euler(a=amat0[j], echo=False)
            grains_new.append([phi1,phi,phi2,grains[i][-1]/len(amat0)])
    return grains_new

def inplanerot(rot=0,grains=None):
    import numpy as np
    if rot==0: return grains
    from euler import euler
    arot = euler(ph=rot, th=0, tm=0, echo=False)
    for i in range(len(grains)):
        phi1, phi, phi2, wgt = grains[i]
        ag = euler(ph=phi1, th= phi, tm=phi2, echo=False) # ca<-sa
        agg = np.dot(ag, arot)
        phi1, phi, phi2 = euler(a=agg, echo=False)
        grains[i][0] = phi1
        grains[i][1] = phi
        grains[i][2] = phi2
    return grains
