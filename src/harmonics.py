def extPF(PF):
    import numpy as np
    mn, nn = PF.shape
    dbeta  = 360. / (mn-1) # rotating
    dalpha =  90. / (nn-1) # tilting
    beta = np.linspace(0,360., mn) # 0.~360.
    alpha = np.linspace(0.,90., nn) # 0.~90.

    nn2 = 180/dalpha + 1
    alpha2 = np.linspace(0.,180., nn2)

    ## I want to fill alaph from 90~ 180
    # PF[rot, tilt]
    PF_ext = np.zeros((mn, nn2))

    for i in range(len(PF_ext)):
        b = beta[i]
        for j in range(len(PF_ext[i])):
            a = alpha2[j]
            if a>90.:
                a0 = 180 - a    # alpha for PF
                j0 = a0/dalpha  #
                b0 = b + 180
                if b0> 360: b0 = b0 - 360
                i0 = b0/dbeta
            else: i0 = i; j0 = j
            PF_ext[i,j] = PF[i0,j0]
    return PF_ext

def redPF(PF):
    import numpy as np
    mn, nn = PF.shape
    dbeta  = 360. / (mn-1) # rotating
    dalpha = 180. / (nn-1) # tilting
    beta   = np.linspace(0.,360.,mn) # 0.~360.
    alpha  = np.linspace(0.,180.,nn) # 0.~180.

    nn2 = 90/dalpha + 1
    alpha2 = np.linspace(0.,90.,nn2)
    PF_red = np.zeros((mn, nn2))

    for i in range(len(PF_red)):
        b = beta[i]
        for j in range(len(PF_red[i])):
            PF_red[i,j] = PF[i,j]
    return PF_red

#----------------------------------------------------------------------#
# Fitting grid pole figure with the Harmonics series...                #
def Qlm(PF=None, dalpha=None, dbeta=None, alpha=None,
        beta=None, l=23, m=3):
    """
    Returns Q_lm, the coefficients

    m shouldn't be larger than l
    l shouldn't be a negative value

    This def returns Q_lm, the harminics coefficients.

    Arguments
    ========
    PF     = Pole figure intensity grid[m, n]
    dalpha = delta alpha
    dbeta  = delta beta
    alpha  = array of alpha   (The tilting angle from the centeral pole)
    beta   = array of thetas  (Rotating angle)
    l      = 23
    m      = 3
    """
    import numpy as np
    im = complex(0.,1.)
    Q_lm = complex(0.,0.) # complex number

    ALP = []
    for i in range(len(alpha)):
        ALP.append(AssociatedLegendreP(l,m,alpha[i])[-1])

    for j in range(len(alpha)):
        for i in range(len(beta)):
            P_lm =  ALP[j]
            Q_lm = Q_lm +\
                PF[i,j] *\
                P_lm *\
                np.exp(-im * m * beta[i]) *\
                np.sin(alpha[j]) *\
                dbeta *\
                dalpha

    #print 'Qlm', Q_lm
    return Q_lm

def Qlm2(PF, Plm, dalpha, dbeta, alpha, beta, l, m):
    import numpy as np
    im = complex(0.,1.)
    Q_lm = complex(0.,0.) # complex number

    for j in range(len(alpha)):
        for i in range(len(beta)):
            Q_lm = Q_lm + PF[i,j] * Plm[j] *\
                np.exp(-im * m * beta[i]) *\
                np.sin(alpha[j]) *\
                dbeta *\
                dalpha
    return Q_lm


def AssociatedLegendreP(l,m,alpha):
    """
    Returns polynomial P_l^m for (l=0, ... l) for the given m.
    """
    from scipy.special import lpn
    import numpy as np
    from numpy import cos, sin

    dat = lpn(l, cos(alpha))  # Pl, and d(Pl)/dx
    LegPol = dat[0]
    dLdx = dat[1]  # 0,1,2,... l

    # dLdx is in dimension: 0,1,2,...l
    np.zeros((l+1))
    #return (-1)**m * sin(alpha) * dLdx # 0,1,2,...l for the given m.
    return (-1)**m * (sin(alpha))**abs(m) * dLdx

def harm_PF(PF, l=23):
    """
    PF is assumed to be a 'node' as shown below:

    --------------------------------------------------
    ** nodes = np.zeros((mgrid+1,ngrid+1)) , m: rot, n: tilt
     The mesh nodes is on the central point, 'o', of
    each grid as depicted below:

    90o---o---o
      |   |   |
      .
      .
      .

      o---o---o
      |   |   |
    10o---o---o
      |   |   |   . . .
     5o---o---o
      |   |   |
     0o---o---o
      0   5  10 ....  360

    beta (rotating)
    alpha(tilting)

    Arguments
    =========
    PF : Pole figure nodes from polefigure class
    l  : Order of harmonics
    """
    import numpy as np
    import time
    mn, nn = PF.shape
    dbeta  = 360. / (mn-1) # rotating
    dalpha =  180. / (nn-1) # tilting
    beta = np.linspace(0,360., mn) # 0.~360.
    alpha = np.linspace(0.,180., nn) # 0.~90.

    #print 'Rot dbeta', dbeta
    #print 'Tilt dalpha', dalpha

    beta = beta     * np.pi / 180.
    alpha = alpha   * np.pi / 180.
    dbeta = dbeta   * np.pi / 180.
    dalpha = dalpha * np.pi / 180.

    im = complex(0, 1)
    PF_ab = np.zeros(np.array(PF).shape, dtype='complex')


    t0 = time.time()

    ALP = np.zeros((len(alpha), l+1, len(np.arange(-l, l+1))  ))
    for a0 in range(len(alpha)):
        ms = np.arange(-l, l+1)
        for m0 in range(len(ms)): #  -l ... + l
            polynomial = AssociatedLegendreP( # dum : 0 ... m
                l=l, m=ms[m0],
                alpha = alpha[a0]) # ... 0... l
            for l0 in range(l+1):
                ALP[a0,l0,m0] = polynomial[abs(ms[m0])]

    print '%5.2f seconds elasped to calculate'%(time.time()-t0)+\
        'Associated Legendre Polynomilal'


    for l0 in range(l):
        print 'l=',l0
        print 'm=',
        ms = np.arange(-l0, l0+1)
        for m0 in range(len(ms)):
            print ms[m0],
            # Plm should be [m0]
            # ALP is [a0, l0, m0]
            # ALP.swapaxes[0,2] -> [m0, l0, a0]
            Q_lm = Qlm2(PF=PF,
                        Plm=ALP.swapaxes(0,2)[m0, l0],
                        dalpha=dalpha, dbeta = dbeta,
                        beta=beta, alpha = alpha, l=l0, m=ms[m0])
            for b0 in range(len(beta)):
                bet = beta[b0]
                for a0 in range(len(alpha)):
                    alph = alpha[a0]
                    P_lm = ALP[a0,l0,m0]
                    dum = Q_lm * P_lm * \
                        np.exp(im * ms[m0] * bet)
                    PF_ab[b0,a0] = PF_ab[b0,a0] + dum
        print '\n'

    return PF_ab

def harm_pf(grains=None, filename=None, l=2, dm=7.5, dn=7.5,
            pole=[[1,0,0]], csym='cubic'):
    """
    Something is wrong and I couldn't figure out.
    Further development is defferred.

    Arguments
    =========
    grains   = None
    filename = None
    l        = 10
    dm       = 7.5
    dn       = 7.5
    pole     = [[1, 0, 0]]
    csym     = 'cubic'
    """
    import numpy as np
    import cmb
    from upf import polefigure, circle
    pi = np.pi

    if grains==None and filename==None:
        gr = cmb.random(filename='dum.tex',ngrain=1000,
                        phi1=360., phi=90., phi2=180.)
        mypf = polefigure(grains=grains, csym=csym)
    else: mypf = polefigure(grains=grains, filename=filename, csym=csym)

    mypf.pf(pole=pole, dm=dm, dn=dn)
    hpf = []

    import matplotlib.pyplot as plt
    fact = 2. #size factor
    figsize = (len(pole)*2.*fact, 1.*2.*fact)
    fig = plt.figure(33, figsize=figsize)
    ax = plt.gca()

    for ip in range(len(pole)):
        extended_PF = extPF(mypf.pfnodes[ip])
        print 'extended_PF.shape', extended_PF.shape
        #dum = harm_PF(PF=mypf.pfnodes[ip], l=l)
        dum = harm_PF(PF=extended_PF, l=l)
        reduced_PF = redPF(dum)
        hpf.append(reduced_PF)

#        print 'max:', np.real(max(hpf[ip].flatten()))
#        print 'min:', np.real(min(hpf[ip].flatten()))

    theta = np.linspace(pi, pi/2., (90.)/dn + 1)
    phi = np.linspace(0.,2.*pi,  (360.)/dm+1)
    r = np.sin(theta)/(1-np.cos(theta))
    R, PHI = np.meshgrid(r,phi)
    PHI = PHI + pi/2.
    x = R*np.cos(PHI); y = R*np.sin(PHI)

    #    return x, y, mypf.pfnodes[0]
    cnt = ax.contourf(x,y, hpf[0])
    ax.set_frame_on(False)
    ax.set_axis_off()
    ax.set_aspect('equal')
    rx, ry = circle()
    ax.plot(rx, ry, 'k')
    ax.set_xlim(-1.2,1.5)
    ax.set_ylim(-1.2,1.5)

    tcolors = cnt.tcolors
    clev = cnt._levels
    for i in range(len(tcolors)):
        cc = tcolors[i][0][0:3]
                    #if levels==None:
        if ip==len(pole)-1:
                        ## level line
            ax.plot([1.28, 1.35],
                    [1. - i * 0.2, 1. - i * 0.2],
                    color=cc)
                        ## level text
            ax.text(x=1.40, y= 1. - i*0.2 - 0.05,
                    s='%3.2f'%(clev[i]),
                    fontsize=4*fact)

    return hpf, mypf.pfnodes
