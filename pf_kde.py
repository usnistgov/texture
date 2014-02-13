"""
Kernel density estimate module for pole figure plotting
"""
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import upf
def ex(fact=100, res=1.0,filename=None,ngrain=None):

    kde = stats.gaussian_kde

    mypf = upf.polefigure(ngrain=ngrain, filename=filename, csym='cubic')
    mypf.pf(pole=[[1,0,0]],mode='contourf')

    # PF[# of pole, mgrid(rot=azmthl), ngrid(tilt= declination)]
    PF = mypf.pfnodes[0]

    mn, nn = PF.shape
    dbeta  = 360./(mn -1) # Rotating
    dalpha =  90./(nn-1)  # tilting
    beta   = np.linspace(0.,360., mn) # 0.~360.
    alpha  = np.linspace(0., 90., nn) # 0.~90.

    #dum1 = np.arange(0.,360. + res/10., res)
    #dum2 = np.arange(0.,90.  + res/10., res)

    x = []
    y = []
    #len(dum2)))
    # |---------------|
    # |   |   |   |   |
    # |---------------|
    # |   |   |   |   |
    # |---------------|
    # |   |   |   |   |
    # |---------------|
    # |   |   |   |   |
    # |---------------|

    maxint = max(PF.flatten())
    for i in range(len(PF)):
        for j in range(len(PF[i])):
            nsample = int(PF[i][j]/maxint*fact)
            for k in range(nsample):
                x.append(beta[i])
                y.append(alpha[j])


    # For a new resolution
    beta  = np.arange(0.,360.+res/2.,res)
    alpha = np.arange(0.,90.+res/2.,res)
    A, B = np.meshgrid(alpha, beta)
    positions = np.vstack([B.ravel(), A.ravel()])
    values = np.vstack([x,y]) # beta, alpha
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel.evaluate(positions).T, A.shape)

    plt.figure(3)
    ax = plt.gca()
    ax.set_frame_on(False)
    ax.set_aspect('equal')
    ax.set_axis_off()
    rx, ry = upf.circle()
    ax.plot(rx,ry,'k')
    pi = np.pi

    nm = (360.0-0.)/res; nn = (180.-90.)/res
    theta = np.linspace(pi, pi/2., nn+1)
    phi = np.linspace(0.,2.*pi, nm+1)
    r = np.sin(theta)/(1-np.cos(theta))
    R, PHI = np.meshgrid(r,phi)
    PHI = PHI + pi/2.
    x = R*np.cos(PHI); y= R*np.sin(PHI)

    ax.contourf(x, y, Z)

def hist2eqsmp(x,pf,fact=100):
    """
    Histogram type PF node to a group of equally sampled points

    inten
    |
    |           |
    |     |     |  |  |
    |  |  |  |  |  |  |
    |--------------------->   x

    |
    |           
    |                 
    |  | ||  | ||| || ||
    |--------------------->
    """
    datpoints = []
    maxpf = []
    maxpf = max(pf)
    for i in range(len(pf)):
        nsample = int(pf[i] / maxpf * fact)
        for j in range(nsample):
            datpoints.append(x[i])
    return np.array(datpoints)
