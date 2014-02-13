"""
ver: 2011-13-Mar

Features:

  crystal symmetries : cubic, hexagonal (as of 2011-14-March)
    This module can easily exnted to other symmetries but have not
    been completed yet.

Pole figure and Inverse pole figure plotting by stereographic
projection. Both contour and dot types are available for pole
figures, whereas only dot type is available for inverse pole figure, yet.

It should be machine independent. Make the way treating the data uniform and
intiutive and simple.

--------
examples
--------
 >>> import upf
 >>> mypf = upf.polefigure(ngrain=8000, # or filename ='texture/08000.cmb'
                           csym = 'cubic', cdim=[1.,1.,1.], cang =[90.,90.,90.])
 >>> cnts = mypf.pf(pole=[ [1,0,0],[1,1,0],[1,1,1]], mode='contourf', ifig=2,
                    dm=7.5, dn=7.5, levels=None, cmode='gray_r')

        * Pole figure plotting upon the polycrystal aggregate.
          --> refer to 'class polefigure'

        * Experimental pole figure plotting software where
         the range of khi (tiliting angle) is often limited.
         However, the plotting scheme will be very similiar to the
         pole figure for polycrystalline aggregates, in that
         experimental pole figure often comes as a grid set of data
         in the polar coordinate system (theta as angle, and khi as
         radius)
          --> First of all, I should come up with a uniform format
          having a grid format consisting of theta and khi axes.
          the resolution of grid should be decided after getting the
          data.

          range of data:
          theta : 0. ~ 2.*pi
          khi : 0. ~ pi/2.

          detecting the np.shape and decide the incremental resolution
          of each axis.

          Let's say it will be something like below:

          >>> m, n = np.shape(data)
               # m is the number of grid along theta axis
               # while n is the number of grid along khi

          --> Often, as far as I know, experimental pole figure should be
          normalized since the incremental scanning area is varying with
          increasing khi(tilting angle). The issues related to background
          and defocus effect should be dealt with independently, as well.

          The format controller is also necessary.
          Current cases of interests:
              1. Bruker D8 Discovery (at GIFT POSTECH) (UXD)
              2. Unknown system at Dr. Steglich's institute
              3. Experimental PF (popLA format)
"""
#print __doc__
"""
Updates
  # 1
  2011-16-Mar
  Contour plot is done on the given axes if necessary.

  # 2
  2011-17-April
  The symmetry module is cythonized.
  The __equiv__ method has been taken out from the loop over grains,
  which saves some computational waste.

  # 3
  2011-18-April
  The hardwired smoothing rim is expanded to the 2nd rim.
  Thus now it is n=0, n=1 axes along which the intensities are averaged.

  # 4
  2011-14-Sept
  The pixel view superimposition to contour pole figure is initiated.
  Cython part having been used for 'sym.py' has been removed.

"""
##########################################################################
# For developing the ultimate pole figure plotting software              #
# which can be alternative to the pole8 program by Dr. Tome.             #
# After been tried so many times to wrap the code with f2py              #
# decided to make my own rather than being desperated by the             #
# trivial problems raised during making the legacy code more             #
# objective-oriented.                                                    #
#                                                                        #
# Need a exemplary texture file along with the consideration that the    #
# equivalent information can be given by arrays of grains.               #
#                                                                        #
# EITHER FILE OR ARRAY -- This is going to be start.                     #
# can be given as experimental pole figure or texture file style.        #
# I want to start with texture file style. Let's go for it! (2011-02-23) #
#                                                                        #
# A stable and decent version now completed. A contour type of pole      #
# figure is plottable, whereas inverse pole figure does not come along   #
# with the contour type. (2011-03-13)                                    #
#                                                                        #
#                                                                        #
# inverse and pole figure seems to be working properly. The symmetry     #
# operations is not completed yet. As of now (2011-05-Mar) only cubic    #
# crystal symmetry is considered.                                        #
#                                                                        #
# Symmetry operation has been working fine.                              #
#                                                                        #
##########################################################################
## import library blocks

import numpy as np
import matplotlib.pyplot as plt
import matplotlib #matplotlib as raw
import os, glob, math
from randomEuler import randomEuler as re
from euler import euler #in euler module def euler : A-matrix and Euler angles
import time
import random
# import pp ## parallel

def pfnorm(data):
    """
    experimental incomplete pole figure preliminary normalization

    data format should be correct
    grid (phi, khi)

    resolution should be 5.0 degress in both phi and khi
    phi range: 0~ 355
    khi range: 0~ 5*(nn-1)
    """
    # All angles are in radian
    if len(data)!=72:
        print 'number of phi grid:  %i'%len(data)
        raise IOError, 'Unexpected resolution along phi axis'

    dphi = 360. / len(data)
    dphi = dphi * np.pi / 180. # dphi
    dkhi = 5. * np.pi / 180. # dkhi
    print 'dkhi, dphi', dphi*180./np.pi, dkhi*180./np.pi

    nkhi = len(data[0])
    phi_i = 0.
    phi_f = np.pi * 2.
    khi_i = 0.
    khi_f = dkhi * (nkhi - 1)
    print 'khi range', khi_i, khi_f*180./np.pi
    print 'phi range', phi_i, phi_f*180./np.pi
    # spanned area, i.e., the area of incomplete hemisphere:
    # area = (np.cos(khi_f) - np.cos(khi_i)) * (phi_f - phi_i)

    Nf = 0.

    a = 0.
    b = 0.
    for i in range(len(data)):
        for j in range(len(data[0])):
            a = a + np.sin(dkhi*j)

    for i in range(len(data)):
        for j in range(len(data[i])):
            b = b + data[i,j] * np.sin(dkhi*j)

    # for k in range(len(data)):
    #     for l in range(len(data[k])):
    #         b = b + data[k, l] * np.sin(dkhi*l)

    for i in range(len(data)): #phi
        for j in range(len(data[i])): #khi
            data[i,j] = data[i,j] * a / b
            pass
        pass
    return data

def epfformat(mode=None, filename=None):
    """
    Experimental pole figure format controller
    mode:
      "steglich"
      "bruker"  *.uxd file
      "epf" (2011-Oct-6) epf popLA experimental pole figure format

    Returns the pole figure data as
    the standard format (m x n numpy array)
    each of axes stands for rotating (phi) and tilting (khi)
    angle in the laboratory space.
    These angles will be angle and radius in the
    space onto which a pole figure is projected.

    conventions:
     tilting angle: khi
     rotating angle: phi
     dk: incremental khi angle
     dp: incremental phi angle
     nk: number of points along khi axis
     np: number of points along phi axis

     angles are converted into radian whenever possible
     """

    ## steglich's format
    nskip_calc = 13
    nskip_raw = 28
    ##

    if mode=='steglich':
        # Calculated or raw pole figure
        i = 0
        print 'filename=', filename
        if os.path.isfile(filename)!=True:
            raise IOError, 'file is not available'

        while True:
            try:
                data = np.loadtxt(
                    filename, skiprows=i)
            except: i = i + 1
            else:
                print 'number of skipped rows: %i'%i
                break
            if i>1000: raise IOError, 'something is wrong'
            pass

        ## raw pole figure format
        if i==nskip_raw:
            # axes: (phi, khi)
            data = data.T # (khi, phi)

            # upon each axis
            f = open(filename)
            temp = f.readlines()
            khi = map(float, temp[nskip_raw - 1].split()[1:])
            khi = np.array(khi)
            khi = khi * np.pi/ 180.
            phi = data[0] #first column is for phi
            phi = phi * np.pi/ 180.
            data = data[1:] #remove phi
            data = data.T # (phi, khi) axis

            ## dp and dk
            dp = phi[1] - phi[0]
            dk = khi[1] - khi[0]

            ## shift phi back-modification
            phi, data = shiftphi(data, dp, phi)
            ## shift phi back-modification ends

            ## check if there is phi=0, 360 at the same time
            ## if that's the case remove data along phi=360
            phi, data = duplicate0_360(phi, data)
            ##

            ## check if this is incomplete pole figure
            tiny = 0.000001
            isincomplete=False
            if np.pi/2. - khi[-1] > tiny:
                print 'Incomplete pole figure'
                print 'khi range: %3.2f ~%3.2f'%(
                    khi[0]*180./np.pi, khi[-1]*180./np.pi)
                isincomplete=True

                ## normalization
                if raw_input('y(norm), or n(no)>>>')=='y':
                    data = pfnorm(data)
                    pass

                dum = np.zeros((np.pi*2./dp, np.pi/2./dk+1))
                for i in range(len(data)):
                    for j in range(len(data[i])):
                        dum[i, j] = data[i, j]
                data = dum.copy()
                del dum
                pass
            pass

        ## calculated pole figure format
        elif i==nskip_calc: #He had two formats: raw and calculated.
            print 'Calculated pole figure format'
            # axes: (phi, khi)
            data = data.T #(khi, phi)
            f = open(filename)
            temp = f.readlines()
            khi = map(float, temp[nskip_calc - 1].split()[1:])
            khi = np.array(khi)
            khi = khi * np.pi / 180.
            phi = data[0]
            phi = phi * np.pi / 180.
            data = data[1:]
            data = data.T #(phi, khi)

            ## dp and dk
            dp = phi[1] - phi[0]
            dk = khi[1] - khi[0]

            ## shift phi back-modification
            phi, data = shiftphi(data, dp, phi)
            ##

            ## check if there is phi=0, 360 at the same time
            phi, data = duplicate0_360(phi, data)
            ##
            pass

    elif mode=='bruker':
        ## make full use of existing uxd.py script
        ## --> normalization is missing...
        ## This must be completed!
        import uxd
        print 'You are now in the bruker mode under epfformat'
        print 'given file name is %s'%filename
        myuxd = uxd.pf(filename=filename, mode='pf')
        if len(myuxd.polefigures)>1:
            print 'multiple pole figures are found'
            raw_input()
        for i in range(len(myuxd.polefigures)):
            pf = myuxd.polefigures[i]
            pf = pfnorm(pf) ##normalize
            pass
        pass

    elif mode=='epf':
        """
        ready made popLA epf format reader
        consider the possibility of multiple number of polefigure
        ## phi must be 0~355, khi must be 0~90 with 5 as an ang resolution
        """
        print 'You are now reading %s'%filename
        blocks = open(filename, 'rU').read().split('(')[1:]

        npf = len(blocks)
        datasets = []
        max_khi = []
        for i in range(len(blocks)):
            #d = __epffiletrimmer__(blocks[i]) #only for popLA epf format
            ## epffiletrimmer should return a block of intensity grid map
            ## as well as maximum khi.
            ## --> modification (2011 NOV 8)
            hkl = blocks[i][0:3] #hkl
            hkl = map(int, [hkl[0],hkl[1],hkl[2]])
            if blocks[i][3]!=')':
                print 'Caution: unexpected hkl labeling format'
                pass
            d, mxk = __epffiletrimmer__(blocks[i]) #only for popLA epf format
            datasets.append(d)
            max_khi.append(mxk)
            pass
        print "number of pole figures:", len(datasets)

        ## presumption of using 19x72 should be deprecated...
        data = np.zeros((len(datasets), 19, 72)) #npf, nphi, nkhi
        dum = data.copy()
        for i in range(len(datasets)):
            for j in range(len(datasets[i])):
                for k in range(len(datasets[i][j])):
                    data[i,j,k] = datasets[i][j][k]
        ## Swap the axis
        data = data.swapaxes(1,2) # E.g. 72 x 19 than 19 x 72

        ## psuedo-normalization
        for i in range(len(data)):
            data[i] = pfnorm(data[i].copy())
        ## --

        return data, max_khi, hkl
    else: raise IOError, 'Unexpected mode is given'
    return data

def __epffiletrimmer__(block):
    """
    EPF (experimental pole figure) trimmer for a single block of data
    EPF: The popLA format
    """
    pfdata = block.split('\n')
    dkhi = float(pfdata[0][5:9])   #khi incremental angle
    fkhi = float(pfdata[0][9:14])  #khi end angle
    dphi = float(pfdata[0][14:19]) #phi incremental angle
    fphi = float(pfdata[0][19:24]) - dphi #written as 360 but it is only upto 355

    ## now it is enfored to have 5 deg resolution, which
    ## results in 19x72
    data = np.zeros((19,72)) ## intensities along khi and phi axes coordinate
    ## the rest information is neglected from the header

    pfdata = pfdata[1:]
    lines = []
    for i in range(len(pfdata)):
        if len(pfdata[i])==72:
            lines.append(pfdata[i][:])
            pass
        elif len(pfdata[i])==73:
            lines.append(pfdata[i][1:])
            pass
        pass
    pfdata = lines     #       (maxkhi/dkhi + 1) * 4
    if len(pfdata)!=76:# 76 =  ((90 / 5) + 1) * 4 (4lines for one khi level)
        print 'len(pfdata) =', len(pfdata)
        print pfdata
        raise IOError, 'Unexpected pfdata format or type'

    if True:
        for j in range(19): #90 / 5 + 1 #number of khi threads
            kh = pfdata[j*4: (j+1)*4] #along the same kh level
            khline = ''
            for k in range(4):
                khline = khline + kh[k]
                pass
            khline # all data in a thread of string every 4digits corresponds to a datum
            kp = 0 #each point along phi
            for k in range(72): #72 = 288 / 4 (4digits) : 0~355 5.0
                datum = khline[k*4:(k+1)*4]
                data[j,k] = int(datum)
                pass
            pass
        pass

    # block of mesh grid pole figure map and khi final angle.
    return data, fkhi


def shiftphi(data, dp, phi):
    """ shifted phi modification """
    tiny = 0.0000001
    if abs(phi[0] - 0.)> tiny and abs((dp-phi[0])-dp/2.)< tiny:
        dum = data.copy()
        for i in range(len(data)): #phi axis
            for j in range(len(data[i])): #khi axis
                dum[i, j] = (data[i-1, j] + data[i, j]) / 2.
        data = dum.copy()
    return phi, data

def duplicate0_360(phi, data):
    """
    If both phi=0 and phi=360 exist,
    delete phi=360
    """
    tiny = 0.00000000001
    isallsame = True
    for i in range(len(data[0])):
        print data[0,i]
        print data[-1,i]
    if any(abs(data[0,i]-data[-1,i])>tiny for i in range(len(data[0]))):
        print data[0,i]
        print data[-1,i]
        isallsame = False

    isphioverlapped=False
    print abs(phi[0] - phi[-1] - 2. * np.pi)
    if abs(phi[0] - phi[-1] + 2. * np.pi)<tiny:
        isphioverlapped=True

    if isphioverlapped and isallsame:
        newdata = np.zeros((data.shape[0]-1, data.shape[1]))
        for i in range(len(newdata)):
            for j in range(len(newdata[i])):
                newdata[i,j] = data[i,j]
        return phi[:-1], newdata
    elif isphioverlapped and isallsame==False:
        print "conflict results!"
        print "phi=0, and 360 turned out to coexist but"
        print "Values along these two axes are not equivalent"
        raise IOError
    elif not(isphioverlapped) and isallsame:
        print "Conflict results!"
        print "Phi is overlapped but phi[0] and phi[-1] is same"
        raise IOError
    else:
        print "No duplicated axis is found"
        return phi, data

def cart2polar(x,y):
    """ cartesian to polar coordinate """
    r = np.sqrt(x**2+y**2)
    theta = np.arctan2(y,x)
    return r, theta

def circle(center=[0,0], r=1.):
    """
    Draw a circle around the given center point with a radius of r(as given)
    The default settings are as below.
    --------
    Arugment
    center = [0,0]
    r = 1.
    """
    ang = np.linspace(0,2*np.pi,1000)
    #unit circle * radius
    x = np.cos(ang)*r
    y = np.sin(ang)*r
    #circle transloation
    x = x + center[0]
    y = y + center[0]
    return x,y

def __isunique__(a, b):
    """
    Is a(3) in b(m, 3)
    """
    for i in range(len(b)):
        ## is a//b[i]? either + or -
        diff0 = abs(b[i][0] - a[0]) + abs(b[i][1] - a[1]) + abs(b[i][2] - a[2])
        diff1 = abs(b[i][0] + a[0]) + abs(b[i][1] + a[1]) + abs(b[i][2] + a[2])
        if diff0 < 0.1**4 or diff1 < 0.1**4:
            return False
    return True

def __circle__(center=[0,0], r=1.):
    """
    Draw a circle around the given center point with a radius of r(as given)
    The default settings are as below.

    --------
    Arugment
    center = [0,0]
    r = 1.
    """
    ang = np.linspace(0,2*np.pi,1000)
    #unit circle * radius
    x = np.cos(ang)*r
    y = np.sin(ang)*r
    #circle transloation
    x = x + center[0]
    y = y + center[0]
    return x,y

def projection(pole=None, agrain=None):
    """
    Projects the a pole to the projection plane.
    (default is stereography projection)

    pole = [1,1,1] or [1,1,0] something like this.

    ---------
    Arguments
    pole = None
    agrain = [ph1, phi, phi2, vf]
    """
    #normalization of the miller indices
    norm = np.sqrt(pole[0]**2 + pole[1]**2 + pole[2]**2)
    pole = pole / norm

    a = pole[0]; b = pole[1]; c = pole[2]
    ###  mid-plane projection (z=0)
    if c==1:
        # pole[0] = 0; pole[1]=0; pole[2] = 1
        X=0; Y=0
    else:
        X = a/(c-1)
        Y = b/(c-1)
    return X,Y

def invproj(x=None,y=None):
    """
    Converts projected point to the pole
    """
    X = 2*x/(1+x**2+y**2)
    Y = 2*y/(1+x**2+y**2)
    Z = (-1+x**2+y**2)/(1+x**2+y**2)
    return np.array([X,Y,Z])

def vect2sphe(pole):
    """
    cartensian vector to spherical cooridnate system
    """
    seca = math.sqrt(pole[0]**2 + pole[1]**2)
    if seca < 0.1**6 :
        x = 0.0
        y = cos(seca)
    else:
        x = math.atan2(pole[1], pole[0])
        y = pole[2]
    return x, y

def cart2sph(pole):
    """
    argument: pole
    Returns phi and theta in such order indeed.

    """
    pole = np.array(pole)
    r = np.sqrt((pole**2).sum())
    theta = math.acos(pole[2]/r)
    phi = math.atan2(pole[1],pole[0])
    return np.array([phi, theta])


def agr2pol(agrain=None, miller=None, proj=None):
    """
    -- for pole figure projection (proj='pf')
    For the given grain, crystallographic direction
    is mapped onto the sample axes by doing tensor
    product of the pole vector with rotation matrix.

    -- for inverse pole figure projection (proj='ipf')
    For the given grain, its principal axis mapped
    onto the crystal axes.

    ---------
    Arguments
    agrain = None
    miller = None
    proj   = None
    """
    if proj==None: print "argument proj should be given"; raise IOError
    elif proj!='pf' and proj!='ipf':
        print " proj should be either 'pf' or 'ipf'"
        raise IOError
    if type(miller).__name__=='list': miller = np.array(miller)

    # a-matrix between the grain's coordinate and sampl coordinate
    phi1 = agrain[0]; phi = agrain[1]; phi2 = agrain[2] #VPSC convention
    amat = euler(ph=phi1, th=phi, tm=phi2, echo=None) #rotation matrix (sa->ca)

    #normalize miller
    norm = math.sqrt( miller[0]**2 + miller[1]**2 + miller[2]**2)
    miller = miller / norm

    if proj=='pf':
        #Returns the dot product of the transposed A-matrix and miller vector
        "A^{T}_{ij} * V_{j}"
        return np.dot(amat.transpose(), miller) #ca to sa

    elif proj=='ipf':
        #for inverse projection,
        #the sample axis should be one of principal axes.(x or y or z)
        if   miller[0]==1: p_ca=0
        elif miller[1]==1: p_ca=1
        elif miller[2]==1: p_ca=2
        else: print"something is wrong"; raise IOError
        #map the pole in the sample axes into the crystal axes
        "A_{ij} * V_{j}"
        return np.dot(amat, miller) #sa to ca returns the pole in ca
    else:
        print "projection should be pf or ipf"
        raise IOError

def ipfline(center=[0,0],csym='cubic'):
    """
    the boundary line for inverse pole figure
    ---------
    Arguments
    center = [0,0]
    csym = 'cubic'
    """
    xc = []; yc = []
    if csym!='cubic': print "Only Cubic!"; raise IOError
    xc.append( center[0])
    yc.append( center[1])

    for i in np.linspace(0.,1/math.sqrt(3.)):
        yaux = i
        xaux = math.sqrt((1. - yaux**2)/2.)
        zaux = xaux
        t1 = math.sqrt(1. - zaux) / math.sqrt(xaux**2 + yaux**2)
        t2 = t1/math.sqrt(1. + zaux)
        ## equal area
        # xc.append(xaux*t1)
        # yc.append(yaux*t1)
        ## stereo
        xc.append(xaux*t2)
        yc.append(yaux*t2)

    xc.append(center[0])
    yc.append(center[1])
    return np.array([xc,yc])

"""
Sample symmetry application is performed over RVE calculation
refer to the RVE class in cmb.py module.
"""
# def planar_sym(gr, nrot):
#     """
#     multiplication of the given grain considering sample symmetry
#     """
#     gr = np.array(gr)
#     old = gr.copy()
#     new = [ ]

#     for i in range(len(gr)):
#         dang = 360/nrot
#         for n in range(nrot):
#             current_grain = gr[i]
#             new.append([current_grain[0] + n*dang, current_grain[1],
#                         current_grain[2], current_grain[3]/nrot])
#             pass
#         pass

#     return np.array(new)

class polefigure:
    # decides if the given set is in the texture file form or array
    def __init__(self, grains=None, filename=None, csym=None, ngrain=100,
                 cdim=[1.,1.,1.], cang=[90.,90.,90.], ssym=False, epf=None):
        """
        ----------------
        class polefigure
        ----------------
        Makes or accepts grains and saves them into the global
        variable self.gr. As a local method, there is 'def core'
        in which miller indices of poles (pole), are returned for
        a grain (agrain). The crystallographically equivalent poles
        are calculated if symmetry operation is allowed (isym).
        ---------
        Arguments
        grains = None
        filename = None
        csym = 'cubic' or'hexag'  #crystal symmetry
        ngrain = 100
        cdim=[1.,1.,1.]
        cang=[90.,90.,90.]
        ssym=False : Sample symmetry: not even started...
        epf = None : experimental pole figure file
        """
        # The grain aggregte can be given either through a file or #
        # passing an array of them to the class directly.          #
        # either grains or filename                                #
        # if none of them is given, a 500-grains file is generated #
        # and returns its grains to the global gr variable.        #

        if grains==None and filename==None and epf==None:
            print " ****************************** "
            print " Since no argument is passed,"
            print " 100 random grains are created"
            print " ****************************** \n"
            a = re(ngrain=ngrain)
            gr = np.array(a.euler).transpose()
            gr = np.array([gr[1],gr[2],gr[3]]).transpose()
            temp = []
            for i in range(len(gr)):
                temp.append([gr[i][0],gr[i][1],gr[i][2],0.01])
            self.gr = np.array(temp)

        self.epf = epf # global

        if grains!=None:
            self.gr = np.array(grains)
        elif filename!=None:
            self.gr = np.genfromtxt(fname=filename,skiprows=4)
            pass
        elif epf!=None: # None is the default for epf
            """
            experimental pole figures..
             # available format:
                 - UXD
                 - steglich
                 - bruker
                 - epf*
            """
            if type(epf).__name__=='list': self.epf_fn = epf
            elif type(epf).__name__=='str': self.epf_fn = [epf]
            elif epf==True:
                fn = [] # list of file names
                print 'type the experimental pole figure file names'
                print "To finish input, press enter"
                while True:
                    dum = raw_input(">>> ")
                    if len(dum)==0: break
                    fn.append(dum)
                    pass
                self.epf_fn = fn
                pass
            else: raise IOError, 'Unexpected epf type found'

            ## check if the file name is correct ##
            for i in range(len(self.epf_fn)):
                if not(os.path.isfile(self.epf_fn[i])):
                    raise IOError, "Could not find %s"%self.epf_fn[i]
                pass
            ## --------------------------------- ##

            ## POLE FIGURE MODE --------------------------------------
            print "Type the experimental polfe figure mode"
            print "Available options:", #continuation
            print "bruker, steglich, epf (default: %s)"%'epf'
            epf_mode = raw_input(" >>>" )
            if len(epf_mode)==0:
                epf_mode='steglich'
                pass
            ##---------------------------------------------------------

            self.grid = []; self.hkl = []
            ## more than one pf can be included.
            npole_per_file = []
            if epf_mode=='epf': self.max_khi = [] #Available only for epf_mode yet.

            for i in range(len(self.epf_fn)):
                if epf_mode=='epf':
                    data, maxk, hkl = epfformat(
                        mode=epf_mode,
                        filename=self.epf_fn[i]
                        )
                    # one file may include multiple poles
                    for i in range(len(data)):
                        self.grid.append(data[i])
                        self.max_khi.append(maxk[i])
                        self.hkl.append(hkl)
                    npole_per_file.append(len(data)) # of pole per a file

                else:
                    data = epfformat(
                        mode=epf_mode,
                        filename=self.epf_fn[i]
                        )
                    self.grid.append(
                        data
                        )
                    self.hkl.append(None)
            self.grid = np.array(self.grid)
            self.epf_mode=epf_mode

        ## EXPERIMENTAL POLE FIGURE
        ## ---------------------------------------------------------- ##
        ## POLE FIGURES BINNINGED FROM THE POLYCRYSTALLINE AGGREGATES ##

        if epf==None:
            dat = self.gr.transpose()
            phi1 = dat[0]; phi = dat[1]; phi2 = dat[2]

            print 'phi1: %i ~ %i'%(
                int(round(min(dat[0]/90.)))*90, int(
                    round(max(dat[0]/90.)))*90)
            print 'phi:  %i ~ %i'%(
                int(round(min(dat[1]/90.)))*90, int(
                    round(max(dat[1]/90.)))*90)
            print 'phi2: %i ~ %i'%(
                int(round(min(dat[2]/90.)))*90, int(
                    round(max(dat[2]/90.)))*90)
            ph1min, ph1max= int(
                round(min(dat[0]/90.)))*90, int(
                round(max(dat[0]/90.)))*90
            phmin, phmax  = int(
                round(min(dat[1]/90.)))*90, int(
                round(max(dat[1]/90.)))*90
            ph2min, ph2max= int(
                round(min(dat[2]/90.)))*90, int(
                round(max(dat[2]/90.)))*90

            ## symmetric multiplication over self.gr is performed unless ph1max==360
            """
            Sample symmetry application is pending,
            because it is done over rve (refer to cmb.py)
            """
            # nrot = int(round(360./ph1max))
            # if ssym==True:
            #     if nrot==4: self.gr = planar_sym(gr=self.gr, nrot=2)
            #     else: raise IOError, "not ready for other nrot than 4"
            #     pass

            ### environments global variables
            #1 symmetry
            self.csym = csym
            self.ngr = len(self.gr)
            self.cdim = cdim
            self.cang = cang
            pass
        pass

    def epfplot(self, ifig, mode, cmode, levels, rot):
        """
        if self.epf!=None:, pf is redirected to this defunc.

        * pole figure master data: self.grid
          [i]: i-th pole
          [i][j]: i-th pole's j-th phi segment's khi
          [i][j][k]: An intensity (count) at i-th pole's j-th phi sement at the k-th khi
        """
        print 'List of files:'
        for f in self.epf_fn: print '%s '%f
        #print 'dimension of self.grid:', self.grid.shape
        print 'ifig=',ifig

        fact = 2.
        nrow = len(self.grid)

        figs = []
        for i in range(nrow):
            figs.append(plt.figure(ifig+21+i))
            pass

        ## loop over each of pole figures
        for ip in range(len(self.grid)): #upon each of self.eps_fn
            # each pole figure set
            pf = np.zeros((self.grid[ip].shape[0]+1, self.grid[ip].shape[1]))
            for i in range(len(self.grid[ip])):
                for j in range(len(self.grid[ip][i])):
                    #0, 5, 10 ... 350, 355, 360 : phi
                    pf[i,j] = self.grid[ip][i][j]
                    pf[-1,j] = self.grid[ip][0][j]
                    pass
                pass

            nm = len(pf); nn = len(pf[0]) #phi, khi
            dp = 360. / nm; dk = 90. / nn
            phi = np.linspace(0., 2.*np.pi, nm)
            khi = np.linspace(np.pi, np.pi/2., nn)

            #print self.grid.shape
            try: max_khi = self.max_khi[ip] # if popLA
            except:
                #upto 80 degree is convetional to me...
                max_khi = 80.
                pass

            ## trimmed out the unmeasured khi rims ---- ##
            iidk = 0; tiny = 10**-4
            #print max_khi
            #print abs(np.rad2deg(khi[::-1] - np.pi))
            while True:
                if abs(abs(np.rad2deg(khi[::-1][iidk]-np.pi)) - max_khi) < tiny:
                    break
                iidk = iidk + 1
                pass

            if iidk==0: #full khi available
                khi_r = khi[:]
                pass
            else: khi_r = khi[:-iidk]
            ## ---------------------------------------- ##

            khi_r = khi[0:-2]                     #reduced khi_r

            r = np.sin(khi)/(1-np.cos(khi))
            rr = np.sin(khi_r)/(1-np.cos(khi_r)) #reduced rr..

            ##
            ## prescribed rotation
            phi = phi + rot * np.pi/180.
            R, PHI = np.meshgrid(r, phi)    #meshing radius and rotation angle
            RR, PHIR= np.meshgrid(rr,phi)
            PHI = PHI + np.pi/2. # rotation the pole figure up.
            PHIR = PHIR + np.pi/2. # rotation the pole figure up.
            phi = phi + np.pi/2. # rotation the pole figure up.
            x = R * np.cos(PHI); y = R*np.sin(PHI) #convert the polar coord

            each_ax = figs[ip].add_subplot(111, polar=True)
            each_ax.set_axis_off()
            ## ------------------------------------- ##

            ########################################
            ## border line drawing
            ##
            # polar
            r0 = np.ones(100)
            t0 = np.linspace(0., np.pi*2, 100)
            each_ax.plot(t0, r0, color='gray', alpha=0.5)

            max_khi = (90 - max_khi) + 90
            max_khi = max_khi * np.pi/180.
            # _ max khi outline
            r0 = np.sin(max_khi)/(1-np.cos(max_khi))
            r0 = np.ones(100) * r0
            each_ax.plot( t0, r0, 'k')
            ########################################

            # Trim out the unmeasured rims.
            pfr = pf.T.copy()
            if iidk==0: pass
            else: pfr = pfr[0:-iidk]

            cnt_each = each_ax.contour(
                phi.copy(), rr.copy(), pfr,
                linewidth=1.)

            pcmr_each = each_ax.pcolormesh(phi.copy(), rr.copy(), pfr)

            ## add pole indices or pole figure file name
            x0, y0 = 0.4, -1.18
            r0, t0 = cart2polar(x0, y0)
            if self.epf_mode=='epf':
                #if self.hkl[ip]==None:
                hkl = raw_input(
                    "Type the indices delimiter as a space (e.g. 1 1 1)>>>")
                hkl = map(int, hkl.split())
                #else: hkl = self.hkl[ip]
                # (hkl)
                index = '('
                for hk in hkl: index = index+'%i'%hk
                index = index + ')'
                each_ax.text(
                    x=t0, y=r0, s=index, fontsize=8.*fact)
                pass
            else:
                if self.hkl[ip]==None:
                    hkl = raw_input(
                        "Type the indices delimiter"+\
                            " as a space (e.g. 1 1 1)>>>")
                    hkl = map(int, hkl.split())
                    pass
                # (hkl)
                index = '('
                for hk in hkl: index = index + '%i'%hkl[hk]
                index = index + ')'
                each_ax.text(
                    x=t0, y=r0, s='%s %s'%(index, self.epf_fn[ip]),
                    fontsize=8.*fact)
                pass

            tcolors = cnt_each.tcolors
            clev = cnt_each._levels

            for i in range(len(tcolors)):
                cc = tcolors[i][0][0:3]
                if levels==None or ip==len(pole)-1:
                    ## Colored marker
                    x0, y0 = 1.3, 0.8 - i * 0.2
                    r0, t0 = cart2polar(x0,y0)

                    each_ax.plot(t0, r0, marker='o', mfc=cc, ms=7.,
                                 ls='None', mec='black',
                                 markeredgewidth=0.01)

                    ## Contour level
                    x2, y2 = 1.35, 0.8 - i *0.2 - 0.05
                    r2, t2 = cart2polar(x2, y2)
                    each_ax.text(x=t2, y=r2,s='%4.2f'%(clev[i]),
                                 fontsize=4.*fact)

                    pass

                ## RD and TD indication
                x4, y4 = -0.05, 1.05
                r4, t4 = cart2polar(x4, y4)
                each_ax.text(x=t4, y=r4, s='RD', fontsize = 6.*fact)
                x5, y5 = 1.05, 0.
                r5, t5 = cart2polar(x5, y5)
                each_ax.text(x=t5, y=r5, s='TD', fontsize = 6.*fact)
                pass

            # Save individual pole figure
            figs[ip].savefig('figs_%s.pdf'%str(ip).zfill(2))
            figs[ip].savefig('figs_%s.eps'%str(ip).zfill(2))
            figs[ip].clf()

            pass
        return

    def pf_axis(self, pole=[[1,0,0]], ifig=1):
        """
        Plot each pole without crystal symmetry
        """
        color = ['r','b','g','k','gray']
        #marker =['o','x','+','d','.']
        for ip in range(len(pole)):
            cl = color[ip]
            #mk = marker[i]
            for i in range(len(self.gr)):
                tm = self.dotplot(proj='pf', agrain=self.gr[i],
                                  npole=len(pole), ipole=ip+1,
                                  pole=pole[ip], ifig=ifig,
                                  cdim='None', cang=self.cang,
                                  csym=self.csym, mode=None,
                                  color=cl)


    def pf(self, pole=[[1,0,0],[1,1,0],[1,1,1]], mode='contour',
           ifig=1, dm=7.5, dn=7.5, ssym=None, levels=None,
           axes=None, cmode=None,rot=0., nths=None):
        """
        Plots pole figures, either experimental pole figure or
        polycrystalline aggregate.

        Arguments
          pole = [[1,0,0,], [1,1,1], [1,1,0]]
          mode ='dot', 'contour', 'contourf', 'im'
          ifig = 1
          dm = 7.5
          dn = 7.5
          ssym = None ## --> dummy yet (2011-Sept-13)
          levels = None
          axes = None
          cmode = 'gray_r'  #the color modes under matplotlib.cm

        ** Exemplary cmodes:
        'gray_r, gray, pink_r, pink, summer_r, summer, winter_r, winter,
        Blues, Blues_r, Set1_r, Set1 .... '
        """

        ## PF ploting directly from experimental pole figure
        if self.epf!=None:
            # mode is fixed to be 'im'
            self.epfplot(
                ifig=ifig, mode=mode, cmode=cmode, levels=levels, rot=rot
                )
            return
        ## end of experimental polefigure plotting

        ## if it is not an experimental polefigure,
        ## binning of COD on the spherical surface is carried out.
        if self.csym=='hexag' or self.csym=='trigo':
            p__ = pole  ## for later use on pole index indication
            pole_tmp=[]
            for ip in range(len(pole)):
                p_ = pole[ip]
                p = [0,0,0]
                if len(pole[ip])!=4:
                    raw_input('pole must be four digit for hexag')
                    raise IOError
                p[2] = p_[3]
                p[0] = p_[0] - p_[2]
                p[1] = p_[1] - p_[2]
                p[2] = p_[3]
                #print 'pole as below'
                #print p[0:3]
                pole_tmp.append(p[0:3])
            pole = pole_tmp
        pi = math.pi
        ## pole figure
        temp = []
        if mode=='dot':
            for ip in range(len(pole)):
                color = ['r','b','g','k','gray']
                marker =['o','x','+','d','.']
                for i in range(len(self.gr)):
                    if len(self.gr)<5:
                        cl = color[i]
                        mk = marker[i]
                    else:
                        cl='k'
                        mk='.'
                    tm = self.dotplot(proj='pf', agrain=self.gr[i],
                                      npole=len(pole), ipole=ip+1,
                                      pole=pole[ip], ifig=ifig,
                                      cdim=self.cdim, cang=self.cang,
                                      csym=self.csym, mode=None,
                                      color=cl
                                      )
        elif mode=='trace':
            for ip in range(len(pole)):
                for i in range(len(self.gr)):
                    if i==0:
                        color='red'
                        alpha=1.0
                        marker='o'
                    elif i==-1:
                        color='blue'
                        alpha=1.0
                        marker='o'
                    tm = self.dotplot(proj='pf', agrain=self.gr[i],
                                      alpha=alpha,color=color,
                                      marker=marker,
                                      pole=pole[ip], ifig=ifig,
                                      cdim=self.cdim, cang=self.cang,
                                      csym=self.csym, mode=None)
                    color='black'
                    marker='.'
                    alpha=0.5

        elif mode in ['contour', 'contourf']:
            cnts = [] # cnt container
            ## figure setting
            # figure size
            fact = 2. #size factor
            figsize = (len(pole)*2.*fact, 1.*2.*fact)

            ## If axis is passed to the module, plotting is performed on them.
            if axes!=None:
                if len(axes)!=len(pole): raise IOError
                fig = plt.figure(ifig)
                pass
            elif axes==None:
                fig = plt.figure(ifig, figsize=figsize)
                nrow = len(pole)
                pass
            levs = []


            #### on each of the pole ------------------------------ ####
            N = []
            self.pfnode = []
            start = time.time()
            for ip in range(len(pole)):
                # Polar cell and nodes generation.
                f, nodes = self.cells(
                    pole=pole[ip],
                    ifig=None,
                    dm=dm, dn=dn,
                    cdim=self.cdim,
                    cang=self.cang,
                    csym=self.csym,
                    nths=nths)


                N.append(nodes)
                self.pfnode.append(nodes)
                print 'node shape:', nodes.shape

            print "%5.2f seconds elapsed during calling"\
                " self.cells\n"%(time.time() - start)
            del nodes

            ## resolution and pole figure plotting preferences
            nm = (360.0 - 0.)/dm; nn = (180.-90.)/dn
            # theta and phi and stereographic projection of them.
            theta = np.linspace(pi, pi/2., nn+1) #tilting angle
            phi = np.linspace(0.,2.*pi, nm+1)    #rotation angle
            r = np.sin(theta)/(1-np.cos(theta))  #tilting angle to radius
            R, PHI = np.meshgrid(r,phi)          #meshing radius and rotation angle
            PHI = PHI + pi/2. # rotation the pole figure up.
            x = R*np.cos(PHI); y = R*np.sin(PHI) #convert the polar coord-> cartensian

            ## somehow, I did not recognize that there is an intrinsic
            ## polar projection option in subplot. Therefore, at the time
            ## this defunc was made, I did transform them into the
            ## cartensian coordinate system above. Now I am wishing to
            ## make use of it, over excursion of pixel view of the
            ## binning of pole figures. That is mode=='im'.
            ## I do not like this lousy name again.
            ##
            ## On the top of it, matplotlib contours do not have
            ## polar coordinate system. Thus, if one wants to have
            ## contours with lines, this manual transform is the only
            ## way. (2011-Sept)

            for ip in range(len(pole)):
                nodes = N[ip] ## intensity at the current pole (ip)
                if axes!=None: ax=axes[ip]
                elif axes==None:
                    ax=fig.add_subplot(1,nrow,ip+1)
                    pass

                ax.set_frame_on(False)

                # contour plotting
                if mode=='contour':
                    if levels==None:
                        if cmode!=None:
                            cnt = ax.contour(
                                x, y, nodes,
                                cmap=plt.cm.cmap_d[cmode])
                            pass
                        else: cnt = ax.contour(
                            x, y, nodes)
                        pass
                             #cmap=plt.cm.bone)
                    elif levels!=None:
                        if cmode!=None: cnt = ax.contour(
                            x, y, nodes, cmap=plt.cm.cmap_d[cmode])
                        else: cnt = ax.contour(x, y, nodes, levels)#, cmap=plt.cm.bone);pass
                        pass
                    pass
                elif mode=='contourf':
                    if levels==None:
                        if cmode!=None:
                            cnt = ax.contourf(
                                x, y, nodes,
                                cmap=plt.cm.cmap_d[cmode]
                                )
                            pass
                        else: cnt = ax.contourf(x, y, nodes);pass
                    elif levels!=None:
                        if cmode!=None:
                            cnt = ax.contouf(
                                x,y,nodes,
                                cmap=plt.cm.cmap_d[cmode]
                                )
                        else: cnt = ax.contourf(
                            x, y, nodes,
                            levels)#, cmap=plt.cm.bone)
                        pass
                    pass

                cnts.append(cnt)
                clev = cnt._levels
                levs.append(clev)

                # Contour's details.
                ax.set_axis_off();
                ax.set_aspect('equal')
                rx, ry = circle()
                ax.plot(rx, ry, 'k')

                # misc (decoration of the axes, with information)
                tcolors = cnt.tcolors
                for i in range(len(tcolors)):
                    cc = tcolors[i][0][0:3]
                    #if levels==None:
                    if levels==None or ip==len(pole)-1:
                        ## level line
                        ax.plot(
                            [1.28, 1.35],
                            [1. - i * 0.2, 1. - i * 0.2],
                            color=cc)
                        ## level text
                        ax.text(x=1.47, y= 1. - i*0.2 - 0.05,
                                s='%3.2f'%(clev[i]),
                                fontsize=4.*fact)
                        pass
                    ## pole plane
                    if self.csym=='hexag' or self.csym=='trigo':
                        ax.text(x=0.4,
                                y=-1.18, s='(%1i%1i%1i%1i)'%

                                (p__[ip][0],p__[ip][1],
                                 p__[ip][2],p__[ip][3]),
                                fontsize=6.*fact)
                    else:
                        ax.text(x=0.4, y=-1.18, s='(%1i%1i%1i)'%
                                (pole[ip][0],
                                 pole[ip][1],
                                 pole[ip][2]),
                                fontsize=6.*fact)
                    ## RD and TD indication
                    ax.text(x=-0.05, y = 1.05,
                            s='RD', fontsize = 4.*fact)
                    ax.text(x= 1.05, y = 0.,
                            s='TD', fontsize = 4.*fact)
                # Fixes the frame
                ax.set_xlim(-1.2, 1.5); ax.set_ylim(-1.2, 1.5)
            #### on each of the pole
            ####---------------------------
            ####---------------------------
            #for ip in range(len(pole)):
            return cnts # for mode in ['contour', 'contourf']

        elif mode=='im':
            fact = 2. #plot size factor
            figsize=(len(pole)*2.*fact, 1.*2.*fact)

            if axes!=None:
                raise IOError, "'im' mode does not support"\
                    " imposition of axes"
            fig = plt.figure(ifig, figsize=figsize)
            fig.clf() #clear the figure
            nrow = len(pole)
            Zs = []
            start = time.time()
            print'dm, dn:', dm,dn
            for ip in range(len(pole)):
                f, Z = self.cells(
                    pole=pole[ip], ifig=None,
                    dm=dm, dn=dn,
                    cdim=self.cdim, cang=self.cang,
                    csym=self.csym,
                    nths=nths)
                Zs.append(Z)

            print "%5.2f seconds elapsed during"\
                " calling self.cells\n"%(time.time()-start)

            del Z
            ## resolution and pole figure plotting preferences
            nm = (360.0 - 0.)/dm; nn = (180.-90.)/dn
            # theta and phi and stereographic projection of them.
            theta = np.linspace(pi, pi/2., nn+1) #tilting
            phi = np.linspace(0.,2.*pi, nm+1)    #rotation
            r = np.sin(theta)/(1-np.cos(theta))  #radius
            phi = phi + pi/2. # rotate the RD to the north
                              # pole up (for 'im')
            phi = phi + rot # arbitrary rotation (2011 OCT)

            #phi, r = np.meshgrid(phi,r)
            zmax = 0.
            for ip in range(len(pole)):
                if np.max(Zs[ip])>zmax: zmax=np.max(Zs[ip].copy())
                pass

            for ip in range(len(pole)):
                Z = Zs[ip] ## intensity of the current pole
                nodes = Z.copy()
                Z = Z.transpose()

                ## polar coordinate system
                axp = fig.add_subplot(
                    1, nrow, ip+1, polar=True) #polar

                pcm = axp.pcolormesh(
                    phi.copy(), r.copy(), Z.copy(),
                    #color='red',# alpha=0.75
                    )

                axp.set_axis_off()
                axp.set_aspect('equal')

                cnt = axp.contour(
                    phi.copy(), r.copy(), Z.copy(),
                    #cmap = plt.cm.gray_r,
                    color='red',
                    levels=levels)

                clev = cnt._levels
                #rx, ry = circle()
                tcolors = cnt.tcolors
                for i in range(len(tcolors)):
                    cc = tcolors[i][0][0:3]
                    if levels==None or ip==len(pole)-1:
                        x0, y0 = 1.3, 0.8 - i * 0.2
                        r0, t0 = cart2polar(x0,y0)
                        axp.plot(t0, r0, marker='o',
                                 ls='None', mec='black',
                                 mfc=cc, #color=cc
                                 ms=7./len(pole),
                                 markeredgewidth=0.01/len(pole)
                                 )
                        x2, y2 = 1.40, 0.8 - i *0.2 - 0.05
                        r2, t2 = cart2polar(x2, y2)
                        axp.text(x=t2, y=r2,
                                 s='%4.2f'%(clev[i]),
                                 fontsize=6.*fact/len(pole)
                                 )
                        pass
                    # pole figure indices
                    x3, y3 = 0.4, -1.18
                    r3, t3 = cart2polar(x3, y3)
                    if self.csym=='hexga' or self.csym=='trigo':
                        axp.text(
                            x=t3, y=r3,
                            s='(%1i%1i%1i%1i)'%
                            (
                                p__[ip][0],p__[ip][1],
                                p__[ip][2],p__[ip][3]
                                ) ,
                            fontsize=8.*fact/len(pole)
                            )
                        pass
                    else:
                        axp.text(
                            x=t3, y=r3,
                            s='(%1i%1i%1i)'%
                            (
                                pole[ip][0],
                                pole[ip][1],
                                pole[ip][2]
                                ),
                            fontsize=8.*fact/len(pole)
                            )
                        pass
                    ## RD and TD indication
                    x4, y4 = -0.05, 1.05
                    r4, t4 = cart2polar(x4, y4)
                    axp.text(x=t4, y=r4,
                             s='RD',
                             fontsize = 6.*fact/len(pole))
                    x5, y5 = 1.05, 0.
                    r5, t5 = cart2polar(x5, y5)
                    axp.text(x=t5, y=r5,
                             s='TD',
                             fontsize = 6.*fact/len(pole))
                    axp.set_axis_off()
                    pass
                pass
            ## save figures
            fig.savefig('pcm.pdf')
            fig.savefig('pcm.eps')
            for i in range(len(fig.axes)):
                extent = fig.axes[i].get_window_extent()
                extent = extent.transformed(fig.dpi_scale_trans.inverted())
                fig.savefig(
                    'each_axpcm_%s.pdf'%str(i).zfill(2),
                    bbox_inches=extent.expanded(1.1,1.1)
                    )
                fig.savefig(
                    'each_axpcm_%s.eps'%str(i).zfill(2),
                    bbox_inches=extent.expanded(1.1,1.1)
                    )
                pass
            ##
            print "A figure's been saved to pcm.pdf and .eps"
            print 'each_axpcm_00.pdf and .eps'

            fig.clf()
            #return phi, r, Z
        pass

    def ipf(self,pole=None,color='k', ifig=4):
        """
        Given the pole plot the dot inverse pole figure.
        **The contour version of the inverse pole
        figure is not prepared yet.

        ---------
        Arguments
        ---------
        pole  = [1,0,0]
        color = 'k'
        ifig  = 1
        """
        if pole==None:
            print 'miller index of the pole should be given'
            raise IOError
        temp = []
        for i in range(len(self.gr)):
            tm = self.dotplot(proj='ipf',agrain=self.gr[i],
                              pole=pole, ifig=ifig, color=color)
            temp.append(tm)
            #self.dotplot(proj='ipf',
            #agrain=self.gr[i], pole=[0,1,0], ifig=5)
            #self.dotplot(proj='ipf',agrain=self.gr[i],
            #pole=[0,0,1], ifig=6)
        return temp

    def core(self, pole=None, proj='pf', csym=None, ssym=None,
             agrain=None, isym=True, cdim=[1,1,1],
             cang=[90.,90.,90.],
             equivp=None):
        """
        --------------------------------
        The core of the polefigure class
        --------------------------------

        It is the engine for plotting regular and inverse
        pole figures in that it generates projected
        cartensian coordinate of the 3D vector poles
        onto the pole figure sphere. One can directly plot
        pole figure on a xy plane.

        Provided the miller indices,
        1. Calculates the crystallographically equivalent poles
        2. Maps the miller indices for a grain
           into a vector in sample axes
        3. Returns mapped poles as raw and as projected (xy)

        ---------
        Arguments
          pole = None
          proj = 'pf'
          csym = 'cubic', 'hexag'
          agrain = None
          isym = True
          cdim = [ 1., 1., 1.]
          cang = [90.,90.,90.]
          equivp = None  #crystallographically equivalent pole
        """
        xy = []; POLE = [];
        if csym!='cubic' and csym!='hexag' and csym!='None':
            print "Other symmetries than cubic or hexag nor 'None' "\
                "is not prepared yet"
            raise IOError
        if proj!='pf' and proj!='ipf':
            print "Other modes of projection than pf and "\
                "ipf is not prepared yet"
            raise IOError
        if agrain==None:
            print "A grains must be given to the method"
            raise IOError
        if pole==None:
            print "Pole must be given to core"
            raise IOError
        if type(pole).__name__=='list': pole = np.array(pole)
        elif type(pole).__name__=='ndarray': pass
        else: raise IOError, 'Unexpected type of the pole argument'
        temp = pole.copy()
        del pole; pole = temp

        ##### calculates crystallographically equivalent poles #####

        ## pole figure projection
        if proj=='pf':
            if isym==True:
                npeq = equivp
            else: npeq = [pole]
            nit = len(npeq)
            nppp = []
            for i in range(nit):
                nppp.append(npeq[i]); nppp.append(npeq[i]*-1)
            npeq = np.array(nppp)
            for i in range(len(npeq)):
                for j in range(len(npeq[i])):
                    if abs(npeq[i,j])<0.1**8:
                        npeq[i,j] = 0.

        ## inverse pole figure
        elif proj=='ipf':
            if abs(pole[0]**2+pole[1]**2+pole[2]**2-1)>0.1**6:
                print "The pole must be one of principal axes of"\
                    " the sample"
                print "It should be among [1,0,0], [0,1,0], [0,0,1]"
                print "current pole is as below\n", pole
                raw_input()
                raise IOError
            npeq = [pole] ## Note it is sample axis vector!
            #equivalent pole calculatation is deffered to
            #next for in block
        ## unexpected proj argument
        else: print "It should be either pf or ipf"; raise IOError

        for ip in range(len(npeq)):
            ## 'pf':  converts ca pole to sa pole
            ## 'ipf': converts sa pole to ca pole
            p = agr2pol(agrain=agrain, miller=npeq[ip], proj=proj)
            if proj=='pf':
                ## if a pole is toward the north pole of the unit circle,
                #if p[2]>0: pass
                #else:
                temp = projection(pole=p)
                POLE.append(p)
                xy.append(temp)
            elif proj=='ipf':
                ## calculates equivalent p by
                ## applying symmetry operations
                if isym==True:
                    """
                    Sample axis is now referred to crystal
                    coordinate system. That vector has to
                    be mutiplicated by symmetry operations.
                    """
                    npoles = self.__equiv__(
                        miller=p, csym=csym,
                        cdim=cdim, cang=cang)
                    temp = []

                    for i in range(len(npoles)):
                        temp.append(npoles[i])
                        temp.append(npoles[i]*-1)
                        pass

                    temp = np.array(temp)
                    npoles = temp

                else: npoles=[p]
                for npp in range(len(npoles)):
                    prj_xy = projection(pole=npoles[npp])
                    xy.append(prj_xy)
                    POLE.append(npoles[npp])
                    pass
                pass # if over 'pf' or 'ipf'

            pass # End of for loop over ipf
        return xy, POLE

    def cells(self,pole=[1,0,0], ifig=None, dm=15., dn = 15.,
              csym=None, cang=[90.,90.,90.], cdim=[1.,1.,1.],
              nths=None):
        """
        Creates cells whose resolutioin is mgrid * ngrid.
        Given the delta m and delta n (dm, dn), each pole's
        weight is assigned to a cell which braces it.
        Plots the cell's weight and returns the cell in array.

        ---------
        Arguments

        pole = [1,0,0]
        ifig = None
        dm = 7.5.
        dn = 7.5.
        csym = None
        cang = [90.,90.,90.]
        cdim = [1.,1.,1.,]
        nths = number of line to be smooths in the central region
              of polefigures (default: None)
        """
        ## Frequently used local functions or methods
        pi   = math.pi;  npi  = np.pi
        cos  = math.cos; sin  = math.sin
        ncos = np.cos;   nsin = np.sin

        pol = [] # a group of pole points.

        #start = time.time()

        ## npeq calculations: equivalent x-tal poles
        npeq = self.__equiv__(
            miller=pole, csym=csym,
            cdim=cdim, cang=cang)
        for i in range(len(self.gr)):
            xy, p = self.core(
                pole=pole, proj='pf', csym=csym,
                agrain=self.gr[i], cang=cang,
                cdim=cdim, equivp=npeq)
            # p is n equivalent poles
            p = np.array(p)
            for j in range(len(p)):
                if p[j][2] <0: p[j] = p[j] * - 1  #make it postive.
                ## phi, and theta in radian
                x, y = cart2sph(p[j])
                pol.append([x, y, self.gr[i][3]]) #phi, theta, intensity

        #print '%5.3f seconds elapsed during calling self.core\n'%
        #(time.time() - start)

        ## weight normalization: total sum enforced to be 1.
        pol = np.array(pol)
        pol = pol.transpose()
        wgts = pol[2].sum()
        pol[2] = pol[2]/wgts
        pol = pol.transpose()
        del wgts
        ##

        # grid number along phi   axis (x) is referred to as m (rotating)
        # grid number along theta axis (y) is referred to as n (tilting)
        # m x n grid
        # m and n's range
        # m is in [-pi.,pi], (rotating)
        # n is in [0, pi/2.] (tilting)

        mgrid = 360./dm; ngrid = 90./dn # Number of grid points

        # phi_angle = np.arange(-pi, pi)/dm + dm/2
        # theta_angle = np.arange( 0, pi/2.)/dn + dn/2

        if np.mod(360., dm)!=0: raise IOError, \
                'dm should be a devisor of 360.'
        if np.mod(90., dn)!=0: raise IOError, \
                'dn should be a devisor of 90.', dn

        f = np.zeros(map(int,[mgrid, ngrid])) # central
        #  The mesh f is on the central point, 'x', of each grid as
        # depicted below:
        #
        # 90|---|---|
        #   .
        #   .
        #   .

        #   |---|---|
        #   | x | x |
        # 10|---|---|
        #   | x | x |   . . .
        #  5|---|---|
        #   | x | x |
        #  0|---|---|
        #   0   5  10 ....  360

        nodes = np.zeros((mgrid+1,ngrid+1))
        #  The mesh nodes is on the central point, 'o', of each grid as
        # depicted below:
        #
        # 90o---o---o
        #   |   |   |
        #   .
        #   .
        #   .

        #   o---o---o
        #   |   |   |
        # 10o---o---o
        #   |   |   |   . . .
        #  5o---o---o
        #   |   |   |
        #  0o---o---o
        #   0   5  10 ....  360


        for i in range(len(pol)):
            phi = pol[i][0]; theta = pol[i][1]
            mi = int((phi+pi)/(dm*pi/180.)-0.1**6) #subtract tiny
            ni = int(theta/(dn*pi/180.)   -0.1**6) #subtract tiny
            if mi<0 or ni<0 :
                raw_input('Negative index'); raise IOError
            # print 'phi, theta', phi/np.pi*180., theta/np.pi*180.
            # print 'mi, ni ', mi, ni
            elif mi>mgrid or ni>ngrid:
                raw_input("Unexpected Error")
                raise IOError
            try: f[mi,ni] = f[mi,ni] + pol[i][2] #add pole intensity
                 #accumulate the intensity
            except:
                print phi*180./np.pi, theta*180./np.pi
                raw_input('Error  raised >> ')
                raise IOError, "Something wrong in the index for f[mi,ni]"

        ## ----------------------------------------
        """
        Symmetrization over the cell (f)
        """
        #1. Symmetrization abount the x-axis

        #2. Symmetrization about the y-axis

        ## ----------------------------------------

        # ismooth = False
        # if ismooth==True:
        #     fpole = 0.
        #     for m in range(int(mgrid)):
        #         fpole = fpole + f[m,0]
        #     #fnorm = (1.-cos(7.5 * pi/180. ))*2.*pi
        #     fnorm = 1.*pi
        #     fpole = fpole /fnorm
        #     for m in range(mgrid):
        #         f[m,0] = fpole
        # else:  pass

        z = np.zeros((ngrid+1,))
        deltz = (pi/2.)/float(ngrid)
        for i in range(int(ngrid)+1): z[i] = deltz*i

        deltx = 2.*pi/mgrid
        for m in range(int(mgrid)):
            for n in range(int(ngrid)):
                fnorm = (cos(z[n]) - cos(z[n+1])) * deltx/(2*pi)
                #print 'fnorm =', fnorm;raw_input()
                f[m,n] = f[m,n] / fnorm

        if ifig!=None:
            fig = plt.figure(ifig)
            ax = fig.add_subplot(111)
            ax.plot(f)
            ax.set_ylim(0.,)


        ## Assigns the same intensity for the first ring
        # South-pole : along  n=0 axis ---
        # f0 = 0.
        # for m in range(len(f)): f0 = f0 + f[m,0]
        # f0avg = f0 / len(f)
        # for m in range(len(f)+1): nodes[m,0] = f0avg

        # ## ------------------------------

        # # # along n=1 axis --------------
        # f1 = 0.
        # for m in range(len(f)):
        #     f1 = f1 + f[m,1]
        # f1avg = f1 / len(f)
        # for m in range(len(f)+1): nodes[m,1] = f1avg
        # ## ------------------------------

        # # # along n=2 axis --------------
        # f2 = 0.
        # for m in range(len(f)):
        #     f2 = f2 + f[m,2]
        # f2avg = f2 / len(f)
        # for m in range(len(f)+1): nodes[m,2] = f2avg
        # ## ------------------------------



        # Node's phi(rot) and theta(tilt)
        ph = np.linspace(-pi, pi   , mgrid+1)
        th = np.linspace( 0., pi/2., ngrid+1)

        self.ph = ph
        self.th = th

        regions = np.zeros((mgrid+1,ngrid+1,4,2)) # deprecated??
        for m in range(int(mgrid+1)):
            for nn in range(int(ngrid)):
                n = nn + 1
                nodes[m,n] # <--- currently interesting node
                ## find the nodes' phi and theta
                ## that's phi[m] and theta[n]

                reg_ph = np.zeros((2,)); reg_th = np.zeros((2,))

                reg_ph[0] = ph[m] - (dm*pi/180.)/2.
                reg_ph[1] = ph[m] + (dm*pi/180.)/2.
                reg_th[0] = th[n] - (dn*pi/180.)/2.
                reg_th[1] = th[n] + (dn*pi/180.)/2.

                reg = []
                reg.append([reg_ph[0], reg_th[0]])
                reg.append([reg_ph[0], reg_th[1]])
                reg.append([reg_ph[1], reg_th[0]])
                reg.append([reg_ph[1], reg_th[1]])

                for i in range(len(reg)):
                    p = reg[i][0]
                    t = reg[i][1]
                    if p>pi: p = p - 2.*pi
                    elif p<-pi: p = p + 2.*pi
                    elif p==pi or p==-pi:
                        print 'Unexpected..'; raise IOError
                        pass
                    if t>pi/2.:
                        t = t - (dn*pi/180.)
                        p = p + pi
                        if   p> pi: p = p - 2.*pi
                        elif p<-pi: p = p + 2.*pi
                        pass
                    reg[i][0] = p
                    reg[i][1] = t

                ## for each node, find the 4 adjacent
                ## regions' intensites
                ## and average them out and assign to the node.
                inten = 0.
                for i in range(4):
                    p = reg[i][0]
                    t = reg[i][1]
                    mi = int((p+pi)/(dm*pi/180.)-0.1**6)
                    ni = int(t/(dn*pi/180.) -0.1**6)
                    if mi<0 or ni<0 :
                        raw_input('Negative index')
                        raise IOError
                    inten = inten + f[mi,ni]
                nodes[m,n] = inten/4.


        ## along n=i-th axis
        #nths = 2
        # if nths!=None:
        #     iths = np.arange(nths)
        #     for ith in iths:
        #         n_i = 0.
        #         for m in range(len(nodes)):
        #             n_i = n_i + nodes[m,ith]
        #         niavg = n_i / float(len(nodes))
        #         print 'average:', niavg
        #         for m in range(len(nodes)): nodes[m,ith] = niavg


        return f, nodes

    def dotplot(self, pole=None, ifig=None, npole=1,
                ipole=1,
                alpha=1.0, color='k', marker='.',
                proj='ipf', csym='cubic',
                ssym='tric', agrain=None,
                isym=True, irdc=True,
                cdim=[1.,1.,1.],
                cang=[90.,90.,90.], mode='trace' # or None
                ):
        """
        Plots the individual poles. No contour.

        Arguments:
        pole = None
        ifig = None
        proj = 'pf' or 'ipf'
        csym = 'cubic' Crystal symmetry
        ssym = 'tric'  Sample symmetry
        isym = True:   Symmetry operation to the pole
        irdc = True:   Reduction of inverse pole figure region
        """
        if pole==None:
            print "Pole must be given"
            raise IOError

        ## If a grain is assigned,
        ## Operations are carried on it, instead of on
        ## Global grains.
        if color==None:
            print 'Default color is black'; color='k'
        if agrain!=None:
            "* a grains' euler angle is given"
            "* dot plotting is performed on a grain"
            agr = np.array([agrain]) #p1,p,p2, inten
        else:
            print 'Agrain must be asigned'
            raise IOError
        #print "------------------------------------"
        XY = []; polz = []
        ## Draws the boundary of inverse and normal pole figures
        if ifig==None: pass
        else:
            fact = 3.
            figsize=(npole*fact, 1*fact)
            fig = plt.figure(ifig, figsize=figsize)
            ax = fig.add_subplot(1,npole,ipole)
            ax.set_axis_off(); ax.set_aspect('equal')
            ax.text(x=-0.05, y= 1.05, s='RD', fontsize=4.*fact)
            ax.text(x= 1.05, y= 0.  , s='TD', fontsize=4.*fact)
            ax.text(x=0.5, y=-1.05, s='(%i%i%i)'%
                    (pole[0],pole[1],pole[2]), fontsize=3. * fact)
            if proj=='pf':
                rx, ry = __circle__(center=[0,0], r=1)
                ax.plot(rx, ry, color='grey')
                ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2)
            elif proj=='ipf':
                # in ipfline module, projection can be
                # changed into equal area type
                cxy = ipfline(center=[0,0], csym=csym)
                ax.plot(cxy[0], cxy[1], color='grey', alpha=0.1)
                pass

        ## add poles
        for i in range(len(agr)):
            ### crystallographically equivalent poles are calculated.
            agr[i][0] = agr[i][0] - 90.
            npeq = self.__equiv__(
                miller=pole, csym=csym, cdim=cdim, cang=cang)
            ### -----------------------------------------------------

            xy, POLE = self.core(
                pole=pole, proj=proj, csym=csym,
                agrain=agr[i], isym=isym, cdim=cdim, cang=cang,
                equivp = npeq)
            ##### POLE FIGURE PROJECTIN #####
            if proj=='pf':
                for j in range(len(xy)):
                    if POLE[j][2]<=0: #north pole is [0,0,1]
                        XY.append(
                            [xy[j][0], xy[j][1], agr[i][3]]
                            ) #x,y, intensity
                    else: pass

            ##### INVERSE POLE FIGURE PROJECTION #####
            elif proj=='ipf':
                for j in range(len(xy)):
                    r = np.sqrt(xy[j][0]**2 + xy[j][1]**2)
                    if r>1.0000000000: pass
                    else:
                        ### reduced region filter
                        ### must be applied here.
                        #y must be positive.
                        tiny = 0.
                        phi = math.atan2(xy[j][1],xy[j][0])
                        phi = phi * 180.0 / math.pi
                        if phi > 0. - tiny and phi < 45.0 +tiny:
                            ## another filter
                            ## 2nd and 3rd triangles
                            #a = POLE[j][0]; c = POLE[j][2]
                            a,b,c = invproj(x=xy[j][0], y=xy[j][1])
                            #print 'a,b,c= ',a,b,c
                            #print 'atan2(a,c)',
                            #math.atan2(a,-c)*180./math.pi
                            if math.atan2(
                                a,-c)*180./math.pi < 45.0 + tiny:
                                XY.append(
                                    [xy[j][0], xy[j][1],
                                     agr[i][3]]
                                    ) #x, y, intensity of the grain
                        else: pass
        if ifig!=None:
            # plotting --------------------
            try:
                XY = np.array(XY)
                xxx = XY.copy().transpose()[0]
                yyy = XY.copy().transpose()[1]
                ax.plot(xxx,yyy,ls='None', marker=marker,color=color,
                        ms=10)
                if proj=='pf':
                    ax.set_xlim(-1.1,1.1)
                    ax.set_ylim(-1.1,1.1)
                elif proj=='ipf':
                    ax.set_xlim(-0.02, 0.5)
                    ax.set_ylim(-0.02,0.5)
            except:
                if len(XY)==0: print 'Empty XY is returned'
                else: print "Unexpected Error"; raw_input()
            #-------------------------------
        return np.array(XY)

    def __equiv__(self, miller=None, csym=None,
                  cdim=[1.,1.,1.], cang=[90.,90.,90.]):
        """
        Provided the miller indices,
        Crystallographically equivalent and only unique
        vectors are returned.

        ---------
        Arguments
        ---------

        miller = None  , e.g. [1,1,1]
        csym ='cubic'
        """
        start = time.time()
        from sym import cv
        from sym import cubic, hexag
        #from sym_cy import cubic, hexag
        import sym    #python compiled
        #import sym_cy #cython compiled
        #from sym.py cvec, cubic, and hexgonal modules are brought in
        if miller==None:
            print "Miller index should be given"
            raise IOError
        vect = np.array(miller)
        norm = 0.; sneq = []
        temp = vect.copy()
        #norm = vect[0]**2 + vect[1]**2 + vect[2]**2
        #norm = np.sqrt(norm)
        #vect = vect/ norm
        #print 'elapsed time before v calculation: %8.6f'%
        #(time.time()-start)

        ##---------------------------------
        ##---------------------------------
        #start = time.time()
        if csym=='cubic':
            #H = sym_cy.cubic(1)  #cython compiled
            H = sym.cubic()  #operators
            for i in range(len(H)):
                sneq.append(np.dot(H[i], vect))
                pass
            pass
        elif csym=='hexag':
            #H = sym_cy.hexag(1) #cython compiled
            H = sym.hexag() #operators
            v = cv(pole=vect, cdim=cdim, cang=cang)
            for i in range(len(H)):
                sneq.append(np.dot(H[i], v))
                pass
            pass
        elif csym=='None':
            #H = [np.identity(3)]
            sneq = [vect]
        else:
            print 'Given symmetry, %s is not prepared'%csym
            raw_input('Enter to raise an error and quits the job');
            raise IOError

        #print 'elapsed time during v calculation: %8.6f'%
        #(time.time()-start)
        #####-------------------------------
        #####--------------------------------

        start = time.time()
        stacked = [] #empty unique vectors
                # is cH in the already existing stacked list?
                # yes: pass
                # no : add

        ## filtering the sneq under whether or not it is unique
        for i in range(len(sneq)):
            cH = sneq[i].copy()  #current vector
            if __isunique__(a=cH, b=stacked):
                stacked.append(cH)
            else: pass
            pass

        ## if v[2] is minus, mutiply minus sign to the vector.
        for i in range(len(stacked)):
            if stacked[i][2]<0:
                stacked[i] = stacked[i]*-1
                pass
        #print 'elapsed time during the rest: %8.6f'%
        #(time.time()-start)
        return np.array(stacked)
    pass # end of class polefigure

### Excutes the module in the command line with arguments and options
def main(filename, pfmode, gr, csym):
    """
    plots the (100),(110) and (111) pole figures of cubic crystal.
    (to be including hexagonal)

    Arugments
      filename = '00500.cmb' for instance.
      csym : crystal symmetry 'cubic', 'hexag'
      pfmode ='contour','contourf', 'dot'
      """
    if gr!=None: a = polefigure(grains=gr, csym=csym)
    else: a = polefigure(filename=filename, csym=csym)
    a.pf(mode=pfmode)
    pass


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import getopt, sys


    ## arguments ------------------------------- ##
    try: opts, args = getopt.getopt(sys.argv[1:],
                                    'm:i:o:c:s')#, 'output=', 'csym='])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
        pass
    ## ----------------------------------------- ##

    ## default options
    ishow = False
    mode = 'contourf'
    csym = 'cubic'
    outputfile = 'temp.pdf'

    for o, a in opts:
        if o in ('-i'): inputfile = a
        elif o in ('-o'): outputfile = a
        elif o in ('-c'): csym = a
        elif o in ('-m'): mode = a
        elif o in ('-s'): ishow = True
        else: assert False, 'unhandled option'
        pass

    main(filename=inputfile, csym=csym, pfmode=mode)
    plt.gcf().savefig(outputfile)
    if ishow==True: plt.show()
    pass
