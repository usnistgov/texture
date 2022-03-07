"""
This module started as a hobby project while learning Python back in 2011

Features:
  crystal symmetries : cubic and hexagonal
    This module can be easily extened for other crystal symmetries
    but have not been thoroughly checked yet.

Pole figure and Inverse pole figure plotting by stereographic
projection. Both contour and dot types are available for pole
figures, whereas only dot type is available for inverse pole figure, yet.

It should be machine independent. Make the way treating the data
uniform and intiutive and simple.

--------
examples
--------
import upf
mypf = upf.polefigure(ngrain=8000,
                           # or filename ='texture/08000.cmb'
                           csym = 'cubic', cdim=[1.,1.,1.],
                           cang =[90.,90.,90.])
cnts = mypf.pf(pole=[ [1,0,0],[1,1,0],[1,1,1]],
                    mode='contourf', ifig=2,
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

         m, n = np.shape(data)
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

import upf
mypf = upf.polefigure(epf='???.epf')
..
..
cnts = mypf.pf()
"""

# print __doc__
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
# ----------------------------------------------------------------------c
## import library blocks
import warnings
import logging
# warnings.filterwarnings("ignore")
warnings.filterwarnings("error")

import numpy as np
import matplotlib.pyplot as plt
import matplotlib  # matplotlib as raw
import os, glob, math
from .randomEuler import randomEuler as re
from .euler import euler  # in euler module def euler:
# A-matrix and Euler angles
import time, random
import fortranformat as ff

try:
    import MP
except:
    logging.debug('-' * 50)
    logging.debug('MP was not installed')
    logging.debug('You may clone it and install via:')
    logging.debug('git@github.com:youngung/mpl-lib.git')
    logging.debug('-' * 50)
else:
    from MP import progress_bar

    t2s = progress_bar.convert_sec_to_string
    uet = progress_bar.update_elapsed_time

## Compiling fortran modules has been very very tricky
# If cannot import for_lib correctly, upf will use
# corresponding python modules/functions/methods
try:
    import pf_for_lib

    agr2pol_f = pf_for_lib.agr2pol
    proj_f = pf_for_lib.projection
    euler_f = pf_for_lib.euler
    gr2psa = pf_for_lib.grain2pole_sa
    i_for = True
except:
    i_for = False
    logging.debug('----------------------------------------------------')
    logging.debug('     pf_for_lib is not available in the system')
    logging.debug('   pf_for_lib is a fortran binary that is loaded')
    logging.debug('   in Python using f2py, which may be helpful to ')
    logging.debug('         speed-up calculations in upf.py')
    logging.debug('----------------------------------------------------')

try:
    import joblib
except:
    is_joblib = False
    logging.debug('** joblib was not found - will not be used in TX.upf')
    # print '-'*60
    # print 'One might improve the speed by installing joblib'
    # print 'JOBLIB is to run embarrasingly parallel runs for multiple poles'
    # print "Find about joblib in  https://github.com/joblib/joblib"
    # print '-'*60
else:
    is_joblib = True
    from joblib import Parallel, delayed

pi = math.pi
cos = math.cos;
sin = math.sin

warnings.resetwarnings()
# import pp ## parallel

def cub(filename=None, gr=None, ifig=3, **kwargs):
    """
    arguments
    =========
    filename = None
    gr       = None
    ifig     = 3
    pole     = [[1,0,0],[1,1,0],[1,1,1]],
    cmode    = None):
    """
    if filename is not None:
        mypf = polefigure(filename=filename, csym='cubic')
    elif gr is not None:
        mypf = polefigure(grains=gr, csym='cubic')
    fig = mypf.pf_new(**kwargs)
    return fig


def cubgr(gr=None, ifig=3, poles=[[1, 0, 0], [1, 1, 0], [1, 1, 1]]):
    """
    Arguments
    =========
    gr
    ifig
    poles
    """
    mypf = polefigure(grains=gr, csym='cubic')
    fig = mypf.pf_new(poles=poles, cmap='jet')
    return fig


def pfnorm(data):
    """
    experimental incomplete pole figure preliminary normalization

    data format should be correct
    grid (phi, khi)

    resolution should be 5.0 degress in both phi and khi
    phi range: 0~ 355
    khi range: 0~ 5*(nn-1)

    Argument
    ========
    data
    """
    # All angles are in radian
    if len(data) != 72:
        logging.debug('number of phi grid:  %i' % len(data))
        raise IOError('Unexpected resolution along phi axis')

    dphi = 360. / len(data)
    dphi = dphi * np.pi / 180.  # dphi
    dkhi = 5. * np.pi / 180.  # dkhi
    logging.debug('dkhi, dphi', dphi * 180. / np.pi, dkhi * 180. / np.pi)

    nkhi = len(data[0])
    phi_i = 0.
    phi_f = np.pi * 2.
    khi_i = 0.
    khi_f = dkhi * (nkhi - 1)
    logging.debug('khi range', khi_i, khi_f * 180. / np.pi)
    logging.debug('phi range', phi_i, phi_f * 180. / np.pi)
    # spanned area, i.e., the area of incomplete hemisphere:
    # area = (np.cos(khi_f) - np.cos(khi_i)) * (phi_f - phi_i)

    Nf = 0.

    a = 0.
    b = 0.
    ## below needs vectorization
    for i in range(len(data)):
        for j in range(len(data[0])):
            a = a + np.sin(dkhi * j)

    for i in range(len(data)):
        for j in range(len(data[i])):
            b = b + data[i, j] * np.sin(dkhi * j)

    # for k in xrange(len(data)):
    #     for l in xrange(len(data[k])):
    #         b = b + data[k, l] * np.sin(dkhi*l)

    for i in range(len(data)):  # phi
        for j in range(len(data[i])):  # khi
            data[i, j] = data[i, j] * a / b

    return data


def epfformat(mode=None, filename=None):
    """
    Experimental pole figure format controller
    mode:
      "steglich"
      "bruker"  *.uxd file
      "epf" (2011-Oct-6) epf popLA experimental pole figure format
      "xpc" (2017-Feb) xcp format compliant with MAUD

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

    Arguments
    =========
    mode     = None
    filename = None

    Returns
    =======
    data
    max_khi
    hkl
    """

    ## steglich's format
    nskip_calc = 13
    nskip_raw = 28
    ##

    if mode == 'steglich':
        # Calculated or raw pole figure
        i = 0
        logging.debug('filename=', filename)
        if not os.path.isfile(filename):
            raise IOError('file is not available')

        while True:
            try:
                data = np.loadtxt(
                    filename, skiprows=i)
            except:
                i = i + 1
            else:
                logging.debug('number of skipped rows: %i' % i)
                break
            if i > 1000: raise IOError('something is wrong')

        ## raw pole figure format
        if i == nskip_raw:
            # axes: (phi, khi)
            data = data.T  # (khi, phi)

            # upon each axis
            f = open(filename)
            temp = f.readlines()
            khi = list(map(float, temp[nskip_raw - 1].split()[1:]))
            khi = np.array(khi)
            khi = khi * np.pi / 180.
            phi = data[0]  # first column is for phi
            phi = phi * np.pi / 180.
            data = data[1:]  # remove phi
            data = data.T  # (phi, khi) axis

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
            isincomplete = False
            if np.pi / 2. - khi[-1] > tiny:
                logging.debug('Incomplete pole figure')
                logging.debug('khi range: %3.2f ~%3.2f' % (
                    khi[0] * 180. / np.pi, khi[-1] * 180. / np.pi))
                isincomplete = True

                ## normalization
                if input('y(norm), or n(no)>>>') == 'y':
                    data = pfnorm(data)

                dum = np.zeros((np.pi * 2. / dp, np.pi / 2. / dk + 1))
                for i in range(len(data)):
                    for j in range(len(data[i])):
                        dum[i, j] = data[i, j]
                data = dum.copy()
                del dum

        ## calculated pole figure format
        elif i == nskip_calc:  # He had two formats: raw and calculated.
            logging.debug('Calculated pole figure format')
            # axes: (phi, khi)
            data = data.T  # (khi, phi)
            f = open(filename)
            temp = f.readlines()
            khi = list(map(float, temp[nskip_calc - 1].split()[1:]))
            khi = np.array(khi)
            khi = khi * np.pi / 180.
            phi = data[0]
            phi = phi * np.pi / 180.
            data = data[1:]
            data = data.T  # (phi, khi)

            ## dp and dk
            dp = phi[1] - phi[0]
            dk = khi[1] - khi[0]

            ## shift phi back-modification
            phi, data = shiftphi(data, dp, phi)
            ##

            ## check if there is phi=0, 360 at the same time
            phi, data = duplicate0_360(phi, data)
            ##

    elif mode == 'bruker':
        ## make full use of existing uxd.py script
        ## --> normalization is missing...
        ## This must be completed!
        from . import uxd
        logging.debug('You are now in the bruker mode under epfformat')
        logging.debug('given file name is %s' % filename)
        myuxd = uxd.pf(filename=filename, mode='pf')
        if len(myuxd.polefigures) > 1:
            logging.debug('multiple pole figures are found')
            input()
        for i in range(len(myuxd.polefigures)):
            pf = myuxd.polefigures[i]
            pf = pfnorm(pf)  ##normalize

    elif mode == 'epf':
        """
        ready made popLA epf format parser
        consider the possibility of multiple number of polefigure
        ## phi must be 0~355, khi must be 0~90 with 5 as an ang
        resolution
        """
        logging.debug('You are now reading experimental pole figure(s) :%s' % filename)
        blocks = open(filename, 'rU').read().split('(')[1:]
        if len(blocks) == 0:
            msg1 = 'epf parser in upf assumes that hkl is embraced by paratheses'
            msg2 = ' %s has no paranthesis that embraces hkl was found' % filename
            msg = '%s \n %s' % (msg1, msg2)
            logging.debug('-' * 52)
            logging.debug(msg)
            # raise IOError, msg
            logging.debug('upf will keep proceding with using upf.parse_epf method with n_unit=79')
            logging.debug('-' * 52)
            blocks = parse_epf(filename)

        npf = len(blocks)
        if npf == 0: raise IOError('No pf block found.')

        datasets = [];
        max_khi = []
        if npf > 1: hkls = []  ## multiple number of pole figures in a file
        for i in range(npf):
            # d = __epffiletrimmer__(blocks[i]) #only for popLA epf format
            ## epffiletrimmer should return a block of intensity grid map
            ## as well as maximum khi.
            ## --> modification (2011 NOV 8)
            hkl = blocks[i][0:3]  # hkl
            hkl = list(map(int, [hkl[0], hkl[1], hkl[2]]))

            if npf > 1: hkls.append(hkl)

            if blocks[i][3] != ')':
                logging.debug('Caution: unexpected hkl labeling format')

            d, mxk = __epffiletrimmer__(blocks[i])
            # only for popLA epf format
            datasets.append(d)
            max_khi.append(mxk)

        logging.debug("number of pole figures:", len(datasets))

        ## presumption of using 19x72 should be deprecated...
        data = np.zeros((len(datasets), 19, 72))  # npf, nphi, nkhi
        dum = data.copy()
        for i in range(len(datasets)):
            for j in range(len(datasets[i])):
                for k in range(len(datasets[i][j])):
                    data[i, j, k] = datasets[i][j][k]
        ## Swap the axis
        data = data.swapaxes(1, 2)  # from 72 x 19 to 19 x 72

        ## psuedo-normalization
        for i in range(len(data)):
            data[i] = pfnorm(data[i].copy())
        ## --
        if npf > 1:
            hkl = hkls[::]
        else:
            hkl = [hkl]
        return data, max_khi, hkl

    elif mode == 'xpc':
        """
        Adapted .xpc format parser, from
        https://github.com/usnistgov/texture
        commit 9c0ac85
        based on ready made popLA epf format parser
        """
        import pandas as pd
        logging.debug('You are now reading experimental pole figure(s) :%s' % filename)
        blocks = open(filename, 'rU').read().split('\n\n\n\n')[1:]
        logging.debug('There are %s blocks of data found' % len(blocks))
        if len(blocks) == 0:
            msg1 = 'xpc parser in upf assumes that pole figures are separated by 4 new lines'
            msg2 = ' searching %s finds no set of 4 new lines in ' % filename
            msg = '%s \n %s' % (msg1, msg2)
            raise IOError(msg)
            # blocks = parse_epf(filename)

        npf = len(blocks)
        if npf == 0: raise IOError('No pf block found.')

        datasets = [];
        max_khi = []
        if npf > 1: hkls = ["HKL"]  ## multiple number of pole figures in a file

        for part in blocks:
            line = part.split('\n')
            # print len(line)

            # header lines
            structureline = ff.FortranRecordReader('(6f10.4,1x,i4,1x,i4)')  # the head-line in each data block
            [a, b, c, alpha, beta, gamma, crystalclass, something] \
                = structureline.read(line[1])
            pfDefline = ff.FortranRecordReader('(1x,3i3,6f5.1,2i2)')
            [h, k, l, unknown1, tilt, tiltinc, unknown2, rotation, \
             rotationinc, unknown3, unknown4] = pfDefline.read(line[2])

            # for the rest of the lines, do the following
            dataline = ff.FortranRecordReader('(1x,18i4)')

            # Pretty ugly code, but works...
            grouping = [[3, 4, 5, 6], [7, 8, 9, 10], [11, 12, 13, 14], [15, 16, 17, 18], [19, 20, 21, 22],
                        [23, 24, 25, 26], \
                        [27, 28, 29, 30], [31, 32, 33, 34], [35, 36, 37, 38], [39, 40, 41, 42], [43, 44, 45, 46],
                        [47, 48, 49, 50], \
                        [51, 52, 53, 54], [55, 56, 57, 58], [59, 60, 61, 62], [63, 64, 65, 66], [67, 68, 69, 70],
                        [71, 72, 73, 74], \
                        [75, 76, 77, 78]]

            # dataset=np.array([])
            dataset = []
            for item in grouping:
                # print item[0],item[1],item[2],item[3]
                parsed = dataline.read(line[item[0]])
                parsed.extend(dataline.read(line[item[1]]))
                parsed.extend(dataline.read(line[item[2]]))
                parsed.extend(dataline.read(line[item[3]]))
                dataset.append(parsed)

            # print dataset

            # now saves as a dataframe, and wraps 360 to 0 degrees
            # row and column indexes are by degrees
            df = pd.DataFrame(dataset, index=np.arange(0, 91, 5))
            df.columns = [np.arange(0, 360, 5)]
            df[360] = df.ix[:, 0]
            # print df.ix[:,0]

            # df["Tilt Angle"] =np.arange(0,92,5)
            # df.index.names=[np.transpose(np.arange(0,90,5))]
            # print df

            hkl = [h, k, l]  # hkl
            # hkl = map(int, [hkl[0],hkl[1],hkl[2]])
            # print hkl
            hkls.append(hkl)
            datasets.append(df)
            # max_khi.append(mxk)

        # print hkls
        logging.debug("number of pole figures:", len(datasets))

        ### convert the panda data frame to numpy
        npf = len(datasets)
        datasets_compliant = np.zeros((npf, 72, 19))
        for i in range(npf):
            array = datasets[i].values
            arrayt = array.T
            ## ignore the last phi axis as it is repeated (phi360=phi0)
            datasets_compliant[i, :, :] = arrayt[:72, :]

        return datasets_compliant, hkls

    else:
        raise IOError('Unexpected mode is given')


def __epffiletrimmer__(block):
    """
    EPF (experimental pole figure) trimmer for a single block of data
    EPF: The popLA format

    Arguments
    =========
    block
    """
    pfdata = block.split('\n')
    dkhi = float(pfdata[0][5:9])  # khi incremental angle
    fkhi = float(pfdata[0][9:14])  # khi end angle
    dphi = float(pfdata[0][14:19])  # phi incremental angle

    # written as 360 but it is only upto 355
    fphi = float(pfdata[0][19:24]) - dphi

    ## now it is enfored to have 5 deg resolution, which
    ## results in 19x72
    ## intensities along khi and phi axes coordinate
    data = np.zeros((19, 72))
    ## the rest information is neglected from the header

    pfdata = pfdata[1:]
    lines = []
    for i in range(len(pfdata)):
        if len(pfdata[i]) == 72:
            lines.append(pfdata[i][:])
        elif len(pfdata[i]) == 73:
            lines.append(pfdata[i][1:])

    pfdata = lines  # (maxkhi/dkhi + 1) * 4
    if len(pfdata) != 76:  # 76 =  ((90 / 5) + 1) * 4
        # (4lines for one khi level)
        logging.debug('len(pfdata) =', len(pfdata))
        logging.debug(pfdata)
        raise IOError('Unexpected pfdata format or type')

    if True:
        for j in range(19):  # 90 / 5 + 1 #number of khi threads
            kh = pfdata[j * 4: (j + 1) * 4]  # along the same kh level
            khline = ''
            for k in range(4):
                khline = khline + kh[k]
            # khline  # all data in a thread of
            # string every 4digits corresponds to a datum
            kp = 0  # each point along phi
            for k in range(72):  # 72 = 288 / 4 (4digits) : 0~355 5.0
                datum = khline[k * 4:(k + 1) * 4]
                data[j, k] = int(datum)

    # block of mesh grid pole figure map and khi final angle.
    return data, fkhi


def shiftphi(data, dp, phi):
    """
    shifted phi modification

    Arguments
    =========
    data
    dp
    phi
    """
    tiny = 0.0000001
    if abs(phi[0] - 0.) > tiny > abs((dp - phi[0]) - dp / 2.):
        dum = data.copy()
        for i in range(len(data)):  # phi axis
            for j in range(len(data[i])):  # khi axis
                dum[i, j] = (data[i - 1, j] + data[i, j]) / 2.
        data = dum.copy()
    return phi, data


def duplicate0_360(phi, data):
    """
    If both phi=0 and phi=360 exist,
    delete phi=360

    Arguments
    =========
    phi
    data
    """
    tiny = 0.00000000001
    isallsame = True
    for i in range(len(data[0])):
        logging.debug(data[0, i])
        logging.debug(data[-1, i])
    if any(abs(data[0, i] - data[-1, i]) > tiny for i in range(len(data[0]))):
        logging.debug(data[0, i])
        logging.debug(data[-1, i])
        isallsame = False

    isphioverlapped = False
    logging.debug(abs(phi[0] - phi[-1] - 2. * np.pi))
    if abs(phi[0] - phi[-1] + 2. * np.pi) < tiny:
        isphioverlapped = True

    if isphioverlapped and isallsame:
        newdata = np.zeros((data.shape[0] - 1, data.shape[1]))
        for i in range(len(newdata)):
            for j in range(len(newdata[i])):
                newdata[i, j] = data[i, j]
        return phi[:-1], newdata
    elif isphioverlapped and isallsame == False:
        logging.debug("conflict results!")
        logging.debug("phi=0, and 360 turned out to coexist but")
        logging.debug("Values along these two axes are not equivalent")
        raise IOError
    elif not isphioverlapped and isallsame:
        logging.debug("Conflict results!")
        logging.debug("Phi is overlapped but phi[0] and phi[-1] is same")
        raise IOError
    else:
        logging.debug("No duplicated axis is found")
        return phi, data


def pole2f(poles_sa, poles_wgt, dth, dph, f):
    tiny = 1e-9
    for i in range(len(poles_sa)):
        theta, phi = cart2sph(poles_sa[i])
        ix = int((theta * 180. / np.pi + 180) / dth - tiny)
        iy = int((phi * 180. / np.pi) / dph - tiny)
        f[ix, iy] = f[ix, iy] + poles_wgt[i]
    return f


def cart2polar(x, y):
    """
    cartesian to polar coordinate

    Arguments
    =========
    x
    y
    """
    r = np.sqrt(x ** 2 + y ** 2)
    theta = np.arctan2(y, x)
    return r, theta


def circle(center=[0, 0], r=1.):
    """
    Draw a circle around the given center point with
    a radius of r(as given)
    The default settings are as below.

    Arugments
    =========
    center = [0,0]
    r = 1.
    """
    ang = np.linspace(0, 2 * np.pi, 1000)
    # unit circle * radius
    x = np.cos(ang) * r
    y = np.sin(ang) * r
    # circle transloation
    x = x + center[0]
    y = y + center[0]
    return x, y


def basic_triangle():
    """
    provide the boundary of basic triangle
    """
    v100_ca = np.array([1, 0, 0])
    v110_ca = np.array([1, 1, 0]) / np.sqrt(2.)
    v111_ca = np.array([1, 1, 1]) / np.sqrt(3.)

    delt_thet_100_110 = np.pi / 2.0
    delt_thet_111_110 = \
        np.arccos(np.dot(v110_ca, v111_ca))
    delt_thet_100_111 = \
        np.arccos(np.dot(v100_ca, v111_ca))

    # line1  # 100_110
    # line2  # 110_111
    # line3  # 111_100


def trace(thet_f, n, v1, v2):
    """
    Trance of vectors that runs on the plane made
    by the two given vectors

    ---------
    Arguments
    ---------
    thet_f (in Radian)
    n
    v1
    v2
    """
    v1 = v1 / np.sqrt(v1[0] ** 2 + v1[1] ** 2 + v1[2] ** 2)
    v2 = v2 / np.sqrt(v2[0] ** 2 + v2[1] ** 2 + v2[2] ** 2)
    v3 = np.cross(v1, v2)
    v3 = v3 / np.sqrt(v3[0] ** 2 + v3[1] ** 2 + v3[2] ** 2)

    dthet = np.linspace(0., thet_f, n)
    # around v3 rotate v1 by incremental angle upto thet_f
    # rotate v1 about v3 axis.

    trace = []
    for iang in range(n):
        th = dthet[iang]
        r = vector_ang(v3, th)
        vnew = np.dot(r, v1)
        trace.append(vnew)
    return trace


def xytrace(thet_f, n, v1, v2):
    """
    Returned projected trace
    """
    poles = trace(thet_f, n, v1, v2)
    X, Y = [], []
    for i in range(len(poles)):
        x, y = projection(poles[i], agrain=[0, 0, 0, 1])
        X.append(x)
        Y.append(y)
    return X, Y


def vector_ang(u, th):
    """
    transformation matrix that rotates a vector about an axis u
    by the angle, th.

    ---------
    Arguments
    ---------
    u
    th (in Radian)
    """
    pi = np.pi
    idx = np.zeros((3, 3))
    r = np.zeros((3, 3))
    for i in range(3):
        idx[i, i] = 1.

    ct = np.cos(th)
    st = np.sin(th)
    cm = crossop(u)

    for i in range(3):
        for j in range(3):
            r[i, j] = idx[i, j] * ct + st * cm[i, j] + \
                      (1 - ct) * u[i] * u[j]

    return r


def crossop(u):
    m = np.zeros((3, 3))
    m[0, 1] = -u[2]
    m[0, 2] = u[1]
    m[1, 0] = u[2]
    m[1, 2] = -u[0]
    m[2, 0] = -u[1]
    m[2, 1] = u[0]
    return m


def __isunique__(a, b):
    """
    Is a(3) in b(m, 3)

    Arguments
    =========
    a
    b
    """
    for i in range(len(b)):
        ## is a//b[i]? either + or -
        diff0 = abs(b[i][0] - a[0]) + \
                abs(b[i][1] - a[1]) + abs(b[i][2] - a[2])
        diff1 = abs(b[i][0] + a[0]) + abs(b[i][1] + a[1]) \
                + abs(b[i][2] + a[2])
        if diff0 < 0.1 ** 4 or diff1 < 0.1 ** 4:
            return False
    return True


def __circle__(center=[0, 0], r=1.):
    """
    Draw a circle around the given center point with a
    radius of r(as given)
    The default settings are as below.

    Arugments
    =========
    center = [0,0]
    r = 1.
    """
    ang = np.linspace(0, 2 * np.pi, 1000)
    # unit circle * radius
    x = np.cos(ang) * r
    y = np.sin(ang) * r
    # circle transloation
    x = x + center[0]
    y = y + center[0]
    return x, y


def deco_pf(ax, cnt, miller=[0, 0, 0],
            iopt=0, iskip_last=False,
            ix='1', iy='2', mode='line'):
    """
    Decorate matplotlib.pyplot.axes used for plotting pole figures

    iopt==1: skip level lines

    Arguments
    ---------
    ax     (matplotlib axes)
    cnt    (contour object generated by ax.contour)
    miller (miller indices)
    iopt
    iskip_last (whether or not to skip the maximum
                iso-contour line or level)
    ix     (xlabel, i.e., horizontal axis label)
    iy     (ylabel, i.e., vertial axis label)
    mode   pole figure plotting mode (line or fill)
    """
    ## fontsize of appended text to pole figures will be 4*fact
    fact = 2.5
    # --

    clev = cnt._levels
    tcolors = cnt.tcolors
    if iskip_last:
        nlev = len(tcolors) - 1
    else:
        nlev = len(tcolors)

    # print 'len(tcolors):',len(tcolors)
    # print 'nlev in deco_pf:',nlev

    if iopt == 1:
        pass
    elif iopt == 0:
        for i in range(nlev):
            cc = tcolors[i][0][0:3]

            x = [1.30, 1.37]
            y = [1. - i * 0.2, 1. - i * 0.2]

            if not iskip_last and i == nlev - 1 and mode == 'line':
                ax.plot((x[0] + x[1]) / 2., (y[0] + y[1]) / 2.,
                        '+', mew=2., color=cc)
            else:
                ax.plot(x, y, color=cc)

            ## level text
            if clev[i] < 10:
                s = '  %4.2f' % clev[i]
            else:
                s = '%5.2f' % clev[i]

            ax.text(x=1.44, y=1. - i * 0.2 - 0.05,
                    s=s, fontsize=4. * fact)

    ## axis label/    ## Ticks
    ax.text(1.14, 0., ix, va='center', ha='center')
    ax.text(0., 1.14, iy, va='center', ha='center')
    ax.plot([0.0, 0.0], [0.97, 1.00], 'k-')
    ax.plot([0.97, 1.00], [0.0, 0.0], 'k-')

    ## pole
    s = '('
    for i in range(len(miller)):
        if miller[i] < 0:
            h = r'\bar{%s}' % str(-1 * miller[i])
        else:
            h = '%s' % str(miller[i])
        s = '%s%s' % (s, h)
    s = '%s)' % s
    s = r'$\mathbf{%s}$' % s
    ax.text(0.6, -0.95, s, fontsize=12)

    ## circle
    x, y = __circle__()
    ax.plot(x, y, 'k-')
    ax.set_axis_off()
    ax.set_xlim(-1.1, 1.4)
    ax.set_ylim(-1.1, 1.4)
    ax.set_aspect('equal')


def projection(pole=None, agrain=None):
    """
    Projects a pole (vector) to projection plane.
    (default is stereography projection)

    pole = [1,1,1] or [1,1,0] something like this.

    Arguments
    ---------
    pole = None
    agrain = [ph1, phi, phi2, vf]
    """
    # normalization of the miller indices
    pole = pole / np.sqrt(pole[0] ** 2 + pole[1] ** 2 + pole[2] ** 2)
    a, b, c = pole[0:3]
    ###  mid-plane projection (z=0)
    if c == 1:
        # pole[0] = 0; pole[1]=0; pole[2] = 1
        X = 0;
        Y = 0
    else:
        X = a / (c - 1)
        Y = b / (c - 1)
    return X, Y


def invproj(x=None, y=None):
    """
    Converts projected point to the pole

    Arguments
    =========
    x = None
    y = None
    """
    X = 2 * x / (1 + x ** 2 + y ** 2)
    Y = 2 * y / (1 + x ** 2 + y ** 2)
    Z = (-1 + x ** 2 + y ** 2) / (1 + x ** 2 + y ** 2)
    return np.array([X, Y, Z])


def vect2sphe(pole):
    """
    cartensian vector to spherical cooridnate system

    Argument
    ========
    pole
    """
    seca = math.sqrt(pole[0] ** 2 + pole[1] ** 2)
    if seca < 0.1 ** 6:
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

    Argument
    ========
    pole

    Returns
    -------
    [phi,theta]
    """
    pole = np.array(pole)
    r = np.sqrt((pole ** 2).sum())
    theta = math.acos(pole[2] / r)
    phi = math.atan2(pole[1], pole[0])
    return np.array([phi, theta])


def agr2pol(agrain=None, miller=None, proj=None):
    """
    -- for pole figure projection (proj='pf')
    For the given grain, crystallographic direction
    is mapped onto the sample axes by doing tensor
    product of the pole vector with rotation matrix.

    -- for inverse pole figure projection (proj='ipf')
    For the given grain, a vector in sa transformed to
    one referred in the crystal axes.

    Arguments
    =========
    agrain = None
    miller = None
    proj   = None
    """
    if proj is None:
        logging.debug("argument proj should be given"); raise IOError
    elif proj != 'pf' and proj != 'ipf':
        logging.debug(" proj should be either 'pf' or 'ipf'")
        raise IOError
    if type(miller).__name__ == 'list': miller = np.array(miller)

    # a-matrix between the grain's coordinate and sampl coordinate
    phi1 = agrain[0];
    phi = agrain[1];
    phi2 = agrain[2]  # VPSC convention
    # rotation matrix (sa->ca)
    amat = euler(ph=phi1, th=phi, tm=phi2, echo=None)

    # normalize miller
    norm = math.sqrt(miller[0] ** 2 + miller[1] ** 2 + miller[2] ** 2)
    miller = miller / norm

    if proj == 'pf':
        # Returns the dot product of the
        # transposed A-matrix and miller vector
        "A^{T}_{ij} * V_{j}"
        return np.dot(amat.transpose(), miller)  # ca to sa

    elif proj == 'ipf':
        # #for inverse projection,
        # #the sample axis should be one of principal axes.(x or y or z)
        # if   miller[0]==1: p_ca=0
        # elif miller[1]==1: p_ca=1
        # elif miller[2]==1: p_ca=2
        # # else: print"something is wrong"; raise IOError
        # # map the pole in the sample axes into the crystal axes
        "A_{ij} * V_{j}"
        return np.dot(amat, miller)  # sa to ca returns the pole in ca
    else:
        logging.debug("projection should be pf or ipf")
        raise IOError


def ipfline(center=[0, 0], csym='cubic'):
    """
    The boundary line for inverse pole figure
    ---------
    Arguments
    ---------
    center = [0,0]
    csym   = 'cubic'
    """
    xc = [];
    yc = []
    if csym != 'cubic': logging.debug("Only Cubic!"); raise IOError
    xc.append(center[0])
    yc.append(center[1])

    for i in np.linspace(0., 1 / math.sqrt(3.)):
        yaux = i
        xaux = math.sqrt((1. - yaux ** 2) / 2.)
        zaux = xaux
        t1 = math.sqrt(1. - zaux) / math.sqrt(xaux ** 2 + yaux ** 2)
        t2 = t1 / math.sqrt(1. + zaux)
        ## equal area
        # xc.append(xaux*t1)
        # yc.append(yaux*t1)
        ## stereo
        xc.append(xaux * t2)
        yc.append(yaux * t2)

    xc.append(center[0])
    yc.append(center[1])
    return np.array([xc, yc])


"""
Sample symmetry application is performed over RVE calculation
refer to the RVE class in cmb.py module.
"""


class polefigure:
    # decides if the given set is in the texture file form or array
    def __init__(self, grains=None, filename=None, csym=None,
                 ngrain=100, cdim=[1., 1., 1.], cang=[90., 90., 90.],
                 ssym=False, epf=None, epf_mode=None):
        """
        cdim=[1.,1.,1.6235] ## AZ31

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
        ---------
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

        # Turn interactive plotting off
        plt.ioff()

        if type(grains) == type(None) and type(filename) == type(None) \
                and type(epf) == type(None):
            logging.debug(" ****************************** ")
            logging.debug(" Since no argument is passed,   ")
            logging.debug(" 1000 random grains are created ")
            logging.debug(" ****************************** \n")
            from .cmb import random
            self.gr = random(phi1=360, phi2=360, phi=90, ngrain=1000, iplot=False)

        self.epf = epf  # global

        if type(grains) != type(None):
            self.gr = np.array(grains)
        elif type(filename) != type(None):
            self.gr = np.loadtxt(fname=filename, skiprows=4)
        elif type(epf) != type(None):  # None is the default for epf
            """
            experimental pole figures..
             # available formats:
                 - UXD
                 - steglich
                 - bruker
                 - epf*
                 - xpc
            """
            if type(epf).__name__ == 'list':
                self.epf_fn = epf
            elif type(epf).__name__ == 'str':
                self.epf_fn = [epf]
            elif type(epf) == type(True):
                fn = []  # list of file names
                logging.debug('type the experimental pole figure file names')
                logging.debug("To finish input, press enter")
                while True:
                    dum = input(">>> ")
                    if len(dum) == 0: break
                    fn.append(dum)
                self.epf_fn = fn
            else:
                raise IOError('Unexpected epf type found')

            ## check if the file name is correct ##
            for i in range(len(self.epf_fn)):
                if not (os.path.isfile(self.epf_fn[i])):
                    raise IOError("Could not find %s" % self.epf_fn[i])
            ## --------------------------------- ##

            ## POLE FIGURE MODE --------------------------------------
            if type(epf_mode) == type(None):
                logging.debug("Type the experimental polfe figure mode")
                logging.debug("Available options:", end=' ')  # continuation
                logging.debug("bruker, steglich, epf, xpc (default: %s)" % 'epf')
                epf_mode = input(" >>>")
                if len(epf_mode) == 0: epf_mode = 'epf'  # default
            ##---------------------------------------------------------

            self.grid = [];
            self.hkl = []
            ## more than one pf can be included.
            npole_per_file = []
            if epf_mode == 'epf': self.max_khi = []  # Available only for epf_mode yet.

            for i in range(len(self.epf_fn)):
                if epf_mode == 'epf':
                    data, maxk, hkl = epfformat(
                        mode=epf_mode,
                        filename=self.epf_fn[i])
                    # one file may include multiple poles
                    for j in range(len(data)):
                        self.grid.append(data[j])
                        self.max_khi.append(maxk[j])
                        self.hkl.append(hkl[j])
                    npole_per_file.append(len(data))  # of pole per a file
                elif epf_mode == 'xpc':
                    data, hkl = epfformat(
                        mode=epf_mode,
                        filename=self.epf_fn[i])
                    for j in range(len(data)):
                        self.grid.append(data[j])
                        # self.max_khi.append(90.)
                        self.hkl.append(hkl[j])
                else:
                    data = epfformat(
                        mode=epf_mode,
                        filename=self.epf_fn[i])
                    self.grid.append(data)

                    self.hkl.append(None)
            self.grid = np.array(self.grid)
            self.epf_mode = epf_mode

        ## EXPERIMENTAL POLE FIGURE
        ## ------------------------------------------------------- ##
        ## POLE FIGURES BINNED FROM THE POLYCRYSTALLINE AGGREGATES ##

        if epf is None:
            dat = self.gr.transpose()
            phi1 = dat[0];
            phi = dat[1];
            phi2 = dat[2]
            ph1min, ph1max = int(
                round(min(dat[0] / 90.))) * 90, int(
                round(max(dat[0] / 90.))) * 90
            phmin, phmax = int(
                round(min(dat[1] / 90.))) * 90, int(
                round(max(dat[1] / 90.))) * 90
            ph2min, ph2max = int(
                round(min(dat[2] / 90.))) * 90, int(
                round(max(dat[2] / 90.))) * 90

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
            # 1 symmetry
            self.csym = csym
            self.ngr = len(self.gr)
            self.cdim = cdim
            self.cang = cang

    def epfplot(self, ifig, cmap, nlev, mn, mx, ix, iy, rot, iline_khi80):
        """
        This function is expected to be called
        within self.pf or self.pf_new

        if type(self.epf).__name__!='NoneType':
            ## This function is called within self.pf or self.pf_new

        * pole figure master data: self.grid
          [i]: i-th pole
          [i][j]: i-th pole's j-th phi segment's khi
          [i][j][k]: An intensity (count) at i-th pole's j-th phi sement at the k-th khi

        Arguments
        =========
        ifig
        cmap    ## color map, e.g., 'jet' or 'gray_r'
        nlev
        mn
        mx
        ix
        iy
        rot
        iline_khi80
        """
        from matplotlib.colors import LogNorm
        logging.debug('List of files:')
        for f in self.epf_fn: logging.debug('%s ' % f)
        logging.debug('dimension of self.grid:', self.grid.shape)
        logging.debug('ifig=', ifig)

        fact = 2.
        nrow = len(self.grid)
        fig = plt.figure(figsize=(3.3 * nrow, 3.0))

        logging.debug('self.grid.shape:', self.grid.shape)
        mns, mxs, indices_mx = self.calcMXN(self.grid, mx, mn, 'line', 1)

        # print 'mns:',mns
        # print 'mxs:',mxs

        ## loop over each of pole figures
        for ip in range(len(self.grid)):  # upon each of self.eps_fn
            ax = fig.add_subplot(1, nrow, ip + 1)
            pf = np.zeros((self.grid[ip].shape[0] + 1,
                           self.grid[ip].shape[1]))
            for i in range(len(self.grid[ip])):
                for j in range(len(self.grid[ip][i])):
                    pf[i, j] = self.grid[ip][i][j]
                    pf[-1, j] = self.grid[ip][0][j]

            # if type(mx).__name__=='NoneType':
            #     mx = np.array(pf).flatten().max()
            # if type(mn).__name__=='NoneType':
            #     mn = np.array(pf).flatten().min()

            # if mn==0: mn=0.5
            # if mx>100: mx=99.
            # levels = np.logspace(
            #     np.log10(mn),np.log10(mx),nlev)

            if mns[ip] == 0: mns[ip] = 0.5
            levels = np.logspace(
                np.log10(mns[ip]), np.log10(mxs[ip]), nlev)

            norm = LogNorm()

            nm = len(pf);
            nn = len(pf[0])  # phi, khi
            dp = 360. / nm;
            dk = 90. / nn
            khi = np.linspace(np.pi, np.pi / 2., nn)
            phi = np.linspace(0., 2. * np.pi, nm)
            r = np.sin(khi) / (1 - np.cos(khi))
            R, PHI = np.meshgrid(r, phi)
            PHI = PHI + rot * np.pi / 180.  # default=0.

            x = R * np.cos(PHI);
            y = R * np.sin(PHI)

            logging.debug('levels in epfplot:')
            logging.debug(levels)

            # cnt=ax.contour(
            cnt = ax.contourf(
                x, y, pf, levels=levels,
                cmap=cmap, norm=norm)
            deco_pf(ax=ax, cnt=cnt, miller=self.hkl[ip], ix=ix, iy=iy)

            ## dots like in pf_new.
            xs = [];
            ys = []
            for j in range(len(x) - 1):
                for k in range(len(x[j])):
                    if pf[j, k] < levels[0]:
                        if k == 0 and j > 1:
                            pass
                        else:
                            xs.append(x[j][k])
                            ys.append(y[j][k])
            if len(xs) > 0:
                ax.plot(xs, ys, 'k.',
                        alpha=0.17, markersize=2.0)

            if iline_khi80:
                max_khi = 80.
                r_khi = 2 - np.sin(max_khi * np.pi / 180.) / (1 - np.cos(max_khi * np.pi / 180.))
                rx, ry = __circle__(center=[0, 0], r=r_khi)
                ax.plot(rx, ry, '--', color='gray')

        return fig

    def pf_axis(self, pole=[[1, 0, 0]], ifig=1):
        """
        Plot each pole without crystal symmetry
        """
        color = ['r', 'b', 'g', 'k', 'gray']
        # marker =['o','x','+','d','.']
        for ip in range(len(pole)):
            cl = color[ip]
            # mk = marker[i]
            for i in range(len(self.gr)):
                tm = self.dotplot(proj='pf', agrain=self.gr[i],
                                  npole=len(pole), ipole=ip + 1,
                                  pole=pole[ip], ifig=ifig,
                                  cdim='None', cang=self.cang,
                                  csym=self.csym, mode=None,
                                  color=cl)

    def pf2xyw(self, pole=[1, 0, 0], csym='cubic', cdim=[1., 1., 1.],
               cang=[90., 90., 90.], fn='dat.xyz'):
        """
        Read pole and write xyw to a file

        Arguments
        =========
        pole
        """
        f = open(fn, 'w')
        xyzw = []
        for i in range(len(self.gr)):
            gr = self.gr[i][::]
            phi1, phi, phi2 = gr[:3:]
            phi1 = phi1 - 90.

            npeq = __equiv__(
                miller=pole, csym=csym, cdim=cdim, cang=cang)

            xy, POLE = self.core(
                pole=pole, proj='pf', csym=csym,
                agrain=gr, isym=True,
                cdim=cdim, cang=cang, equivp=npeq)

            w = gr[-1]
            # for j in xrange(len(xy)):
            #     x,y = xy[j]
            #     z = 0
            #     f.write('%4.2f %4.2f %4.2f %11.4e\n'%(x,y,z,w))
            #     xyzw.append([x,y,z,w])


            for j in range(len(POLE)):
                xyz = POLE[j]
                x, y, z = xyz
                f.write('%4.2f %4.2f %4.2f %11.4e\n' % (x, y, z, w))
                xyzw.append([x, y, z, w])

        f.close()
        return np.array(xyzw).T

    def pf(self, pole=[[1, 0, 0], [1, 1, 0], [1, 1, 1]], mode='contour',
           ifig=1, dm=7.5, dn=7.5, ssym=None, levels=None,
           axes=None, cmode='gray_r', rot=0., proj='pf', pole_mode='sys',
           poles_gr=None, ix='2', iy='1'):
        """
        Plots pole figures, either experimental pole figure or
        polycrystalline aggregate.

        ** Exemplary cmodes:
        'gray_r, gray, pink_r, pink, summer_r, summer, winter_r, winter,
        Blues, Blues_r, Set1_r, Set1 .... '

        ---------
        Arguments
        ---------
          pole = [[1,0,0,], [1,1,1], [1,1,0]]
          mode ='dot', 'contour', 'contourf', 'im'
          ifig = 1
          dm = 7.5
          dn = 7.5
          ssym = None ## --> dummy yet (2011-Sept-13)
          levels = None
          axes = None
          cmode = 'gray_r'  #the color modes under matplotlib.cm
          rot = 0.
          proj = 'pf'
          pole_mode='sys'
          ix = '2'  (label for horinzontal axis)
          iy = '1'  (label for vertical axis)

        ** Exemplary cmodes:
        gray_r, gray, pink_r, pink, summer_r, summer, winter_r, winter,
        Blues, Blues_r, Set1_r, Set1 .... '
        """

        logging.debug('pf is deprecated. Considering using pf_new')

        ## PF ploting directly from experimental pole figure
        if type(self.epf) != type(None):
            # mode is fixed to be 'im'
            self.epfplot(
                ifig=ifig,
                cmap=cmode,
                rot=rot, nlev=7, mn=None, mx=None,
                iline_khi80=False)
            return
        ## end of experimental polefigure plotting

        ## if it is not an experimental polefigure,
        ## binning of COD on the spherical surface is carried out.
        if self.csym == 'hexag' or self.csym == 'trigo':
            p__ = pole  ## for later use on pole index indication
            pole_tmp = []
            for ip in range(len(pole)):
                p_ = pole[ip]
                p = [0, 0, 0]
                if len(pole[ip]) != 4:
                    input('pole must be four digit for hexag')
                    raise IOError
                p[2] = p_[3]
                p[0] = p_[0] - p_[2]
                p[1] = p_[1] - p_[2]
                p[2] = p_[3]
                # print 'pole as below'
                # print p[0:3]
                pole_tmp.append(p[0:3])
            pole = pole_tmp
        pi = math.pi
        ## pole figure
        temp = []
        if mode == 'dot':
            for ip in range(len(pole)):
                color = ['r', 'b', 'g', 'k', 'gray']
                marker = ['o', 'x', '+', 'd', '.']
                for i in range(len(self.gr)):
                    if len(self.gr) < 5:
                        cl = color[i]
                        mk = marker[i]
                    else:
                        cl = 'k'
                        mk = '.'
                    tm = self.dotplot(
                        proj='pf', agrain=self.gr[i],
                        npole=len(pole), ipole=ip + 1,
                        pole=pole[ip], ifig=ifig,
                        cdim=self.cdim, cang=self.cang,
                        csym=self.csym, mode=None,
                        color=cl)

        elif mode == 'trace':
            for ip in range(len(pole)):
                for i in range(len(self.gr)):
                    if i == 0:
                        color = 'red'
                        alpha = 1.0
                        marker = 'o'
                    elif i == -1:
                        color = 'blue'
                        alpha = 1.0
                        marker = 'o'
                    tm = self.dotplot(
                        proj='pf', agrain=self.gr[i],
                        alpha=alpha, color=color,
                        marker=marker,
                        pole=pole[ip], ifig=ifig,
                        cdim=self.cdim, cang=self.cang,
                        csym=self.csym, mode=None)
                    color = 'black'
                    marker = '.'
                    alpha = 0.5

        elif mode in ['contour', 'contourf']:
            cnts = []  # cnt container
            ## figure setting
            ## figure size
            fact = 2.  # size factor
            figsize = (len(pole) * 2. * fact, 1. * 2. * fact)

            ## If axis is passed to the module, plotting is performed on them.
            if axes is not None:
                if len(axes) != len(pole): raise IOError
                fig = plt.figure(ifig)
            elif axes is None:
                fig = plt.figure(ifig, figsize=figsize)
                nrow = len(pole)
            levs = []

            #### on each of the pole --------------------------------- #
            N = []
            start = time.time()
            self.pfnodes = []
            for ip in range(len(pole)):
                # Polar cell and nodes generation.
                f, nodes = self.cells(
                    pole=pole[ip], ifig=None, dm=dm, dn=dn,
                    cdim=self.cdim, cang=self.cang,
                    csym=self.csym, proj=proj, pole_mode=pole_mode,
                    poles_gr=poles_gr)

                N.append(nodes)
                logging.debug('node shape:', nodes.shape)
                self.pfnodes.append(nodes)

            logging.debug("%5.2f seconds elapsed during calling" \
                  " self.cells\n" % (time.time() - start))
            del nodes

            ## resolution and pole figure plotting preferences
            nm = (360.0 - 0.) / dm;
            nn = (180. - 90.) / dn
            # theta and phi and stereographic projection of them.
            theta = np.linspace(pi, pi / 2., nn + 1)  # tilting angle
            phi = np.linspace(0., 2. * pi, nm + 1)  # rotation angle
            r = np.sin(theta) / (1 - np.cos(theta))  # tilting angle to radius
            R, PHI = np.meshgrid(r, phi)  # meshing radius and rotation angle
            PHI = PHI + pi / 2.  # rotation the pole figure up.
            x = R * np.cos(PHI);
            y = R * np.sin(PHI)  # convert the polar coord-> cartensian

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
                nodes = N[ip]  ## intensity at the current pole (ip)
                if axes is not None:
                    ax = axes[ip]
                elif axes is None:
                    ax = fig.add_subplot(1, nrow, ip + 1)
                ax.set_frame_on(False)
                # contour plotting
                if mode == 'contour':
                    if type(levels) == type(None):
                        if type(cmode) != type(None):
                            cnt = ax.contour(
                                x, y, nodes, cmap=plt.cm.cmap_d[cmode])
                        else:
                            cnt = ax.contour(x, y, nodes)
                    elif type(levels) != type(None):
                        if type(cmode) != type(None):
                            cnt = ax.contour(
                                x, y, nodes, levels=levels, cmap=plt.cm.cmap_d[cmode])
                        else:
                            cnt = ax.contour(x, y, nodes, levels)

                    if proj == 'ipf':
                        v100_ca = np.array([0, 0, -1])
                        v110_ca = np.array([-1, 0, -1]) / np.sqrt(2.)
                        v111_ca = np.array([-1, -1, -1]) / np.sqrt(3.)

                        delt_thet_100_110 = np.pi / 4
                        delt_thet_111_110 = \
                            np.arccos(np.dot(v110_ca, v111_ca))
                        delt_thet_100_111 = \
                            np.arccos(np.dot(v100_ca, v111_ca))

                        X1, Y1 = xytrace(delt_thet_100_110, 20,
                                         v100_ca, v110_ca)  # 100_110
                        X2, Y2 = xytrace(delt_thet_111_110, 20,
                                         v111_ca, v110_ca)  # 110_111
                        X3, Y3 = xytrace(delt_thet_100_111, 20,
                                         v100_ca, v111_ca)  # 111_100

                        ax.plot(X1, Y1, 'k-')  # , alpha=0.5)
                        ax.plot(X2, Y2, 'k-')  # , alpha=0.5)
                        ax.plot(X3, Y3, 'k-')  # , alpha=0.5)

                elif mode == 'contourf':
                    if type(levels) == type(None):
                        if type(cmode) != type(None):
                            cnt = ax.contourf(
                                x, y, nodes, cmap=plt.cm.cmap_d[cmode])
                        else:
                            cnt = ax.contourf(x, y, nodes);pass
                    elif type(levels) != type(None):
                        if type(cmode) != type(None):
                            cnt = ax.contourf(
                                x, y, nodes, levels, cmap=plt.cm.cmap_d[cmode])
                        else:
                            cnt = ax.contourf(
                                x, y, nodes, levels)  # , cmap=plt.cm.bone)

                cnts.append(cnt)
                clev = cnt._levels
                levs.append(clev)

                # Contour's details.
                ax.set_axis_off()
                ax.set_aspect('equal')
                rx, ry = circle()
                ax.plot(rx, ry, 'k')

                # misc (decoration of the axes, with information)
                tcolors = cnt.tcolors
                for i in range(len(tcolors)):
                    cc = tcolors[i][0][0:3]
                    # if levels==None:
                    if type(levels) == type(None) or ip == len(pole) - 1:
                        ## level line
                        ax.plot([1.28, 1.35],
                                [1. - i * 0.2, 1. - i * 0.2],
                                color=cc)
                        ## level text
                        ax.text(x=1.47, y=1. - i * 0.2 - 0.05,
                                s='%3.2f' % (clev[i]),
                                fontsize=4. * fact)
                    ## pole plane
                    if self.csym == 'hexag' or self.csym == 'trigo':
                        ax.text(x=0.4, y=-1.18, s='(%1i%1i%1i%1i)' %
                                                  (p__[ip][0], p__[ip][1],
                                                   p__[ip][2], p__[ip][3]),
                                fontsize=6. * fact)
                    else:
                        if proj == 'pf':
                            s = '(%1i%1i%1i)' % (
                                pole[ip][0], pole[ip][1], pole[ip][2])
                        elif proj == 'ipf':
                            ifac = 1  # 0
                            s = '~(%1i%1i%1i)' % (
                                pole[ip][0] * ifac, pole[ip][1] * ifac,
                                pole[ip][2] * ifac)

                        ax.text(x=0.4, y=-1.18, s=s,
                                fontsize=6. * fact)
                    if proj == 'pf':
                        ## RD and TD indication
                        ax.text(x=-0.05, y=1.05,
                                s=iy, fontsize=4. * fact)
                        ax.text(x=1.05, y=0.,
                                s=ix, fontsize=4. * fact)
                # Fixes the frame
                # if proj!='ipf':
                ax.set_xlim(-1.2, 1.5);
                ax.set_ylim(-1.2, 1.5)
            #### on each of the pole
            ####---------------------------
            ####---------------------------
            # for ip in xrange(len(pole)):
            return cnts  # for mode in ['contour', 'contourf']

        elif mode == 'im':
            fact = 2.  # plot size factor
            figsize = (len(pole) * 2. * fact, 1. * 2. * fact)

            if axes is not None: raise IOError("'im' mode does not support imposition of axes")
            fig = plt.figure(ifig, figsize=figsize)
            fig.clf()  # clear the figure
            nrow = len(pole)
            Zs = []
            start = time.time()
            logging.debug('dm, dn:', dm, dn)
            for ip in range(len(pole)):
                f, Z = self.cells(
                    pole=pole[ip], ifig=None,
                    dm=dm, dn=dn,
                    cdim=self.cdim, cang=self.cang,
                    csym=self.csym)
                Zs.append(Z)
            logging.debug("%5.2f seconds elapsed during calling self.cells\n" % (
                time.time() - start))
            del Z
            ## resolution and pole figure plotting preferences
            nm = (360.0 - 0.) / dm;
            nn = (180. - 90.) / dn
            # theta and phi and stereographic projection of them.
            theta = np.linspace(pi, pi / 2., nn + 1)  # tilting
            phi = np.linspace(0., 2. * pi, nm + 1)  # rotation
            r = np.sin(theta) / (1 - np.cos(theta))  # radius
            phi = phi + pi / 2.  # rotate the RD to the north pole up (for 'im')
            phi = phi + rot  # arbitrary rotation (2011 OCT)

            # phi, r = np.meshgrid(phi,r)
            zmax = 0.
            for ip in range(len(pole)):
                if np.max(Zs[ip]) > zmax: zmax = np.max(Zs[ip].copy())

            for ip in range(len(pole)):
                Z = Zs[ip]  ## intensity of the current pole
                nodes = Z.copy()
                Z = Z.transpose()
                ## polar coordinate system
                axp = fig.add_subplot(
                    1, nrow, ip + 1, polar=True)  # polar
                pcm = axp.pcolormesh(
                    phi.copy(), r.copy(), Z.copy())
                # color='red',# alpha=0.75
                axp.set_axis_off()
                axp.set_aspect('equal')
                cnt = axp.contour(
                    phi.copy(), r.copy(), Z.copy(),
                    # cmap = plt.cm.gray_r,
                    color='red',
                    levels=levels)

                clev = cnt._levels
                # rx, ry = circle()
                tcolors = cnt.tcolors
                for i in range(len(tcolors)):
                    cc = tcolors[i][0][0:3]
                    if levels is None or ip == len(pole) - 1:
                        x0, y0 = 1.3, 0.8 - i * 0.2
                        r0, t0 = cart2polar(x0, y0)
                        axp.plot(t0, r0, marker='o',
                                 ls='None', mec='black',
                                 mfc=cc,  # color=cc
                                 ms=7. / len(pole),
                                 markeredgewidth=0.01 / len(pole))
                        x2, y2 = 1.40, 0.8 - i * 0.2 - 0.05
                        r2, t2 = cart2polar(x2, y2)
                        axp.text(x=t2, y=r2,
                                 s='%4.2f' % (clev[i]),
                                 fontsize=6. * fact / len(pole))

                    # pole figure indices
                    x3, y3 = 0.4, -1.18
                    r3, t3 = cart2polar(x3, y3)
                    if self.csym == 'hexga' or self.csym == 'trigo':
                        axp.text(
                            x=t3, y=r3,
                            s='(%1i%1i%1i%1i)' %
                              (p__[ip][0], p__[ip][1],
                               p__[ip][2], p__[ip][3]),
                            fontsize=8. * fact / len(pole))
                    else:
                        axp.text(
                            x=t3, y=r3,
                            s='(%1i%1i%1i)' %
                              (pole[ip][0],
                               pole[ip][1],
                               pole[ip][2]),
                            fontsize=8. * fact / len(pole))

                    ## RD and TD indication
                    x4, y4 = -0.05, 1.05
                    r4, t4 = cart2polar(x4, y4)
                    axp.text(x=t4, y=r4, s=iy,
                             fontsize=6. * fact / len(pole))
                    x5, y5 = 1.05, 0.
                    r5, t5 = cart2polar(x5, y5)
                    axp.text(x=t5, y=r5, s=ix,
                             fontsize=6. * fact / len(pole))
                    axp.set_axis_off()

            ## save figures
            fig.savefig('pcm.pdf')
            fig.savefig('pcm.eps')
            for i in range(len(fig.axes)):
                extent = fig.axes[i].get_window_extent()
                extent = extent.transformed(
                    fig.dpi_scale_trans.inverted())
                fig.savefig(
                    'each_axpcm_%s.pdf' % str(i).zfill(2),
                    bbox_inches=extent.expanded(1.1, 1.1))
                fig.savefig(
                    'each_axpcm_%s.eps' % str(i).zfill(2),
                    bbox_inches=extent.expanded(1.1, 1.1))

            ##
            logging.debug("A figure's been saved to pcm.pdf and .eps")
            logging.debug('each_axpcm_00.pdf and .eps')

            fig.clf()
            # return phi, r, Z
        pass

    def ipf(self, pole=None, ifig=4, mode='dot', deco=True, **kwargs):
        """
        Given the pole plot the dot inverse pole figure.
        **The contour version of the inverse pole
        figure is not prepared yet.

        ---------
        Arguments
        ---------
        pole  = [1,0,0]
        ifig  = 1
        mode  = 'dot', 'contour'
        deco  = True
        **kwargs - key-worded argument passed to plots.
        """
        if type(pole) == type(None):
            raise IOError('miller index of the pole should be given')

        if mode == 'dot':
            temp = []
            for i in range(len(self.gr)):
                if deco == True and i == 0:
                    _deco_ = True
                else:
                    _deco_ = False

                tm, fig = self.dotplot(proj='ipf', csym=self.csym,
                                       agrain=self.gr[i],
                                       pole=pole, ifig=ifig, deco=_deco_,
                                       **kwargs)
                temp.append(tm)
            # self.dotplot(proj='ipf',
            # agrain=self.gr[i], pole=[0,1,0], ifig=5)
            # self.dotplot(proj='ipf',agrain=self.gr[i],
            # pole=[0,0,1], ifig=6)
            return temp
        elif mode == 'contour':
            fig = plt.figure(ifig)
            ipf_ax = fig.add_subplot(111)
            for i in range(len(self.gr)):
                pass
            pass
        else:
            raise IOError('Unexpected model for ipf')

    def core(self, pole=None, proj='pf', csym=None, ssym=None,
             agrain=None, isym=True, cdim=[1, 1, 1],
             cang=[90., 90., 90.],
             equivp=None):
        """
        --------------------------------
        The core of the polefigure class
        --------------------------------

        It is the engine for plotting regular and inverse
        pole figures by generating projected
        cartensian coordinates for the 3D vectors
        onto the pole figure sphere (stereographic sphere).
        One can directly plot pole figure on a xy plane.

        Provided the miller indices,
        1. Calculates the crystallographically equivalent poles
        2. Maps the miller indices for a grain
           into a vector in sample axes
        3. Returns mapped poles as raw and as projected (xy)

        ---------
        Arguments
        ---------
          pole = None
          proj = 'pf'
          csym = 'cubic', 'hexag'
          ssym = None,
          agrain = None
          isym = True
          cdim = [ 1., 1., 1.]
          cang = [90.,90.,90.]
          equivp = None  #crystallographically equivalent pole
        """
        xy = [];
        POLE = []
        if csym != 'cubic' and csym != 'hexag' and csym != 'None' \
                and csym != 'centro':
            raise IOError("Other symmetries than cubic" + \
                           " or hexag or 'None'" + " nor 'centro' is" + \
                           " not prepared yet")
        if proj != 'pf' and proj != 'ipf':
            raise IOError("Other modes of projection than pf and" + \
                           "ipf is not prepared yet")

        if type(agrain) == type(None):
            raise IOError("A grains must be given to the method")
        if type(pole) == type(None):
            raise IOError("Pole must be given to core")

        if type(pole).__name__ == 'ndarray':
            pass
        elif type(pole).__name__ == 'list':
            pole = np.array(pole)

        else:
            raise IOError('Unexpected type of the pole argument')

        temp = pole.copy()
        del pole;
        pole = temp

        ##### calculates crystallographically equivalent poles #####
        ## pole figure projection
        if proj == 'pf':
            if isym:
                npeq = equivp  # xtallographically equiv. poles
            else:
                npeq = [pole]
            nit = len(npeq)
            nppp = []
            for i in range(nit):
                nppp.append(npeq[i])
                nppp.append(npeq[i] * -1)
            npeq = np.array(nppp)
            for i in range(len(npeq)):
                for j in range(len(npeq[i])):
                    if abs(npeq[i, j]) < 1e-9:
                        npeq[i, j] = 0.

        ## inverse pole figure
        elif proj == 'ipf':
            # if abs(pole[0]**2+pole[1]**2+pole[2]**2-1)>0.1**6:
            #     print "The pole must be one of principal axes of"\
            #         " the sample"
            #     print "It should be among [1,0,0], [0,1,0], [0,0,1]"
            #     print "current pole is as below\n", pole
            #     raw_input()
            #     raise IOError
            npeq = [pole]  ## Note it is sample axis vector!
            # equivalent pole calculatation is deffered to
            # next for in block
        ## unexpected proj argument
        else:
            logging.debug("It should be either pf or ipf"); raise IOError

        t0 = time.time()
        t_agr2pol = 0.
        t_proj = 0.

        for ip in range(len(npeq)):
            ## 'pf':  converts ca pole to sa pole
            ## 'ipf': converts sa pole to ca pole
            t_1 = time.time()
            if proj == 'ipf':
                p = agr2pol(agrain=agrain, miller=npeq[ip], proj=proj)
            elif proj == 'pf':
                ag = agrain[:3].copy()
                _np_eq_ = npeq[ip]
                if i_for:
                    p = agr2pol_f(ag, _np_eq_)
                else:
                    p = agr2pol(ag, _np_eq_, proj=proj)
            # p = agr2pol(agrain=agrain, miller=npeq[ip], proj=proj)

            t_agr2pol = t_agr2pol + time.time() - t_1
            if proj == 'pf':  # p is in sa
                ## if a pole is toward the north pole of the unit circle,
                # if p[2]>0: pass
                # else:
                POLE.append(p)
                t_1 = time.time()

                if i_for:
                    xy.append(proj_f(p[:3]))
                else:
                    xy.append(projection(pole=p))

                t_proj = t_proj + time.time() - t_1
            elif proj == 'ipf':  # p is in ca
                ## calculates equivalent p by
                ## applying symmetry operations
                if isym:
                    """
                    Sample axis is now referred to crystal
                    coordinate system. That vector has to
                    be mutiplicated by symmetry operations.
                    """
                    npoles = __equiv__(
                        miller=p, csym=csym,
                        cdim=cdim, cang=cang)
                    temp = []
                    for i in range(len(npoles)):
                        temp.append(npoles[i])
                        temp.append(npoles[i] * -1)
                    temp = np.array(temp)
                    npoles = temp

                else:
                    npoles = [p]
                for npp in range(len(npoles)):

                    if i_for:
                        prj_xy = proj_f(npoles[npp][:3])
                    else:
                        prj_xy = projection(pole=npoles[npp])
                    xy.append(prj_xy)
                    POLE.append(npoles[npp])
                pass  # if over 'pf' or 'ipf'
            pass  # End of for loop over ipf

        # print 'Elapsed time for agr2pol ', t2s(t_agr2pol)
        # print 'Elapsed time for proj ', t2s(t_proj)

        return xy, POLE

    def cells(self, pole=[1, 0, 0], ifig=None, dm=15., dn=15.,
              csym=None, cang=[90., 90., 90.], cdim=[1., 1., 1.],
              proj='pf', pole_mode='sys', poles_gr=None):
        """
        Creates cells whose resolutioin is mgrid * ngrid.
        Given the delta m and delta n (dm, dn), each pole's
        weight is assigned to a cell which braces it.
        Plots the cell's weight and returns the cell in array.

        ---------
        Arguments
        ---------
        pole = [1,0,0]
        ifig = None
        dm   = 7.5.
        dn   = 7.5.
        csym = None
        cang = [90.,90.,90.]
        cdim = [1.,1.,1.,]
        proj = 'pf'
        pole_mode='sys', 'indv'
        poles_gr = None, [] array shape: (ngr, 3)
        """
        ## Frequently used local functions or methods
        pol = []  # a group of pole points.

        ## npeq calculations
        if pole_mode == 'sys':
            npeq = __equiv__(
                miller=pole, csym=csym,
                cdim=cdim, cang=cang)

        t0 = time.time()

        # is_joblib=False
        if is_joblib:
            ## This loop is quite extensive and slow.
            rst = Parallel(n_jobs=3)(delayed(core)(
                self, pole=pole, proj=proj, csym=csym,
                cang=cang, cdim=cdim, equivp=npeq,
                agrain=self.gr[i]) for i in range(len(self.gr)))

            for i in range(len(rst)):
                xy, p = rst[i]
                # p is n equivalent poles
                p = np.array(p)
                for j in range(len(p)):
                    # make it postive.
                    if p[j][2] < 0: p[j] = p[j] * - 1
                    ## phi, and theta in radian
                    x, y = cart2sph(p[j])
                    pol.append([x, y, self.gr[i][3]])  # phi, theta, intensity

        elif not is_joblib:
            for i in range(len(self.gr)):
                ## Either systematic representative
                ## poles or individual pole for
                ## each and every grain.
                if pole_mode == 'sys':
                    xpole = pole
                elif pole_mode == 'indv':
                    xpole = poles_gr[i]
                    npeq = [xpole, xpole * -1]  ## two coaxial poles
                xy, p = self.core(
                    pole=xpole, proj=proj, csym=csym,
                    agrain=self.gr[i], cang=cang,
                    cdim=cdim, equivp=npeq)
                # p is n equivalent poles
                p = np.array(p)
                for j in range(len(p)):
                    # make it postive.
                    if p[j][2] < 0: p[j] = p[j] * - 1
                    ## phi, and theta in radian
                    x, y = cart2sph(p[j])
                    pol.append([x, y, self.gr[i][3]])  # phi, theta, intensity

        logging.debug('Elapsed time for self.core:', t2s(time.time() - t0))

        t0 = time.time()

        ## weight normalization
        pol = np.array(pol).T
        pol[2] = pol[2] / pol[2].sum()  ## for normalization that follows.
        pol = pol.transpose()

        # grid number along phi   axis (x) is
        # referred to as m
        # grid number along theta axis (y) is
        # referred to as n

        # m x n grid;  m and n's range
        # m is in [-pi.,pi];  n is in [0, pi/2.]


        """
        (ngrid) (0,pi/2.) (theta, tilting)
        ^
        |
        |
        |
        |
        L_______________>  (mgrid) (-pi,+pi) (phi, rotation)

        phi = (-pi,+pi) + dm/2.
        theta = (-pi,+pi) + dm/2.
        """
        dmp = dm * pi / 180.
        dnp = dn * pi / 180.

        mgrid = int(360. / dm);
        ngrid = int(90. / dn)
        phi_angle = np.arange(-pi, pi) / dm + dm / 2
        theta_angle = np.arange(0., pi / 2.) / dn + dn / 2
        f = np.zeros((mgrid, ngrid))
        for i in range(len(pol)):
            phi = pol[i][0];
            theta = pol[i][1]
            mi = int((phi + pi) / dmp - 1e-9)  # subtract tiny
            ni = int(theta / dnp - 1e-9)  # subtract tiny
            f[mi, ni] = f[mi, ni] + pol[i][2]

            # if mi<0 or ni<0 :
            #     raise IOError, 'Negative index'
            # elif mi>mgrid or ni>ngrid:
            #     raise IOError, 'Unexpected Error'
            # try:
            # f[mi,ni] = f[mi,ni] + pol[i][2]
            # except:
            #     print phi*180./np.pi, theta*180./np.pi
            #     raise IOError, "Something wrong in "+\
            #         "the index for f[mi,ni]"

        ## ----------------------------------------
        """
        Symmetrization over the cell (f)
        """
        # 1. Symmetrization abount the x-axis
        # 2. Symmetrization about the y-axis
        ## ----------------------------------------

        # #ismooth = False
        # ismooth = True
        # if ismooth:
        #     fpole = 0.
        #     for m in xrange(int(mgrid)):
        #         fpole = fpole + f[m,0]
        #     fnorm = (1.-cos(dn * pi/180. ))*2.*pi
        #     #fnorm = 1.*pi
        #     fpole = fpole / fnorm
        #     for m in xrange(int(mgrid)):
        #         f[m,0] = fpole
        # else:  pass

        # fnorm = dcos(z) * deltx / 360

        """
        Spherical coordinates result in an area element
        dA = sin(n)*dm*dn
        This area element dA is dependent on the tilting (n) grid

        Normalization of pole figure intensity, i.e., f(theta, phi) is:
        1/(2pi) \int f(theta,phi) sin(theta),dphi,dtheta = 1
        """

        z = np.zeros((ngrid + 1,))
        deltz = (pi / 2.) / float(ngrid)
        for i in range(ngrid + 1):
            z[i] = deltz * i

        deltx = 2. * pi / mgrid
        for m in range(mgrid):
            for n in range(int(ngrid)):
                fnorm = (cos(z[n]) - cos(z[n + 1])) * deltx / (2 * pi)
                f[m, n] = f[m, n] / fnorm

        if ifig is not None:
            fig = plt.figure(ifig)
            ax = fig.add_subplot(111)
            ax.plot(f)
            ax.set_ylim(0., )

        nodes = np.zeros((mgrid + 1, ngrid + 1))  # rot, tilting = (azimuth, declination) ...

        ## Assigns the same intensity for the first n rings
        n = 1
        for i in range(n):
            # South-pole : along  n=0 axis ---
            for m in range(len(f) + 1):
                nodes[m, i] = f[:, i].sum() / len(f)
        ## ------------------------------

        # # along n=1 axis ----------------
        # f1 = 0.
        # for m in xrange(len(f)):
        #     f1 = f1 + f[m,1]
        # f1avg = f1 / len(f)
        # for m in xrange(len(f)+1):
        #     nodes[m,1] = f1avg
        # ## ------------------------------

        ph = np.linspace(-pi, pi, mgrid + 1)
        th = np.linspace(0., pi / 2., ngrid + 1)

        regions = np.zeros((mgrid + 1, ngrid + 1, 4, 2))
        for m in range(int(mgrid + 1)):
            for nn in range(int(ngrid)):
                n = nn + 1
                # nodes[m, n]  # <--- currently interesting node
                ## find the nodes' phi and theta
                ## that's phi[m] and theta[n]

                reg_ph = np.zeros((2,));
                reg_th = np.zeros((2,))

                reg_ph[0] = ph[m] - (dm * pi / 180.) / 2.
                reg_ph[1] = ph[m] + (dm * pi / 180.) / 2.
                reg_th[0] = th[n] - (dn * pi / 180.) / 2.
                reg_th[1] = th[n] + (dn * pi / 180.) / 2.

                reg = [[reg_ph[0], reg_th[0]], [reg_ph[0], reg_th[1]], [reg_ph[1], reg_th[0]], [reg_ph[1], reg_th[1]]]

                for i in range(len(reg)):
                    p = reg[i][0]
                    t = reg[i][1]
                    if p > pi:
                        p = p - 2. * pi
                    elif p < -pi:
                        p = p + 2. * pi
                    elif p == pi or p == -pi:
                        logging.debug('Unexpected..');
                        raise IOError
                    if t > pi / 2.:
                        t = t - (dn * pi / 180.)
                        p = p + pi
                        if p > pi:
                            p = p - 2. * pi
                        elif p < -pi:
                            p = p + 2. * pi
                    reg[i][0] = p
                    reg[i][1] = t

                ## for each node, find the 4 adjacent
                ## regions' intensites
                ## and average them out and assign to the node.
                inten = 0.
                for i in range(4):
                    p = reg[i][0]
                    t = reg[i][1]
                    mi = int((p + pi) / (dm * pi / 180.) - 0.1 ** 6)
                    ni = int(t / (dn * pi / 180.) - 0.1 ** 6)
                    if mi < 0 or ni < 0:
                        input('Negative index')
                        raise IOError
                    inten = inten + f[mi, ni]
                nodes[m, n] = inten / 4.

        logging.debug('Elapsed time in self.core after cells:', t2s(time.time() - t0))
        return f, nodes

    def define_axes(self):
        """
        convert three bases vectors (i.e., [1,0,0], [0,1,0], [0,0,1]
        And find to which direction they are pointing at in the final
        pole figure projection.
        """
        self.bases = np.identity(3)
        x, y, z = self.bases[0, :], self.bases[1, :], self.bases[2, :]
        ## Very often, x/y/z are aligned with RD/TD/ND.
        ## Note that self.bases are referenced in the laboratory axes.

    def pf_new(
            self, ifig=None, axs=None,
            poles=[[1, 0, 0], [1, 1, 0]], ix='1', iy='2',
            mode='line',
            dth=10, dph=10, n_rim=2, cdim=None, ires=True, mn=None, mx=None,
            lev_norm_log=True, nlev=7, ilev=1, cmap='magma',
            rot=0., iline_khi80=False):
        """
        New version of pf that will succeed upf.polefigure.pf
        Note that upf.polefigure.pf is deprecated and will be deleted soon.

        Arguments
        ---------
        <ifig> or <axs>
            <ifig> and <axs> should be mutually exclusive.
            It is acceptable for both to be *not* specified.
            However, it is unacceptable for both to be specified.

        <poles>
           For cubics, three digits; for hexagonals four digits
        <ix>, <iy>
           x and y tick labels appended to each pole figure
        <dph>:
            (tilting angle : semi-sphere 0, +90 or full-sphere 0, +180)
        <dth>:
            (rotation angle: -180,+180)

        <n_rim>:
             The number of 'central' rims to be *averaged*.
             For better understandings, see the algorithm notebook
             located in
                ./ipynb/UPF_Algorithm.ipynb

        <cdim>:  crystal dimension
           For cubic, [1,1,1]
           For a particularly AZ31 sheet, it is [1,1,1.6235]
           Users should know what is the lattice dimension for the crystal
           structure of which he/she plots the pole figures.

        <ires>  = True;
           If True, indicate the grid
           If <mode> is 'fill' and ires is True, overlay the resolution
              all over the pole.
           if <mode> is 'line' and ires is True, only the spots lower than
              minimum level of contour is plotted.

        <mn>,<mx>
          minimum and maximum levels of contour
          If not specified, mn and mx is determined using levels of
          calculated contours
          if <mode> is 'fill', <mn> is overriden by the levels of
          calculated contours.

        <lev_norm_log>
           If True, use logarithmic scales. If False, linear scale.
        <nlev> = 7
           Level of iso contour bins.
           The total number of lines will be nlev+1
        <cmap>
           Color map used to color-code the contour levels.
           Refer to 'http://matplotlib.org/users/colormaps.html'
           for more color-map options.
        <iline_khi80> = False
           Whether or not to draw a line of chi=80 in experimental
           pole figure plot. - I usually obtain incomplete pole figure
           upto a tilting <chi> of 80.
        <mode>
           Contour model: 'line' or 'fill'
        <ilev>
           level option: 0 commonly contour levels for all poles generated
                         1 individual levels applied for individual poles

        Returns
        -------
        fig: matplotlib.figure.Figure
        """
        import MP.lib.mpl_lib

        # Turn interactive plotting off
        plt.ioff()

        ## check mutually exclusive arguments (ifig and axs)
        if type(ifig) != type(None) and type(axs) != type(None):
            raise IOError('** Err: ifig and axs are mutually exclusive')

        ##################################################
        ## PF plots for experimental pole figure is
        ## separately conducted by epfplot function
        if type(self.epf).__name__ != 'NoneType':
            logging.debug('Writing Experimental pole figures')
            return self.epfplot(
                ifig=ifig, cmap=cmap, nlev=nlev, mn=mn, mx=mx,
                ix=ix, iy=iy, rot=rot, iline_khi80=iline_khi80)
        ##################################################

        nlev = nlev + 1  ##
        from matplotlib.colors import LogNorm
        miller = poles[::]

        if type(cdim) != type(None): self.cdim = cdim
        ## 4 digits miller indices are used for hexagon and trigo
        if self.csym == 'hexag' or self.csym == 'trigo':
            pole_ = []
            for i in range(len(poles)):
                p = [0, 0, 0]
                p_ = poles[i]
                if len(p_) != 4: raise IOError('4 digits should be given')
                p[2] = p_[3]
                p[0] = p_[0] - p_[2]
                p[1] = p_[1] - p_[2]
                pole_.append(p)
            poles = pole_[::]

        tiny = 1.e-9
        N = []
        t0 = time.time()
        # is_joblib=False ## Debug
        if is_joblib and len(poles) > 1:
            rst = Parallel(n_jobs=len(poles))(
                delayed(cells_pf)(
                    pole_ca=poles[i], dth=dth, dph=dph,
                    csym=self.csym, cang=self.cang, cdim=self.cdim,
                    grains=self.gr, n_rim=n_rim) for i in range(len(poles)))
            for i in range(len(rst)):
                N.append(rst[i])
        else:
            for i in range(len(poles)):
                N.append(cells_pf(
                    pole_ca=poles[i], dth=dth, dph=dph,
                    csym=self.csym, cang=self.cang,
                    cdim=self.cdim, grains=self.gr,
                    n_rim=n_rim))

        et = time.time() - t0
        # uet(et, head='Elapsed time for calling cells_pf')
        # logging.debug()

        x_node = np.arange(-180., 180. + tiny, dth)
        y_node = np.arange(0., 90. + tiny, dph)
        XN, YN = np.meshgrid(x_node, y_node)

        # --------------------------------------------------#
        ## plotting / resolution
        nm = int((360.0 - 0.) / dth)
        nn = int((180. - 90.) / dph)
        theta = np.linspace(pi, pi / 2., nn + 1)
        phi = np.linspace(0., 2. * pi, nm + 1)
        r = np.sin(theta) / (1 - np.cos(theta))
        R, PHI = np.meshgrid(r, phi)
        PHI = PHI + rot  ## default: rot=0.
        x = R * np.cos(PHI);
        y = R * np.sin(PHI)

        nArray = np.array(N)
        xyCoords = np.array([x, y])

        # print 'nArray.shape:',nArray.shape

        mns, mxs, indices_mx = self.calcMXN(nArray, mx, mn, mode, ilev)

        # ## normalized color map
        # mpl_cmap, cm_func, cm_norm = MP.lib.mpl_lib.norm_cmap(
        #     mx=mx,mn=mn,cm_name=cmap,inorm=True)

        if type(ifig) == type(None):
            fig = plt.figure(figsize=(3.3 * len(poles), 3.0))
        else:
            fig = plt.figure(ifig, figsize=(3.3 * len(poles), 3.0))

        ##
        axs = []
        for i in range(len(poles)):
            _ax_ = fig.add_subplot(1, len(poles), i + 1)
            axs.append(_ax_)
        plt.subplots_adjust(left=0, right=0.8)

        for i in range(len(poles)):
            if lev_norm_log:
                ## To prevent log (0) -> np.nan
                ## hardwire minimum value
                if mns[i] == 0: mns[i] = 0.5
                levels = np.logspace(
                    np.log10(mns[i]), np.log10(mxs[i]), nlev)
                norm = LogNorm()
            else:
                levels = np.linspace(mns[i], mxs[i], nlev)
                norm = None

            import matplotlib.cm
            cmap_mpl = matplotlib.cm.get_cmap(cmap)
            color_mapping = matplotlib.cm.ScalarMappable(
                norm=norm, cmap=cmap_mpl)

            if mode == 'line':
                func = axs[i].contour
            elif mode == 'fill':
                func = axs[i].contourf

            ## contour plot
            nArray[i][np.isnan(nArray[i])] = 0.
            nArray[i][nArray[i] <= 0] = 1e-4
            cnts = func(x, y, nArray[i], levels=levels,
                        cmap=cmap, norm=norm, zorder=10)

            ## x, y coordinates of maximum intensity in grid
            i0, j0 = indices_mx[i]
            mx_coord_x = x[i0, j0]
            mx_coord_y = y[i0, j0]

            if mode == 'line':
                axs[i].plot(mx_coord_x, mx_coord_y, '+', mew=2,
                            color=color_mapping.to_rgba(levels[-1]))

            if ires:  # and mode!='fill':
                xs = [];
                ys = []
                filt = nArray[i, :, :] < levels[0]
                filt[0, 1:] = False
                filt[1:, 0] = False
                xs = x[filt];
                ys = y[filt]
                if len(xs) > 0:
                    axs[i].plot(
                        xs, ys, 'k.',
                        alpha=0.17 * len(poles),
                        markersize=2.0)
            # if ires and mode=='fill': ## overlay the resolution
            #     axs[i].plot(x,y,'k+',
            #                 alpha=0.17*len(poles),
            #                 markersize=2.0,zorder=100)

            deco_pf(axs[i], cnts, miller[i], 0,
                    iskip_last=False, ix=ix, iy=iy, mode=mode)
        return fig
        # --------------------------------------------------#

    def calcMXN(self, nArray=None, mx=None, mn=None, mode='line', ilev=0):
        """
        Arguments
        ---------
        nArray: ndarray that contains all pole figure nodes.
        mx
        mn
        mode
        ilev  : option (0: common mx and mn, 1: individual mx and mn)

        Returns
        -------
        mns, mxs, indices_mx
        """
        npole = nArray.shape[0]
        mxs = np.zeros(npole)
        mns = np.zeros(npole)

        if ilev == 0:
            ## determine maximum and minimum levels.
            if type(mx) == type(None): mx = nArray.flatten().max()
            if mx > 100: mx = 99.

            if type(mn) == type(None) and mode != 'fill':
                mn = 0.5
            elif type(mn).__name__ == 'float':
                pass
            else:
                mn = nArray.flatten().min()
            ## commonly assigned
            mxs[:] = mx
            mns[:] = mn
        elif ilev == 1:
            for ipole in range(npole):
                if type(mx) == type(None):
                    mx_ = nArray[ipole].flatten().max()
                else:
                    mx_ = mx * 1.

                if mx_ > 100: mx_ = 99.

                if type(mn) == type(None) and mode != 'fill':
                    mn_ = nArray[ipole].flatten().min()
                elif type(mn) == type(None) and mode == 'fill':
                    mn_ = nArray[ipole].flatten().min()
                else:
                    mn_ = mn * 1.

                mxs[ipole] = mx_
                mns[ipole] = mn_
        else:
            raise IOError

        ## Find coordinates of 'maximum' intensity
        indices_mx = []
        for ipole in range(npole):
            indx_mx = np.unravel_index(nArray[ipole, :, :].argmax(),
                                       nArray[ipole, :, :].shape)
            ## indx_mx = np.argmax(nArray[ipole,:,:],axis=1)
            indices_mx.append(indx_mx)
        indices_mx = np.array(indices_mx)
        return mns, mxs, indices_mx

    def dotplot(self, pole=None, ifig=None, npole=1,
                ipole=1,
                proj='ipf', csym='cubic',
                ssym='tric', agrain=None,
                isym=True, irdc=True,
                cdim=[1., 1., 1.],
                cang=[90., 90., 90.], mode='trace',  # or None
                deco=True,
                **kwargs
                # alpha=1.0, color='k',
                ):
        """
        Plots the individual poles. No contour.

        ---------
        Arguments
        ---------
        pole  = None
        ifig  = None
        npole = 1
        ipole = 1
        proj = 'pf' or 'ipf'
        csym = 'cubic' Crystal symmetry
        ssym = 'tric'  Sample symmetry
        agrain = None
        isym = True:   Symmetry operation to the pole
        irdc = True:   Reduction of inverse pole figure region
        cdim = [1., 1., 1.]
        cang = [90.,90.,90.]
        mode = 'trace'
        deco = True
        **kwargs : matplotlib pyplot key-worded arguments
        """
        if pole is None:
            logging.debug("Pole must be given")
            raise IOError

        ## If a grain is assigned,
        ## Operations are carried on it, instead of on
        ## Global grains.
        # if color==None:
        #     print 'Default color is black'; color='k'
        if type(agrain) != type(None):
            "* a grains' euler angle is given"
            "* dot plotting is performed on a grain"
            agr = np.array([agrain])  # p1,p,p2, inten
        else:
            logging.debug('Agrain must be asigned')
            raise IOError
        # print "------------------------------------"
        XY = [];
        polz = []
        ## Draws the boundary of inverse and normal pole figures
        if type(ifig) == type(None):
            fig = None
            pass
        else:
            ## if plt.figure(ifig) is present, use it
            ## otherwise, create one.

            fact = 3.
            figsize = (npole * fact, 1. * fact)
            if plt.fignum_exists(ifig):
                fig = plt.figure(ifig)
            else:
                fig = plt.figure(ifig, figsize=figsize)

            ax = fig.add_subplot(1, npole, ipole)

            if deco:
                ax.set_axis_off();
                ax.set_aspect('equal')
                ax.text(x=-0.08, y=-0.07, s='(100)', fontsize=4. * fact, transform=ax.transAxes)
                ax.text(x=0.7, y=-0.07, s='(110)', fontsize=4. * fact, transform=ax.transAxes)
                ax.text(x=0.65, y=0.8, s='(111)', fontsize=4. * fact, transform=ax.transAxes)
                # ax.text(x=0.5, y=-1.05, s='(%i%i%i)'%
                #         (pole[0],pole[1],pole[2]), fontsize=3. * fact,transform=ax.transAxes)
            if proj == 'pf' and deco:
                rx, ry = __circle__(center=[0, 0], r=1)
                ax.plot(rx, ry, color='grey')
                ax.set_xlim(-1.2, 1.2);
                ax.set_ylim(-1.2, 1.2)
            elif proj == 'ipf' and deco:
                # in ipfline module, projection can be
                # changed into equal area type
                cxy = ipfline(center=[0, 0], csym=csym)
                ax.plot(cxy[0], cxy[1], color='grey', alpha=0.1)

        ## add poles
        for i in range(len(agr)):
            ### crystallographically equivalent poles are calculated.
            agr[i][0] = agr[i][0] - 90.
            npeq = __equiv__(
                miller=pole, csym=csym, cdim=cdim, cang=cang)
            ### -----------------------------------------------------
            xy, POLE = self.core(
                pole=pole, proj=proj, csym=csym,
                agrain=agr[i], isym=isym, cdim=cdim,
                cang=cang, equivp=npeq)

            ##### POLE FIGURE PROJECTIN #####
            if proj == 'pf':
                for j in range(len(xy)):
                    if POLE[j][2] <= 0:  # north pole is [0,0,1]
                        XY.append(
                            [xy[j][0], xy[j][1], agr[i][3]]
                        )  # x,y, intensity

            ##### INVERSE POLE FIGURE PROJECTION #####
            elif proj == 'ipf':
                for j in range(len(xy)):
                    r = np.sqrt(xy[j][0] ** 2 + xy[j][1] ** 2)
                    if r > 1.0:
                        pass
                    else:
                        ### reduced region filter
                        ### must be applied here.
                        # y must be positive.
                        tiny = 0.
                        phi = math.atan2(xy[j][1], xy[j][0])
                        phi = phi * 180.0 / math.pi
                        if 0. - tiny < phi < 45.0 + tiny:
                            ## another filter
                            ## 2nd and 3rd triangles
                            # a = POLE[j][0]; c = POLE[j][2]
                            a, b, c = invproj(x=xy[j][0], y=xy[j][1])
                            # print 'a,b,c= ',a,b,c
                            # print 'atan2(a,c)',
                            # math.atan2(a,-c)*180./math.pi
                            if math.atan2(
                                    a, -c) * 180. / math.pi < 45.0 + tiny:
                                XY.append([xy[j][0], xy[j][1],
                                           agr[i][3]])
                                # x, y, intensity of the grain

        if type(ifig) != type(None):
            # plotting --------------------
            try:
                XY = np.array(XY)
                xxx = XY.copy().transpose()[0]
                yyy = XY.copy().transpose()[1]
                ax.plot(xxx, yyy, ls='None', **kwargs)
                if proj == 'pf':
                    ax.set_xlim(-1.1, 1.1)
                    ax.set_ylim(-1.1, 1.1)
                elif proj == 'ipf':
                    ax.set_xlim(-0.02, 0.5)
                    ax.set_ylim(-0.02, 0.5)
            except:
                if len(XY) == 0:
                    pass
                    #   print 'Empty XY is returned'
                else:
                    raise IOError("Unexpected Error")
                    # raw_input()
                    # -------------------------------
        return np.array(XY), fig


def cells_pf(
        pole_ca=[1, 0, 0], dph=7.5,
        dth=7.5, csym=None, cang=[90., 90., 90.],
        cdim=[1., 1., 1.], grains=None, n_rim=2):
    """
    Creates cells whose resolutioin is nphi, ntheta
    Given the delta x and delt y (dm, dn), each pole's
    weight is assigned to a cell which braces it.
    Plots the cell's weight and returns the cell in array.

    ---------
    Arguments
    ---------
    pole_ca = [1,0,0]
    dph  = 7.5. (tilting angle : semi-sphere 0, +90 or full-sphere 0, +180)
    dth  = 7.5. (rotation angle: -180,+180)
    csym = None
    cang = [90.,90.,90.]
    grains = None, [] array shape: (ngr, 3)
    n_rim=2
    """
    tiny = 1e-9
    ## Find poles in sa [pole_sa]
    # phs, ths, wgts = [], [], []

    ## Set of equivalent vectors based on crystal symmetry
    p0 = __equiv__(miller=pole_ca, csym=csym,
                   cdim=cdim, cang=cang)
    poles_ca = np.zeros((len(p0) * 2, 3))
    poles_ca[:len(p0), :] = p0[:, :]
    poles_ca[len(p0):, :] = -p0[:, :]
    poles_ca = poles_ca / np.sqrt((poles_ca ** 2).sum())

    nx, ny = len(grains), len(poles_ca)
    poles_sa = np.zeros((nx, ny, 3))
    poles_wgt = np.zeros((nx, ny))

    # i_for=False # debug
    if i_for:
        poles_sa, poles_wgt = gr2psa(
            ngr=len(grains), grains=grains,
            npol=len(poles_ca), poles_ca=poles_ca)  # np.array(poles_ca))
    else:
        for i in range(len(grains)):
            phi1, phi, phi2, wgt = grains[i]
            ## arg = euler_f(2,phi1,phi,phi2,np.zeros((3,3))) ## ca<-sa
            ## amat = arg[-1]
            amat = euler(phi1, phi, phi2, a=None, echo=False)
            for j in range(len(poles_ca)):
                poles_sa[i, j, :] = np.dot(amat.T, poles_ca[j])
                poles_wgt[i, j] = wgt

    poles_sa = poles_sa.reshape((len(grains) * len(poles_ca), 3))
    poles_wgt = poles_wgt.reshape((len(grains) * len(poles_ca)))

    ## Full Sphere
    x = np.arange(-180., 180. + tiny, dth)
    y = np.arange(0., 180. + tiny, dph)
    nx, ny = int(360. / dth), int(180. / dph)
    f = np.zeros((nx, ny))

    ## Semi Sphere
    x_node = np.arange(-180., 180. + tiny, dth)
    y_node = np.arange(0., 90. + tiny, dph)
    nx_node = len(x_node);
    ny_node = len(y_node)
    nodes = np.zeros((nx_node, ny_node))

    f = pole2f(poles_sa, poles_wgt, dth, dph, f.copy())

    ## Normalization (m.u.r)
    fsum = f[:, :int(ny / 2)].flatten().sum()
    z = np.zeros((ny + 1))
    for i in range(ny + 1): z[i] = np.pi / float(ny) * i
    dx_ = 2. * np.pi / nx
    dcosz = -np.diff(np.cos(z))
    fnorm = dcosz * dx_ / (2 * np.pi)
    f = f / fnorm / fsum

    ## Extension of f_bounds - see algorithm ipynb
    f_bounds = np.zeros((nx + 2, ny + 2))
    f_bounds[1:-1, 1:-1] = f[:, :]
    f_bounds[-1, 1:-1] = f[0, :]
    f_bounds[0, 1:-1] = f[-1, :]
    f_bounds[1:-1, 0] = f[:, -1]
    f_bounds[0, 0] = f_bounds[0, -2]
    f_bounds[-1, 0] = f_bounds[-1, -2]
    f_bounds[:, -1] = f_bounds[:, 1]

    ## Use average of the four adjacent neighbouring nodes
    for i in range(len(nodes)):
        for j in range(len(nodes[i])):
            nodes[i, j] = (f_bounds[i:i + 2, j:j + 2]).sum() / 4.

    ## Centeral region is using an avergage around the rim
    for i in range(n_rim):
        nodes[:, i] = (nodes[:, i].sum()) / len(nodes[:, i])

    XN, YN = np.meshgrid(x_node, y_node)
    return nodes


def __equiv__(miller=None, csym=None,
              cdim=[1., 1., 1.], cang=[90., 90., 90.]):
    """
    Provided the miller indices,
    Crystallographically equivalent and only unique
    vectors are returned.

    ---------
    Arguments
    ---------
    miller = None  , e.g. [1,1,1]
    csym   = 'cubic'
    cdim   = [ 1, 1, 1]
    cang   = [90,90,90]
    """
    start = time.time()
    from .sym import cv
    from .sym import cubic, hexag
    # from sym_cy import cubic, hexag
    from . import sym  # python compiled
    # import sym_cy #cython compiled
    # from sym.py cvec, cubic, and hexgonal modules are brought in
    if type(miller) == type(None): raise IOError("Miller index should be given")

    vect = np.array(miller)
    norm = 0.;
    sneq = []
    temp = vect.copy()
    # start = time.time()
    if csym == 'cubic':
        # H = sym_cy.cubic(1)  #cython compiled
        H = sym.cubic()  # operators
        for i in range(len(H)):
            sneq.append(np.dot(H[i], vect))
    elif csym == 'hexag':
        # H = sym_cy.hexag(1) #cython compiled
        H = sym.hexag()  # operators
        v = cv(pole=vect, cdim=cdim, cang=cang)
        sneq = np.tensordot(H, vect, axes=[-1, 0])
        # for i in xrange(len(H)):
        #     sneq.append(np.dot(H[i], v))
    elif csym == 'None':
        # H = [np.identity(3)]
        sneq = [vect]
    elif csym == 'centro':
        sneq = [vect, -vect]
    else:
        logging.debug('Given symmetry, %s is not prepared' % csym)
        input('Enter to raise an error and quits the job')
        raise IOError

    # print 'elapsed time during v calculation: %8.6f'%
    # (time.time()-start)
    #####-------------------------------
    # start = time.time()
    stacked = []  # empty unique vectors
    # is cH in the already existing stacked list?
    # yes: pass
    # no : add

    ## filtering the sneq under whether or not it is unique
    for i in range(len(sneq)):
        cH = sneq[i].copy()  # current vector
        if __isunique__(a=cH, b=stacked):
            stacked.append(cH)

    ## if v[2] is minus, mutiply minus sign to the vector.
    for i in range(len(stacked)):
        if stacked[i][2] < 0:
            stacked[i] = stacked[i] * -1
    # print 'elapsed time during the rest: %8.6f'%
    # (time.time()-start)
    return np.array(stacked)


### Excutes the module in the command line with arguments and options
def main(filename, pfmode, gr, csym):
    """
    plots the (100),(110) and (111) pole figures of cubic crystal.
    (to be including hexagonal)

    Arugments
    =========
      filename = '00500.cmb' for instance.
      csym : crystal symmetry 'cubic', 'hexag'
      pfmode ='contour','contourf', 'dot'
      """
    if gr is not None:
        a = polefigure(grains=gr, csym=csym)
    else:
        a = polefigure(filename=filename, csym=csym)
    a.pf(mode=pfmode)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import getopt, sys

    ## arguments ------------------------------- ##
    try:
        opts, args = getopt.getopt(
            sys.argv[1:], 'm:i:o:c:s')  # , 'output=', 'csym='])

    except getopt.GetoptError as err:
        logging.debug(str(err))
        sys.exit(2)

    ## ----------------------------------------- ##
    ## default options
    ishow = False
    mode = 'contourf'
    csym = 'cubic'
    outputfile = 'temp.pdf'

    for o, a in opts:
        if o in '-i':
            inputfile = a
        elif o in '-o':
            outputfile = a
        elif o in '-c':
            csym = a
        elif o in '-m':
            mode = a
        elif o in '-s':
            ishow = True
        else:
            assert False, 'unhandled option'

    main(filename=inputfile, csym=csym, pfmode=mode)
    plt.gcf().savefig(outputfile)
    if ishow: plt.show()


def parse_epf(fn, n_unit=79):
    """
    Parse popLA-generated epf files

    Arguments
    ---------
    fn
    n_unit (the number of lines in an epf file)
    """
    nline_each_block = n_unit
    with open(fn, 'r') as fo:
        string = fo.read()
        lines = string.split('\n')[:-1]
        nl = len(lines)

    logging.debug('# of lines : %i' % nl)
    logging.debug('# of blocks: %i' % (float(nl) / float(n_unit)))

    nb = int(float(nl) / float(n_unit))

    blocks = []
    for i in range(nb):
        i0 = n_unit * i + 1
        i1 = n_unit * (i + 1)

        l = lines[i0:i1]

        l[0] = l[0][1:]
        # print '*'*10
        # print l[0]
        # print l[-1]

        b = ''
        # print 'number of lines:', len(l)
        for j in range(len(l)):
            b = '%s%s\n' % (b, l[j])
        blocks.append(b)

    return blocks
