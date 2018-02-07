"""
Sampling a set of discrete orientations that represent the
crystallographic texture given as the crystallographic
orientation distribution.


Youngung Jeong
youngung.jeong@gmail.com
"""

###
## Representative volume element maker
###

"""
cmb module calculates represetnative volume element, also
known as the population of discrete grains.

check if update is followed by submodule?
check if update is followed by cloning?
"""

import matplotlib.pyplot as plt
import numpy as np
import math
rand = np.random.rand
randi = np.random.random_integers

def steglich_format(filename=None):
    """
    2011-Sept-21
    Convert Dr. Steglich's COD format into LABOTEX's

    Reference file: 'ODF numerisch.txt'

    Argument
    --------
    filename
    """
    f = open(filename, 'r')
    contents = f.read()

    ## Assumes that descrete COD is provided by slice of sections
    ## that are perpendicular to phi axis
    blocks = contents.split('Phi1=')
    header = blocks[0]
    planes = blocks[1:]

    axis_p1 = []
    axis_P = []
    axis_p2 = []
    cod = []

    for i in xrange(len(planes)): #each block of phi=constant plane
        clines = planes[i].split('\n')
        block = clines[1:][:-1:] #tail off
        block = np.array(block)
        dum = []
        for i in xrange(len(block)): #phi2
            if i!=0 and len(block[i]) > 3: #PHI
                dum.append(
                    map(float,block[i].split()[1:])) #remove the first row
        dum = np.array(dum) # dum: (phi2, PHI)
        dum = dum.T         # dum: (PHI, phi2)
        # dum = dum[0:]
        dum = dum.tolist()  # make numpy array into list type
        cod.append(dum) # cod: (phi1, PHI, phi2)
        pass

    rst = np.zeros((len(cod), len(cod[0]), len(cod[0][0])))
    for i in xrange(len(cod)): #phi1
        for j in xrange(len(cod[i])): #PHI
            for k in xrange(len(cod[i][j])): #phi2
                rst[i][j][k] = cod[i][j][k]

    print 'rst shape:', rst.shape

    ## write this into LABOTEX descrete COD format file
    ##  phi1 phi phi2 COD
    ##   0   0    0    0.002
    ##   5   0    0    0.012
    ##   ..
    ## 360   0    0    0.023
    ##   0   5    0    0.100
    ##   5   5    0    0.123
    ##   ..
    ##   0   0    5    0.603

    # permute the rst(phi1, phi, phi2) -> temp(phi, phi2, phi1)
    temp = np.transpose(rst, (1,2,0))
    print 'temp shape:', temp.shape
    fout = open('%s_labo.txt'%filename.split('.')[0], 'w')
    fout.writelines('%s %s %s %s \n'%('PHI1','PHI2','PHI', 'COD'))
    for i in range(len(temp)): #phi
        for j in range(len(temp[i])): #phi2
            for k in range(len(temp[i][j])): #phi1
                fout.writelines(
                    '  %6.2f  %6.2f  %6.2f  %12.7e\n'%(
                        k*5., j*5., i*5., temp[i][j][k]))
    return rst


def randomGrain(phi1,phi2,phi):
    """
    Generate a single 'randomly' oriented
    grain and return its 'rotation' matrix

    Argument
    --------
    ## maximum bounds of the thee euler angles.
    phi1
    phi2
    phi
    """
    from euler import euler
    cp1 = rand() * phi1  #phi1
    cp2 = rand() * phi2  #phi2
    cp = rand() # 0~1
    cp = math.acos(cp) * 180./math.pi
    return euler(ph=cp1,th=cp,tm=cp2,echo=False)

def random(phi1=90, phi2=90, phi=90,
           ngrain=1000, iplot=False,
           filename=None):
    """
    Random iostropic texture maker
    Bunge convention

    gr consists of phi1, phi, phi2, and the volume fraction in order.

    Arguments
    =========
    phi1     = 90
    phi2     = 90
    phi      = 90
    ngrain   = 1000
    iplot    = False
    filename = None

    Return
    ------
    gr
    """
    print 'phi1, phi, phi2', phi1, phi, phi2
    gr = np.zeros((ngrain,4))
    for i in range(ngrain):
        cp1 = rand() * phi1  #phi1
        cp2 = rand() * phi2  #phi2
        cp = rand() # 0~1
        cp = math.acos(cp) * 180./math.pi
        # if phi==180:
        #     if randi(0,1)==0: cp = cp
        #     else: cp = 180 - cp
        # elif phi==90: pass
        # else:
        #     #raise IOError, "Unexpected phi range is given"
        #     pass
        gr[i,:]=cp1, cp, cp2, 1./ngrain
        # gr.append([cp1, cp, cp2, 1./ngrain])

    # gr = np.array(gr)
    p1 = gr.transpose()[0]
    p = gr.transpose()[1]
    p2 = gr.transpose()[2]

    if filename!=None and type(filename).__name__=='str':
        FILE = open(filename, 'w')
        FILE.writelines(
            'dummy\ndummy\ndummy\n%s %i\n'%
            ('B', ngrain))
        for i in range(ngrain):
            FILE.writelines(
                "%10.4e %10.4e %10.4e %10.4e\n"%
                (p1[i], p[i], p2[i], 1./ngrain))

    if iplot==True:
        fig = plt.figure(figsize=(7,10))
        ax1 = fig.add_subplot(311); ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)
        ax1.plot(p1, p, '.')
        ax2.plot(p1, p2, '.')
        ax3.plot(p2, p, '.')
        ax1.set_xlabel(r'$\phi_1$',dict(fontsize=28))
        ax1.set_ylabel(r'$\Phi$',dict(fontsize=28))
        ax2.set_xlabel(r'$\phi_1$',dict(fontsize=28))
        ax2.set_ylabel(r'$\phi_2$',dict(fontsize=28))
        ax3.set_xlabel(r'$\phi_2$',dict(fontsize=28))
        ax3.set_ylabel(r'$\Phi$',dict(fontsize=28))
        plt.tight_layout()

    return gr


def rve_ortho(cod, rve):
    """
    Apply the orthonormal sample symmetry to the given discrete
    crystallographic orientation distribution (COD)

    Orthonormal (phi1 = 90)
    Monoclinic (phi1=180)
    None (Triclinic) (phi1=360)
    """
    from euler import euler

    codt = cod.transpose()
    ## information ------------------
    p1max = max(codt[0])  #phi1
    print 'p1max: %4.1f'%p1max
    # phi1 = codt[0]
    # phi2 = codt[1]
    # phi = cot[2]
    ## ------------------------------

    if p1max==90: ssym="Orth"
    elif p1max==180: ssym="Mono"
    elif p1max==360: ssym="Tric"
    else: raise IOError, "Unexpected maximum phi1 anlge"
    print 'symmetry: %s'%ssym

    new_rve = [ ]
    for igr in range(len(rve)):
        ## Phi1, Phi, Phi2 angles and volume fraction
        p1 = rve[igr][0]; p = rve[igr][1]
        p2 = rve[igr][2]; vf = rve[igr][3]

        ## rotation matrix of the current grain
        amat = euler(p1, p, p2, echo=False)
        amat_t = amat.transpose()
        amat_new = []
        if ssym=="Orth": # rolling sample symmetry: mmm
            ## multiplication of the matrix according to the symmetry

            # x-mirror
            oldt = amat_t.copy()
            oldt[1] = oldt[1]*-1
            oldt[2] = oldt[2]*-1
            amat_new.append(oldt.transpose())

            # y-mirror
            oldt = amat_t.copy()
            oldt[0] = oldt[0]*-1
            oldt[2] = oldt[2]*-1
            amat_new.append(oldt.transpose())

            # x and y-mirror
            oldt = amat_t.copy()
            oldt[0] = oldt[0]*-1
            oldt[1] = oldt[1]*-1
            amat_new.append(oldt.transpose())

            nvol = 4
            pass

        elif ssym=="Mono":
            # x-mirror (along TD)
            oldt = amat_t.copy()
            oldt[1] = oldt[1]*-1
            oldt[2] = oldt[2]*-1
            amat_new.append(oldt.transpose())
            nvol = 2

            pass

        elif ssym=="Tric":
            nvol=1
            #no mirror axis
            pass

        ## assigns the newly multiplied A-matrix to the new_rve
        temp = rve[igr].copy(); temp[3] = vf/nvol
        new_rve.append(temp)
        for i in range(len(amat_new)):
            ph1, ph, ph2 = euler(a=amat_new[i],echo=False)
            new_rve.append([ph1,ph,ph2,vf/nvol])
            pass
        pass
    return np.array(new_rve)


class RVE:
    """
    Generates the random iostropic aggregate
    and combines with the discrete COD.

    It assumes that the give COD is the reduced space
    in that the boundary of the Euler space of the COD
    indicates the sample symmetry

    Arguments:
       ngrain : number of grains in the RVE
       odf : Orientation Distribution File (Labotex format)
       cmbfile : output file of the created RVE
       ssym : None (0, 1) sample symmetry
       fmt = 'labo' : labotex Discrete ODF
             'mtex' : mtex discrete ODF
    """
    def __init__(self, ngrain=100, odf=None, cmbfile='temp.cmb',
                 ssym=None, fmt='labo'):
        """
        Arguments
        ---------
        """
        ## globally shared info ----------------------------- ##
        if fmt=='labo':
            self.cod = np.loadtxt(odf,skiprows=1)
            self.codt = self.cod.T
            self.resolution = self.codt[0][1] - self.codt[0][0]
            self.inc = self.resolution
            self.p1max = max(self.codt[0])  #phi1
            self.p2max = max(self.codt[1])  #phi2
            self.p3max = max(self.codt[2])  #phi
        elif fmt=='mtex':
            ## determine the number of head lines:
            with open(odf,'r') as fo:
                ibreak=False
                nhead=0
                while(not(ibreak)):
                    try:
                        map(float,fo.readline().split())
                    except:
                        nhead=nhead+1
                    else:
                        ibreak=True

            # nhead = 4
            F     = np.loadtxt(odf,skiprows=nhead)
            Ft    = F.T
            p1mx  = Ft[0,-1]#max(Ft[0]) # phi1
            pmx   = Ft[1,-1]#max(Ft[1]) # PHI
            p2mx  = Ft[2,-1]#max(Ft[2]) # phi2

            if p1mx==355 and pmx==90 and p2mx==85:
                self.p1max=360.
                self.p2max=90.
                self.p3max=90.
                self.resolution = 5.
                self.inc=self.resolution
            else:
                print 'phi1:',p1mx
                print 'PHI :',pmx
                print 'phi2:',p2mx
                raise IOError, 'unexpected format'+\
                    ' in the given mtex odf file'

            res   = 5
            nphi1 = 360 / res + 1
            nphi  =  90 / res + 1
            nphi2 =  90 / res + 1
            cod_  = np.zeros((nphi1,nphi,nphi2))

            n=0
            for i in xrange(nphi2-1):
                for j in xrange(nphi):
                    for k in xrange(nphi1-1):
                        cod_[k,j,i] = F[n][3]
                        n=n+1

            cod_[-1,:,:] = cod_[0,:,:]
            cod_[:,:,-1] = cod_[:,:,0]

            ## swap 2nd and 3rd axes to comply with labo format
            cod_ = cod_.swapaxes(1,2) ## [phi1,phi2,phi]
            self.cod = np.zeros((nphi1*nphi*nphi2,4))
            n = 0
            for i in xrange(nphi):
                for j in xrange(nphi2):
                    for k in xrange(nphi1):
                        self.cod[n,:3] = k*5.,j*5,i*5
                        self.cod[n,3]  = cod_[k,j,i]
                        n=n+1
            self.codt = self.cod.T
        else:
            raise IOError, 'Unexpected file format requested'


        self.ngrain = ngrain
        ## -------------------------------------------------- ##

        ## Whether the sample symmetry is applied ----------- ##
        #nrot = int(360./self.p1max)
        #print "The max phi1 angle in the COD is %f"%self.p1max
        #ngrain = ngrain / nrot
        ## the isotropic random grain number is set to be 1 of nrot
        ## since the number of grain will be multiplied during
        ## sample symmetry application (see def rve_orth above)
        ## -------------------------------------------------- ##

        print 'max values %5.1f %5.1f %5.1f'%(
            self.p1max, self.p3max, self.p2max)
        self.cmbfile = cmbfile

        ## isotropic sampling, to which the COD is mapped.
        print 'maximum angles along each axis'
        print 'phi1: %6.2f'%max(self.codt[0])
        print 'phi2: %6.2f'%max(self.codt[1])
        print 'phi : %6.2f'%max(self.codt[2])

        print 'minimum angles along each axis'
        print 'phi1: %6.2f'%min(self.codt[0])
        print 'phi2: %6.2f'%min(self.codt[1])
        print 'phi : %6.2f'%min(self.codt[2])

        # if ssym==0: ngrain = ngrain/2
        # elif ssym==1: ngrain = ngrain/4
        # elif ssym==2: ngrain = ngrain/8
        if ssym=='mmm': ngrain = ngrain / 4
        elif ssym==None: ngrain=ngrain
        else: raise IOError, "Unexpected sample symmetry is given"

        print 'ssym:', ssym, 'ngrain=',ngrain

        self.gr = random(
            phi1=max(self.codt[0]), phi2=max(self.codt[1]),
            phi=max(self.codt[2]), ngrain=ngrain, iplot=False)

        ## Calculates the RVE
        self.rve = self.cmb()

        ## Sample symmetry option
        if ssym=='mmm': self.sample_sym_ortho()
        elif ssym==None: pass
        else: raise IOError, "Symm should be either 'mmm' or None"


        ## file output
        self.rve = np.array(self.rve)
        FILE = open(self.cmbfile, 'w')


        import time
        FILE.writelines('%s\n%s\n%s\n %s %i\n'%(
                time.asctime(),'Current texture file was made by cmb.py',
                'contact: youngung.jeong@gmail.com',
                'B', self.ngrain))
        for i in range(len(self.rve)):
            FILE.writelines(
                '%13.7f %13.7f %13.7f   %13.7e \n'%
                (self.rve[i][0], self.rve[i][1], self.rve[i][2],\
                     self.rve[i][3]))
        FILE.close()

    def rot(self,r):
        """
        Given rotation matrix that transform the sample axes to another sample axes,
        apply to the all constituent grains in the aggregate.

        Argument
        --------
        r
        """
        from euler import euler
        for i in xrange(len(self.rve)):
            a1,a2,a3 = self.rve[i][:3]
            rot = euler(a1,a2,a3,echo=False) ## CA<-SA(old)
            rt = rot.T        ## SA<-CA
            rt = np.dot(r,rt) ## SA`<-SA<-CA
            rot = rt.T         ## CA<-SA`
            a1,a2,a3 = euler(a=rot,echo=False)
            self.rve[i][:3]=a1,a2,a3

    def write(self,grs,fn):
        import time
        ngrain = len(grs)
        with open(fn, 'w') as FILE:
            FILE.writelines('%s\n%s\n%s\n %s %i\n'%(
                time.asctime(),'Current texture file was made by cmb.py',
                'contact: youngung.jeong@gmail.com',
                'B', ngrain))
            for i in range(len(grs)):
                FILE.writelines(
                    '%13.7f %13.7f %13.7f   %13.7e \n'%
                    (grs[i][0], grs[i][1], grs[i][2],\
                     grs[i][3]))

    def plot(self,csym='cubic',**kwargs):
        """
        Plot pole figure
        Arguments
        ---------
        csym : crystal symmetry
        **kwargs: kw arguments to be passed to pf_new
        """
        import upf,time
        mypf = upf.polefigure(grains =  self.rve, csym=csym)
        fig=mypf.pf_new(**kwargs)
        plt.show()
        return fig

    def reflect(self,th):
        import math
        pi = math.acos(0.) * 2
        c2t = math.cos(2. * th * pi / 180.)
        s2t = math.sin(2. * th * pi / 180.)
        r = np.zeros((3,3))

        # reflection matrix: (reflection by xz plane with th=0)
        #                    (reflection by yz plane with th=2/pi)

        r[0,0] = c2t
        r[0,1] = s2t
        r[1,0] = s2t
        r[1,1] = c2t * (-1.)
        r[2,2] = 1.0
        return r

    def sample_sym_ortho(self):
        """
        Application of sample symmetry

        'mmm' operation (orthorhombic sample symmetry operation
        """
        from euler import euler

        m0 = [[  1,  0,  0],  [  0,  1,  0],  [  0,  0,  1]]
        m1 = [[  1,  0,  0],  [  0, -1,  0],  [  0,  0, -1]]
        m2 = [[ -1,  0,  0],  [  0,  1,  0],  [  0,  0, -1]]
        m3 = [[ -1,  0,  0],  [  0, -1,  0],  [  0,  0,  1]]
        mmm = [m1,m2,m3]

        # # Rmir0 = self.reflect(th=0.) # RD-ND
        # # Rmir1 = self.reflect(th=90.) # TD-ND
        # # Rmir2 = np.zeros((3,3))
        # # Rmir2[0,0] = 1; Rmir2[1,1] = 1; Rmir2[2,2] = -1

        # print Rmir0
        # print Rmir1
        # print Rmir2
        # # raw_input()
        # # Rmir2  rotate Rmir0 by x axis 90 degree

        # if iopt==0:
        #     Rm = [Rmir0]
        # elif iopt==1:
        #     Rm = [Rmir0,Rmir1]
        # elif iopt==2:
        #     Rm = [Rmir0,Rmir1,Rmir2]
        # else: raise IOError, 'Wrong option'

        gr = [] # new grains
        for i in range(len(self.rve)):
            tempgr=[]
            # itself
            tempgr.append([self.rve[i][0],self.rve[i][1],self.rve[i][2],
                       self.rve[i][3]])
            cRs = euler(self.rve[i][0],self.rve[i][1],self.rve[i][2],
                        echo=False) # R (ca<-sa)
            for j in range(len(mmm)):
                R0 = np.dot(cRs, mmm[j]) # (ca<-sa<-sa)
                phi1, phi, phi2 = euler(a=R0, echo=False)
                tempgr.append([phi1,phi,phi2,self.rve[i][3]])

            for j in range(len(tempgr)):
                gr.append(tempgr[j])


            # R0 = np.dot(cRs, Rm[0]) # RD-TD
            # phi1, phi, phi2 = euler(a=R0, echo=False)
            # tempgr.append([phi1,phi,phi2,self.rve[i][3]])

            # #
            # for i in range(len(Rm)):
            #     np.dot(cRs, Rm[i]
            # if iopt>0:
            #     n = len(tempgr)
            #     for j in range(n):
            #         cRs = euler(tempgr[j][0],tempgr[j][1],tempgr[j][2],
            #                     echo=False)
            #         R0 = np.dot(cRs, Rm[1])
            #         phi1,phi,phi2 = euler(a=R0, echo=False)
            #         tempgr.append([phi1,phi,phi2,self.rve[i][3]])

            # if iopt>1:
            #     n = len(tempgr)
            #     for j in range(n):
            #         cRs = euler(tempgr[j][0], tempgr[j][1], tempgr[j][2],
            #                     echo=False)
            #         R0 = np.dot(cRs, Rm[2])
            #         phi1,phi,phi2 = euler(a=R0, echo=False)
            #         tempgr.append([phi1,phi,phi2,self.rve[i][3]])

            # for j in range(len(tempgr)):
            #     gr.append(tempgr[j])


        self.rve = np.array(gr)

        wgt = self.rve.T[-1]
        twgt = wgt.sum()
        self.rve[:,3] = self.rve[:,3]/twgt


    def index_odf(self, i, j, k):
        """
        Returns the value of the cod
        corresponding index i,j,k for
        phi1, phi, and phi2, respectively.
        """
        isize= int(self.p1max/self.resolution) + 1
        jsize= int(self.p2max/self.resolution) + 1
        ksize= int(self.p3max/self.resolution) + 1

        a = k * isize * jsize
        b = j * isize
        c = i
        ind = int(a) + int(b) + int(c)
        return self.cod[ind][3]

    def search(self,  phi1, phi2, phi):
        """
        Provided the phi1, phi2, and phi,
        find the Euler cube which contains the point.
        Then, estimiates the volume fraction...
        """
        if phi1>self.p1max or phi2>self.p2max or phi>self.p3max:
            print 'phi1, phi2, phi: %6.3f  %6.3f  %6.3f'%(
                phi1, phi2, phi)
            raise IOError, "The given angle is not available"

        i = int(phi1/self.resolution)
        j = int(phi2/self.resolution)
        k = int(phi/self.resolution)
        return self.index_odf(i,j,k)

    def interpolate(self, phi1, phi2, phi):
        """
        Interpolates the volume fraction
        """
        r = [ [], [], [] ]
        r[0].append(phi1 - phi1%self.inc)             #phi
        r[0].append(phi1 - phi1%self.inc  +self.inc)
        r[1].append(phi2 - phi2%self.inc)             #phi2
        r[1].append(phi2 - phi2%self.inc  +self.inc)
        r[2].append(phi  - phi %self.inc)             #phi
        r[2].append(phi  - phi %self.inc  +self.inc)
        value = 0
        x, y, z = 0, 0, 0
        for j in range(2):              #phi1
            for k in range(2):          #phi2
                for l in range(2):      #phi
                    x = self.inc - abs(phi1 - r[0][j])
                    y = self.inc - abs(phi2 - r[1][k])
                    z = self.inc - abs(phi  - r[2][l])
                    value = value + self.search(
                        phi1=r[0][j], phi2=r[1][k],
                        phi=r[2][l]) * x * y * z / (self.inc**3)
                    pass
                pass
            pass
        return value

    def cmb(self):
        """
        Interpolates the individual grains' volume fraction
        from the discrete crystallographic orientation distribution file
        and returns the representative volume element.
        """
        RVE = [ ]
        tiny = 10**-6
        for i in range(len(self.gr)):
            ## Trick to bypass the case in which the random grain is at the edge
            ## of the given Euler space of the COD.
            if self.gr[i][0]==self.p1max: self.gr[i][0]=self.gr[i][0]-tiny#phi1
            if self.gr[i][1]==self.p3max: self.gr[i][1]=self.gr[i][1]-tiny#phi
            if self.gr[i][2]==self.p2max: self.gr[i][2]=self.gr[i][2]-tiny#phi2
            ## -------------------------------------------------------------- ##

            ## interpolation of the COD value for each of the grain
            value = self.interpolate(phi1=self.gr[i][0],
                                     phi2=self.gr[i][2],
                                     phi=self.gr[i][1])
            RVE.append([self.gr[i][0], self.gr[i][1],
                        self.gr[i][2], value])
            ## ----------------------------------------------------
            pass


        return RVE
    pass # end of the class RVE

def main(odf, ngrain, outputfile, iplot=True, irandom=False):
    """
    Arguments
    =========

      odf: od filename
      ngrain : number of grains in the RVE
      outputfile : combined file
      iplot =True : flag for plotting
      irandom=True
    """
    import upf,time
    import matplotlib.pyplot as plt
    temp = RVE(ngrain=ngrain, odf=odf, cmbfile=outputfile)

    if iplot==True:
        mypf = upf.polefigure(grains =  temp.rve, csym='cubic')
        mypf.pf(pole=[[1,0,0],[1,1,0],[1,1,1]], mode='contourf', ifig=2)
        plt.show()
        pass

    if irandom==True:
        filename='iso.cmb'
        FILE = open(filename, 'w')
        FILE.writelines('%s\n%s\n%s\n %s %i\n'%(
                time.asctime(),'Current texture file was made by cmb.py',
                'contact: youngung.jeong@gmail.com',
                'B', int(ngrain)))
        for i in range(len(temp.gr)):
            FILE.writelines('%11.7f %11.7f %11.7f  %10.4e \n'%
                            (temp.gr[i][0], temp.gr[i][1],
                             temp.gr[i][2], temp.gr[i][3]))
            pass
        FILE.close()
        #np.savetxt(filename, temp.gr)
        print 'The random isotropic aggregate is written to %s'%filename
        pass
    pass

if __name__=='__main__':
    import getopt, sys
    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                    'i:n:o:sr')
        pass
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
        pass

    ##  default options  ##
    iplot = False
    ngrain = 2000
    outputfile ='temp.cmb'
    irandom = False
    ## ----------------- ##

    for o, a in opts:
        if o in ('-i'): inputfile = a
        elif o in ('-n'): ngrain = int(a)
        elif o in ('-o'): outputfile = a
        elif o in ('-s'): iplot=True
        elif o in ('-r'): irandom=True
        else: assert False, 'Unhandled option'
        pass

    main(odf=inputfile, ngrain=ngrain,
         outputfile=outputfile, iplot=iplot,
         irandom=irandom
         )

    pass


def random_gen(ngrain=100,mmm=False,phi1=360,phi2=360,phi=360):
    """
    Generate isotropic random polycrystal texture files in case of
    cubic crystal symmetry

    Arguments
    =========
    ngrain = 100
    mmm = False
    phi1
    phi2
    phi
    """
    import upf
    # if mmm:
    #     phi1 = 90.
    #     phi2 = 90.
    #     phi  = 90.
    #     ngrain = int(ngrain/4.)
    # if mmm==False:
    #     phi1 = 360.
    #     phi2 = 360.
    #     phi  = 360.
    #     ngrain = ngrain / 1
    # apply mmm sample symmetry...

    filename='iso_%s.cmb'%str(ngrain).zfill(5)
    random(phi1,phi2,phi,ngrain,iplot=False,filename=filename)
    print '%s has been created'%filename
    upf.cub(filename,ifig=10)
