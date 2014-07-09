try: np
except: import numpy as np
from numpy import linalg as LA
import math
import euler
eul = euler.euler


def __inv__(a):
    """
    Returns the inverse matrix of A
    
    *Note:
       Matrix a should be np.array
    """
    if type(a).__name__ =='ndarray': pass
    else:
        print 'Error: Unexpected Matrix a type'
        print 'It should be numpy.ndarry!'
        raise IOError
    ainv = LA.inv(a)
    if np.allclose(np.dot(ainv,a),np.eye(3)): return ainv
    else:
        print'Error: problem in getting an inverse matrix'
        raise IOError

def miso(a, b):
    """
    misorientation between given A and B rotation matrices
    """
    types = [type(a).__name__, type(b).__name__]
    if any(types[i]!='ndarray' for i in range(2)):
        print 'Error: Unexpected Matrix-a type'
        print 'It should be numpy.ndarry!'
        raise IOError
    dg = np.dot(a,__inv__(b))
    cth = (dg[0,0] + dg[1,1] + dg[2,2] - 1)/2.
    r1 = dg[1][2] - dg[2][1]
    r2 = dg[2][0] - dg[0][2]
    r3 = dg[0][1] - dg[1][0]
    return math.acos(cth)*180./math.pi  #Return as an angle

def del_A(p1,p2):
    """
    Misorientation between given p1 and p2 euler angle
    """
    #a= eul(ph=phi[0], th=phi[1], tm=phi[2],echo=False)
    a = eul(ph=p1[0],  th=p1[1],  tm=p1[2], echo=False)
    b = eul(ph=p2[0],  th=p2[1],  tm=p2[2], echo=False)
    a = np.array(a)
    b = np.array(b)
    m = miso(a,b)
    return m

def fibr(filename='texture/500.tex'):
    """
    Provided the file object, getting a certain euler space sliced.
    To do so, interpolation of a grid of that slice is performed.
    The factor of contribution of a each point is decayed with
    increasing the Eulerian distance between the grid point.
    The number of point in the Euler space is prescribed.
    """
    FILE = open(filename,'r')
    lines = FILE.readlines()
    nblock = 0  # total number of texture file blocks
    iblock = 0 # id for each block
    iline = 0   # id for each line
    grains = []  # raw data

    while True:  #loop over of blocks

        try: len(lines[iline])
        except IndexError: print 'EOF';  break
        else: grains.append([])

        if len(lines[iline])<2: print 'EOF'; break
        
        for i in range(4):
            cwl = lines[iline]
            iline = iline + 1
            pass
        ngr = int(cwl.split()[1])

        for i in range(ngr):  #loop over grains
            cwl = lines[iline]
            grains[iblock].append([])
            grains[iblock][i] = (map(float, cwl.split()))
            iline = iline + 1
            pass
        iblock = iblock + 1
                
    nblock = iblock
    return grains
    #return lines


def interpolation(gr, origin):
    """
    Provided:
        a grain's info
        a origin's euler angle

    returns vol_fration, misorientation angle
    """
    phi = gr[0:3]
    vf = gr[3]
    o = origin[0:3]

    """ Calculate the misorientation """
    a = eul(ph=phi[0], th=phi[1], tm=phi[2],echo=False)
    b = eul(ph=o[0], th=phi[1], ntm=phi[2],echo=False)
    mis = miso(a,b)
    if mis<0:
        print "Warning: negative misorientation happens"
        raw_input('Enter if you want to keep going on')

    return vf, mis


def slice(filename='texture/500.tex'):
    """
    to be filled
    """
    FILE = open(filename, 'r')
    pass


def grid(phi1 = None, phi = None, phi2 = None, pmx = [180,90,90],
         inc = 5 ):  # Provides a certain grid points
                     # around which interpolation is performed.
    """
    Provided one of the three Euler angles,
    grid of that plane is made and returned.
    """
    phi_ = [phi1, phi, phi2]
    if all(phi_[i]==None for i in range(3)):
        print 'At least one of the Euler angles must be given'
    else:pass

    inon = 0
    for i in range(3):
        if phi_[i] == None:
            inon = inon + 1
    
    if inon!=2:
        print 'Only two of angles should be None'
        print 'You must give one of three angles'
        raise IOError

    #xyz = np.resize((),(pmx[0]/inc+1,pmx[1]/inc+1,pmx[2]/inc+1))
    xyz = []

    inctyp = []

    for i in range(3):
        inctyp.append(type(pmx[i]/inc))
    
    if any(inctyp[i].__name__!='int' for i in range(len(inctyp))):
        print 'increment does not produce a integer grid'
        print 'consider to change this'

    if phi1!=None:
        " x: phi2, y: phi"
        for i in range(pmx[1]/inc + 1): 
            for j in range(pmx[2]/inc +1): 
                xyz.append([phi1, i*inc, j*inc]) #phi1, phi, phi2
        return xyz

    if phi2!=None:
        " x: phi1, y: phi"
        for i in range(pmx[0]/inc + 1):
            for j in range(pmx[1]/inc + 1):
                xyz.append([i*inc, j*inc, phi2])
        return xyz

    if phi!=None:
        " x: phi1, y: phi2"
        for i in range(pmx[0]/inc + 1):
            for j in range(pmx[2]/inc + 1):
                xyz.append([i*inc, phi, j*inc])
        return xyz




def dg(phi1,phi,phi2):
    """
    dg = 1/ 8pi^2 * sin(phi) dp1 dp dp2
    """
    pass




def neighbours(mis, grains, xyz):
    """
    Provided the coordinate of a point,
    returns certain grains which are within
    the misorientation width
    
    mis = misorientation 
    grains = grains
    xyz = origion
    """
    rstg = [] #eul(ph=phi[0], th=phi[1], tm=phi[2])
    a_xyz = eul(ph=xyz[0],th=xyz[1],tm=xyz[2],echo=False)
    for i in range(len(grains)):
        a_gr = eul(ph=grains[i][0],
                   th=grains[i][1],
                   tm=grains[i][2],
                   echo=False)
        mo = abs(miso(a_xyz,a_gr))
        if mo <= mis:
            temp = grains[i][0:4]
            temp.append(mo)
            rstg.append(temp)
    return rstg



def __a2i__(ang=5.,inc = 5,pm=90.):
    """
    Convert angle to index
    """
    ang = int(ang)
    inc = int(inc)
    if ang>pm: 
        print 'Your angle is beyond the maximum'
        raise IOError
#    if ang==pm:
#        ang = 0
    return int(int(ang)/int(inc))



def as2i(p=[5,10,30], inc= 5):
    a=[]
    a.append(__a2i__(ang=p[0], inc=inc))
    a.append(__a2i__(ang=p[1], inc=inc))
    a.append(__a2i__(ang=p[2], inc=inc))
    return a


def corner8(ind):
    c = np.resize((),(8,3))
    c[0] = ind
    c[1] = [ind[0]+1, ind[1]  , ind[2]]
    c[2] = [ind[0]  , ind[1]+1, ind[2]]
    c[3] = [ind[0]+1, ind[1]+1, ind[2]]
    c[4] = [ind[0]  , ind[1]  , ind[2]+1]
    c[5] = [ind[0]+1, ind[1]  , ind[2]+1]
    c[6] = [ind[0]  , ind[1]+1, ind[2]+1]
    c[7] = [ind[0]+1, ind[1]+1, ind[2]+1]
    return np.array(c,int)


class distr:
    def __init__(self,grains, inc = 5):
        """
        grains must be a single block of texture file
        """
        ng = len(grains)
        pmx = [360., 90., 90.]
        
        self.ODF = np.resize((0.),(__a2i__(pmx[0],inc=inc,pm=pmx[0])+1,
                              __a2i__(pmx[1],inc=inc,pm=pmx[1])+1,
                              __a2i__(pmx[2],inc=inc,pm=pmx[2])+1))
        
        print len(self.ODF)
        print len(self.ODF[0])
        print len(self.ODF[0][0])
        #raw_input()

        phi1 = [] ; phi =[]; phi2 = []; odf =[]
        for i in range(len(grains)):
            phi1.append(grains[i][0])
            phi.append (grains[i][1])
            phi2.append(grains[i][2])

        for i in range(ng):
            inte = grains[i][3]
            cph = [grains[i][0],grains[i][1],grains[i][2]]
            cph_ind =[ __a2i__(grains[i][0], inc=inc, pm=pmx[0]+inc),
                       __a2i__(grains[i][1], inc=inc, pm=pmx[1]+inc),
                       __a2i__(grains[i][2], inc=inc, pm=pmx[2]+inc)]

            corners = corner8(cph_ind)
            
            if grains[i][0]==pmx[0]: grains[i][0]=grains[i][0]-0.00001
            if grains[i][1]==pmx[1]: grains[i][1]=grains[i][1]-0.00001
            if grains[i][2]==pmx[2]: grains[i][2]=grains[i][2]-0.00001
            self.interpolate(grain=grains[i],corners =corners)



            
    def write(self, phi1=None, phi=None, phi2=None, filename='od_section.out', inc=5):
        """
        Provided one of the Euler angles,
        section of OD is given and written into the file 
        whose name is given in the current working directory
        """
        FILE = open(filename, 'w')
        FILE.writelines('** header \n')
        #FILE.writelines('*
        p = [phi1, phi, phi2]
        nnone = 0
        if all(p[i]==None for i in range(len(p))):
            print 'Wrong! one of the angles must be given'
            raise IOError
        else:
            for i in range(len(p)):
                if p[i]!=None: nnone=nnone+1
            if nnone==1:
                if phi1!=None:
                    FILE.writelines(' ** OD section of phi1 = %4i \n'%(phi1))
                    """
                    x=phi2, y=phi
                    """
                    ind = phi1/inc
                    for i in range(len(self.ODF[ind])):
                        for j in range(len(self.ODF[ind][i])):
                            FILE.writelines('%8.3f '%(self.ODF[ind][i][j]))
                        FILE.writelines('\n')
                if phi !=None:
                    FILE.writelines(' ** OD section of phi  = %4i \n'%(phi))
                    """
                    x=phi1, y= phi2
                    """
                    ind = phi/inc
                    for i in range(len(self.ODF)):
                        for j in range(len(self.ODF[i][ind])):
                            FILE.writelines('%8.3f '%(self.ODF[i][ind][j]))
                        FILE.writelines('\n')
                if phi2!=None:
                    FILE.writelines(' ** OD section of phi2 = %4i \n'%(phi2))                    
                    """
                    x=phi1, y=phi
                    """
                    ind = phi2/inc
                    for i in range(len(self.ODF)):
                        for j in range(len(self.ODF[i][j])):
                            FILE.writelines('%8.3f '%(self.ODF[i][j][ind]))
            else: print'Only one of the euler angles must be given'; raise IOError

        # temp = ODF[0][0][0] + ODF[-1][0][0]
        # ODF[0][0][0] = temp
        # ODF[-1][0][0] = temp
        
        # temp = ODF[0][1][0] + ODF[-1][1][0]
        # ODF[0][1][0] = temp
        # ODF[-1][0][0] = temp

        # temp = ODF[0][..][loop] + ..

        # temp = ODF[0][-1][0] +ODF[0][0][0]
        # ODF[0][-1][0] = temp
        # ODF[0][0][0] = temp 



    def interpolate(self, grain, corners):
        inten = grain[3]
        for i in range(len(corners)):
            c = corners[i]
            #print c[0],c[1],c[2]
            try:self.ODF[c[0]][c[1]][c[2]]
            except:
                print ''
                print 'Error Happened when c is below'
                print c[0],c[1],c[2]
                print ''
                print len(self.ODF)
                print len(self.ODF[0])
                print len(self.ODF[0][0])
                print grain
                print ' '
            else:
                temp = grain[3]
                self.ODF[c[0]][c[1]][c[2]] = self.ODF[c[0]][c[1]][c[2]] + temp

    def dist(self,p1,p2):
        x1=p1[0]
        y1=p1[1]
        z1=p1[2]

        x2=p2[0]
        y2=p2[1]
        z2=p2[2]

        return math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

class sect:
    """
    Provide the sectional grid info of ODF calculated from VPSC
    texture type file (having four lines of header on the top which
    is, then, followed by grains, with each consisted of three
    Euler angles and corresponding volume fraction

    Provided the texture and the list of arguments,
    the grided ODF section is written down to a file.
    """
    def __init__(self, filename='texture/500.tex', 
                 iblock=None, phi1=None, phi=None, phi2=None, 
                 misrad=15., grid_ang=5, pmx = [180,90,90], fout='od_section.out'):
        """
        Arguments:
            phi1, phi, phi2
            iblock : id of block in case you have multi blocks in 
                     the TEXTURE file
            misrad : radius of misorientation
                     The radius of the circle in which the interpolation is performed
                     in the Euler space
            grid_ang = 5,
            pmx = [180,90,90]
            fout = 'od_section.out'

        Note: empty
        """
        FILE = open(filename,'r')
        self.fn = filename
        self.pmx = pmx
        self.ang = [phi1, phi, phi2]
        self.grid_ang = grid_ang
        mrad = misrad * math.pi / 180.

        self.grains = fibr(FILE.name)
        if iblock == None: iblock = 0
        self.grains = self.grains[iblock]
        self.xyz = grid(phi1=phi1, phi=phi, phi2=phi2, inc=grid_ang, pmx=self.pmx)

        """ e.g, xyz= [   [0., phi, 0.,]  ... [180., phi., 90.] ]  """

        intensity = []
        for i in range(len(self.xyz)):   #xyz: Cooridnates for the grid of the section
            n = neighbours(mis=misrad, grains=self.grains, xyz=self.xyz[i])
            # n[j] = [phi1, phi, phi2, ODF, misorientation]
            inten = 0.
            #print n
            #raw_input()
            for j in range(len(n)):
                inten = inten + n[j][3]/n[j][4]
            intensity.append(inten)
            #intensity.append(inten  /(mrad**3))
        isum = sum(intensity)
        for i in range(len(self.xyz)):
            self.xyz[i].append(intensity[i]/isum)
            pass

        self.print_xyz(filename=fout)


    def print_xyz(self,filename='od_section.out'):
        """
        Calcaulted xyz into a Grid file
        """
        FILE = open(filename,'w')
        grid = self.grid_ang
        xyz  = self.xyz
        phi1 = self.ang[0]
        phi  = self.ang[1]
        phi2 = self.ang[2]
        p1   = self.pmx[0]
        P    = self.pmx[1]
        p2   = self.pmx[2]

        FILE.writelines('** filename :'+self.fn+'\n')

        if phi1!=None:
            FILE.writelines('** Phi1 = %5.1f section ** \n'%(phi1))
            FILE.writelines(' x- axis : phi2, y-axis : phi \n')
            
            igrid = 0 
            for i in range(p2/grid+1):
                for j in range(P/grid +1):
                    FILE.writelines('   %11.3e'%(xyz[igrid][3]))
                    igrid = igrid + 1
                FILE.writelines('\n')
            """
            axis:
                x:phi2, y:phi
            """
            pass

        elif phi!=None:
            FILE.writelines('** Phi = %5.1f section ** \n'%(phi))
            FILE.writelines(' x- axis : phi1, y-axis : phi2 \n')

            igrid = 0
            for i in range(p1/grid + 1):
                for j in range(p2/grid +1):
                    FILE.writelines('   %11.3e'%(xyz[igrid][3]))
                    igrid = igrid + 1
                FILE.writelines('\n')
            """
            axis
                x:phi1, y:phi2
            """
            pass

        elif phi2!=None:
            FILE.writelines('** Phi2 = %5.1f section ** \n'%(phi2))
            FILE.writelines(' x- axis : phi1, y-axis : phi \n')

            igrid = 0
            for i in range(p1/grid + 1):
                for j in range(P/grid + 1):
                    FILE.writelines('   %11.3e'%(xyz[igrid][3]))
                    igrid = igrid + 1
                FILE.writelines('\n')
            """
            axis:
                x:phi1, y=phi
            """
            pass
        FILE.close()
