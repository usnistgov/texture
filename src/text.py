"""
Gaussian texture components.

With given components (final aim, will start from given Euler angle set)


Materials Mechanics Lab,
Graduate Institute of Ferrous Technology
POSTECH, Korea

Youngung Jeong

youngung.jeong@gmail.com
"""
import random,math,os,upf
from euler import euler as eul
gauss = random.gauss
sin = math.sin
cos = math.cos
atan = math.atan
atan2 = math.atan2
asin = math.asin
acos = math.acos
pi = math.pi
sqrt = math.sqrt
# coef= math.atan(1.0)/45. #*math.pi/180.

def pn(n,b):
    from cs import uniqpnset# f2py-wrapped binary module.
    """
    normal and direction
    """
    import numpy as np
    uniqset, nn = uniqpnset(n=n, b=b, isym=1, ipr=False)
    uniqset = uniqset.swapaxes(0,2)
    rst = []
    for i in range(nn):
        pn = uniqset[i][0]
        di = uniqset[i][1]
        rst.append([pn,di])

    return np.array(rst)

def orthogonal():
    """
    Returns orthogonal symmetry operators

    (100)_sa mirroring  - RD
    (010)_sa mirroring  - TD


    RD = [[ 1.,  0.,  0.],
          [ 0., -1.,  0.],
          [ 0.,  0.,  1.]]
    TD = [[-1.,  0.,  0.],
          [ 0.,  1.,  0.],
          [ 0.,  0.,  1.]]

    v`[sa] = RD[sa,sa] * v[sa]
    """
    from cs import reflect, vector_ang

    import numpy as np
    rxz = reflect(th = 0) # rxz xz-plane reflection, v`<- [rxz]v
    r = vector_ang(u=[0,0,1], th=90)

    # matrix that rotates about z-axis by 90degree
    ryz = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            ryz[i,j] = 0.
            for k in range(3):
                for l in range(3):
                    ryz[i,j] = ryz[i,j] + r[k,i] * rxz[k,l] * r[l,j]

    RD = rxz
    TD = ryz
    return RD, TD

def c2e(component='cube', ngrain=100, w0=10, ifig=1, mode='contour'):
    """
    Provided (hkl)// ND and [uvw]//RD,

    arguments
    =========
    component = 'cube'
    ngrain    = 100
    w0        = 10
    ifig      = 1
    mode      = 'contour'
    """
    import numpy as np
    from euler import euler
    import upf
    comp = component.lower()
    print comp
    if comp=='cube':
        hkl=[1,0,0]
        uvw=[0,1,0]
    elif comp=='copper':
        hkl=[1,1,2]
        uvw=[-1,-1,1]
    elif comp=='goss':
        hkl=[1,1,0]
        uvw=[0,0,1]
    elif comp=='rbrass':
        hkl=[-1,1,0]
        uvw=[1,1,1]

    print 'hkl',hkl
    print 'uvw',uvw,'\n'

    # mat0: trasformation matrix that converts
    # sample axes to crystal ones [ca<-sa]
    mat0 = miller2mat(hkl=map(round, hkl), uvw=map(round, uvw))

    print 'mat0'
    for i in range(len(mat0)):
        for j in range(len(mat0[i])):
            print '%5.1f '%mat0[i][j],
        print ''

    # Mirror reflection by RD and TD in the sample coordinate system.
    # mat0[ca<-sa]
    RD, TD = orthogonal()
    mat1 = np.dot(mat0, RD)
    mat2 = np.dot(mat0, TD)
    mat3 = np.dot(np.dot(mat0, RD), TD)
    print 'mat1'
    for i in range(len(mat1)):
        for j in range(len(mat1[i])):
            print '%5.1f '%mat1[i][j],
        print ''
    print 'mat2'
    for i in range(len(mat2)):
        for j in range(len(mat2[i])):
            print '%5.1f '%mat2[i][j],
        print ''
    print 'mat3'
    for i in range(len(mat3)):
        for j in range(len(mat3[i])):
            print '%5.1f '%mat3[i][j],
        print ''


    ang0 = euler(a=mat0, echo=False)
    ang1 = euler(a=mat1, echo=False)
    ang2 = euler(a=mat2, echo=False)
    ang3 = euler(a=mat3, echo=False)

    angs = [ang0, ang1, ang2, ang3]

    for i in range(len(angs)):
        for j in iter(angs[i]):
            print '%5.1f '%j,
        print ''

    xtal = []
    for i in range(4):
        phi1, p, phi2 = angs[i]
        myx = textur(p1=phi1,p=p, p2=phi2, w0=w0, ngrain = ngrain/4,
                     filename='dum', iplot=False)

        for j in range(len(myx.gr)):
            xtal.append(myx.gr[j])

    mypf = upf.polefigure(grains=xtal,csym='cubic')
    mypf.pf(pole=[[1,0,0],[1,1,0],[1,1,1]], mode=mode, ifig=ifig)

    ## Write
    gr2f(gr=xtal, header='module: %s,   w0: %i'%('ce2', w0),
         filename='%i%s%i_%s.cmb'%(len(xtal), component, w0, 'ce2'))
    return xtal

def gr2f(gr=None, header=None, filename=None):
    """ Writes grains in the texture file format"""
    import numpy as np
    gr = np.array(gr)
    f = open(filename, 'w')
    if header==None: f.write('dum\ndum\ndum\nB %i\n'%len(gr))
    else: f.write('%s\ndum\ndum\nB %i\n'%(header, len(gr)))
        
    for i in range(len(gr)):
        f.write('%9.2f %9.2f %9.2f  %s\n'%(
                gr[i][0], gr[i][1], gr[i][2], str(1./len(gr))))
    f.close()

def comp2euler(component='cube'):
    """ Provided (hkl)// ND and [uvw]//RD, """
    import numpy as np
    comp = component.lower()
    print comp
    if comp=='cube':
        hkl=[1,0,0]
        uvw=[0,1,0]
    elif comp=='copper':
        hkl=[1,1,2]
        uvw=[-1,-1,1]
    elif comp=='goss':
        hkl=[1,1,0]
        uvw=[0,0,1]
    elif comp=='rbrass':
        hkl=[-1,1,0]
        uvw=[1,1,1]

    print 'hkl',hkl
    print 'uvw',uvw,'\n'

    rst = pn(n=hkl, b=uvw)

    eul = []
    for i in range(len(rst)):
        print 'hkl: [%3.1f,%3.1f,%3.1f]'%(
            rst[i][0][0],rst[i][0][1],rst[i][0][2])
        print 'uvw: [%3.1f,%3.1f,%3.1f]'%(
            rst[i][1][0],rst[i][1][1],rst[i][1][2])
        #hkl = map(int,rst[i][0])
        #uvw = map(int,rst[i][1])
        p1, phi, p2 = miller2euler(hkl=map(round,rst[i][0]),
                                   uvw=map(round,rst[i][1]))
        print 'p1, phi2, p2: %3.1f %3.1f %3.1f'%(p1, phi, p2)
        eul.append([p1, phi, p2])

    return np.array(eul)


def main(component='cube', ngrain=100, w0=10, ifig=1):
    """
    arguments
    =========
    component = 'cube', 'copper', 'rbrass', 'goss'
    ngrain    = 100
    w0        = 10
    ifig      = 1
    """
    import upf
    angs = comp2euler(component)
    xtal = []
    for i in range(len(angs)):
        myx = textur(p1=angs[i][0], p=angs[i][1], p2=angs[i][2], w0=w0,
               ngrain = ngrain/len(angs), filename='dum',
               iplot=False)
        myx.gr

        for j in range(len(myx.gr)):
            xtal.append(myx.gr[j])

    mypf = upf.polefigure(grains = xtal, csym='cubic')
    mypf.pf(pole=[[1,0,0],[1,1,0],[1,1,1]], ifig=ifig, mode='contourf')

    f = open('%s%s.cmb'%(str(len(xtal)),component), 'w')
    f.write('dum\ndum\ndum\nB %i\n'%len(xtal))
    for i in range(len(xtal)):
        f.write('%9.2f %9.2f %9.2f  %s\n'%(
                xtal[i][0], xtal[i][1], xtal[i][2], str(1./len(xtal))))

    return xtal

class textur:
    """
    Design an aggregate for the given distribution parameters
    Arguments :

    Euler angles
    p1         : phi1
    p2         : phi2
    p          : PHI
    w0         : width...
    ngrain     : Number of to be created grains
    filename = 'temp.txt'
    dist = 'g','n','e','l'

    iplot=True
    mode = 'contour'
    pole=[[1,0,0]]
    ifig = 1

    Cube: {100}<100>
    Rotated Brass: {111}<110>
    copper: {112}<-1-11>
    Goss: {110}<100>
    """
    def __init__(self, p1, p2, p, w0, ngrain,
                 filename='temp.txt', dist='g',
                 iplot=True, mode='contour',
                 pole=[[1,0,0],[1,1,0],[1,1,1]], ifig=1):
        f = open(filename, 'w')
        f.writelines('Designed texture '\
                         'using probability distributions \n')
        f.writelines('Given Euler angles are as below \n')
        f.writelines('ph1,   phi2,   PHI = ' + str(p1) + str(p2)+ str(p) )

        f.write('\nB   '+ str(ngrain))
        self.gr = []
        for i in range(ngrain):
            txt = text(p1=p1, p2=p2, p=p, w0=w0, dist=dist)
            angle = txt.angle
            f.write('\n %9.2f %9.2f %9.2f %9.3f'%(
                    angle[0],angle[1],angle[2], 0.1))
            self.gr.append([angle[0],angle[1],angle[2], 0.1])

        f.close()
        if iplot==True:
            mypf = upf.polefigure(filename=filename, csym='cubic')
            mypf.pf(pole=pole, mode=mode, ifig=ifig)


class text:
    """
    Gaussian texture for the given angle
    Arguments:

    p1
    p2
    p

    w0
    dist = 'g','n','e','l'
    """
    def __init__(self, p1=45., p2=45., p=45., w0=15., dist='g'):
        """
        Arguments :
        p1
        p2
        p
        w0
        """
        import numpy as np

        ## generates a random rotation axis, about which
        ## the given grain rotates. Axis is represented by 
        ## two angles: delta and phi. 
        # delta, phi = self.rot_axis()
        # print "Rotation axis's delta, phi ",\
        #     delta*180./np.pi, phi*180./np.pi

        # if dist=='g':
        #     B = self.transform_matrix(delta=delta,
        #                               phi=phi,
        #                               w=self.gaussian(w0))
        # if dist=='e':
        #     B = self.transform_matrix(delta=delta,
        #                               phi=phi,
        #                               w=self.expo(w0))
        # if dist=='l':
        #     B = self.transform_matrix(delta=delta,
        #                               phi=phi,
        #                               w=self.lognorm(w0))
        # if dist=='n':
        #     B = self.transform_matrix(delta = delta,
        #                             phi = phi,
        #                             w=self.normal(w0))

            
        # A transforms from SA to CA [ca<-sa]
        A = eul(ph = p1, th = p, tm = p2, echo=False)

        if dist=='g':
            w = self.gaussian(w0)
        elif dist=='e':
            w = self.expo(w0)
        elif dist=='l':
            w = self.lognorm(w0)
        elif dist=='n':
            w = self.normal(w0)

        C = self.rot_vectang(th = w, r = A)

        p1, p, p2 = eul(a=C, echo=False)

        self.angle=np.array([p1,p,p2])

    def transform_matrix(self, delta, phi, w):
        """
        delta, phi, w =
             delta*math.pi/180. , phi*math.pi/180. , w*math.pi/180.
        Generates the rotation matrix
        """
        print 'deprecated def'
        w = w * math.pi/180. # to radian
        d1 = sin(delta) * cos(phi)
        d2 = sin(delta) * sin(phi)
        d3 = cos(delta)
        cw = cos(w)
        sw = sin(w)
        p=[[None,None,None],[None,None,None],[None,None,None]]
        p[0][0] = ( 1 - d1**2 ) * cw + d1**2
        p[0][1] = d1 * d2 * ( 1 - cw ) + d3 * sw
        p[0][2] = d1 * d3 * ( 1 - cw ) - d2 * sw
        p[1][0] = d1 * d2 * ( 1 - cw ) - d3 * sw
        p[1][1] = ( 1 - d2**2 ) * cw + d2**2
        p[1][2] = d2 * d3 * ( 1 - cw ) + d1 * sw
        p[2][0] = d1 * d3 * ( 1 - cw ) + d2 * sw
        p[2][1] = d2 * d3 * ( 1 - cw ) - d1 * sw
        p[2][2] = ( 1 - d3**2 ) * cw + d3**2
        return p

    def rot_vectang(self, th, r):
        """ Rotate the given rotation matrix r [ca<-sa] by a 
        random axis with th degree.
        """
        from cs import vector_ang
        import numpy as np
        delta, phi = self.rot_axis()
        v = self.polar2vect(delta, phi)
        rot = vector_ang(u=v, th=th)
        newr = np.dot(rot, r)
        return newr

    def rot_axis(self):
        """
        Random rotation axis generator
        """
        delta = random.uniform(-1.,1.) 
        delta = acos(delta)                 # arc
        phi = random.uniform(-math.pi, math.pi) # rotation
        return delta, phi   #radians...

    def polar2vect(self, delta, phi):
        """ Given delta and phi returns a vector"""
        import numpy as np
        x = cos(phi)
        y = sin(phi)
        z = sin(delta)
        vector = np.array([x,y,z])
        vector = vector / np.linalg.norm(vector)
        return vector
        
    def gaussian(self, w0):
        dp = gauss(mu=0., sigma=w0)
        return dp

    def expo(self, w0):
        dp = random.expovariate(w0)
        return dp

    def lognorm(self, w0):
        dp = random.lognormvariate(mu=0, sigma=w0)
        return dp

    def normal(self,w0):
        dp = random.normalvariate(mu=0, sigma=w0)
        return dp


def miller2euler(hkl, uvw):
    """
    RD // uvw  // 1
    ND // hkl  // 3

    3 x 1 = 2
    """
    from euler import euler
    # mat [ca<-sa]
    mat = miller2mat(hkl, uvw)
    phi1, phi, phi2 = euler(a=mat, echo=False)
    return phi1, phi, phi2


def miller2mat(hkl, uvw):
    """
    RD // uvw  // 1
    ND // hkl  // 3

    3 x 1 = 2
    """
    import numpy as np
    from euler import euler
    uvw = map(int, uvw)
    hkl = map(int, hkl)

    uvw = np.array(uvw)
    hkl = np.array(hkl)

    # x, y, z bases vectors.
    cx = uvw/np.linalg.norm(uvw)
    cz = hkl/np.linalg.norm(hkl)
    cy = np.cross(cz, cx)

    mat = np.zeros((3,3))
    for i in range(3):
        mat[i][0] = cx[i]
        mat[i][1] = cy[i]
        mat[i][2] = cz[i]

    # mat [ca<-sa]
    return mat

def miller2mat_RT(uvw,xyz):
    """
    RD // uvw  // 1
    TD // xyz  // 2
    1 x 2 = 3
    """
    import numpy as np
    uvw = np.array(map(int,uvw))
    xyz = np.array(map(int,xyz))

    # x, y, z bases vectors
    cx = uvw/np.linalg.norm(uvw)
    cy = xyz/np.linalg.norm(xyz)
    cz = np.cross(cx,cy)

    mat = np.zeros((3,3))
    for i in range(3):
        mat[i][0] = cx[i]
        mat[i][1] = cy[i]
        mat[i][2] = cz[i]

    # mat [ca<-sa]
    return mat
    


# def miller2euler(hkl=None,uvw=None):
#     """
#     Provided (hkl)// ND and [uvw]//RD,
#     calculates and returns the Euler angles.

#     * Euler angles are using Bunge nomenclature
#     """
#     import numpy as np
#     if hkl==None or uvw==None:
#         print "both hkl and uvw must be given"
#         raise IOError
#     hkl = np.array(hkl)
#     uvw = np.array(uvw)
#     #hkl = hkl/sqrt(((hkl.copy())**2).sum())
#     #uvw = uvw/sqrt(((uvw.copy())**2).sum())
#     print 'uvw = ', uvw
#     print 'hkl = ' , hkl

#     h = hkl.copy()[0]; k = hkl.copy()[1]; l = hkl.copy()[2]
#     u = uvw.copy()[0]; v = uvw.copy()[1]; w = uvw.copy()[2]


#     # phi
#     value = math.sqrt(h**2. + k**2. + l**2.)
#     if value>1 or value<-1:
#         if abs(abs(value)-1)<1.e-10:
#             if value>0: value =1
#             else: value = -1
#         else: raise IOError, 'Unexpected case'

#     phi= acos(value)

#     # p2
#     value = h/math.sqrt(h**2. + k**2.)
#     if value>1 or value<-1:
#         if abs(abs(value)-1)<1.e-10:
#             if value>0: value =1
#             else: value = -1
#         else: raise IOError, 'Unexpected case'
#     p2 = asin(value)

#     p2_= acos(k/math.sqrt(h**2. + k**2.))


#     # p1
#     value = w     / sqrt( u**2. + v**2. + w**2.)
#     value = value * sqrt( h**2. + k**2. + l**2. )
#     value = value / sqrt( h**2. + k**2.)
#     if value>1 or value<-1:
#         if abs(abs(value)-1)<1.e-10:
#             if value>0: value =1
#             else: value = -1
#         else: raise IOError, 'Unexpected case'
#     p1 = asin(value)


#     return np.array([p1,phi,p2])*180./pi
