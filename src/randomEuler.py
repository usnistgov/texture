"""
Generates random texture and combines with ODF descritized

Practical extension(s) are located down below the class declaration

"""

from math import *
from random import *
import numpy as np

# what's wrong?
class randomEuler:
    """
    A random euler angle generator
    Arguments :
    d_ang : increment of a particular angle among three Euler angles
    p1 : Maximum Phi1 , (usually 360)
    p2 : Maximum Phi2,  (usually 90)
    
    Available number of grains
    100, 300, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000
    """
    def __init__(self, d_ang=10., p1=360., p2=90., ngrain = None, echo=True):
        temp = []
	if ngrain==None: 
            if d_ang==None: dang = 10.
	    else: dang = d_ang
        elif ngrain==100: dang = 26.5
        elif ngrain==200: dang = 21.140625
        elif ngrain==300: dang = 18.35
	elif ngrain==500:  dang = 15.5
	elif ngrain==1000: dang = 12.3
	elif ngrain==2000: dang = 9.753
	elif ngrain==3000: dang = 8.525
	elif ngrain==4000: dang = 7.742
	elif ngrain==5000: dang = 7.187
	elif ngrain==6000: dang = 6.765
	elif ngrain==7000: dang = 6.425
	elif ngrain==8000: dang = 6.146
	elif ngrain==9000: dang = 5.908
	elif ngrain==10000: dang =  5.705
	else:
            print 'Wrongly designated number of grains'
            print "## Using an optimization scheme can be beneficial"
            raise IOError
	i = 0
        while True:
            self.__initial__(d_ang=dang, p1=p1, p2=p2, echo=echo)
            if ngrain == None:  break
            else:
                if self.ngr==ngrain:
                    temp = np.array(temp)
		    temp.std(); temp.mean()
                    print i,'th iteration'
		    if i==0: pass
		    else:
                        print 'Standard deviaition: %9.3f'%temp.std()
			print 'Mean: %9.3f'%temp.mean()
			print 'Variance: %9.3f'%temp.var()
			print 'Max Ngr : %9.3f'%temp.max()
			print 'Min Ngr : %9.3f'%temp.min()
			#raw_input(' press enter ... >>> ')
		    break
		else:
		    temp.append(self.ngr)
            i = i + 1
                
    def __initial__(self, d_ang = 10, p1=360, p2=90, echo=True):
        """
        """
        self.ngr = 0
        self.d_ang = d_ang
        self.p1, self.p2,self.p3 = p1,p2,90
        self.euler = []
        self.nomenclature = 'Bunge'
        self.coef = atan(1.0)/45.
        coef_ = self.coef
        istep = (p1  - 0) /   self.d_ang  #phi1
        jstep = (p2  - 0) /   self.d_ang  #phi2
	kstep = (1.0 - 0) / (coef_*d_ang) #PHI
	itot = 0
	for i in range(int(istep)+3):
            for j in range(int(jstep)+3):
                for k in range(int(kstep)+3):
                    phi1=self.d_ang*i+random()*self.d_ang
                    if phi1<=self.p1:
                        phi2= self.d_ang*j+random()*self.d_ang
                        if phi2<=self.p2:
                            phi = self.d_ang*coef_*k+random()*self.d_ang*coef_
                            if phi<=1:
                                itot=itot+1
                                self.euler.append([itot,
                                                   phi1,
                                                   180./pi*acos(phi),
                                                   phi2])
                                 #self.euler.append([itot,phi1,180./pi*acos(phi),phi,phi2])
        self.ngr = itot
        if echo==True:
            print '*********************************'
            print '* Random Euler angle generator  *'
            print '* ',self.ngr,' grains are generated', '  *'
            print '*********************************'
            pass

    def chg_nomen(self, name):
        """
        chg_nomen is for changing the nomenclature of the Euler angle /n
	*************************************************** /n
	*     The initial nomenclature given is Bunge     * /n
	*  One may change this either to 'Kocks' or 'Roe  * /n
	*************************************************** /n
	"""
	print '*  the current nomenclature is', name
	self.nomenclature = name

    def write(self, filename='random.tex'):
	"""
        write def for writing euler angles into a file 
	in which a name is given as the argument
        Argument:
        filename ='random.tex'
        """
	f = open(filename, 'w')
	header = 'Randomeuler angle writer from randomeuler.py script\n'
	header = header + 'youngung.jeong@gmail.com'+'\n'
	header = header + 'Materials Mechanics Laboratory, GIFT, POSTECH\n'
	header = header + self.nomenclature[0]
	header = header + '   '+ str(self.ngr) + '\n'
	f.write(header)
	for i in range(len(self.euler)):
            for j in range(len(self.euler[0])-1):
		f.write('%8.5f '%(self.euler[i][j+1]))
	    f.write('%8.1f\n'%(10.))
        f.close()
			
    def odf_reading(self, filename_='C:\Python26\odf.txt'):
        """
        Discritized ODF output from LABOTEX, the texture analyser
        Arguments : filename_
        """
        f = open(filename_,'r')
        self.filename = filename_
        dummy=f.readline()
        self.data=[]
	while True:
            s=f.readline()
	    if len(s) < 5:
                print 'an EOF is reached'
                break #EOF
            self.data.append(map(float,s.split()))
        self.maxdata=[]
        self.maxdata.append(self.data[-1][0])  #phi1
        self.maxdata.append(self.data[-1][1])  #phi2
        self.maxdata.append(self.data[-1][2])  #phi

        print 'self.maxdata',
        print self.maxdata
        raw_input()
        
        self.inc = self.data[1][0]-self.data[0][0]
                
    def odf(self, phi1=0, phi2=0, phi=0):
        """
        Arugments : phi1, phi2, phi
        Returns odf(phi1,phi2,phi)
        """
        return self.ODF(phi1,phi2,phi)
            
    def ODF(self, phi1=0, phi2=0, phi=0):
        """
        Returns the ODF of the given Euler angle set.
	It serves as a search engine for a given Euler angle set
        """
        randEul=[]
        randEul.append(int(self.p1))
        randEul.append(int(self.p2))
        randEul.append(int(self.p3))

        if randEul == self.maxdata:
            value = 0
            try:
                self.data[0]
            except:
                print 'ODF does not exist. Please try odf_reading'
                value = -1
            else:
                inc  = self.inc
                pmax = self.maxdata
                isize = int(pmax[0]/inc)+1   #72  for phi1
                jsize = int(pmax[1]/inc)+1   #19  for phi2
                ksize = int(pmax[2]/inc)+1   #19  for phi

                i = phi1/inc
                j = phi2/inc
                k = phi /inc

                a = k * isize * jsize 
                b = j * isize
                c = i
                ind = int(a)+ int(b) + int(c)
                value = self.data[ind][3]
            return value
        else:
            #print 'Euler angle range of the random texture',
            #print ' and ODF is not matched'
            value = 0
            try:
                self.data[0]
            except:
                print 'ODF does not exist. Please try odf_reading'
                value = -1
            else:
                inc  = self.inc
                pmax = self.maxdata
                isize = int(pmax[0]/inc)+1   #72  for phi1
                jsize = int(pmax[1]/inc)+1   #19  for phi2
                ksize = int(pmax[2]/inc)+1   #19  for phi
                while(phi1 > pmax[0]):
                    phi1 = phi1 - pmax[0]
                while(phi2 > pmax[1]):
                    phi2 = phi2 - pmax[1]
                while(phi  > pmax[2]):
                    phi  = phi  - pmax[2]
                i = phi1/inc
                j = phi2/inc
                k = phi /inc
                a = k * isize * jsize 
                b = j * isize
                c = i
                ind = int(a)+ int(b) + int(c)
                value = self.data[ind][3]
	    return value
        return -1
    
    def combine(self, filename='example.cmb'):
        """
	Combines the generated random file and the read ODF file.
	Arguments
	filename = 'example.cmb'
	"""
	try:
	    self.data[0]
	except:
            print 'ODF does not exist. Please try odf_reading'
	    return -1
	else:
            """
	    temp  =[]
	    temp_ =[]
	    for i in range(3):
	    temp.append(int(self.maxdata[i]))
	    temp_.append(int(self.p1))
	    temp_.append(int(self.p2))
	    temp_.append(int(self.p3))
	    """
	    print self.inc
	    e = self.euler
	    data = []
	    f = open(filename,'w')
	    f.write('discretized grain file\n')
	    f.write('dummy\n')
	    f.write('dummy\n')
	    f.write(self.nomenclature[0]+'   '+str(self.ngr)+'\n')            
	    for i in range(len(e)):
                phi1 = []
		phi2 = []
		phi  = []
		down = []
		up   = []
		r    = [ [],[],[] ]
		p    = [ [],[],[] ]
		p[0] = e[i][1]  #phi1
		p[1] = e[i][2]  #phi
		p[2] = e[i][3]  #phi2
		for j in range(3):
                    r[j].append(e[i][j+1] - e[i][j+1]%self.inc)
		    r[j].append(e[i][j+1] - e[i][j+1]%
                                self.inc + self.inc)
		value = self.interpolate(p[0],p[2],p[1])
		data.append(value)
		f.write('%8.5f  '%(p[0]))       #phi1
		f.write('%8.5f  '%(p[1]))       #phi
		f.write('%8.5f  '%(p[2]))       #phi2
		f.write('  %13.8e\n'%(data[i]))
                pass
            pass
	pass
        
    def interpolate(self,phi1=0,phi2=0,phi=0):
        """
        Interpolates the ODF intensity for the given Euler angle set
        """
        r=[ [],[],[] ]
        r[0].append(phi1 - phi1%self.inc)             #phi1
        r[0].append(phi1 - phi1%self.inc  +self.inc)
	r[1].append(phi2 - phi2%self.inc)             #phi2
	r[1].append(phi2 - phi2%self.inc  +self.inc)
	r[2].append(phi  - phi %self.inc)             #phi
	r[2].append(phi  - phi %self.inc  +self.inc)
	value = 0
	x,y,z = 0,0,0
	for j in range(2):              #phi1
	    for k in range(2):          #phi2
                for l in range(2):      #phi
                    x = self.inc - abs(phi1 - r[0][j])
                    y = self.inc - abs(phi2 - r[1][k])
		    z = self.inc - abs(phi  - r[2][l])
		    value = value + self.ODF(r[0][j],
					     r[1][k],
					     r[2][l])*x*y*z/(self.inc**3)
                    pass
                pass
            pass
	return value

##
## practical extensions of using the above class
## ex1) sampling over many typified number of grains
def sampling(odf=None, p1=360., p2=90., *args):
    """
    Sampling several representative grain sets.
    
    ---------
    Arguments
     odf = None : discrete OD file name
     *args: number of grains. (500,~)

    -------
    Example
    In order to obtain 500, 1000, 2000, 4000, RVEs,
    follow the below:
    >>> import randomEuler
    >>> randomEuler.sampling('304_surface.TXT', 500,1000,2000,4000,...)
              . . .
    """
    import upf
    files = []
    print 'p1=', p1, 'p2=', p2
    for i in range(len(args)):
        re = randomEuler(ngrain=args[i], p1=p1, p2=p2)
	re.odf_reading(odf)
        filename = "%s_%s%s"%(odf.split('.TXT')[0],str(args[i]).zfill(5),'.cmb')
	re.combine(filename)   #"%s_%s%s"%(odf.split('.TXT')[0],str(args[i]).zfill(5),'.cmb'))
        files.append(filename)
        pass

    ## pole figure plotting
    flag = raw_input('Polefigure? (bcc or fcc)')
    if flag=='bcc' or 'fcc':
        for i in range(len(files)):
            pf = upf.polefigure(filename=files[i], csym='cubic',
                                cdim=[1.,1.,1.], cang=[90.,90.,90.])
            pf.pf(pole=[[1,1,1]], mode='contourf', ifig=1+i, cmode=None)
            pass
        pass
    else: pass
    pass
## print document of def sampling, for it is frequnetly used.
# print sampling.__doc__

## Finds a relavant incremental angle for a given ngrain
## using an optimization scheme.
def __refunc__(dang=[10,0], ngrain=200, p1=360., p2=90.):
    """
    The function for which the parameters are optimized 
    """
    dang = dang[0] #the rest be the dummy
    a = randomEuler(d_ang=dang, p1=p1, p2=p2, echo=False)
    a.ngr # number of grain
    a.euler = np.array(a.euler)
    return abs(a.ngr - ngrain)

def __callback__(dang=None):
    """
    A callback function
    """
    # import matplotlib.pyplot as plt
    dang = dang[0]
    # a = randomEuler(dang, echo=False)
    print 'Incremental angle =', dang
    # print 'number of grain =', a.ngr
    # ax = plt.gca()
    # ax.plot(dang,'.')
    # plt.show()

    pass
    
def sampling_simplex(odf=None, iang = 10.,
                     ngrain=100,maxiter=40,xtol=10,
                     header = None, p1=360., p2=90.,
                     iplot=False,ifig=30):
    """
    RVE sampling over an OD file using an optimization scheme

    argument
    odf = grided crystallograhpic orientation distribution file
    iang = initial anglular increment
    ngrain = 100, # the targetted number of grains
    maxiter = 40, maximum number of iteration for which the optimization loop
                  iterates
    xtol = 10, the tolerance
    header = None
    p1 = 360,  (maximum phi1 angle)
    p2 = 90 ,  (maximum phi2 angle)
    """
    import matplotlib.pyplot as plt
    #import upf
    from scipy.optimize import fmin #optimization scheme
    rst = fmin(__refunc__,                #function to be optimized
               [iang,0],                  # the parameters
               args=(ngrain,p1,p2),           # *args
               callback = __callback__, # callback function
               xtol = xtol, ##changeable...
               maxiter=maxiter
               )

    # raw_input('press enter to close >>>')
    # plt.close()
   
    print 'angle =', rst

    grains = []
    while True:
        tmp = randomEuler(d_ang=rst[0], p1=p1, p2=p2, echo=False)
        grains.append(tmp.ngr)
        if tmp.ngr==ngrain:
            break
        pass

    print "-----------------------"
    print "Interation: %i"%len(grains)
    print "mean number: %6.2f"%np.array(grains).mean()
    print "-----------------------"
    
    tmp.odf_reading(odf)
    if header==None: filename= '%s.cmb'%str(ngrain).zfill(5)
    else: filename= '%s_%s.cmb'%(header, str(ngrain).zfill(5))
    tmp.combine(filename)
    print "written to '%s'"%filename


    if iplot==True:
        import upf
        pf = upf.polefigure(filename=filename, csym='cubic',
                            cdim=[1.,1.,1.], cang=[90.,90.,90.])
        pf.pf(pole=[[1,0,0],[1,1,0],[1,1,1]], mode='contourf', ifig=ifig)
        pass
        
    return rst


def sampling_simplex_grains(odf=None,ngrains = []):
    for ng in ngrains:
        rst = sampling_simplex(odf=odf, iang=10.,
                               ngrain=ng, maxiter=40, xtol=10)


