"""
UXD (BRUKER D8 DISCOVER) to epf(popLA compatible) converter 
developed on (2010-Dec-29)

YOUNGUNG JEONG
   Materials Mechanics Laboratory,
   GIFT, POSTECH

Editor used : GNU Emacs

planned feature 
 1. Background subtraction 
 2. Defocusing correction 
 3. Convert UXD file into EPF 
 4. Make and save defocuse data 
 5. All features need be interactive for later modifications
 6. pole figure monitoring with respect to corrections (bg, dfc)



UXD 2THETA - INTENSITY PROFILE plotting software
  developed on (2011-Jan-25)

  ** UXD 2theta scan plotter with different khis.. 
  ** This job may require 3D plotting with one of three axes being the khi.


***************************************
 Experimental polefigure post-processor
***************************************
Exemplery use,
>> import uxd>> myPF = uxd.pf(filename ='as_R.UXD',,
 echo = False, mode='pf', bgmode='auto',,
 monitoring = False) 
 uxd.pf class has arguments,
 : filename, echo(False, True), mode('pf','df'),,
 bgmode('manual','auto', None)
 >> polefigures = myPF.INTENSITIES
 >> pole0 = np.array(polefigures[0])
 >> uxd.pfsplot(pfs=pole0)     ---> This will show up the pole figure 

***************************************
   2theta(fixed khi)-intensity plotting 
***************************************
>> x,y,khi = uxd.prof(filename='304.uxd')
>> plot(x,y) ... 
"""
#### prints out the documentation to the module.####
print __doc__
####################################################

import numpy as np
import glob
import os
try: from mpl_toolkits.mplot3d import Axes3D
except: print'Axes3D toolkit is not loaded '
try:
    import pylab
    from pylab import *
except: 
    print 'pylab module could not be imported'
sort = np.sort
mean = np.mean


from upf import pfnorm  #pf normalize data

def prof(filename=None):
    """
    Returns x,y, and khi
    """
    FILE = open(filename,'r')
    contents = FILE.read()
    lines = contents.split('_2THETACOUNTS')[1].split('\n')
    info = contents.split('_2THETACOUNTS')[0].split('\n')
    lines = lines[1:len(lines)-1]

    ### 2theta-intensity profile
    for i in range(len(lines)): lines[i] = map(float, lines[i].split())
    lines = np.array(lines)
    lines = lines.transpose()
    x = lines[0]; y = lines[1]
    for i in range(len(info)):
        if info[i][0:5]=='_KHI=':
            khi = float(info[i].split('_KHI=')[1])
            #print 'khi angle is %f'%khi
            break
    return x,y,khi

def profall(preppend):
    files = glob.glob('%s*'%preppend)
    intensities = []
    x = []; y = []; zs = []
    for i in files: 
        tx, ty, kh = prof(filename=i)
        intensities.append(zip(tx,ty))
        zs.append(kh)

    #return intensities, zs
    # make it numpyed
    intensities= np.array(intensities)
    zs=np.array(zs)

    #plotting
    fig = pylab.plt.figure()
    ax = fig.gca(projection='3d')
    fig = pylab.plt.figure()
    ax1 = fig.gca()
    for i in range(len(intensities)):
        th2 = intensities[i].transpose()[0]
        inten = intensities[i].transpose()[1]
        khi = np.resize(zs[i],(th2.shape[0],))
        ax.plot3D(xs=th2, ys=khi, zs=inten)
        ax1.plot(th2, inten)

    ax.set_xlabel(r'$2\theta$', dict(fontsize=20))
    ax.set_ylabel(r'$khi$',dict(fontsize=20))
    ax.set_zlabel('INTENSITY COUNTS',dict(fontsize=12))

    return intensities, zs


def pfsplot(pfs=None, ifig=6, idot=False,x0=0,y0=0): 
    """plot pole figures using matplotlib package """
    try: plt.figure(ifig)
    except :
        print "TclError happened. This is probably because your terminal\n"
        print "does not support DISPLAY as expected. ",
        print "Please try fig.figure(3) and see if is true"
        return -1
    try: plt
    except: 
        print "plt module is not existed"
    print "\n********** pole figure plot internal module ************\n"
    if pf==None: 
        print "\n*****************************************************"
        print " Wrong argument type found. It must be a polefigure "
        print "Will return -1"
        print "*****************************************************\n"
        return -1

    ## loop over individual pole figure in the pfs block
    """ x and y must be a list considering the number of polefigures """
    khi = np.linspace(0,80.,80/5+1)
    phi = np.linspace(0.,360,360/5+1)
    r = khi/90.
    R, PHI = np.meshgrid(r,phi)
    """
    for i in range(len(pfs)):  pass
    pole = pfs[i]
    """
    pole = pfs
    TEMP = np.resize((),(17,73))
    #pole = pole.transpose()       #pole.shape = [17,72]
    print 'pole shape =', pole.shape
    print 'temp shape =', TEMP.shape
    for i in range(len(pole)):
        for j in range(len(pole[i])):
            TEMP[i][j] = pole[i][j]

    for i in range(len(pole)):
        TEMP[i][72]=pole[i][0]
    TEMP = TEMP.transpose()

    x = R*np.cos(PHI*np.pi/180.)
    y = R*np.sin(PHI*np.pi/180.)
    x = x + x0
    y = y + y0

    fig = plt.figure(ifig); tempfig=plt.figure(ifig+1)
    ax = fig.add_subplot(111, aspect='equal')
    tempax = tempfig.add_subplot(111, aspect='equal')
    ax.set_axis_off()

    ### contour     
    if idot==False:
        ## plot and delete the temperal axis
        cnt = tempax.contour(x,y,TEMP)
        maxlev = max(cnt._levels)
        minlev = min(cnt._levels)
        plt.close(tempfig)
        nlev = 10
        lev = np.linspace(minlev,maxlev,nlev)
        ax = fig.add_axes([0.05,0.05,0.89,0.89])
        ax.set_aspect('equal')
        ax.set_axis_off()
        ax.contour(x,y,TEMP,lev)
    elif idot==True:
        ### scatter
        ax.scatter(x,y,c=TEMP)
        pass

    theta = np.linspace(0.,360.,1000)
    Radius = 1.0
    x = Radius * np.cos(theta*pi/180.)
    y = Radius * np.sin(theta*pi/180.)
    x = x + x0
    y = y + y0
    ax.plot(x,y,color='black')
    Radius = 80./90.
    x = Radius * np.cos(theta*pi/180.)
    y = Radius * np.sin(theta*pi/180.)
    x = x + x0
    y = y + y0
    ax.plot(x,y,color='gray',ls='--',alpha = 0.3)


def __pf_plot__(intensity=None, 
                ifig=6, 
                x0=0, y0=0):   #dispose arguments :delt_alpha = 5, delt_khi = 5
    """pole figure plot using matplotlib package """
    try: plt.figure(ifig)
    except :
        print "TclError happened. This is probably because your terminal\n"
        print "does not support DISPLAY as expected. ",
        print "Please try fig.figure(3) and see if is true"
        return -1
    try: plt
    except: 
        print "plt module is not existed"

    print "\n********** pole figure plot internal module ************\n"
    if pf==None: 
        print "\n*****************************************************"
        print " Wrong argument type found. It must be a polefigure "
        print "Will return -1"
        print "*****************************************************\n"
        return -1

    #########################################################
    # Decided to fix the condition to write the code easily.#
    #########################################################

    khi = np.linspace(0,80,80/5+1)
    phi = np.linspace(0.,355,355/5+1)

    #r = np.sin(khi*np.pi/180.)
    r = khi/90.
    R,PHI = np.meshgrid(r,phi)

    x = R*np.cos(PHI*np.pi/180.)
    y = R*np.sin(PHI*np.pi/180.)

    pole = intensity
    pole = np.array(pole)


    ##### 
    phi_temp = np.linspace(0.,360,360/5+1)
    R_TEMP,PHI_TEMP = np.meshgrid(r,phi_temp)
    TEMP = np.resize((),(pole.shape[1],pole.shape[0]))
    print 'pole shape =', pole.shape
    print 'temp shape =', TEMP.shape
    for i in range(len(pole)):
        for j in range(len(pole[i])):
            TEMP[i][j] = pole[i][j]

    for i in range(len(pole)):
        TEMP[i][72]=pole[i][0]

    TEMP = TEMP.transpose()

    x = R_TEMP*np.cos(PHI_TEMP*np.pi/180.)
    y = R_TEMP*np.sin(PHI_TEMP*np.pi/180.)

    x = x + x0
    y = y + y0
    fig = plt.figure(ifig)
    ax = fig.add_subplot(111, aspect='equal')
    
    contour(x,y,TEMP)
    #contourf(x,y,TEMP)  #cotour -filled mode
    plt.show()
    
    ### Circular outer rim drawing
    theta = np.linspace(0.,360.,1000)
    Radius = 1.0
    x = Radius * np.cos(theta*pi/180.)
    y = Radius * np.sin(theta*pi/180.)
    x = x + x0
    y = y + y0
    plt.plot(x,y,color='black')
    Radius = 80./90.
    x = Radius * np.cos(theta*pi/180.)
    y = Radius * np.sin(theta*pi/180.)
    x = x + x0
    y = y + y0
    plt.plot(x,y,color='gray',marker = '--',alpha = 0.3)


def xyz(intensity,ix,iy):
    """ 
    xyz
    Planned to be used, but looks disposed now. 
    """
    if len(intensity)==17:
        xgrid = 5
    else: 
        print 'Unexpected Khi step'
        raise IOError


    if len(intensity[0])==72:
        ygrid = 5
    else:
        print 'Unexpected phi step'
        raise IOError
    

    z = intensity[ix][iy]
    return ix*xgrid,y*ygrid,z


def __near__(a=None, b=[ ]):
    """
    Among elements in the b list who's closed to the value a
    Left-closest, right-closest
    """
    diff = []
    abs_diff = set()
    rst = []
    positive, negative = [], []
    pos , neg  = [], []
    pos_, neg_ = [], []

    for i in range(len(b)):
        if a>b[i]: positive.append(b[i])
        if a<b[i]: negative.append(b[i])

    for i in range(len(positive)):
        pos.append(abs(positive[i]-a))
        pos_.append(abs(positive[i]-a))

    for i in range(len(negative)):
        neg.append(abs(negative[i]-a))
        neg_.append(abs(negative[i]-a))

    pi = pos_.index(sort(pos)[0])
    ni = neg_.index(sort(neg)[0])

    close_pos = positive[pi]
    close_neg = negative[ni]


    rst.append(b.index(close_pos))
    rst.append(b.index(close_neg))

    return rst


def _info(block,name):
    """
    Finds corresponding environmental variable,
    and returns it as string
    """
    for i in range(len(block)):
        try: block[i].split('=')
        except: pass
        else:
            cwl = block[i].split('=')
            #print 'cwl = ', cwl
            if cwl[0][0:len(name)]==name:
                return cwl[1]
    print 'Could not fine ', name,' argument in the block'
    raise IOError


def make_blocks(filename='as_R.UXD', sep = '; (Data for Range number'):
    """
    Fix the separating string as the given sep,
    returns the contents of the file into several blocks
    """
    FILE = open(filename,'r')
    source = FILE.read()
    # for i in range(len(lines)):
    #     if lines[i][0:24]=='; (Data for Range number':
    #         return lines[i+1:len(lines)]
    blocks = source.split(sep)
    return blocks



def block_in(dum):
    """
    splits the bulky dum into lined dum and returns it
    """
    blockl = []
    for i in range(len(dum)):
        blockl.append([])
        blockl[i]=dum[i].split('\n')
    return blockl[0:len(blockl)]


def th2count(block, block_indicator = '_2THETACOUNTS'):
    """
    RETURNS data block's 2th and intensity respectively.
    """
    intensities,phi2 =[],[]
    i = 0
    while True:
        if block[i][0:len(block_indicator)]==block_indicator:
            istart = i; break
        i = i + 1

    while True:
        i = i + 1
        try : block[i]
        except IndexError : return intensities, phi2
        else: 
            if len(block[i])< 3: return intensities, phi2
            try:
                current_intensity = map(float, block[i].split())
            except ValueError:
                #print block[i].split()
                #raw_input()
                pass
            else:
                intensities.append(current_intensity[1])
                phi2.append(current_intensity[0])

            
class pf:
    """
    pole figures 
    Arguments : 
         filename = 'as_R.UXD'
         echo = False
         mode = 'pf'  : flag for defocusing data or pole figure aquisition 
                     if mode =='pf': pole figure aquisition
                     elif mode =='df': defocusing correction
                      -> defocusing data is written down to a files

    note :
       Due to use of 'set' type variable,
       the order of blocks are messed up.
       I've tried to bypass this but failed.
       
       Since non-ordered case is more general,
       in applicability-wise, it is good.
       While, it is bad in coding-wise.

    Reference:
       TEXTURE AND ANISOTROPY, Cambridge University press
        - U.F. KOCKS, C. N. TOME AND H.R. WENK
          ** Chapter 4. Pole Figure Measurements with Diffraction Techniques (Wenk)

    Standard procedure:
       >>> import uxd
       >>> myclass=uxd.pf(filename= 'as_R.uxd', mode='pf')

         uxd.pf(arguments...)
            argument:
              filename,
              echo = FALSE
              mode = 'df, 
                     'df' for defocus curve making,
                     'pf' for post-process the experimental pole figure
              bgmode = 'manual', 'auto'
              sep = "; (Data for Range number)"
     
            *argument sep was added on 2011-01-10 
             sep is a separator in between data blocks.
             This can be different even between uxd files.
             Care must be taken for this.


     *** pole figure projections ***
       >>> mypfs = myclass.INTENSITIES
       >>> uxd.__pf_plot__(intensity = mypfs[0], ifig=9, x0 = 0, y0=0)
       >>> 
              
    """
    
    def __init__(self, filename='as_R.UXD', echo=False, 
                 mode='pf', bgmode='manual', monitoring = False,
                 sep = "; (Data for Range number" ):
        if mode == 'pf':
            print ''
            print '---------------------------'
            print 'POLE FIGURE AQUISITION MODE'
            print '---------------------------\n'
        elif mode =='df':
            print ''
            print '---------------------------'
            print 'DEFOCUS DATA AQUSITION MODE'
            print '---------------------------\n'
        
        blocks = make_blocks(filename=filename, sep = sep)
        self.blocks = blocks
        blocks = block_in(blocks)
        header = blocks[0]
        self.data_block = blocks[1:len(blocks)]
   
        # Trim the data_block neatly.
        for i in range(len(self.data_block)):
            self.data_block[i]=self.data_block[i][0:len(self.data_block[i])-1]
           
        print '** Total number of data blocks: ', len(self.data_block)

        set_2thet = set()
        for i in range(len(self.data_block)):
            cb = self.data_block[i] #cb : current block
            info = self.block_info(block=cb, echo=echo)
            # _2theta, _khi, _steptime, _stepsize
            set_2thet.add(round(float(info[0]), 3))
        th2 = set_2thet.copy()
        #print th2"
        #if raw_input()=='q': raise IOError"

        print 'Kinds of _2theta are printed out'
        while True:
            try:
                print set_2thet.pop(),' ',
            except: break

        self.listh2 = []
        while True:
            try: self.listh2.append(th2.pop())
            except: break

        self.pfs = []
        for i in range(len(self.listh2)):
            self.pfs.append(self.lookfor(th2=self.listh2[i]))

        ###################################
        #   polefigures and backgrounds   #
        ###################################
        self.polefigures=[]
        self.backgrounds=[]
        for i in range(len(self.pfs)):
            if self.bg_or_pf(self.pfs[i], condition = 'digits')=='pf':
                self.polefigures.append(self.pfs[i])
            elif self.bg_or_pf(self.pfs[i], condition = 'digits') =='bg':
                self.backgrounds.append(self.pfs[i])


        print '\n'
        for i in range(len(self.polefigures)):
            print 'PF #',i+1
            _2th, _st, _sz, d_alpha, d_khi, _khis = self.pf_info(self.polefigures[i])
            print 'peak at Bragg the 2theta of ', round(_2th,3)
            print 'delta alpha  = ', d_alpha,' delta khi = ', d_khi,
            print '    step time :', _st,
            print '    step size :', _sz


        print '\n'
        for i in range(len(self.backgrounds)):
            print 'BG #', i+1
            _2th, _st, _sz, d_alpha, d_khis, _khis = self.pf_info(self.backgrounds[i])
            print 'peak at Bragg the 2theta of ', round(_2th,3)
            if d_alpha !='unknown': print 'delta alpha = ', d_alpha, 
            print ' delta khi   = ', d_khi,
            print '    step time :', _st,
            print '    step size :', _sz

        raw_input(' Press Enter  >>> ')

        self.__pf_selection__()

        raw_input("Press enter if you'd like to proceed >> ")

        if os.name=='nt': os.system('cls')
        elif os.name=='posix': os.system('clear')
            
        print "\n\n***************************************************************"
        print "d_alpha is given to backgrounds, that's most probably because"
        print "the backgrounds are measured only at along a certain phi angle"
        print "In other words, it's been partially measured."
        print "***************************************************************"

        "Combines a certain set of pole figure and backgrounds"

        #----------------------------------------------------------
        # Recommends sets of polefigure and its background measure
        #----------------------------------------------------------
        self.__pf_bg_sets__(bgmode = bgmode)  
            #access the combination set by self.combi_pf_bg
        """ 
        Note that if bgmode is None, background is not subtracted
        in the final end. This Ipf = Ipf - dI is not performed in the
        next big loop below
        """

        #----------------------------------------------------------
        #  Core of this scripts
        #      Outer-Loop over polefifures
        #      Inner-Loop over its corresponding two backgrounds
        #----------------------------------------------------------

        INTENSITIES =[]
        
        if bgmode==None:
            for i in range(len(self.polefigures)):
                INTENSITIES.append([])
                for j in range(len(self.polefigures[i])):
                    INTENSITIES[i].append([])
                    C_pf = self.polefigures[i][j]      
                    Ipf = th2count(block = C_pf )[0]    #intensity
                    info_temp = self.block_info(C_pf, echo=False)
                    #normalization by the step time
                    for k in range(len(Ipf)): 
                        Ipf[k] = Ipf[k]/float(info_temp[2])
                        INTENSITIES[i][j].append(Ipf[k])
            
            #monitoring upon polefigure
            if monitoring ==True:
                if os.name=='posix': os.system('clear')
                elif os.name=='nt': os.system('cls')
                print " You chose not to subtract the backgrounds"
                print " Now you have %i polefigure(s) "%(len(self.polefigures))
                print " Do you like to plot pole figures ? "
                if raw_input( 'yes/no') == 'yes':
                    for i in range(len(self.polefigures)):
                        __pf_plot__(intensity = polefigures[i], ifig=6+i)
                    else : pass
            else: pass
                    
        elif bgmode!=None:
            for i in range(len(self.polefigures)):
                INTENSITIES.append([])            
                "On each polefigure set"
                iask = True
                for j in range(len(self.polefigures[i])):
                    INTENSITIES[i].append([])                 
                    " On each chi "
                    R_bg = self.backgrounds[self.combi_pf_bg[i][0]][j]
                    L_bg = self.backgrounds[self.combi_pf_bg[i][1]][j]
                    C_pf = self.polefigures[i][j]
                    Ipf  = th2count(block = C_pf)[0]
                    Ibgl = th2count(block = R_bg)[0]
                    Ibgr = th2count(block = L_bg)[0]
                    info_temp = self.block_info(C_pf, echo=False)
                    info_temp = map(float,info_temp[0:4])
                    pf_2th = info_temp[0]
                    pf_khi = info_temp[1]
                    pf_steptime = info_temp[2]
                    pf_stepsize = info_temp[3]
                    info_temp = self.block_info(L_bg, echo=False)
                    info_temp = map(float,info_temp[0:4])
                    L_bg_2th = info_temp[0]
                    L_bg_khi = info_temp[1]
                    L_bg_steptime = info_temp[2]
                    L_bg_stepsize = info_temp[3]
                    info_temp = self.block_info(R_bg, echo=False)
                    info_temp = map(float,info_temp[0:4])
                    R_bg_2th = info_temp[0]
                    R_bg_khi = info_temp[1]
                    R_bg_steptime = info_temp[2]
                    R_bg_stepsize = info_temp[3]
                    "Normalize the intensity by its steptime"
                    for k in range(len(Ipf)):  Ipf[k] =Ipf[k]/float(pf_steptime)
                    for k in range(len(Ibgl)): Ibgl[k]=Ibgl[k]/float(L_bg_steptime)
                    for k in range(len(Ibgr)): Ibgr[k]=Ibgr[k]/float(L_bg_steptime)
                    
                    if bgmode!=None:
                        bglr_len = [len(Ibgl),len(Ibgr)]
                        if any(bglr_len[k] !=len(Ipf) for k in range(2)):
                            print '** Partial background measured **'
                            pass

                    elif bgmode==None: print '** No Background subtraction **'

                    for k in range(len(Ipf)):
                        """
                        If Ibgl and Ibgr were measured at a certain phi, 
                        at different phi then
                        that certain phi is assumed to be that of the only measured 
                        """
                        try: cibgl = Ibgl[k]
                        except IndexError: cibgl = Ibgl[0]
                        try: cibgr = Ibgr[k]
                        except IndexError: cibgr = Ibgr[0]
                        slope = (cibgr - cibgl)/(R_bg_2th-L_bg_2th)
                        dI = slope * (pf_2th - L_bg_2th) + cibgl
                        Ipf[k] = Ipf[k] - dI
                        if Ipf[k] < 0:
                            if iask==True:
                                print 'Caution) Negative value from prior BG subtraction: ',
                                print 'value = ', Ipf[k]
                                print 'Do you want to keep on going(yes)?, or',
                                print "Don't ask this to the end(ignore)"
                                ans = raw_input('Type answer (yes/ignore) >>  ')
                                if len(ans)==0: pass
                                elif ans =='yes': pass
                                elif ans =='ignore': iask = False; pass
                                elif ans =='no': 
                                    print "\n******************************************"
                                    print "There's no 'no' answer here"
                                    print "Negative intensity is physically non-sense"
                                    print "The code raises an error"
                                    print "******************************************"
                                    raw_input()
                                    raise IOError
                                else: raise IOError
                                Ipf[k] = 1
                            elif iask == False:
                                print "If negative intensity will be returned to be 1."
                                Ipf[k] = 1
                        INTENSITIES[i][j].append(Ipf[k])


        self.INTENSITIES = INTENSITIES
        if mode =='pf':
            dff = glob.glob('*.dfc')
            if len(dff)==0:
                print 'You do not have any *.dfc file'
            else:
                print '*********************************************'
                print '%15s\n'%('Available FILES and its _2THETA')
                print '%3s  %15s  %6s  %6s'%('ID','FILENAME','PEAK_AT','COMMENT')
                for i in range(len(dff)):
                    ff = open(dff[i],'r')
                    lines = ff.readlines()
                    try: float(lines[3].split('=')[1])
                    except:
                        print 'Could not get %s file rightly.\n'%(dff[i])
                        print '**************************************'
                        print "Contents of the file is shown as below\n"
                        if os.name =='posix': os.system('cat %s'%(dff[i]))
                        elif os.name =='nt': os.system('type %s'%(dff[i]))
                        print '\n Please type enter to proceed'; raw_input()
                        pass
                    else:
                        _2th = float(lines[3].split('=')[1])
                        comment = lines[2]
                        ff.close()
                        print '%3i  %15s  %5.3f  %s'%(i, dff[i], _2th, comment)

                for i in range(len(self.INTENSITIES)):
                    print '  Type the defocus correction file id (from 0)'
                    print '  minus value (e.g. -1) will turn down ',
                    print 'the defoucssing correction'
                    id_dfc = raw_input(' >>>    ')
                    id_dfc = int(id_dfc)
                    if id_dfc<0: pass
                    else:  self.__dfc_correct__(filename=dff[id_dfc], pf_id=i)

            #--------------------------------------
            # Normalize the intensity so that max
            # intensity in the all pole figures in
            # the given file equals to 9999 (modified from 999 to 9999)
            #--------------------------------------
            self.__normalize___()

            #--------------------------------------
            # WRITING ACTIVITY TO POLE FIGURE FILE
            #--------------------------------------
            if os.name=='nt': os.system('cls')
            elif os.name=='posix': os.system('clear')
            print "############################"
            print "      WRITING ACTIVITY     "
            print "############################\n"
            print " available formats: epf(0), list(1), No file writing(-1)"
            fileformat = raw_input("Type format flag(0,1,-1)(default=0) >>    ")
            if len(fileformat)==0: fileformat='epf'
            else:
                if fileformat =='0': fileformat='epf'
                elif fileformat =='1': fileformat='list'
                elif fileformat =='-1': fileformat=None
                else: print 'Wrong fileformat input'; raise IOError
            if fileformat == None: pass
            else:
                print " Type the file name"
                filename = raw_input(" >> ")
                self.write(filename=filename, mode=fileformat)
                
        elif mode=='df':  #defocus correction curve making mode
            #--------------------------------
            # DEFOCUS CORRECTION CURVE FILE 
            #--------------------------------
            print "\n\n ****************************************"
            print " * Defocus correction curve file maker *"
            print " ****************************************\n"
            for i in range(len(self.polefigures)):
                self.defc(pfid=i, filename='dfc_'+str(i)+'.dfc', mode='avg')


    def __dfc_correct__(self, filename, pf_id):
        """
        Correct the defucussed error for given self.INTENSITIES[pf_id]
        """
        current_int = self.INTENSITIES[pf_id]
        print 'Defocusing file name is =', filename
        #df_factor = 
        FILE = open(filename,'r')
        
        lines = FILE.readlines()
        lines = lines[5:len(lines)]
        FILE.close()

        temp = []
        for i in range(len(lines)):
            try:  temp.append(float(lines[i].split()[-1]))
            except: pass
        df_factor = temp

        if len(df_factor)==len(current_int): pass
        else: 
            print 'len(df_factor), len(current_int)', len(df_factor),len(current_int)
            raw_input()
            print "defocus file and given polefigure have different khi range"
            raise IOError

        for i in range(len(current_int)): #khi
            for j in range(len(current_int[i])): #phi
                self.INTENSITIES[pf_id][i][j] = current_int[i][j] * df_factor[i]


    def defc(self, pfid=None, 
             filename='dfc_file.dfc', 
             mode='avg'):
        """
        Save file for defocusing correction 
        from the random texture sample.

        >> import uxd
        >> pf = uxd.pf(filename='random.UXD')
               .
               .
               .
        >> pf.defc(pfid=0, filename='my_DFC.dfc')


        Arguments:
           pfid = None (starts from 0)
           filename = 'defocusing_file'
           mode = 'avg', 'all', 'RD'*
                  *If is 'RD' then records only at along phi=0 

        Preliminary:
           You should subtract the background in advance,
           which indicates that you must measure background level as well.

           This method takes the global variable, self.INTENSITIES,
           as the ground state. Among possible pole figures, takes
           one, in that self.INTENSITIES[pfid] is taken.
        """

        curve = []
        if pfid==None:
            print 'You should input pole figure id',
            print ' for your self.polefigures[i]'
            return -1

        FILE = open(filename,'w')
        print '\nWrite your own comment which will be prepended '
        print 'to your defocus file *.dfc'
        comments = raw_input( ' >> ')
        I_stand = self.INTENSITIES[pfid]

        info = self.pf_info(self.polefigures[pfid])
        _2theta = info[0]
        _d_khi = float(info[4])

        FILE.writelines('** XRD DEFOCUS CORRECTION FILE **\n')
        FILE.writelines('** mode: %s\n '%(mode))
        FILE.writelines(comments+'\n')
        FILE.writelines('_2THETA=%8.3f\n'%(_2theta))
        FILE.writelines('  Khi  Intensity    DF_factor\n')
        
        for i in range(len(I_stand)):  # khi
            if mode== 'RD':
                # conventional case
                FILE.writelines(' %5i %5i\n'%(_d_khi*i, I_stand[i][0])) #saves only at along phi=0 (conventionally along RD)

            elif mode=='all':
                # along every and each one of phi axes
                FILE.writelines(' %5i'%(_d_khi*i))
                for j in range(len(I_stand[i])):
                    FILE.writelines(' %5i'%(I_stand[i][j]))
                FILE.writelines('\n')

            elif mode =='avg':
                cwl = []
                FILE.writelines(' %5i'%(_d_khi*i))
                for j in range(len(I_stand[i])):
                    cwl.append(I_stand[i][j])
                curve.append(mean(cwl))
                FILE.writelines(' %5i      %7.3f\n'%
                                (mean(cwl), float(curve[0])/float(mean(cwl))))
                
        pass

                

    def __pf_selection__(self):
        print " \n\n"
        print " ****************************************"
        print " *        POLE FIGURE SELECTION         *"
        print " * Select only some of the pole figures *"
        print " ****************************************\n\n"

        print "**********************************"
        print "*     list of pole figures       *"
        print "**********************************\n"


        print "%8s %7s\n"%('PF id', '2theta')
        pf_slc = []
        for i in range(len(self.polefigures)):
            _2th, _st, _sz, d_alpha, d_khi, _khis = self.pf_info(self.polefigures[i])
            print "%8i %7.3f\n"%(i,_2th)
        
        print "\nPlease type the id of your polefigure\n"
        print "with delimiter as ',' e.g. 0,1,2\n"
        pf_slc_id =map(int, raw_input(' >  ').split(','))

        temp_pf = []
        for ID in pf_slc_id:
            temp_pf.append(self.polefigures[ID])
    
        self.polefigures = temp_pf


    def __normalize___(self):
        """
        Averages out the intensities so that the bigest intensity can be below 9999

        Normalized when the maximum intensity measured does exceed 9999, otherwise
        it does not do .
        """
        #tmx = max(self.INTENSITIES)
        tmx = []
        for i in range(len(self.INTENSITIES)):
            for j in range(len(self.INTENSITIES[i])):
                tmx.append(max(self.INTENSITIES[i][j]))
        mx = max(tmx)
        
        if mx > 9999:
            for i in range(len(self.INTENSITIES)):
                for j in range(len(self.INTENSITIES[i])):
                    for k in range(len(self.INTENSITIES[i][j])):
                        self.INTENSITIES[i][j][k] = self.INTENSITIES[i][j][k]*9999/mx
        else: pass
        
        tmx = []
        for i in range(len(self.INTENSITIES)):
            for j in range(len(self.INTENSITIES[i])):
                tmx.append(max(self.INTENSITIES[i][j]))
        mx = max(tmx)
        
        if mx > 9999.499: 
            print 'unexpected result in def _avg_out_'
            print mx
            raw_input()
            #raise IOError

        
    def write(self, filename='temp.pf', mode='list'):
        """
        Writes intensities
        
        filename ='temp.pf'
        mode = 'list', 'epf', ... and so on (only if necessary later on)
        """
        # Writes to a file having
        # one column listing all the raw intensities 
        # as of given at the current stream through.
        if mode=='list':
            FILE = open(filename,'w')
            FILE.writelines('** Header \n')
            for i in range(len(self.INTENSITIES)):
                #FILE.writelines(' %3i th pole figure \n'%(i+1))
                FILE.writelines('--\n--\n--\n')
                for j in range(len(self.INTENSITIES[i])):
                    #FILE.writelines(' %3i th khi \n'%(j+1))
                    FILE.writelines('--\n--\n')
                    for k in range(len(self.INTENSITIES[i][j])):
                        FILE.writelines('%8i\n'%(self.INTENSITIES[i][j][k]))
            FILE.close()

        elif mode=='epf':  
            'Experimental Pole Figure format (in abbreviation epf)'
            'compatible to popLA, LABOTEX'
            try: ext = filename.split('.')[-1]
            except IndexError: filename=filename+'.epf'
            else:
                if ext =='epf': pass
                else : 
                    filename = filename.split('.')[0]
                    filename = filename+'.epf'
                    
            FILE = open(filename,'w')

            header ='** Needs header here'
            for i in range(len(self.INTENSITIES)):
                FILE.writelines('%s\n'%(header))
                # i: ID of pole figure self.polefigures[i]
                info = self.pf_info(self.polefigures[i])
                print 'The Bragg 2theta of pf',i+1,'#',
                print info[0]
                print 'Type the indices of the current polefigure',
                print 'with spacing'
                index = raw_input(' >>>   ')
                index = map(int, index.split())

                khi_inc = info[4]
                khi_max = max(map(float,info[5]))
                phi_inc = info[3]
                phi_max = 355.
                FILE.writelines('(%1i%1i%1i)'%(index[0],index[1],index[2]))
                FILE.writelines('%5.1f%5.1f'%(khi_inc,khi_max))
                FILE.writelines('%5.1f%5.1f'%(phi_inc,360.))
                FILE.writelines('%2i%2i%2i%2i%2i'%(1,1,2,-1,3))
                FILE.writelines('%5i%5i\n'%(99,1))
                
                #for j in range(len(self.INTENSITIES[i])):
                for j in range(int(90/khi_inc) + 1):
                    # j : ID of each khi scan in self.polefigures[i]
                    # corresponding to self.polefigures[i][j]
                    FILE.writelines(' ')
                    for k in range(int(phi_max/phi_inc) + 1): 
                        try:
                            stream = self.INTENSITIES[i][j][k]
                            if int(stream) == 0: stream = 1  #This is for preventing a measured intensity to be zero
                        except IndexError: stream = 0
                        FILE.writelines('%4i'%(stream))
                        if abs(float(k+1)/18.0 - int(float(k+1)/18.)) < 0.0001:
                            if k == phi_max/phi_inc  : pass
                            else: FILE.writelines('\n ')
                        pass
                    FILE.writelines('\n')
                FILE.writelines('\n')
        else:
            print "You typed wrong type flag"; return -1
                

    def __pf_bg_sets__(self, bgmode='manual'):
        """
        Make self.combi_pf_bg list whose dimension will be like [npf * nbg]
        """
        self.combi_pf_bg = []
        bg2ths = []
        th2_diff = []
        print "\n\n ***  POLEFIGURES"
        for i in range(len(self.polefigures)):
            print 'PF #',i,"'s info"
            print 'The Bragg 2theta = ', self.pf_info(self.polefigures[i])[0]

        print "\n ***  BACKGROUNDS "
        for i in range(len(self.backgrounds)):
            print 'BG #',i,"'s info"
            bg2ths.append(self.pf_info(self.backgrounds[i])[0])
            print 'The Bragg 2theta = ', bg2ths[i]
        print '\n\n' 
        
        if bgmode=='manual':
            for i in range(len(self.polefigures)):
                answer = raw_input('Type the '+str(i)+'th PFs bg ids (eg. 0,2)')
                self.combi_pf_bg.append(map(int,answer.split(',')))
        else:
            if bgmode=='auto':
                for i in range(len(self.polefigures)):
                    crr_2theta = self.pf_info(self.polefigures[i])[0]
                    temp = __near__(a=crr_2theta, b=bg2ths)
                    self.combi_pf_bg.append(temp)
            elif bgmode==None: 
                # Even bgmode is give as None behaves as if it were 'auto'
                # Background will not be subtracted in the end, however.
                for i in range(len(self.polefigures)):
                    crr_2theta = self.pf_info(self.polefigures[i])[0]
                    #This is for checking 
                    """
                    temp = __near__(a=crr_2theta, b=bg2ths)
                    self.combi_pf_bg.append(temp)
                    """
                    return None
        return self.combi_pf_bg


    def pf_info(self, pf):
        """
        print information of pole figure block

        phi 
        range of phi, khi
        grid of phi(0~355), khi(0~80)
        stepsize
        steptime
        _khis
        """
        peak_at= []
        _khis  = []
        _stepsize = []
        _steptime = []
        _d_alpha = []
        for i in range(len(pf)):
            info        = self.block_info(block = pf[i])
            _2theta     = float(info[0])
            _khi        = float(info[1])
            if i == 0: temp = _khi
            if i == 1: _d_khi = abs(_khi-temp)
                
            _steptime.append(float(info[2]))
            _stepsize.append(float(info[3]))
            try:  _delta_alpha =float(info[4])
            except: _delta_alpha=info[4]

            _khis.append(_khi)
            peak_at.append(_2theta)
            _d_alpha.append(_delta_alpha)
        
        if self.__isuniq__(peak_at):pass
        else:
            print 'Positions of bragg peaks are not unique ',
            print 'within the given pole figure'
            raise IOError
        if self.__isuniq__(_steptime): pass
        else:
            print '_STEPTIME is not unique ',
            print 'within the given pole figure'
            raise IOError            
        if self.__isuniq__(_stepsize): pass
        else:
            print '_STEPSIZE is not unique ',
            print 'within the given pole figure'
            raise IOError                        
        if self.__isuniq__(_d_alpha):pass
        else:
            print '_d_alpha is not unique ',
            print 'within the given pole figure'
            raise IOError
        return peak_at[0], _steptime[0], _stepsize[0], _d_alpha[0], _d_khi, _khis

            
    def __isuniq__(self,alist):
        """
        Sees if all the elements in a tuple variable is the same or not.
        Returns False: if they are not all the same
                True: if they are all the same
        """
        a = set()
        if len(alist)>1:
            for i in range(len(alist)):
                a.add(alist[i])
            itrial = 0
            
            while True:
                try:
                    a.pop()
                except:
                    break
                else: itrial = itrial + 1

            if itrial == 1: return True
            else: return False
                
        else:
            print "length of given tuple variable must be exceeding 1"
            raise IOError
        

    def bg_or_pf(self, pfs, condition = 'digits'):
        """
        Provides a good guess if a given "pf-like blocks" is a background
        or pole figure

        pfs[i] : i-th block

        condition ='digits':
           when 'digitis' is the condition, those pseudo pole figure blocks 
           having background with 0 subzero digts are determined to be background
              e.g.   if _2theta == 56.000000 -> BACKGROUND
                     if _2THETA == 56.3432   -> POLEFIGURE  

        condition ='chi0_measure':
           when 'chi0_measure' is the condition,those among pseduo pole figure blocks,
           measured only chi=0 are determined to be background
        """
        if condition == 'digits':
            set_2thet = set()
            for i in range(len(pfs)):
                info = self.block_info(block=pfs[i], echo=False)
                set_2thet.add(info[0])
                
                th2 = float(set_2thet.pop())
                if abs(th2 - int(th2)) > 0.0001: return 'pf'
                else: return 'bg'
        if condition == 'short':
            nint = []
            for i in range(len(pfs)):
                cpfs = pfs[i]
                intensities = th2count(block = cpfs)
                nint.append(len(intensities))
            if all(nint[i]==1 for i in range(len(nint[i]))): return 'bg'
            if all(nint[i]>10 for i in range(len(nint[i]))): return 'pf'


    def lookfor(self, th2 = None, echo = False):
        """
        the provided block is given
        echo = False  :flag for if echo on the screen
        th2 = None    :2theta (at which the Bragg condition is satisfied)
        """
        rst = []
        if th2!=None: 
            for i in range(len(self.data_block)):
                cb = self.data_block[i] #cb : current block                
                if abs( float(_info(cb,'_2THETA')) - th2) < 0.1:
                    if echo ==True:
                        self.block_info(cb)
                    rst.append(cb)
            return rst

        else: print 'th2 must be given' ; raise IOError

    def block_info(self, block, echo = False):
        """
        Print information of the block on the screen

        block = block
        echo = False   :flag for if echo on the screen
        """
        _2theta   = _info(block,'_2THETA')
        _khi      = _info(block,'_KHI')
        _steptime = _info(block,'_STEPTIME')
        _stepsize = _info(block,'_STEPSIZE')

        # information on intensities
        intensities, alpha = th2count(block=block)
        if len(intensities)>1: delta_alpha = alpha[1] - alpha[0]
        else: delta_alpha = 'unknown'
        if echo == True:
            print "current block info"
            print "** _2THETA   = ", _info(block,'_2THETA')
            print "** _KHI      = ", _info(block,'_KHI')
            print "** _STEPTIME = ", _info(block,'_STEPTIME')
            print "** _STEPSIZE = ", _info(block,'_STEPSIZE')
            print "** length of 2phi intensities = ", len(th2count(block))
            print "_2theta, _khi, _steptime, _stepsize"
            print "len(intensities) = ", len(intensities)
            print "intensities[0] = ", intensities[0]
            print "alpha[0] = ", alpha[0]
            print "delta_alpha = ", delta_alpha
        else: pass

        return _2theta, _khi, _steptime, _stepsize, str(delta_alpha)
    pass # end of class pf

