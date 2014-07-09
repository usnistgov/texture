"""
popLA epf format to defocus correction


# defocus correction curve generation
  1. read pole figures of isotropic powder samples
  2. Each of the files may include a multiple pole figures
  3. For each of pole, average the intensities along the same khi rims
  4. Correction curve = 1/(averaged_intensity)
"""

def __gen__(ext='*.epf'):
    """ Defocus correction curve generator """
    import matplotlib.pyplot as plt
    import numpy as np
    import glob
    files = glob.glob(ext)
    masterdata = []
    for f in files:
        datasets, max_khi = __read__(f)
        # need to trim up until max_khi though ...
        for i in range(len(datasets)):
            if len(datasets[i])!=19:
                print 'The resolution along chi axis is not 5degree.'
                print 'There is %i elements along the chi axis.'%len(datasets)
                raw_input('>>>')
                raise IOError
            else: dkhi = 5.
            tiny = 0.001
            ang = np.arange(0, 90+tiny, dkhi)

            # maxi khi trim the datasets[i]

            avg = __avg__(datasets[i])
            DefocusCurve = __NormReci__(avg) #avg: a list type variable
            plt.plot(ang, DefocusCurve)
            ax=plt.gca(); ax.set_xlabel(r'$\chi$')
            #print DefocusCurve
    pass

def __read__(fn):
    """ read a epf file """
    from upf import __epffiletrimmer__ as epfTrim
    # Read blocks
    blocks= open(fn, 'rU').read().split('(')[1:]
    npf = len(blocks)
    datasets = []
    max_khi = []
    for i in range(len(blocks)):
        hkl = blocks[i][0:3] #hkl
        hkl = map(int, [hkl[0],hkl[1],hkl[2]])
        if blocks[i][3]!=')':
            print 'Caution: unexpected hkl labeling format'
            pass
        d, mxk = epfTrim(blocks[i]) #only for popLA epf format
        datasets.append(d)
        max_khi.append(mxk)
        pass
    return datasets, max_khi

def __avg__(SingleDataSet):
    """ Returns a list including average intensities
    along the phi axis
    """
    import numpy as np

    AvgThroughPhi = []
    for i in range(len(SingleDataSet)): #each khi
        avg = np.average(SingleDataSet[i])
        AvgThroughPhi.append(avg)
        pass
    return AvgThroughPhi

def __NormReci__(SingleList):
    """ Given the single list of intensity drop lines,
    normalize it with respect to the intensity at khi=0.
    Then 1/ curve is returned.
    """
    import numpy as np
    
    dat = np.array(SingleList)/ SingleList[0]
    rst = [] 
    for i in range(len(dat)):
        if dat[i]!=0: rst.append(1./dat[i])
        else: rst.append(np.NAN)
        pass
    rst = np.array(rst)
    return rst

    

# def __DefocussedIntensities__(datasets):
#     """ 
#     """
#     data = np.zeros((len(datasets), 19, 72))
#     dum = data.copy()
 
#     for i in range(len(datasets)): #each pf
#         for j in range(len(datasets[i])): #each khi
#             for k in range(len(datasets[i][j])):
#                 data[i,j,k] = datasets[i][j][k]
#         ## Swap the axis
#     data = data.swapaxes(1,2) # E.g. 72 x 19 than 19 x 72
    
                


    
    
    
    
