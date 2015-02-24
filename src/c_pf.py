"""
Continuous pole figure for a grain.
"""

def reader(fn='TEX_PH1.OUT', step=100, igr=1):
    """ Returns a particular grain's orientations during deformation """
    conts = open(fn, 'r').read()
    blocks = conts.split('TEXTURE AT STRAIN =')
    blocks = blocks[1::]
    grains = []
    eps = []
    for i in range(len(blocks)):
        lines = blocks[i].split('\n')
        e = float(lines[0].split()[-1])
        eps.append(e)
        gr = map(float, lines[4+igr-1].split())
        grains.append(gr)

    return grains[::step], eps

def readerstep(fn='TEX_PH1.OUT', step=0, mode='texture'):
    """ Returns a snapshot of orientation at a certain deformation step """
    import numpy as np
    conts = open(fn, 'r').read()
    if mode=='texture':
        blocks = conts.split('TEXTURE AT STRAIN =')
    elif mode=='morph':
        blocks = conts.split('MORPHOLOGY AT STRAIN =')
    blocks = blocks[1::]
    blocks[step]
    lines = blocks[step].split('\n')
    e = float(lines[0].split()[-1])

    print 'strains:', e
    lines = lines[4:-1:]

    grains = []
    for i in range(len(lines)):
        gr = map(float, lines[i].split())
        grains.append(gr)
    return np.array(grains)

def main(fn='TEX_PH1.OUT', step=100, igr=1,ifig=3):
    """
    Given the texture file, plots a continuous trace of
    texture evolution in the pole figure representation.

    Arguments
    =========
    fn   = 'TEX_PH1.OUT'
    step = 100
    igr  = 1
    ifig = 3
    """
    import upf
    grains, eps = reader(fn, step, igr) # particular grain
    mypf = upf.polefigure(csym='cubic', grains=grains)
    mypf.pf(pole=[[0,0,1]], mode='trace',ifig=ifig)


def progressive_texture(steps=[0,1]):
    import upf
    import matplotlib.pyplot as plt
    plt.ioff()
    for i in range(len(steps)):
        step = steps[i]
        gr=readerstep(fn='TEX_PH1.OUT', step=step, mode='texture')
        mypf=upf.polefigure(csym='cubic',grains=gr)
        mypf.pf(pole=[[0,0,1]],ifig=i+1)
        plt.figure(i+1).savefig('pf_%i.pdf'%(steps[i]+1))
        plt.close(i+1)
        print 'i th pf saved'
