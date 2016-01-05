## parsing texture file consisting of multiple blocks.
import numpy as np
from TX import upf
import matplotlib.pyplot as plt

def read(fn='TEX_PH1.OUT'):
    """
    Read blocks from the given file <string>

    Argument
    --------
    fn='TEX_PH1.OUT'

    Returns
    -------
    blocks
    """
    with open(fn,'r') as fo:
        blocks=fo.read().split('TEXTURE')[1:]
    return blocks

def block2gr(block):
    """
    Convert a tex_ph block to grain arrays
    """
    lines = block.split('\n')
    ngr = int(lines[3].split()[1])
    l = lines[4:4+ngr]
    grs = np.zeros((ngr,4))
    for i in xrange(len(l)):
        grs[i,:]=map(float,l[i].split())
    return grs

def main(fn='TEX_PH1.OUT',csym='cubic',
         poles=[[1,0,0],[1,1,0],[1,1,1]],
         mn=None,mx=None,nlev=5):
    """
    Arguments
    =========
    fn    = 'TEX_PH1.OUT'
    csym  = 'cubic'
    poles = [[1,0,0],[1,1,0],[1,1,1]]
    mn    = None
    mx    = None
    nlev  = 5

    Save pole figures to a series of files <pf_%i.pdf>
    """
    blocks = read(fn)
    figs=[]
    for i in xrange(len(blocks)):
        gr=block2gr(blocks[i])
        mypf=upf.polefigure(grains=gr,csym=csym)
        fig=mypf.pf_new(ifig=None,poles=poles,mn=None,mx=None,nlev=nlev)
        figs.append(fig)
    return figs

if __name__=='__main__':
    import argparse    
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--fn', type=str, help='Texture File name',
        default='TEX_PH1.OUT')
    parser.add_argument(
        '--csym',type=str, help='Crystal symmetry',
        default='cubic')
    parser.add_argument(
        '--mn', type=float, help='minimum intensity level',
        default=None)
    parser.add_argument(
        '--mx', type=float, help='maximum intensity level',
        default=None)
    parser.add_argument(
        '--nlev', type=int, help='Number of levels',
        default=9)

    args = parser.parse_args()
    fn   = args.fn
    csym = args.csym
    mn   = args.mn
    mx   = args.mx
    nlev = args.nlev

    main(fn=fn,csym=csym,
         poles=[[1,0,0],[1,1,0],[1,1,1]],
         mn=mn,mx=mx,nlev=5)
