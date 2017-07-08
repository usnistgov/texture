"""

"""
## parsing texture file consisting of multiple blocks.
import numpy as np
from TX import upf
import matplotlib.pyplot as plt

def read(fn='TEX_PH1.OUT'):
    """
    Read blocks from the given file <string>

    Argument
    --------
    fn='TEX_PH1.OUT' ; or tarfile object

    Returns
    -------
    blocks
    """

    if type(fn).__name__ in ['TarFile','ExFileObject']:
        return read_fo(fn)
    elif type(fn).__name__=='str':
        with open(fn,'r') as fo:
            blocks = read_fo(fo)
        return blocks
    else:
        raise IOError, 'unexpected type of fn given;'

def read_fo(fo):
    """
    Read from object.
    """
    return fo.read().split('TEXTURE')[1:]

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

def main(fn='TEX_PH1.OUT',csym='cubic',**kwargs):
    """
    Save pole figures to a series of files <pf_%i.pdf>

    Arguments
    =========
    fn    = 'TEX_PH1.OUT' or tarfile object
    csym  = 'cubic'
    **kwargs: key-worded arguments to upf.pf_new
    """
    blocks = read(fn)
    figs=[]
    for i in xrange(len(blocks)):
        gr=block2gr(blocks[i])
        mypf=upf.polefigure(grains=gr,csym=csym)
        fig=mypf.pf_new(**kwargs)
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
