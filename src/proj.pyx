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
    import numpy as np
    #normalization of the miller indices
    pole = pole[0:3].copy()
    pole = pole /  (pole**2).sum()
    a,b,c = pole[:]
    ###  mid-plane projection (z=0)
    if c==1:
        # pole[0] = 0; pole[1]=0; pole[2] = 1
        X=0; Y=0
    else:
        X = a/(c-1)
        Y = b/(c-1)
    return X,Y
