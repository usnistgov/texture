## test..

def test1():
    import upf;reload(upf)
    import matplotlib.pyplot as plt
    mypf =upf.polefigure(
        filename='/Users/yj/repo/vpsc/vpsc-dev-fld/examples/ex15_FLD/magnesium/az31/az31dirk.cmb',
        csym='hexag',cdim=[1,1,1.6235])
    fig  = mypf.pf_new(poles=[[0,0,0,2],[1,0,-1,1],[1,1,-2,1]],nlev=5,
                       mode='line',cmap='jet',ires=True,mn=0.5,
                       dth=15,dph=15)


if __name__=='__main__':
    test1()
