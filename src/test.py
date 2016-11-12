## test..
import MP.lib.whichcomp
path_to_vpsc = MP.lib.whichcomp.find_vpsc_repo()

def test1():
    import upf;reload(upf)
    import matplotlib.pyplot as plt
    import os

    ex_text_filename = os.path.join(
        path_to_vpsc,'examples',
        'ex15_FLD','magnesium','az31','az31dirk.cmb')

    mypf =upf.polefigure(
        filename=ex_text_filename,
        csym='hexag',cdim=[1,1,1.6235])

    fig  = mypf.pf_new(
        poles=[[0,0,0,2],[1,0,-1,1],[1,1,-2,1]],nlev=5,
        mode='line',cmap='jet',ires=True,mn=0.5,
        dth=15,dph=15)

    fig  = mypf.pf_new(
        poles=[[0,0,0,2],[1,0,-1,1],[1,1,-2,1]],nlev=5,
        mode='fill',cmap='jet',ires=True,mn=0.5,
        dth=15,dph=15)

    fig  = mypf.pf_new(
        poles=[[0,0,0,2],[1,0,-1,1],[1,1,-2,1]],nlev=5,
        mode='fill',cmap='rainbow',ires=True,mn=0.5,
        dth=15,dph=15)

if __name__=='__main__':
    test1()
