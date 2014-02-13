"""
popLA epf format khi angle (tilting) rim trimmer
"""
from upf import __epffiletrimmer__ as epfTrim

def __open__(wc='*.epf'):
    """ Open files with the indicated wild card (wc) """
    import glob
    files = glob.glob(extension)
    return files

def __print__(fn):
    """ show the file whose name is fn """
    import os
    if not(os.path.isfile(nf)):
        raise IOError, "Could not find the file"
    if os.name=='posix':
        clf = 'clear'
        fshow = 'cat'
    else: raise IOError, 'unsupported command'
    os.system('%s %s'%(fshow, fn))
    
