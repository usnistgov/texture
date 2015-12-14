from distutils.core import setup


## Libs written in Fortran
#from __future__ import division, absolute_import, print_function
from numpy.distutils.core import Extension
ext1 = Extension(name = 'for_lib',sources = ['src/for_lib.f'])

setup(name='upf',
      version='1.0',
      py_modules=['cmb','upf','sym','randomEuler','epf2dfc','epf_khi_trimmer'],
      description='The ultimate pole figure plotting system on python',
      author ='Youngung Jeong',
      author_email='youngung.jeong@gmail.com',
      ext_modules=[ext1]
      )
