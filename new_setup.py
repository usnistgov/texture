from setuptools import setup, Extension

setup(name='TEXTURE_personal_lib',
      version='0.0',
      description='Collection of texture related modules',
      author='Youngung Jeong',
      author_email='youngung.jeong@gmail.com',
      packages=['TX'],
      package_dir={'TX': 'src'},
      install_requires=['fortranformat'])

from numpy.distutils.core import setup as setup_numpy

## Fortran subroutines with f2py
ext_modules = []
ext_modules += [
    Extension(
        name="pf_for_lib",
        sources=['src/for_lib.for', 'src/cs.f'],
        ## extra_compile_args=['-O3']
    )]
setup_numpy(ext_modules=ext_modules)
