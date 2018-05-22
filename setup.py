## Dependents
from distutils.core import setup
from numpy.distutils.core import setup as setup_numpy
from distutils.extension import Extension
# from Cython.Build import cythonize
# from Cython.Distutils import build_ext
from numpy.distutils.core import Extension

cmdclass = {}
ext_modules = []
# ext_modules += [
#     Extension("proj_cy",["src/proj.pyx"] ),]
# cmdclass.update({'build_ext': build_ext})

setup(name='TEXTURElib',
      version='0.0',
      description='Collection of texture related modules',
      author='Youngung Jeong',
      author_email='youngung.jeong@gmail.com',
      packages=['TX'],
      package_dir={'TX':'src'},
      install_requires=['fortranformat','pandas','matplotlib','numpy'],
      cmdclass=cmdclass,
      ext_modules=ext_modules)

## Fortran subroutines with f2py
ext_modules = []
ext_modules += [
    Extension(
        name="pf_for_lib",
        sources=['src/for_lib.for','src/cs.f'],
        ## extra_compile_args=['-O3']
    )]
setup_numpy(ext_modules=ext_modules)
