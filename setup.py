## Dependents
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext


cmdclass = {}
ext_modules=[]

ext_modules += [
    Extension("proj_cy",["src/proj.pyx"]),]
cmdclass.update({'build_ext': build_ext})

setup(name='TEXTURE_personal_lib',
      version='0.0',
      description='collection of texture related modules',
      author='Youngung Jeong',
      author_email='youngung.jeong@gmail.com',
      packages=['TX'],
      package_dir={'TX':'src'},
      cmdclass=cmdclass,
      ext_modules=ext_modules
      )
