""" 
Use distutils to make a module

1) python setup_linux.py build
2) python setup_linux.py install

In the python file, import collision

"""

from distutils.core import setup, Extension

module1 = Extension('krbcollision', sources=['cpp/krbcollision.cpp'])

setup(name = 'KRbCollisionModule',
	version = '1.0',
	description = 'C extension for calculating collision pairs for Monte Carlo gas dynamics',
	ext_modules = [module1])