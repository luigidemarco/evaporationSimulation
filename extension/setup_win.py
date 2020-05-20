""" 
Use distutils to make a module
See https://docs.python.org/2.7/extending/windows.html

1) python setup_win.py build
2) python setup_win.py install

In the python file, import collision

"""
from setuptools import setup
from distutils.core import Extension

module1 = Extension('krbcollision', sources=['cpp/krbcollision.cpp'])

setup(name = 'KRbCollisionModule',
	version = '1.0',
	description = 'C extension for calculating collision pairs for Monte Carlo gas dynamics',
	ext_modules = [module1])