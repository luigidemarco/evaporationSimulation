""" 
Use distutils to make a module

1) python setup.py build
2) python setup.py install

In the python file, import collision

"""
from setuptools import setup
from distutils.core import Extension

module1 = Extension('collision', sources=['collisionmodule.cpp'])

setup(name = 'CollisionModule',
	version = '1.0',
	description = 'C extension for calculating collision pairs',
	ext_modules = [module1])