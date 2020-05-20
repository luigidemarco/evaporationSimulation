import platform
import os

plt = platform.system()

if plt == "Linux" or plt == "Darwin":
    os.system("python setup_linux.py build")
    os.system("python setup_linux.py install")
elif plt == "Windows":
    os.system("python setup_win.py build")
    os.system("python setup_win.py install")
