import platform
import os

plt = platform.system()

if plt == "Linux":
	os.system("cd setup/")
	os.system("python setup_linux.py build")
	os.system("python setup_linux.py install")
elif plt == "Windows":
	os.system("cd setup/")
	os.system("python setup_win.py build")
	os.system("python setup_win.py install")
elif plt == "Darwin":
	print("Sorry, MacOS not yet supported.")