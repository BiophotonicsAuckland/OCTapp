from cx_Freeze import setup, Executable
import os
PYTHON_INSTALL_DIR = os.path.dirname(os.path.dirname(os.__file__))
os.environ['TCL_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tcl8.6')
os.environ['TK_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tk8.6')

additional_mods = ['numpy.core._methods', 'numpy.lib.format', 'numpy.matlib']
import numpy.core._methods
import numpy.lib.format

import scipy
includefiles = [os.path.dirname(scipy.__file__)]
excludes = ["tkinter", "PyQt4.QtSql", "sqlite3", "scipy.lib.lapack.flapack", "PyQt4", "tensorflow", "bokeh", "pandas" , "notebook", "babel", "botocore", "mpl_data", "theano", "netCDF4", "IPython"]


setup( 
name = "script" , 
version = "0.1" ,
description = "test" , 
options = {'build_exe':{'includes':additional_mods, 'include_files':includefiles, 'excludes':excludes, 'optimize':2} },
executables = [Executable("OCTApp.py")] , 
)