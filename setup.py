#!/usr/bin/env python

"""
setup.py file for SWIG wrapping of cosmocalc
"""

import os
from distutils.core import setup, Extension
import glob

csrc = glob.glob("*.c")
srcs = []
for src in csrc:
    if src != "main.c" and src != "fftlog.c" and src != "test_code.c":
        srcs.append(src)
srcs.append("cosmocalc.i")        
cosmocalc_module = Extension('_cosmocalc',
                             sources=srcs,
                             extra_compile_args = [os.path.expandvars("-I${SLAC_GSL_INC}"),"-I/opt/local/include"],
                             extra_link_args = [os.path.expandvars("-L${SLAC_GSL_LIB}"),"-L/opt/local/lib","-lm","-lgsl","-lgslcblas"],
                             )

setup (name = 'cosmocalc',
       version = '0.1',
       author      = "Matthew R. Becker",
       description = """cosmocalc wrapped by SWIG""",
       ext_modules = [cosmocalc_module],
       py_modules = ["cosmocalc"],
       )

