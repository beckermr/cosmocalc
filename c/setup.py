#!/usr/bin/env python

"""
setup.py file for SWIG wrapping of cosmocalc
"""

from distutils.core import setup, Extension
import glob

csrc = glob.glob("*.c")
srcs = []
for src in csrc:
    if src != "main.c" and src != "fftlog.c" and src != "test_code.c":
        srcs.append(src)
        
cosmocalc_module = Extension('_cosmocalc',
                             sources=srcs,
                             extra_compile_args = ["-I/afs/slac.stanford.edu/g/ki/software/gsl/1.15/@sys/include"],
                             extra_link_args = ["-L/afs/slac.stanford.edu/g/ki/software/gsl/1.15/@sys/lib","-lm","-lgsl","-lgslcblas"],
                             )

setup (name = 'cosmocalc',
       version = '0.1',
       author      = "Matthew R. Becker",
       description = """cosmocalc wrapped by SWIG""",
       ext_modules = [cosmocalc_module],
       py_modules = ["cosmocalc"],
       )

