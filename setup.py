from setuptools import setup, Extension, find_packages
import glob

csrc = glob.glob("src/*.c")
srcs = []
for src in csrc:
    if src not in ["src/main.c", "test_code.c"]:
        srcs.append(src)
srcs += ['cosmocalc/cosmocalc_wrap.c']
print(srcs)
cosmocalc_module = Extension(
    'cosmocalc._cosmocalc',
    sources=srcs,
    extra_link_args=["-lfftw3", "-lm", "-lgsl", "-lgslcblas"],
)

setup(
    name='cosmocalc',
    version='0.2.1',
    author="Matthew R. Becker",
    description="cosmocalc wrapped by SWIG",
    ext_modules=[cosmocalc_module],
    packages=find_packages(),
)
