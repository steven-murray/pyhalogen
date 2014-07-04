from setuptools import setup, Extension

import os
import sys

version = '0.1.0'

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

# # GET THE DEFINE MACROS
defines = [("ONLYBIG", None),
           ("VERB", None),
           ("NDENS", None)]  # These are defaults...

for define in ["DEBUG", "ULTRADEBUG",
               "NO_EXCLUSION", "NO_MASS_CONSERVATION", "RANKED",
               "MASS_OF_PARTS"]:
    if define in sys.argv:
        defines += [(sys.argv.pop(sys.argv.index(define)), None)]

for define in ["ONLYBIG", "VERB", "NDENS"]:
    if define in sys.argv:
        defines.remove((sys.argv.pop(sys.argv.index(define)), None))

place = Extension('pyhalogen.place_halos',
                    sources=['pyhalogen/source/place_halos.c'],
                    libraries=['m', 'gomp'],
                    extra_compile_args=["-fopenmp", "-O2"],
                    define_macros=defines)

pop_mf = Extension("pyhalogen.pop_mf",
                   sources=["pyhalogen/source/populate_mass_function.c"],
                   libraries=["m", "gomp"],
                   extra_compile_args=["-fopenmp", "-O2"],
                   define_macros=defines)

setup(
    name="pyhalogen",
    version=version,
    packages=['pyhalogen'],
    install_requires=["hmf"],
    scripts=["scripts/pyhalogen"],  # , "scripts/analyse"],
    author="Steven Murray and Santiago Avila Perez",
    author_email="steven.murray@uwa.edu.au",
    description="Fast Synthetic Statistical Galaxy Catalogues",
    # long_description=read('README.rst'),
    license="MIT",
    keywords="halo mass function 2LPT nbody simulations",
    url="https://github.com/steven-murray/pyhalogen",
    ext_modules=[place, pop_mf]
    # could also include long_description, download_url, classifiers, etc.
)
