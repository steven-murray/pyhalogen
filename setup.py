from setuptools import setup, Extension

import os
import sys

version = '0.0.3'

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

# if sys.argv[-1] == "publish":
#     os.system("python setup.py sdist upload")
#     os.system("python setup.py bdist_egg upload")
#     sys.exit()

module1 = Extension('halogen.place_halos',
                    sources=['halogen/place_halos.c'],
                    libraries=['m'],
                    extra_compile_args=["-fopenmp", "-O2"])

setup(
    name="halogen",
    version=version,
    packages=['halogen'],
    install_requires=["hmf"],
    scripts=["scripts/halogen", "scripts/analyse"],
    author="Steven Murray and Santiago Avila Perez",
    author_email="steven.murray@uwa.edu.au",
    description="FAst Synthetic galaxy caTalogues",
    # long_description=read('README.rst'),
    license="MIT",
    keywords="halo mass function 2LPT nbody simulations",
    # url="https://github.com/steven-murray/hmf",
    ext_modules=[module1]
    # could also include long_description, download_url, classifiers, etc.
)
