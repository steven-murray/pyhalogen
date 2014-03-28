from setuptools import setup, find_packages, Extension

import os
import sys

version = '0.0.2'

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

# if sys.argv[-1] == "publish":
#     os.system("python setup.py sdist upload")
#     os.system("python setup.py bdist_egg upload")
#     sys.exit()

module1 = Extension('pyfast.place_halos',
                    sources=['pyfast/place_halos.c'],
                    libraries=['m', 'gomp'],
                    extra_compile_args=["-fopenmp", "-O2", "-std=c99"])

setup(
    name="pyfast",
    version=version,
    packages=['pyfast'],
    install_requires=["hmf"],
    scripts=["scripts/runfast", "scripts/analyse"],
    author="Steven Murray and Santiago Avila Perez",
    author_email="steven.murray@uwa.edu.au",
    description="FAst Synthetic galaxy caTalogues",
    # long_description=read('README.rst'),
    license="MIT",
    keywords="halo mass function",
    url="https://github.com/steven-murray/hmf",
    ext_modules=[module1]
    # could also include long_description, download_url, classifiers, etc.
)
