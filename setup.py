from __future__ import print_function

import codecs
import os

from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))
about = {}

with codecs.open("README.md", "r") as fh:
    long_description = fh.read()

with codecs.open(os.path.join(here, 'pymerra2', '__init__.py'), 'r', 'utf-8') as f:
    exec(f.read(), about)

INSTALL_REQUIRES = [line.strip() for line in open('requirements.txt')]

KEYWORDS = "nasa merra2 netcdf climate forecast"

setup(
    # -- meta information --------------------------------------------------
    name=about["__title__"],
    version=str(about["__version__"]),
    author=about["__author__"],
    author_email=about["__author_email__"],
    description=about["__description__"],
    long_description=long_description,
    url=about["__url__"],
    license=about["__license__"],
    platforms="all",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.6"
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
    ],
    keywords=KEYWORDS,

    # -- Package structure -------------------------------------------------

    packages=find_packages(exclude=['tests', 'templates']),
    include_package_data=None,
    python_requires='>=2.6, !=3.0.*, !=3.1.*, !=3.2.*, <4',
    install_requires=INSTALL_REQUIRES

)
