import codecs
import os

from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))
about = dict(
    __version__=0.3,
    __title__="pymerra2",
    __description__="A tool for downloading and repackaging NASA MERRA-2 Data",
    __url__="https://github.com/Ouranosinc/pymerra2",
    __author__="Trevor James Smith",
    __author_email__="smith.trevorj@ouranos.ca",
    __license__="Apache Software License 2.0",
    __copyright__="Copyright 2018 Ouranos Inc.",
)

with codecs.open("README.md", "r") as fh:
    long_description = fh.read()

INSTALL_REQUIRES = [line.strip() for line in open("requirements.txt")]

KEYWORDS = "nasa merra2 netcdf climate forecast reanalysis"

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
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "Topic ::Utilities",
    ],
    keywords=KEYWORDS,
    packages=find_packages(exclude=["tests", "templates"]),
    include_package_data=None,
    python_requires=">=3.5, <4",
    install_requires=INSTALL_REQUIRES,
)
