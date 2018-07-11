# pymerra2

[![Build Status](https://travis-ci.org/Ouranosinc/pymerra2.svg?branch=master)](https://travis-ci.org/Ouranosinc/pymerra2)

## Usage

Examples are found in the template files.

Before beginning, be sure to read the GES-DISC how-to which explains the method of entering your credentials for data access!

https://disc.gsfc.nasa.gov/information/howto?title=How%20to%20Download%20Data%20Files%20from%20HTTP%20Service%20with%20wget

The second step is to find the variables of interest in:

https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf

Note that there are often both instantaneous fields and time-averaged fields for the same variable.

The collection name and the esdt are in the name of each section of that document, e.g. collection = tavg1_2d_flx_Nx and esdt = M2T1NXFLX.
However, the actual esdt directory on the server also has a suffix, that can be verified on the servers themself:

https://goldsmr4.gesdisc.eosdis.nasa.gov/data/
https://goldsmr5.sci.gsfc.nasa.gov/data/

so we set esdt_dir = M2T1NXFLX.5.12.4

These information can be given to the various functions as a dictionary,
or some default values can be set in merra2_variables.py.
CF Standard Names can be found at http://cfconventions.org/standard-names.html

## Limitations

- Chunk sizes are currently set rather arbitrarily.
- Naming files is forced to some common cases, not fully customizable.
- Constant fields require manually downloading a file first.
- Daily NetCDF functions are slightly modified copies of subdaily functions, should be combined into one with calls to time frequency dependent subfunctions.
