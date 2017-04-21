import os
import glob
import datetime
from calendar import monthrange

import numpy as np
import numpy.ma as ma
import netCDF4

# Aliases for default fill values
defi2 = netCDF4.default_fillvals['i2']
defi4 = netCDF4.default_fillvals['i4']
deff4 = netCDF4.default_fillvals['f4']

subdaily_vars = {'uas': {'esdt_dir': 'M2I1NXASM.5.12.4',
                         'collection': 'inst1_2d_asm_Nx',
                         'merra_name': 'V10M',
                         'standard_name': 'northward_wind'}}


class NetCDFError(Exception):
    pass


def _time_vectors_int(time_vectors, force=False, raise_exception=False,
                      allow_masked=True, dtype='int32'):
    # tries to make the time vectors array an integer array if possible.
    if allow_masked:
        time_vectors_int = ma.array(time_vectors, dtype=dtype)
    else:
        time_vectors_int = np.array(time_vectors, dtype=dtype)
    if force or ((time_vectors_int - time_vectors).sum() == 0):
        return time_vectors_int
    else:
        if raise_exception:
            raise NetCDFError("Floats in time vectors.")
        else:
            return time_vectors


def _datetimes_to_time_vectors(datetimes):
    """Convert list of datetimes to Nx6 matrix.

    Parameters
    ----------
    datetimes - list of datetime

    Returns
    -------
    out - numpy array
        Nx6 matrix of time vectors

    """

    # Does not support microsecond (datetime shape of length 7)
    def datetime_timetuple(one_datetime):
        if one_datetime is None:
            return ma.masked_all([6], dtype='int32')
        return one_datetime.timetuple()[0:6]

    try:
        time_tuples = map(datetime_timetuple, datetimes)
        return _time_vectors_int(ma.array(time_tuples))
    except (AttributeError, TypeError):
        time_tuples = datetimes.timetuple()
        return _time_vectors_int(ma.array(time_tuples))


def subdaily_download(dataset_esdt, merra2_collection, initial_year,
                      final_year, initial_month=1, final_month=12,
                      initial_day=1, final_day=None, output_directory=None):
    # This will not work for constant fields
    # This will not work for monthly data
    if output_directory is None:
        add_output_dir = ''
    else:
        add_output_dir = "--directory-prefix={0} ".format(output_directory)
    merra_cmd = ("wget -c {0}--load-cookies ~/.urs_cookies "
                 "--save-cookies ~/.urs_cookies --keep-session-cookies "
                 "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/")
    merra_cmd = merra_cmd.format(add_output_dir)
    data_path = ("MERRA2/{4}/{0}/{1}/"
                 "MERRA2_{3}.{5}.{0}{1}{2}.nc4")
    for yyyy in range(initial_year, final_year+1):
        if yyyy < 1992:
            merra_stream = '100'
        elif yyyy < 2001:
            merra_stream = '200'
        elif yyyy < 2011:
            merra_stream = '300'
        else:
            merra_stream = '400'
        if yyyy == initial_year:
            mi = initial_month
        else:
            mi = 1
        if yyyy == final_year:
            mf = final_month
        else:
            mf = 12
        for mm in range(mi, mf+1):
            if (yyyy == initial_year) and (mm == mi):
                di = initial_day
            else:
                di = 1
            if final_day and (yyyy == final_year) and (mm == mf):
                df = final_day
            else:
                mrange = monthrange(yyyy, mm)
                df = mrange[1]
            for dd in range(di, df+1):
                cdp = data_path.format(str(yyyy), str(mm).zfill(2),
                                       str(dd).zfill(2), merra_stream,
                                       dataset_esdt, merra2_collection)
                os.system(merra_cmd + cdp)


def subdaily_netcdf(path_data, output_file, var_name, initial_year,
                    final_year):
    # Currently does not deal with the level aspect (and therefore 3d fields)
    search_str = "*{0}*.nc4".format(subdaily_vars[var_name]['collection'])
    nc_files = glob.glob(os.path.join(path_data, search_str))
    nc_files.sort()

    relevant_files = []
    divided_files = [[]]
    nt_division = [0]
    nt = 0
    nmb = 0
    for nc_file in nc_files:
        yyyy = int(nc_file.split('.')[-2][0:4])
        if (yyyy >= initial_year) and (yyyy <= final_year):
            relevant_files.append(nc_file)
            nc = netCDF4.Dataset(nc_file, 'r')
            ncvar = nc.variables[subdaily_vars[var_name]['merra_name']]
            nmb += ncvar.size*4/(1024.*1024.)
            if nmb > 500:
                divided_files.append([nc_file])
                nt_division.append(0)
                nmb = ncvar.size*4/(1024.*1024.)
            else:
                divided_files[-1].append(nc_file)
            nt_division[-1] += len(nc.dimensions['time'])
            nt += len(nc.dimensions['time'])
            nc.close()

    nc_reference = netCDF4.Dataset(relevant_files[0], 'r')
    var_ref = nc_reference.variables[subdaily_vars[var_name]['merra_name']]

    # 2.1 Filename
    #     NetCDF files should have the file name extension ".nc".
    nc_file = output_file

    now = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    nc1 = netCDF4.Dataset(nc_file, 'w', format='NETCDF4_CLASSIC')

    # 2.6.1 Identification of Conventions
    nc1.Conventions = 'CF-1.6'

    # 2.6.2. Description of file contents
    nc1.title = ('Modern-Era Retrospective analysis for Research and '
                 'Applications, Version 2')
    nc1.history = "%s: Convert from original format to NetCDF" % (now,)
    nc1.institution = 'NASA'
    nc1.source = 'Reanalysis'
    nc1.references = 'https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/'

    # Create netCDF dimensions
    nc1.createDimension('time', nt)
    nc1.createDimension('ts', 6)
    # nc1.createDimension('level', k)
    nc1.createDimension('lat', len(nc_reference.dimensions['lat']))
    nc1.createDimension('lon', len(nc_reference.dimensions['lon']))

    # Create netCDF variables
    # Compression parameters include:
    # zlib=True,complevel=9,least_significant_digit=1
    # Set the fill value (shown with the 'f4' default value here) using:
    # fill_value=netCDF4.default_fillvals['f4']
    # In order to also follow COARDS convention, it is suggested to enforce the
    # following rule (this is used, for example, in nctoolbox for MATLAB):
    #     Coordinate Variables:
    #     1-dimensional netCDF variables whose dimension names are identical to
    #     their variable names are regarded as "coordinate variables"
    #     (axes of the underlying grid structure of other variables defined on
    #     this dimension).

    # COADS requirements
    nc1.createVariable('ts', 'i2', ('ts',), zlib=True, fill_value=defi2)

    # 4.4. Time Coordinate
    time = nc1.createVariable('time', 'i4', ('time',), zlib=True)
    time.axis = 'T'
    time.units = "hours since 1980-01-01 00:00:00"
    time.long_name = 'time'
    time.standard_name = 'time'
    # 4.4.1. Calendar
    time.calendar = 'gregorian'

    time_vectors = nc1.createVariable('time_vectors', 'f4', ('time', 'ts'),
                                      zlib=True)

    # 4.3. Vertical (Height or Depth) Coordinate
    # level = nc1.createVariable('level','f4',('level',),zlib=True)
    # level.axis = 'Z'
    # level.units = 'Pa'
    # level.positive = 'up'
    # level.long_name = 'air_pressure'
    # level.standard_name = 'air_pressure'

    # 4.1. Latitude Coordinate
    lat = nc1.createVariable('lat', 'f4', ('lat',), zlib=True)
    lat.axis = 'Y'
    lat.units = 'degrees_north'
    lat.long_name = 'latitude'
    lat.standard_name = 'latitude'
    lat[:] = nc_reference.variables['lat'][:]

    # 4.2. Longitude Coordinate
    lon = nc1.createVariable('lon', 'f4', ('lon',), zlib=True)
    lon.axis = 'X'
    lon.units = 'degrees_east'
    lon.long_name = 'longitude'
    lon.standard_name = 'longitude'
    lon[:] = nc_reference.variables['lon'][:]

    var1 = nc1.createVariable(var_name, 'f4', ('time', 'lat', 'lon'),
                              zlib=True, chunksizes=(720, 30, 30),
                              fill_value=deff4)
    # 3.1. Units
    var1.units = var_ref.units
    # 3.2. Long Name
    var1.long_name = var_ref.long_name
    # 3.3. Standard Name
    var1.standard_name = subdaily_vars[var_name]['standard_name']

    nc_reference.close()

    t = 0
    for i, division in enumerate(divided_files):
        tmp_data = ma.masked_all([nt_division[i], len(nc1.dimensions['lat']),
                                  len(nc1.dimensions['lon'])])
        tmp_time = ma.masked_all([nt_division[i]])
        ttmp = 0
        for nc_file in division:
            print(nc_file)
            nc = netCDF4.Dataset(nc_file, 'r')
            ncvar = nc.variables[subdaily_vars[var_name]['merra_name']]
            nctime = nc.variables['time']
            ncdatetime = netCDF4.num2date(nctime[:], nctime.units)
            nctime_1980 = np.round(netCDF4.date2num(ncdatetime, time.units))
            tmp_data[ttmp:ttmp+ncvar.shape[0],:,:] = ncvar[:,:,:]
            tmp_time[ttmp:ttmp+ncvar.shape[0]] = nctime_1980[:]
            ttmp += ncvar.shape[0]
            nc.close()
        var1[t:t+tmp_data.shape[0],:,:] = tmp_data[:,:,:]
        time[t:t+tmp_data.shape[0]] = tmp_time[:]
        t += tmp_data.shape[0]

    # Methods to fill time and time_vectors variables:
    datetimes = netCDF4.num2date(time[:], time.units, time.calendar)
    time_vectors[:,:] = _datetimes_to_time_vectors(datetimes)

    nc1.close()
