import os
import glob
import datetime
from calendar import monthrange
import tempfile
import shutil

import numpy as np
import numpy.ma as ma
import netCDF4

from merra2_variables import subdaily_vars

# Aliases for default fill values
defi2 = netCDF4.default_fillvals['i2']
defi4 = netCDF4.default_fillvals['i4']
deff4 = netCDF4.default_fillvals['f4']


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


def fixed_netcdf(path_data, output_file, var_name, merra2_var_dict=None):
    """MERRA2 invariant NetCDF.

    Parameters
    ----------
    path_data : str
    output_file : str
    var_name : str
    merra2_var_dict : dict
        Dictionary containing the following keys:
        esdt_dir, collection, merra_name, standard_name,
        see the Bosilovich paper for details.
    verbose : bool

    """

    if not merra2_var_dict:
        merra2_var_dict = subdaily_vars[var_name]

    search_str = "*{0}*.nc4".format(merra2_var_dict['collection'])
    nc_files = glob.glob(os.path.join(path_data, search_str))

    nc_reference = netCDF4.Dataset(nc_files[0], 'r')
    var_ref = nc_reference.variables[merra2_var_dict['merra_name']]

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
    nc1.history = "%s: Extract variable." % (now,)
    nc1.institution = nc_reference.Institution
    nc1.source = 'Reanalysis'
    nc1.references = nc_reference.References
    for attr in nc_reference.ncattrs():
        setattr(nc1, 'original_file_' + attr, getattr(nc_reference, attr))

    # Create netCDF dimensions
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

    var1 = nc1.createVariable(var_name, 'f4', ('lat', 'lon'), zlib=True,
                              fill_value=deff4)
    # 3.1. Units
    var1.units = var_ref.units
    # 3.2. Long Name
    var1.long_name = var_ref.long_name
    # 3.3. Standard Name
    var1.standard_name = merra2_var_dict['standard_name']
    var_ref = nc_reference.variables[merra2_var_dict['merra_name']]
    var1[:,:] = var_ref[0,:,:]

    nc_reference.close()
    nc1.close()


def subdaily_download(dataset_esdt, merra2_collection, initial_year,
                      final_year, initial_month=1, final_month=12,
                      initial_day=1, final_day=None, output_directory=None):
    """MERRA2 subdaily download.

    Parameters
    ----------
    dataset_esdt : string
        See the Bosilovich paper for details.
    merra2_collection : string
        See the Bosilovich paper for details.
    initial_year : int
    final_year : int
    initial_month : int
    final_month : int
    initial_day : int
    final_day : int
    output_directory : string

    """

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
                    final_year, merra2_var_dict=None, verbose=False):
    """MERRA2 subdaily NetCDF.

    Parameters
    ----------
    path_data : str
    output_file : str
    var_name : str
    initial_year : int
    final_year : int
    merra2_var_dict : dict
        Dictionary containing the following keys:
        esdt_dir, collection, merra_name, standard_name,
        see the Bosilovich paper for details.
    verbose : bool

    Notes
    -----
    Currently does not deal with the level aspect (and therefore 3d fields).

    """

    if not merra2_var_dict:
        merra2_var_dict = subdaily_vars[var_name]

    search_str = "*{0}*.nc4".format(merra2_var_dict['collection'])
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
            ncvar = nc.variables[merra2_var_dict['merra_name']]
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
    var_ref = nc_reference.variables[merra2_var_dict['merra_name']]

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
    nc1.history = "%s: Extract variable and merge in time." % (now,)
    nc1.institution = nc_reference.Institution
    nc1.source = 'Reanalysis'
    nc1.references = nc_reference.References
    for attr in nc_reference.ncattrs():
        setattr(nc1, 'first_file_' + attr, getattr(nc_reference, attr))

    # Create netCDF dimensions
    nc1.createDimension('time', nt)
    # nc1.createDimension('ts', 6)
    # nc1.createDimension('level', k)
    nc1.createDimension('lat', len(nc_reference.dimensions['lat']))
    nc1.createDimension('lon', len(nc_reference.dimensions['lon']))
    if merra2_var_dict['cell_methods']:
        nc1.createDimension('nv', 2)

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
    # nc1.createVariable('ts', 'i2', ('ts',), zlib=True, fill_value=defi2)

    # 4.4. Time Coordinate
    if merra2_var_dict['cell_methods']:
        time = nc1.createVariable('time', 'f4', ('time',))
    else:
        time = nc1.createVariable('time', 'i4', ('time',))
    time.axis = 'T'
    time.units = "hours since 1980-01-01 00:00:00"
    time.long_name = 'time'
    time.standard_name = 'time'
    # 4.4.1. Calendar
    time.calendar = 'gregorian'

    if merra2_var_dict['cell_methods']:
        tbounds = nc1.createVariable('time_bounds', 'f4', ('time', 'nv'))

    # time_vectors = nc1.createVariable('time_vectors', 'i2', ('time', 'ts'),
    #                                   zlib=True)

    # 4.3. Vertical (Height or Depth) Coordinate
    # level = nc1.createVariable('level','f4',('level',),zlib=True)
    # level.axis = 'Z'
    # level.units = 'Pa'
    # level.positive = 'up'
    # level.long_name = 'air_pressure'
    # level.standard_name = 'air_pressure'

    # 4.1. Latitude Coordinate
    lat = nc1.createVariable('lat', 'f4', ('lat',))
    lat.axis = 'Y'
    lat.units = 'degrees_north'
    lat.long_name = 'latitude'
    lat.standard_name = 'latitude'
    lat[:] = nc_reference.variables['lat'][:]

    # 4.2. Longitude Coordinate
    lon = nc1.createVariable('lon', 'f4', ('lon',))
    lon.axis = 'X'
    lon.units = 'degrees_east'
    lon.long_name = 'longitude'
    lon.standard_name = 'longitude'
    lon[:] = nc_reference.variables['lon'][:]

    var1 = nc1.createVariable(var_name, 'f4', ('time', 'lat', 'lon'),
                              zlib=True, chunksizes=(360, 30, 30),
                              fill_value=deff4)
    # 3.1. Units
    # Force kg kg-1 to 1
    if var_ref.units == 'kg kg-1':
        var1.units = '1'
    else:
        var1.units = var_ref.units
    # 3.2. Long Name
    var1.long_name = var_ref.long_name
    # 3.3. Standard Name
    var1.standard_name = merra2_var_dict['standard_name']
    if merra2_var_dict['cell_methods']:
        var1.cell_methods = merra2_var_dict['cell_methods']

    nc_reference.close()

    t = 0
    for i, division in enumerate(divided_files):
        tmp_data = ma.masked_all([nt_division[i], len(nc1.dimensions['lat']),
                                  len(nc1.dimensions['lon'])])
        tmp_time = ma.masked_all([nt_division[i]])
        ttmp = 0
        for nc_file in division:
            if verbose:
                print(nc_file)
            nc = netCDF4.Dataset(nc_file, 'r')
            ncvar = nc.variables[subdaily_vars[var_name]['merra_name']]
            nctime = nc.variables['time']
            ncdatetime = netCDF4.num2date(nctime[:], nctime.units)
            if merra2_var_dict['cell_methods']:
                nctime_1980 = netCDF4.date2num(ncdatetime, time.units)
            else:
                nctime_1980 = np.round(
                    netCDF4.date2num(ncdatetime, time.units))
            tmp_data[ttmp:ttmp+ncvar.shape[0],:,:] = ncvar[:,:,:]
            tmp_time[ttmp:ttmp+ncvar.shape[0]] = nctime_1980[:]
            ttmp += ncvar.shape[0]
            nc.close()
        var1[t:t+tmp_data.shape[0],:,:] = tmp_data[:,:,:]
        time[t:t+tmp_data.shape[0]] = tmp_time[:]
        if merra2_var_dict['cell_methods']:
            if tmp_time[1]-tmp_time[0] == 1.0:
                tbounds[t:t+tmp_data.shape[0],0] = tmp_time[:]-0.5
                tbounds[t:t+tmp_data.shape[0],1] = tmp_time[:]+0.5
            elif tmp_time[1]-tmp_time[0] == 3.0:
                tbounds[t:t+tmp_data.shape[0],0] = tmp_time[:]-1.5
                tbounds[t:t+tmp_data.shape[0],1] = tmp_time[:]+1.5
        t += tmp_data.shape[0]

    # Methods to fill time and time_vectors variables:
    # datetimes = netCDF4.num2date(time[:], time.units, time.calendar)
    # time_vectors[:,:] = _datetimes_to_time_vectors(datetimes)

    nc1.close()


def subdaily_download_and_convert(var_names, initial_year, final_year,
                                  initial_month=1, final_month=12,
                                  initial_day=1, final_day=None,
                                  merra2_var_dicts=None, output_dir=None,
                                  delete_temp_dir=True, verbose=True):
    """MERRA2 subdaily download and conversion.

    Parameters
    ----------
    var_names : list of string
        Variable short names, must be defined in merra2_variables.py
        if merra2_var_dict is not provided. If more than one variable,
        they are assumed to have the same original files and those will only
        be downloaded once.
    initial_year : int
    final_year : int
    initial_month : int
    final_month : int
    initial_day : int
    final_day : int
    merra2_var_dicts : list of dict
        Dictionary containing the following keys:
        esdt_dir, collection, merra_name, standard_name,
        see the Bosilovich paper for details. Same order as var_names.
    output_dir : string
    delete_temp_dir : bool
    verbose : bool

    """

    if output_dir is None:
        output_dir = os.getcwd()
    temp_dir_download = tempfile.mkdtemp(dir=output_dir)
    for i, var_name in enumerate(var_names):
        if not merra2_var_dicts:
            merra2_var_dict = subdaily_vars[var_name]
        else:
            merra2_var_dict = merra2_var_dicts[i]
        # Download subdaily files
        if i == 0:
            subdaily_download(
                merra2_var_dict['esdt_dir'], merra2_var_dict['collection'],
                initial_year, final_year, initial_month=initial_month,
                final_month=final_month, initial_day=initial_day,
                final_day=final_day, output_directory=temp_dir_download)
        # Name the output file
        if (initial_year == final_year) and (initial_month == final_month):
            file_name_str = "{0}_1hr_merra2_reanalysis_{1}{2}.nc"
            out_file_name = file_name_str.format(
                var_name, str(initial_year), str(initial_month).zfill(2))
        else:
            file_name_str = "{0}_1hr_merra2_reanalysis_{1}{2}-{3}{4}.nc"
            out_file_name = file_name_str.format(
                var_name, str(initial_year), str(initial_month).zfill(2),
                str(final_year), str(final_month).zfill(2))
        out_file = os.path.join(output_dir, out_file_name)
        # Extract variable
        subdaily_netcdf(
            temp_dir_download, out_file, var_name, initial_year, final_year,
            verbose=verbose)
    if delete_temp_dir:
        shutil.rmtree(temp_dir_download)


def daily_netcdf(path_data, output_file, var_name, initial_year, final_year,
                 merra2_var_dict=None, verbose=False):
    """MERRA2 daily NetCDF.

    Parameters
    ----------
    path_data : str
    output_file : str
    var_name : str
    initial_year : int
    final_year : int
    merra2_var_dict : dict
        Dictionary containing the following keys:
        esdt_dir, collection, merra_name, standard_name,
        see the Bosilovich paper for details.
    verbose : bool

    Notes
    -----
    Currently does not deal with the level aspect (and therefore 3d fields).
    This is almost identical to subdaily_netcdf but with a different chunksize,
    perhaps not an ideal separation as this is duplicated code...

    """

    if not merra2_var_dict:
        merra2_var_dict = subdaily_vars[var_name]

    search_str = "*{0}*.nc4".format(merra2_var_dict['collection'])
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
            ncvar = nc.variables[merra2_var_dict['merra_name']]
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
    var_ref = nc_reference.variables[merra2_var_dict['merra_name']]

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
    nc1.history = "%s: Extract variable and merge in time." % (now,)
    nc1.institution = nc_reference.Institution
    nc1.source = 'Reanalysis'
    nc1.references = nc_reference.References
    for attr in nc_reference.ncattrs():
        setattr(nc1, 'first_file_' + attr, getattr(nc_reference, attr))

    # Create netCDF dimensions
    nc1.createDimension('time', nt)
    # nc1.createDimension('ts', 6)
    # nc1.createDimension('level', k)
    nc1.createDimension('lat', len(nc_reference.dimensions['lat']))
    nc1.createDimension('lon', len(nc_reference.dimensions['lon']))
    if merra2_var_dict['cell_methods']:
        nc1.createDimension('nv', 2)

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
    # nc1.createVariable('ts', 'i2', ('ts',), zlib=True, fill_value=defi2)

    # 4.4. Time Coordinate
    time = nc1.createVariable('time', 'i4', ('time',), zlib=True)
    time.axis = 'T'
    time.units = "hours since 1980-01-01 00:00:00"
    time.long_name = 'time'
    time.standard_name = 'time'
    # 4.4.1. Calendar
    time.calendar = 'gregorian'

    if merra2_var_dict['cell_methods']:
        tbounds = nc1.createVariable('time_bounds', 'i4', ('time', 'nv'))

    #time_vectors = nc1.createVariable('time_vectors', 'i2', ('time', 'ts'),
    #                                  zlib=True)

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
                              zlib=True, chunksizes=(30, 120, 120),
                              fill_value=deff4)
    # 3.1. Units
    var1.units = var_ref.units
    # 3.2. Long Name
    var1.long_name = var_ref.long_name
    # 3.3. Standard Name
    var1.standard_name = merra2_var_dict['standard_name']
    if merra2_var_dict['cell_methods']:
        var1.cell_methods = merra2_var_dict['cell_methods']

    nc_reference.close()

    t = 0
    for i, division in enumerate(divided_files):
        tmp_data = ma.masked_all([nt_division[i], len(nc1.dimensions['lat']),
                                  len(nc1.dimensions['lon'])])
        tmp_time = ma.masked_all([nt_division[i]])
        ttmp = 0
        for nc_file in division:
            if verbose:
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
        if merra2_var_dict['cell_methods']:
            tbounds[t:t+tmp_data.shape[0],0] = tmp_time[:]-12
            tbounds[t:t+tmp_data.shape[0],1] = tmp_time[:]+12
        t += tmp_data.shape[0]

    # Methods to fill time and time_vectors variables:
    #datetimes = netCDF4.num2date(time[:], time.units, time.calendar)
    #time_vectors[:,:] = _datetimes_to_time_vectors(datetimes)

    nc1.close()


def daily_download_and_convert(var_names, initial_year, final_year,
                               initial_month=1, final_month=12,
                               initial_day=1, final_day=None,
                               merra2_var_dicts=None, output_dir=None,
                               delete_temp_dir=True, verbose=True):
    """MERRA2 daily download and conversion.

    Parameters
    ----------
    var_names : list of string
        Variable short names, must be defined in merra2_variables.py
        if merra2_var_dict is not provided. If more than one variable,
        they are assumed to have the same original files and those will only
        be downloaded once.
    initial_year : int
    final_year : int
    initial_month : int
    final_month : int
    initial_day : int
    final_day : int
    merra2_var_dicts : list of dict
        Dictionary containing the following keys:
        esdt_dir, collection, merra_name, standard_name,
        see the Bosilovich paper for details. Same order as var_names.
    output_dir : string
    delete_temp_dir : bool
    verbose : bool

    Notes
    -----
    This is almost identical to subdaily_download_and_convert,
    perhaps not an ideal separation as this is duplicated code...

    """

    if output_dir is None:
        output_dir = os.getcwd()
    temp_dir_download = tempfile.mkdtemp(dir=output_dir)
    for i, var_name in enumerate(var_names):
        if not merra2_var_dicts:
            merra2_var_dict = subdaily_vars[var_name]
        else:
            merra2_var_dict = merra2_var_dicts[i]
        # Download subdaily files
        if i == 0:
            subdaily_download(
                merra2_var_dict['esdt_dir'], merra2_var_dict['collection'],
                initial_year, final_year, initial_month=initial_month,
                final_month=final_month, initial_day=initial_day,
                final_day=final_day, output_directory=temp_dir_download)
        # Name the output file
        if (initial_year == final_year):
            file_name_str = "{0}_day_merra2_reanalysis_{1}.nc"
            out_file_name = file_name_str.format(var_name, str(initial_year))
        else:
            file_name_str = "{0}_day_merra2_reanalysis_{1}-{2}.nc"
            out_file_name = file_name_str.format(
                var_name, str(initial_year), str(final_year))
        out_file = os.path.join(output_dir, out_file_name)
        # Extract variable
        daily_netcdf(
            temp_dir_download, out_file, var_name, initial_year, final_year,
            verbose=verbose)
    if delete_temp_dir:
        shutil.rmtree(temp_dir_download)
