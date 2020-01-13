import datetime
import shutil
import subprocess
import sys
import glob
import tempfile
import logging
from calendar import monthrange
from pathlib import Path
from typing import List
from typing import Optional
from typing import Union

import netCDF4
import numpy as np
import numpy.ma as ma

from pymerra2 import __version__
from pymerra2.variables import var_list

# Aliases for default fill values
defi2 = netCDF4.default_fillvals["i2"]
defi4 = netCDF4.default_fillvals["i4"]
deff4 = netCDF4.default_fillvals["f4"]

# Aliases for different sizes
KiB = 2 ** 10
MiB = 2 ** 20
GiB = 2 ** 30


class NetCDFError(Exception):
    pass


# TODO: Integrate more opportunities for tests and simplify function calls
# TODO: Create a class for data requests that can be passed to functions rather than script calls.
# TODO: Refactor keywords and variables to remove ambiguity/confusion
# TODO: Split functions between downloads and subsets
# TODO: Integrate data downloads methods for similar data access types (Not just MERRA/MERRA2)
# TODO: Write docstrings for ALL methods and classes from now on


def _time_vectors_int(
    time_vectors, force=False, raise_exception=False, allow_masked=True, dtype="int32"
):
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


def _datetimes_to_time_vectors(
    datetimes: Union[datetime.datetime, List[datetime.datetime]]
):
    """Convert list of datetimes to Nx6 matrix.

    Parameters
    ----------
    datetimes : Union[datetime.datetime, List[datetime.datetime]]

    Returns
    -------
    out : numpy.array
        Nx6 matrix of time vectors

    """

    # Does not support microsecond (datetime shape of length 7)
    def datetime_timetuple(one_datetime):
        if one_datetime is None:
            return ma.masked_all([6], dtype="int32")
        return one_datetime.timetuple()[0:6]

    try:
        time_tuples = map(datetime_timetuple, datetimes)
        return _time_vectors_int(ma.array(time_tuples))
    except (AttributeError, TypeError):
        time_tuples = datetimes.timetuple()
        return _time_vectors_int(ma.array(time_tuples))


def get_nc_attr(nc, name, default=None):
    """Non-error raising netCDF attribute getter"""
    try:
        return nc.getncattr(name)
    except AttributeError:
        return default


def fixed_netcdf(
    path_data: Union[str, Path],
    output_file: Union[str, Path],
    var_name: str,
    merra2_var_dict: Optional[dict] = None,
):
    """MERRA2 invariant NetCDF.

    Parameters
    ----------
    path_data : Union[str, Path]
    output_file : Union[str, Path]
    var_name : str
    merra2_var_dict : Optional[dict]
        Dictionary containing the following keys:
        esdt_dir, collection, merra_name, standard_name,
        see the Bosilovich paper for details.

    """
    if not isinstance(path_data, Path):
        path_data = Path(path_data)

    if not merra2_var_dict:
        merra2_var_dict = var_list[var_name]

    search_str = "*{0}*.nc4".format(merra2_var_dict["collection"])
    nc_files = [f for f in path_data.rglob(search_str)]

    nc_reference = netCDF4.Dataset(nc_files[0], "r")
    var_ref = nc_reference.variables[merra2_var_dict["merra_name"]]

    # 2.1 Filename
    #     NetCDF files should have the file name extension ".nc".
    nc_file = output_file

    now = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    nc1 = netCDF4.Dataset(nc_file, "w", format="NETCDF4_CLASSIC")

    # 2.6.1 Identification of Conventions
    nc1.Conventions = "CF-1.7"

    # 2.6.2. Description of file contents
    nc1.title = (
        "Modern-Era Retrospective analysis for Research and " "Applications, Version 2"
    )
    nc1.history = (
        "{0}\n{1} (pymerra2-{2}): " "Reformat to CF-1.7 & " "Extract variable."
    ).format(nc_reference.History, now, __version__)
    nc1.institution = nc_reference.Institution
    nc1.source = "Reanalysis"
    nc1.references = nc_reference.References
    # Using lower case c for conventions because lower() is used below...
    attr_overwrite = ["conventions", "title", "institution", "source", "references"]
    ordered_attr = {}
    for attr in nc_reference.ncattrs():
        if attr == "History":
            continue
        if attr.lower() in attr_overwrite:
            ordered_attr["original_file_{}".format(attr)] = getattr(nc_reference, attr)
        else:
            ordered_attr[attr] = getattr(nc_reference, attr)
    for attr in sorted(ordered_attr.keys(), key=lambda v: v.lower()):
        setattr(nc1, attr, ordered_attr[attr])

    # Create netCDF dimensions
    nc1.createDimension("lat", len(nc_reference.dimensions["lat"]))
    nc1.createDimension("lon", len(nc_reference.dimensions["lon"]))

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
    lat = nc1.createVariable("lat", "f4", ("lat",), zlib=True)
    lat.axis = "Y"
    lat.units = "degrees_north"
    lat.long_name = "latitude"
    lat.standard_name = "latitude"
    lat[:] = nc_reference.variables["lat"][:]

    # 4.2. Longitude Coordinate
    lon = nc1.createVariable("lon", "f4", ("lon",), zlib=True)
    lon.axis = "X"
    lon.units = "degrees_east"
    lon.long_name = "longitude"
    lon.standard_name = "longitude"
    lon[:] = nc_reference.variables["lon"][:]

    least_digit = merra2_var_dict.get("least_significant_digit", None)
    var1 = nc1.createVariable(
        var_name,
        "f4",
        ("lat", "lon"),
        zlib=True,
        fill_value=deff4,
        least_significant_digit=least_digit,
    )
    # 3.1. Units
    var1.units = var_ref.units
    # 3.2. Long Name
    var1.long_name = var_ref.long_name
    # 3.3. Standard Name
    var1.standard_name = merra2_var_dict["standard_name"]
    var_ref = nc_reference.variables[merra2_var_dict["merra_name"]]
    var1[:, :] = var_ref[0, :, :]

    nc_reference.close()
    nc1.close()


def subdaily_download(
    merra2_server: str,
    dataset_esdt: str,
    merra2_collection: str,
    initial_year: int,
    final_year: int,
    initial_month: int = 1,
    final_month: int = 12,
    initial_day: int = 1,
    final_day: Optional[int] = None,
    output_directory: Union[str, Path] = None,
):
    """MERRA2 subdaily download.

    Parameters
    ----------
    merra2_server : str
        Must contain trailing slash.
        e.g. https://goldsmr4.gesdisc.eosdis.nasa.gov/data/
    dataset_esdt : str
        See the Bosilovich paper for details.
    merra2_collection : str
        See the Bosilovich paper for details.
    initial_year : int
    final_year : int
    initial_month : int
    final_month : int
    initial_day : int
    final_day : Optional[int]
    output_directory : Union[str, Path]
    """

    if output_directory is None:
        add_output_dir = ""
    else:
        add_output_dir = "--directory-prefix={0}".format(output_directory)

    data_path = "MERRA2/{4}/{0}/{1}/" "MERRA2_{3}.{5}.{0}{1}{2}.nc4"
    for yyyy in range(initial_year, final_year + 1):
        if yyyy < 1992:
            merra_stream = "100"
        elif yyyy < 2001:
            merra_stream = "200"
        elif yyyy < 2011:
            merra_stream = "300"
        else:
            merra_stream = "400"
        if yyyy == initial_year:
            mi = initial_month
        else:
            mi = 1
        if yyyy == final_year:
            mf = final_month
        else:
            mf = 12
        for mm in range(mi, mf + 1):
            if (yyyy == initial_year) and (mm == mi):
                di = initial_day
            else:
                di = 1
            if final_day and (yyyy == final_year) and (mm == mf):
                df = final_day
            else:
                mrange = monthrange(yyyy, mm)
                df = mrange[1]
            for dd in range(di, df + 1):
                cdp = data_path.format(
                    str(yyyy),
                    str(mm).zfill(2),
                    str(dd).zfill(2),
                    merra_stream,
                    dataset_esdt,
                    merra2_collection,
                )
                tries = 0
                while True:
                    try:
                        subprocess.check_call(
                            [
                                "wget",
                                "-c",
                                add_output_dir,
                                "--load-cookies",
                                str(Path("~/.urs_cookies").expanduser()),
                                "--save-cookies",
                                str(Path("~/.urs_cookies").expanduser()),
                                "--keep-session-cookies",
                                merra2_server + cdp,
                            ]
                        )
                        logging.info("File `{}` successfully downloaded.".format(cdp))
                        break
                    except subprocess.CalledProcessError:
                        tries += 1
                        if tries > 3:
                            logging.exception("Failed to download file: {}".format(cdp))
                            raise
                        continue


def subdaily_netcdf(
    path_data: Union[str, Path],
    output_file: Union[str, Path],
    var_name: str,
    initial_year: int,
    final_year: int,
    initial_month: int = 1,
    final_month: int = 12,
    merra2_var_dict: Optional[dict] = None,
    verbose: bool = False,
):
    """MERRA2 subdaily NetCDF.

    Parameters
    ----------
    path_data : Union[str, Path]
    output_file : Union[str, Path]
    var_name : str
    initial_year : int
    final_year : int
    initial_month : int
    final_month : int
    merra2_var_dict : Optional[dict]
        Dictionary containing the following keys:
        esdt_dir, collection, merra_name, standard_name,
        see the Bosilovich paper for details.
    verbose : bool

    """
    if not isinstance(path_data, Path):
        path_data = Path(path_data)

    if not merra2_var_dict:
        merra2_var_dict = var_list[var_name]

    search_str = "*{0}*.nc4".format(merra2_var_dict["collection"])
    nc_files = [str(f) for f in path_data.rglob(search_str)]
    nc_files.sort()

    relevant_files = []
    divided_files = [[]]
    nt_division = [0]
    nt = 0
    nmb = 0
    for nc_file in nc_files:
        yyyy = int(nc_file.split(".")[-2][0:4])
        mm = int(nc_file.split(".")[-2][4:6])
        if (
            (yyyy >= initial_year)
            and (yyyy <= final_year)
            and (mm >= initial_month)
            and (mm <= final_month)
        ):
            relevant_files.append(nc_file)
            nc = netCDF4.Dataset(nc_file, "r")
            ncvar = nc.variables[merra2_var_dict["merra_name"]]
            nmb += (ncvar.size * 4) / MiB
            if nmb > 512 * 2:
                divided_files.append([nc_file])
                nt_division.append(0)
                nmb = (ncvar.size * 4) / MiB
            else:
                divided_files[-1].append(nc_file)
            nt_division[-1] += len(nc.dimensions["time"])
            nt += len(nc.dimensions["time"])
            nc.close()

    nc_reference = netCDF4.Dataset(relevant_files[0], "r")
    var_ref = nc_reference.variables[merra2_var_dict["merra_name"]]

    # 2.1 Filename
    #     NetCDF files should have the file name extension ".nc".
    nc_file = output_file

    now = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    nc1 = netCDF4.Dataset(nc_file, "w", format="NETCDF4_CLASSIC")

    # 2.6.1 Identification of Conventions
    nc1.Conventions = "CF-1.7"

    # 2.6.2. Description of file contents
    nc1.title = (
        "Modern-Era Retrospective analysis for Research and " "Applications, Version 2"
    )
    if (len(divided_files) == 1) and (len(divided_files[0]) == 1):
        nc1.history = (
            "{0}\n{1} (pymerra2-{2}): " "Reformat to CF-1.7 & " "Extract variable."
        ).format(nc_reference.History, now, __version__)
    else:
        nc1.history = (
            "{0}\n{1} (pymerra2-{2}): "
            "Reformat to CF-1.7 & "
            "Extract variable & "
            "Merge in time."
        ).format(nc_reference.History, now, __version__)
    nc1.institution = nc_reference.Institution
    nc1.source = "Reanalysis"
    nc1.references = nc_reference.References

    # Using lower case c for conventions because lower() is used below...
    attr_overwrite = ["conventions", "title", "institution", "source", "references"]
    ordered_attr = {}
    for attr in nc_reference.ncattrs():
        if attr == "History":
            continue
        if attr.lower() in attr_overwrite:
            ordered_attr["original_file_" + attr] = getattr(nc_reference, attr)
        else:
            ordered_attr[attr] = getattr(nc_reference, attr)
    for attr in sorted(ordered_attr.keys(), key=lambda v: v.lower()):
        setattr(nc1, attr, ordered_attr[attr])

    # Create netCDF dimensions
    nc1.createDimension("time", nt)
    # nc1.createDimension('ts', 6)
    if "lev" in nc_reference.dimensions:
        nc1.createDimension("level", nc_reference.dimensions["lev"].size)
    nc1.createDimension("lat", len(nc_reference.dimensions["lat"]))
    nc1.createDimension("lon", len(nc_reference.dimensions["lon"]))
    if merra2_var_dict["cell_methods"]:
        nc1.createDimension("nv", 2)

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
    if merra2_var_dict["cell_methods"]:
        time = nc1.createVariable("time", "f4", ("time",))
    else:
        time = nc1.createVariable("time", "i4", ("time",))
    time.axis = "T"
    time.units = "hours since 1980-01-01 00:00:00"
    time.long_name = "time"
    time.standard_name = "time"
    # 4.4.1. Calendar
    time.calendar = "gregorian"

    tbounds = None
    if merra2_var_dict["cell_methods"]:
        # WARNING: This breaks time. Caveat emptor.
        # time.bounds = "time_bnds"
        tbounds = nc1.createVariable("time_bnds", "f4", ("time", "nv"))

    # time_vectors = nc1.createVariable('time_vectors', 'i2', ('time', 'ts'),
    #                                   zlib=True)

    # 4.3. Vertical (Height or Depth) Coordinate
    if "lev" in nc_reference.dimensions:
        level = nc1.createVariable("level", "i4", ("level",), zlib=True)
        level.axis = "Z"
        level.units = "Pa"
        level.positive = "down"
        level.long_name = "air_pressure"
        level.standard_name = "air_pressure"
        level_ref = nc_reference.variables["lev"]
        if level_ref.units == "hPa":
            level[:] = np.round(level_ref[:] * 100)
        else:
            raise NotImplementedError()

    # 4.1. Latitude Coordinate
    lat = nc1.createVariable("lat", "f4", ("lat",))
    lat.axis = "Y"
    lat.units = "degrees_north"
    lat.long_name = "latitude"
    lat.standard_name = "latitude"
    lat[:] = nc_reference.variables["lat"][:]

    # 4.2. Longitude Coordinate
    lon = nc1.createVariable("lon", "f4", ("lon",))
    lon.axis = "X"
    lon.units = "degrees_east"
    lon.long_name = "longitude"
    lon.standard_name = "longitude"
    lon[:] = nc_reference.variables["lon"][:]

    if "lev" in nc_reference.dimensions:
        var_dims = ("time", "level", "lat", "lon")
    else:
        var_dims = ("time", "lat", "lon")
    # Here, the chunksize is set rather arbitrarily...
    if "lev" in nc_reference.dimensions:
        var1_chunksizes = (8, 6, 121, 144)
    else:
        var1_chunksizes = (360, 30, 30)
    least_digit = merra2_var_dict.get("least_significant_digit", None)
    var1 = nc1.createVariable(
        var_name,
        "f4",
        var_dims,
        zlib=True,
        chunksizes=var1_chunksizes,
        fill_value=deff4,
        least_significant_digit=least_digit,
    )
    # 3.1. Units
    var1.units = var_ref.units
    # 3.2. Long Name
    var1.long_name = var_ref.long_name
    # 3.3. Standard Name
    var1.standard_name = merra2_var_dict["standard_name"]
    if merra2_var_dict["cell_methods"]:
        var1.cell_methods = merra2_var_dict["cell_methods"]

    nc_reference.close()

    t = 0
    for i, division in enumerate(divided_files):
        if "lev" in nc_reference.dimensions:
            tmp_data = ma.masked_all(
                [
                    nt_division[i],
                    len(nc1.dimensions["level"]),
                    len(nc1.dimensions["lat"]),
                    len(nc1.dimensions["lon"]),
                ]
            )
        else:
            tmp_data = ma.masked_all(
                [nt_division[i], len(nc1.dimensions["lat"]), len(nc1.dimensions["lon"])]
            )
        tmp_time = ma.masked_all([nt_division[i]])
        ttmp = 0
        for nc_file in division:
            if verbose:
                print(nc_file)
            nc = netCDF4.Dataset(nc_file, "r")
            ncvar = nc.variables[merra2_var_dict["merra_name"]]
            nctime = nc.variables["time"]

            ncdatetime = netCDF4.num2date(nctime[:], nctime.units)

            if merra2_var_dict["cell_methods"]:
                nctime_1980 = netCDF4.date2num(ncdatetime, time.units)
            else:
                nctime_1980 = np.round(netCDF4.date2num(ncdatetime, time.units))
            if "lev" in nc_reference.dimensions:
                tmp_data[ttmp : ttmp + ncvar.shape[0], :, :, :] = ncvar[:, :, :, :]
            else:
                tmp_data[ttmp : ttmp + ncvar.shape[0], :, :] = ncvar[:, :]
            tmp_time[ttmp : ttmp + ncvar.shape[0]] = nctime_1980[:]
            ttmp += ncvar.shape[0]
            nc.close()
        if "lev" in nc_reference.dimensions:
            var1[t : t + tmp_data.shape[0], :, :, :] = tmp_data[:, :, :, :]
        else:
            var1[t : t + tmp_data.shape[0], :, :] = tmp_data[:, :, :]
        time[t : t + tmp_data.shape[0]] = tmp_time[:]
        if merra2_var_dict["cell_methods"]:
            # TODO: There is something not working here. 29 Oct 2019.
            if tmp_time[1] - tmp_time[0] == 1.0:
                tbounds[t : t + tmp_data.shape[0], 0] = tmp_time[:] - 0.5
                tbounds[t : t + tmp_data.shape[0], 1] = tmp_time[:] + 0.5
            elif tmp_time[1] - tmp_time[0] == 3.0:
                tbounds[t : t + tmp_data.shape[0], 0] = tmp_time[:] - 1.5
                tbounds[t : t + tmp_data.shape[0], 1] = tmp_time[:] + 1.5
        t += tmp_data.shape[0]

    nc1.close()


def file_namer(
    var_name,
    time_frequency: str,
    initial_year,
    final_year,
    initial_month,
    final_month,
    initial_day: int = None,
    final_day: int = None,
):
    # Name the output file
    if (initial_year == final_year) and (initial_month == final_month):
        if initial_day == final_day and (None not in [initial_day, final_day]):
            file_name_str = "{0}_{1}_merra2_reanalysis_{2}{3}{4}.nc"
            out_file_name = file_name_str.format(
                var_name,
                time_frequency,
                str(initial_year),
                str(initial_month).zfill(2),
                str(initial_day).zfill(2),
            )
        else:
            file_name_str = "{0}_{1}_merra2_reanalysis_{2}{3}.nc"
            out_file_name = file_name_str.format(
                var_name, time_frequency, str(initial_year), str(initial_month).zfill(2)
            )
    else:
        file_name_str = "{0}_{1}_merra2_reanalysis_{2}{3}-{4}{5}.nc"
        out_file_name = file_name_str.format(
            var_name,
            time_frequency,
            str(initial_year),
            str(initial_month).zfill(2),
            str(final_year),
            str(final_month).zfill(2),
        )

    return out_file_name


def subdaily_download_and_convert(
    merra2_server: str,
    var_names: Union[str, List[str]],
    initial_year: int,
    final_year: int,
    initial_month: int = 1,
    final_month: int = 12,
    initial_day: int = 1,
    final_day: Optional[int] = None,
    merra2_var_dicts: Optional[List[dict]] = None,
    output_dir: Union[str, Path] = None,
    delete_temp_dir: bool = False,
    verbose: bool = False,
    time_frequency: str = "1hr",
):
    """MERRA2 subdaily download and conversion.

    Parameters
    ----------
    merra2_server : str
        Must contain trailing slash.
        e.g. https://goldsmr4.gesdisc.eosdis.nasa.gov/data/
    var_names : Union[str, List[str]]
        Variable short names, must be defined in pymerra2_variables.py
        if merra2_var_dict is not provided. If more than one variable,
        they are assumed to have the same original files and those will only
        be downloaded once.
    initial_year : int
    final_year : int
    initial_month : int
    final_month : int
    initial_day : int
    final_day : Optional[int]
    merra2_var_dicts : Optional[List[dict]]
        Dictionary containing the following keys:
        esdt_dir, collection, merra_name, standard_name,
        see the Bosilovich paper for details. Same order as var_names.
    output_dir : Union[str, Path]
    delete_temp_dir : bool
    verbose : bool
    time_frequency: str
    """
    if isinstance(output_dir, Path):
        output_dir = Path(output_dir)
    if output_dir is None:
        output_dir = Path.cwd()
    output_dir.mkdir(exist_ok=True)

    if (2, 7) < sys.version_info < (3, 6):
        output_dir = str(output_dir)

    if isinstance(var_names, str):
        var_names = [var_names]

    temp_dir_download = tempfile.mkdtemp(dir=output_dir)
    for i, var_name in enumerate(var_names):
        if not merra2_var_dicts:
            merra2_var_dict = var_list[var_name]
        else:
            merra2_var_dict = merra2_var_dicts[i]
        # Download subdaily files
        subdaily_download(
            merra2_server,
            merra2_var_dict["esdt_dir"],
            merra2_var_dict["collection"],
            initial_year,
            final_year,
            initial_month=initial_month,
            final_month=final_month,
            initial_day=initial_day,
            final_day=final_day,
            output_directory=temp_dir_download,
        )

        # Name the output file
        out_file_name = file_namer(
            var_name,
            time_frequency,
            initial_year,
            final_year,
            initial_month,
            final_month,
            initial_day,
            final_day,
        )
        out_file = Path(output_dir).joinpath(out_file_name)

        # Extract variable
        subdaily_netcdf(
            temp_dir_download,
            out_file,
            var_name,
            initial_year,
            final_year,
            final_month,
            final_year,
            verbose=verbose,
        )

    if delete_temp_dir:
        shutil.rmtree(temp_dir_download)


def daily_netcdf(
    path_data: Union[str, Path, List[str]],
    output_file: Union[str, Path],
    var_name: str,
    initial_year: int,
    final_year: int,
    time_from_filename: bool = False,
    merra2_var_dict: Optional[dict] = None,
    verbose: bool = False,
):
    """MERRA2 daily NetCDF.

    Parameters
    ----------
    path_data : Union[str, Path, List[str]]
        If path_data is a list, skips the glob search.
    output_file : Union[str, Path]
    var_name : str
    initial_year : int
    final_year : int
    time_from_filename : bool
        Whether to take the timestamp from the file or not.
    merra2_var_dict : Optional[dict]
        Dictionary containing the following keys:
        esdt_dir, collection, merra_name, standard_name,
        see the Bosilovich paper for details.
        May also contain the following netCDF attributes:
        Title, Institution, Reference, Source
    verbose : bool

    Notes
    -----
    Currently does not deal with the level aspect (and therefore 3d fields).
    This is almost identical to subdaily_netcdf but with a different chunksize,
    perhaps not an ideal separation as this is duplicated code...

    """
    if isinstance(path_data, list):
        nc_files = path_data
    else:
        if not isinstance(path_data, Path):
            path_data = Path(path_data)

        if not merra2_var_dict:
            merra2_var_dict = var_list[var_name]

        search_str = "*{0}*.nc4".format(merra2_var_dict["collection"])
        nc_files = [str(f) for f in path_data.rglob(search_str)]
    nc_files.sort()

    relevant_files = []
    divided_files = [[]]
    nt_division = [0]
    nt = 0
    nmb = 0
    for nc_file in nc_files:
        yyyy = int(nc_file.split(".")[-2][0:4])
        if (yyyy >= initial_year) and (yyyy <= final_year):
            relevant_files.append(nc_file)
            nc = netCDF4.Dataset(nc_file, "r")
            ncvar = nc.variables[merra2_var_dict["merra_name"]]
            nmb += (ncvar.size * 4) / MiB
            if nmb > 512:
                divided_files.append([nc_file])
                nt_division.append(0)
                nmb = (ncvar.size * 4) / MiB
            else:
                divided_files[-1].append(nc_file)
            nt_division[-1] += len(nc.dimensions["time"])
            nt += len(nc.dimensions["time"])
            nc.close()

    nc_reference = netCDF4.Dataset(relevant_files[0], "r")
    var_ref = nc_reference.variables[merra2_var_dict["merra_name"]]

    # 2.1 Filename
    #     NetCDF files should have the file name extension ".nc".
    nc_file = output_file

    now = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    nc1 = netCDF4.Dataset(nc_file, "w", format="NETCDF4_CLASSIC")

    # 2.6.1 Identification of Conventions
    nc1.Conventions = "CF-1.7"

    # 2.6.2. Description of file contents
    nc1.title = merra2_var_dict.get(
        "Title",
        "Modern-Era Retrospective analysis for Research and Applications, Version 2",
    )
    if (len(divided_files) == 1) and (len(divided_files[0]) == 1):
        nc1.history = (
            "{0}\n{1} (pymerra2-{2}): " "Reformat to CF-1.7 & " "Extract variable."
        ).format(nc_reference.History, now, __version__)
    else:
        nc1.history = (
            "{0}\n{1} (pymerra2-{2}): "
            "Reformat to CF-1.7 & "
            "Extract variable & "
            "Merge in time."
        ).format(nc_reference.History, now, __version__)

    nc1.institution = get_nc_attr(
        nc_reference,
        "Institution",
        get_nc_attr(
            nc_reference,
            "Center",
            merra2_var_dict.get(
                "Institution", "NASA Global Modeling and Assimilation Office"
            ),
        ),
    )
    nc1.source = get_nc_attr(
        nc_reference, "Source", merra2_var_dict.get("Source", "Reanalysis")
    )
    nc1.references = get_nc_attr(
        nc_reference,
        "References",
        merra2_var_dict.get("References", "http://gmao.gsfc.nasa.gov"),
    )

    # Using lower case c for conventions because lower() is used below...
    attr_overwrite = ["conventions", "title", "institution", "source", "references"]
    ordered_attr = {}
    for attr in nc_reference.ncattrs():
        if attr == "History":
            continue
        if attr.lower() in attr_overwrite:
            ordered_attr["original_file_" + attr] = getattr(nc_reference, attr)
        else:
            ordered_attr[attr] = getattr(nc_reference, attr)
    for attr in sorted(ordered_attr.keys(), key=lambda v: v.lower()):
        setattr(nc1, attr, ordered_attr[attr])

    # Create netCDF dimensions
    nc1.createDimension("time", nt)
    # nc1.createDimension('ts', 6)
    # nc1.createDimension('level', k)
    nc1.createDimension("lat", len(nc_reference.dimensions["lat"]))
    nc1.createDimension("lon", len(nc_reference.dimensions["lon"]))
    if merra2_var_dict["cell_methods"]:
        nc1.createDimension("nv", 2)

    # TODO: Remove these unused comments from source code to somewhere else?
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
    time = nc1.createVariable("time", "i4", ("time",), zlib=True)
    time.axis = "T"
    time.units = "hours since 1980-01-01 00:00:00"
    time.long_name = "time"
    time.standard_name = "time"
    # 4.4.1. Calendar
    time.calendar = "gregorian"

    tbounds = None
    if merra2_var_dict["cell_methods"]:
        time.bounds = "time_bnds"
        tbounds = nc1.createVariable("time_bnds", "i4", ("time", "nv"))

    # TODO: Remove these unused comments from source code to somewhere else?
    # 4.3. Vertical (Height or Depth) Coordinate
    # level = nc1.createVariable('level','f4',('level',),zlib=True)
    # level.axis = 'Z'
    # level.units = 'Pa'
    # level.positive = 'up'
    # level.long_name = 'air_pressure'
    # level.standard_name = 'air_pressure'

    # 4.1. Latitude Coordinate
    lat = nc1.createVariable("lat", "f4", ("lat",), zlib=True)
    lat.axis = "Y"
    lat.units = "degrees_north"
    lat.long_name = "latitude"
    lat.standard_name = "latitude"
    lat[:] = nc_reference.variables["lat"][:]

    # 4.2. Longitude Coordinate
    lon = nc1.createVariable("lon", "f4", ("lon",), zlib=True)
    lon.axis = "X"
    lon.units = "degrees_east"
    lon.long_name = "longitude"
    lon.standard_name = "longitude"
    lon[:] = nc_reference.variables["lon"][:]

    least_digit = merra2_var_dict.get("least_significant_digit", None)
    var1 = nc1.createVariable(
        var_name,
        "f4",
        ("time", "lat", "lon"),
        zlib=True,
        chunksizes=(30, 120, 120),
        fill_value=deff4,
        least_significant_digit=least_digit,
    )
    # 3.1. Units
    var1.units = var_ref.units
    # 3.2. Long Name
    var1.long_name = var_ref.long_name
    # 3.3. Standard Name
    var1.standard_name = merra2_var_dict["standard_name"]
    if merra2_var_dict["cell_methods"]:
        var1.cell_methods = merra2_var_dict["cell_methods"]

    nc_reference.close()

    t = 0
    for i, division in enumerate(divided_files):
        tmp_data = ma.masked_all(
            [nt_division[i], len(nc1.dimensions["lat"]), len(nc1.dimensions["lon"])]
        )
        tmp_time = ma.masked_all([nt_division[i]])
        ttmp = 0
        for nc_file in division:
            if verbose:
                print(nc_file)
            nc = netCDF4.Dataset(nc_file, "r")
            ncvar = nc.variables[merra2_var_dict["merra_name"]]
            if time_from_filename:
                dt_str = nc_file.split(".")[-2]
                if len(dt_str) == 6:
                    dt_str += "01"
                ncdatetime = [
                    datetime.datetime(
                        int(dt_str[:4]), int(dt_str[4:6]), int(dt_str[6:])
                    )
                ]
            else:
                nctime = nc.variables["time"]
                ncdatetime = netCDF4.num2date(nctime[:], nctime.units)
            nctime_1980 = np.round(netCDF4.date2num(ncdatetime, time.units))
            tmp_data[ttmp : ttmp + ncvar.shape[0], :, :] = ncvar[:, :, :]
            tmp_time[ttmp : ttmp + ncvar.shape[0]] = nctime_1980[:]
            ttmp += ncvar.shape[0]
            nc.close()
        var1[t : t + tmp_data.shape[0], :, :] = tmp_data[:, :, :]
        time[t : t + tmp_data.shape[0]] = tmp_time[:]
        if merra2_var_dict["cell_methods"]:
            tbounds[t : t + tmp_data.shape[0], 0] = tmp_time[:] - 12
            tbounds[t : t + tmp_data.shape[0], 1] = tmp_time[:] + 12
        t += tmp_data.shape[0]

    nc1.close()


def daily_download_and_convert(
    merra2_server: str,
    var_names: Union[str, List[str]],
    initial_year: int,
    final_year: int,
    initial_month: int = 1,
    final_month: int = 12,
    initial_day: int = 1,
    final_day: Optional[int] = None,
    merra2_var_dicts: Optional[List[dict]] = None,
    output_dir: Union[str, Path] = None,
    delete_temp_dir: bool = True,
    verbose: bool = True,
):
    """MERRA2 daily download and conversion.

    Parameters
    ----------
    merra2_server : str
        Must contain trailing slash.
        e.g. https://goldsmr4.gesdisc.eosdis.nasa.gov/data/
    var_names : Union[str, List[str]]
        Variable short names, must be defined in pymerra2_variables.py
        if merra2_var_dict is not provided. If more than one variable,
        they are assumed to have the same original files and those will only
        be downloaded once.
    initial_year : int
    final_year : int
    initial_month : int
    final_month : int
    initial_day : int
    final_day : Optional[int]
        The range to download. final_day defaults to the end of the month.
    merra2_var_dicts : Optional[List[dict]]
        Dictionary containing the following keys:
        esdt_dir, collection, merra_name, standard_name,
        see the Bosilovich paper for details. Same order as var_names.
    output_dir : Union[str, Path]
    delete_temp_dir : bool
    verbose : bool

    Notes
    -----
    This is almost identical to subdaily_download_and_convert,
    perhaps not an ideal separation as this is duplicated code...

    """
    if isinstance(output_dir, Path):
        output_dir = Path(output_dir)
    if output_dir is None:
        output_dir = Path.cwd()

    if (2, 7) < sys.version_info < (3, 6):
        output_dir = str(output_dir)

    time_frequency = "day"

    if isinstance(var_names, str):
        var_names = [var_names]

    temp_dir_download = tempfile.mkdtemp(dir=output_dir)
    for i, var_name in enumerate(var_names):
        if not merra2_var_dicts:
            merra2_var_dict = var_list[var_name]
        else:
            merra2_var_dict = merra2_var_dicts[i]
        # Download subdaily files
        subdaily_download(
            merra2_server,
            merra2_var_dict["esdt_dir"],
            merra2_var_dict["collection"],
            initial_year,
            final_year,
            initial_month=initial_month,
            final_month=final_month,
            initial_day=initial_day,
            final_day=final_day,
            output_directory=temp_dir_download,
        )

        # Name the output file
        out_file_name = file_namer(
            var_name,
            time_frequency,
            initial_year,
            final_year,
            initial_month,
            final_month,
            initial_day,
            final_day,
        )

        out_file = Path(output_dir).joinpath(out_file_name)
        # Extract variable
        daily_netcdf(
            temp_dir_download,
            out_file,
            var_name,
            initial_year,
            final_year,
            verbose=verbose,
        )
    if delete_temp_dir:
        shutil.rmtree(temp_dir_download)


def _date_range_gen(
    init: datetime.datetime, end: datetime.datetime, freq: str = "daily"
):
    """Generator yielding dates (formatted as "YYYYMM[DD]"") and years in a range
    Iterates over days if freq=='daily', over months otherwise.
    """
    for year in range(init.year, end.year + 1):
        for month in range(
            init.month if year == init.year else 1,
            end.month + 1 if year == end.year else 13,
        ):
            if freq == "daily":
                init_day = (
                    init.day if (year == init.year) and (month == init.month) else 1
                )
                end_day = (
                    end.day
                    if (year == end.year) and (month == end.month)
                    else monthrange(year, month)[1]
                )
                for day in range(init_day, end_day + 1):
                    yield "{yyyy}{mm:02d}{dd:02d}".format(
                        yyyy=year, mm=month, dd=day
                    ), year
            elif freq == "monthly":
                yield "{yyyy}{mm:02d}".format(yyyy=year, mm=month), year


def download_from_url(
    url_template: str,
    freq: str,
    initial_year: int,
    final_year: int,
    initial_month: int = 1,
    final_month: int = 12,
    initial_day: int = 1,
    final_day: Optional[int] = None,
    merra2_names_map: Optional[dict] = None,
    output_dir: Union[str, Path] = None,
    delete_temp_dir: bool = True,
    verbose: bool = True,
    **kwargs
):
    """Download netCDF files from urls and merges them.

    This is for datasets that are outside of the raw merra2 output, but still related
    and similarly built.

    Parameters
    ----------
    url_template : str
        A string template for the download urls. Will be formatted with:
            date : a YYYYMM[DD] str, year (int), freq (same as below)
            and all kwargs (below)
    freq : str
        Either "daily" or "monthly", frequency of the files.
    initial_year : int
    final_year : int
    initial_month : int
    final_month : int
    initial_day : int
    final_day : Optional[int]
        The range to download. final_day defaults to the end of the month.
    merra2_names_map: Optional[dict]
        A mapping between names in the download files and names in the saved file.
        If given, restricts the merged file to those names.
    output_dir : Union[str, Path]
        The temporary directory and he final files will be put there.
    delete_temp_dir: bool
        Whether to delete the temporary download directory or not.
    verbose : bool
        Whether to print statuses to the screen, or not.

    Notes
    -----
    All other kwargs are passed to the url formatting and to the merging function.
    In the latter, they might be attributes missing in the files restricted to:

        Title, Institution, Source, Reference

    Returns
    -------
    List[Path]
        A list of the Paths of all netCDF files generated.
    """
    if isinstance(output_dir, Path):
        output_dir = Path(output_dir)
    if output_dir is None:
        output_dir = Path.cwd()

    if (2, 7) < sys.version_info < (3, 6):
        output_dir = str(output_dir)

    init_date = datetime.datetime(initial_year, initial_month, initial_day)
    end_date = datetime.datetime(
        final_year, final_month, final_day or monthrange(initial_year, initial_month)[1]
    )

    temp_dir_download = tempfile.mkdtemp(dir=output_dir)

    for date, year in _date_range_gen(init_date, end_date, freq=freq):
        subprocess.call(
            [
                "wget",
                "-c",
                "--directory-prefix={0}".format(temp_dir_download),
                "--load-cookies",
                str(Path("~/.urs_cookies").expanduser()),
                "--save-cookies",
                str(Path("~/.urs_cookies").expanduser()),
                "--keep-session-cookies",
                url_template.format(date=date, year=year, freq=freq, **kwargs),
            ]
        )

    nc_files = glob.glob(str(Path(temp_dir_download) / "*.nc"))
    nc_reference = netCDF4.Dataset(nc_files[0], "r")

    merra2_names_map = merra2_names_map or {}
    var_names = []
    merra2_var_dicts = []
    for var_name, variable in nc_reference.variables.items():
        if (merra2_names_map and var_name in merra2_names_map) or (
            var_name not in ["lat", "lon", "time"]
        ):
            var_names.append(merra2_names_map.get(var_name, var_name))
            merra2_var_dicts.append(
                dict(
                    merra_name=var_name,
                    standard_name=merra2_names_map.get(var_name, var_name),
                    **kwargs,
                )
            )

    try:
        netCDF4.num2date(nc_reference["time"][:], nc_reference["time"].units)
    except ValueError:
        time_from_filename = True
    else:
        time_from_filename = False

    nc_reference.close()

    out_files = []
    for var_name, merra2_var_dict in zip(var_names, merra2_var_dicts):
        out_file_name = file_namer(
            var_name,
            "day" if freq == "daily" else "month",
            initial_year,
            final_year,
            initial_month,
            final_month,
            initial_day,
            final_day,
        )
        out_file = Path(output_dir).joinpath(out_file_name)
        out_files.append(out_file)

        daily_netcdf(
            nc_files,
            output_file=out_file,
            var_name=var_name,
            initial_year=initial_year,
            final_year=final_year,
            merra2_var_dict=merra2_var_dict,
            time_from_filename=time_from_filename,
            verbose=verbose,
        )

    if delete_temp_dir:
        shutil.rmtree(temp_dir_download)

    return out_files
