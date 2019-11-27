import logging
from pathlib import Path
from pymerra2 import download

# Here we process multiple variables at a time to avoid downloading
# original data twice (all these variables are in the same files).
# These variables names are user choices, their merra-2 equivalent are
# specified below or in the default pymerra2_variables.py
var_names = ["evspsbl", "huss", "prbc", "tas", "sic", "snw", "uas", "vas", "ps"]
# Frequency of the given variables (1hr, 3hr, 6hr)
time_frequency = "1hr"
delete_temp_dir = False
download_dir = Path.cwd().joinpath("downloaded")
merra2_server = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/"

# The variables specification is in the same order as var_names above.
# esdt_dir, collection and merra_name can be found from
# https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf
# https://goldsmr4.gesdisc.eosdis.nasa.gov/data/
# standard_name comes from
# http://cfconventions.org/standard-names.html
# Optionally, if all the variables are already in the default
# pymerra2_variables.py, this can be set to None.

# This loop will create monthly files of hourly MERRA2 data
for yyyy in range(2017, 2019):
    for mm in range(1, 13):
        try:
            download.subdaily_download_and_convert(
                merra2_server,
                var_names,
                merra2_var_dicts=None,
                initial_year=yyyy,
                final_year=yyyy,
                initial_month=mm,
                final_month=mm,
                initial_day=1,
                final_day=None,
                output_dir=download_dir,
                delete_temp_dir=delete_temp_dir,
                time_frequency=time_frequency,
            )
        except Exception as e:
            msg = "{}: File not found".format(e)
            logging.error(msg)
            continue


def file_namer(
    var_name,
    initial_year,
    final_year,
    initial_month,
    final_month,
    initial_day: int = None,
    final_day: int = None,
):
    # Name the output file
    if (initial_year == final_year) and (initial_month == final_month):
        if initial_day == final_day:
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
                var_name,
                time_frequency,
                str(initial_year),
                str(initial_month).zfill(2),
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


for yyyy in range(2017, 2019):
    for mm in range(1, 13):
        try:

            download.subdaily_netcdf(
                path_data=download_dir,
                output_file="toto1.nc",
                var_name="pr",
                initial_year=yyyy,
                final_year=yyyy,
                merra2_var_dict=None,
                verbose=False,
            )
        except Exception as e:
            msg = "{}: File not found".format(e)
            logging.error(msg)
            continue
