import logging
from pathlib import Path
from pymerra2 import download

# Here we process multiple variables at a time to avoid downloading
# original data twice (all these variables are in the same files).
# These variables names are user choices, their merra-2 equivalent are
# specified below or in the default pymerra2_variables.py
var_names = ["evspsbl", "huss", "prbc", "tas", "sic", "snw", "uas", "vas", "ps"]
var_names = ["hur"]
delete_temp_dir = False
download_dir = Path.cwd().joinpath("downloaded")
merra2_server = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/"
merra2_server = "https://goldsmr5.gesdisc.eosdis.nasa.gov/data/"

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
            )
        except Exception as e:
            msg = "{}: File not found".format(e)
            logging.error(msg)
            continue
