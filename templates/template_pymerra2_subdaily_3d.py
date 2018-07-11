from calendar import monthrange

from pymerra2 import download

# Here we process multiple variables at a time to avoid downloading
# original data twice (all these variables are in the same files).
# These variables names are user choices, their merra-2 equivalent are
# specified below or in the default pymerra2_variables.py
var_names = ['hur']
# Frequency of the given variables (1hr, 3hr, 6hr)
time_frequency = '3hr'
delete_temp_dir = True
download_dir = '/path/to/output'
merra2_server = 'https://goldsmr5.sci.gsfc.nasa.gov/data/'

# The variables specification is in the same order as var_names above.
# esdt_dir, collection and merra_name can be found from
# https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf
# https://goldsmr5.sci.gsfc.nasa.gov/data/
# standard_name comes from
# http://cfconventions.org/standard-names.html
# Optionally, if all the variables are already in the default
# pymerra2_variables.py, this can be set to None.
merra2_var_dicts = [{'esdt_dir': 'M2T3NPCLD.5.12.4',
                     'collection': 'tavg3_3d_cld_Np',
                     'merra_name': 'RH',
                     'standard_name': 'relative_humidity',
                     'cell_methods': 'time: mean'}]

# This loop will create daily files of 3 hourly MERRA2 data
for yyyy in range(1980, 2017):
    for mm in range(1, 13):
        for dd in range(1, monthrange(yyyy, mm)[1] + 1):
            download.subdaily_download_and_convert(
                merra2_server, var_names, merra2_var_dicts=merra2_var_dicts,
                initial_year=yyyy, final_year=yyyy, initial_month=mm,
                final_month=mm, initial_day=dd, final_day=dd,
                output_dir=download_dir, delete_temp_dir=delete_temp_dir,
                time_frequency=time_frequency)
