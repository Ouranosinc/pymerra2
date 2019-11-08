import os

from pymerra2 import download

# Need to download the constant files on disk from
# https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/
path_data = "/path/to/constant_files"
output_path = "/path/to/output"

# These variables names are user choices, their merra-2 equivalent are
# specified below or in the default pymerra2_variables.py
var_names = ["phis", "sftgif", "sftlf", "sftof"]

# The variables specification is in the same order as var_names above.
# esdt_dir, collection and merra_name can be found from
# https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf
# https://goldsmr4.gesdisc.eosdis.nasa.gov/data/
# standard_name comes from
# http://cfconventions.org/standard-names.html
# Optionally, if all the variables are already in the default
# pymerra2_variables.py, this can be set to None.
merra2_var_dicts = [
    {
        "esdt_dir": "M2C0NXASM.5.12.4",
        "collection": "const_2d_asm_Nx",
        "merra_name": "PHIS",
        "standard_name": "surface_geopotential",
        "cell_methods": None,
    },
    {
        "esdt_dir": "M2C0NXASM.5.12.4",
        "collection": "const_2d_asm_Nx",
        "merra_name": "FRLANDICE",
        "standard_name": "land_ice_area_fraction",
        "cell_methods": None,
    },
    {
        "esdt_dir": "M2C0NXASM.5.12.4",
        "collection": "const_2d_asm_Nx",
        "merra_name": "FRLAND",
        "standard_name": "land_area_fraction",
        "cell_methods": None,
    },
    {
        "esdt_dir": "M2C0NXASM.5.12.4",
        "collection": "const_2d_asm_Nx",
        "merra_name": "FROCEAN",
        "standard_name": "sea_area_fraction",
        "cell_methods": None,
    },
]

for i, var_name in enumerate(var_names):
    output_file = os.path.join(
        output_path, "{0}_fx_merra2_reanalysis.nc".format(var_name)
    )
    if merra2_var_dicts:
        merra2_var_dict = merra2_var_dicts[i]
    else:
        merra2_var_dict = None
    download.fixed_netcdf(
        path_data, output_file, var_name, merra2_var_dict=merra2_var_dict
    )
