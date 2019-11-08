import glob
import os
import shutil
import tempfile

from pymerra2 import download

# This assumes that for each year, every month is available and that the
# subset files downloaded have the .SUB.nc4 suffix and are otherwise named
# exactly as the raw files.
# This has only been tested with variable subsetting (i.e. the files still
# contain the full domain)

path_data = "/path/to/downloaded/subset/files"
path_output = "/path/to/output/results"

year_ini = 1980
year_fin = 2017
var_name = "tas"  # This is not the merra2 variable name, see var_info below.
freq = "1hr"

# Look at pymerra2_variables.py for examples of that this should be:
var_info = {
    "cell_methods": None,
    "collection": "inst1_2d_asm_Nx",
    "esdt_dir": None,  # Not required here...
    "least_significant_digit": 3,
    "merra_name": "T2M",
    "standard_name": "air_temperature",
}

output_file_template = "{0}_{1}_merra2_reanalysis_{2}{3}.nc"

for yyyy in range(year_ini, year_fin + 1):
    for mm in range(1, 13):
        # list files and rename them (actually a copy) without the .SUB.
        search_str = "*{0}{1}*.nc4".format(str(yyyy), str(mm).zfill(2))
        nc_files = glob.glob(os.path.join(path_data, search_str))
        path_tmp = tempfile.mkdtemp(dir=path_output)
        for nc_file in nc_files:
            renamed_file = os.path.basename(nc_file).replace(".SUB.", ".")
            shutil.copy(nc_file, os.path.join(path_tmp, renamed_file))

        output_file = output_file_template.format(
            var_name, freq, str(yyyy), str(mm).zfill(2)
        )
        download.subdaily_netcdf(
            path_tmp,
            os.path.join(path_output, output_file),
            var_name,
            yyyy,
            yyyy,
            var_info,
            True,
        )

        shutil.rmtree(path_tmp)
