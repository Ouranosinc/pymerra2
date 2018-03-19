using PyCall
unshift!(PyVector(pyimport("sys")["path"]), "")
@pyimport merra2

# Here we process multiple variables at a time to avoid downloading the
# original data twice (all these variables are in the same files).
# These variables names are user choices, their merra-2 equivalent are
# specified below or in the default merra2_variables.py
var_names = ["tasmin"]#, "tasmax","prmax"]
delete_temp_dir = true
download_dir = "/path/to/output"
merra2_server = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/"

# The variables specification is in the same order as var_names above.
# esdt_dir, collection and merra_name can be found from
# https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf
# https://goldsmr4.gesdisc.eosdis.nasa.gov/data/
# standard_name comes from
# http://cfconventions.org/standard-names.html
# Optionally, if all the variables are already in the default
# merra2_variables.py, this can be set to None.
global merra2_var_dicts
merra2_var_dicts = [Dict("esdt_dir"=> "M2SDNXSLV.5.12.4",
                     "collection"=> "statD_2d_slv_Nx",
                     "merra_name"=> "T2MMAX",
                     "standard_name"=> "air_temperature",
                     "cell_methods"=> "time: max"),
                    Dict("esdt_dir"=> "M2SDNXSLV.5.12.4",
                     "collection"=> "statD_2d_slv_Nx",
                     "merra_name"=> "T2MMIN",
                     "standard_name"=> "air_temperature",
                     "cell_methods"=> "time: min"),
                    Dict("esdt_dir"=> "M2SDNXSLV.5.12.4",
                     "collection"=> "statD_2d_slv_Nx",
                     "merra_name"=> "TPRECMAX",
                     "standard_name"=> "precipitation_flux",
                     "cell_methods"=> "time: max")]

# This loop will create annual files of daily MERRA2 data
for yyyy in range(2017,2017)
    merra2.daily_download_and_convert(
        merra2_server, var_names, merra2_var_dicts=nothing, initial_year=yyyy,
        final_year=yyyy, initial_month=1, final_month=12, initial_day=1,
        final_day=nothing, download_dir=download_dir,
        delete_temp_dir=delete_temp_dir)
end
