import merra2

# Here we process 2 variables at a time to avoid downloading the original
# data twice (both are in the same file)
var_names = ['tasmin', 'tasmax','prmax']
delete_temp_dir = True
download_dir = '/scen3/stdenis/projects/climate_datasets/merra2/data'

# This loop will create annual files of hourly MERRA2 data
for yyyy in range(1980, 2017):
    merra2.daily_download_and_convert(
        var_names, merra2_var_dicts=None, initial_year=yyyy,
        final_year=yyyy, initial_month=1, final_month=12, initial_day=1,
        final_day=None, output_dir=download_dir,
        delete_temp_dir=delete_temp_dir)
