import merra2

var_name = 'uas'
delete_temp_dir = True
download_dir = '/scen3/stdenis/projects/climate_datasets/merra2/data'

# This loop will create monthly files of hourly MERRA2 data
for yyyy in range(1980, 1981):
    for mm in range(1, 13):
        merra2.subdaily_download_and_convert(
            var_name, merra2_var_dict=None, initial_year=yyyy, final_year=yyyy,
            initial_month=mm, final_month=mm, initial_day=1, final_day=None,
            output_dir=download_dir, delete_temp_dir=delete_temp_dir)
