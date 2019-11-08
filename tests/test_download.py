# import fnmatch
# import glob
# import os
# import shutil
# import sys
# import tempfile
# import pymerra2.download
#

# TODO: Figure out a reasonable way to test download capabilities?
# TODO: Add an example MERRA2 dataset to perfrom conversions on?

#
# def test_daily_download_convert():
#     var_names = ['tasmin', 'tasmax', 'prmax']
#     delete_temp_dir = True
#     download_dir = tempfile.mkdtemp()
#     merra2_server = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/'
#
#     pymerra2.download.daily_download_and_convert(
#         merra2_server, var_names, merra2_var_dicts=None, initial_year=2015,
#         final_year=2015, initial_month=1, final_month=1, initial_day=1,
#         final_day=None, output_dir=download_dir,
#         delete_temp_dir=delete_temp_dir)
#
#     if sys.version_info >= (3, 5):
#         daily_files = glob.glob(os.path.join(download_dir, '**', '*.nc'), recursive=True)
#         print(len(daily_files))
#         try:
#             assert len(daily_files) == 3
#         finally:
#             shutil.rmtree(download_dir)
#
#     if sys.version_info < (3, 5):
#         daily_files = []
#         for root, dirnames, filenames in os.walk(download_dir):
#             for filename in fnmatch.filter(filenames, '*.nc'):
#                 daily_files.append(os.path.join(root, filename))
#         try:
#             assert len(daily_files) == 3
#         finally:
#             shutil.rmtree(download_dir)
