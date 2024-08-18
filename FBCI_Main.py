import arcpy
from arcpy.sa import *
import requests as r
from datetime import date, datetime, timedelta
import time
import os
import pandas as pd
from arcpy import env
import numpy as np
from tqdm import tqdm
import re
from collections import defaultdict
from TVDI_calculator import tvdi
import warnings
env.overwriteOutput = True
##### Input
username = ""  # Update this line with the Appeears username
password = "" # Update this line with the Appeears password

name_set = False #If false the name of the task will be the date of the day
AHTN = True # Already Have Task Number
DownloadPresent = True #Already Have downloaded the task
task_id_input = "" #The task ID if a task already exists
Task_Name = "FBCITest" #The name of the task
DateStart = "12-07-24" #The start date of the task
DateEnd = "12-08-24" #The end date of the task
shapefile_path = r"F:\MemoireTest\SA_24Y\SA20Y.shp" #The shapefile that covers the study area
FuelMap = "" #The Fuel map over the study area, leave empty if not available
FuelMapEU = r'C:\Users\fouca\Documents\DocumentsDATA\GEO Master 2\Q1 Master 2\Memoire\FuelMap_LAEA\FuelMap2000_NFFL_LAEA.tif' #The fuel map at Europe scale

if name_set == False:
    today = datetime.now().date().strftime("%m-%d-%Y")
else:
    today = Task_Name
Dir = fr"F:\MemoireTest\{today}"
if not os.path.exists(Dir): os.makedirs(Dir)
inDir = os.path.join(Dir, "Input")  # Set up output directory using input directory and task name
outDir = os.path.join(Dir, "Output")
TVDIDir = os.path.join(outDir, "TVDI")
FBCIDir = os.path.join(outDir, "FBCI_FireMap")
CleanDir = os.path.join(outDir, "Cleaned")
if not os.path.exists(inDir): os.makedirs(inDir)  # Create the output directory
if not os.path.exists(outDir): os.makedirs(outDir)
if not os.path.exists(TVDIDir): os.makedirs(TVDIDir)
if not os.path.exists(FBCIDir): os.makedirs(FBCIDir)
if not os.path.exists(CleanDir): os.makedirs(CleanDir)
  # Get the input shapefile path
desc_c = arcpy.Describe(shapefile_path)
current_spatial_ref = desc_c.spatialReference
# Define the desired spatial reference (WGS84)
wgs84 = arcpy.SpatialReference(4326)  # EPSG code for WGS84

# Check if the current spatial reference is not WGS84
if current_spatial_ref.name != wgs84.name:
    # Define the output path for the reprojected shapefile
    output_shapefile_path = os.path.join(Dir, f"Extent_WGS84_{today}.shp")

    # Reproject the shapefile to WGS84
    arcpy.Project_management(shapefile_path, output_shapefile_path, wgs84)
    desc_new_shape = arcpy.Describe(output_shapefile_path)
    arcpy.AddMessage('Shapefile reprojected to WGS84 and saved as: {}'.format(output_shapefile_path))

else:
    arcpy.AddMessage('Shapefile is already in WGS84.')
    desc_new_shape = arcpy.Describe(shapefile_path)

ext=desc_new_shape.extent
ext.XMin, ext.YMin, ext.XMax, ext.YMax
bbox = {"west": ext.XMin, "south": ext.YMin, "east": ext.XMax, "north": ext.YMax, "crs": "EPSG:4326"}


###Date

startDateTime = datetime.strptime(DateStart, "%d-%m-%y")
startDatestr = startDateTime.strftime("%m-%d-%Y")

Year_start = startDateTime.year

endDateTime = datetime.strptime(DateEnd, "%d-%m-%y")
endDatestr = endDateTime.strftime("%m-%d-%Y")

Year_End = endDateTime.year
from datetime import datetime, timedelta
#Coordonnées de la zone d'étude

coordinates = [
    [
        [bbox["west"], bbox["south"]],
        [bbox["west"], bbox["north"]],
        [bbox["east"], bbox["north"]],
        [bbox["east"], bbox["south"]]
    ]
]

arcpy.env.extent = arcpy.Extent(bbox["west"], bbox["south"], bbox["east"],bbox["north"])
arcpy.env.workspace = inDir

###### APPEARS : Telechargement du LST

###LOGIN
if DownloadPresent == False:
    response = r.post('https://appeears.earthdatacloud.nasa.gov/api/login', auth=(username, password))
    token_response = response.json()
    print(token_response)
    token = token_response['token']                      # Save login token to a variable
    head = {'Authorization': 'Bearer {}'.format(token)}  # Create a header to store token information, needed to submit a request
    api = 'https://appeears.earthdatacloud.nasa.gov/api/'

# Création du task request

# Define common variables
Product1 = "MOD13A2.061"
Layer1_1 = "_1_km_16_days_NDVI"
Product2 = "MOD21A2.061"
Layer1_2 = "QC_Day"
Layer2_2 = "LST_Day_1KM"

# Create combined task JSON
task = {
    "params": {
        "geo": {
            "type": "FeatureCollection",
            "features": [{
                "type": "Feature",
                "geometry": {
                    "type": "Polygon",
                    "coordinates": coordinates
                },
                "properties": {}
            }]
        },
        "dates": [{
            "endDate": endDatestr,
            "startDate": startDatestr,
            "yearRange": [2000, 2050]
        }],
        "layers": [{
            "layer": Layer1_1,
            "product": Product1
        }, {
            "layer": Layer1_2,
            "product": Product2
        }, {
            "layer": Layer2_2,
            "product": Product2
        }],
        "output": {
            "format": {
                "type": "geotiff"
            },
            "projection": "geographic"
        }
    },
    "task_name": f"Task_{today}",
    "task_type": "area"
}
if DownloadPresent == False:
    if AHTN == False:
        token = token_response['token']
        response = r.post(
            'https://appeears.earthdatacloud.nasa.gov/api/task',
            json=task,
            headers={'Authorization': 'Bearer {0}'.format(token)})
        task_response = response.json()
        print(task_response)
    else:
        arcpy.AddMessage('Task already exists')

    starttime = time.time()
    if AHTN == False:
        task_id = task_response['task_id'] #récupération de l'ID propre à la task effectuée
    else:
        task_id = task_id_input
    while r.get('{}task/{}'.format(api, task_id), headers=head).json()['status'] != 'done':
        arcpy.AddMessage(r.get('{}task/{}'.format(api, task_id), headers=head).json()['status'])
        #arcpy.AddMessage(print(r.get('{}task/{}'.format(api, task_id), headers=head).json()['status']))
        time.sleep(20.0 - ((time.time() - starttime) % 20.0))
    arcpy.AddMessage(r.get('{}task/{}'.format(api, task_id), headers=head).json()['status'])
    #Permet de savoir quand la tâche est terminée

    ###Telechargement des données : sources

    bundle = r.get('{}bundle/{}'.format(api, task_id),headers=head).json()  # Call API and return bundle contents for the task_id as json
    #bundle  # Print bundle contents

    files = {}  # Create empty dictionary
    for f in bundle['files']: files[f['file_id']] = f['file_name']  # Fill dictionary with file_id as keys and file_name as values
    #files

    ###### DOWNLOAD
    total_files = len(files)
    progress_counter = 0

    with tqdm(total=total_files, desc="Downloading files") as pbar:
        for f in files:
            dl = r.get('{}bundle/{}/{}'.format(api, task_id, f), headers=head, stream=True, allow_redirects = 'True')                                # Get a stream to the bundle file
            if files[f].endswith('.tif'):
                filename = files[f].split('/')[1]
            else:
                filename = files[f]
            filepath = os.path.join(inDir, filename)                                                       # Create output file path
            with open(filepath, 'wb') as f:                                                                  # Write file to dest dir
                for data in dl.iter_content(chunk_size=8192): f.write(data)

            progress_counter += 1
            pbar.update(1)  # Update the progress bar
    print('Downloaded files can be found at: {}'.format(inDir))
else:
    print('Task has already been downloaded')
# List of LST_Day raster files
all_files = arcpy.ListRasters()

available_dates = []

from datetime import datetime, timedelta
if DownloadPresent == False:
    def julian_to_date(julian_day, base_year=2000):
        base_date = datetime(base_year, 1, 1)
        date = base_date + timedelta(days=julian_day - 1)
        return date.strftime('%m-%d')

    # Function to differentiate file type
    def differentiate_file_type(file_name):
        if "QC" in file_name:
            return "QC"
        elif "LST" in file_name:
            return "LST"
        elif "NDVI" in file_name:
            return "NDVI"
        elif "EVI" in file_name:
            return "EVI"
        elif "VI_Quality" in file_name:
            return "VI_Quality"
        else:
            return "Unknown"

    available_dates = []
    ndvi_dates = set()
    pattern = re.compile(r'doy(\d{4})(\d{3})')
    # List all rasters in the input directory
    all_files = arcpy.ListRasters()

    # Set the workspace environment to match the first raster's properties
    if all_files:
        first_raster = arcpy.Raster(os.path.join(inDir, all_files[0]))
        arcpy.env.extent = first_raster.extent
        arcpy.env.outputCoordinateSystem = arcpy.Describe(shapefile_path).spatialReference
        arcpy.env.cellSize = first_raster.meanCellWidth

    # First pass: Collect dates for NDVI files
    for i in all_files:
        match = pattern.search(i)
        if match:
            year = int(match.group(1))
            julian_day = int(match.group(2))
            if Year_start <= year <= Year_End:
                file_type = differentiate_file_type(i)
                if file_type == "NDVI":
                    new_date = f'{year}-{julian_to_date(julian_day, year)}'
                    ndvi_dates.add(new_date)

    # Second pass: Process files, only if there's a corresponding NDVI file
    for i in all_files:
        match = pattern.search(i)
        if match:
            year = int(match.group(1))
            julian_day = int(match.group(2))

            if Year_start <= year <= Year_End:
                file_type = differentiate_file_type(i)
                new_date = julian_to_date(julian_day, year)

                # Skip processing LST files without corresponding NDVI files
                #if file_type == "LST" and f'{year}-{new_date}' not in ndvi_dates:
                #    continue ==> #No because we need to calculate the average value across th epast month

                # Append the new date to the available_dates list
                available_dates.append(f'{year}-{new_date}')

                # Create the new file name
                new_file_name = f'{file_type}_Day_{year}-{new_date}.tif'
                old_path = os.path.join(inDir, i)
                new_path = os.path.join(inDir, new_file_name)

                if file_type == "LST":
                    # Read raster data for LST and apply the transformation
                    lst_data = (arcpy.Raster(old_path) * 0.02) - 273.15
                    lst_data.save(new_path)
                    arcpy.management.Delete(old_path)
                    print("Renamed {} to {}".format(i, new_file_name))
                    del lst_data
                elif file_type == "NDVI":
                    # Read raster data for NDVI and apply the transformation
                    ndvi_data = arcpy.Raster(old_path) * 0.0001
                    ndvi_data.save(new_path)
                    arcpy.management.Delete(old_path)
                    print("Renamed {} to {}".format(i, new_file_name))
                    del ndvi_data
                else:
                    # For QC files, check if the file exists before renaming
                    if os.path.exists(old_path):
                        os.rename(old_path, new_path)
                        print("Renamed {} to {}".format(i, new_file_name))
                    else:
                        print(f"File not found: {old_path}")
        else:
            print(f"No valid date found in file {i}")

    # Remove duplicates by converting the list to a set and back to a list
    unique_dates = list(set(available_dates))
    with open(os.path.join(Dir, 'valid_available_dates.txt'), 'w') as date_file:
        for date in unique_dates:
            date_file.write(date + '\n')
else:
    print("Importing dates from files")
    print(Dir)
    with open(os.path.join(Dir, 'valid_available_dates.txt'), 'r') as date_file:
        unique_dates = [line.strip() for line in date_file]
#re-assign the first raster
all_files = arcpy.ListRasters()
first_raster = arcpy.Raster(os.path.join(inDir, all_files[0]))
arcpy.env.extent = first_raster.extent
arcpy.env.outputCoordinateSystem = arcpy.Describe(shapefile_path).spatialReference
arcpy.env.cellSize = first_raster.meanCellWidth
valid_dates = []
valid_lst = []
print("Number of observations : {}".format(len(unique_dates)))

bad_lst_qc_values = {2, 5, 21, 53, 85, 116, 117, 164, 165, 181, 228, 229, 245, 177} #{2,5,16,17,21,33,48,49,53, 85, 116, 117, 164, 165, 180, 181, 228, 229, 244, 245}
#244 left out despite being labelled as "missing pixel", but does not seems to not be a missing pixel
# List of good VI Quality values for NDVI
good_vi_quality_values = {2112, 2116, 2120, 2372, 2376, 2380, 4160, 4164, 4168, 6208, 6212, 6216, 18504, 18508, 20548, 20552, 35144, 35148, 35152, 35156, 36936, 36940, 36944, 38984, 38988, 51540, 53328}
# Statistics list to collect data
statistics = []
arcpy.env.outputCoordinateSystem = wgs84

def remove_bad_pixels(lst_raster_path, qc_raster_path):
    try:
        # Read LST raster
        lst_raster = arcpy.Raster(lst_raster_path)
        lst_array = arcpy.RasterToNumPyArray(lst_raster, nodata_to_value=np.nan)

        # Read QC raster
        qc_raster = arcpy.Raster(qc_raster_path)
        qc_array = arcpy.RasterToNumPyArray(qc_raster)

        # Mask bad pixels in LST array
        mask = np.isin(qc_array, list(bad_lst_qc_values))
        lst_array[mask] = np.nan

        return lst_array

    except Exception as e:
        tqdm.write(f"Error processing LST raster: {e}")
        return None

def filter_vi_quality(ndvi_raster_path, vi_quality_raster_path):
    try:
        # Read NDVI raster
        ndvi_raster = arcpy.Raster(ndvi_raster_path)
        ndvi_array = arcpy.RasterToNumPyArray(ndvi_raster, nodata_to_value=np.nan)

        # Read VI Quality raster
        vi_quality_raster = arcpy.Raster(vi_quality_raster_path)
        vi_quality_array = arcpy.RasterToNumPyArray(vi_quality_raster)

        # Mask bad pixels in NDVI array
        mask = ~np.isin(vi_quality_array, list(good_vi_quality_values))
        ndvi_array[mask] = np.nan

        return ndvi_array

    except Exception as e:
        tqdm.write(f"Error processing NDVI raster: {e}")
        return None

def calculate_valid_percentage(intersection_array):
    valid_pixels = np.sum(~np.isnan(intersection_array))
    total_pixels = intersection_array.size
    return (valid_pixels / total_pixels) * 100

def save_raster(array, template_raster_path, output_raster_path):
    try:
        template_raster = arcpy.Raster(template_raster_path)
        modified_raster = arcpy.NumPyArrayToRaster(array, template_raster.extent.lowerLeft, template_raster.meanCellWidth, template_raster.meanCellHeight, np.nan)
        modified_raster.save(output_raster_path)
        return True
    except Exception as e:
        tqdm.write(f"Error saving raster {output_raster_path}: {e}")
        return False

statistics = []
valid_dates = []
valid_lst = []

# Clean and save all LST rasters
with tqdm(total=len(unique_dates), desc="Processing dates") as pbar:
    for date in unique_dates:
        lst_file = f'LST_Day_{date}.tif'
        qc_file = f'QC_Day_{date}.tif'
        ndvi_file = f'NDVI_Day_{date}.tif'
        vi_quality_file = f'VI_Quality_Day_{date}.tif'

        lst_path = os.path.join(inDir, lst_file)
        qc_path = os.path.join(inDir, qc_file)

        output_lst_path = os.path.join(CleanDir, f'LST_Day_{date}_Cleaned.tif')

        if os.path.exists(lst_path) and os.path.exists(qc_path):
            lst_array = remove_bad_pixels(lst_path, qc_path)
            valid_percentage_lst = calculate_valid_percentage(lst_array)
            if valid_percentage_lst > 10:
                lst_saved = save_raster(lst_array, lst_path, output_lst_path)
                if lst_saved:
                    valid_lst.append(date)
                    statistics.append([date, 'Accepted', valid_percentage_lst, 'LST Cleaned'])
                else:
                    statistics.append([date, 'Rejected', valid_percentage_lst, 'LST Save Failed'])
            else:
                statistics.append([date, 'Rejected', valid_percentage_lst, 'LST Bad Quality'])
        else:
            statistics.append([date, 'Rejected', 0, 'LST Missing'])

        pbar.update(1)

#write valid_lst to a file
with open(os.path.join(outDir, 'valid_lst_dates.txt'), 'w') as file:
    for date in valid_lst:
        file.write(f"{date}\n")

# Convert statistics to DataFrame and save as CSV
stats_df = pd.DataFrame(statistics, columns=['Date', 'Status', 'Valid_Percentage', 'Processing'])
stats_df.to_csv(os.path.join(outDir, 'valid_dates_statistics.csv'), index=False)

print(f"Processing completed for valid dates with more than 10% good pixels: {valid_lst}")
print(f"Statistics saved to {os.path.join(outDir, 'valid_dates_statistics.csv')}")

# Print or save the valid dates list
print("Valid dates for LST monthly average calculation:", len(valid_lst))

### Average LST calculation
arcpy.env.workspace = CleanDir
MAXIMUM_ALL = 0
MINIMUM_ALL = 100

# Dictionary to store rasters by date
rasters_by_prev_period = defaultdict(list)

for raster in valid_lst:
    # Extract the date from the raster file name
    date_obj = datetime.strptime(raster, "%Y-%m-%d")
    start_date = date_obj - timedelta(days=32)
    rasters_by_prev_period[raster] = [raster]

    # Iterate through all rasters to find those within the 32-day period
    for other_raster in valid_lst:
        other_date_obj = datetime.strptime(other_raster, "%Y-%m-%d")
        if start_date <= other_date_obj < date_obj:
            rasters_by_prev_period[raster].append(other_raster)

valid_dates = []  # Reset valid dates for the new criteria
statistics = []  # Reset statistics for the new criteria

# Calculate the average LST for the previous 32 days and check intersection with NDVI
with tqdm(total=len(rasters_by_prev_period), desc="Calculating avg LST and checking NDVI") as pbar:
    for target_date, rasters in rasters_by_prev_period.items():
        if rasters:
            raster_objects = [arcpy.Raster("LST_Day_" + raster + "_Cleaned.tif") for raster in rasters]
            avg_lst = arcpy.sa.CellStatistics(raster_objects, "MEAN")
            avg_output_path = os.path.join(CleanDir, f"avg_lst_32days_until_{target_date}.tif")
            avg_lst.save(avg_output_path)

            min_temp = arcpy.GetRasterProperties_management(avg_output_path, "MINIMUM").getOutput(0)
            min_temp = min_temp.replace(",", ".")
            if float(min_temp) < MINIMUM_ALL:
                MINIMUM_ALL = float(min_temp)

            max_temp = arcpy.GetRasterProperties_management(avg_output_path, "MAXIMUM").getOutput(0)
            max_temp = max_temp.replace(",", ".")
            if float(max_temp) > MAXIMUM_ALL:
                MAXIMUM_ALL = float(max_temp)

            ndvi_file = f'NDVI_Day_{target_date}.tif'
            vi_quality_file = f'VI_Quality_Day_{target_date}.tif'
            ndvi_path = os.path.join(inDir, ndvi_file)
            vi_quality_path = os.path.join(inDir, vi_quality_file)
            output_ndvi_path = os.path.join(CleanDir, f'NDVI_Day_{target_date}_Cleaned.tif')

            if os.path.exists(ndvi_path) and os.path.exists(vi_quality_path):
                ndvi_array = filter_vi_quality(ndvi_path, vi_quality_path)
                avg_lst_array = arcpy.RasterToNumPyArray(avg_output_path)
                intersection_array = np.where(~np.isnan(avg_lst_array) & ~np.isnan(ndvi_array), avg_lst_array, np.nan)
                valid_percentage = calculate_valid_percentage(intersection_array)

                if valid_percentage > 10:
                    ndvi_saved = save_raster(ndvi_array, ndvi_path, output_ndvi_path)
                    if ndvi_saved:
                        valid_dates.append(target_date)
                        statistics.append([target_date, 'Accepted', valid_percentage, 'Avg LST and NDVI Cleaned'])
                    else:
                        statistics.append([target_date, 'Rejected', valid_percentage, 'NDVI Save Failed'])
                else:
                    statistics.append([target_date, 'Rejected', valid_percentage, 'Bad Intersection Quality'])
            else:
                statistics.append([target_date, 'Rejected', 0, 'NDVI Missing'])
        pbar.update(1)

# Write final valid dates to a file
with open(os.path.join(outDir, 'filtered_valid_dates.txt'), 'w') as file:
    for date in valid_dates:
        file.write(f"{date}\n")

# Convert final statistics to DataFrame and save as CSV
final_stats_df = pd.DataFrame(statistics, columns=['Date', 'Status', 'Valid_Percentage', 'Processing'])
final_stats_df.to_csv(os.path.join(outDir, 'filters_statistics.csv'), index=False)

print(f"Final processing completed for valid dates with good avg LST and NDVI: {valid_dates}")
print(f"Final statistics saved to {os.path.join(outDir, 'filters_statistics.csv')}")

print("Valid available dates:", len(valid_dates))
print("Minimum temperature across all rasters:", MINIMUM_ALL)
print("Maximum temperature across all rasters:", MAXIMUM_ALL)

# Set the workspace to the destination directory for LST rasters
arcpy.env.workspace = CleanDir
list_new = arcpy.ListRasters()  # YYYY-MM-DD
# Example usage

if arcpy.Exists(FuelMap):
    FM = arcpy.Raster(FuelMap)
else:
    # Input raster
    ras = arcpy.Raster(FuelMapEU)

    # Define the remap ranges

    remap = RemapValue([
            (1, 74),
            (2, 400),
            (3, 300),
            (4, 1300),
            (5, 350),
            (6, 600),
            (7, 490),
            (8, 500),
            (9, 350),
            (10, 1200)
        ])
    # Reclassify the raster
    VT = Reclassify(ras, "Value", remap, "NODATA")
    #normalisation
    FuelMap = VT / 1300
    #Extract by mask over the extent of the shapefile
    FuelMap = ExtractByMask(FuelMap, shapefile_path)
    env.workspace = Dir
    FuelMap.save("FuelMap.tif")

env.extent = first_raster.extent
env.cellSize = first_raster.meanCellWidth
arcpy.env.snapRaster = first_raster

def compute_tvdi_for_dates(matching_dates, ndvi_rasters, lst_rasters):
    CLM_OK = 1
    warning_rasters = []
    dates_to_remove = []

    for selected_date in tqdm(matching_dates, desc="Processing TVDI"):
        try:
            arcpy.env.workspace = CleanDir
            sdatestr = selected_date
            selected_ndvi = [raster for raster in ndvi_rasters if sdatestr in raster][0]
            selected_lst = [raster for raster in lst_rasters if sdatestr in raster][0]
            NDVI_s = arcpy.Raster(selected_ndvi)

            io_inf = {'ndvi_file': selected_ndvi,
                      'ts_file': selected_lst,
                      'CLM_file': '',
                      'delta_file': '',
                      'output_dir': TVDIDir,
                      'ndvi_mult': 1,
                      'ts_mult': 1}

            roi_inf = {'geoflag': 2,
                       'x_pix': 0, 'y_pix': 0,
                       'dim_cols': 0, 'dim_rows': 0,
                       'moving_window': 0,
                       'window_size_x': 0, 'window_size_y': 0}

            alg_inf = {'dry_edge': 0,
                       'ts_min': 0,
                       'output': 0,
                       'dry_edge_params': [0.01, 0.1],  # [ndvi_step, ndvi_lower_limit]
                       'ts_min_params': 10,
                       'ts_min_file': "",
                       'output_params': [0.1, 0.9, 1.26,
                                         1]}  # [min_ndvi, max_ndvi, max_fi, use_1/delat_for_max_fi, interpolation]

            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                tvdi(io_inf, roi_inf, alg_inf)

                for warning in w:
                    if "Mean of empty slice" in str(warning.message):
                        warning_rasters.append(selected_ndvi)
                    else:
                        print(warning.message)
                        warning_rasters.append(selected_ndvi)
                        break

            arcpy.env.workspace = TVDIDir
            tvdi_fini = arcpy.Raster(f"NDVI_Day_{sdatestr}_Cleaned_TVDI.tif")

            # Apply the filter to remove pixels not between 0 and 1
            tvdi_filtered = arcpy.sa.SetNull((tvdi_fini < 0) | (tvdi_fini > 1), tvdi_fini)
            arcpy.management.DefineProjection(tvdi_filtered, NDVI_s.spatialReference)

            # Check if the filtered raster still contains negative values
            min_value = float(arcpy.GetRasterProperties_management(tvdi_filtered, "MINIMUM").getOutput(0))
            if min_value < 0:
                dates_to_remove.append(selected_date)
                print(f"Warning: Negative values found in TVDI raster for {selected_date} after filtering")
                continue

            tvdi_filtered.save(f"TVDI_{sdatestr}.tif")
            arcpy.management.Delete(f"NDVI_Day_{sdatestr}_Cleaned_TVDI.tif")
        except Exception as e:
            tqdm.write(f"Error processing {selected_date}: {e}")
            dates_to_remove.append(selected_date)
            continue

    return warning_rasters, dates_to_remove


# Example usage
ndvi_rasters = []
lst_rasters = []
avg_lst_rasters = []

for raster_name in list_new:
    if "NDVI_Day_" in raster_name:
        ndvi_rasters.append(raster_name)
    elif "LST_Day_" in raster_name:
        lst_rasters.append(raster_name)
    elif "avg_lst_32days_until_" in raster_name:
        avg_lst_rasters.append(raster_name)

warning_rasters, dates_to_remove = compute_tvdi_for_dates(valid_dates, ndvi_rasters, lst_rasters)

#Adding manually the dates that got an error
dates_to_remove = ['2013-03-06', '2022-05-09', '2024-01-17', '2005-03-06', '2015-08-13', '2003-02-02', '2012-02-18']
#remove the dates with exceptions from the valid_dates list
for date in dates_to_remove:
    if date in valid_dates:
        valid_dates.remove(date)

# Save the list of rasters with warnings
with open(os.path.join(outDir, 'rasters_with_warnings.txt'), 'w') as file:
    for raster in warning_rasters:
        file.write(f"{raster}\n")

print(f"Rasters with warnings: {len(warning_rasters)}")
print(f"Rasters that got an error and has been deleted: {len(dates_to_remove)} : {dates_to_remove}")
print(f"List of rasters with warnings saved to {os.path.join(outDir, 'rasters_with_warnings.txt')}")

arcpy.env.workspace = TVDIDir

# List all TVDI rasters .tif in the directory
all_TVDI = arcpy.ListRasters("*", "TIF")

# Define the mask using the shapefile path
mask_shapefile = r"C:\Users\fouca\Documents\ArcGIS\Projects\M_QA\RemoveArtefacts.shp"

# Loop through each TVDI raster and apply the ExtractByMask function
for tvdi_raster in all_TVDI:
    try:
        # Apply ExtractByMask
        extracted_raster = ExtractByMask(tvdi_raster, mask_shapefile)

        # Define the output path for the extracted raster
        output_raster_path = os.path.join(TVDIDir, f"Mask\Extracted_{tvdi_raster}")

        # Save the extracted raster
        extracted_raster.save(output_raster_path)

        print(f"Extracted and saved: {output_raster_path}")
    except Exception as e:
        print(f"Error processing {tvdi_raster}: {e}")

# Set the extent of the environment to the extent of the mask
arcpy.env.extent = arcpy.Describe(mask_shapefile).extent

import arcpy
from arcpy.sa import *
import os

def calculate_FBCI(selected_ndvi, selected_lst, tvdi_path, VT, avg_lst, MAXIMUM_ALL, MINIMUM_ALL):
    # Load NDVI and LST rasters
    print(f"Calculating FBCI for {selected_ndvi} and {selected_lst}")
    env.workspace = CleanDir
    NDVI = arcpy.Raster(selected_ndvi)
    # Load TVDI raster
    tvdi_fini = arcpy.Raster(tvdi_path)
    avg_lst_r = arcpy.Raster(avg_lst)

    #Calculating normalised LST values:

    normalized_LST = (avg_lst_r - MINIMUM_ALL) / (MAXIMUM_ALL - MINIMUM_ALL)

    # Calculate FBCI using TVDI and normalized LST
    FBCI = ((normalized_LST + tvdi_fini - NDVI + (VT)) + 1) / 4
    print(f"FBCI range: {arcpy.GetRasterProperties_management(FBCI, 'MINIMUM').getOutput(0)} - {arcpy.GetRasterProperties_management(FBCI, 'MAXIMUM').getOutput(0)}")
    #FBCI_Scale = (FBCI/FBCI.min)* 100

    # Save the output raster
    output_path = os.path.join(FBCIDir, f"FBCI_{os.path.splitext(os.path.basename(tvdi_path))[0]}.tif")
    FBCI.save(output_path)

from datetime import datetime

# Folder containing TVDI rasters
tvdifolder = TVDIDir
# Calculate FBCI for each selected date
for selected_date in tqdm(valid_dates, desc="Creating FBCI rasters"):
    selected_ndvi = [raster for raster in ndvi_rasters if selected_date in raster][0]
    selected_lst = [raster for raster in lst_rasters if selected_date in raster][0]
    tvdi_path = os.path.join(tvdifolder, f"TVDI_{selected_date}.tif")
    select_avg_lst = [raster for raster in avg_lst_rasters if selected_date in raster][0]
    calculate_FBCI(selected_ndvi, selected_lst, tvdi_path, FuelMap, select_avg_lst, MAXIMUM_ALL, MINIMUM_ALL)

#Create an average FBCI out of every FBCI available
arcpy.env.workspace = FBCIDir
all_FBCI = arcpy.ListRasters()
FBCI_rasters = []
for raster in all_FBCI:
    FBCI_rasters.append(arcpy.Raster(raster))
FBCI_avg = CellStatistics(FBCI_rasters, "MEAN")
FBCI_avg.save(os.path.join(outDir, "FBCI_Average.tif"))
### Average for each year
#Create a list of years
years = []
for date in valid_dates:
    year = datetime.strptime(date, "%Y-%m-%d").year
    if year not in years:
        years.append(year)

for year in years:
    FBCI_rasters = []
    for raster in all_FBCI:
        if f"{year}-" in raster:
            FBCI_rasters.append(arcpy.Raster(raster))
    FBCI_avg = CellStatistics(FBCI_rasters, "MEAN")
    FBCI_avg.save(os.path.join(outDir, f"FBCI_Average_{year}.tif"))

arcpy.env.workspace = TVDIDir
all_TVDI = arcpy.ListRasters()
TVDI_rasters = []
for raster in all_TVDI:
    TVDI_rasters.append(arcpy.Raster(raster))
TVDI_avg = CellStatistics(TVDI_rasters, "MEAN")
TVDI_avg.save(os.path.join(outDir, "TVDI_Average.tif"))


#Using extract by mask to crop out the artefacts on the TVDI and FBCI rasters
arcpy.env.workspace = Dir
arcpy.env.extent = arcpy.Describe(shapefile_path).extent
arcpy.env.cellSize = first_raster.meanCellWidth
shapefile_ExtractArt = r"C:\Users\fouca\Documents\ArcGIS\Projects\M_QA\RemoveArtefacts.shp"""
#Extract by mask for TVDI

TVDI = arcpy.Raster(os.path.join(outDir, "TVDI_Average.tif"))
TVDI = ExtractByMask(TVDI, shapefile_path)
TVDI.save(os.path.join(outDir, "TVDI_Average.tif"))

#Extract by mask for FBCI

FBCI = arcpy.Raster(os.path.join(outDir, "FBCI_Average.tif"))
FBCI = ExtractByMask(FBCI, shapefile_path)
FBCI.save(os.path.join(outDir, "FBCI_Average.tif"))

##LOGOUT -> Disconnect from Appeears server
if DownloadPresent == False:
    token = token_response['token']
    response = r.post(
        'https://appeears.earthdatacloud.nasa.gov/api/logout',
        headers={'Authorization': 'Bearer {0}'.format(token)})
    print(response.status_code)

