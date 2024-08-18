# FBCI
Arcpy-based calculation of a Fire Behaviour Composite Index

The code is composed of two Python scripts for processing and analyzing raster data, with
the final output being the FBCI. Below is a detailed explanation of the code, including the
input parameters and expected outputs.

# Input Parameters
name_set: Boolean flag to determine if a custom task name is set.
AHTN: Boolean flag indicating if the task number is already available.
DonwloadPresent: Boolean flag indicating if the data has already been downloaded.
task_id_input: String representing the task ID if the task is already created.
Task_Name: String representing the custom task name.
DateStart: String representing the start date in the format dd-mm-yy.
DateEnd: String representing the end date in the format dd-mm-yy.
shapefile_path: String representing the path to the input shapefile.
roads_shapefile: String representing the path to the roads shapefile.
FuelMap: String representing the path to the fuel map (optional).
FuelMapEU: String representing the path to the European fuel map raster.
Source of the Fuel Map for EU;
=> EFFIS Fuel Map
https://forest-fire.emergency.copernicus.eu/applications/data-and-services

# Setting Up the Environment
Directories: The script sets up several directories for input, output, and intermediate results.
Spatial Reference: The input shapefile is reprojected to WGS84 if it is not already in that coordinate system.
Date Handling: The script parses the start and end dates and sets up the bounding box for the study area.
# Data Download and Processing
Login and Task Creation: If DonwloadPresent is False, the script logs into the NASA API, creates a task, and downloads the data.
Data Cleaning: The script processes the downloaded LST and NDVI rasters, removing bad pixels based on quality control (QC) values.
Valid Dates: The script identifies valid dates for LST and NDVI data and saves them to text files.
# Calculations and Analysis
Average LST Calculation: The script calculates the average LST for the previous 32 days for each valid date.
TVDI Calculation: The script calculates the Temperature Vegetation Dryness Index (TVDI) for each valid date.
FBCI Calculation: The script calculates the Fire Danger Condition Index (FBCI) using the TVDI, normalized LST, and NDVI.
# Outputs
Reprojected Shapefile: If the input shapefile is not in WGS84, it is reprojected and saved with a new name.
Downloaded Data: Various raster files (e.g., LST, NDVI) are downloaded from the NASA API.
Cleaned LST Rasters: LST rasters are processed to remove bad pixels and saved in a cleaned format.
Valid Dates Files: Lists of valid dates for LST and NDVI are saved to text files (valid_lst_dates.txt and final_valid_dates.txt).
Statistics CSVs: Processing statistics are saved to CSV files (raster_processing_statistics.csv and final_raster_processing_statistics.csv).
Average LST and NDVI Rasters: Average LST and NDVI rasters are calculated and saved.
TVDI and FBCI Rasters: TVDI and FBCI rasters are calculated and saved, including an average FBCI raster.
Cleaned TVDI Raster: The TVDI raster is cleaned using a mask and saved.
Exported Layouts: Layouts are exported to PNG files for each fire date.
Log Files: Various log messages and status updates are printed to the console or saved to files.
# Summary
Before running the script, ensure that the input parameters are correctly set, including paths to shapefiles, task names, and date ranges. The script will process the data, calculate the TVDI & FBCI, then save the results in the specified directories. The main outputs include cleaned and processed rasters, calculated indices, overall average of the indices and .csv with the statistics of the filtered data. The TVDI is based on https://github.com/hectornieto/pyTVDI/tree/master, but adapted to use Arcpy. 

Important : If the script fails in any way before finishing calculating the outputs, be sure to delete all the files created in the "Output" folder, then run the code again with the correct input parameters.
