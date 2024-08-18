import arcpy
from arcpy.sa import *
import os

# Set environment settings
arcpy.env.workspace = "F:/MemoireTest/SA_24Y/Output/FDCI_FireMap"
arcpy.env.overwriteOutput = True

# List of fire dates
fire_dates = ["2011-04-26", "2023-05-30", "2015-04-23", "2015-04-21", "2018-04-19", "2018-04-19",
              "2019-02-27", "2019-02-27", "2019-02-26", "2022-03-09"]

# Get all FDCI rasters from directory
all_FDCI = arcpy.ListRasters("FDCI_TVDI*.tif")

# Extract dates from filtered FDCI rasters and sort them
fdci_dates = sorted([raster.split('_')[2].replace('.tif', '') for raster in all_FDCI])

# Determine global min and max values for consistent legend
global_min = 0.18
global_max = 0.67

# Path to fire shapefile
fire_shapefile_path = 'C:/Users/fouca/Documents/DocumentsDATA/GEO Master 2/Q1 Master 2/Memoire/CEMS/f1e14c6597114f5ebf68360b27fd6c86.shp'

# Define a reference raster for consistent extent and resolution
ref_raster = arcpy.Raster(all_FDCI[0])
ref_extent = ref_raster.extent
ref_cellsize = ref_raster.meanCellWidth

# Set up project and layout
aprx = arcpy.mp.ArcGISProject(r"C:\Users\fouca\Documents\ArcGIS\Projects\M_QA\M_QA.aprx")
layouts = [layout.name for layout in aprx.listLayouts()]
#select the layout PlotFBCI
layout = aprx.listLayouts("PlotFBCI")[0]


# Create a new map frame for each fire date
for fire_date in fire_dates:
    fire_date_obj = fire_date

    # Find the four previous observations
    previous_observations = [date for date in fdci_dates if date < fire_date_obj][-4:]

    # Check for matching FDCI rasters
    found_rasters = [raster for date in previous_observations for raster in all_FDCI if date in raster]

    # Align and resample data from the four previous FDCI rasters
    aligned_rasters = []
    for raster_path in found_rasters:
        data = arcpy.Raster(raster_path)
        aligned_data = arcpy.sa.Resample(data, arcpy.env.workspace + "/aligned_" + os.path.basename(raster_path),
                                         ref_cellsize, "BILINEAR")
        aligned_rasters.append(aligned_data)

    # Add the rasters and fire shapefile to the map
    for i, aligned_raster in enumerate(aligned_rasters):
        map_frame = layout.listElements("MAPFRAME_ELEMENT")[0]
        new_layer = map_frame.map.addDataFromPath(aligned_raster)

        # Apply symbology
        symbology = new_layer.symbology
        symbology.updateRenderer('RasterStretchRenderer')
        symbology.renderer.colorRamp = aprx.listColorRamps('Viridis')[0]
        symbology.renderer.stretchType = 'MinMax'
        symbology.renderer.minMaxValues = [(global_min, global_max)]
        new_layer.symbology = symbology

        # Add fire shapefile
        fire_layer = map_frame.map.addDataFromPath(fire_shapefile_path)

        # Set symbology for fire layer
        fire_symbology = fire_layer.symbology
        fire_symbology.updateRenderer('SimpleRenderer')
        fire_symbology.renderer.symbol.color = {'RGB': [255, 0, 0, 100]}
        fire_layer.symbology = fire_symbology

        # Set extent
        map_frame.camera.setExtent(ref_extent)

        # Add title
        title_text = "FBCI values before the fire of " + fire_date
        title_element = layout.listElements("TEXT_ELEMENT", "Title")[0]
        title_element.text = title_text

        # Export to PNG
        output_path = os.path.join("C:/path/to/output_directory", f"FBCI_Before_Fire_{fire_date}.png")
        layout.exportToPNG(output_path, resolution=300)

        # Remove layers for next iteration
        map_frame.map.removeLayer(new_layer)
        map_frame.map.removeLayer(fire_layer)

print("Process complete")

import arcpy
import os
import datetime
# Function to set symbology of the raster layers
def set_raster_symbology(layer, min_val, max_val):
    symbology = layer.symbology
    symbology.updateColorizer('RasterStretchColorizer')

    # Get the color ramp from the project
    aprx = arcpy.mp.ArcGISProject(project_path)
    color_ramp = aprx.listColorRamps('Inferno')[0]

    symbology.colorizer.colorRamp = color_ramp
    symbology.colorizer.minValue = min_val
    symbology.colorizer.maxValue = max_val
    symbology.transparency = 50  # Set transparency to 50%
    layer.symbology = symbology

# Paths and variables
project_path = r'C:\Users\fouca\Documents\ArcGIS\Projects\PlotsFBCI\PlotsFBCI.aprx'
layout_name = 'Layout'
map_frame_names = ['Map Frame', 'Map Frame 1', 'Map Frame 2', 'Map Frame 3']
raster_folder = r'F:\MemoireTest\SA_24Y\Output\FDCI_FireMap'
output_folder = r'C:\Users\fouca\Documents\DocumentsDATA\GEO Master 2\Q1 Master 2\Memoire\Figures\Cartes\FBCI'
global_min = 0.18
global_max = 0.67

# List of fire dates
fire_dates = ['2011-04-26', '2023-05-30', '2015-04-23', '2015-04-21', '2018-04-19',
              '2018-04-19', '2019-02-27', '2019-02-27', '2019-02-26', '2022-03-09']

# Get all FDCI rasters from directory
all_fdci_rasters = arcpy.ListRasters("FDCI_TVDI*.tif")



# Extract dates from filtered FDCI rasters and sort them
fdci_dates = sorted([os.path.basename(f).split('_')[2].split('.tif')[0] for f in all_fdci_rasters])

# Open the existing project
aprx = arcpy.mp.ArcGISProject(project_path)
layout = aprx.listLayouts(layout_name)[0]

for fire_date in fire_dates:
    fire_date_obj = datetime.datetime.strptime(fire_date, '%Y-%m-%d')

    # Find the four previous observations
    previous_observations = sorted([date for date in fdci_dates if datetime.datetime.strptime(date, '%Y-%m-%d') < fire_date_obj], reverse=True)[:4]

    # Check for matching FDCI rasters
    found_rasters = [raster for date in previous_observations for raster in all_fdci_rasters if date in raster]
    found_rasters = [os.path.join(raster_folder, raster) for raster in found_rasters]

    # Debugging: Print the names of found rasters
    print(found_rasters)

    # Clear existing layers in map frames
    for map_frame_name in map_frame_names:
        map_frame = layout.listElements('MAPFRAME_ELEMENT', map_frame_name)[0]
        for lyr in map_frame.map.listLayers(): #except if its a base map
            if lyr.name != 'World Topographic Map':
                map_frame.map.removeLayer(lyr)


    # Add rasters to the map frames
    for i, raster_path in enumerate(found_rasters):
        map_frame = layout.listElements('MAPFRAME_ELEMENT', map_frame_names[i])[0]
        raster_layer = map_frame.map.addDataFromPath(raster_path)
        set_raster_symbology(raster_layer, global_min, global_max)
        raster_layer.transparency = 50  # Set transparency to 50%

    # Update the title text element
    title_text_element = layout.listElements('TEXT_ELEMENT')[1]
    title_text_element.text = f"FBCI values before the fire of {fire_date}"

    # Export layout to PNG
    output_path = os.path.join(output_folder, f"FBCI_before_fire_{fire_date.replace('-', '')}.png")
    layout.exportToPNG(output_path, 300)

# Save the project
aprx.save()
del aprx

import arcpy
import os
import datetime

# Function to set symbology of the raster layers
def set_raster_symbology(layer, min_val, max_val):
    symbology = layer.symbology
    symbology.updateColorizer('RasterStretchColorizer')

    # Get the color ramp from the project
    aprx = arcpy.mp.ArcGISProject(project_path)
    color_ramp = aprx.listColorRamps('Inferno')[0]

    symbology.colorizer.colorRamp = color_ramp
    symbology.colorizer.minValue = min_val
    symbology.colorizer.maxValue = max_val
    symbology.transparency = 30  # Set transparency to 50%
    layer.symbology = symbology

# Paths and variables
project_path = r'C:\Users\fouca\Documents\ArcGIS\Projects\PlotsFBCI\PlotsFBCI.aprx'
layout_name = 'Layout'
map_frame_names = ['Map Frame', 'Map Frame 1', 'Map Frame 2', 'Map Frame 3']
raster_folder = r'F:\MemoireTest\SA_24Y\Output\FDCI_FireMap'
output_folder = r'C:\Users\fouca\Documents\DocumentsDATA\GEO Master 2\Q1 Master 2\Memoire\Figures\Cartes\FBCI'
global_min = 0.18
global_max = 0.67
fire_shapefile_path = r'C:/Users/fouca/Documents/DocumentsDATA/GEO Master 2/Q1 Master 2/Memoire/CEMS/f1e14c6597114f5ebf68360b27fd6c86.shp'

# List of fire dates
fire_dates = ['2011-04-26', '2023-05-30', '2015-04-23', '2015-04-21', '2018-04-19',
              '2018-04-19', '2019-02-27', '2019-02-27', '2019-02-26', '2022-03-09']

# Get all FDCI rasters from directory
all_fdci_rasters = arcpy.ListRasters("FDCI_TVDI*.tif")

# Extract dates from filtered FDCI rasters and sort them
fdci_dates = sorted([os.path.basename(f).split('_')[2].split('.tif')[0] for f in all_fdci_rasters])

# Open the existing project
aprx = arcpy.mp.ArcGISProject(project_path)
layout = aprx.listLayouts(layout_name)[0]

# Add fire shapefile to the map frames
for map_frame_name in map_frame_names:
    map_frame = layout.listElements('MAPFRAME_ELEMENT', map_frame_name)[0]
    map_frame.map.addDataFromPath(fire_shapefile_path)

for fire_date in fire_dates:
    fire_date_obj = datetime.datetime.strptime(fire_date, '%Y-%m-%d')

    # Find the four previous observations
    previous_observations = sorted([date for date in fdci_dates if datetime.datetime.strptime(date, '%Y-%m-%d') < fire_date_obj], reverse=True)[:4]

    # Check for matching FDCI rasters
    found_rasters = [raster for date in previous_observations for raster in all_fdci_rasters if date in raster]
    found_rasters = [os.path.join(raster_folder, raster) for raster in found_rasters]

    # Debugging: Print the names of found rasters
    print(found_rasters)

    # Clear existing layers in map frames (except base maps and fire shapefile)
    for map_frame_name in map_frame_names:
        map_frame = layout.listElements('MAPFRAME_ELEMENT', map_frame_name)[0]
        for lyr in map_frame.map.listLayers():
            if lyr.name != 'World Topographic Map' and lyr.name != "f1e14c6597114f5ebf68360b27fd6c86" or lyr.name == "Burned Area":
                map_frame.map.removeLayer(lyr)

    # Add rasters to the map frames
    for i, raster_path in enumerate(found_rasters):
        map_frame = layout.listElements('MAPFRAME_ELEMENT', map_frame_names[i])[0]
        raster_layer = map_frame.map.addDataFromPath(raster_path)
        set_raster_symbology(raster_layer, global_min, global_max)
        raster_layer.transparency = 30  # Set transparency to 50%

    # Apply definition query to the fire shapefile layer
    fire_layer = None
    for lyr in map_frame.map.listLayers():
        print(lyr.name)
        if lyr.name == "f1e14c6597114f5ebf68360b27fd6c86" or lyr.name == "Burned Area":
            fire_layer = lyr
            break

    if fire_layer:
        fire_layer.definitionQuery = f"initialdat = timestamp '{fire_date}'"
    #make sure the layer is on top of the raster

    # Move the fire layer to the top of the map frame
    for lyr in map_frame.map.listLayers():
        if lyr.name == "Burned Area":
            map_frame.map.moveLayer(lyr, 'TOP')


    # Update the title text element
    title_text_element = layout.listElements('TEXT_ELEMENT')[1]
    title_text_element.text = f"FBCI values before the fire of {fire_date}"
    # Loop through the map frames and update the date text elements
    for i in range(4):
        title_element = layout.listElements("TEXT_ELEMENT", f"Date {i + 1}")[0]
        # Set the date of the corresponding map frame raster
        title_element.text = f"Date : {previous_observations[i]}"

    # Export layout to PNG
    output_path = os.path.join(output_folder, f"FBCI_before_fire_{fire_date.replace('-', '')}.png")
    layout.exportToPNG(output_path, 300)
    #clear the definition query

    if fire_layer:
        fire_layer.definitionQuery = ""

# Save the project
#aprx.save()
#del aprx
import arcpy
import os
import datetime

# Function to set symbology of the raster layers
def set_raster_symbology(layer, min_val, max_val):
    symbology = layer.symbology
    symbology.updateColorizer('RasterStretchColorizer')

    # Get the color ramp from the project
    aprx = arcpy.mp.ArcGISProject(project_path)
    color_ramp = aprx.listColorRamps('Inferno')[0] #inverse the color ramp


    symbology.colorizer.colorRamp = color_ramp
    symbology.colorizer.minValue = min_val
    symbology.colorizer.maxValue = max_val
    symbology.transparency = 50  # Set transparency to 50%
    layer.symbology = symbology

# Paths and variables
project_path = r'C:\Users\fouca\Documents\ArcGIS\Projects\PlotsFBCI\PlotsFBCI2.aprx'
layout_name = 'Layout'
map_frame_names = ['Map Frame', 'Map Frame 1', 'Map Frame 2', 'Map Frame 3']
raster_folder = r'F:\MemoireTest\SA_24Y\Output\FDCI_FireMap'
output_folder = r'C:\Users\fouca\Documents\DocumentsDATA\GEO Master 2\Q1 Master 2\Memoire\Figures\Cartes\FBCI'
global_min = 0.18
global_max = 0.67
fire_layer_path = r'C:/Users/fouca/Documents/DocumentsDATA/GEO Master 2/Q1 Master 2/Memoire/CEMS/Burned Area.lyrx'

# List of fire dates
fire_dates = ['2011-04-26', '2023-05-30','2015-04-21', '2018-04-19',
              '2018-04-19', '2019-02-27', '2019-02-27', '2022-03-09', '2018-07-30']

# Open the existing project
aprx = arcpy.mp.ArcGISProject(project_path)
layout = aprx.listLayouts(layout_name)[0]

# Add fire shapefile to the map frames if not already added
for map_frame_name in map_frame_names:
    map_frame = layout.listElements('MAPFRAME_ELEMENT', map_frame_name)[0]
    if not any(layer.name == "Burned Area" for layer in map_frame.map.listLayers()):
        fire_layer_template = arcpy.mp.LayerFile(fire_layer_path)
        map_frame.map.addLayer(fire_layer_template, "TOP")

for fire_date in fire_dates:
    fire_date_obj = datetime.datetime.strptime(fire_date, '%Y-%m-%d')

    # Get all FDCI rasters from directory
    arcpy.env.workspace = raster_folder
    all_fdci_rasters = arcpy.ListRasters("FDCI_TVDI*.tif")

    # Extract dates from filtered FDCI rasters and sort them
    fdci_dates = sorted([os.path.basename(f).split('_')[2].split('.tif')[0] for f in all_fdci_rasters])
    #Remove the dates with too few pixels; 2015-02-02, 2015-01-01, 2019-02-02, 2018-12-19
    fdci_dates = [date for date in fdci_dates if date not in ['2015-02-02', '2015-01-01', '2019-02-02', '2018-12-19']]

    # Find the four previous observations
    previous_observations = sorted([date for date in fdci_dates if datetime.datetime.strptime(date, '%Y-%m-%d') < fire_date_obj], reverse=True)[:4]

    # Check for matching FDCI rasters
    found_rasters = [raster for date in previous_observations for raster in all_fdci_rasters if date in raster]
    found_rasters = [os.path.join(raster_folder, raster) for raster in found_rasters]

    # Debugging: Print the names of found rasters
    print(found_rasters)

    # Clear existing layers in map frames (except base maps and fire shapefile)
    for map_frame_name in map_frame_names:
        map_frame = layout.listElements('MAPFRAME_ELEMENT', map_frame_name)[0]
        layers_to_remove = [lyr for lyr in map_frame.map.listLayers() if lyr.name != 'World Topographic Map' and lyr.name != "Burned Area"]
        for lyr in layers_to_remove:
            map_frame.map.removeLayer(lyr)

    # Add rasters to the map frames
    for i, raster_path in enumerate(found_rasters):
        map_frame = layout.listElements('MAPFRAME_ELEMENT', map_frame_names[i])[0]
        raster_layer = map_frame.map.addDataFromPath(raster_path)
        set_raster_symbology(raster_layer, global_min, global_max)
        raster_layer.transparency = 50  # Set transparency to 50%

    # Apply definition query to the fire shapefile layer
    for map_frame_name in map_frame_names:
        map_frame = layout.listElements('MAPFRAME_ELEMENT', map_frame_name)[0]
        fire_layer = next((lyr for lyr in map_frame.map.listLayers() if lyr.name == "Burned Area"), None)
        if fire_layer:
            fire_layer.definitionQuery = f"initialdat = timestamp '{fire_date}'"

    # Make sure the layer is on top of the raster
    for map_frame_name in map_frame_names:
        map_frame = layout.listElements('MAPFRAME_ELEMENT', map_frame_name)[0]
        fire_layer = next((lyr for lyr in map_frame.map.listLayers() if lyr.name == "Burned Area"), None)
        if fire_layer:
            for lyr in map_frame.map.listLayers():
                if lyr.name != 'World Topographic Map' and lyr.name != "Burned Area":
                    reference_layer = lyr
                    map_frame.map.moveLayer(reference_layer, fire_layer, "BEFORE")
                    break

    # Update the title text element
    title_text_element = layout.listElements('TEXT_ELEMENT', "Text")[0]
    title_text_element.text = f"FBCI values before the fire of {fire_date}"

    # Loop through the map frames and update the date text elements
    for i in range(4):
        title_element = layout.listElements("TEXT_ELEMENT", f"Date {i + 1}")[0]
        # Set the date of the corresponding map frame raster
        title_element.text = f"{previous_observations[i]}"

    # Export layout to PNG
    output_path = os.path.join(output_folder, f"FBCI_before_fire_{fire_date.replace('-', '')}.png")
    layout.exportToPNG(output_path, 300)

    # Clear the definition query
    for map_frame_name in map_frame_names:
        map_frame = layout.listElements('MAPFRAME_ELEMENT', map_frame_name)[0]
        fire_layer = next((lyr for lyr in map_frame.map.listLayers() if lyr.name == "Burned Area"), None)
        if fire_layer:
            fire_layer.definitionQuery = ""

# Save the project
aprx.save()
del aprx

out_raster = arcpy.sa.Reclassify("FDCI_Raster", "VALUE", "0,299658 0,400000 1;0,400000 0,500000 2;0,500000 0,600000 3;0,600000 7 4", "DATA"); out_raster.save(r"C:\Users\fouca\Documents\ArcGIS\Projects\M_QA\M_QA.gdb\Reclass_FDCI_AllValues")

out_raster = arcpy.sa.Reclassify("wui_sa", "VALUE", "0 0;1 10;2 10;3 10;4 10;5 0;6 0;8 0", "DATA"); out_raster.save(r"C:\Users\fouca\Documents\ArcGIS\Projects\M_QA\M_QA.gdb\Reclass_WUI_1")

out_raster = arcpy.sa.Reclassify("FDCI_Moy_Clean", "VALUE", "0,299658 0,368530 1;0,368530 0,395407 2;0,395407 0,423124 3;0,423124 0,455040 4;0,455040 0,513833 5", "DATA"); out_raster.save(r"C:\Users\fouca\Documents\ArcGIS\Projects\M_QA\M_QA.gdb\Reclass_FDCI2")

with arcpy.EnvManager(outputCoordinateSystem='GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]', snapRaster="FDCI_Moy_Clean", extent='5.8666666661411 50.4166666621501 6.29999999943561 50.6916666621254 GEOGCS["GCS_WGS_1984".DATUM["D_WGS_1984".SPHEROID["WGS_1984".6378137.0.298.257223563]].PRIMEM["Greenwich".0.0].UNIT["Degree".0.0174532925199433]]', cellSize="DangerMap", mask="FDCI_Moy_Clean"):
    output_raster = arcpy.ia.RasterCalculator(' "Reclass_WUI_1" + "Reclass_FDCI2"'); output_raster.save(r"c:\Users\fouca\documents\ArcGIS\Projects\M_QA\M_QA.gdb\DangerMap")

out_raster = arcpy.sa.Reclassify("FDCI_Moy_Clean", "VALUE", "0,299658 0,368530 1;0,368530 0,395407 2;0,395407 0,423124 3;0,423124 0,455040 4;0,455040 0,513833 5", "DATA"); out_raster.save(r"C:\Users\fouca\Documents\ArcGIS\Projects\M_QA\M_QA.gdb\Reclass_FDCI2")

import arcpy
import os
import datetime
arcpy.env.overwriteOutput = True
from arcpy.sa import *

# Function to set symbology of the raster layers
def set_raster_symbology(layer, color_list):
    sym = layer.symbology
    sym.updateColorizer('RasterUniqueValueColorizer')
    colorizer = sym.colorizer
    #colorizer.colorRamp = None  # Clear any existing color ramp
    if colorizer.type == "RasterUniqueValueColorizer":
        colorizer.field = 'Value'
        items = colorizer.groups[0].items
        if len(items) > 0:
            items[0].label = "≤ 0.39"
            items[0].color = {'RGB': [253, 255, 165, 255]}
        if len(items) > 1:
            items[1].label = "0.4 - 0.49"
            items[1].color = {'RGB': [250, 141, 10, 255]}
        if len(items) > 2:
            items[2].label = "0.5 - 0.59"
            items[2].color = {'RGB': [188, 55, 85, 255]}
        if len(items) > 3:
            items[3].label = "0.6 - 0.69"
            items[3].color = {'RGB': [87, 16, 110, 255]}
        if len(items) > 4:
            items[4].label = "≥ 0.7"
            items[4].color = {'RGB': [0, 0, 4, 255]}
    layer.symbology = sym

# Paths and variables
project_path = r'C:\Users\fouca\Documents\ArcGIS\Projects\PlotsFBCI\PlotsFBCI4.aprx'
layout_name = 'Layout'
map_frame_names = ['Map Frame', 'Map Frame 1', 'Map Frame 2', 'Map Frame 3']
raster_folder = r'F:\MemoireTest\SA_24Y\Output\FDCI_FireMap'
output_folder = r'C:\Users\fouca\Documents\DocumentsDATA\GEO Master 2\Q1 Master 2\Memoire\Figures\Cartes\FBCI'
reclass_raster_folder = r'F:\MemoireTest\SA_24Y\Output\Reclass_FDCI_FireMap'
fire_layer_path = r'C:/Users/fouca/Documents/DocumentsDATA/GEO Master 2/Q1 Master 2/Memoire/CEMS/Burned Area.lyrx'
color_list = ["#007206", "#7DB810", "#F2FE1E", "#FFAC12", "#FC3B09"]

# Create the reclass raster folder if it doesn't exist
if not os.path.exists(reclass_raster_folder):
    os.makedirs(reclass_raster_folder)

# List of fire dates
fire_dates = ['2011-04-26', '2023-05-30', '2015-04-21', '2018-04-19',
              '2018-04-19', '2019-02-27', '2019-02-27', '2022-03-09', '2018-07-30']

# Open the existing project
aprx = arcpy.mp.ArcGISProject(project_path)
layout = aprx.listLayouts(layout_name)[0]

# Add fire shapefile to the map frames if not already added
for map_frame_name in map_frame_names:
    map_frame = layout.listElements('MAPFRAME_ELEMENT', map_frame_name)[0]
    if not any(layer.name == "Burned Area" for layer in map_frame.map.listLayers()):
        fire_layer_template = arcpy.mp.LayerFile(fire_layer_path)
        map_frame.map.addLayer(fire_layer_template, "TOP")

for fire_date in fire_dates:
    fire_date_obj = datetime.datetime.strptime(fire_date, '%Y-%m-%d')

    # Get all FDCI rasters from directory
    arcpy.env.workspace = raster_folder
    all_fdci_rasters = arcpy.ListRasters("FDCI_TVDI*.tif")

    # Extract dates from filtered FDCI rasters and sort them
    fdci_dates = sorted([os.path.basename(f).split('_')[2].split('.tif')[0] for f in all_fdci_rasters])
    fdci_dates = [date for date in fdci_dates if date not in ['2015-02-02', '2015-01-01', '2019-02-02', '2018-12-19']]

    # Find the four previous observations
    previous_observations = sorted(
        [date for date in fdci_dates if datetime.datetime.strptime(date, '%Y-%m-%d') < fire_date_obj], reverse=True)[:4]

    # Check for matching FDCI rasters
    found_rasters = [raster for date in previous_observations for raster in all_fdci_rasters if date in raster]
    found_rasters = [os.path.join(raster_folder, raster) for raster in found_rasters]

    # Debugging: Print the names of found rasters
    print(found_rasters)

    # Clear existing layers in map frames (except base maps and fire shapefile)
    for map_frame_name in map_frame_names:
        map_frame = layout.listElements('MAPFRAME_ELEMENT', map_frame_name)[0]
        layers_to_remove = [lyr for lyr in map_frame.map.listLayers() if
                            lyr.name != 'World Topographic Map' and lyr.name != "Burned Area"]
        for lyr in layers_to_remove:
            map_frame.map.removeLayer(lyr)

    # Reclassify and add rasters to the map frames
    for i, raster_path in enumerate(found_rasters):
        reclass_map = "0 0,400000 1;0,400000 0,500000 2;0,500000 0,600000 3;0,600000 0,700000 4;0,700000 1 5"
        reclass_raster = arcpy.sa.Reclassify(raster_path, "VALUE", reclass_map, "DATA")
        reclass_raster_path = os.path.join(reclass_raster_folder, f"reclass_{os.path.basename(raster_path)}")
        reclass_raster.save(reclass_raster_path)

        map_frame = layout.listElements('MAPFRAME_ELEMENT', map_frame_names[i])[0]
        raster_layer = map_frame.map.addDataFromPath(reclass_raster_path)
        set_raster_symbology(raster_layer, color_list)
        raster_layer.transparency = 50  # Set transparency to 50%

    # Apply definition query to the fire shapefile layer
    for map_frame_name in map_frame_names:
        map_frame = layout.listElements('MAPFRAME_ELEMENT', map_frame_name)[0]
        fire_layer = next((lyr for lyr in map_frame.map.listLayers() if lyr.name == "Burned Area"), None)
        if fire_layer:
            fire_layer.definitionQuery = f"initialdat = timestamp '{fire_date}'"

    # Make sure the layer is on top of the raster
    for map_frame_name in map_frame_names:
        map_frame = layout.listElements('MAPFRAME_ELEMENT', map_frame_name)[0]
        fire_layer = next((lyr for lyr in map_frame.map.listLayers() if lyr.name == "Burned Area"), None)
        if fire_layer:
            for lyr in map_frame.map.listLayers():
                if lyr.name != 'World Topographic Map' and lyr.name != "Burned Area":
                    reference_layer = lyr
                    map_frame.map.moveLayer(reference_layer, fire_layer, "BEFORE")
                    break

    # Update the title text element
    title_text_element = layout.listElements('TEXT_ELEMENT', "Text")[0]
    title_text_element.text = f"FBCI values\n before the fire\n of {fire_date}"


    # Loop through the map frames and update the date text elements
    for i in range(4):
        title_element = layout.listElements("TEXT_ELEMENT", f"Date {i + 1}")[0]
        title_element.text = f"{previous_observations[i]}"

    # Export layout to PNG
    output_path = os.path.join(output_folder, f"FBCI_before_fire_only2vert_{fire_date.replace('-', '')}.png")
    layout.exportToPNG(output_path, 300)

    # Clear the definition query
    for map_frame_name in map_frame_names:
        map_frame = layout.listElements('MAPFRAME_ELEMENT', map_frame_name)[0]
        fire_layer = next((lyr for lyr in map_frame.map.listLayers() if lyr.name == "Burned Area"), None)
        if fire_layer:
            fire_layer.definitionQuery = ""

print("Process complete")

del aprx