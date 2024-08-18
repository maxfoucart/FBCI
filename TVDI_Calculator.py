# run the TVDI
import arcpy
import os
import numpy as np
# methods for calcualting ts min
cMEAN = 0
cMEDIAN = 1
cVAR_MAX_NDVI = 2
cCONST_TS_MIN = 3
cDAILY_TS_MIN = 4

# methods for calculating dry edge
cSIMPLE = 0
cTANG = 1

# output variable
cTVDI = 0
cEF = 1

# interpolation
cSQUARE = 0
cLINEAR = 1

# geoflag values
cGEOG_COORD = 0
cPIXEL_COORD = 1
cWHOLE_IMG = 2

# moving window or tiles
cNO_WINDOW = 0
cTILES = 1
cMOVING_WINDOW = 2
CLM_OK = 1
def tvdi(io_inf, roi_inf, alg_inf):
    # Open files
    cEF = 1
    if not arcpy.Exists(io_inf['ndvi_file']):
        print(io_inf['ndvi_file'] + ' not found')
        return
    if not arcpy.Exists(io_inf['ts_file']):
        print(io_inf['ts_file'] + ' not found')
        return
    if not arcpy.Exists(io_inf['CLM_file']):
        clm_fid = -1
    if alg_inf['output'] == cEF and arcpy.Exists(io_inf['delta_file']):
        delta = arcpy.Raster(io_inf['delta_file']).to_numpy()
    else:
        delta = 0
    if not arcpy.Exists(io_inf['output_dir']):
        os.mkdir(io_inf['output_dir'])
    ndvi_fid = arcpy.Raster(io_inf['ndvi_file'])
    ts_fid = arcpy.Raster(io_inf['ts_file'])
    #clm_fid = arcpy.Raster(io_inf['CLM_file'])
    # Get Geotransform
    geo = ts_fid.extent
    prj = ts_fid.spatialReference
    # Read and process data
    if roi_inf['geoflag'] == 2:
        ndvi = arcpy.RasterToNumPyArray(ndvi_fid, nodata_to_value= np.nan)
        ts = arcpy.RasterToNumPyArray(ts_fid, nodata_to_value=np.nan)
        if clm_fid:
            cld_mask = arcpy.sa.Con(arcpy.sa.IsNull(ndvi_fid), 0, 1)
        else:
            cld_mask = arcpy.sa.Con(arcpy.sa.IsNull(ndvi_fid), 0, 1)
        dims = [0, 0, ndvi.shape[0] - 1 , ndvi.shape[1] - 1 ] # -1 not needed with arpcy
    else:
        dims = [roi_inf['y_pix'], roi_inf['x_pix'],
                roi_inf['y_pix'] + roi_inf['dim_rows'], roi_inf['x_pix'] + roi_inf['dim_cols']]
        ndvi = ndvi_fid[dims[0]:dims[2], dims[1]:dims[3]].to_numpy()
        ts = ts_fid[dims[0]:dims[2], dims[1]:dims[3]].to_numpy()
        if clm_fid:
            cld_mask = clm_fid[dims[0]:dims[2], dims[1]:dims[3]].to_numpy()
        else:
            cld_mask = arcpy.sa.Con(arcpy.sa.IsNull(ndvi), 0, 1)
        if alg_inf['output'] == cEF and arcpy.Exists(io_inf['delta_file']):
            delta = delta[dims[0]:dims[2], dims[1]:dims[3]]
    # Closing NDVI and LST files
    del ndvi_fid, clm_fid
    ndvi = ndvi / float(io_inf['ndvi_mult'])
    ts = ts / float(io_inf['ts_mult'])
    nan_mask = np.isnan(ndvi) | np.isnan(ts)
    # Bad pixel if ndvi being lower than the lower limit set in algorithm options or higher than 1 or temperature is negative
    ndvi[ndvi < alg_inf['dry_edge_params'][1]] = np.nan
    ts[ndvi < alg_inf['dry_edge_params'][1]] = np.nan
    ndvi[ts < 0] = np.nan
    ts[ts < 0] = np.nan
    #cld_mask_array = arcpy.RasterToNumPyArray(cld_mask)
#bro
    # Now you can perform operations on the NumPy array
    CLM_OK = 1
    #ndvi[cld_mask_array != CLM_OK] = 0.0
    #ts[cld_mask_array != CLM_OK] = 0.0
    # Create a text file where the edge line equations will be saved
    filename = os.path.splitext(os.path.basename(io_inf['ndvi_file']))[0]
    tvdi_out_file = os.path.join(io_inf['output_dir'], filename + "_TVDI")
    with open(tvdi_out_file + '_line_equations.txt', 'w') as lun:
        string = 'image \t dry_edge_intercept \t dry_edge_slope \t wet_edge_intercept\t wet_edge_slope \t r^2 \t dry_edge_points \t total_points \t points_per_bin \t ndvi_range_ratio\n'
        lun.write(string)
        lun.flush()
    #Mask

    # Calculate TVDI using the method specified
    tvdi_array = triangle_window(ndvi, ts, delta, dims, roi_inf, alg_inf, tvdi_out_file)
    tvdi_array[nan_mask] = np.nan
    # Get rid of bad or cloudy pixels in the output
    #tvdi[(ndvi < alg_inf['dry_edge_params'][1]) | (ts < 0)] = -1
    #tvdi[cld_mask_array != CLM_OK] = -1
    # Write GDAL tvdi file
    # Set output name
    out_name = tvdi_out_file + '.tif'
    # Create the output GeoTIFF file
    width = ts_fid.meanCellWidth
    height = ts_fid.meanCellHeight
    lowerLeft = arcpy.Point(ts_fid.extent.XMin, ts_fid.extent.YMin)
    arcpy.NumPyArrayToRaster(tvdi_array, lowerLeft, width, height).save(out_name)
    print('Done: TVDI saved in ' + out_name)
    #print("File used; pyTVDI_Backup_OSJ")


def triangle_window(ndvi, ts, delta, dims, roi_inf, alg_inf, tvdi_out_file):
    ''' Calculates the TVDI/EF from numpy arrays and returns the TVDI as another array.
    It also writes an ASCII file with the TVDI statistics

    Parameters
    ----------
    ndvi : array_like
        Vegetation Index array.
    ts : array_like
        Surface Temperature array.
    delta : array_like
        Delta/Deta+psicr array
    dims : list
        Dimensions of the Region of Interest: [row min, column min, row max, column max]
    roi_inf : dict
        specifies the Region of Interest.
        - geoflag = cGEOG_COORD if you use geographic coordinates.
        - geoflag = cPIXEL_COORD if you use ns, nl.
        - geoflag = cWHOLE_IMG if you want to use the whole picture.
        x_pix, y_pix : float
            the coordinate of the top left corner of ROI, units depend on geoflag.
        dim_rows, dim_cols : int
            the size, in pixels, of the ROI.
        moving_window : bool
            set to 1 if moving window is to be used, 0 if tiles are to be used.
        window_size_x, window_size_y : int
            size of the window/tiles, in pixels, if 0 then the whole ROI is one window.
    alg_inf : dict
        specifies the algorithm to use to calcualte the triangle
        Set dry_edge = cSIMPLE or cTANG depending on which dry edge algorithm is to be used
        Set ts_min = cMEAN, cMEDIAN, cVAR_MAX_NDVI or cCONST_TS_MIN depending on which ts min algorithm is to be used
        dry_edge_params - array containing paramenters required for the selected dry edge algorithm:
            - [ndvi step, ndvi lower limit]
        ts_min_params - array containing parameters required for the selected ts min algorithm:
            [[number of max NDVI values to use, N/A]              - cMEAN
            [number of max NDVI values to use, N/A]              - cMEDIAN
            [given value of NDVI, scaling factor for max NDVI]   - cVAR_MAX_NDVI
            [given value of TS min, N/A]]                        - cCONST_TS_MIN
            [daily Ts min value read from file, N/A]]             - cDAILY_TS_MIN
        output_params - array containing paremeters required for calculating EF:
            [min NDVI, max NDVI, max FI, use 1/delta for max FI, interpolation type]
    tvdi_out_file : str
        Name of the the output ASCII file
    Returns
    -------
    output : array_like
        TVDI  or Evaporative Fraction array

    '''
    from numpy import zeros, amin
    from re import search
    from os.path import dirname
    # Creat the output array
    output = zeros(ndvi.shape)
    ts_min = 0
    # Read algorithm options for alg_inf struct
    dry_edge_method = alg_inf['dry_edge']
    ts_min_method = alg_inf['ts_min']
    ts_min_params = alg_inf['ts_min_params']
    ndvi_step = alg_inf['dry_edge_params'][0]
    ndvi_lower_limit = alg_inf['dry_edge_params'][1]
    output_var = alg_inf['output']
    # EF options
    ndvi_min = alg_inf['output_params'][0]
    ndvi_max = alg_inf['output_params'][1]
    fi_max = alg_inf['output_params'][2]
    cSQUARE = 0
    # use square or linear interpolation
    if alg_inf['output_params'][3] == cSQUARE:
        power = 2
    else:
        power = 1
    # Initialize counting variables for loops in sub-images
    row = dims[0]
    col = dims[1]
    # Use the whole image if window_size_y and window_size_y are set to 0
    if roi_inf['window_size_y'] == 0:
        roi_inf['window_size_y'] = ndvi.shape[0]
    if roi_inf['window_size_x'] == 0:
        roi_inf['window_size_x'] = ndvi.shape[1]
    # Vertical loop in sub-images
    while row < dims[2]:
        # Horizontal loop in sub-images
        while col < dims[3]:
            window_size_y = np.nanmin([roi_inf['window_size_y'], dims[2] - row])
            window_size_x = np.nanmin([roi_inf['window_size_x'], dims[3] - col])
            # Select the window data
            window_col = col - dims[1]
            window_row = row - dims[0]
            print('col: ' + str(window_col) + '   row: ' + str(window_row))
            ndvi_window = ndvi[window_row:window_row + window_size_y, window_col:window_col + window_size_x]
            ts_window = ts[window_row:window_row + window_size_y, window_col:window_col + window_size_x]
            output_window = zeros(ndvi_window.shape)
            if delta:
                delta_window = delta[window_row:window_row + window_size_y, window_col:window_col + window_size_x]
            # ==============================================================================
            #             Calculate the dry edge and ts_min and TVDI
            # ==============================================================================
            # set plot output file name
            out_name = tvdi_out_file + '_cols_' + str(col) + '_' + str(col + window_size_x - 1) + '_rows_' + str(
                row) + '_' + str(row + window_size_y - 1) + '.png'
            if dry_edge_method == cSIMPLE:
                [lin_fit, ts_min, fit_stats] = calc_triangle_simple(ndvi_window,
                                                                    ts_window, ndvi_step, ndvi_lower_limit,
                                                                    ts_min_method, ts_min_params,
                                                                    plot_out_file=out_name)
            if dry_edge_method == cTANG:
                [lin_fit, ts_min, fit_stats] = calc_triangle_tang(ndvi_window,
                                                                  ts_window, ndvi_step, ndvi_lower_limit, ts_min_method,
                                                                  ts_min_params, plot_out_file=out_name)
            # if dry edge couldn't be found just output error value
            if lin_fit[0] == 0 and lin_fit[1] == 0:
                print('Could not calcualte dry edge for ' + out_name)
                output_window = zeros(ndvi_window.shape)
                output_window = -1
            # otherwise calculate the required output
            else:
                # write edge line equations to a text file
                lun = open(tvdi_out_file + '_line_equations.txt', 'a')
                string = out_name + '\t' + str(lin_fit[1]) + '\t' + str(lin_fit[0]) + '\t' + str(
                    ts_min) + '\t' + '0' + '\t' + str(fit_stats['r']) + '\t' + str(
                    fit_stats['dry_edge_points']) + '\t' + str(fit_stats['total_points']) + '\t' + str(
                    fit_stats['points_per_bin']) + '\t' + str(fit_stats['ndvi_range_ratio']) + '\n'
                lun.write(string)
                lun.flush()
                lun.close()
                # if the output variable is to be TVDI then calcualte it
                if output_var == cTVDI:
                    tvdi_window = (ts_window - ts_min) / (lin_fit[1] + lin_fit[0] * ndvi_window - ts_min)
                    tvdi_window[tvdi_window > 1.0] = 1.0
                    tvdi_window[tvdi_window < 0.0] = 0.0
                    output_window = tvdi_window
                # else calculate FI or EF
                else:
                    fi_min = fi_max * ((ndvi_window - ndvi_min) / (ndvi_max - ndvi_min)) ** power
                    ts_max = lin_fit[1] + lin_fit[0] * ndvi_window
                    fi_window = (ts_max - ts_window) / (ts_max - ts_min) * (fi_max - fi_min) + fi_min
                    # if delta file is present output EF, ohterwise output fi
                    if delta:
                        output_window = fi_window * delta_window
                    else:
                        output_window = fi_window
            # incorporate the tvdi_window into the appropriate place in the large tvdi matrix or just use the central pixel if moving window is used
            if roi_inf['moving_window'] == cMOVING_WINDOW:
                output[row + window_size_y / 2 + 1 - dims[0],
                       col + window_size_x / 2 + 1 - dims[1]] = output_window[window_size_y / 2 + 1,
                                                                              window_size_x / 2 + 1]
                col = col + 1;
            else:
                output[row - dims[0]:row - dims[0] + window_size_y,
                col - dims[1]:col - dims[1] + window_size_x] = output_window
                col = col + roi_inf['window_size_x']
        col = dims[0]
        # increase row by 1 if moving window otherwise by window size
        if roi_inf['moving_window'] == cMOVING_WINDOW:
            row = row + 1;
        else:
            row = row + roi_inf['window_size_y']
    return output

def calc_triangle_simple(ndvi, ts, ndvi_step, ndvi_lower_limit, ts_min_method,
                         ts_min_params, plot_out_file=None):
    ''' Estimates the dry and wet edges of the LST-VI triangle using the
    simple method for the dry edge.

    Parameters
    ----------
    ndvi, ts, ndvi_step, ndvi_lower_limit, ts_min_method,
                         ts_min_params, plot_out_dir, plot_out_file
    '''
    import numpy as np
    from numpy import amax, amin, floor, zeros, isfinite, mean, where, arange, size, logical_and
    from scipy.stats import linregress as linfit
    # Dry edge value plot
    max_ndvi = np.nanmax(ndvi)
    range_ndvi = max_ndvi - ndvi_lower_limit
    lin_fit = [0, 0]
    ts_min = 0
    fit_stats = dict()
    if range_ndvi < 2 * ndvi_step:
        print('Range in VI sample is too small, skipping Ts/VI space')
        return lin_fit, ts_min, fit_stats
    steps = int(floor(range_ndvi / ndvi_step))
    ndvi_val_old = ndvi_lower_limit
    ts_max_arr = zeros(steps)
    ts_min_arr = zeros(steps)
    ndvi_arr = zeros(steps)
    # fill into one dimensional array
    for j in arange(steps):
        ndvi_val_new = ndvi_step + ndvi_val_old
        ndvi_ind = logical_and(ndvi >= ndvi_val_old, ndvi < ndvi_val_new)
        # need at least two points in each bin
        if np.sum(ndvi_ind) >= 2:
            ts_max_arr[j] = np.nanmax(ts[ndvi_ind])
            ts_min_arr[j] = np.nanmin(ts[ndvi_ind])
        else:
            ts_max_arr[j] = 0
            ts_min_arr[j] = 0
        ndvi_arr[j] = ndvi_val_new
        ndvi_val_old = ndvi_val_new
    zeroindex = ts_max_arr != 0
    if np.sum(zeroindex) == 0: return lin_fit, ts_min, fit_stats
    ts_nozero = ts_max_arr[zeroindex]
    ndvi_nozero = ndvi_arr[zeroindex]
    # ==============================================================================
    # Cutting away points to the left of the maximum LST value and points lying lower
    # than the mean minimum temperature
    # ==============================================================================
    ts_first = ts_nozero
    ndvi_first = ndvi_nozero
    if size(ts_first) < 2: return lin_fit, ts_min, fit_stats
    # determine temporary Ts_min
    ts_min = mean(ts_min_arr[ts_min_arr > 0])
    # Determine size of reduced arrays
    array_size = int(size(ts_first))
    # Finding maximum LST cut-off
    max_plot_ind = np.argmax(ts_first)
    # Filling arrays after maximum LST value cut-off
    ts_second = ts_first[max_plot_ind:array_size]
    ndvi_second = ndvi_first[max_plot_ind:array_size]
    # Cutting off points lower than the mean minimum temperature
    ts_min_index = ts_second > ts_min
    if np.any(ts_min_index):
        ts_third = ts_second[ts_min_index]
        ndvi_third = ndvi_second[ts_min_index]
    else:
        ndvi_third = ndvi_second
        ts_third = ts_second
    # ==============================================================================
    #     Making linear  fit
    # ==============================================================================
    lin_fit = linfit(ndvi_third, ts_third)
    if not isfinite(lin_fit[0]):
        lin_fit = [0, 0]
        return lin_fit, ts_min, fit_stats
    fit_stats = calc_fit_stats(lin_fit, ts_third, ndvi_third, ts, ndvi, steps, range_ndvi)
    # calculate ts_min
    ts_min = calc_ts_min(ts_min_method, ts_min_arr, ts_min_params, lin_fit, max_ndvi)
    if plot_out_file:
        # call plot routine . Create Scatter plot
        yrange = [np.nanmin(ts[ts > 0]), np.nanmax(ts)]
        plot_pro(plot_out_file, ndvi, ts, ndvi_third, ts_third, ts_min, lin_fit, yrange)
    return lin_fit, ts_min, fit_stats

def calc_triangle_tang(ndvi, ts, ndvi_step, ndvi_lower_limit, ts_min_method,
                       ts_min_params, plot_out_file=None):
    #==============================================================================
    # for algorithm details see Tang et al., Remote Sensing of Environment 114 (2010)
    # "An application of the Ts-VI triangle method with enhaced edges determination for evapotranspiration estimation ..."
    #==============================================================================
    #import constants
    from numpy import amax,amin, floor,zeros, nan, isfinite, mean ,std, round, sqrt, where, arange,size, logical_and
    from scipy.stats import linregress as linfit
    lin_fit=[0,0]
    ts_min=0
    fit_stats=dict()
    max_ndvi = np.nanmax(ndvi)
    range_ndvi = max_ndvi - ndvi_lower_limit
    if range_ndvi < 2*ndvi_step:
        print('Range in VI sample is too small, skipping file')
        return lin_fit,ts_min,fit_stats
    intervals = int(floor(range_ndvi/ndvi_step))
    subintervals = 5
    ts_max_arr = zeros(intervals)
    ts_min_arr = zeros(intervals)
    ndvi_arr = zeros(intervals)
    ndvi_val_old = ndvi_lower_limit
    #for steps descriptions see the referenced paper, page 543
    #step (i)
    for j in arange(intervals):
        subint_max_ts = zeros(subintervals)
        for i in arange(subintervals):
            ndvi_val_new = float(ndvi_step/subintervals) + ndvi_val_old
            #step (ii)
            #need at least two points in each bin
            ndvi_ind = logical_and(ndvi >= ndvi_val_old,ndvi < ndvi_val_new)
            if np.sum(ndvi_ind)>2:
                subint_max_ts[i] = np.nanmax(ts[ndvi_ind])
            else:
                subint_max_ts[i] = nan
            ndvi_val_old = ndvi_val_new
        #step (iii)
        not_NaN = isfinite(subint_max_ts)
        if np.sum(not_NaN)==0:
            mean_max_ts = nan
            std_max_ts = nan
        else:
            subint_max_ts = subint_max_ts[not_NaN]
            mean_max_ts = mean(subint_max_ts)
            std_max_ts = std(subint_max_ts)
        while True:
            if not isfinite(std_max_ts): break
            #step (iv)
            num_el = size(subint_max_ts)
            subint_max_ts = subint_max_ts[subint_max_ts >= mean_max_ts - std_max_ts]
            if size(subint_max_ts) == num_el: break
            #step (v)
            mean_max_ts = mean(subint_max_ts)
            std_max_ts = std(subint_max_ts)
            #step (vi)
            if size(subint_max_ts) <= round(subintervals/2.0) or std_max_ts <= 4.0: break
        #step (vii)
        if isfinite(mean_max_ts): ts_max_arr[j] = mean_max_ts #otherwise ts_max_arr[j] stays as 0
        ndvi_arr[j] = ndvi_lower_limit + j*ndvi_step
        ndvi_val_old = ndvi_lower_limit + (j+1.0)*ndvi_step
        #finding ts_min_arr is not in the Tang algorithm
        #only do this if ts_min is set using the mean or median methods
        if ts_min_method == cMEAN or ts_min_method == cMEDIAN:
            ndvi_ind=np.logical_and(ndvi >= ndvi_val_old-ndvi_step, ndvi < ndvi_val_old)
            if np.any(ndvi_ind):
                ts_min_arr[j] = np.nanmin(ts[ndvi_ind])
    #the following is not part of the original Tang algorithm
    #remove all the invalid values before fitting the line
    valid_pix = ts_max_arr > 0
    # need at least two pixels to fit a straight line
    if size(valid_pix) <= 1: return lin_fit,ts_min,fit_stats
    ts_max_arr = ts_max_arr[valid_pix]
    ndvi_arr = ndvi_arr[valid_pix]
    #Finding maximum LST cut-off
    max_plot_ind = np.argmax(ts_max_arr)
    #Filling arrays after maximum LST value cut-off
    ts_max_arr = ts_max_arr[max_plot_ind:]
    ndvi_arr = ndvi_arr[max_plot_ind:]
    #need at least to pixels to fit a straight line
    valid_pix = ts_max_arr > 0
    if np.sum(valid_pix) <= 1 : return lin_fit,ts_min,fit_stats
    while True:
        #step (viii)
        lin_fit= linfit(ndvi_arr,ts_max_arr)
        fitted_ts = lin_fit[1] + lin_fit[0]*ndvi_arr
        rmse = sqrt(sum((ts_max_arr - fitted_ts)**2)/size(ndvi_arr))
        #step (ix)
        good_points = where(abs(fitted_ts-ts_max_arr) <= 2*rmse)
        if size(good_points) == size(ndvi_arr) : break
        ndvi_arr = ndvi_arr[good_points]
        ts_max_arr = ts_max_arr[good_points]
        if size(ndvi_arr) < 5: break
    #step (x)
    lin_fit= linfit(ndvi_arr,ts_max_arr)
    fit_stats = calc_fit_stats(lin_fit, ts_max_arr, ndvi_arr, ts, ndvi, intervals, range_ndvi)
    #calculate ts_min
    ts_min = calc_ts_min(ts_min_method, ts_min_arr, ts_min_params, lin_fit, max_ndvi)
    #set output name
    if plot_out_file:
        #call plot routine . Create Scatter plot
        ts_nonzero = ts > 0
        yrange = [np.nanmin(ts[ts_nonzero]), np.nanmax(ts)]
        plot_pro(plot_out_file, ndvi, ts, ndvi_arr, ts_max_arr, ts_min, lin_fit, yrange)
    return lin_fit, ts_min, fit_stats
def calc_ts_min(ts_min_method, ts_min_arr, ts_min_params, lin_fit, max_ndvi):
    #import constants
    from numpy import where, mean, median, amin, size
    a = lin_fit[1]
    b = lin_fit[0]
    if ts_min_method == cVAR_MAX_NDVI:
        return a + b*np.nanmin([ts_min_params, max_ndvi])
    if ts_min_method == cCONST_TS_MIN :
        return ts_min_params
    if ts_min_method == cDAILY_TS_MIN :
        return ts_min_params
    if ts_min_method == cMEAN or ts_min_method == cMEDIAN:
        #take the mean or median of the min Ts values of specified number of largest NDVIs
        point_num = int(ts_min_params)
        temp = where(ts_min_arr > 0)[0]
        temp_size=size(temp)
        if temp_size<=0: return 0
        if temp_size > point_num : temp = temp[temp_size-point_num:temp_size]
        if ts_min_method == cMEAN : return mean(ts_min_arr[temp])
        if ts_min_method == cMEDIAN : return median(ts_min_arr[temp])
def calc_fit_stats(lin_fit, ts_select, ndvi_select, ts_all, ndvi_all, bins, range_ndvi):
    from numpy import where, size, amax, amin
    fit_stats=dict()
    fit_stats['r'] =float(lin_fit[2])
    fit_stats['dry_edge_points'] = float(size(ts_select))
    fit_stats['total_points'] = float(size(where((ndvi_all != 0) & (ts_all != 0))))
    fit_stats['points_per_bin'] = float(fit_stats['total_points']/bins)
    fit_stats['ndvi_range_ratio'] = float((np.nanmax(ndvi_select)-np.nanmin(ndvi_select))/range_ndvi)
    print('Linear fit parameters, y=A+Bx')
    print('A= '+str(lin_fit[1])+ ' B= '+str(lin_fit[0])+ ' r= '+str(fit_stats['r'])
          + ' dep= '+str(fit_stats['dry_edge_points'])+ ' tp= ' +str(fit_stats['total_points'])
          + ' ppb= '+str(fit_stats['points_per_bin'])+ ' nrr= '+str(fit_stats['ndvi_range_ratio']))
    return fit_stats

def plot_pro(out_name, x1array, y1array, x2array, y2array, ts_min, linefit, yrange):
    import matplotlib.pylab as plt
    from numpy import amin, arange,zeros, size
    # Turn interactive plotting off
    plt.ioff()
    fig=plt.figure()
    #set up the plot and do a scatterplot of x1array, y1array
    # do scatterplotts of the other arrays if they are present
    plt.plot(x1array,y1array, 'o',markersize=1, alpha=0.5, markeredgecolor='black', color='white')
    if size(x2array) >0 and size(y2array) > 0: plt.plot(x2array, y2array,'^', markersize=6, color='red')
    plt.xlabel('NDVI')
    plt.ylabel('TS')
    plt.title('Triangle')
    plt.ylim(yrange)
    # plot Dry edge
    x_dummy=arange(20)/19.0
    # scale the line to go through the whole range of x
    x_min=np.nanmin(x1array)
    x_dummy = x_dummy*(1. - x_min) + x_min
    y = linefit[1] + x_dummy*linefit[0]
    plt.plot(x_dummy,y,'k-', lw=2)
    if linefit[0]>0:
        string='y='+str(linefit[1])+ '+' +str(linefit[0])+'*x'
    else:
        string='y='+str(linefit[1])+ str(linefit[0])+'*x'
    label_pos_x=0.2
    label_pos_y=yrange[0]+0.9*(yrange[1]-yrange[0])
    # plot ts_min
    x_dummy=[x_min, 1]
    if ts_min > 0:
        plt.plot(x_dummy, zeros(2)+ts_min,'k-', lw=2)
        string=string+'\n TS min = '+str(ts_min)
    plt.text(label_pos_x,label_pos_y,string, color='red')
    # Save the file
    plt.savefig(out_name)
    plt.close(fig)