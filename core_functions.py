# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 15:01:19 2024

@author: Sebastian
"""

import rasterio
from rasterio import features
from rasterio.transform import from_bounds
import geopandas as gpd
import geowombat as gw
import xarray as xr
import os
from osgeo import gdal
import time

######################################### Polyogon to raster functionality
def polygons_to_raster(shapefile_path, raster_out_path, raster_resolution, nodata_value=0, custom_boundary= None):
    """
    Convert polygons from a vector to a raster and save the raster to the specified path.

    Args:
    shapefile_path (str): Path to the input shapefile.
    raster_out_path (str): Path to save the output raster file.
    raster_resolution (float): Resolution of the output raster in the same units as the shapefile's coordinate system.
    nodata_value (int, optional): Value to assign to cells outside the polygons. Defaults to 0.
    custom_boundary (tuple, optional): Custom boundary extent (minx, miny, maxx, maxy). If not provided, the extent is derived from the shapefile.

    Returns:
    None

    Raises:
    ValueError: Raised if an error occurs during rasterization. Consider increasing raster_resolution.

    Example:
    # Example usage
    shapefile_path = "Path/to/Input_shapefile"
    raster_out_path = "Path/to/output_tiff_file"
    raster_resolution = 1.0  # Default resolution of the raster

    After successful execution, a raster TIFF file is saved at the specified raster_out_path.
    """
    # Read the shapefile using geopandas:
    gdf = gpd.read_file(shapefile_path)
    
    # calculate bounds based on custom_boundary if provided
    if custom_boundary:
        minx, miny, maxx, maxy = custom_boundary
    else:
        minx, miny, maxx, maxy = gdf.total_bounds
        
    # calculate the width and hight of the raster:
    width = int((maxx - minx) / raster_resolution)
    height = int((maxy - miny) / raster_resolution)
    
    # Define the transformation:
    transform = from_bounds(minx, miny, maxx, maxy, width, height)
    
    # create an empty raster array
    raster_array = features.rasterize(
        [(geometry, 1) for geometry in gdf.geometry],
        out_shape = (height, width),
        transform = transform,
        fill = nodata_value,
        all_touched=True,
        dtype=rasterio.uint8
    )
    
    # Write the raster array to a GeoTIFF file:
    with rasterio.open(
        raster_out_path,
        'w',
        driver= 'GTiff',
        width = width,
        height= height,
        count= 1,
        dtype=rasterio.uint8,
        crs = gdf.crs,
        transform=transform,
    
    ) as dst:
        dst.write(raster_array, 1)
        
        
######################################### High resolution Raster Processing Functions
def open_raster_get_FORCE_bounds(raster, chunksize):
    """
    Opens a raster file, processes values < zero to zero, and returns its data and attributes along with bounds.

    Args:
        raster (str): Path to the raster file.
        chunksize (int): Chunk size for raster processing.

    Returns:
        tuple: A tuple containing the raster data, its attributes,
            its bounds, its nodata value, and its resolution.
    """
    with gw.open(raster, chunks=chunksize) as veg_h:
        veg_h = veg_h.gw.set_nodata(src_nodata=veg_h.attrs['_FillValue'], dst_nodata=-9999) #get original nodata value and replace with new one
        veg_h_attrs = veg_h.attrs #get attributes
        veg_h_nodata = veg_h.attrs['_FillValue'] #get new nodata value
        veg_h_res = veg_h.attrs['res'][0] #get resolution
        veg_h = xr.where((veg_h < 0) & (veg_h != veg_h_nodata), 0, veg_h).assign_attrs(veg_h_attrs) #set veg_h raster < zero to zero and reassign attrs

    # Get bounds of vegetation height raster in FORCE CRS
    with gw.config.update(ref_crs='EPSG:3035'):
        with gw.open(raster) as veg_h_rep:
            veg_h_rep_bounds = veg_h_rep.gw.bounds

    return veg_h, veg_h_attrs, veg_h_rep_bounds, veg_h_nodata, veg_h_res


def transform_to_gv(veg_h):
    """
    Transforms a raster to a green volume (GV) raster based on predefined thresholds.

    Args:
        raster (xarray.DataArray): The input raster data.

    Returns:
        xarray.DataArray: The transformed green volume data.
    """
    veg_h_nodata = veg_h.attrs['_FillValue']
    multip_fact = (veg_h.res[1] ** 2)
    gv_low = xr.where(veg_h < 5, multip_fact * veg_h, 0) #low vegetation 
    gv_med = xr.where((veg_h >= 5) & (veg_h <= 9), 0.9 * (multip_fact * veg_h), 0) #medium vegetation
    gv_high = xr.where(veg_h > 9, 0.75 * (multip_fact * veg_h), 0) #high vegetation
    gv = gv_low + gv_med + gv_high #add 
    gv = xr.where(veg_h == veg_h_nodata, veg_h_nodata, gv) #re-set nodata value
    
    return gv


def transform_to_canopy_cover(veg_h, canopy_threshold):
    """
    Transforms a raster to a binary canopy cover raster based on a specified threshold.

    Args:
        raster (xarray.DataArray): The input raster data.
        canopy_threshold (float): The threshold value for canopy cover.

    Returns:
        xarray.DataArray: The binary canopy cover data.
    """
    veg_h_nodata = veg_h.attrs['_FillValue']
    canopy_binary = xr.where(veg_h > canopy_threshold, 1, 0)
    canopy_binary = xr.where(veg_h == veg_h_nodata, veg_h_nodata, canopy_binary) #re-set nodata value
    return canopy_binary


######################################### Bounding Box and Aggregation Functions
    
def fit_box(box1, box2):
    """
    Adjusts the first bounding box until it fits inside the second bounding box.

    Args:
        box1 (tuple): Tuple of coordinates (left, bottom, right, top) for the first bounding box.
        box2 (tuple): Tuple of coordinates (left, bottom, right, top) for the second bounding box.

    Returns:
        tuple: The adjusted bounding box as a tuple of coordinates (left, bottom, right, top).
    """
    step = 10
    left = box1[0]
    bottom = box1[1]
    right = box1[2]
    top = box1[3]

    while left < box2[0]:
        left += step

    while bottom < box2[1]:
        bottom += step

    while right > box2[2]:
        right -= step

    while top > box2[3]:
        top -= step

    return (left, bottom, right, top)


def average_aggregate_raster_10m_force(input_raster, adjusted_bbox):
    """Aggregates a raster to 10m resolution using average resampling.

    Args:
        input_raster (str): Path to the input raster file.
        adjusted_bbox (tuple): The adjusted bounding box to use for aggregation.

    Returns:
        tuple: A tuple containing the aggregated raster data and its attributes.
    """
    with gw.config.update(ref_crs='EPSG:3035', ref_res=(10, 10), ref_bounds=adjusted_bbox, nodata=-9999):
        with gw.open(input_raster, resampling='average') as input_aggregated:
            input_aggregated_attrs = input_aggregated.attrs
            return input_aggregated, input_aggregated_attrs


def sum_aggregate_raster_10m_force(input_raster, adjusted_bbox, veg_h_res):
    """Aggregates a raster to 10m resolution using sum aggregation.

    Args:
        input_raster (str): Path to the input raster file.
        adjusted_bbox (tuple): The adjusted bounding box to use for aggregation.
        veg_h_res (float): The resolution of the vegetation height raster.

    Returns:
        tuple: A tuple containing the aggregated raster data and its attributes.
    """
    with gw.config.update(ref_crs='EPSG:3035', ref_res=(10, 10), ref_bounds=adjusted_bbox, nodata=-9999):
        with gw.open(input_raster, resampling='average') as input_aggregated:
            input_aggregated_attrs = input_aggregated.attrs
            # Get n of pixels within 10m pixel to go from average to sum
            # (Area of one pixel of raster a / Area of one pixel of raster b)
            conversion_factor = ((10**2)/(veg_h_res**2))
            # Calculate gv sum in 10m pixel where there is no nodata
            input_aggregated_sum = xr.where(input_aggregated != -9999,
                                            input_aggregated * conversion_factor,
                                            input_aggregated)
            return input_aggregated_sum, input_aggregated_attrs
        
def std_aggregate_raster_10m_force(input_raster, adjusted_bbox, veg_h_res):
    """Aggregates a raster to 10m resolution using standard deviation.

    Args:
        input_raster (str): Path to the input raster file.
        adjusted_bbox (tuple): The adjusted bounding box to use for aggregation.
        veg_h_res (float): The resolution of the vegetation height raster.

    Returns:
        tuple: A tuple containing the aggregated raster data and its attributes.
    """
    
    with gw.config.update(ref_crs='EPSG:3035', ref_bounds=adjusted_bbox, ref_res=10,nodata=-9999):
        with gw.open(input_raster) as dummy_raster:
            #get attributes to be used later for exporting coarsened raster 
            input_dummy_attrs = dummy_raster.attrs

    with gw.config.update(ref_crs='EPSG:3035', ref_bounds=adjusted_bbox, ref_res=0.5, nodata=0):
        with gw.open(input_raster) as input_raster:

            # Get the number of pixels to group by in each dimension
            factor = round(10 / veg_h_res)
        
    # Coarsen Array along x and y dimension and take the SD for each coarsened Pixel - 
    input_std_agg = input_raster.coarsen(x=factor, y=factor, boundary='pad').std().assign_attrs(input_dummy_attrs) #.squeeze('time')
            
            # res set nodata using numpy 
            #input_std_agg_adj = xr.where(dummy_raster.values == -9999, -9999, input_std_agg.values) #fix dimension error
            
            #input_std_agg_adj = xr.DataArray(input_std_agg_adj, coords=dummy_raster.coords, dims=dummy_raster.dims, attrs=dummy_raster.attrs)
            
    return input_std_agg, input_dummy_attrs


def mask_nodata(input_raster, mask_raster):
    """
    Opens two rasters with the same dimensions and sets the input raster nodata where the mask raster is nodata.

    Parameters
    ----------
    input_raster : str
        The path to the input raster file.
    mask_raster : str
        The path to the mask raster file.

    Returns
    -------
    xarray.DataArray
        The masked input raster as an xarray.DataArray object.
    """
    # update the configuration using the mask raster as the reference
    with gw.config.update(ref_image  =  mask_raster):
        # open both rasters using the same configuration
        with gw.open(input_raster) as input_ras:
            with gw.open(mask_raster) as mask_ras:
                #get attributes of input raster for export later 
                export_attrs = mask_ras.attrs
                # get the nodata value from the mask raster
                mask_nodata = mask_ras.attrs['_FillValue']
                
    # use xarray.where to mask the input raster where the mask raster is nodata
    masked_input_ras = xr.where(mask_ras == mask_nodata, mask_nodata, input_ras) #re-set nodata value
    

    # return the masked input raster
    return masked_input_ras, export_attrs
        
        

######################################### Raster Export Function
def export_raster(export_raster, attrs, dest_path, chunksize):
    """Exports a raster to a file.

    Args:
        export_raster (xarray.DataArray): The raster data to export.
        attrs (dict): The attributes of the raster data.
        dest_path (str): The path to the destination file.
    """
    start_time = time.time()
    
    export_raster = export_raster.assign_attrs(attrs).astype('float32').chunk(chunks=chunksize)

    export_raster.gw.save(dest_path, 
                          num_workers = 16,
                          nodata = -9999)
    
    # Compress the exported raster using gdal
    print('Compressing raster...')
    compression_start_time = time.time()
    
    ds = gdal.Open(dest_path, gdal.GA_Update)
    filename, file_extension = os.path.splitext(dest_path)
    temp_path = filename + '_temp' + file_extension
    gdal.Translate(temp_path, ds, creationOptions=['COMPRESS=LZW', 'BIGTIFF=YES'])
    ds = None
    
    # Replace original raster with compressed raster
    os.replace(temp_path, dest_path)
    
    compression_end_time = time.time()
    compression_elapsed_time = compression_end_time - compression_start_time
    print(f'Finished compressing raster in {compression_elapsed_time:.2f} seconds ({compression_elapsed_time/60:.2f} minutes).')
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f'Finished exporting raster in {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes).')