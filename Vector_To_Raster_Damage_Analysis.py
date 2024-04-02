# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 09:00:44 2024

@author: LUP
"""

# Import Library
from core_functions import polygons_to_raster,open_raster_get_FORCE_bounds, export_raster, average_aggregate_raster_10m_force, fit_box
import os


# Please Define the variables to rasterize your vector

# insert the path for the vector file
shapefile_path = "path/to/shapefile"

# define the place to save the output raster
raster_out_path = "path/output/ratser"

# Define the resolution of the raster file, maybe need to be increased for large images to avoid Error
raster_resolution = 1.0  # Resolution of the raster


# Define the custom boundary to change the default extent (Optional)
custom_boundary = None

################## Rasterize the Vector
raster = polygons_to_raster(shapefile_path, raster_out_path, raster_resolution, custom_boundary)

######################## Calculate the Average in every 10 m pixel
# vegetation damaged rasters variables:       
veg_damaged_raster = raster_out_path

# advanced config - location of FORCE 10m raster for germany, canopy threshold and xarray chunksize
force_germany_bounds = (4016026.363042, 2654919.607965, 4676026.363042001, 3554919.607965) #bound of FORCE GRID for Germany for EPSG:3035 - needs to be expanded if working outside germany 
chunksize = 64 #chunksize for geowombat export - may need to be decreased for small images to avoid RasterBlockError

# The main function to claculate the Average value 
def main():
    """Main function."""

    # Get result directory for current dataset
    result_dir_i = os.path.dirname(raster_out_path)

    # Generate result paths for current dataset
    result_veg_d_noneg_i = os.path.join(result_dir_i, r"_vegetation_damaged_noneg_50cm.tif")
    result_veg_d_10m_i = os.path.join(result_dir_i, r"_vegetation_damaged_average_10m_FORCE.tif")


    # Open vegetation height raster
    veg_d, veg_d_attrs, veg_d_bounds, veg_d_nodata, veg_d_res = open_raster_get_FORCE_bounds(veg_damaged_raster)

    # 10m aggregation bounding box
    adjusted_bbox = fit_box(force_germany_bounds, veg_d_bounds)
        

    # export adjusted veg_h raster 
    print(f'exporting veg_h_no_negative raster ...')
    export_raster(veg_d, veg_d_attrs, result_veg_d_noneg_i)
            
    # Open the raster in 10m resolution with adjusted bounding box and export
    veg_d_aggregated, veg_d_aggregated_attrs = average_aggregate_raster_10m_force(result_veg_d_noneg_i, adjusted_bbox)
    print(f'exporting aggregated veg_h raster...')
    export_raster(veg_d_aggregated, veg_d_aggregated_attrs, result_veg_d_10m_i)


if __name__ == '__main__':
    main()
