# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 15:23:07 2024

@author: Sebastian Lehmler
"""
import os

from core_functions import open_raster_get_FORCE_bounds, fit_box, export_raster, transform_to_gv, transform_to_canopy_cover 
from core_functions import sum_aggregate_raster_10m_force, average_aggregate_raster_10m_force, std_aggregate_raster_10m_force, mask_nodata

########################################## Configuration
# list of city names
city = [
        'Berlin'
        ]

# list of corresponding years       
year = [
        '2024'
        ]

# list of corresponding vegetation height raster paths       
veg_height_raster = [
                    r'PATH\TO\veg_height.tif'
                     ]

# Switches - which functionalities should be run
veg_h_aggregate = True #set veg_h values < 0 to 0 and and aggregate to 10m average height
veg_h_aggregate_SD = True
gv = True # Set to True if you want to perform the transform_to_gv function and the associated aggregation to 10m - all negative values will be set to 0
canopy = True # Set to True if you want to perform the transform_to_canopy_cover function and the associated aggregation to 10m
only_upscale = False  # Set to True if you only want to apply upscaling without transformation


# advanced config - location of FORCE 10m raster for germany, canopy threshold and xarray chunksize
force_germany_bounds = (4016026.363042, 2654919.607965, 4676026.363042001, 3554919.607965) #bound of FORCE GRID for Germany for EPSG:3035 - needs to be expanded if working outside germany 
canopy_threshold = 2.5 #vegetation height cutoff for being considered as canopy
chunksize = 512 #chunksize for geowombat export - may need to be decreased for small images to avoid RasterBlockError

######################################### loop over Datasets in Main Function

def main():
    """Main function."""
    # Iterate over each dataset using zip()
    for city_i, year_i, veg_height_raster_i in zip(city, year, veg_height_raster):
        # Get result directory for current dataset
        result_dir_i = os.path.dirname(veg_height_raster_i)

        # Generate result paths for current dataset
        result_veg_h_noneg_i = os.path.join(result_dir_i, f"{city_i}_{year_i}_vegetation_height_noneg_50cm.tif")
        result_gv_i = os.path.join(result_dir_i, f"{city_i}_{year_i}_green_volume_50cm.tif")
        result_canopy_i = os.path.join(result_dir_i, f"{city_i}_{year_i}_canopy_cover_binary_50cm.tif")
        result_veg_h_10m_i = os.path.join(result_dir_i, f"{city_i}_{year_i}_vegetation_height_average_10m_FORCE.tif")
        result_veg_h_10m_std_i = os.path.join(result_dir_i, f"delete_{city_i}_{year_i}_vegetation_height_standard_deviation_10m_FORCE_unmasked.tif")
        result_veg_h_10m_std_masked_i = os.path.join(result_dir_i, f"{city_i}_{year_i}_vegetation_height_standard_deviation_10m_FORCE.tif")
        result_canopy_10m_i = os.path.join(result_dir_i, f"{city_i}_{year_i}_canopy_cover_fraction_10m_FORCE.tif")
        result_gv_10m_i = os.path.join(result_dir_i, f"{city_i}_{year_i}_green_volume_sum_10m_FORCE.tif")

        # Open vegetation height raster
        veg_h, veg_h_attrs, veg_h_bounds, veg_h_nodata, veg_h_res, processing_time = open_raster_get_FORCE_bounds(veg_height_raster_i, chunksize)

        # 10m aggregation bounding box
        adjusted_bbox = fit_box(force_germany_bounds, veg_h_bounds)
        
        
        # Export vegetation_height raster with values below zero set to zero
        if veg_h_aggregate:
            # export adjusted veg_h raster 
            print(f'exporting veg_h_no_negative raster for {city_i}...')
            export_raster(veg_h, veg_h_attrs, result_veg_h_noneg_i, chunksize)
            
            # Open vegetation height raster in 10m resolution with adjusted bounding box and export
            veg_h_aggregated, veg_h_aggregated_attrs = average_aggregate_raster_10m_force(result_veg_h_noneg_i, adjusted_bbox)
            print(f'exporting aggregated veg_h raster for {city_i}...')
            export_raster(veg_h_aggregated, veg_h_aggregated_attrs, result_veg_h_10m_i, chunksize)
        
        if veg_h_aggregate_SD:
            
            # Open vegetation height raster in 10m resolution with adjusted bounding box and export
            print('performing SD aggregation...')
            veg_h_aggregated_std_unmasked, veg_h_aggregated_attrs = std_aggregate_raster_10m_force(result_veg_h_noneg_i,adjusted_bbox, veg_h_res)
            print(f'exporting Standard Deviation for veg_h raster for {city_i}...')
            export_raster(veg_h_aggregated_std_unmasked, veg_h_aggregated_attrs, result_veg_h_10m_std_i, chunksize)
            
            
            #mask and export SD raster 
            print(f'masking Standard Deviation for veg_h raster for {city_i}...')
            veg_h_aggregated_std_masked, veg_h_aggregated_attrs = mask_nodata(result_veg_h_10m_std_i, result_veg_h_10m_i)
            export_raster(veg_h_aggregated_std_masked, veg_h_aggregated_attrs, result_veg_h_10m_std_masked_i, chunksize)
            
                

        # Check if gv switch is on
        if gv:
            # Transform vegetation height to GV
            gv_data = transform_to_gv(veg_h)

            # Export GV raster
            print(f'exporting gv raster for {city_i}...')
            export_raster(gv_data, veg_h_attrs, result_gv_i, chunksize)

            # Open gv raster in 10m resolution with sum aggregation and adjusted bounding box and export
            gv_aggregated, gv_aggregated_attrs = sum_aggregate_raster_10m_force(result_gv_i, adjusted_bbox, veg_h_res)
            print(f'exporting aggregated gv raster for {city_i}...')
            export_raster(gv_aggregated, gv_aggregated_attrs, result_gv_10m_i, chunksize)

        # Check if canopy switch is on
        if canopy:
            # Transform vegetation height to canopy cover
            canopy_binary = transform_to_canopy_cover(veg_h, canopy_threshold)

            # Export Canopy cover raster
            print(f'exporting canopy raster for {city_i}...')
            export_raster(canopy_binary, veg_h_attrs, result_canopy_i, chunksize)

            # Open canopy raster in 10m resolution with adjusted bounding box and export
            canopy_aggregated, canopy_aggregated_attrs = average_aggregate_raster_10m_force(result_canopy_i, adjusted_bbox)
            print(f'exporting aggregated canopy raster for {city_i}...')
            export_raster(canopy_aggregated, canopy_aggregated_attrs, result_canopy_10m_i, chunksize)

        if only_upscale:
            print(f'Only upscaling mode enabled for {city_i}...')
            # Open canopy raster in 10m resolution with adjusted bounding box and export
            canopy_aggregated, canopy_aggregated_attrs = average_aggregate_raster_10m_force(result_canopy_i, adjusted_bbox)
            print(f'exporting aggregated canopy raster for {city_i}...')
            export_raster(canopy_aggregated, canopy_aggregated_attrs, result_canopy_10m_i, chunksize)



if __name__ == '__main__':
    main()
