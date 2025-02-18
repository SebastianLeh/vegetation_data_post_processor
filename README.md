
# Vegetation Data Post-Processor

## Description

A toolkit for transforming, aggregating, and upscaling various types of vegetation data. Currently used to transform CNN 
derived canopy height predictions into green volume and upscale them to [FORCE Datacube Grid](https://force-eo.readthedocs.io/en/latest/howto/datacube.html) for further use in satellite deep learning models.

## Key Functionalities:

### 1. High-Resolution Data Transformation 
- **Input**: High-resolution canopy height data from CNN predictions
- **Outputs**: Green volume data, binary canopy data (>2.5m = canopy), original canopy height data
- **Code**: canopy_height_upscaling.py


### 2. Data Aggregation to FORCE Grid 
- **Input**: Three raster products from step 1
- **Output**: For each input, mean, SD, and sum per 10m pixel aligned to FORCE Grid
- **Code**: canopy_height_upscaling.py

### 3. Flexible Upscaling (e.g., Damage Analysis)
- **Input**: Vector data on damaged vegetation
- **Output**: Raster data of vegetation damage at 10m resolution, FORCE grid-aligned
- **Code**: vector_to_raster_damage_analysis.py

This versatile toolkit handles data type transformations and resolution changes for various ecological and forestry applications.

## Green Volume calculation
The following table describes how green volume data is derived from the canopy height prediction. The constants for sealed surface, grassland and cropland in m³/m² are the same in m for the CNN based vegetation height predictions.
The constants of 10 % and 25 % substracted for high vegetation represent estimations of missing vegetation volume around tree stems. 

| Land Cover Type | Green Volume |
|-----------------|--------------|
| Sealed surface & Water | 0 m³/m² |
| Grassland | 0.5 m³/m² |
| Cropland | 1.0 m³/m² |
| Shrubs (< 5m) | Pixelsize x canopy height |
| Shrubs and Trees (5 - 7m) | Pixelsize x canopy height - 10% |
| Trees (> 7m) | Pixelsize x canopy height - 25% |

## Getting Started

### Dependencies

* GDAL, Geopandas, Geowombat, xarray
* Anaconda [download](https://www.anaconda.com/download) 
* developed on Windows

### Installation
* The Python GDAL package version must match the GDAL binaries version. Be sure to use the same version printed from:
```python I'm A tab
gdalinfo --version
```
* Create a virtual Conda environment with the required Python version and requirements file:
```python I'm A tab
conda create --name gwenv python=3.8
conda activate gwenv
conda install -c conda-forge gdal==3.4.3
conda config --env --add channels conda-forge
conda config --env --set channel_priority strict
cd ../environment
pip install -r requirements.txt
```
## Authors

 
 - [Sebastian Lehmler](https://github.com/SebastianLeh)- Core functionality and canopy height upscaling use case
 - [Shadi Ghantous](https://github.com/LUP-LuftbildUmweltPlanung) - Code improvements and damage analysis use case


## Reference

* [GeoWombat documentation](https://geowombat.readthedocs.io/en/latest/).
* [FORCE Framework](https://force-eo.readthedocs.io/en/latest/index.html)
* [Geopandas](https://geopandas.org/en/stable/)
* [GDAL](https://gdal.org/index.html)




