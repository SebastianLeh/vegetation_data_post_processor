
# FORCE_GRID_SPATIAL_DATA_AGGREGATION



## Description

"The FORCE_GRID_SPATIAL_DATA_AGGREGATION contains functions to aggregate geospatial data to a lower resolution Grid. We use the functionalities to aggregate raster or vector data of different environmental indicators to the [FORCE Grid](https://force-eo.readthedocs.io/en/latest/howto/datacube.html) at 10m resolution to facilitate usage in Sentinel-2 based machine learning models. Aggregation functions include e.g. Mean, Standard Deviation and Sum. 

Additionally, some functionality for high resolution raster transformation is included. For our projects, we for example transform high resolution canopy height rasters to green volume (applying some rules based transformations) and into binary canopy rasters (using a specified threshold value). These are then upscaled to the 10*10m FORCE Grid using various aggregation statistics such as mean, sum and standard deviation." 


## Getting Started

### Dependencies

* GDAL, Geopandas, Geowombat, xarray
* Anaconda [download](https://www.anaconda.com/download) 
* developed on Windows

### Installation

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




