#!/usr/bin/env python3

import xarray as xr
import rioxarray
import geopandas as gpd
import numpy as np
from pathlib import Path
import argparse
from enum import Enum

class AggregationMethod(Enum):
    MEAN = 'mean'
    ANNUAL_MEAN = 'annual_mean'
    ANNUAL_SUM_MEAN = 'annual_sum_mean'  # For precipitation: sum within years, then mean across years

def get_time_grouping(da, time_dim):
    """
    Extract year from time coordinate for grouping
    
    Parameters:
    -----------
    da : xarray.DataArray
        Input data array
    time_dim : str
        Name of time dimension
    
    Returns:
    --------
    xarray.DataArray.GroupBy
        DataArray grouped by year
    """
    # xarray automatically decodes CF-compliant time coordinates to datetime64
    try:
        return da.groupby(f"{time_dim}.year")
    except AttributeError:
        raise ValueError(
            f"Could not extract years from time coordinate. "
            f"Ensure the NetCDF file has CF-compliant time values and encoding."
        )

def temporal_aggregate(da, time_dim, method: AggregationMethod):
    """
    Aggregate data temporally according to specified method
    
    Parameters:
    -----------
    da : xarray.DataArray
        Input data array
    time_dim : str
        Name of time dimension
    method : AggregationMethod
        Method to use for temporal aggregation
    
    Returns:
    --------
    xarray.DataArray
        Temporally aggregated data
    """
    if method == AggregationMethod.MEAN:
        return da.mean(dim=time_dim)
    
    # Group by year
    grouped = get_time_grouping(da, time_dim)
    
    if method == AggregationMethod.ANNUAL_MEAN:
        # Take mean within each year, then mean across years
        # This avoids bias from leap years in daily data
        return grouped.mean().mean(dim='year')
    
    elif method == AggregationMethod.ANNUAL_SUM_MEAN:
        # Sum within each year, then mean across years
        return grouped.sum().mean(dim='year')
    
    raise ValueError(f"Unknown aggregation method: {method}")

def process_netcdf(netcdf_path, shapefile_path, output_path, agg_method: AggregationMethod):
    """
    Process a 3D NetCDF file by cropping it to a shapefile extent and computing temporal aggregation.
    
    Parameters:
    -----------
    netcdf_path : str
        Path to input NetCDF file
    shapefile_path : str
        Path to shapefile for cropping
    output_path : str
        Path for output GeoTIFF file
    agg_method : AggregationMethod
        Method to use for temporal aggregation
    """
    # Read the NetCDF file - xarray will automatically decode time coordinates
    ds = xr.open_dataset(netcdf_path)
    
    # Get the single data variable (assuming there's only one)
    data_var = list(ds.data_vars)[0]
    da = ds[data_var]
    
    # Ensure the data has spatial coordinates
    if not all(coord in da.coords for coord in ['latitude', 'longitude']):
        # Try alternate common coordinate names
        coord_map = {
            'lat': 'latitude',
            'lon': 'longitude',
            'y': 'latitude',
            'x': 'longitude'
        }
        for old, new in coord_map.items():
            if old in da.coords:
                da = da.rename({old: new})
    
    # Add CRS if not present (assuming WGS84)
    da.rio.write_crs("EPSG:4326", inplace=True)
    
    # Read and reproject shapefile
    shape = gpd.read_file(shapefile_path)
    if shape.crs is None:
        shape.set_crs("EPSG:4326", inplace=True)
    elif shape.crs != "EPSG:4326":
        shape = shape.to_crs("EPSG:4326")
    
    # Clip data to shapefile extent
    clipped = da.rio.clip(shape.geometry, shape.crs)
    
    # Find the time dimension dynamically
    time_dim = None
    for dim in clipped.dims:
        if dim not in ['latitude', 'longitude']:
            time_dim = dim
            break
    
    if time_dim is None:
        raise ValueError("Could not find time dimension in the dataset")
    
    # Calculate temporal aggregation
    result = temporal_aggregate(clipped, time_dim, agg_method)
    
    # Save as GeoTIFF
    result.rio.to_raster(output_path)

def main():
    parser = argparse.ArgumentParser(description='Process NetCDF file with shapefile masking and temporal averaging')
    parser.add_argument('netcdf_file', help='Input NetCDF file path')
    parser.add_argument('shapefile', help='Input shapefile path')
    parser.add_argument('output_file', help='Output GeoTIFF file path')
    parser.add_argument('--method', choices=[m.value for m in AggregationMethod], 
                      default='annual_mean', help='Temporal aggregation method')
    
    args = parser.parse_args()
    
    agg_method = AggregationMethod(args.method)
    process_netcdf(args.netcdf_file, args.shapefile, args.output_file, agg_method)

if __name__ == "__main__":
    main()
