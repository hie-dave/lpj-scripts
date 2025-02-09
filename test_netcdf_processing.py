#!/usr/bin/env python3

import numpy as np
import xarray as xr
import rioxarray
import geopandas as gpd
from shapely.geometry import box
import pandas as pd
from datetime import datetime
import os
from process_netcdf import process_netcdf, AggregationMethod

def create_test_data():
    """Create synthetic NetCDF datasets for testing"""
    
    # Create a small spatial grid (10x10 degrees over NSW)
    lats = np.linspace(-36, -32, 20)
    lons = np.linspace(145, 149, 20)
    
    # Create daily timestamps for 2 years
    dates = pd.date_range(start='2020-01-01', periods=731, freq='D')  # 2 years including leap year
    
    # Create temperature-like data
    # Annual cycle with noise and a warming trend
    time_index = np.arange(len(dates))
    annual_cycle = 15 + 10 * np.sin(2 * np.pi * time_index / 365.25)
    warming_trend = 0.002 * time_index
    
    temp_data = np.zeros((len(dates), len(lats), len(lons)))
    for i in range(len(lats)):
        lat_effect = -0.5 * (lats[i] + 34)  # Temperature decreases with latitude
        for j in range(len(lons)):
            lon_effect = 0.2 * (lons[j] - 147)  # Slight longitudinal gradient
            base_temp = annual_cycle + warming_trend + lat_effect + lon_effect
            temp_data[:, i, j] = base_temp + np.random.normal(0, 1, len(dates))
    
    # Create precipitation-like data
    # Seasonal pattern with many zeros and occasional heavy rainfall
    precip_data = np.zeros_like(temp_data)
    seasonal_prob = 0.3 + 0.2 * np.sin(2 * np.pi * time_index / 365.25)  # More rain in winter
    
    for i in range(len(lats)):
        for j in range(len(lons)):
            rain_days = np.random.random(len(dates)) < seasonal_prob
            intensity = np.random.exponential(5, len(dates))  # Some days with heavy rain
            precip_data[:, i, j] = rain_days * intensity
    
    # Create datasets with proper CF attributes
    temp_ds = xr.Dataset(
        data_vars=dict(
            temperature=(['time', 'latitude', 'longitude'], temp_data, {
                'units': 'celsius',
                'long_name': 'Air Temperature'
            })
        ),
        coords=dict(
            longitude=(['longitude'], lons, {
                'units': 'degrees_east',
                'long_name': 'Longitude'
            }),
            latitude=(['latitude'], lats, {
                'units': 'degrees_north',
                'long_name': 'Latitude'
            }),
            time=dates
        )
    )
    
    precip_ds = xr.Dataset(
        data_vars=dict(
            precipitation=(['time', 'latitude', 'longitude'], precip_data, {
                'units': 'mm/day',
                'long_name': 'Daily Precipitation'
            })
        ),
        coords=dict(
            longitude=(['longitude'], lons, {
                'units': 'degrees_east',
                'long_name': 'Longitude'
            }),
            latitude=(['latitude'], lats, {
                'units': 'degrees_north',
                'long_name': 'Latitude'
            }),
            time=dates
        )
    )
    
    # Save datasets
    temp_ds.to_netcdf('test_temperature.nc')
    precip_ds.to_netcdf('test_precipitation.nc')
    
    return 'test_temperature.nc', 'test_precipitation.nc'

def create_test_shapefile():
    """Create a test shapefile covering part of the domain"""
    # Create a smaller box within our domain
    geometry = box(146, -35, 148, -33)
    gdf = gpd.GeoDataFrame({'geometry': [geometry]}, crs="EPSG:4326")
    
    # Save to file
    gdf.to_file('test_region.shp')
    return 'test_region.shp'

def run_tests():
    """Run tests with different aggregation methods"""
    print("Creating test data...")
    temp_file, precip_file = create_test_data()
    shape_file = create_test_shapefile()
    
    print("\nTesting temperature data...")
    # Test both mean methods on temperature
    process_netcdf(temp_file, shape_file, 'temp_mean.tif', AggregationMethod.MEAN)
    process_netcdf(temp_file, shape_file, 'temp_annual_mean.tif', AggregationMethod.ANNUAL_MEAN)
    
    # Load and compare results
    temp_mean = rioxarray.open_rasterio('temp_mean.tif')
    temp_annual = rioxarray.open_rasterio('temp_annual_mean.tif')
    print(f"Temperature difference (mean - annual_mean): {(temp_mean - temp_annual).mean().values:.4f}°C")
    
    print("\nTesting precipitation data...")
    # Test all methods on precipitation
    process_netcdf(precip_file, shape_file, 'precip_mean.tif', AggregationMethod.MEAN)
    process_netcdf(precip_file, shape_file, 'precip_annual_mean.tif', AggregationMethod.ANNUAL_MEAN)
    process_netcdf(precip_file, shape_file, 'precip_annual_sum_mean.tif', AggregationMethod.ANNUAL_SUM_MEAN)
    
    # Load and compare results
    precip_mean = rioxarray.open_rasterio('precip_mean.tif')
    precip_annual = rioxarray.open_rasterio('precip_annual_mean.tif')
    precip_sum_mean = rioxarray.open_rasterio('precip_annual_sum_mean.tif')
    
    print(f"Mean daily precipitation (simple mean): {precip_mean.mean().values:.2f} mm/day")
    print(f"Mean daily precipitation (annual mean): {precip_annual.mean().values:.2f} mm/day")
    print(f"Mean annual precipitation: {precip_sum_mean.mean().values:.2f} mm/year")
    
    # Clean up
    for f in [temp_file, precip_file, shape_file, 
              'temp_mean.tif', 'temp_annual_mean.tif',
              'precip_mean.tif', 'precip_annual_mean.tif', 'precip_annual_sum_mean.tif',
              'test_region.dbf', 'test_region.prj', 'test_region.shx']:
        if os.path.exists(f):
            os.remove(f)

if __name__ == "__main__":
    run_tests()
