#!/usr/bin/env python3
import netCDF4
import numpy as np
from sys import argv, exit

def calculate_grid_corners(lons, lats):
    """Calculate four corners for each grid cell."""
    # Create arrays for the corners
    nlon, nlat = len(lons), len(lats)
    lon_corners = np.zeros((nlon, nlat, 4))
    lat_corners = np.zeros((nlon, nlat, 4))
    
    # Calculate coordinate differences
    dlon = np.diff(lons)
    dlat = np.diff(lats)
    
    # Handle wrapping for longitude
    dlon = np.where(dlon > 180, dlon - 360, dlon)
    dlon = np.where(dlon < -180, dlon + 360, dlon)
    
    # Extend differences to match grid size
    dlon = np.append(dlon, dlon[-1])
    dlat = np.append(dlat, dlat[-1])
    
    # Create meshgrid for vectorized operations
    lon_grid, lat_grid = np.meshgrid(lons, lats, indexing='ij')
    dlon_grid, dlat_grid = np.meshgrid(dlon, dlat, indexing='ij')
    
    # Calculate corners for each cell (counterclockwise from bottom left)
    for i in range(nlon):
        for j in range(nlat):
            # Get half-steps
            half_dlon = dlon_grid[i,j] / 2
            half_dlat = dlat_grid[i,j] / 2
            
            # Bottom left
            lon_corners[i,j,0] = lon_grid[i,j] - half_dlon
            lat_corners[i,j,0] = lat_grid[i,j] - half_dlat
            
            # Bottom right
            lon_corners[i,j,1] = lon_grid[i,j] + half_dlon
            lat_corners[i,j,1] = lat_grid[i,j] - half_dlat
            
            # Top right
            lon_corners[i,j,2] = lon_grid[i,j] + half_dlon
            lat_corners[i,j,2] = lat_grid[i,j] + half_dlat
            
            # Top left
            lon_corners[i,j,3] = lon_grid[i,j] - half_dlon
            lat_corners[i,j,3] = lat_grid[i,j] + half_dlat
    
    return lon_corners, lat_corners

def main(file_path):
    with netCDF4.Dataset(file_path, mode="r+") as nc:
        # Extract the center coordinates
        rlon = nc.variables['rlon'][:]
        rlat = nc.variables['rlat'][:]
        
        # Calculate corners
        rlon_corners, rlat_corners = calculate_grid_corners(rlon, rlat)
        
        # Add vertices dimension
        if 'vertices' not in nc.dimensions:
            nc.createDimension('vertices', 4)
        
        # Create or update bounds variables
        for var_name, corners, dim_name in [
            ('rlon_vertices', rlon_corners, 'rlon'),
            ('rlat_vertices', rlat_corners, 'rlat')
        ]:
            if var_name in nc.variables:
                bounds_var = nc.variables[var_name]
            else:
                bounds_var = nc.createVariable(var_name, 'f8', (dim_name, dim_name, 'vertices'))
                nc.variables[dim_name].bounds = var_name
            
            bounds_var[:] = corners
            bounds_var.units = nc.variables[dim_name].units
            bounds_var.long_name = f'{dim_name} cell vertices'

if __name__ == '__main__':
    if len(argv) == 1:
        print(f"Usage: {argv[0]} <infile.nc>")
        exit(1)
    main(argv[1])
