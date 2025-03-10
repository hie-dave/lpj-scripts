#!/usr/bin/env python3
import netCDF4
import numpy as np
from sys import argv, exit

def calculate_1d_bounds(coords):
    """Calculate bounds for a 1D coordinate array."""
    n = len(coords)
    bounds = np.zeros((n, 4))
    
    # Calculate differences between adjacent coordinates
    diffs = np.diff(coords)
    # Handle wrapping for longitude
    diffs = np.where(diffs > 180, diffs - 360, diffs)
    diffs = np.where(diffs < -180, diffs + 360, diffs)
    # Extend diffs to match original array size
    diffs = np.append(diffs, diffs[-1])
    
    for i in range(n):
        half_diff = diffs[i] / 2
        # Set corners in counterclockwise order
        bounds[i] = [coords[i] - half_diff,  # left/bottom
                    coords[i] + half_diff,   # right/bottom
                    coords[i] + half_diff,   # right/top
                    coords[i] - half_diff]   # left/top
    
    return bounds

def calculate_2d_bounds(lons, lats):
    """Calculate corner coordinates for a 2D lat/lon grid."""
    ny, nx = lats.shape
    lon_bounds = np.zeros((ny, nx, 4))
    lat_bounds = np.zeros((ny, nx, 4))
    
    # Calculate differences
    dlon = np.diff(lons, axis=1)
    dlat = np.diff(lats, axis=0)
    
    # Handle longitude wrapping
    dlon = np.where(dlon > 180, dlon - 360, dlon)
    dlon = np.where(dlon < -180, dlon + 360, dlon)
    
    # Extend differences to match grid size
    dlon = np.pad(dlon, ((0,0), (0,1)), mode='edge')
    dlat = np.pad(dlat, ((0,1), (0,0)), mode='edge')
    
    # Calculate corners for each cell
    for j in range(ny):
        for i in range(nx):
            half_dlon = dlon[j,i] / 2
            half_dlat = dlat[j,i] / 2
            
            # Set corners in counterclockwise order
            lon_bounds[j,i] = [lons[j,i] - half_dlon,  # bottom left
                             lons[j,i] + half_dlon,   # bottom right
                             lons[j,i] + half_dlon,   # top right
                             lons[j,i] - half_dlon]   # top left
            
            lat_bounds[j,i] = [lats[j,i] - half_dlat,  # bottom left
                             lats[j,i] - half_dlat,   # bottom right
                             lats[j,i] + half_dlat,   # top right
                             lats[j,i] + half_dlat]   # top left
    
    return lon_bounds, lat_bounds

def main(file_path):
    with netCDF4.Dataset(file_path, mode="r+") as nc:
        # Extract coordinates
        rlon = nc.variables["rlon"][:]
        rlat = nc.variables["rlat"][:]
        lon = nc.variables["lon"][:]  # 2D array
        lat = nc.variables["lat"][:]  # 2D array
        
        # Calculate bounds for rotated coordinates (1D)
        rlon_bounds = calculate_1d_bounds(rlon)
        rlat_bounds = calculate_1d_bounds(rlat)
        
        # Calculate bounds for true coordinates (2D)
        lon_bounds, lat_bounds = calculate_2d_bounds(lon, lat)
        
        # Add vertices dimension if it doesn't exist
        if "vertices" not in nc.dimensions:
            nc.createDimension("vertices", 4)
        
        # Create or update bounds variables
        # 1D bounds for rotated coordinates
        for var_name, bounds, dim_name in [
            ("rlon_bounds", rlon_bounds, "rlon"),
            ("rlat_bounds", rlat_bounds, "rlat")
        ]:
            if var_name in nc.variables:
                bounds_var = nc.variables[var_name]
            else:
                bounds_var = nc.createVariable(var_name, "f8", (dim_name, "vertices"))
                nc.variables[dim_name].bounds = var_name
            bounds_var[:] = bounds
            bounds_var.units = nc.variables[dim_name].units
            bounds_var.long_name = f"{dim_name} cell bounds"
        
        # 2D bounds for true coordinates
        for var_name, bounds, coord_name in [
            ("lon_bounds", lon_bounds, "lon"),
            ("lat_bounds", lat_bounds, "lat")
        ]:
            if var_name in nc.variables:
                bounds_var = nc.variables[var_name]
            else:
                bounds_var = nc.createVariable(var_name, "f8", ("rlat", "rlon", "vertices"))
                nc.variables[coord_name].bounds = var_name
            bounds_var[:] = bounds
            bounds_var.units = nc.variables[coord_name].units
            bounds_var.long_name = f"{coord_name} cell bounds"

if __name__ == "__main__":
    if len(argv) == 1:
        print(f"Usage: {argv[0]} <infile.nc>")
        exit(1)
    main(argv[1])
