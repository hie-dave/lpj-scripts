#!/usr/bin/env python3
import netCDF4, numpy
from sys import argv, exit


NAME_LON = "lon"
NAME_LAT = "lat"
NAME_RLON = "rlon"
NAME_RLAT = "rlat"

STD_LAT = "latitude"
STD_LON = "longitude"

STD_RLON = "grid_longitude"
STD_RLAT = "grid_latitude"

LONG_RLON = "longitude in rotated pole grid"
LONG_RLAT = "latitude in rotated pole grid"

LONG_LON = "longitude"
LONG_LAT = "latitude"

# Lookup table of standard names.
std_names = {
    NAME_LON: STD_LON,
    NAME_LAT: STD_LAT,
    NAME_RLON: STD_RLON,
    NAME_RLAT: STD_RLAT,
}

long_names = {
    NAME_LON: LONG_LON,
    NAME_LAT: LONG_LAT,
    NAME_RLON: LONG_RLON,
    NAME_RLAT: LONG_RLAT,
}

def get_std_name(var_name: str) -> str:
    """
    Get the standard name for a variable.
    """
    # var_name will be "rlon_bounds", "rlat_bounds", "lon_bounds", etc.
    # Remove "_bounds" suffix.
    var_name = var_name.replace("_bounds", "")

    return std_names.get(var_name, var_name)

def get_long_name(var_name: str) -> str:
    """
    Get the long name for a variable.
    """
    # var_name will be "rlon_bounds", "rlat_bounds", "lon_bounds", etc.
    # Remove "_bounds" suffix.
    var_name = var_name.replace("_bounds", "")

    return long_names.get(var_name, var_name)

def get_var(nc: netCDF4.Dataset, var_name: str) -> netCDF4.Variable:
    """
    Get a bounds variable from a NetCDF file or create it if it doesn't exist.
    """
    if var_name in nc.variables:
        return nc.variables[var_name]
    
    # Create variable with appropriate dimensions
    base_var_name = var_name.replace("_bounds", "")
    base_var = nc.variables[base_var_name]
    
    if len(base_var.dimensions) == 1:
        # 1D variables like rlon, rlat
        dims = (base_var.dimensions[0], "vertices")
    else:
        # 2D variables like lon, lat
        dims = (base_var.dimensions[0], base_var.dimensions[1], "vertices")
    
    var = nc.createVariable(var_name, "f8", dims)
    var.units = "degrees"
    var.standard_name = get_std_name(var_name)
    var.long_name = get_long_name(var_name)
    return var

# Function to calculate bounds for wrapped coordinates
def calculate_wrapped_bounds(coords):
    bounds = numpy.zeros((len(coords), 2))
    diffs = numpy.diff(coords)  # Differences between adjacent coordinates
    diffs = numpy.where(diffs > 180, diffs - 360, diffs)  # Handle wrapping
    diffs = numpy.where(diffs < -180, diffs + 360, diffs)  # Handle wrapping

    # Calculate left bounds
    bounds[1:, 0] = coords[:-1] + diffs / 2  # Midpoints between adjacent coordinates
    bounds[0, 0] = coords[0] - diffs[0] / 2  # Extrapolate for the first point

    # Calculate right bounds
    bounds[:-1, 1] = bounds[1:, 0]  # Right bounds are the next point's left bound
    bounds[-1, 1] = coords[-1] + diffs[-1] / 2  # Extrapolate for the last point

    return bounds

def get_dim(nc: netCDF4.Dataset, name: str, size: int) -> netCDF4.Dimension:
    """
    Get a dimension from a NetCDF file or create it if it doesn't exist.
    """
    if name in nc.dimensions:
        return nc.dimensions[name]
    # Create dimension.
    dim = nc.createDimension(name, size)
    return dim

def main(file_path):
    # Open the file in "r+" mode for in-place editing
    with netCDF4.Dataset(file_path, mode="r+") as nc:
        # Extract the center coordinates (rlon, rlat)
        rlon = nc.variables["rlon"][:]
        rlat = nc.variables["rlat"][:]
        lon = nc.variables["lon"][:]
        lat = nc.variables["lat"][:]

        # Calculate bounds for rlon and rlat
        rlon_bounds = calculate_wrapped_bounds(rlon)
        rlat_bounds = calculate_wrapped_bounds(rlat)

        # Add the bounds variables to the file
        vertices_dim = get_dim(nc, "vertices", 2)

        # Set bounds for rotated coordinates (1D)
        rlon_bounds_var = get_var(nc, "rlon_bounds")
        rlon_bounds_var[:] = rlon_bounds

        rlat_bounds_var = get_var(nc, "rlat_bounds")
        rlat_bounds_var[:] = rlat_bounds

        # Set bounds for geographic coordinates (2D)
        # For each point in the rotated grid, calculate geographic bounds
        lon_bounds = numpy.zeros(lon.shape + (2,)) # (i, j, 2)
        lat_bounds = numpy.zeros(lat.shape + (2,)) # (i, j, 2)

        lon_diffs = numpy.diff(lon)
        lat_diffs = numpy.diff(lat)
        
        # Use the rotated bounds to calculate geographic bounds
        for i in range(lon.shape[0]):
            for j in range(lon.shape[1]):
                # Get the rotated coordinates for this point
                rlon_pt = rlon[j]  # Note: lon[i,j] corresponds to rlat[i], rlon[j]
                rlat_pt = rlat[i]
                
                # Calculate bounds in rotated coordinates
                rlon_bounds_pt = [rlon_bounds[j,0], rlon_bounds[j,1]]
                rlat_bounds_pt = [rlat_bounds[i,0], rlat_bounds[i,1]]
                
                # Store the bounds
                lon_bounds[i,j,:] = [lon[i,j] - abs(rlon_bounds_pt[1] - rlon_bounds_pt[0])/2,
                                   lon[i,j] + abs(rlon_bounds_pt[1] - rlon_bounds_pt[0])/2]
                lat_bounds[i,j,:] = [lat[i,j] - abs(rlat_bounds_pt[1] - rlat_bounds_pt[0])/2,
                                   lat[i,j] + abs(rlat_bounds_pt[1] - rlat_bounds_pt[0])/2]

        # Create and set the bounds variables
        lon_bounds_var = get_var(nc, "lon_bounds")
        lon_bounds_var[:] = lon_bounds
        
        lat_bounds_var = get_var(nc, "lat_bounds")
        lat_bounds_var[:] = lat_bounds

        # Add attributes to describe the bounds
        nc.variables["rlon"].bounds = "rlon_bounds"
        nc.variables["rlat"].bounds = "rlat_bounds"
        nc.variables["lon"].bounds = "lon_bounds"
        nc.variables["lat"].bounds = "lat_bounds"

if __name__ == '__main__':
    if len(argv) == 1:
        print(f"Usage: {argv[0]} <infile.nc>")
        exit(1)
    main(argv[1])
