#!/usr/bin/env python3
from argparse import ArgumentParser
import netCDF4
import numpy as np
from sys import argv, exit
from shutil import copy2

class Options:
    """
    CLI Options.
    """
    def __init__(self, in_file: str, out_file: str):
        """
        Constructor.

        Args:
            in_file: Input file.
            out_file: Output file.
        """
        self.in_file = in_file
        self.out_file = out_file

def parse_args(argv: list[str]) -> Options:
    """
    Parse the command line arguments.

    Args:
        argv: The command line arguments.

    Returns:
        The parsed options.
    """
    parser = ArgumentParser(prog=argv[0], description="Set the contents of the "
        "specified variable in a netcdf file to the values given in a plaintext"
        " file.", exit_on_error=True)
    parser.add_argument(
        "-i", "--in-file", required=True, type=str, help="Input file.")
    parser.add_argument(
        "-o", "--out-file", type=str, help="Output file. If not specified, the "
        "input file will be modified in place.")
    parsed = parser.parse_args(argv[1:])
    return Options(**vars(parsed))

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

def set_bounds(nc: netCDF4.Dataset):
    """
    Set bounds variables in a NetCDF file.

    This function calculates bounds for rotated coordinates (rlon, rlat)
    and true coordinates (lon, lat) and adds them to the file.
    """
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

def main(opts: Options):
    if opts.out_file is None:
        # Edit in-place.
        opts.out_file = opts.in_file
    elif opts.in_file != opts.out_file:
        # Copy input file to output file, then edit output file in-place.
        copy2(opts.in_file, opts.out_file)
    # else output file is the same as input file, so edit in-place.

    # Open output file for in-place editing.
    with netCDF4.Dataset(opts.out_file, mode="r+") as nc:
        set_bounds(nc)

if __name__ == "__main__":
    opts = parse_args(argv)
    try:
        main(opts)
    except BaseException as error:
        print_exc()
        exit(1)
