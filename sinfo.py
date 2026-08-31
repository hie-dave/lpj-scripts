#!/usr/bin/env python3
import netCDF4
import numpy
import sys

STD_LON = "longitude"
STD_LAT = "latitude"

def var_from_std_name(nc: netCDF4.Dataset, std_name: str):
    """
    Retrieve a variable from a netCDF dataset based on its standard name.

    Parameters:
        nc (netCDF4.Dataset): The netCDF dataset object.
        std_name (str): The standard name of the variable to retrieve.

    Returns:
        netCDF4.Variable: The variable with the specified standard name, or None if not found.
    """
    for _, var in nc.variables.items():
        if getattr(var, "standard_name", None) == std_name:
            return var
    raise RuntimeError(f"No variable with standard name '{std_name}' found in the dataset.")

def print_sinfo(nc: netCDF4.Dataset):
    """
    Print information about the spatial variables (latitude and longitude) in a netCDF file.

    Parameters:
        nc (netCDF4.Dataset): The netCDF dataset object.
    """

    lat = var_from_std_name(nc, STD_LAT)
    lon = var_from_std_name(nc, STD_LON)

    if len(lat.dimensions) != 1:
        raise RuntimeError("Latitude variable is not 1-dimensional.")

    if len(lon.dimensions) != 1:
        raise RuntimeError("Longitude variable is not 1-dimensional.")

    lats = lat[:]
    lons = lon[:]

    dlat = numpy.unique(numpy.diff(lats))
    dlon = numpy.unique(numpy.diff(lons))

    # Print latitude and longitude coverage and resolution.
    print(f"Latitude coverage: {lats[0]} to {lats[-1]}")
    print(f"Longitude coverage: {lons[0]} to {lons[-1]}")
    if len(dlat) == 1:
        print(f"Latitude resolution: {dlat[0]}")
    else:
        print(f"Latitude resolution: {dlat} (multiple unique resolutions found)")
    if len(dlon) == 1:
        print(f"Longitude resolution: {dlon[0]}")
    else:
        print(f"Longitude resolution: {dlon} (multiple unique resolutions found)")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: {} <file>".format(sys.argv[0]))
        sys.exit(1)

    filename = sys.argv[1]
    with netCDF4.Dataset(filename, 'r') as nc:
        print_sinfo(nc)
