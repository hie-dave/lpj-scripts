#!/usr/bin/env python3
#
# Generate a gridlist file readable by LPJ-Guess, which contains all gridcells
# in the input netcdf file which are also within the bounds defined by a
# landmask shape file.
#
# Usage: generate_gridlist.py <infile.nc> <landmask.shp> <outfile.txt>
#

from sys import argv
from argparse import ArgumentParser

class Options:
    """
    Class for storing CLI arguments from the user.
    """
    def __init__(self, in_file: str, out_file: str, shape_file: str,
                 indices: bool):
        self.in_file = in_file
        self.out_file = out_file
        self.shape_file = shape_file
        self.indices = indices

def parse_args(argv: list[str]) -> Options:
    """
    Parse CLI arguments, return a parsed options object.

    @param argv: Raw CLI arguments.

    @return Parsed Options object.
    """
    parser = ArgumentParser(prog=argv[0], description = "Generate a gridlist from a netcdf file, filtering out any ocean gridcells based on a given shapefile")
    parser.add_argument("-i", "--in-file", required = True, help = "Path to a netcdf file on the target grid")
    parser.add_argument("-o", "--out-file", required = True, help = "Path to the output gridlist file")
    parser.add_argument("-s", "--shape-file", required = True, help = "Path to the shapefile")
    parser.add_argument("--indices", action = "store_true", default = False, help = "If set, the output file will be in X Y format, suitable for the CF input module (rather than LON LAT format)")

    parsed = parser.parse_args(argv[1:])
    return Options(parsed.in_file, parsed.out_file, parsed.shape_file,
                   parsed.indices)

if __name__ == "__main__":
    opts = parse_args(argv)

from geopandas import GeoDataFrame, points_from_xy, read_file, sjoin
from numpy import meshgrid
from io import TextIOWrapper
from netCDF4 import Dataset
from traceback import format_exc

LAT_NAMES: list[str] = ["lat", "latitude"]
LON_NAMES: list[str] = ["lon", "longitude"]

def find_variable(nc: Dataset, names: list[str]):
    for name in names:
        if name in nc.variables:
            return nc.variables[name]
    raise ValueError(f"Failed to find variable. Tried names: {names}")

def index_of(needle, haystack) -> int:
    for i in range(len(haystack)):
        if haystack[i] == needle:
            return i
    return -1

def gen_gridlist(shapefile: str, nc: Dataset, gridlist: TextIOWrapper, indices: bool = False):
    # Get the longitude/latitude variables from the netcdf file.
    lat = find_variable(nc, LAT_NAMES)
    lon = find_variable(nc, LON_NAMES)

    # Read coordinates from the netcdf file.
    latitudes = lat[:]
    longitudes = lon[:]

    lon_grid, lat_grid = meshgrid(longitudes, latitudes)

    # Flatten the grid arrays
    lon_flat = lon_grid.flatten()
    lat_flat = lat_grid.flatten()

    # Create a GeoDataFrame of the grid points.
    df = GeoDataFrame(geometry = points_from_xy(lon_flat, lat_flat))

    # Load the land mask shapefile.
    land = read_file(shapefile)

    # Perform a spatial join to find points within land areas.
    points = sjoin(df, land, op='within')

	# Lookup tables for finding an index for a given coordinate value.
    lx = {}
    ly = {}

    # Write to gridlist.
    for point in points.geometry:
        x = point.x
        y = point.y
        if indices:
            if x in lx.keys():
                x = lx[x]
            else:
                x = index_of(x, longitudes)
                lx[point.x] = x

            if y in ly.keys():
                y = ly[y]
            else:
                y = index_of(y, latitudes)
                ly[point.y] = y

            if x == -1:
                raise ValueError(f"Failed to find longitude value: {point.x}")
            if y == -1:
                raise ValueError(f"Failed to find latitude value: {point.y}")
        gridlist.writelines(f"{x} {y}\n")

if __name__ == "__main__":
    try:
        with open(opts.out_file, "w") as gridlist:
            with Dataset(opts.in_file, "r") as nc:
                gen_gridlist(opts.shape_file, nc, gridlist, opts.indices)
    except Exception as error:
        # Basic error handling.
        print(format_exc())
        exit(1)
