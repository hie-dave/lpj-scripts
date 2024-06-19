#!/usr/bin/env python3
#
# Generate a gridlist file readable by LPJ-Guess, which contains all gridcells
# in the input netcdf file which are also within the bounds defined by a
# landmask shape file.
#
# Usage: generate_gridlist.py <infile.nc> <landmask.shp> <outfile.txt>
#

from sys import argv, exit

if __name__ == "__main__":
	err = not ("-h" in argv or "--help" in argv)
	if len(argv) != 4:
		print(f"Usage: {argv[0]} <infile.nc> <landmask.shp> <outfile.txt>")
		exit(1 if err else 0)

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

def gen_gridlist(shapefile: str, nc: Dataset, gridlist: TextIOWrapper):
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

	# Write to gridlist.
	gridlist.writelines([f"{point.x} {point.y}\n" for point in points.geometry])

if __name__ == "__main__":
	# Parse CLI args
	in_file = argv[1]
	shapefile = argv[2]
	out_file = argv[3]

	try:
		with open(out_file, "w") as gridlist:
			with Dataset(in_file, "r") as nc:
				gen_gridlist(shapefile, nc, gridlist)
	except Exception as error:
		# Basic error handling.
		print(format_exc())
		exit(1)
