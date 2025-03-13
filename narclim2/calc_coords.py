#!/usr/bin/env python3
#
# Calculate interpretable lon/lat coordinates from a rotated-pole NetCDF file.
#

from sys import argv, exit
import xarray as xr
import numpy as np
import os, sys
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import pyproj

def read_data(file_name):
    """
    Read netcdf file and return arrays: lat, lon, var, pole_longitude, pole_latitude.
    """
    ds = xr.open_dataset(file_name)
    var  = ds.huss[0,:,:]
    rlat = ds.rlat[:]
    rlon = ds.rlon[:]
    pole = ds.crs

    try:
        if hasattr(pole,"grid_north_pole_longitude"):
            px = pole.attrs["grid_north_pole_longitude"]
        if hasattr(pole,"grid_north_pole_latitude"):
            py = pole.attrs["grid_north_pole_latitude"]
    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise

    return rlon, rlat, var, px, py

# Function to calculate bounds for wrapped coordinates
def calculate_wrapped_bounds(coords):
    bounds = np.zeros((len(coords), 2))
    diffs = np.diff(coords)  # Differences between adjacent coordinates
    diffs = np.where(diffs > 180, diffs - 360, diffs)  # Handle wrapping
    diffs = np.where(diffs < -180, diffs + 360, diffs)  # Handle wrapping

    # Calculate left bounds
    bounds[1:, 0] = coords[:-1] + diffs / 2  # Midpoints between adjacent coordinates
    bounds[0, 0] = coords[0] - diffs[0] / 2  # Extrapolate for the first point

    # Calculate right bounds
    bounds[:-1, 1] = bounds[1:, 0]  # Right bounds are the next point"s left bound
    bounds[-1, 1] = coords[-1] + diffs[-1] / 2  # Extrapolate for the last point

    return bounds

def plot_coords(rlon: np.ndarray, rlat: np.ndarray, var: np.ndarray, pole_lon: float, pole_lat: float, file_name: str):
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()
    ax.set_extent([110, 150, -60, 12], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.OCEAN, color="white", zorder=0)
    ax.add_feature(cartopy.feature.LAND, color="lightgray",zorder=0,
                   linewidth=0.5, edgecolor="black")
    ax.gridlines(draw_labels=True, linewidth=0.5, color="gray",
                 xlocs=range(-180,180,15), ylocs=range(-90,90,15))
    ax.coastlines(resolution="50m", linewidth=0.3, color="black")
    ax.set_title("Py: de-rotated grid (cartopy)", fontsize=10, fontweight="bold")
    colormap = "RdYlBu_r"
    crs = ccrs.RotatedPole(pole_longitude=pole_lon, pole_latitude=pole_lat)
    ax.contourf(rlon, rlat, var, levels=15, cmap=colormap, transform=crs)
    plt.savefig(file_name, bbox_inches="tight", dpi=200)

def main(fname: str):
    """
    Main function.
    """
    # read file content and return relevant variables
    rlon, rlat, var, pole_lon, pole_lat = read_data(fname)

    # Create a meshgrid of rotated coordinates
    rlon2D, rlat2D = np.meshgrid(rlon, rlat)

    # Define the transformation using pyproj
    transform = pyproj.Transformer.from_crs(
        f"+proj=ob_tran +o_proj=latlon +o_lon_p={pole_lon} +o_lat_p={pole_lat} +lon_wrap=180",
        "EPSG:4326",  # Convert to standard lat/lon
        always_xy=True
    )

    # Convert to geographic coordinates
    lon2D, lat2D = transform.transform(rlon2D, rlat2D)
    print(f"lon2D = [{np.min(lon2D)}..{np.max(lon2D)}] {np.shape(lon2D)}")
    print(f"lat2D = [{np.min(lat2D)}..{np.max(lat2D)}] {np.shape(lat2D)}")

    # plot_file = "plot_Python_curvilinear_grid_derotated_1.png"
    # plot_coords(lon2D, lat2D, var, pole_lon, pole_lat, plot_file)

if __name__ == "__main__":
    if len(argv) != 2:
        print(f"Usage: {argv[0]} <infile.nc>")
        exit(1)
    main(argv[1])
