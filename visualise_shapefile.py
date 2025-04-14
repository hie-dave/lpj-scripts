#!/usr/bin/env python

from geopandas import read_file
import matplotlib.pyplot as plt
from sys import argv, exit
from os.path import basename

def die(msg: str):
    print(msg)
    exit(1)

def main(shapefile: str):
    # Read the shapefile
    gdf = read_file(shapefile)

    # Print some info about the shapefile
    print(gdf.head())
    print(gdf.crs)

    # Plot the shapefile
    gdf.plot(edgecolor='black', figsize=(10, 10))
    plt.title(basename(shapefile))
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    # Path to your shapefile (replace this with the actual path)
    if len(argv) != 2:
        die(f"Usage: {argv[0]} <shapefile.shp>")
    main(argv[1])
