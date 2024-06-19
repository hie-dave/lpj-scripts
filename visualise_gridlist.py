#!/usr/bin/env python3
#
# Visualise a gridlist file. This should be a plaintext file containing
# one or more rows of space-separated lon-lat pairs.
#

from sys import argv, exit

if len(argv) != 2:
    print(f"Usage: {argv[0]} <gridlist.txt>")
    exit(1)

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point

# Load the gridlist file
gridlist_path = argv[1]
grid_data = pd.read_csv(gridlist_path, delim_whitespace=True, header=None, names=['longitude', 'latitude'])

# Convert the dataframe to a GeoDataFrame
geometry = [Point(xy) for xy in zip(grid_data['longitude'], grid_data['latitude'])]
gdf = gpd.GeoDataFrame(grid_data, geometry=geometry)

# Load a basemap to provide context (e.g., a world map)
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

# Plot the basemap and the grid points
fig, ax = plt.subplots(figsize=(10, 10))
world.plot(ax=ax, color='lightgrey')
gdf.plot(ax=ax, marker='o', color='red', markersize=5)

# Set the map's title and labels
ax.set_title('Grid Coordinates Visualization')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# Set axis limits to focus on the desired region (e.g., Australia)
ax.set_xlim([110, 155])
ax.set_ylim([-45, -10])

# Show the plot
plt.show()
