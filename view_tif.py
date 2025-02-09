#!/usr/bin/env python3

import rioxarray
import matplotlib.pyplot as plt
import sys

def plot_tif(tif_file):
    """Quick visualization of a GeoTIFF file"""
    # Read the file
    da = rioxarray.open_rasterio(tif_file)
    
    # Create figure
    plt.figure(figsize=(10, 8))
    
    # Plot with a colorbar
    im = plt.imshow(da[0], cmap='viridis')
    plt.colorbar(im, label=f'Value ({getattr(da, "units", "unknown units")})')
    
    # Add title
    plt.title(f'Visualization of {tif_file}')
    
    # Show lat/lon grid
    plt.grid(True, linestyle='--', alpha=0.5)
    
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python view_tif.py <tif_file>")
        sys.exit(1)
    
    plot_tif(sys.argv[1])
