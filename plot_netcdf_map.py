#!/usr/bin/env python3
"""
Script to plot a spatial map of the 's0' variable at the first timestep from a NetCDF file.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import netCDF4 as nc
import os

def plot_netcdf_map(netcdf_file, variable='s0', timestep=0, output_file=None):
    """
    Plot a spatial map of a variable at a specific timestep from a NetCDF file.
    
    Parameters:
    -----------
    netcdf_file : str
        Path to the NetCDF file
    variable : str
        Name of the variable to plot (default: 's0')
    timestep : int
        Timestep index to plot (default: 0, which is the first timestep)
    output_file : str, optional
        Path to save the output plot. If None, the plot will be displayed.
    """
    # Open the NetCDF file
    print(f"Opening NetCDF file: {netcdf_file}")
    try:
        ds = nc.Dataset(netcdf_file, 'r')
    except Exception as e:
        print(f"Error opening NetCDF file: {e}")
        return
    
    # Print information about the file
    print("\nNetCDF file information:")
    print(f"Dimensions: {list(ds.dimensions.keys())}")
    print(f"Variables: {list(ds.variables.keys())}")
    
    # Check if the variable exists
    if variable not in ds.variables:
        print(f"Error: Variable '{variable}' not found in the NetCDF file.")
        print(f"Available variables: {list(ds.variables.keys())}")
        ds.close()
        return
    
    # Get the variable data
    var_data = ds.variables[variable]
    print(f"\nVariable '{variable}' information:")
    print(f"Dimensions: {var_data.dimensions}")
    print(f"Shape: {var_data.shape}")
    
    # Check if the variable has a time dimension
    if 'time' in var_data.dimensions:
        time_index = var_data.dimensions.index('time')
        if timestep >= var_data.shape[time_index]:
            print(f"Error: Timestep {timestep} is out of range. The variable has {var_data.shape[time_index]} timesteps.")
            ds.close()
            return
        
        # Extract the data for the specified timestep
        if time_index == 0:
            data = var_data[timestep, :, :]
        elif time_index == 1:
            data = var_data[:, timestep, :]
        elif time_index == 2:
            data = var_data[:, :, timestep]
        else:
            print(f"Error: Unsupported dimension order for variable '{variable}'.")
            ds.close()
            return
    else:
        print(f"Warning: Variable '{variable}' does not have a time dimension. Using all data.")
        data = var_data[:]
    
    # Get longitude and latitude data
    if 'lon' in ds.variables and 'lat' in ds.variables:
        lons = ds.variables['lon'][:]
        lats = ds.variables['lat'][:]
    else:
        print("Warning: 'lon' or 'lat' variables not found. Using array indices instead.")
        lons = np.arange(data.shape[1])
        lats = np.arange(data.shape[0])
    
    # Create a meshgrid for plotting
    lon_mesh, lat_mesh = np.meshgrid(lons, lats)
    
    # Get time information if available
    time_str = "first timestep"
    if 'time' in ds.variables:
        time_var = ds.variables['time']
        if hasattr(time_var, 'units'):
            time_units = time_var.units
            if timestep < len(time_var):
                time_value = time_var[timestep]
                calendar = time_var.calendar if hasattr(time_var, 'calendar') else 'standard'
                # Parse the time value.
                dt = nc.num2date(time_value, time_units, calendar)
                time_str = f"{dt.strftime('%Y-%m-%d %H:%M:%S')} (timestep {timestep})"
    
    # Create the plot
    plt.figure(figsize=(10, 8))
    
    # Use Cartopy for map projection if the coordinates are in degrees
    if np.max(lons) <= 360 and np.min(lons) >= -180 and np.max(lats) <= 90 and np.min(lats) >= -90:
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()
        ax.gridlines(draw_labels=True)
        
        # Adjust longitude values if they're in 0-360 range to -180 to 180 range
        if np.min(lons) >= 0 and np.max(lons) > 180:
            lons_plot = np.where(lons > 180, lons - 360, lons)
            lon_mesh, lat_mesh = np.meshgrid(lons_plot, lats)
        
        # Create the plot
        im = ax.pcolormesh(lon_mesh, lat_mesh, data, transform=ccrs.PlateCarree(), cmap='viridis')
    else:
        # For non-geographic coordinates, use a regular plot
        ax = plt.axes()
        im = ax.pcolormesh(lon_mesh, lat_mesh, data, cmap='viridis')
        ax.set_xlabel('Longitude Index')
        ax.set_ylabel('Latitude Index')
    
    # Add colorbar and labels
    cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.02)
    if hasattr(var_data, 'units'):
        cbar.set_label(f"{variable} ({var_data.units})")
    else:
        cbar.set_label(variable)
    
    # Set title
    plt.title(f"{variable} at {time_str}")
    
    # Save or display the plot
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_file}")
    else:
        plt.tight_layout()
        plt.show()
    
    # Close the NetCDF file
    ds.close()

def main():
    parser = argparse.ArgumentParser(description='Plot a spatial map from a NetCDF file.')
    parser.add_argument('netcdf_file', help='Path to the NetCDF file')
    parser.add_argument('--variable', '-v', default='s0', help='Variable to plot (default: s0)')
    parser.add_argument('--timestep', '-t', type=int, default=0, help='Timestep index to plot (default: 0)')
    parser.add_argument('--output', '-o', help='Output file path (optional)')
    
    args = parser.parse_args()
    
    plot_netcdf_map(args.netcdf_file, args.variable, args.timestep, args.output)

if __name__ == '__main__':
    main()
