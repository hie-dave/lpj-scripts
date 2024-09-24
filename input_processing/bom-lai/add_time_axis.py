#!/usr/bin/env python3

import datetime
import netCDF4 as nc
import numpy as np
from sys import argv, exit
import re
import os

# Open the existing NetCDF file
input_file = argv[1]
output_file = argv[2]

# Get timestamp from file name.
date_fmt = "%Y-%m-%d"
pattern = r"(\d{4}-\d{2}-\d{2})"
rx = re.compile(pattern)
match = rx.search(input_file)
if match is None:
    print(f"Failed to parse valid timestamp from file name '{input_file}'")
    exit(1)
timestr = match[0]
timestamp = datetime.datetime.strptime(timestr, date_fmt)

data_var_name = "Band1"
new_data_var_name = data_var_name

# Open the input file for reading
with nc.Dataset(input_file, 'r') as src:

    # Create a new NetCDF file to write to
    with nc.Dataset(output_file, 'w') as dst:

        # Copy the existing dimensions (lat, lon)
        for name, dimension in src.dimensions.items():
            dst.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
        
        # Add the time dimension
        dst.createDimension('time', 1)
        
        # Copy the existing variables (lat, lon)
        for name, variable in src.variables.items():
            if name != data_var_name:
                out_var = dst.createVariable(name, variable.datatype, variable.dimensions)
                out_var.setncatts({k: variable.getncattr(k) for k in variable.ncattrs()})
                out_var[:] = variable[:]
        
        # Create the time variable
        time_var = dst.createVariable('time', 'f8', ('time',))
        time_var.units = 'hours since 1900-01-01 00:00:00'
        time_var.calendar = 'gregorian'
        
        # Set the time variable value
        time_point = nc.date2num([timestamp], units=time_var.units, calendar=time_var.calendar)
        time_var[:] = time_point
        
        # Add a new variable with the added time dimension
        data_var = src.variables[data_var_name]  # Replace with your data variable name
        new_var = dst.createVariable(new_data_var_name, data_var.datatype, ('time', 'lat', 'lon'))  # Replace with desired variable name
        new_var.setncatts({k: data_var.getncattr(k) for k in data_var.ncattrs()})
        
        # Copy the data with an added time dimension
        new_var[0, :, :] = data_var[:]
