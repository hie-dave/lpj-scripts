#!/usr/bin/env python3

"""
Calculate Thome (long-term mean maximum temperature of the warmest month) from air temperature data.

This script processes gridded air temperature data from a NetCDF file to calculate Thome.
The input file should be CF-compliant and contain air temperature data with dimensions (time, lat, lon)
or similar. The script will automatically detect dimension order using standard_name attributes.
"""

import argparse
import datetime
import numpy
from netCDF4 import Dataset
import os
from sys import argv, exit
from traceback import format_exc

from ozflux_common import *
from ozflux_logging import *
from ozflux_netcdf import *

# Standard name for air temperature variable.
STD_AIR_TEMP = "air_temperature"

class Options:
    """Command line options for calc_thome.py."""

    def __init__(self, input_file: str, output_file: str, log_level: LogLevel,
                 show_progress: bool):
        """Initialize with default values."""
        self.input_file = input_file
        self.output_file = output_file
        self.log_level = log_level
        self.show_progress = show_progress

def parse_args(argv: list[str]) -> Options:
    """Parse command line arguments and return Options object."""
    parser = argparse.ArgumentParser(description = "Calculate Thome from air temperature data")
    parser.add_argument("-v", "--verbosity", type = int, help = f"Logging verbosity ({LogLevel.ERROR}-{LogLevel.DEBUG}, default: {LogLevel.INFORMATION})", default = LogLevel.INFORMATION)
    parser.add_argument("-p", "--show-progress", action = "store_true", help = "Report progress")
    parser.add_argument("-i", "--in-file", help = "Input NetCDF file containing air temperature data")
    parser.add_argument("-o", "--out-file", help = "Output NetCDF file to write Thome values")

    args = parser.parse_args(argv[1:])
    opts = Options(args.in_file, args.out_file, args.verbosity, args.show_progress)
    return opts

def get_time_info(nc: Dataset) -> tuple[Dimension, Variable, list[datetime.datetime]]:
    """
    Get time dimension and variable information.

    @return: Time dimension, time variable, and time values.
    """
    # Find time dimension using standard_name.
    time_dim = get_dim_from_std_name(nc, STD_TIME)
    time_var = get_dimension_variable(nc, time_dim)

    # Read and parse time values.
    times = num2date(time_var[:], time_var.units, time_var.calendar,
                     only_use_cftime_datetimes = False,
                     only_use_python_datetimes = True)
    return time_dim, time_var, times

def read_temp_data(var: Variable) -> numpy.ndarray:
    """
    Read temperature data from a NetCDF variable in ℃.

    @param var: Variable to read data from.
    @return: Numpy array of temperature data.
    """
    # Read data from input file.
    data = var[:]

    # Convert to ℃ if necessary.
    if var.units != "degC":
        conversion = find_units_conversion(var.units, "degC")
        log_information(f"Converting temperature from {var.units} to ℃...")
        # Timestep is irrelevant here, so pass dummy value of 1.
        data = conversion(data, 1)
    else:
        log_diagnostic("Temperature already in ℃.")
    return data

def group_by_month(times):
    """Group time indices by month."""
    month_indices = {}
    pandas.groupby(times, lambda t: t.month)
    for i, t in enumerate(times):
        month_indices.setdefault(t.month, []).append(i)
    return month_indices

def calculate_thome_gridcell(temps: numpy.ndarray, times: list[datetime.datetime]):
    """
    Calculate Thome (long-term mean maximum temperature of warmest month) for a single gridcell.

    The input arrays should have the same shape, and should contain only the
    data for a single gridcell.

    @param temps: 1-dimensional numpy array of temperature data.
    @param times: List of datetime objects corresponding to each temperature data point.
    @return: Numpy array of Thome values.
    """
    log_diagnostic("Calculating monthly maximums...")

    df = pandas.DataFrame({"temp": temps, "time": times})
    df.set_index("time", inplace=True)
    df.groupby(lambda t: t.month)

def calculate_thome(nc: Dataset,temps: numpy.ndarray, times: list[datetime.datetime]):
    """
    Calculate Thome (long-term mean maximum temperature of warmest month).

    @param temps: Numpy array of temperature data.
    @param times: List of datetime objects corresponding to each temperature data point.
    @return: Numpy array of Thome values, with same dimensionality as the original temperature data.
    """
    # Find temperature variable
    temp_var = get_var_from_std_name(nc_in, STD_AIR_TEMP)
    log_information(f"Found temperature variable: {temp_var.name}")
    log_information(f"├── Shape: {temp_var.shape}")
    log_information(f"└── Units: {temp_var.units}")

    # Get time information.
    time_dim, time_var, times = get_time_info(nc_in)
    log_information(f"Found time dimension: {time_dim.name}")
    log_information(f"└── Time coverage: {times[0]} to {times[-1]} in {len(times)} timesteps")

    dim_lat = get_dim_from_std_name(nc_in, STD_LAT)
    var_lat = get_dimension_variable(nc_in, dim_lat)
    log_information(f"Found latitude dimension: {dim_lat.name}")
    log_information(f"├── Shape: {var_lat.shape}")
    log_information(f"└── Units: {var_lat.units}")

    dim_lon = get_dim_from_std_name(nc_in, STD_LON)
    var_lon = get_dimension_variable(nc_in, dim_lon)
    log_information(f"Found longitude dimension: {dim_lon.name}")
    log_information(f"├── Shape: {var_lon.shape}")
    log_information(f"└── Units: {var_lon.units}")

    # Read temperature data.
    log_information("Reading temperature data...")
    temp_data = read_temp_data(temp_var)

    if len(temp_var.dimensions) == 1:
        # Flatten temperature data to 1-dimensional array.
        # TODO: test this
        temp_data = temp_data.reshape(-1)
        return calculate_thome_gridcell(temp_data, times)

    # Get the index of the latitude and longitude dimensions in the temperature
    # variable's dimensions.
    index_latitude = index_of_throw(dim_lat, temp_var.dimensions,
    lambda: "Latitude dimension not found in temperature variable's dimensions")
    index_longitude = index_of_throw(dim_lon, temp_var.dimensions,
    lambda: "Longitude dimension not found in temperature variable's dimensions")

    # Iterate through gridcells.
    for i in range(temps.shape[index_latitude]):
        for j in range(temps.shape[index_longitude]):
            thome = calculate_thome_gridcell(i, j, temp_data, times)
            write_output(nc_in, opts.output_file, thome)
            log_information("Done!")

    # Group data by month.
    month_indices = group_by_month(times)

    # Calculate maximum temperature for each month
    monthly_maxes = numpy.zeros((12,) + temps.shape[1:])
    for month in range(1, 13):
        if month in month_indices:
            indices = month_indices[month]
            monthly_maxes[month-1] = numpy.nanmax(temps[indices], axis=0)

    # Calculate long-term mean of monthly maximums
    mean_monthly_maxes = numpy.nanmean(monthly_maxes, axis=1)

    # Find warmest month's mean maximum temperature
    thome = numpy.nanmax(mean_monthly_maxes, axis=0)

    return thome

def write_output(nc_in: Dataset, output_file: str, thome: numpy.ndarray):
    """Write Thome values to output NetCDF file."""
    log_information(f"Writing output to {output_file}")

    # Create output file
    nc_out = Dataset(output_file, 'w', format=NC_FORMAT)

    # Copy lat/lon dimensions and variables
    lat_dim = get_dim_with_attr(nc_in, 'standard_name', 'latitude',
                               'Failed to find latitude dimension')
    lon_dim = get_dim_with_attr(nc_in, 'standard_name', 'longitude',
                               'Failed to find longitude dimension')

    # Create dimensions
    nc_out.createDimension(lat_dim.name, len(lat_dim))
    nc_out.createDimension(lon_dim.name, len(lon_dim))

    # Copy coordinate variables
    lat_var = get_dimension_variable(nc_in, lat_dim)
    lon_var = get_dimension_variable(nc_in, lon_dim)

    nc_out.createVariable(lat_dim.name, lat_var.dtype, (lat_dim.name,))
    nc_out.createVariable(lon_dim.name, lon_var.dtype, (lon_dim.name,))

    nc_out.variables[lat_dim.name][:] = lat_var[:]
    nc_out.variables[lon_dim.name][:] = lon_var[:]

    # Copy attributes
    copy_attributes(lat_var, nc_out.variables[lat_dim.name])
    copy_attributes(lon_var, nc_out.variables[lon_dim.name])

    # Create Thome variable
    thome_var = nc_out.createVariable('thome', 'f8', (lat_dim.name, lon_dim.name))
    thome_var[:] = thome

    # Add metadata
    thome_var.long_name = 'Long-term mean maximum temperature of warmest month'
    thome_var.units = 'K'  # Assuming input is in Kelvin, as per CF convention
    thome_var.standard_name = 'air_temperature'  # Using same standard_name as input

    # Add global attributes
    nc_out.description = 'Thome (long-term mean maximum temperature of warmest month)'
    nc_out.history = f'Created {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'
    nc_out.source = f'Calculated from {os.path.basename(nc_in.filepath())}'

    nc_out.close()

def main(opts: Options):
    """
    Main program.

    @param opts: Parsed CLI options.
    """
    log_information(f"Processing {opts.input_file}")

    # Open input file.
    with open_netcdf(opts.input_file) as nc_in:
        # Calculate Thome.
        log_information("Calculating Thome...")
        thome = calculate_thome(nc_in, temp_data, times)

        # Write output.
        log_information("Writing output...")
        write_output(nc_in, opts.output_file, thome)
        log_information("Done!")

if __name__ == "__main__":
	# Parse CLI args
	opts = parse_args(argv)

	set_log_level(opts.log_level)
	set_show_progress(opts.show_progress)

	try:
		main(opts)
	except BaseException as error:
		# Basic error handling.
		log_error(format_exc())
		exit(1)
