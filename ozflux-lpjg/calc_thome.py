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
from calendar import month_name

from ozflux_common import *
from ozflux_logging import *
from ozflux_netcdf import *

# Standard name for air temperature variable.
STD_AIR_TEMP = "air_temperature"

# Name of the Thome variable in the NetCDF output file.
NAME_THOME = "thome"

# Long name of the Thome variable in the NetCDF output file.
LONG_NAME_THOME = "Long-term mean maximum temperature of warmest month"

# Name of the latitude column in the output file, when writing in CSV format.
COL_LAT = "latitude"

# Name of the longitude column in the output file, when writing in CSV format.
COL_LON = "longitude"

# Name of the Thome column in the output file, when writing in CSV format.
COL_THOME = "thome"

class OutputFormat(enum.Enum):
    """Output file format."""
    CSV = "csv"
    NETCDF = "nc"

def parse_output_format(s: str) -> OutputFormat:
    if s == OutputFormat.CSV.value:
        return OutputFormat.CSV
    elif s == OutputFormat.NETCDF.value:
        return OutputFormat.NETCDF
    else:
        raise ValueError(f"Unknown output format: {s}")

class Options:
    """Command line options for calc_thome.py."""

    def __init__(self, input_file: str, output_file: str, log_level: LogLevel,
                 show_progress: bool, format: OutputFormat):
        """Initialize with default values."""
        self.input_file = input_file
        self.output_file = output_file
        self.log_level = log_level
        self.show_progress = show_progress
        self.format = format

class DataPoint:
    """
    Represents a Thome value for a single gridcell.
    """
    def __init__(self, thome: float, latitude: float, longitude: float):
        self.thome = thome
        self.latitude = latitude
        self.longitude = longitude

def parse_args(argv: list[str]) -> Options:
    """Parse command line arguments and return Options object."""
    parser = argparse.ArgumentParser(description = "Calculate Thome from air temperature data")
    parser.add_argument("-v", "--verbosity", type = int, help = f"Logging verbosity ({LogLevel.ERROR}-{LogLevel.DEBUG}, default: {LogLevel.INFORMATION})", default = LogLevel.INFORMATION)
    parser.add_argument("-p", "--show-progress", action = "store_true", help = "Report progress")
    parser.add_argument("-i", "--in-file", help = "Input NetCDF file containing air temperature data")
    parser.add_argument("-o", "--out-file", help = "Output NetCDF file to write Thome values")
    parser.add_argument("-f", "--format", default = "csv", help = "Output file format (csv|nc)")

    args = parser.parse_args(argv[1:])
    opts = Options(args.in_file, args.out_file, args.verbosity,
                   args.show_progress, parse_output_format(args.format))
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

def calculate_thome_gridcell(temps: numpy.ndarray, times: list[datetime.datetime]) -> float:
    """
    Calculate Thome (long-term mean maximum temperature of warmest month) for a single gridcell.

    The input arrays should have the same shape, and should contain only the
    data for a single gridcell.

    @param temps: 1-dimensional numpy array of temperature data.
    @param times: List of datetime objects corresponding to each temperature data point.
    @return: The Thome value for the gridcell - ie the long-term mean maximum daily temperature of the warmest month.
    """
    log_diagnostic("Calculating monthly maximums...")

    # Convert times to month and day arrays for efficient grouping
    months = numpy.array([t.month for t in times])
    dates = numpy.array([t.date() for t in times])

    # Initialize array for monthly mean maximums
    max_temperatures = numpy.zeros(12)

    # Calculate for each month
    for month in range(1, 13):
        # Get data for this month
        month_mask = months == month
        month_temps = temps[month_mask]
        month_dates = dates[month_mask]

        if len(month_temps) == 0:
            log_warning(f"No data found for month {month_name[month]}")
            continue

        # Get daily maximums
        unique_dates = numpy.unique(month_dates)
        daily_maxes = numpy.array([
            numpy.max(month_temps[month_dates == date])
            for date in unique_dates
        ])

        # Calculate mean of daily maximums for this month
        max_temperatures[month - 1] = numpy.mean(daily_maxes)

    # Get the month with the highest mean maximum daily temperature
    return numpy.max(max_temperatures)

def calculate_thome(nc: Dataset) -> list[DataPoint]:
    """
    Calculate Thome (long-term mean maximum temperature of warmest month).

    @param nc: NetCDF file containing temperature data.
    @return: Numpy array of Thome values, with same dimensionality as the original temperature data.
    """
    # Find temperature variable
    temp_var = get_var_from_std_name(nc, STD_AIR_TEMP)
    log_information(f"Found temperature variable: {temp_var.name}")
    log_information(f"├── Shape: {temp_var.shape}")
    log_information(f"└── Units: {temp_var.units}")

    # Get time information.
    time_dim, time_var, times = get_time_info(nc)
    log_information(f"Found time dimension: {time_dim.name}")
    log_information(f"└── Time coverage: {times[0]} to {times[-1]} in {len(times)} timesteps")

    dim_lat = get_dim_from_std_name(nc, STD_LAT)
    var_lat = get_dimension_variable(nc, dim_lat)
    log_information(f"Found latitude dimension: {dim_lat.name}")
    log_information(f"├── Shape: {var_lat.shape}")
    log_information(f"└── Units: {var_lat.units}")

    dim_lon = get_dim_from_std_name(nc, STD_LON)
    var_lon = get_dimension_variable(nc, dim_lon)
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
        lat = 0.0
        lon = 0.0
        # Try reading coordinates from global attributes.
        if hasattr(nc, "latitude"):
            lat = float(nc.latitude)
            log_diagnostic(f"Using latitude {lat} from global attribute")
        if hasattr(nc, "longitude"):
            lon = float(nc.longitude)
            log_diagnostic(f"Using longitude {lon} from global attribute")
        return [DataPoint(calculate_thome_gridcell(temp_data, times), lat, lon)]

    # Get the index of the latitude and longitude dimensions in the temperature
    # variable's dimensions.
    index_latitude = index_of_throw(dim_lat.name, temp_var.dimensions,
    lambda: "Latitude dimension not found in temperature variable's dimensions")
    index_longitude = index_of_throw(dim_lon.name, temp_var.dimensions,
    lambda: "Longitude dimension not found in temperature variable's dimensions")
    index_time = index_of_throw(time_dim.name, temp_var.dimensions,
    lambda: "Time dimension not found in temperature variable's dimensions")

    # Iterate through gridcells.
    data: list[DataPoint] = []
    for i in range(temp_data.shape[index_latitude]):
        for j in range(temp_data.shape[index_longitude]):
            thome = calculate_thome_gridcell(temp_data[i, j], times)
            lat = float(var_lat[i])
            lon = float(var_lon[j])
            log_diagnostic(f"thome at ({lon}, {lat}): {thome}")
            data.append(DataPoint(thome, lat, lon))

    return data

def write_netcdf(nc_in: Dataset, out_file: str, data: list[DataPoint]):
    log_information(f"Writing output to {out_file} in NetCDF format")

    if os.path.exists(out_file):
        log_information(f"Clobbering existing output file: '{out_file}'")
        os.remove(out_file)

    # Create output file
    with open_netcdf(out_file, True) as nc_out:
        # Get coordinate dimensions and variables.
        lat_dim = get_dim_from_std_name(nc_in, STD_LAT)
        lon_dim = get_dim_from_std_name(nc_in, STD_LON)

        lat_var = get_dimension_variable(nc_in, lat_dim)
        lon_var = get_dimension_variable(nc_in, lon_dim)

        # Create dimensions and variables in output file.
        create_dim_if_not_exists(nc_out, lat_dim.name, len(lat_dim))
        create_dim_if_not_exists(nc_out, lon_dim.name, len(lon_dim))

        lat_type = get_nc_datatype(str(lat_var.dtype))
        create_var_if_not_exists(nc_out, lat_var.name, lat_type, (lat_dim.name,))

        lon_type = get_nc_datatype(str(lon_var.dtype))
        create_var_if_not_exists(nc_out, lon_var.name, lon_type, (lon_dim.name,))

        # Copy coordinates.
        nc_out.variables[lat_var.name][:] = lat_var[:]
        nc_out.variables[lon_var.name][:] = lon_var[:]

        # Copy attributes
        copy_attributes(lat_var, nc_out.variables[lat_var.name])
        copy_attributes(lon_var, nc_out.variables[lon_var.name])

        # Create Thome variable
        create_var_if_not_exists(nc_out, NAME_THOME, FORMAT_FLOAT, (lat_dim.name, lon_dim.name))
        thome_var = nc_out.variables[NAME_THOME]
        for i in range(lat_var.shape[0]):
            for j in range(lon_var.shape[0]):
                lat = float(lat_var[i])
                lon = float(lon_var[j])
                # Find DataPoint in data with matching latitude and longitude.
                values = [d.thome for d in data if floats_equal(d.latitude, lat) and floats_equal(d.longitude, lon)]
                if len(values) < 1:
                    raise ValueError(f"No DataPoint found with latitude={lat} and longitude={lon}")
                if len(values) > 1:
                    raise ValueError(f"Multiple DataPoints found with latitude={lat} and longitude={lon}")
                thome_var[i, j] = values[0]

        # Add metadata
        thome_var.long_name = LONG_NAME_THOME
        thome_var.units = "degC"
        # Doesn't really make sense to add a standard_name to Thome.

        # Add global attributes
        nc_out.description = 'Thome (long-term mean maximum temperature of warmest month)'
        nc_out.history = f'Created {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'
        nc_out.source = f'Calculated from {os.path.basename(nc_in.filepath())}'

def write_csv(out_file: str, data: list[DataPoint]):
    """
    Write the Thome values to a CSV file.

    @param out_file: Output file path.
    @param data: Data to write.
    """
    log_information(f"Writing output to {out_file} in CSV format")
    with open(out_file, "w") as file:
        file.write(f"{COL_LAT},{COL_LON},{COL_THOME}\n")
        for datum in data:
            file.write(f"{datum.latitude},{datum.longitude},{datum.thome}\n")

def write_output(opts: Options, data: list[DataPoint]):
    """
    Write Thome values to an output file.

    @param data: Thome values to write.
    @param output_file: Output file to write to.
    """
    if opts.format == OutputFormat.CSV:
        write_csv(opts.output_file, data)
    elif opts.format == OutputFormat.NETCDF:
        with open_netcdf(opts.input_file) as nc_in:
            write_netcdf(nc_in, opts.output_file, data)

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
        thome = calculate_thome(nc_in)

        # Write output.
        log_information("Writing output...")
        write_output(opts, thome)
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
