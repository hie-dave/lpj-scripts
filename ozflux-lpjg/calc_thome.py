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
from typing import Optional

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

    def __init__(
        self,
        input_files: list[str],
        output_file: str,
        log_level: LogLevel,
        show_progress: bool,
        format: OutputFormat,
        gridlist: Optional[str] = None,
    ):
        """Initialize with default values."""
        self.input_files = input_files
        self.output_file = output_file
        self.log_level = log_level
        self.show_progress = show_progress
        self.format = format
        self.gridlist = gridlist

class DataPoint:
    """
    Represents a Thome value for a single gridcell.
    """
    def __init__(self, thome: float, latitude: float, longitude: float):
        self.thome = thome
        self.latitude = latitude
        self.longitude = longitude

class Coordinate:
    """
    Represents a coordinate in a gridlist.
    """
    def __init__(self, lon: float, lat: float):
        self.lon = lon
        self.lat = lat
    def __str__(self):
        return f"({self.lon}, {self.lat})"

def parse_args(argv: list[str]) -> Options:
    """Parse command line arguments and return Options object."""
    parser = argparse.ArgumentParser(description = "Calculate Thome from air temperature data")
    parser.add_argument("-v", "--verbosity", type = int, help = f"Logging verbosity ({LogLevel.ERROR}-{LogLevel.DEBUG}, default: {LogLevel.INFORMATION})", default = LogLevel.INFORMATION)
    parser.add_argument("-p", "--show-progress", action = "store_true", help = "Report progress")
    parser.add_argument("-o", "--out-file", required = True, help = "Output NetCDF file to write Thome values")
    parser.add_argument("-f", "--format", default = "csv", help = "Output file format (csv|nc)")
    parser.add_argument("-g", "--gridlist", help = "Optional file containing longitude/latitude pairs to filter gridcells")
    parser.add_argument("files", nargs = "+", help = "Input NetCDF files containing air temperature data to be processed")

    args = parser.parse_args(argv[1:])
    opts = Options(args.files, args.out_file, args.verbosity,
                   args.show_progress, parse_output_format(args.format),
                   args.gridlist)
    return opts

def fail(msg):
    """
    Abort execution by throwing an exception with the specified message.
    """
    raise ValueError(msg)

def get_time_info(nc: Dataset) -> tuple[Dimension, Variable, list[datetime.datetime]]:
    """
    Get time dimension and variable information.

    @return: Time dimension, time variable, and time values.
    """
    # Find time dimension using standard_name.
    time_dim = get_dim_from_std_name(nc, STD_TIME)
    time_var = get_dimension_variable(nc, time_dim)

    # Read and parse time values.
    # Note: using cftime values because unusal calendars such as 365_day don't
    # play well with python datetime objects.
    times = num2date(time_var[:], time_var.units, time_var.calendar)
    return time_dim, time_var, times

def calculate_thome_gridcell(temps: numpy.ndarray, times: list[datetime.datetime], pcb: Callable[[float], None]) -> float:
    """
    Calculate Thome (long-term mean maximum temperature of warmest month) for a single gridcell.

    The input arrays should have the same shape, and should contain only the
    data for a single gridcell.

    @param temps: 1-dimensional numpy array of temperature data.
    @param times: List of datetime objects corresponding to each temperature data point.
    @param pcb: Callback function to report progress.
    @return: The Thome value for the gridcell - ie the long-term mean maximum daily temperature of the warmest month.
    """
    log_diagnostic("Calculating monthly maximums...")

    # Convert times to month and day arrays for efficient grouping
    months = numpy.array([t.month for t in times])
    # dates = numpy.array([t.date() for t in times])
    # cftime dates don't have a .date() method, so we use an integer encoding.
    dates = numpy.array([t.year * 10000 + t.month * 100 + t.day for t in times])

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

def calculate_thome(nc: Dataset, pcb: Callable[[float], None], gridlist: list[Coordinate]) -> list[DataPoint]:
    """
    Calculate Thome (long-term mean maximum temperature of warmest month).

    @param nc: NetCDF file containing temperature data.
    @param pcb: Callback function to report progress.
    @param gridlist: Optional path to file containing gridcells to process.
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

    # Convert to ℃ if necessary.
    conversion = None
    if temp_var.units != "degC":
        _conversion = find_units_conversion(temp_var.units, "degC")
        # Timestep is irrelevant here, so pass dummy value of 1.
        conversion = lambda x: _conversion(x, 1)
        log_information(f"Temperature will be converted from {temp_var.units} to ℃...")
    else:
        log_diagnostic("Temperature already in ℃.")

    if len(temp_var.dimensions) == 1:
        # Flatten temperature data to 1-dimensional array.
        # TODO: test this
        temp_data = temp_var[:].reshape(-1)
        if conversion is not None:
            temp_data = conversion(temp_data)

        lat = 0.0
        lon = 0.0
        # Try reading coordinates from global attributes.
        if hasattr(nc, "latitude"):
            lat = float(nc.latitude)
            log_diagnostic(f"Using latitude {lat} from global attribute")
        if hasattr(nc, "longitude"):
            lon = float(nc.longitude)
            log_diagnostic(f"Using longitude {lon} from global attribute")
        return [DataPoint(calculate_thome_gridcell(temp_data, times, pcb), lat, lon)]

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
    subset = gridlist is not None and len(gridlist) > 0
    nlon = 1 if subset else temp_var.shape[index_longitude]
    nlat = len(gridlist) if subset else temp_var.shape[index_latitude]

    lats = var_lat[:]
    lons = var_lon[:]

    step_start = 0.0
    step_size = 1 / (nlat * nlon)

    # Full slice for time dimension
    indices = [slice(None)] * 3

    for i in range(nlat):
        for j in range(nlon):
            indices[index_latitude] = index_of(lats, gridlist[i].lat) if subset else i
            # Note: j is always 0 for subset mode, so use i for gridlist lookup.
            indices[index_longitude] = index_of(lons, gridlist[i].lon) if subset else j
            temperature = temp_var[tuple(indices)]
            if conversion is not None:
                temperature = conversion(temperature)

            pcb(step_start)
            thome = calculate_thome_gridcell(temperature, times, lambda p: pcb(step_start + step_size * p))
            step_start += step_size
            lat = gridlist[i].lat if subset else lats[i]
            lon = gridlist[i].lon if subset else lons[j]
            log_diagnostic(f"thome at ({lon}, {lat}): {thome}")
            pcb(step_start)
            data.append(DataPoint(thome, lat, lon))

    pcb(1.0)
    return data

def write_netcdf(opts: Options, data: list[DataPoint]):
    log_information(f"Writing output to {opts.output_file} in NetCDF format")

    if os.path.exists(opts.output_file):
        log_information(f"Clobbering existing output file: '{opts.output_file}'")
        os.remove(opts.output_file)

    lats = numpy.unique(numpy.sort([d.latitude for d in data]))
    lons = numpy.unique(numpy.sort([d.longitude for d in data]))

    # Create output file
    with open_netcdf(opts.output_file, True) as nc_out:
        # Create coordinate dimensions.
        create_dim_if_not_exists(nc_out, DIM_LAT, len(lats))
        create_dim_if_not_exists(nc_out, DIM_LON, len(lons))

        # Create coordinate variables.
        create_var_if_not_exists(nc_out, DIM_LAT, FORMAT_FLOAT, (DIM_LAT,))
        create_var_if_not_exists(nc_out, DIM_LON, FORMAT_FLOAT, (DIM_LON,))

        # Write coordinates.
        nc_out.variables[DIM_LAT][:] = lats
        nc_out.variables[DIM_LON][:] = lons

        # Write coordinate metadata.
        setattr(nc_out.variables[DIM_LAT], ATTR_UNITS, UNITS_LAT)
        setattr(nc_out.variables[DIM_LAT], ATTR_STD_NAME, STD_LAT)
        setattr(nc_out.variables[DIM_LAT], ATTR_LONG_NAME, LONG_LAT)

        setattr(nc_out.variables[DIM_LON], ATTR_UNITS, UNITS_LON)
        setattr(nc_out.variables[DIM_LON], ATTR_STD_NAME, STD_LON)
        setattr(nc_out.variables[DIM_LON], ATTR_LONG_NAME, LONG_LON)

        # Create data variable.
        create_var_if_not_exists(nc_out, NAME_THOME, FORMAT_FLOAT, (DIM_LAT, DIM_LON))
        thome_var = nc_out.variables[NAME_THOME]

        # Write data.
        # Technically, the data should be correctly ordered, but just to be safe
        # we will still use a linear search on each iteration. This could be
        # optimised if needed by having calc_thome return a 2D array or similar.
        for i in range(len(data)):
            i = index_of_throw(data[i].latitude, lats, lambda: f"Failed to get latitude index for latitude {data[i].latitude}")
            j = index_of_throw(data[i].longitude, lons, lambda: f"Failed to get longitude index for longitude {data[i].longitude}")
            thome_var[i, j] = data[i].thome

        # Write meatadata to data variable.
        thome_var.long_name = LONG_NAME_THOME
        thome_var.units = "degC"
        # No standard name - t_home is not a CF-recognised variable.

        # Write global metadata.
        nc_out.description = 'Thome (long-term mean maximum temperature of warmest month)'
        nc_out.history = f'Created {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'
        nc_out.source = f'Calculated from {", ".join([os.path.basename(f) for f in opts.input_files])}'

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
        write_netcdf(opts, data)

def parse_gridlist(gridlist: str) -> list[Coordinate]:
    """
    Parse the gridlist file and return a list of Coordinate objects.
    @param gridlist: Path to the gridlist file.
    """
    with open(gridlist) as file:
        i = 1
        coords: list[Coordinate] = []
        for line in file.readlines():
            line = line.strip()
            parts = [p.strip() for p in line.split(" ")]
            if len(parts) != 2:
                fail(f"parse_gridlist(): line {i} contains {len(parts)} elements (should be 2: <lon> <lat>)")
            lon = float(parts[0])
            lat = float(parts[1])
            if lon < -180 or lon > 180:
                fail(f"parse_gridlist(): line {i} contains invalid longitude {lon}: must be in range [-180, 180]")
            if lat < -90 or lat > 90:
                fail(f"parse_gridlist(): line {i} contains invalid latitude {lat}: must be in range [-90, 90]")
            coords.append(Coordinate(lon, lat))
            i += 1
        return coords

def main(opts: Options):
    """
    Main program.

    @param opts: Parsed CLI options.
    """
    # Open input file.
    thome: list[DataPoint] = []
    step_start = 0.0
    step_size = 1 / len(opts.input_files)
    gridlist: list[Coordinate] = []
    if opts.gridlist is not None:
        gridlist = parse_gridlist(opts.gridlist)

    for in_file in opts.input_files:
        log_information(f"Processing {in_file}")
        with open_netcdf(in_file) as nc_in:
            # Calculate Thome.
            log_information("Calculating Thome...")
            thome.extend(calculate_thome(nc_in, lambda p: log_progress(step_start + step_size * p), gridlist))
            step_start += step_size
            log_progress(step_start)

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
