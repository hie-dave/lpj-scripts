#!/usr/bin/env python3
#
# Convert a netcdf file to a csv file. The file must be 3-dimensional (latitude,
# longitude, time).
#

import argparse, gzip, math, numpy, os, pandas
from sys import argv
from ozflux_common import *
from ozflux_netcdf import *
from ozflux_logging import *
from netCDF4 import num2date
from typing import TextIO

DEFAULT_MAX_CHUNK_SIZE = 8192

class Options:
    """
    Options object for the netcdf2csv script.
    """
    def __init__(self,
                 input_file: str,
                 output_file: str,
                 variable: str,
                 compress_output: bool,
                 lat_column: str,
                 lon_column: str,
                 time_column: str,
                 max_chunk_size: int,
                 log_level: LogLevel,
                 show_progress: bool):
        self.input_file = input_file
        self.output_file = output_file
        self.variable = variable
        self.compress_output = compress_output
        self.lat_column = lat_column
        self.lon_column = lon_column
        self.time_column = time_column
        self.max_chunk_size = max_chunk_size
        self.log_level = log_level
        self.show_progress = show_progress

def parse_cli(argv: list[str]) -> Options:
    """
    Parse CLI arguments, return a parsed options object.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--in-file", "-i", required=True, help="Input NetCDF file")
    parser.add_argument("--out-file", "-o", required=True, help="Output CSV file")
    parser.add_argument("--variable", "-v", required=True, help="Variable to extract from the NetCDF file")
    parser.add_argument("--compress-output", "-c", action="store_true", help="Compress the output file (.gz will be appended to output file name if it's not already there)")
    parser.add_argument("--lat-column", default="Lat", help="Name of the latitude column in the output file")
    parser.add_argument("--lon-column", default="Lon", help="Name of the longitude column in the output file")
    parser.add_argument("--time-column", default="Time", help="Name of the time column in the output file")
    parser.add_argument("--max-chunk-size", type=int, default=DEFAULT_MAX_CHUNK_SIZE, help="Maximum chunk size in number of elements")
    parser.add_argument("--verbosity", type=int, default=LogLevel.WARNING, choices=[LogLevel.NONE, LogLevel.ERROR, LogLevel.WARNING, LogLevel.INFORMATION, LogLevel.DIAGNOSTIC, LogLevel.DEBUG], help="Log level")
    parser.add_argument("--show-progress", action="store_true", default=False, help="Show progress")
    parsed = parser.parse_args(argv)
    return Options(parsed.in_file,
                   parsed.out_file,
                   parsed.variable,
                   parsed.compress_output,
                   parsed.lat_column,
                   parsed.lon_column,
                   parsed.time_column,
                   parsed.max_chunk_size,
                   parsed.verbosity,
                   parsed.show_progress)

def choose_chunk_sizes(nc: Dataset, var: Variable, max_chunk_size: int) -> list[int]:
    """
    Choose chunk sizes for a variable.

    If the variable is chunked, we'll use the chunk sizes as a starting point.
    Otherwise, we'll use the dimension sizes.

    We'll then try to find a set of chunk sizes that will result in a total
    chunk size of at most max_chunk_size while optimising read performance.
    """
    chunk_sizes = var.chunking()
    if (chunk_sizes is None):
        chunk_sizes = [min(max_chunk_size, nc.dimensions[dim].size) for dim in var.dimensions]
    else:
        # Don't allow any chunk size to exceed the size of its dimension. Not
        # sure if the NetCDF API even allows for the creation of such a file.
        chunk_sizes = [min(nc.dimensions[dim].size, c) for dim, c in zip(var.dimensions, chunk_sizes)]

    for i in range(len(chunk_sizes)):
        total_chunk_size = numpy.prod(chunk_sizes)
        if total_chunk_size <= max_chunk_size:
            break
        if i < len(chunk_sizes) - 1:
            # First N dimensions set to chunksize=1
            chunk_sizes[i] = 1
        else:
            # Last dimension set to chunksize=8192
            chunk_sizes[i] = min(var.dimensions[i].size, max_chunk_size)
    return chunk_sizes

def copy_to_csv(nc: Dataset, file: TextIO, opts: Options, progress_callback: Callable[[float], None]):
    """
    Copy a variable from a NetCDF file to a CSV file.
    """
    # Get the dimension variables.
    lat = get_var_from_std_name(nc, STD_LAT)
    lon = get_var_from_std_name(nc, STD_LON)
    time = get_var_from_std_name(nc, STD_TIME)

    # Get the variable to extract.
    if not opts.variable in nc.variables:
        raise Exception(f"Variable '{opts.variable}' not found in the NetCDF file.")
    var = nc.variables[opts.variable]

    # Copy the variable to the output file.
    # Get dimensions
    dimensions = var.dimensions
    shape = var.shape

    # Determine if we need to handle a 3D variable (lat, lon, time) or something else
    if len(dimensions) != 3:
        raise Exception(f"Variable '{opts.variable}' has {len(dimensions)} dimensions, but only 3D variables (lat, lon, time) are supported.")

    time_index = index_of_throw(DIM_TIME, dimensions, lambda: f"Variable {opts.variable} has no time dimension")
    lat_index = index_of_throw(DIM_LAT, dimensions, lambda: f"Variable {opts.variable} has no latitude dimension")
    lon_index = index_of_throw(DIM_LON, dimensions, lambda: f"Variable {opts.variable} has no longitude dimension")

    # Get dimension values
    lat_values = lat[:]
    lon_values = lon[:]
    time_values = time[:]

    # Convert time values to datetime objects if they're numeric
    if hasattr(time, 'units'):
        try:
            time_dates = num2date(time_values, time.units, time.calendar if hasattr(time, 'calendar') else 'gregorian')
        except:
            log_warning("Could not convert time values to datetime objects. Using raw values.")
            time_dates = time_values
    else:
        time_dates = time_values

    dimension_values = [[], [], []]
    dimension_values[lat_index] = lat_values
    dimension_values[lon_index] = lon_values
    dimension_values[time_index] = time_dates

    # Determine chunk sizes for processing.
    chunk_sizes = choose_chunk_sizes(nc, var, opts.max_chunk_size)

    # Determine number of iterations over each dimension.
    niter = [math.ceil(s / c) for (s, c) in zip(shape, chunk_sizes)]

    # Total number of iterations.
    it_max = niter[0] * niter[1] * niter[2]

    # Current iteration.
    it = 0

    # Process data in chunks to avoid memory issues
    for i in range(niter[0]):
        ilow = i * chunk_sizes[0]
        ihigh = min(shape[0], ilow + chunk_sizes[0])
        ir = range(ilow, ihigh)
        for j in range(niter[1]):
            jlow = j * chunk_sizes[1]
            jhigh = min(shape[1], jlow + chunk_sizes[1])
            jr = range(jlow, jhigh)
            for k in range(niter[2]):
                klow = k * chunk_sizes[2]
                khigh = min(shape[2], klow + chunk_sizes[2])
                kr = range(klow, khigh)

                # Read this chunk from the input file.
                # NetCDF4 library automatically applies scale_factor and
                # add_offset attributes when reading data, so we get the
                # actual values, not the raw packed values.
                chunk = var[ir, jr, kr]

                # Write progress update.
                it += 1
                progress_callback(it / it_max)

                # Create rows for this chunk.
                rows = []

                # Get non-NaN values and their indices
                valid_mask = ~numpy.isnan(chunk)
                if numpy.any(valid_mask):
                    # Get the indices of non-NaN values
                    ii, jj, kk = numpy.nonzero(valid_mask)
                    values = chunk[valid_mask]

                    # Convert ranges to numpy arrays for array indexing
                    ir_array = numpy.array(list(ir))
                    jr_array = numpy.array(list(jr))
                    kr_array = numpy.array(list(kr))

                    # Create arrays for the actual indices in the full dataset
                    i_indices = ir_array[ii]
                    j_indices = jr_array[jj]
                    k_indices = kr_array[kk]

                    # Create arrays for each dimension's indices
                    indices = numpy.zeros((len(values), 3), dtype=int)
                    indices[:, 0] = i_indices
                    indices[:, 1] = j_indices
                    indices[:, 2] = k_indices

                    # Get the dimension values using vectorized operations
                    lat_vals = dimension_values[lat_index][indices[:, lat_index]]
                    lon_vals = dimension_values[lon_index][indices[:, lon_index]]
                    time_vals = dimension_values[time_index][indices[:, time_index]]

                    # Create a dictionary of rows directly
                    rows = [
                        {
                            opts.lat_column: lat_vals[i],
                            opts.lon_column: lon_vals[i],
                            opts.time_column: time_vals[i],
                            opts.variable: values[i]
                        }
                        for i in range(len(values))
                    ]

                # Convert to DataFrame and add to the list of chunks to write
                if rows:
                    df = pandas.DataFrame(rows)
                    # Write this chunk to the CSV file
                    df.to_csv(file, index=False, header=False)

def main(opts: Options, progress_callback: Callable[[float], None]):
    """
    Main function for the netcdf2csv script.
    """
    # Open the NetCDF file
    with open_netcdf(opts.input_file) as nc:
        # Determine if we need to compress the output.
        output_file = opts.output_file
        if opts.compress_output and not output_file.endswith(".gz"):
            output_file += ".gz"

        # Create the output directory if it doesn't exist.
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)

        # Open the output file for writing.
        open_func = gzip.open if output_file.endswith(".gz") else open

        # Write the header row first
        with open_func(output_file, "w") as f:
            header = f"{opts.lat_column},{opts.lon_column},{opts.time_column},{opts.variable}\n"
            f.write(header.encode() if output_file.endswith(".gz") else header)

            copy_to_csv(nc, f, opts, progress_callback)

        log_information(f"Successfully wrote data to {output_file}")

if __name__ == "__main__":
    try:
        opts = parse_cli(argv[1:])
        set_log_level(opts.log_level)
        set_show_progress(opts.show_progress)
        main(opts, log_progress)
    except Exception as e:
        print(f"Error: {e}")
        exit(1)
