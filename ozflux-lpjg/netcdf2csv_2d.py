#!/usr/bin/env python3
#
# Convert a netcdf file to a csv file. The file must be 2-dimensional (latitude,
# longitude).
#

import argparse, gzip, math, numpy, os, pandas
from sys import argv
from ozflux_common import *
from ozflux_netcdf import *
from ozflux_logging import *
from netCDF4 import num2date
from typing import TextIO

_DEFAULT_MAX_CHUNK_SIZE_BYTES = 8192

class Options:
    """
    Options object for the netcdf2csv script.

    Attributes:
        input_file (str): Input NetCDF file
        output_file (str): Output CSV file
        variables (list[str]): Variables to extract from the NetCDF file
        compress_output (bool): Compress the output file (.gz will be appended to output file name if it's not already there)
        lat_column (str): Name of the latitude column in the output file
        lon_column (str): Name of the longitude column in the output file
        max_chunk_size (int): Maximum chunk size in bytes
        coord_decimals (int): Number of decimal places for latitude/longitude in CSV
        value_decimals (int): Number of decimal places for values in CSV
        log_level (LogLevel): Log level
        show_progress (bool): Show progress
        write_tsv (bool): Write tab-separated values instead of CSV
    """
    def __init__(self,
                 input_file: str,
                 output_file: str,
                 variables: list[str],
                 compress_output: bool,
                 lat_column: str,
                 lon_column: str,
                 max_chunk_size: int,
                 coord_decimals: int,
                 value_decimals: int,
                 log_level: LogLevel,
                 show_progress: bool,
                 write_tsv: bool):
        self.input_file = input_file
        self.output_file = output_file
        self.variables = variables
        self.compress_output = compress_output
        self.lat_column = lat_column
        self.lon_column = lon_column
        self.max_chunk_size = max_chunk_size
        self.coord_decimals = coord_decimals
        self.value_decimals = value_decimals
        self.log_level = log_level
        self.show_progress = show_progress
        self.write_tsv = write_tsv

def parse_cli(argv: list[str]) -> Options:
    """
    Parse CLI arguments, return a parsed options object.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--in-file", "-i", required=True, help="Input NetCDF file")
    parser.add_argument("--out-file", "-o", required=True, help="Output CSV file")
    parser.add_argument("variables", nargs="+", help="Variables to extract from the NetCDF file")
    parser.add_argument("--compress-output", "-c", action="store_true", help="Compress the output file (.gz will be appended to output file name if it's not already there)")
    parser.add_argument("--lat-column", default="Lat", help="Name of the latitude column in the output file")
    parser.add_argument("--lon-column", default="Lon", help="Name of the longitude column in the output file")
    parser.add_argument("--log-level", "-l", type=int, default=LogLevel.WARNING, choices=[LogLevel.NONE, LogLevel.ERROR, LogLevel.WARNING, LogLevel.INFORMATION, LogLevel.DIAGNOSTIC, LogLevel.DEBUG], help="Log level")
    parser.add_argument("--show-progress", "-p", action="store_true", default=False, help="Show progress")
    parser.add_argument("--coord-decimals", type=int, default=6, help="Number of decimal places for latitude/longitude in CSV (default: %(default)s)")
    parser.add_argument("--value-decimals", type=int, default=-1, help="Number of decimal places for values in CSV. Negative means no rounding (default: %(default)s)")
    parser.add_argument("--write-tsv", action="store_true", default=False, help="Write tab-separated values instead of CSV.")

    # Mutually exclusive group for setting max chunk size in bytes, KiB, etc.
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--max-chunk-size-bytes", type=int, default=_DEFAULT_MAX_CHUNK_SIZE_BYTES, help="Maximum chunk size in bytes (default: %(default)s)")
    group.add_argument("--max-chunk-size-kib", type=int, default=None, help=f"Maximum chunk size in KiB (default: {_DEFAULT_MAX_CHUNK_SIZE_BYTES} bytes)")
    group.add_argument("--max-chunk-size-mib", type=int, default=None, help=f"Maximum chunk size in MiB (default: {_DEFAULT_MAX_CHUNK_SIZE_BYTES} bytes)")
    group.add_argument("--max-chunk-size-gib", type=int, default=None, help=f"Maximum chunk size in GiB (default: {_DEFAULT_MAX_CHUNK_SIZE_BYTES} bytes)")

    parsed = parser.parse_args(argv)

    max_chunk_size = parsed.max_chunk_size_bytes
    if parsed.max_chunk_size_kib is not None:
        max_chunk_size = parsed.max_chunk_size_kib * 1024
    if parsed.max_chunk_size_mib is not None:
        max_chunk_size = parsed.max_chunk_size_mib * 1024 * 1024
    if parsed.max_chunk_size_gib is not None:
        max_chunk_size = parsed.max_chunk_size_gib * 1024 * 1024 * 1024

    return Options(parsed.in_file,
                   parsed.out_file,
                   parsed.variables,
                   parsed.compress_output,
                   parsed.lat_column,
                   parsed.lon_column,
                   max_chunk_size,
                   parsed.coord_decimals,
                   parsed.value_decimals,
                   parsed.log_level,
                   parsed.show_progress,
                   parsed.write_tsv)

def choose_chunk_sizes(nc: Dataset, var: Variable,
                       max_chunk_size_bytes: int) -> list[int]:
    """
    Choose chunk sizes for a variable.

    If the variable is chunked, we'll use the chunk sizes as a starting point.
    Otherwise, we'll use the dimension sizes.

    We'll then try to find a set of chunk sizes that will result in a total
    chunk size of at most max_chunk_size while optimising read performance.

    Args:
        nc (Dataset): NetCDF dataset
        var (Variable): Variable to choose chunk sizes for
        max_chunk_size_bytes (int): Maximum chunk size in bytes

    Returns:
        list[int]: Chunk sizes
    """
    chunk_sizes = var.chunking()
    # Get maximum chunk size in number of elements.
    max_chunk_size = max_chunk_size_bytes // var.dtype.itemsize
    if (chunk_sizes is None):
        chunk_sizes = [min(max_chunk_size, nc.dimensions[dim].size) for dim in var.dimensions]
    else:
        # Don't allow any chunk size to exceed the size of its dimension. Not
        # sure if the NetCDF API even allows for the creation of such a file.
        chunk_sizes = [min(nc.dimensions[dim].size, c) for dim, c in zip(var.dimensions, chunk_sizes)]

    for i in range(len(chunk_sizes)):
        total_chunk_size = numpy.prod(chunk_sizes)
        dim = nc.dimensions[var.dimensions[i]]
        if total_chunk_size <= max_chunk_size:
            break
        if i < len(chunk_sizes) - 1:
            # First N dimensions set to chunksize=1
            dest_chunk_size = 1
            if total_chunk_size / chunk_sizes[i] < max_chunk_size:
                # We can reduce this chunk to larger than 1 and still be under
                # the max chunk size.
                dest_chunk_size = int(max_chunk_size / (total_chunk_size / chunk_sizes[i]))
                dest_chunk_size = min(dest_chunk_size, dim.size)
                dest_chunk_size = max(1, dest_chunk_size)

            if chunk_sizes[i] != dest_chunk_size:
                chunk_bytes = total_chunk_size * var.dtype.itemsize
                log_diagnostic(f"Reducing chunk size for dimension {i} from {chunk_sizes[i]} to {dest_chunk_size}, because total chunk size of {total_chunk_size} exceeds max chunk size of {max_chunk_size}. Increasing your max chunk size (via --max-chunk-size) would likely improve performance, if the computer has enough memory ({chunk_bytes} bytes would be enough to hold the chunk size before this particular reduction in size).")
            chunk_sizes[i] = dest_chunk_size
        else:
            # Last dimension set to chunksize=8192
            chunk_sizes[i] = min(dim.size, max_chunk_size)
    return chunk_sizes

def copy_to_csv(nc: Dataset, file: TextIO, opts: Options,
                progress_callback: Callable[[float], None]):
    """
    Copy a variable from a NetCDF file to a CSV file.

    Args:
        nc (Dataset): NetCDF dataset
        file (TextIO): Output file
        opts (Options): Options object
        progress_callback (Callable[[float], None]): Progress callback
    """
    if len(opts.variables) == 0:
        raise Exception("No variables specified.")

    # Validate variables.
    dims_prev = None
    var_prev = None
    for v in opts.variables:
        if v not in nc.variables:
            raise Exception(f"Variable '{v}' not found in the NetCDF file.")

        dimensions = nc.variables[v].dimensions

        # Determine if we need to handle a 2D variable (lat, lon) or something else
        if len(dimensions) != 2:
            raise Exception(f"Variable '{v}' has {len(dimensions)} dimensions, but only 2D variables (lat, lon) are supported.")

        if dims_prev is not None:
            if dimensions != dims_prev:
                raise Exception(f"Variable '{v}' has different dimensions to variable '{var_prev}'.")
        dims_prev = dimensions
        var_prev = v

    # Get the coordinate variables.
    lat = find_var(nc, STD_LAT, LONG_LAT, UNITS_LAT)
    lon = find_var(nc, STD_LON, LONG_LON, UNITS_LON)

    # Copy the variable to the output file.
    # Get dimensions
    var = nc.variables[var_prev]
    dimensions = var.dimensions
    shape = var.shape

    if len(lat.dimensions) != 1:
        raise Exception(f"Variable '{lat}' has {len(lat.dimensions)} dimensions, but only 1D variables (lat) are supported.")

    if len(lon.dimensions) != 1:
        raise Exception(f"Variable '{lon}' has {len(lon.dimensions)} dimensions, but only 1D variables (lon) are supported.")

    dim_lat = lat.dimensions[0]
    dim_lon = lon.dimensions[0]

    lat_index = index_of_throw(dim_lat, dimensions, lambda: f"Variable has no latitude dimension")
    lon_index = index_of_throw(dim_lon, dimensions, lambda: f"Variable has no longitude dimension")

    # Get dimension values.
    log_diagnostic(f"Reading coordinates...")
    lat_values = lat[:]
    lon_values = lon[:]
    log_debug(f"Successfully read all coordinates")

    # Round coordinates to a human-friendly number of decimals for CSV output
    if opts.coord_decimals >= 0:
        if isinstance(lat_values, numpy.ndarray):
            lat_values = numpy.round(lat_values.astype(float), opts.coord_decimals)
        if isinstance(lon_values, numpy.ndarray):
            lon_values = numpy.round(lon_values.astype(float), opts.coord_decimals)

    # Determine chunk sizes for processing.
    chunk_sizes = choose_chunk_sizes(nc, var, opts.max_chunk_size)
    log_diagnostic(f"Using chunk sizes: {chunk_sizes}")

    # Determine number of iterations over each dimension.
    niter = [math.ceil(s / c) for (s, c) in zip(shape, chunk_sizes)]

    # Total number of iterations.
    it_max = niter[0] * niter[1]

    # Current iteration.
    it = 0

    # Process data in chunks to avoid memory issues
    log_diagnostic(f"Processing data in chunks. Total iterations: {it_max}")
    for i in range(niter[0]):
        ilow = i * chunk_sizes[0]
        ihigh = min(shape[0], ilow + chunk_sizes[0])
        ir = range(ilow, ihigh)
        for j in range(niter[1]):
            jlow = j * chunk_sizes[1]
            jhigh = min(shape[1], jlow + chunk_sizes[1])
            jr = range(jlow, jhigh)

            # Read chunks for all requested variables and build a combined valid mask (no NaNs across any variable).
            # NetCDF4 automatically applies scale_factor and add_offset when slicing.
            var_chunks = []
            combined_valid_mask = None
            for vname in opts.variables:
                v = nc.variables[vname]
                v_chunk = v[ir, jr]
                var_chunks.append((vname, v_chunk))

                v_valid = ~numpy.isnan(v_chunk)
                combined_valid_mask = v_valid if combined_valid_mask is None else (combined_valid_mask & v_valid)

            if combined_valid_mask is not None and numpy.any(combined_valid_mask):
                # Indices of elements valid for ALL variables
                ii, jj = numpy.nonzero(combined_valid_mask)

                # Convert ranges to numpy arrays for array indexing.
                ir_array = numpy.array(list(ir))
                jr_array = numpy.array(list(jr))

                # Create arrays for the actual indices in the full dataset.
                i_indices = ir_array[ii]
                j_indices = jr_array[jj]

                # Create arrays for each dimension's indices.
                indices = [i_indices, j_indices]

                nrows = len(ii)
                if nrows > 0:
                    # Build DataFrame columns
                    data = {
                        opts.lat_column: lat_values[indices[lat_index]],
                        opts.lon_column: lon_values[indices[lon_index]]
                    }

                    # Add one column per variable, ensuring rounding if requested
                    for vname, v_chunk in var_chunks:
                        values = v_chunk[ii, jj]
                        if opts.value_decimals >= 0 and isinstance(values, numpy.ndarray):
                            values = numpy.round(values.astype(float), opts.value_decimals)
                        data[vname] = values

                    df = pandas.DataFrame(data)

                    # Write this chunk to the CSV/TSV file
                    header = True if it == 0 else False
                    df.to_csv(file, index=False, header=header, sep="\t" if opts.write_tsv else ",")

            # Write progress update.
            it += 1
            progress_callback(it / it_max)

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
        os.makedirs(os.path.dirname(os.path.abspath(output_file)),
                    exist_ok=True)

        # Open the output file for writing.
        open_func = gzip.open if output_file.endswith(".gz") else open

        # Write the header row first
        with open_func(output_file, "w") as f:
            copy_to_csv(nc, f, opts, progress_callback)

        log_information(f"Successfully wrote data to {output_file}")

if __name__ == "__main__":
    try:
        opts = parse_cli(argv[1:])
        set_log_level(opts.log_level)
        set_show_progress(opts.show_progress)
        main(opts, log_progress)
    except Exception as e:
        print(f"Failed to process file: {e}")
        exit(1)
