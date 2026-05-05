#!/usr/bin/env python3
#
# Convert a netcdf file to a plain-text TSV file. The file must be 2-dimensional
# (latitude, longitude).
#

import argparse, gzip, math, numpy, os, pandas
from sys import argv
from ozflux_common import *
from ozflux_netcdf import *
from ozflux_logging import *
from typing import TextIO

DEFAULT_MAX_CHUNK_SIZE = 8192

class Options:
    """
    Options object for the netcdf2csv script.
    """
    def __init__(self,
                 input_files: str,
                 output_file: str,
                 variables: list[str],
                 lat_column: str,
                 lon_column: str,
                 max_chunk_size: int,
                 log_level: LogLevel):
        self.input_files = input_files
        self.output_file = output_file
        self.variables = variables
        self.lat_column = lat_column
        self.lon_column = lon_column
        self.max_chunk_size = max_chunk_size
        self.log_level = log_level

def parse_cli(argv: list[str]) -> Options:
    """
    Parse CLI arguments, return a parsed options object.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs = "+", help="Input NetCDF file")
    parser.add_argument("--out-file", "-o", required=True, help="Output text file")
    parser.add_argument("-v", "--variables", nargs = "+", required=True, help="Variable to extract from the NetCDF file")
    parser.add_argument("--lat-column", default="Lat", help="Name of the latitude column in the output file")
    parser.add_argument("--lon-column", default="Lon", help="Name of the longitude column in the output file")
    parser.add_argument("--max-chunk-size", type=int, default=DEFAULT_MAX_CHUNK_SIZE, help="Maximum chunk size in number of elements")
    parser.add_argument("--log-level", "-l", type=int, default=LogLevel.WARNING, choices=[LogLevel.NONE, LogLevel.ERROR, LogLevel.WARNING, LogLevel.INFORMATION, LogLevel.DIAGNOSTIC, LogLevel.DEBUG], help="Log level")
    parsed = parser.parse_args(argv)
    return Options(parsed.files,
                   parsed.out_file,
                   parsed.variables,
                   parsed.lat_column,
                   parsed.lon_column,
                   parsed.max_chunk_size,
                   parsed.log_level)

def read_file(in_file: str, variables: list[str]):
    with open_netcdf(in_file) as nc:
        dim_lat = get_dim_from_std_name(nc, STD_LAT)
        var_lat = get_dimension_variable(nc, dim_lat)
        log_diagnostic(f"Found latitude dimension: {dim_lat.name}")
        log_diagnostic(f"├── Shape: {var_lat.shape}")
        log_diagnostic(f"└── Units: {var_lat.units}")

        dim_lon = get_dim_from_std_name(nc, STD_LON)
        var_lon = get_dimension_variable(nc, dim_lon)
        log_diagnostic(f"Found longitude dimension: {dim_lon.name}")
        log_diagnostic(f"├── Shape: {var_lon.shape}")
        log_diagnostic(f"└── Units: {var_lon.units}")

        if len(var_lon) > 1:
            raise ValueError(f"Longitude dimension has length >1")
        if len(var_lat) > 1:
            raise ValueError(f"Latitude dimension has length >1")

        lat = float(var_lat[0])
        lon = float(var_lon[0])

        values = {
            "lon": lon,
            "lat": lat
        }
        for var in variables:
            if var not in nc.variables:
                raise ValueError(f"Input file '{in_file}' does not contain variable: '{var}'")
            ncvar = nc.variables[var]
            if len(ncvar) > 1:
                raise ValueError(f"Variable '{var}' has length >1: {len(ncvar)}")
            values[var] = ncvar[0]
        return values

def main(opts: Options):
    """
    Main CLI entrypoint function
    """
    # Generate DataFrame from all input files
    df = pandas.DataFrame()
    for in_file in opts.input_files:
        log_information(f"Processing {in_file}...")
        # Get dictionary mapping variable names to values.
        values = read_file(in_file, opts.variables)
        # Add a row to the DataFrame.
        df = pandas.concat([df, pandas.DataFrame([values])], ignore_index=True)
    # Write the DataFrame to the TSV output file.
    # df.to_csv(opts.out_file, index=False, sep="\t")
    with open(opts.out_file, "w") as f:
        f.write(df.to_string(index = False))

if __name__ == "__main__":
    opts = parse_cli(argv[1:])
    set_log_level(opts.log_level)
    main(opts)
