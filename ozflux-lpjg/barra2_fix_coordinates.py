#!/usr/bin/env python3
#
# Read longitude/latitude values from input .nc file, and overwrite the values
# in the output file. Need to assert that files have identical grid sizes.

import argparse, sys, netCDF4
from ozflux_logging import *
from ozflux_netcdf import *

class Options:
    """
    Command-line options.

    Attributes:
        input_file (str): Input file path.
        output_file (str): Output file path.
        log_level (LogLevel): Log level.
        dry_run (bool): Whether to perform a dry run (do not overwrite output file).
        warnings_as_errors (bool): Whether to treat all warnings as errors.
    """
    def __init__(self, input_file: str, output_file: str, log_level: LogLevel,
                 dry_run: bool, warnings_as_errors: bool):
        self.input_file = input_file
        self.output_file = output_file
        self.log_level = log_level
        self.dry_run = dry_run
        self.warnings_as_errors = warnings_as_errors

def parse_args(argv: list[str]) -> Options:
    parser = argparse.ArgumentParser(description="Read longitude/latitude values from input .nc file, and overwrite the values in the output file.")
    parser.add_argument("-i", "--in-file", type=str, help="Input file path.")
    parser.add_argument("-o", "--out-file", type=str, help="Output file path.")
    parser.add_argument("-v", "--verbosity", type=int, default=LogLevel.INFORMATION, help=f"Verbosity level (default={LogLevel.INFORMATION}).")
    parser.add_argument("-d", "--dry-run", action="store_true", help="Dry run (do not overwrite output file).")
    parser.add_argument("--warnings-as-errors", action="store_true", help="Treat all warnings as errors.")
    parsed = parser.parse_args(argv)
    return Options(parsed.in_file, parsed.out_file, parsed.verbosity,
                   parsed.dry_run, parsed.warnings_as_errors)

def main(opts: Options):
    with open_netcdf(opts.input_file) as nc_in:
        log_diagnostic("Parsing input file...")

        dim_lon_in = get_dim_from_std_name(nc_in, STD_LON)
        var_lon_in = get_dimension_variable(nc_in, dim_lon_in)

        log_information(f"Found input file longitude dimension: {dim_lon_in.name}")
        log_information(f"├── Shape: {var_lon_in.shape}")
        log_information(f"└── Units: {var_lon_in.units}")

        dim_lat_in = get_dim_from_std_name(nc_in, STD_LAT)
        var_lat_in = get_dimension_variable(nc_in, dim_lat_in)

        log_information(f"Found input file latitude dimension: {dim_lat_in.name}")
        log_information(f"├── Shape: {var_lat_in.shape}")
        log_information(f"└── Units: {var_lat_in.units}")

        with open_netcdf(opts.output_file, True) as nc_out:
            log_diagnostic("Parsing output file...")

            dim_lon_out = get_dim_from_std_name(nc_out, STD_LON)
            var_lon_out = get_dimension_variable(nc_out, dim_lon_out)

            log_information(f"Found output file longitude dimension: {dim_lon_out.name}")
            log_information(f"├── Shape: {var_lon_out.shape}")
            log_information(f"└── Units: {var_lon_out.units}")

            dim_lat_out = get_dim_from_std_name(nc_out, STD_LAT)
            var_lat_out = get_dimension_variable(nc_out, dim_lat_out)

            log_information(f"Found output file latitude dimension: {dim_lat_out.name}")
            log_information(f"├── Shape: {var_lat_out.shape}")
            log_information(f"└── Units: {var_lat_out.units}")

            if var_lon_in.shape != var_lon_out.shape or var_lat_in.shape != var_lat_out.shape:
                raise ValueError("Input and output files have different grid sizes.")

            # Check if values actually differ.
            if numpy.array_equal(var_lon_in[:], var_lon_out[:]):
                log_warning("Input and output files have identical longitude values.")
            if numpy.array_equal(var_lat_in[:], var_lat_out[:]):
                log_warning("Input and output files have identical latitude values.")

            if opts.dry_run:
                log_diagnostic("Dry run: would overwrite longitude/latitude values.")
                return

            log_diagnostic("Overwriting longitude/latitude values...")
            var_lon_out[:] = var_lon_in[:]
            var_lat_out[:] = var_lat_in[:]

            log_information("Done.")

if __name__ == "__main__":
    opts = parse_args(sys.argv[1:])

    set_log_level(opts.log_level)
    set_warnings_as_errors(opts.warnings_as_errors)

    main(opts)
