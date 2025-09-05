#!/usr/bin/env python3
#
# A wrapper script around netcdf2csv.py which converts all AWRA files to CSV.
#
import argparse, glob, netcdf2csv, os
from ozflux_logging import *
from sys import argv, exit

# AWRA is daily, so no need to include timestamp in the time column.
DATE_FMT = "%Y-%m-%d"

class AwraOptions:
    def __init__(self,
                 in_dir: str,
                 out_dir: str,
                 variables: list[str],
                 compress_output: bool,
                 lat_column: str,
                 lon_column: str,
                 time_column: str,
                 max_chunk_size: int,
                 log_level: LogLevel,
                 show_progress: bool):
        self.in_dir = in_dir
        self.out_dir = out_dir
        self.variables = variables
        self.compress_output = compress_output
        self.lat_column = lat_column
        self.lon_column = lon_column
        self.time_column = time_column
        self.max_chunk_size = max_chunk_size
        self.log_level = log_level
        self.show_progress = show_progress

def parse_args(argv: list[str]) -> AwraOptions:
    parser = argparse.ArgumentParser()
    parser.add_argument("--in-dir", "-i", required=True, help="Input directory")
    parser.add_argument("--out-dir", "-o", required=True, help="Output directory")
    parser.add_argument("--variables", "-v", action="append", default=[], help="Variables to extract from the NetCDF file")
    parser.add_argument("--compress-output", "-c", action="store_true", default=False, help="Compress output")
    parser.add_argument("--lat-column", "-lat", default="Lat", help="Latitude column name in the output file")
    parser.add_argument("--lon-column", "-lon", default="Lon", help="Longitude column name in the output file")
    parser.add_argument("--time-column", "-time", default="Time", help="Time column name in the output file")
    parser.add_argument("--max-chunk-size", "-m", type=int, default=netcdf2csv.DEFAULT_MAX_CHUNK_SIZE, help=f"Maximum chunk size (default: {netcdf2csv.DEFAULT_MAX_CHUNK_SIZE})")
    parser.add_argument("--log-level", "-l", type=int, default=LogLevel.WARNING, choices=[LogLevel.NONE, LogLevel.ERROR, LogLevel.WARNING, LogLevel.INFORMATION, LogLevel.DIAGNOSTIC, LogLevel.DEBUG], help="Log level")
    parser.add_argument("--show-progress", "-p", action="store_true", default=False, help="Show progress")
    parsed = parser.parse_args(argv)
    return AwraOptions(parsed.in_dir,
                       parsed.out_dir,
                       parsed.variables,
                       parsed.compress_output,
                       parsed.lat_column,
                       parsed.lon_column,
                       parsed.time_column,
                       parsed.max_chunk_size,
                       parsed.log_level,
                       parsed.show_progress)

def main(opts: AwraOptions):
    log_information(f"Converting AWRA files from {opts.in_dir} to {opts.out_dir}")
    step_start = 0
    big_step_size = 1 / len(opts.variables)
    for var in opts.variables:
        var_dir = os.path.join(opts.in_dir, var)
        if not os.path.exists(var_dir):
            log_error(f"Variable directory {var_dir} does not exist")
            return

        # Create output directory if it doesn't already exist.
        var_out_dir = os.path.join(opts.out_dir, var)
        if not os.path.exists(var_out_dir):
            os.makedirs(var_out_dir)

        in_files = glob.glob(os.path.join(var_dir, "*.nc"))

        # Sort files by basename.
        in_files.sort(key=os.path.basename)

        step_size = big_step_size / len(in_files)
        for in_file in in_files:
            basename = os.path.splitext(os.path.basename(in_file))[0]
            filename = f"{basename}.csv"
            out_file = os.path.join(var_out_dir, filename)
            log_diagnostic(f"Converting: {in_file}...")
            ncopts = netcdf2csv.Options(in_file,
                                        out_file,
                                        var,
                                        opts.compress_output,
                                        opts.lat_column,
                                        opts.lon_column,
                                        opts.time_column,
                                        opts.max_chunk_size,
                                        DATE_FMT,
                                        opts.log_level,
                                        opts.show_progress)
            netcdf2csv.main(ncopts, lambda p: log_progress(step_start + step_size * p))
            step_start += step_size

if __name__ == "__main__":
    opts = parse_args(argv[1:])
    set_log_level(opts.log_level)
    set_show_progress(opts.show_progress)
    main(opts)
