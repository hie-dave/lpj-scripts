#!/usr/bin/env python3
#
# This script will add a time dimension with a single value to a netcdf file.
# A single variable will be modified to have the new time dimension as its first
# dimension. A time dimension must not already exist in the input file.
#
# Run with --help for usage instructions.
#

import datetime, netCDF4, os
from argparse import ArgumentParser
from sys import argv

# Format of date/time values.
TIME_FMT = r"%Y-%m-%d %H:%m:%S"
DATE_FMT = r"%Y-%m-%d"

# Default reference time used for the newly-created time variable.
REFERENCE_TIME = "1900-01-01"

# Calendar to be used for the time variable.
CALENDAR = "gregorian"

class Options:
    """
    Class for storing CLI arguments from the user.
    """
    def __init__(self, infile: str, outfile: str, timestamp: datetime.datetime,
            reftime: datetime.datetime, increment: str, variable_name: str):
        self.in_file = infile
        self.out_file = outfile
        self.timestamp = timestamp
        self.reftime = reftime
        self.increment = increment
        self.var = variable_name

def fail(message):
    print(message)
    quit(1)

def parse_args(argv: list[str]) -> Options:
    """
    Parse CLI arguments, return a parsed options object.

    @param argv: Raw CLI arguments.

    @return Parsed Options object.
    """
    parser = ArgumentParser(prog=argv[0], description = "Add a time dimension to a netcdf file with a single value")
    parser.add_argument("-i", "--in-file", required = True, help = "Path to the input file")
    parser.add_argument("-o", "--out-file", required = True, help = "Path to the output file")
    parser.add_argument("-t", "--time", required = True, help = "The single date/time value to add to the time dimension, in yyyy-MM-dd hh:mm:ss (or yyyy-MM-dd) format")
    parser.add_argument("-v", "--variable-name", required = True, help = "Name of the variable which needs to have a time dimension added")
    parser.add_argument("-r", "--reference-time", default = REFERENCE_TIME, help = "The reference time to be used (ie time units will be X since <reference time>). Should be yyyy-MM-dd format. Default: 1900-01-01, or January 1st on the century before the timestamp")
    parser.add_argument("-I", "--increment", default = "days", help = "Time increment to be used in the time units. Options: days, hours. Default: days")

    parsed = parser.parse_args(argv[1:])
    timestamp = datetime.datetime.strptime(parsed.time, DATE_FMT)
    reference_time = datetime.datetime.strptime(parsed.reference_time, DATE_FMT)

    if parsed.increment != "days" and parsed.increment != "hours":
        fail(f"Time increment must be 'days' or 'hours' but is: '{parsed.increment}'")

    return Options(parsed.in_file, parsed.out_file, timestamp, reference_time,
                   parsed.increment, parsed.variable_name)

def copy_attributes(var_in: netCDF4.Variable, var_out: netCDF4.Variable):
    """
    Copy all attributes from the input variable to the output variable.
    """
    var_out.setncatts({k: var_in.getncattr(k) for k in var_in.ncattrs()})

def main(opts: Options):
    """
    Main CLI entrypoint function.

    @param opts: Parsed CLI options.
    """
    if os.path.exists(opts.out_file):
        print(f"Output file already exists: {opts.out_file}")
    # Open the netcdf files.
    with netCDF4.Dataset(opts.in_file, "r") as nc_in:
        if not opts.var in nc_in.variables:
            fail(f"Variable '{opts.var}' does not exist in input file")

        with netCDF4.Dataset(opts.out_file, "w") as nc_out:
            # Copy the existing dimensions.
            for name, dimension in nc_in.dimensions.items():
                if name == "time":
                    fail(f"Input file already has a time dimension")
                nc_out.createDimension(name, None if dimension.isunlimited() else len(dimension))

            # Create a time dimension with length 1.
            nc_out.createDimension('time', 1)

            # Copy the existing variables which do not vary with time.
            for name, variable in nc_in.variables.items():
                if name != opts.var:
                    out_var = nc_out.createVariable(name, variable.datatype, variable.dimensions)
                    copy_attributes(variable, out_var)
                    out_var[:] = variable[:]

            # Create the time variable.
            time_var = nc_out.createVariable('time', 'f8', ('time',))
            time_var.units = f"{opts.increment} since {opts.reftime}"
            time_var.calendar = CALENDAR

            # Convert the timestamp to its numeric representation.
            time_point = netCDF4.date2num([opts.timestamp],
                                          units = time_var.units,
                                          calendar = time_var.calendar)

            # Write the encoded time value to the time variable.
            time_var[:] = time_point

            # Add a new data variable with the added time dimension.
            data_var = nc_in.variables[opts.var]
            dims = ("time", *data_var.dimensions)
            new_var = nc_out.createVariable(opts.var, data_var.datatype, dims)
            copy_attributes(data_var, new_var)
            
            # Copy the data with an added time dimension
            new_var[0, :, :] = data_var[:]
            copy_attributes(nc_in, nc_out)

if __name__ == "__main__":
    opts = parse_args(argv)
    main(opts)
