#!/usr/bin/env python
#
# Set the values of a variable in a NetCDF file.
#

from netCDF4 import Dataset
from shutil import copyfile
from sys import argv, exit
from argparse import ArgumentParser

class Options:
    """
    Class for storing CLI arguments from the user.

    @param infile: Path to the input file.
    @param outfile: Path to the output file.
    @param var: Name of the variable to set.
    @param values_file: Path to the file containing the values to set the variable to.
    """
    def __init__(self, infile: str, outfile: str, var: str, values_file: str):
        self.in_file = infile
        self.out_file = outfile
        self.var = var
        self.values_file = values_file

def parse_args(argv: list[str]) -> Options:
    """
    Parse CLI arguments, return a parsed options object.

    @param argv: Raw CLI arguments.

    @return Parsed Options object.
    """
    parser = ArgumentParser(prog=argv[0], description = "Add a time dimension to a netcdf file with a single value")
    parser.add_argument("-i", "--in-file", required = True, help = "Path to the input file")
    parser.add_argument("-o", "--out-file", required = True, help = "Path to the output file")
    parser.add_argument("-v", "--variable-name", required = True, help = "Name of the variable which needs to have a time dimension added")
    parser.add_argument("-f", "--values-file", required = True, help = "Path to the file containing the values to set the variable to")

    parsed = parser.parse_args(argv[1:])

    return Options(parsed.in_file, parsed.out_file, parsed.variable_name, parsed.values_file)

def read_values(filename: str) -> list[float]:
    """
    Read values from a text file.

    @param filename: Path to the file containing the values, one value per line.

    @return List of values.
    """
    with open(filename, "r") as f:
        return [float(line) for line in f]

def main(options: Options) -> None:
    """
    Main function.

    @param options: Parsed CLI options provided by the user..
    """

    # Read values from input file.
    values = read_values(options.values_file)

    # Copy input file to output file, then edit the output file in place.
    copyfile(options.in_file, options.out_file)

    with Dataset(options.out_file, "r+") as nc:
        var = nc.variables[options.var]

        # Check that the number of values matches the number of values in the variable.
        if len(values) != var.size:
            raise ValueError(f"Number of values ({len(values)}) does not match number of values in variable ({var.size})")

        # Set the values.
        var[:] = values

if __name__ == "__main__":
    options = parse_args(argv)
    main(options)
