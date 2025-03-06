#!/usr/bin/env python3
#
# Set the contents of the specified variable to the values given in a plaintext
# file.
#
# The main use case here is correcting the erroneous rlon values in the narclim2
# dataset.
#

from argparse import ArgumentParser
from sys import argv, exit
from shutil import copy2
from traceback import print_exc
from netCDF4 import Dataset

class Options:
    """
    CLI Options.
    """
    def __init__(self, in_file: str, out_file: str, values_file: str, var: str):
        """
        Constructor.

        Args:
            in_file: Input file.
            out_file: Output file.
            values_file: File containing the values to set the variable to.
            var: Name of the variable to set.
        """
        self.in_file = in_file
        self.out_file = out_file
        self.values_file = values_file
        self.var = var

def parse_args(argv: list[str]) -> Options:
    """
    Parse the command line arguments.

    Args:
        argv: The command line arguments.

    Returns:
        The parsed options.
    """
    parser = ArgumentParser(prog=argv[0], description="Set the contents of the "
        "specified variable in a netcdf file to the values given in a plaintext"
        " file.", exit_on_error=True)
    parser.add_argument(
        "-i", "--in-file", required=True, type=str, help="Input file.")
    parser.add_argument(
        "-o", "--out-file", type=str, help="Output file. If not specified, the "
        "input file will be modified in place.")
    parser.add_argument(
        "-f", "--values-file", required=True, type=str, help="File containing "
        "the values to set the variable to, one per line.")
    parser.add_argument(
        "-v", "--var", required=True, type=str, help="Name of the variable to "
        "set.")
    parsed = parser.parse_args(argv[1:])
    return Options(**vars(parsed))

def read_values(filename: str) -> list[float]:
    """
    Read the values for the variable from a text file, one per line.

    Args:
        filename: The filename of the text file containing the values.

    Returns:
        The values in the file, as a list.
    """
    with open(filename, "r") as f:
        return [float(line) for line in f]

def main(opts: Options):
    """
    Main function.
    """
    # If no output file was specified, modify the input file in-place.
    if opts.out_file is None:
        opts.out_file = opts.in_file

    # Copy input file to output file and then edit the output file in-place.
    # If the user specifies the same file path, don't copy, and just edit it
    # in-place.
    if opts.in_file != opts.out_file:
        copy2(opts.in_file, opts.out_file)

    values = read_values(opts.values_file)

    # Edit the output file in-place.
    with Dataset(opts.out_file, "r+", format="NETCDF4") as nc:
        var = nc.variables[opts.var]
        if var.size != len(values):
            raise ValueError(f"Variable '{opts.var}' in file '{opts.out_file}'"
                f" has {var.size} values, but {len(values)} were given in file"
                f" '{opts.values_file}'.")
        var[:] = values

if __name__ == "__main__":
    opts = parse_args(argv)
    try:
        main(opts)
    except BaseException as error:
        print_exc()
        exit(1)
