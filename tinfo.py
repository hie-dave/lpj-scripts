#!/usr/bin/env python3
import netCDF4
import numpy
import sys

STD_TIME = "time"

def var_from_std_name(nc: netCDF4.Dataset, std_name: str):
    """
    Retrieve a variable from a netCDF dataset based on its standard name.

    Parameters:
        nc (netCDF4.Dataset): The netCDF dataset object.
        std_name (str): The standard name of the variable to retrieve.

    Returns:
        netCDF4.Variable: The variable with the specified standard name, or None if not found.
    """
    for _, var in nc.variables.items():
        if getattr(var, "standard_name", None) == std_name:
            return var
    raise RuntimeError(f"No variable with standard name '{std_name}' found in the dataset.")

def print_tinfo(nc: netCDF4.Dataset):
    """
    Print information about the time variable in a netCDF file.

    Parameters:
        nc (netCDF4.Dataset): The netCDF dataset object.
    """

    time = var_from_std_name(nc, STD_TIME)

    if len(time.dimensions) != 1:
        raise RuntimeError("Time variable is not 1-dimensional.")

    if not hasattr(time, "units"):
        raise RuntimeError("Time variable does not have 'units' attribute.")

    calendar = getattr(time, "calendar", "standard") if hasattr(time, "calendar") else "standard"
    times = netCDF4.num2date(time[:], time.units, calendar)
    timesteps = numpy.diff(times)
    timestep = numpy.unique(timesteps)

    # Print time coverage and timestep.
    print(f"Time coverage: {times[0]} to {times[-1]}")
    if len(timestep) == 1:
        print(f"Timestep: {timestep[0]}")
    else:
        print(f"Timestep: {timestep} (multiple unique timesteps found)")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: {} <file>".format(sys.argv[0]))
        sys.exit(1)

    filename = sys.argv[1]
    with netCDF4.Dataset(filename, 'r') as nc:
        print_tinfo(nc)
