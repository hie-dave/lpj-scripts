#!/usr/bin/env python3
import argparse

import matplotlib.pyplot as plt
from ozflux_logging import *
from ozflux_netcdf import *

class Options:
    file: str
    variable: str
    output: str | None
    width: int
    height: int
    log_level: LogLevel

class DimensionOrder:
    index_lat: int
    index_lon: int
    index_time: int

def parse_args() -> Options:
    parser = argparse.ArgumentParser(description="Visualise a NetCDF file")
    parser.add_argument("file", help="Path to the NetCDF file")
    parser.add_argument("variable", help="Variable to visualise")
    parser.add_argument("-o", "--output", help="Path to save the plot (optional)", default=None)
    parser.add_argument("--width", type=int, help="Width of the plot in pixels (default: 800)", default=800)
    parser.add_argument("--height", type=int, help="Height of the plot in pixels (default: 600)", default=600)
    parser.add_argument("--log-level", type=int, help="Logging level (0=ERROR, 1=WARNING, 2=INFO, 3=DIAG, 4=DEBUG)", default=2)
    
    args = parser.parse_args()
    
    return Options(
        file=args.file,
        variable=args.variable,
        output=args.output,
        width=args.width,
        height=args.height,
        log_level=LogLevel(args.log_level)
    )

def die(msg: str):
    raise RuntimeError(msg)

def index_of(needle, haystack) -> int:
    for i in range(len(haystack)):
        if haystack[i] == needle:
            return i
    return -1

def index_closest(needle, haystack) -> int:
    closest_index = -1
    closest_distance = float('inf')
    for i in range(len(haystack)):
        distance = abs(haystack[i] - needle)
        if distance < closest_distance:
            closest_distance = distance
            closest_index = i
    return closest_index

def get_dim_order(nc, var) -> DimensionOrder:
    dim_lat = get_dim_from_std_name(nc, STD_LAT)
    dim_lon = get_dim_from_std_name(nc, STD_LON)
    dim_time = get_dim_from_std_name(nc, STD_TIME)

    index_lat = index_of(dim_lat, var.dimensions)
    index_lon = index_of(dim_lon, var.dimensions)
    index_time = index_of(dim_time, var.dimensions)

    if index_lat == -1 or index_lon == -1 or index_time == -1:
        die("Could not find latitude, longitude, or time dimensions in variable")

    return DimensionOrder(index_lat=index_lat, index_lon=index_lon, index_time=index_time)

def read(nc, var, opts: Options):
    var_lat = get_var_from_std_name(nc, STD_LAT)
    var_lon = get_var_from_std_name(nc, STD_LON)

    lats = var_lat[:]
    lons = var_lon[:]

    ilat = index_closest(0, lats)
    ilon = index_closest(0, lons)

    order = get_dim_order(nc, var)

    hyperslab = [slice(None)] * len(var.dimensions)
    hyperslab[order.index_lat] = ilat
    hyperslab[order.index_lon] = ilon
    data = var[tuple(hyperslab)]
    return data

def mkplot(nc, opts: Options):
    if opts.variable not in nc.variables:
        die(f"Variable '{opts.variable}' not found in NetCDF file")
    var = nc.variables[opts.variable]
    units = var.units if "units" in var.ncattrs() else "unknown units"
    data = read(nc, var, opts)

    plt.figure(figsize=(opts.width / 100, opts.height / 100))
    plt.plot(data)
    plt.title(opts.variable)
    plt.xlabel("Time")
    plt.ylabel(f"{opts.variable} ({units})")
    plt.grid()

def main(opts: Options):
    with open_netcdf(opts.file) as nc:
        plt = mkplot(nc, opts)

if __name__ == "__main__":
    opts = parse_args()
    set_log_level(opts.log_level)
    main(opts)

