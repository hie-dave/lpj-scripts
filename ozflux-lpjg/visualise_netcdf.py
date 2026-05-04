#!/usr/bin/env python3
import argparse

import matplotlib.pyplot as plt
from ozflux_logging import *
from ozflux_netcdf import *

@dataclass(frozen=True)
class Options:
    file: str
    variable: str
    output: str | None
    width: int
    height: int
    log_level: LogLevel
    latitude: float
    longitude: float

@dataclass(frozen=True)
class DimensionOrder:
    index_lat: int
    index_lon: int
    index_time: int

def parse_args() -> Options:
    parser = argparse.ArgumentParser(description="Visualise a NetCDF file")
    parser.add_argument("file", help="Path to the NetCDF file")
    parser.add_argument("-v", "--variable", help="Variable to visualise")
    parser.add_argument("-o", "--output", help="Path to save the plot (optional)", default=None)
    parser.add_argument("--latitude", help="Latitude to slice at (default: first latitude)", type=float, default=numpy.nan)
    parser.add_argument("--longitude", help="Longitude to slice at (default: first longitude)", type=float, default=numpy.nan)
    parser.add_argument("--width", type=int, help="Width of the plot in pixels (default: 800)", default=800)
    parser.add_argument("--height", type=int, help="Height of the plot in pixels (default: 600)", default=600)
    parser.add_argument("--log-level", type=int, help="Logging level (0=ERROR, 1=WARNING, 2=INFO, 3=DIAG, 4=DEBUG)", default=int(LogLevel.INFORMATION))

    args = parser.parse_args()
    
    return Options(
        file=args.file,
        variable=args.variable,
        output=args.output,
        width=args.width,
        height=args.height,
        log_level=LogLevel(args.log_level),
        latitude=args.latitude,
        longitude=args.longitude,
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

    index_lat = index_of(dim_lat.name, var.dimensions)
    index_lon = index_of(dim_lon.name, var.dimensions)
    index_time = index_of(dim_time.name, var.dimensions)

    if index_lat == -1 or index_lon == -1 or index_time == -1:
        die(f"The following dimensions could not be found on variable "
            f"'{var.name}': {STD_LAT if index_lat == -1 else ''} "
            f"{STD_LON if index_lon == -1 else ''} "
            f"{STD_TIME if index_time == -1 else ''}")

    return DimensionOrder(index_lat=index_lat, index_lon=index_lon, index_time=index_time)

def read(nc, var, ilat: int, ilon: int):
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

    var_lat = get_var_from_std_name(nc, STD_LAT)
    var_lon = get_var_from_std_name(nc, STD_LON)

    lats = var_lat[:]
    lons = var_lon[:]

    ilat = 0
    if not numpy.isnan(opts.latitude):
        ilat = index_closest(opts.latitude, lats)
    ilon = 0
    if not numpy.isnan(opts.longitude):
        ilon = index_closest(opts.longitude, lons)

    data = read(nc, var, ilat, ilon)

    plt.figure(figsize=(opts.width / 100, opts.height / 100))
    plt.plot(data)
    plt.title(f"{opts.variable} ({lats[ilat]}, {lons[ilon]})")
    plt.xlabel("Time")
    plt.ylabel(f"{opts.variable} ({units})")
    plt.grid()

def main(opts: Options):
    with open_netcdf(opts.file) as nc:
        mkplot(nc, opts)
        if opts.output:
            plt.savefig(opts.output)
            log_information(f"Plot saved to {opts.output}")
        else:
            plt.show()

if __name__ == "__main__":
    opts = parse_args()
    set_log_level(opts.log_level)
    main(opts)

