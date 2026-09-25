#!/usr/bin/env python3
#PBS -l ncpus=1
#PBS -l walltime=48:00:00
#PBS -l mem=4GB
#PBS -q normal
#PBS -l wd
#PBS -m abe
#PBS -M d.holzworth@westernsydney.edu.au
#PBS -l storage=scratch/pt17
#PBS -N extract_many_sites
#PBS -o extract_many_sites.log
#PBS -j oe

from netCDF4 import Dataset, Dimension, Variable
from os.path import exists, join, basename
from os import makedirs
from argparse import ArgumentParser
from dataclasses import dataclass
from typing import List, Optional
from pandas import read_csv, DataFrame

STD_NAME = "standard_name"

STD_LON = "longitude"
STD_LAT = "latitude"
STD_TIME = "time"

@dataclass(frozen = True)
class Options:
    grid_file: str
    input_files: List[str]
    output_directory: str
    latitude_column: Optional[str]
    longitude_column: Optional[str]
    compression_level: int
    compression_type: str

@dataclass(frozen = True)
class Site:
    latitude: float
    longitude: float

@dataclass(frozen = True)
class ResolvedSite:
    latitude: float
    longitude: float
    lat_index: int
    lon_index: int

@dataclass(frozen = True)
class DimensionIndices:
    time: int
    lat: int
    lon: int

def parse_args() -> Options:
    parser = ArgumentParser(description = "Extract a list of sites from a netcdf file")
    parser.add_argument("-g", "--grid-file", required = True, help = "Path to the grid file")
    parser.add_argument("-i", "--input-files", required = True, nargs = "+", help = "List of input files to extract sites from")
    parser.add_argument("-o", "--output-directory", required = True, help = "Directory to write output files to")
    parser.add_argument("--latitude-column", default = None, help = "Name of the column in the grid file which contains latitude values. Default: inferred from the grid file")
    parser.add_argument("--longitude-column", default = None, help = "Name of the column in the grid file which contains longitude values. Default: inferred from the grid file")
    parser.add_argument("--compression-level", type = int, default = 5, help = "Compression level for output netcdf files. Default: 5")
    parser.add_argument("--compression-type", default = "zlib", help = "Compression type for output netcdf files. Default: zlib")

    parsed = parser.parse_args()
    return Options(parsed.grid_file, parsed.input_files,
                   parsed.output_directory, parsed.latitude_column,
                   parsed.longitude_column, parsed.compression_level,
                   parsed.compression_type)

def infer_column_name(default: Optional[str],
                      df: DataFrame,
                      candidates: List[str]) -> str:
    if default is not None:
        if not default in df.columns:
            raise RuntimeError(f"Column {default} does not exist in data")
        return default

    for candidate in candidates:
        for col in df.columns:
            if candidate.lower() == col.lower():
                return col
    detail = str.join(", ", candidates)
    raise RuntimeError(f"Unable to infer column name from candidates: {detail}")

def resolve_dimension(nc: Dataset, std_name: str) -> Dimension:
    for var in nc.variables.values():
        if hasattr(var, STD_NAME) and getattr(var, STD_NAME) == std_name:
            for dim_name in var.dimensions:
                if dim_name == var.name:
                    # Coordinate variables should have one dimension with the
                    # same name as the variable.
                    return nc.dimensions[dim_name]
    raise RuntimeError(f"Unable to resolve dimension with standard_name: {std_name}")

def resolve_variable(nc: Dataset, std_name: str) -> Variable:
    for var in nc.variables.values():
        if hasattr(var, STD_NAME) and getattr(var, STD_NAME) == std_name:
            return var
    raise RuntimeError(f"Unable to resolve variable with standard_name: {std_name}")

def parse_gridlist(path: str,
                   col_lon: Optional[str],
                   col_lat: Optional[str]) -> List[Site]:
    df = read_csv(path)
    col_lon = infer_column_name(col_lon, df, ["lon", "longitude"])
    col_lat = infer_column_name(col_lat, df, ["lat", "latitude"])
    return [Site(float(row[col_lat]), float(row[col_lon])) for row in df.iloc]

def chunk_sizes(var: Variable) -> Optional[List[int]]:
    """Return the chunk sizes of a variable, or None if it is contiguous."""
    if var.chunking() == "contiguous":
        return None
    return list(var.chunking())

def init_outfile(nc_in: Dataset,
                 nc_out: Dataset,
                 sites: List[Site],
                 compression_level: int = 5,
                 compression_type: str = "zlib") -> str:
    print(f"  Initialising output file...")

    # Create 3 dimensions. Bounds are ignored by LPJ-GUESS so don't copy them.
    dim_time = resolve_dimension(nc_in, STD_TIME)
    nc_out.createDimension(dim_time.name, len(dim_time))

    dim_lat = resolve_dimension(nc_in, STD_LAT)
    nc_out.createDimension(dim_lat.name, len(sites))

    dim_lon = resolve_dimension(nc_in, STD_LON)
    nc_out.createDimension(dim_lon.name, len(sites))

    compression = None if compression_level == 0 else compression_type

    # Create variables for these dimensions.
    var_time = resolve_variable(nc_in, STD_TIME)
    var_time_out = nc_out.createVariable(var_time.name, var_time.datatype,
                                         var_time.dimensions,
                                         compression=compression,
                                         complevel=compression_level,
                                         chunksizes=chunk_sizes(var_time))

    # Copy time variable metadata.
    # Note: ncattrs() only gets attributes which don't start with an underscore,
    # which is what we want.
    for attr in var_time.ncattrs():
        setattr(var_time_out, attr, getattr(var_time, attr))

    # Copy the time values.
    var_time_out[:] = var_time[:]

    var_lon = resolve_variable(nc_in, STD_LON)
    var_lon_out = nc_out.createVariable(var_lon.name, var_lon.datatype,
                                        var_lon.dimensions,
                                        compression=compression,
                                        complevel=compression_level,
                                        chunksizes=chunk_sizes(var_lon))
    for attr in var_lon.ncattrs():
        setattr(var_lon_out, attr, getattr(var_lon, attr))
    var_lon_out[:] = [site.longitude for site in sites]

    var_lat = resolve_variable(nc_in, STD_LAT)
    var_lat_out = nc_out.createVariable(var_lat.name, var_lat.datatype,
                                        var_lat.dimensions,
                                        compression=compression,
                                        complevel=compression_level,
                                        chunksizes=chunk_sizes(var_lat))
    for attr in var_lat.ncattrs():
        setattr(var_lat_out, attr, getattr(var_lat, attr))
    var_lat_out[:] = [site.latitude for site in sites]

    # The data variable will be a variable of lat, lon, and time.
    data_vars = [var for var in nc_in.variables.values()
                 if len(var.dimensions) == 3 and
                    dim_time.name in var.dimensions and
                    dim_lat.name in var.dimensions and
                    dim_lon.name in var.dimensions]
    if len(data_vars) == 0:
        raise RuntimeError("No data variables found with dimensions (time, lat, lon)")
    if len(data_vars) > 1:
        raise RuntimeError("Multiple data variables found with dimensions (time, lat, lon). "
                           "This is not supported.")

    var = data_vars[0]
    var_out = nc_out.createVariable(var.name, var.datatype,
                                    var.dimensions,
                                    compression=compression,
                                    complevel=compression_level,
                                    chunksizes=chunk_sizes(var))
    for attr in var.ncattrs():
        setattr(var_out, attr, getattr(var, attr))

    # Don't copy data yet.
    return var.name

def resolve_dimension_indices(nc: Dataset, var: Variable) -> DimensionIndices:
    if len(var.dimensions) != 3:
        raise RuntimeError(f"Variable {var.name} does not have 3 dimensions")

    time_index = None
    lat_index = None
    lon_index = None

    dim_time = resolve_dimension(nc, STD_TIME)
    dim_lat = resolve_dimension(nc, STD_LAT)
    dim_lon = resolve_dimension(nc, STD_LON)

    for i, dim_name in enumerate(var.dimensions):
        if dim_name == dim_time.name:
            time_index = i
        elif dim_name == dim_lat.name:
            lat_index = i
        elif dim_name == dim_lon.name:
            lon_index = i

    if time_index is None or lat_index is None or lon_index is None:
        raise RuntimeError(f"Variable {var.name} does not have dimensions (time, lat, lon)")

    return DimensionIndices(time=time_index, lat=lat_index, lon=lon_index)

def index_closest(needle: float, haystack: List[float]) -> int:
    """Return the index of the closest value in haystack to needle."""
    return min(range(len(haystack)), key=lambda i: abs(haystack[i] - needle))

def resolve_sites(nc: Dataset, sites: List[Site]) -> List[ResolvedSite]:
    var_lon = resolve_variable(nc, STD_LON)
    var_lat = resolve_variable(nc, STD_LAT)

    if len(var_lon.dimensions) != 1 or len(var_lat.dimensions) != 1:
        raise RuntimeError("Longitude and latitude variables must be 1D")

    lons = var_lon[:]
    lats = var_lat[:]

    resolved_sites = []
    for site in sites:
        ilon = index_closest(site.longitude, lons)
        ilat = index_closest(site.latitude, lats)
        resolved = ResolvedSite(latitude=site.latitude,
                                longitude=site.longitude, lon_index=ilon,
                                lat_index=ilat)
        resolved_sites.append(resolved)
    return resolved_sites

def data_ranges(var: Variable,
                indices: DimensionIndices,
                ilon: int,
                ilat: int) -> List[range]:
    ranges = [range(var.shape[i]) for i in range(len(var.dimensions))]

    # Time index stays as-is
    ranges[indices.lat] = range(ilat, ilat + 1)
    ranges[indices.lon] = range(ilon, ilon + 1)

    return ranges[0], ranges[1], ranges[2]

def copy_data(nc_in: Dataset, nc_out: Dataset, sites: List[Site], var_name: str):
    var_in = nc_in.variables[var_name]
    var_out = nc_out.variables[var_name]

    indices_in = resolve_dimension_indices(nc_in, var_in)
    indices_out = resolve_dimension_indices(nc_out, var_out)
    resolved_sites = resolve_sites(nc_in, sites)

    for i, site in enumerate(resolved_sites):
        print(f"  Copying data for site {i + 1}/{len(resolved_sites)}...")

        # Input access range uses the resolved indices.
        iq, jq, kq = data_ranges(var_in, indices_in, site.lon_index,
                                 site.lat_index)

        # Output access range is i because dimension order matches sites.
        ir, jr, kr = data_ranges(var_out, indices_out, i, i)

        # Copy data for this site.
        var_out[ir, jr, kr] = var_in[iq, jq, kq]

def _extract_sites(nc_in: Dataset, nc_out: Dataset, gridlist: List[Site],
                   compression_level: int, compression_type: str):
    var_name = init_outfile(nc_in, nc_out, gridlist, compression_level,
                            compression_type)
    copy_data(nc_in, nc_out, gridlist, var_name)

def extract_sites(input_file: str, output_file: str, gridlist: List[Site],
                  compression_level: int, compression_type: str):
    print(f"Processing {basename(input_file)}...")
    with Dataset(input_file) as nc_in:
        with Dataset(output_file, "w") as nc_out:
            _extract_sites(nc_in, nc_out, gridlist, compression_level,
                           compression_type)

def main(options: Options):
    if not exists(options.output_directory):
        makedirs(options.output_directory)

    gridlist = parse_gridlist(options.grid_file, options.longitude_column,
                              options.latitude_column)
    print(f"Parsed {len(gridlist)} sites from grid file")
    for file in options.input_files:
        out_file = join(options.output_directory, basename(file))
        extract_sites(file, out_file, gridlist, options.compression_level,
                      options.compression_type)

if __name__ == "__main__":
    options = parse_args()
    main(options)
