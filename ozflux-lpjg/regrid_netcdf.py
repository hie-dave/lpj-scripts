#!/usr/bin/env python
#
# Regrid one or more netcdf files which are using a south_north/west_east
# coordinate system. This will produce output files over the same region
# with the same resolution.
#
# Run with --help CLI argument for additional information.
#

from argparse import ArgumentParser
import os, re, time, pandas
import ozflux_dave, ozflux_netcdf, ozflux_parallel
import traceback
from ozflux_logging import *
from netCDF4 import Dataset, date2num
from sys import argv
from typing import Callable
from ozflux_netcdf import *
from ozflux_common import *

SOUTH_NORTH_NAME = "south_north"
WEST_EAST_NAME = "west_east"
LON_NAME = "lon"
LAT_NAME = "lat"
TIME_NAME = "time"

PCB_INT = PROGRESS_CHUNK_SIZE

class Options:
    """
    Class for storing CLI arguments from the user.

    @param loglvl: Log level.
    @param report_progress: True to write progress messages, 0 otherwise.
    @param parallel: True to process files in parallel.
    @param files: Input files.
    @param out_dir: Output directory.
    @param compression_level: Compression quality for output file [0, 9]. 0 = none, 1 = fastest compression, largest filesize, 9 = slowest compression, smallest filesize.
    @param compression_type: Compression algorithm to be used (default 'zlib').
    """
    def __init__(self, loglvl : LogLevel, report_progress: bool, parallel: bool,
                 files: list[str], out_dir: str, compression_level: int
                 , compression_type: str):

        self.log_level = loglvl
        self.report_progress = report_progress
        self.parallel = parallel

        self.files = files
        self.out_dir = out_dir

        self.compression_level = compression_level
        self.compression_type = compression_type

def parse_args(argv: list[str]) -> Options:
    """
    Parse CLI arguments, return a parsed options object.

    @param argv: Raw CLI arguments.

    @return Parsed Options object.
    """
    parser = ArgumentParser(prog=argv[0], description = "Formatting ozflux data into a format suitable for consumption by LPJ-Guess")
    parser.add_argument("-v", "--verbosity", type = int, help = "Logging verbosity (1-5, default 3)", nargs = "?", const = LogLevel.INFORMATION, default = LogLevel.INFORMATION)
    parser.add_argument("files", nargs = "+", help = "Input .nc files to be processed")
    parser.add_argument("-o", "--out-dir", required = True, help = "Path to the output directory. Input files will be written with the same name to this directory.")
    parser.add_argument("-p", "--show-progress", action = "store_true", help = "Report progress")
    parser.add_argument("-P", "--parallel", action = "store_true", help = "Process files in parallel")
    parser.add_argument("--compression-level", type = int, default = -1, help = "Compression level 0-9. 0 means no compression. 9 means highest compression ratio but slowest performance. Default behaviour is to use the same compression level as in the input files.")
    parser.add_argument("--compression-type", default = "zlib", help = "Compression algorithm to be used (default 'zlib')")
    parser.add_argument("--version", action = "version", version = "%(prog)s " + VERSION)

    p = parser.parse_args(argv[1:])

    if p.compression_level > 9:
        raise ValueError(f"Compression level must be between 0-9, but was: {p.compression_level}")

    return Options(p.verbosity, p.show_progress, p.parallel
                   , p.files, p.out_dir, p.compression_level
                   , p.compression_type)

def get_outfile_name(file: str, opts: Options):
    """
    Get the output file path to which the specified input file will be written.

    @param file: The path to the file to be processed.
    @param opts: Parsed CLI arguments.
    """
    return os.path.join(opts.out_dir, os.path.basename(file))

def create_var(name: str, nc_in: Dataset, nc_out: Dataset, cl: int, ct: str):
    """
    Copy the specified coordinate variable into the output file. No data will be
    written, but the variable will be created and all metadata written.

    @param name: Name of the variable.
    @param nc_in: Input netcdf file.
    @param nc_out: Output netcdf file.
    @param cl: Compression level of the output file.
    @param ct: Compression type of the output file.
    """
    var_in = nc_in.variables[name]

    # These will be coordinate variables, with only a single dimension of
    # the same name.
    dims = [name]

    # Get the format of this variable.
    fmt = get_nc_datatype(var_in.datatype.__str__())

    # Compression level
    filters = var_in.filters()
    cl = opts.compression_level
    if cl < 0:
        if filters is not None and "complevel" in filters:
            cl = filters['complevel']
        else:
            cl = 0

    # Compression type
    ct = get_compression_type(filters)

    # Copy chunk sizes from input file.
    chunking = var_in.chunking()
    cs = None if chunking is None or chunking == "contiguous" else tuple(chunking)

    # Create the variable in the output file.
    create_var_if_not_exists(nc_out, name, fmt, dims, cl, ct, cs)

    # Copy all relevant metadata for this variable.
    log_debug("Copying attributes for variable {name}...")
    copy_attributes(var_in, nc_out.variables[name])

def get_output_dim(dim: str) -> str:
    """
    Get the name of the equivalent dimension in the output file.
    """
    if dim == WEST_EAST_NAME:
        return LON_NAME
    if dim == SOUTH_NORTH_NAME:
        return LAT_NAME
    return dim

def get_coord_name(name: str) -> str:
    """
    Get the name of the corresponding coordinate name for the dimension name.
    E.g. for "lon" this will return "west_east".

    @param name: Name of the dimension.
    """
    if name == LON_NAME:
        return WEST_EAST_NAME
    if name == LAT_NAME:
        return SOUTH_NORTH_NAME
    return name

def init_outfile(nc_in: Dataset, nc_out: Dataset, compression_level: int
                 , compression_type: str):
    """
    Initialise the output file by creating dimensions and variables as required.

    @param nc_in: Input netcdf file, opened for reading.
    @param nc_out: Output netcdf file, opened for writing.
    @param compression_level: Compression level in the output file.
    @param compression_type: Compression algorithm to be used in the output file.
    """
    log_information("Initialising output file...")
    log_diagnostic("Creating dimensions in output file...")

    # The input file contains 3 dimensions: west_east, north_south, time.
    # Time is the only dimension we want to copy directly into the output file.
    log_diagnostic(f"Copying time variable into output file...")
    create_dim_if_not_exists(nc_out, TIME_NAME, nc_in.dimensions[DIM_TIME].size)
    create_var(TIME_NAME, nc_in, nc_out, compression_level, compression_type)
    copy_1d(nc_in, nc_out, TIME_NAME, 1, lambda _: ...)
    if not hasattr(nc_out.variables[TIME_NAME], ATTR_CALENDAR):
        setattr(nc_out.variables[TIME_NAME], ATTR_CALENDAR, CALENDAR_GREGORIAN)

    # Create spatial (lon/lat) dimensions and variables in output file.
    log_diagnostic(f"Copying spatial variables into output file...")
    for name in [LON_NAME, LAT_NAME]:
        first = float(nc_in.variables[name][:].min())
        last = float(nc_in.variables[name][:].max())
        size = nc_in.variables[get_coord_name(name)].size
        values = numpy.linspace(first, last, size)
        create_dim_if_not_exists(nc_out, name, size)
        create_var(name, nc_in, nc_out, compression_level, compression_type)
        nc_out.variables[name][:] = values

    # Now we copy all other variables, except west_east and south_north.
    for name in nc_in.variables:
        # Skip this variable if it's one of the coordinates we've already
        # created, or if it's in the list of variables to be skipped.
        if name in nc_out.variables or name in nc_in.dimensions:
            continue

        # Output dimensions will be the same as the input dimensions, except
        # west_east will become lon and south_north will become lat.
        var_in = nc_in.variables[name]
        dims = tuple(map(get_output_dim, var_in.dimensions))

        # Get the format of this variable.
        fmt = get_nc_datatype(var_in.datatype.__str__())

        # Compression level
        filters = var_in.filters()
        cl = opts.compression_level
        if cl < 0:
            if filters is not None and "complevel" in filters:
                cl = filters['complevel']
            else:
                cl = 0

        # Compression type
        ct = get_compression_type(filters)

        # Copy chunking from input file.
        chunking = var_in.chunking()
        cs = None if chunking is None or chunking == "contiguous" else tuple(chunking)
        create_var_if_not_exists(nc_out, name, fmt, dims, cl, ct, cs)
        copy_attributes(var_in, nc_out.variables[name])

    log_diagnostic("Copying global attributes...")
    copy_attributes(nc_in, nc_out)

def copy_var_3d(nc_in: Dataset, nc_out: Dataset, name: str, min_chunk_size: int
        , pcb: Callable[[float], None], check_names = True):
    """
    Copy the contents of the specified variable in the input file into the
    output file.

    @param nc_in: The input .nc file.
    @param nc_out: The output .nc file.
    @param name: The name of the variable to be copied.
    @param min_chunk_size: Minimum chunk size used when copying data.
    @param pcb: Progress callback function.
    """
    var_in = nc_in.variables[name]
    var_out = nc_out.variables[name]

    dim_names_in = var_in.dimensions
    dim_names_out = var_out.dimensions
    nd = len(dim_names_in)

    if nd != len(dim_names_out):
        m = "Inconsistent dimensionality in variable '%s': input file has %d dimensions (%s), output file has %d dimensions (%s)"
        raise ValueError(m % (name, nd, ", ".join(dim_names_in), len(dim_names_out), ", ".join(dim_names_out)))

    dims = [nc_in.dimensions[name] for name in dim_names_in]

    chunk_sizes = var_in.chunking()
    if (chunk_sizes is None):
        chunk_sizes = [min(min_chunk_size, dim.size) for dim in dims]
    else:
        chunk_sizes = [max(min_chunk_size * c, c) for c in chunk_sizes]

    shape = var_in.shape
    niter = [math.ceil(s / c) for (s, c) in zip(shape, chunk_sizes)]

    it_max = niter[0] * niter[1] * niter[2]
    it = 0
    for i in range(niter[0]):
        ilow = i * chunk_sizes[0]
        ihigh = min(shape[0], ilow + chunk_sizes[0])
        ir = range(ilow, ihigh)
        for j in range(niter[1]):
            jlow = j * chunk_sizes[1]
            jhigh = min(shape[1], jlow + chunk_sizes[1])
            jr = range(jlow, jhigh)
            for k in range(niter[2]):
                klow = k * chunk_sizes[2]
                khigh = min(shape[2], klow + chunk_sizes[2])
                kr = range(klow, khigh)
                var_out[ir, jr, kr] = var_in[ir, jr, kr]
                it += 1
                pcb(it / it_max)

def copy_var_data(nc_in: Dataset, nc_out: Dataset, name: str
    , min_chunk_size: int, pcb: Callable[[float], None]):

    # Call the appropriate function based on number of dimensions.
    ndim = len(nc_in.variables[name].dimensions)

    if ndim > 3:
        raise ValueError(f"Unable to copy variable {name}: >3 dimensions not implemented (variable {name} has {ndim} dimensions)")
    elif ndim == 3:
        copy_var_3d(nc_in, nc_out, name, min_chunk_size, pcb)
    else:
        raise ValueError(f"TBI: only 3d variables currently supported, but {name} has {ndim} variables!")
    # elif ndim == 2:
    #     copy_var_2d(nc_in, nc_out, name, min_chunk_size, pcb)
    # elif ndim == 1:
    #     copy_var_1d(nc_in, nc_out, name, min_chunk_size, pcb)

def copy_data(nc_in: Dataset, nc_out: Dataset, pcb: Callable[[float], None]):
    """
    Copy data from input file to output file, regridding as necessary.

    @param nc_in: Input netcdf file.
    @param nc_out: Output netcdf file.
    @param pcb: Progress callback function.
    """
    vars = [name for name in nc_out.variables if name not in nc_out.dimensions]
    total_weight = sum([nc_in.variables[name].size for name in vars])
    start = 0.0
    log_information("Copying variable data...")
    for name in vars:
        if name in nc_out.dimensions:
            continue

        log_diagnostic(f"Copying {name} data...")

        var = nc_out.variables[name]
        wgt = var.size / total_weight

        copy_var_data(nc_in, nc_out, name, 1024, lambda p: pcb(start + p * wgt))
        start += wgt

def process_file(infile: str, outfile: str, compression_level: int
                 , compression_type: str, pcb: Callable[[float], None]):
    """
    Process the specified input file.

    @param infile: Path to the input file to be processed.
    @param outfile: Path to which output will be written.
    @param compression_level: Compression level in the output file.
    @param compression_type: Compression algorithm to be used in the output file.
    @param pcb: Progress callback function.
    """
    log_information(f"Processing {infile}...")
    with open_netcdf(infile) as nc_in:
        with open_netcdf(outfile, True) as nc_out:            
            init_outfile(nc_in, nc_out, compression_level, compression_type)
            copy_data(nc_in, nc_out, pcb)

class PrashantRemappingTask(ozflux_parallel.Task):
    """
    Represents a task which will process one input file, regridding it from
    north_south/east_west to a regular latlon grid.
    """
    def __init__(self, infile: str, outfile: str, compression_level: int
          , compression_type: str):
        self.infile = infile
        self.outfile = outfile
        self.compression_level = compression_level
        self.compression_type = compression_type

    def exec(self, pcb: Callable[[float], None]):
        try:
            process_file(self.infile, self.outfile
                         , self.compression_level, self.compression_type, pcb)
        except BaseException as err:
            raise ValueError(f"Failed to process {self.infile}") from err


def main(opts: Options):
    """
    Main CLI entrypoint function.

    @param opts: Parsed CLI arguments.
    """
    # Create output directory if it doesn't exist.
    if not os.path.exists(opts.out_dir):
        os.makedirs(opts.out_dir)

    job_manager = ozflux_parallel.JobManager()

    for infile in opts.files:
        # Generate an output file path for this input file.
        outfile = get_outfile_name(infile, opts)

        # Clobber an existing output file, if one is found.
        if os.path.exists(outfile):
            os.remove(outfile)

        # Create a processing job for this file and add it to the queue.
        job = PrashantRemappingTask(infile, outfile
                                    , opts.compression_level
                                    , opts.compression_type)
        job_manager.add_job(job, os.path.getsize(infile))

    # Run all jobs as specified by user.
    if opts.parallel:
        job_manager.run_parallel()
    else:
        job_manager.run_single_threaded()

if __name__ == "__main__":
    # Parse CLI args
    opts = parse_args(argv)

    set_log_level(opts.log_level)
    set_show_progress(opts.report_progress)

    try:
        # Actual logic is in main().
        main(opts)
    except BaseException as error:
        # Basic error handling.
        log_error(traceback.format_exc())
        exit(1)
