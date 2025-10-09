#!/usr/bin/env python3
#
# Process GOSIF GPP data. Convert from input .tif.gz to csv format.
#

# Metadata:
# For each time step, GOSIF GPP consists of the mean and standard deviation (SD)
# of GPP; the mean and SD are based on the eight different sets of GPP
# estimates; mean is recommended for most analyses, and SD can be used to
# account for the uncertainty of GPP.
#
# Spatial resolution: 0.05 degree
# Spatial extent: globe
# Temporal resolution: 8 day (and monthly, annual)
# Temporal extent: 2000-2024
# File format: GeoTIFF
# Map projection: Geographic
# Scale factor: 0.001 (8-day GPP); 0.01 (monthly GPP); 0.1 (annual GPP)
# Fill values:
# - 65535 (water bodies)
# - 65534 (lands under snow/ice throughout the year)
# Units:
# - g C m-2 d-1 (8-day GPP)
# - g C m-2 mo-1 (monthly GPP)
# - g C m-2 yr-1 (annual GPP)

from argparse import ArgumentParser
from ozflux_logging import *
from ozflux_sites import *
from sys import argv
import datetime, glob, os, numpy, pandas, rasterio, re

_SCALE_FACTOR_8DAY = 0.001
_SCALE_FACTOR_MONTHLY = 0.01
_SCALE_FACTOR_YEARLY = 0.1

_MASK_WATER = 65535
_MASK_ICE = 65534

class Options:
    """
    Options for the process_gosif script.

    Attributes:
        in_dir (str): Input directory.
        out_file (str): Output file.
        sites (list[str]): List of sites to process.
        log_level (LogLevel): Log level.
        show_progress (bool): Whether to show progress.
    """
    def __init__(self, in_dir: str, out_file: str, sites: list[str],
                 date_col: str, site_col: str, value_col: str,
                 log_level: LogLevel, show_progress: bool):
        self.in_dir = in_dir
        self.out_file = out_file
        self.sites = sites
        self.date_col = date_col
        self.site_col = site_col
        self.value_col = value_col
        self.log_level = log_level
        self.show_progress = show_progress

def parse_args(args: list[str]) -> Options:
    """
    Parse command line arguments.
    """
    parser = ArgumentParser(description="Process GOSIF GPP data.")
    parser.add_argument("--in-dir", "-i", required = True, help = "Input directory")
    parser.add_argument("--out-file", "-o", required = True, help = "Output file")
    parser.add_argument("--date-col", "-d", default = "date", help = "Date column")
    parser.add_argument("--site-col", "-s", default = "site", help = "Site column")
    parser.add_argument("--value-col", "-v", default = "gpp", help = "Value column")
    parser.add_argument("--log-level", "-l", default = LogLevel.INFORMATION, help = "Log level")
    parser.add_argument("--show-progress", "-p", action = "store_true", help = "Show progress")
    p = parser.parse_args(args)

    # TODO: this could be configurable. But is it worth it? Syntactically it
    # would probably be annoying for the user. If we do ever make that change,
    # be sure to call normalise_site_name for each provided site name.
    sites = get_standard_sites()

    return Options(p.in_dir, p.out_file, sites, p.date_col, p.site_col,
                   p.value_col, int(p.log_level), p.show_progress)

def mkdir_p(path: str):
    """
    Make a directory, creating parent directories as needed.
    """
    log_diagnostic(f"Creating directory {path}")
    os.makedirs(path, exist_ok = True)

def _parse_gosif_date_from_name(path: str) -> pandas.Timestamp:
    # Example: GOSIF_GPP_2000057_Mean.tif.gz -> YYYY DDD
    filename = os.path.basename(path)
    m = re.search(r'GOSIF_GPP_(\d{4})(\d{3})_Mean\.tif(?:\.gz)?$', filename)
    if not m:
        # Fallback: try any 4+3 digit YYDDD at end
        m = re.search(r'(\d{4})(\d{3})', filename)
    if not m:
        raise ValueError(f"Cannot parse date from filename: {path}")
    y, doy = int(m.group(1)), int(m.group(2))
    return pandas.Timestamp(datetime.datetime(y, 1, 1)) + pandas.Timedelta(days=doy - 1)

def _read_data(in_file: str, date_col: str, site_col: str,
               value_col: str, site_names: list[str]) -> pandas.DataFrame:
    """
    Read data from a GOSIF .tif.gz file.
    """
    log_diagnostic(f"Reading {in_file}")

    # Use GDAL's VSI to read gzipped GeoTIFF transparently. This requires
    # building a /vsigzip/ path with an absolute filename. For an absolute
    # input like /path/file.tif.gz, the resulting path should be
    # /vsigzip//path/file.tif.gz (note the double slash after vsigzip).
    abs_path = os.path.abspath(in_file)
    vsi_path = f"/vsigzip/{abs_path}"
    log_debug(f"Using VSI path: {vsi_path}")

    # Ensure points in dataset CRS
    sites = [resolve_site(site) for site in site_names]
    pts_lon = numpy.array([s.lon for s in sites], dtype=float)
    pts_lat = numpy.array([s.lat for s in sites], dtype=float)

    with rasterio.open(vsi_path) as src:
        data = src.read(1)
        crs = src.crs
        transform = src.transform
        nodata = src.nodata
        src_crs = "EPSG:4326"
        if src.crs and str(src.crs) != src_crs:
            xs, ys = rasterio.warp.transform(src_crs, src.crs, pts_lon, pts_lat)
        else:
            xs, ys = pts_lon, pts_lat

        # rasterio.sample expects (x, y) in the dataset CRS
        samples = list(src.sample(zip(xs, ys)))
        # Single-band assumed; if multiple bands, pick band index you need
        values = [float(v[0]) if len(v) else np.nan for v in samples]

        # Handle nodata
        nodata = src.nodata
        if nodata is not None:
            values = [np.nan if v == nodata else v for v in values]

        # Mask documented fill values
        fill_values = {_MASK_WATER, _MASK_ICE}
        values = [numpy.nan if (v in fill_values) else v for v in values]

        # Apply scale factor.
        # FIXME: make this configurable somehow.
        scale = _SCALE_FACTOR_8DAY
        values = [numpy.nan if numpy.isnan(v) else (v * scale) for v in values]

    ts = _parse_gosif_date_from_name(in_file)
    rows = []

    for s, v in zip(sites, values):
        rows.append({
            date_col: ts,
            site_col: s.site,
            value_col: v
        })
    return pandas.DataFrame(rows)

def _filter_data(data: pandas.DataFrame, sites: list[str]) -> pandas.DataFrame:
    """
    Filter data to only include the specified sites.
    """
    log_diagnostic(f"Filtering data to sites {sites}")
    raise NotImplementedError("TODO: implement _filter_data")

def main(opts: Options):
    """
    Main function.
    """
    log_information("Processing GOSIF GPP data.")
    mkdir_p(os.path.dirname(opts.out_file))
    # Get all .tif.gz files in the input directory.
    in_files = glob.glob(os.path.join(opts.in_dir, "*.tif.gz"))
    # Order by filename.
    in_files.sort()
    data = pandas.DataFrame()
    i = 0
    for in_file in in_files:
        log_diagnostic(f"Processing {in_file}")
        df = _read_data(in_file, opts.date_col, opts.site_col, opts.value_col, opts.sites)
        data = pandas.concat([data, df], ignore_index = True)
        i += 1
        log_progress(1.0 * i / len(in_files))

    # TODO: if date is index, we will need to write the index.
    # Order by site, then date.
    data.sort_values(by = [opts.site_col, opts.date_col], inplace = True)
    data.to_csv(opts.out_file, index = False)

if __name__ == "__main__":
    options = parse_args(argv[1:])
    set_log_level(options.log_level)
    set_show_progress(options.show_progress)
    main(options)
