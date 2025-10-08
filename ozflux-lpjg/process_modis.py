#!/usr/bin/env python3
#
# Script to process MODIS data for use by the daveanalysis R package.
#
# Input: csv file containing MODIS data.
# Output: csv file containing QC-filtered, standardised MODIS data.
#

from argparse import ArgumentParser
from ozflux_logging import *
from ozflux_sites import *
from typing import Optional
from sys import argv, exit
from traceback import print_exc
import pandas, os

_MODIS_VERSION = "1.0.0"

# Name of the "site name" column in the LAI data file.
_LAI_COL_ID = "ID"

# Name of the date column in the LAI data file.
_LAI_COL_DATE = "Date"

# Name of the LAI column in the LAI data file.
_LAI_COL_LAI = "MYD15A2H_006_Lai_500m"

# Name of the quality control column in the LAI data file.
_LAI_COL_QCFLAG = "MYD15A2H_006_FparLai_QC"

# Name of the greenness column in the greenness input file.
_GREENNESS_COL_GREENNESS = "NDVI_smoothed"

# Name of the date column in the greenness input file.
_GREENNESS_COL_DATE = "xlDateTime"

# Name of the longitude column in the greenness input file.
_GREENNESS_COL_LON = "Lon"

# Name of the latitude column in the greenness input file.
_GREENNESS_COL_LAT = "Lat"

# Date format specifier for the greenness input file.
_GREENNESS_DATE_FMT = "%d/%m/%Y %H:%M"

# Name of the greenness variable in the output file.
_VAR_GREENNESS = "greenness"

# Units written to the greenness variable in the output file.
_GREENNESS_UNITS = "0-1"

# Name of the default sheet name in the LAI data file.
_SHEET_LAI = "OzFlux-sites-LAI-MYD15A2H-006-r"

# Minimum allowed LAI value (m2/m2).
_LAI_MIN = 0

# Maximum allowed LAI value (m2/m2).
_LAI_MAX = 20

class Options:
    """
    Class for storing CLI arguments from the user.

    @param in_file: Input file.
    @param out_file: Output file.
    @param log_level: Log level.
    @param normalise_site_names: Whether to normalise site names.
    @param lai_col: Name of the LAI column in the output file.
    @param site_col: Name of the site column in the output file.
    @param date_col: Name of the date column in the output file.
    """
    def __init__(self, in_file: str, out_file: str, log_level : LogLevel,
                 normalise_site_names: bool, lai_col: str, site_col: str,
                 date_col: str):
        self.log_level = log_level
        self.in_file = in_file
        self.out_file = out_file
        self.normalise_site_names = normalise_site_names
        self.lai_col = lai_col
        self.site_col = site_col
        self.date_col = date_col

def _modis_parse_args(argv: list[str]) -> Options:
    """
    Parse CLI arguments, return a parsed options object.

    @param argv: Raw CLI arguments.

    @return Parsed Options object.
    """
    parser = ArgumentParser(prog=argv[0], description = "Formatting ozflux data into a format suitable for consumption by LPJ-Guess")
    parser.add_argument("-i", "--input-file", required = True, help = "Input file path")
    parser.add_argument("-o", "--output-file", required = True, help = "Output file path")
    parser.add_argument("-v", "--verbosity", type = int, help = "Logging verbosity (1-5, default 3)", nargs = "?", const = LogLevel.INFORMATION, default = LogLevel.INFORMATION)
    parser.add_argument("-n", "--normalise-site-names", action = "store_true", help = "Normalise site names")
    parser.add_argument("--lai-col", type = str, help = "Name of the LAI column in the output file")
    parser.add_argument("--site-col", type = str, help = "Name of the site column in the output file")
    parser.add_argument("--date-col", type = str, help = "Name of the date column in the output file")
    parser.add_argument("--version", action = "version", version = "%(prog)s " + _MODIS_VERSION)

    p = parser.parse_args(argv[1:])
    return Options(p.input_file, p.output_file, p.verbosity,
                   p.normalise_site_names, p.lai_col, p.site_col, p.date_col)

def _modis_read_csv_or_excel(file_path: str, cols: list[str], date_col: str
                           , date_fmt: str) -> pandas.DataFrame:
    """
    Read an excel (.xlsx) or .csv file.
    """
    if not os.path.exists(file_path):
        log_error(f"File does not exist: {file_path}")
    ext = os.path.splitext(file_path)[1]
    log_diagnostic(f"Reading columns: {cols}")
    if ext.lower() == ".csv":
        log_diagnostic(f"Reading CSV file: {file_path}")
        data = pandas.read_csv(file_path, usecols = cols, parse_dates = False)
    else:
        log_diagnostic(f"Reading Excel file: {file_path}")
        data = pandas.read_excel(file_path, _SHEET_LAI, usecols = cols
                                 , parse_dates = False)
    log_diagnostic(f"Successfully read {len(data)} rows")
    log_diagnostic(f"Converting date column to datetime using: {date_fmt}")
    data[date_col] = pandas.to_datetime(data[date_col], format = date_fmt)
    return data

def _bit_set(x: int, n: int) -> bool:
    """
    Check if the N-th bit is set in x.

    @param x
    """
    y = 1 << n
    return x & y == y

def _modis_parse_qc_flag(qc: int) -> bool:
    """
    Check if a QC flag indicates usable data.

    @param qc: A QC flag value.
    """
    # See here for a complete technical guide to the product, which includes
    # documentation of the QC flag:
    # https://lpdaac.usgs.gov/documents/624/MOD15_User_Guide_V6.pdf

    # If bit 2 is set, dead detectors caused >50% adjacent detector retrieval.
    if _bit_set(qc, 2):
        return False

    # Meaning of bit 3 depends on bit 4.
    if _bit_set(qc, 3):
        # If bit 4 is also set, significant clouds were present. Otherwise,
        # cloud state not defined (ie unknown). Either way, we don't trust this
        # value.
        return False

    # If bit 4 set but bit 3 not set, mixed clouds were present. Treat this as
    # acceptable for now.

    # If bit 6 set, RT method failed.
    if _bit_set(qc, 6):
        return False

    # If bit 7 set, pixel not produced at all, value couldn't be retrieved.
    if _bit_set(qc, 7):
        return False

    return True

def _modis_parse_qc_flags(qcflags: pandas.Series) -> list[bool]:
    """
    Filter out invalid QC flags from the input data.

    @param qcflags: Series of QC flags.
    """
    return [_modis_parse_qc_flag(int(qc)) for qc in qcflags]

def modis_read_lai(lai_file: str) -> pandas.DataFrame:
    """
    Read LAI data from a CSV or Excel file.
    """
    log_information("Reading LAI data...")
    cols = [_LAI_COL_ID, _LAI_COL_DATE, _LAI_COL_LAI, _LAI_COL_QCFLAG]
    # _DATE_FMT = "%d/%m/%Y" # 28/07/2002
    _DATE_FMT = "%m/%d/%Y" # 07/28/2002
    return _modis_read_csv_or_excel(lai_file, cols, _LAI_COL_DATE, _DATE_FMT)

def modis_read_greenness(greenness_file: str) -> pandas.DataFrame:
    """
    Read greenness data from a CSV or Excel file.
    """
    log_information("Reading greenness data...")
    cols = [
        _GREENNESS_COL_DATE,
        _GREENNESS_COL_LON,
        _GREENNESS_COL_LAT,
        _GREENNESS_COL_GREENNESS
    ]
    return _modis_read_csv_or_excel(greenness_file, cols, _GREENNESS_COL_DATE,
                                   _GREENNESS_DATE_FMT)

def modis_get_site_data(data: pandas.DataFrame, site: str) -> pandas.DataFrame:
    """
    Get data for a specific site.
    """
    return data[data[_LAI_COL_ID] == site]

def modis_normalise(data: pandas.DataFrame) -> pandas.DataFrame:
    """
    Group data by date and calculate the mean for each group.
    """
    return data.groupby([_LAI_COL_DATE]).mean(numeric_only = True)

def modis_qcfilter(data: pandas.DataFrame) -> pandas.DataFrame:
    """
    Filter out invalid QC flags from the input data, and drop the QCFlag column.
    """
    nrow = len(data)
    data = data.iloc[_modis_parse_qc_flags(data[_LAI_COL_QCFLAG])]
    dropped = nrow - len(data)
    log_diagnostic(f"Filtered out {dropped} rows with invalid QC flags")
    return data.drop(columns = [_LAI_COL_QCFLAG])

def modis_guard(data: pandas.DataFrame) -> pandas.DataFrame:
    """
    Guard clause filtering the LAI data to be within bounds.
    """
    nrow = len(data)
    data = data[data[_LAI_COL_LAI].between(_LAI_MIN, _LAI_MAX)]
    dropped = nrow - len(data)
    log_diagnostic(f"Filtered out {dropped} rows with LAI out of bounds")
    return data

def modis_get_lai_col() -> str:
    """
    Get the name of the column containing LAI data.
    """
    return _LAI_COL_LAI

def modis_normalise_site_names(data: pandas.DataFrame) -> pandas.DataFrame:
    """
    Normalise site names.
    """
    col = _LAI_COL_ID
    # Apply normalisation to the site column in-place and return the full DataFrame
    data[col] = data[col].apply(normalise_site_name)
    return data

def modis_reorder(data: pandas.DataFrame) -> pandas.DataFrame:
    """
    Reorder the columns in the input data.
    """
    # Reorder columns.
    data = data[[_LAI_COL_DATE, _LAI_COL_ID, _LAI_COL_LAI]]
    # Reorder rows: sort by site, then by date.
    data = data.sort_values(by = [_LAI_COL_ID, _LAI_COL_DATE])
    return data

def modis_rename(data: pandas.DataFrame, lai_col: str, site_col: str,
                 date_col: str) -> pandas.DataFrame:
    """
    Rename columns as required.
    """
    if lai_col is not None:
        log_diagnostic(f"Renaming LAI column from {_LAI_COL_LAI} to {lai_col}")
        data = data.rename(columns = {_LAI_COL_LAI: lai_col})
    if site_col is not None:
        log_diagnostic(f"Renaming site column from {_LAI_COL_ID} to {site_col}")
        data = data.rename(columns = {_LAI_COL_ID: site_col})
    if date_col is not None:
        log_diagnostic(f"Renaming date column from {_LAI_COL_DATE} to {date_col}")
        data = data.rename(columns = {_LAI_COL_DATE: date_col})
    return data

def modis_write(data: pandas.DataFrame, out_file: str):
    """
    Write the data to a CSV file.
    """
    log_diagnostic(f"Writing to {out_file}")
    data.to_csv(out_file, index = False)

def main(opts: Options):
    """
    Main function.
    """
    log_information("Processing MODIS data...")

    # Read input data.
    data = modis_read_lai(opts.in_file)

    # Filter/normalise data.
    data = modis_qcfilter(data)
    data = modis_guard(data)
    if opts.normalise_site_names:
        data = modis_normalise_site_names(data)
    data = modis_reorder(data)
    data = modis_rename(data, opts.lai_col, opts.site_col, opts.date_col)

    # Write output data.
    modis_write(data, opts.out_file)
    log_information("Done.")

if __name__ == "__main__":
    try:
        options = _modis_parse_args(argv)
        set_log_level(options.log_level)
        set_show_progress(False) # Progress reporting not supported
        main(options)
    except Exception as e:
        log_error("Error: " + str(e))
        print_exc()
        exit(1)
