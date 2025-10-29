#!/usr/bin/env python3
#
# Script to process data from flux towers (C fluxes, soil moisture, etc) and
# produce per-variable CSV files containing timeseries of this data.
#
# Input: netcdf files containing MODIS data.
# Output: csv files, one per variable, containing timeseries data.
#

from argparse import ArgumentParser
from ozflux_logging import *
from ozflux_sites import *
from ozflux_netcdf import *
from typing import Optional, Callable
from sys import argv, exit
from traceback import print_exc
import pandas, os, numpy
from ozflux_common import find_units_conversion_opt

_VERSION = "0.1.0"

# Aggregation methods. These are applied to the variable in the *input units*.
_AGGREGATORS = {
    "ER_LL": lambda x: numpy.mean(x), # umol/m2/s
    "ER_LT": lambda x: numpy.mean(x), # umol/m2/s
    "ER_SOLO": lambda x: numpy.mean(x), # umol/m2/s
    "GPP_LL": lambda x: numpy.mean(x), # umol/m2/s
    "GPP_LT": lambda x: numpy.mean(x), # umol/m2/s
    "GPP_SOLO": lambda x: numpy.mean(x), # umol/m2/s
    "NEE_LL": lambda x: numpy.mean(x), # umol/m2/s
    "NEE_LT": lambda x: numpy.mean(x), # umol/m2/s
    "NEE_SOLO": lambda x: numpy.mean(x), # umol/m2/s
    "NEP_LL": lambda x: numpy.mean(x), # umol/m2/s
    "NEP_LT": lambda x: numpy.mean(x), # umol/m2/s
    "NEP_SOLO": lambda x: numpy.mean(x), # umol/m2/s
    "ET": lambda x: numpy.mean(x), # kg/m2/s
    "Ta": lambda x: numpy.mean(x), # ℃
    "Fsd": lambda x: numpy.mean(x), # W/m2
    "Precip": lambda x: numpy.sum(x), # mm
    "Sws": lambda x: numpy.mean(x), # m3/m3
}

class Options:
    """
    Options for the process_fluxes script.

    Attributes:
        in_files (list[str]): Input netcdf files.
        out_file (str): Output file.
        variable (str): Variable to process.
        out_units (str): Output units.
        annual (bool): Process annual data.
        site_col (str): Site column name.
        data_col (Optional[str]): Data column name.
        date_col (Optional[str]): Date column name.
        log_level (LogLevel): Log level.
        show_progress (bool): Show progress.
    """
    def __init__(self, in_files: list[str], out_file: str, variable: str,
                 out_units: str, annual: bool, site_col: str,
                 data_col: Optional[str], date_col: Optional[str],
                 log_level: LogLevel, show_progress: bool):
        self.in_files = in_files
        self.out_file = out_file
        self.annual = annual
        self.variable = variable
        self.out_units = out_units
        self.log_level = log_level
        self.show_progress = show_progress
        self.site_col = site_col
        self.data_col = data_col
        self.date_col = date_col

def _parse_args(args: list[str]) -> Options:
    """
    Parse command line arguments.
    """
    parser = ArgumentParser(description = "Process fluxes from netcdf files.")
    parser.add_argument("files", nargs = "+", help = "Input netcdf files")
    parser.add_argument("--out-file", "-o", required = True, help = "Output file")
    parser.add_argument("--variable", "-v", required = True, help = "Variable to process")
    parser.add_argument("--out-units", "-u", required = True, help = "Output units")
    parser.add_argument("--annual", "-a", action = "store_true", help = "Process annual data")
    parser.add_argument("--site-col", default = "site", help = "Site column name")
    parser.add_argument("--data-col", default = None, help = "Data column name")
    parser.add_argument("--date-col", default = "date", help = "Date column name")
    parser.add_argument("--log-level", "-l", type = int, default = LogLevel.INFORMATION, help = "Log level")
    parser.add_argument("--show-progress", "-p", action = "store_true", help = "Show progress")
    parsed = parser.parse_args(args)
    data_col = parsed.data_col if parsed.data_col is not None else parsed.variable
    return Options(parsed.files, parsed.out_file, parsed.variable,
                   parsed.out_units, parsed.annual, parsed.site_col,
                   data_col, parsed.date_col, parsed.log_level,
                   parsed.show_progress)

def mkdir(dir: str):
    """
    Create a directory (and any parent directories) if it does not exist.
    """
    if not os.path.exists(dir):
        log_diagnostic(f"Creating directory {dir}")
        os.makedirs(dir)
    else:
        log_debug(f"Directory {dir} already exists")

def _filter_qc(df: pandas.DataFrame) -> pandas.DataFrame:
    """
    Filter out rows with QC flags.
    """
    # If no QC flags are present, return as-is
    if "QCFlag" not in df.columns:
        return df

    # Vectorized: convert QCFlag to nullable integers, map to validity, drop invalid rows
    qc_int = pandas.to_numeric(df["QCFlag"], errors="coerce").astype("Int64")
    qc_valid_lookup = {k: v[0] for k, v in qc_flag_definitions.items()}
    flag_ok = qc_int.map(qc_valid_lookup).fillna(True)

    bad = ~flag_ok

    # Optional, compact logging by flag value
    if bad.any():
        counts = qc_int[bad].value_counts()
        for flag, n in counts.items():
            msg = qc_flag_definitions.get(int(flag), (False, "Unknown QC flag"))[1]
            log_debug(f"Filtering {n} rows for QC flag {flag}: {msg}")

    # Drop rows with invalid QC flags entirely
    return df.loc[~bad]

def _aggregate_temporally(df: pandas.DataFrame, annual: bool,
                          aggregator: Callable[[list[float]], float]) -> pandas.DataFrame:
    """
    Aggregate the input data to the specified output timestep.

    @param df: Input data frame.
    @param annual: Whether to aggregate annually. If false, aggregate daily.
    @param aggregator: Aggregator function to use.
    """
    # Group by either year or date conditionally.
    if annual:
        df = df.groupby(df.index.year)
    else:
        df = df.groupby(df.index.date)

    # Then aggregate values using the provided aggregator.
    df = df.agg(aggregator)
    return df

def _get_aggregator(variable: str) -> Callable[[list[float]], float]:
    """
    Get the aggregator function for the specified variable.
    """
    if variable not in _AGGREGATORS:
        raise ValueError(f"Variable has no temporal aggregator: {variable}")
    return _AGGREGATORS[variable]

def _read_data(file: str, variable: str, data_col: str) -> tuple[pandas.DataFrame, str]:
    """
    Read data from a netcdf file. Returns a tuple of (data frame, units).
    The data frame will have columns named after the variable.

    @param file: Input netcdf file.
    @param variable: Variable to read.
    @param data_col: Column name for the data.
    """
    with open_netcdf(file) as nc:
        if variable not in nc.variables:
            log_warning(f"Variable {variable} not found in file {file}; skipping")
            return None, None
        var = nc.variables[variable]
        time = get_var_from_std_name(nc, STD_TIME)

        if var.size != time.size:
            raise ValueError(f"Variable {variable} and time variable {time.name} have different sizes ({variable} is probably multi-dimensional; this script currently only supports 1D variables)")

        # Read variable data.
        data = var[:]

        if not hasattr(time, ATTR_UNITS):
            raise ValueError(f"Time variable {time.name} has no units attribute")

        calendar = "standard"
        if hasattr(time, ATTR_CALENDAR):
            calendar = time.calendar

        if not hasattr(var, ATTR_UNITS):
            raise ValueError(f"Variable {variable} has no units attribute")
        units = var.units

        # Read and convert the time variable to dates.
        dates = num2date(time[:], time.units, calendar,
                         only_use_cftime_datetimes = False,
                         only_use_python_datetimes = True)

        # Create a DataFrame with the data and dates; name the column after the variable
        df = pandas.DataFrame({data_col: data.flatten()}, index = dates)
        # Add QC flags if they exist.
        qc_var_name = f"{variable}_QCFlag"
        if qc_var_name in nc.variables:
            qc_var = nc.variables[qc_var_name]
            qc_data = qc_var[:]
            qc_data = qc_data.flatten()
            df["QCFlag"] = qc_data
            df = _filter_qc(df)
            # Remove QCFlag column
            df = df.drop(columns = ["QCFlag"])

        return df, units

def _guess_site_name(filename: str) -> str:
    """
    Guess the site name from the file name.
    """
    # CumberlandPlain_L6_20140101_20231230.nc
    basename = os.path.basename(filename)

    # CumberlandPlain
    site_name = basename.split("_")[0]

    # Normalise sitename. E.g. GWW -> GreatWesternWoodland, etc
    return normalise_site_name(site_name)

def _convert_units(df: pandas.DataFrame, column: str,
                   in_units: str, out_units: str,
                   annual: bool) -> pandas.DataFrame:
    """
    Convert units of the data frame to the specified units.

    @param df: Input data frame.
    @param column: Column to convert units of.
    @param out_units: Output units.
    """
    if in_units == out_units:
        return df

    # Find the conversion function (this will throw if not found).
    conversion = find_units_conversion_opt(in_units, out_units)

    # Compute timestep seconds per row.
    if annual:
        # Index is year (int). Compute exact seconds for each calendar year,
        # e.g., accounting for leap years by using Jan 1 to Jan 1 of next year.
        years = pandas.Index(df.index).astype(int)
        starts = pandas.to_datetime([f"{y}-01-01" for y in years])
        ends = pandas.to_datetime([f"{y+1}-01-01" for y in years])
        delta_s = pandas.Series((ends - starts) / numpy.timedelta64(1, 's'), index=df.index)
    else:
        # Daily (or other) aggregation: derive forward difference from datetime index (or dates).
        idx_series = pandas.Series(pandas.to_datetime(df.index), index=df.index)
        delta_s = (idx_series.shift(-1) - idx_series).dt.total_seconds()
        # Fill the last row's delta with the previous delta if possible, else 0 for a single-row frame.
        delta_s = delta_s.ffill().fillna(0)

    # Apply the conversion elementwise: f(value, timestep_seconds)
    df[column] = df[column].combine(delta_s, lambda v, s: conversion(v, s))
    return df

def main(opts):
    var = opts.variable
    mkdir(os.path.dirname(opts.out_file))
    dfs: list[pandas.DataFrame] = []
    sites: set[str] = set()
    i = 0
    for file in opts.in_files:
        aggregator = _get_aggregator(var)
        df, units = _read_data(file, var, opts.data_col)
        if df is None:
            continue
        df = _aggregate_temporally(df, opts.annual, aggregator)
        df = _convert_units(df, opts.data_col, units, opts.out_units, opts.annual)

        # Add site column from basename of the input file
        site = _guess_site_name(file)
        df[opts.site_col] = site

        # Reorder columns
        df = df[[opts.site_col, opts.data_col]]

        # Check if any df in the list already has this site in the site column.
        if site in sites:
            log_warning(f"Multiple input files contain site {site}; data from input file {file} will be ignored")
            continue
        sites.add(site)

        dfs.append(df)
        i += 1
        log_progress(1.0 * i / len(opts.in_files))

    if len(dfs) == 0:
        log_warning("No input data frames were created")
        return

    # Concatenate all data frames along rows to build a master data frame.
    master_df = pandas.concat(dfs, axis=0)

    # Give the index a name so it becomes a labelled column in CSV output.
    master_df.index.name = opts.date_col

    log_information(f"Built master dataframe with shape {master_df.shape} from {len(dfs)} files")
    master_df.to_csv(opts.out_file)
    log_information(f"Wrote master dataframe to {opts.out_file}")

if __name__ == "__main__":
    try:
        options = _parse_args(argv[1:])
        set_log_level(options.log_level)
        set_show_progress(options.show_progress)
        main(options)
    except Exception as e:
        log_error("Error: " + str(e))
        print_exc()
        exit(1)
