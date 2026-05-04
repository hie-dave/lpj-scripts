#!/usr/bin/env python3
#
# Process NATT stand inventory data into per-variable CSV files.
#

import argparse
import os
from dataclasses import dataclass

import numpy
import pandas

from ozflux_logging import *
from sys import argv

# Input column names.
COL_IN_DATE = "obs_time"
COL_IN_SITE = "site"

# Output key columns.
COL_OUT_DATE = "date"
COL_OUT_SITE = "site"

# Biomass columns in input (kg/ha).
COL_AGB = "agb_drymass_ha"
COL_BGB = "bgb_drymass_ha"
COL_TB = "tb_drymass_ha"

# Optional dead biomass columns in input (kg/ha).
COLS_DEAD_BIOMASS = [
    "dead_biomass_ha",
    "dead_drymass_ha",
    "dead_agb_drymass_ha",
]

# Basal area columns in input (m2/ha).
COL_BA_LIVE = "live_basal_area_ha"
COL_BA_DEAD = "dead_basal_area_ha"

# Height columns in input (m).
COL_H_MEAN = "ht.mean"
COL_H_MAX = "ht.max"
COL_H_QUAD = "ht.quad_mean"
COL_H_P90 = "ht.pc90"
COL_H_SD = "ht.sd"

# Diameter columns in input (cm).
COL_D_MEAN = "diameter.mean"
COL_D_MAX = "diameter.max"
COL_D_QUAD = "diameter.quad_mean"
COL_D_P90 = "diameter.pc90"
COL_D_SD = "diameter.sd"

@dataclass
class Options:
    input_file: str
    output_dir: str
    live_file: str
    dead_file: str
    height_file: str
    diameter_file: str
    ba_file: str
    date_col: str
    site_col: str
    date_in_col: str
    site_in_col: str
    date_infmt: str
    date_fmt: str
    verbosity: LogLevel
    show_progress: bool

def parse_args(args: list[str]) -> Options:
    parser = argparse.ArgumentParser(description="Process NATT stand inventory data.")
    parser.add_argument(
        "-i",
        "--input-file",
        default="data/inventory/raw/NATT_stand_data.csv",
        help="Input NATT CSV file (default: data/inventory/raw/NATT_stand_data.csv)",
    )
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("--live-file", default="live_biomass.csv.gz", help="Live biomass output file")
    parser.add_argument("--dead-file", default="dead_biomass.csv.gz", help="Dead biomass output file")
    parser.add_argument("--height-file", default="height.csv.gz", help="Height output file")
    parser.add_argument("--diameter-file", default="diameter.csv.gz", help="Diameter output file")
    parser.add_argument("--ba-file", default="ba.csv.gz", help="Basal area output file")
    parser.add_argument("--date-col", default=COL_OUT_DATE, help="Date column name in output")
    parser.add_argument("--site-col", default=COL_OUT_SITE, help="Site column name in output")
    parser.add_argument("--date-in-col", default=COL_IN_DATE, help="Date column name in input")
    parser.add_argument("--site-in-col", default=COL_IN_SITE, help="Site column name in input")
    parser.add_argument("--date-infmt", default="%Y-%m-%d", help="Input date format")
    parser.add_argument("--date-fmt", default="%Y-%m-%d", help="Output date format")
    parser.add_argument(
        "-v",
        "--verbosity",
        default=LogLevel.INFORMATION,
        type=int,
        help=f"Log level (0-5, default {LogLevel.INFORMATION})",
    )
    parser.add_argument("-p", "--show-progress", action="store_true", help="Show progress")
    p = parser.parse_args(args)

    return Options(
        input_file=p.input_file,
        output_dir=p.output_dir,
        live_file=p.live_file,
        dead_file=p.dead_file,
        height_file=p.height_file,
        diameter_file=p.diameter_file,
        ba_file=p.ba_file,
        date_col=p.date_col,
        site_col=p.site_col,
        date_in_col=p.date_in_col,
        site_in_col=p.site_in_col,
        date_infmt=p.date_infmt,
        date_fmt=p.date_fmt,
        verbosity=p.verbosity,
        show_progress=p.show_progress,
    )

def get_col(df: pandas.DataFrame, col: str) -> pandas.Series:
    if col in df.columns:
        return df[col]
    log_warning(f"Missing input column '{col}'")
    return pandas.Series(numpy.nan, index=df.index)

def get_first_available_col(df: pandas.DataFrame, cols: list[str]) -> pandas.Series:
    for col in cols:
        if col in df.columns:
            return df[col]
    log_warning(f"Missing all optional input columns: {cols}")
    return pandas.Series(numpy.nan, index=df.index)

def write_variable(
    df: pandas.DataFrame, metric_cols: list[str], date_col: str, site_col: str, path: str, date_fmt: str
):
    filtered = df.dropna(subset=metric_cols, how="all")
    cols = [date_col, site_col] + metric_cols
    log_diagnostic(f"Writing {len(filtered)} rows of [{cols}] to {path}")
    filtered[cols].to_csv(path, index=False, date_format=date_fmt)

def build_output_tables(df: pandas.DataFrame, opts: Options) -> dict[str, pandas.DataFrame]:
    out = pandas.DataFrame()
    out[opts.date_col] = pandas.to_datetime(get_col(df, opts.date_in_col), format=opts.date_infmt)
    out[opts.site_col] = get_col(df, opts.site_in_col)

    # live biomass: kg/ha -> kg/m2
    live = out.copy()
    live["above_ground"] = get_col(df, COL_AGB) / 10_000.0
    live["below_ground"] = get_col(df, COL_BGB) / 10_000.0
    live["total"] = get_col(df, COL_TB) / 10_000.0

    # dead biomass: optional input columns, kg/ha -> kg/m2
    dead = out.copy()
    dead["mean"] = get_first_available_col(df, COLS_DEAD_BIOMASS) / 10_000.0

    # Height metrics (already in m)
    height = out.copy()
    height["mean"] = get_col(df, COL_H_MEAN)
    height["quad_mean"] = get_col(df, COL_H_QUAD)
    height["p90"] = get_col(df, COL_H_P90)
    height["max"] = get_col(df, COL_H_MAX)
    height["sd"] = get_col(df, COL_H_SD)

    # Diameter metrics (cm -> m)
    diameter = out.copy()
    diameter["mean"] = get_col(df, COL_D_MEAN) / 100.0
    diameter["quad_mean"] = get_col(df, COL_D_QUAD) / 100.0
    diameter["p90"] = get_col(df, COL_D_P90) / 100.0
    diameter["max"] = get_col(df, COL_D_MAX) / 100.0
    diameter["sd"] = get_col(df, COL_D_SD) / 100.0

    # Basal area metrics (already in m2/ha)
    ba = out.copy()
    ba["live"] = get_col(df, COL_BA_LIVE)
    ba["dead"] = get_col(df, COL_BA_DEAD)
    ba["total"] = ba["live"] + ba["dead"]

    for table in [live, dead, height, diameter, ba]:
        table.sort_values([opts.site_col, opts.date_col], inplace=True)

    return {
        opts.live_file: live,
        opts.dead_file: dead,
        opts.height_file: height,
        opts.diameter_file: diameter,
        opts.ba_file: ba,
    }

def main(opts: Options):
    log_diagnostic(f"Reading NATT data from {opts.input_file}")
    df = pandas.read_csv(opts.input_file)

    outputs = build_output_tables(df, opts)

    os.makedirs(opts.output_dir, exist_ok=True)
    write_variable(
        outputs[opts.live_file],
        ["above_ground", "below_ground", "total"],
        opts.date_col,
        opts.site_col,
        os.path.join(opts.output_dir, opts.live_file),
        opts.date_fmt,
    )
    write_variable(
        outputs[opts.dead_file],
        ["mean"],
        opts.date_col,
        opts.site_col,
        os.path.join(opts.output_dir, opts.dead_file),
        opts.date_fmt,
    )
    write_variable(
        outputs[opts.height_file],
        ["mean", "quad_mean", "p90", "max", "sd"],
        opts.date_col,
        opts.site_col,
        os.path.join(opts.output_dir, opts.height_file),
        opts.date_fmt,
    )
    write_variable(
        outputs[opts.diameter_file],
        ["mean", "quad_mean", "p90", "max", "sd"],
        opts.date_col,
        opts.site_col,
        os.path.join(opts.output_dir, opts.diameter_file),
        opts.date_fmt,
    )
    write_variable(
        outputs[opts.ba_file],
        ["live", "dead", "total"],
        opts.date_col,
        opts.site_col,
        os.path.join(opts.output_dir, opts.ba_file),
        opts.date_fmt,
    )

if __name__ == "__main__":
    opts = parse_args(argv[1:])
    set_log_level(opts.verbosity)
    set_show_progress(opts.show_progress)
    main(opts)
