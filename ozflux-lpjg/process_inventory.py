#!/usr/bin/env python3
#
# This script processes ozflux flux site forest inventory data.
#
## Biomass
#
# Biomass is calculated as the sum of all biomass for a given start-end visit
# date pair, divided by plot area to get biomass in kg/m2.
#
## Height, Diameter
#
# Height and diamter are outputted as mean and 90th-percentile values for all
# live trees in the plot, measured in a single "visit".
#
## Basal Area
#
# Basal area is outputted in total m2/ha of basal area.
#

import argparse, biom_processing, numpy, os, pandas
from ozflux_sites import *
from ozflux_logging import *
from sys import argv, exit

# Column containing the start date of the visit.
COL_START_DATE = "startVisitDate"

# Column containing the end date of the visit.
COL_END_DATE = "endVisitDate"

# Column containing the date the reading was taken.
COL_DATE = "phenomenonTime"

# Name of the generated column containing plot area (in m2).
COL_PLOT_AREA = "plot_area"

# Name of the generated column containing live biomass (in kg/m2).
COL_LIVE_BIOMASS = "live_biomass"

# Name of the generated column containing dead biomass (in kg/m2).
COL_DEAD_BIOMASS = "dead_biomass"

# Name of the generated column containing mean live tree height (in m).
COL_HEIGHT = "height"

# Name of the generated column containing mean live tree diameter (in m).
COL_DIAMETER = "diameter"

# Name of the generated column containing mean live tree basal area (in m2).
COL_BA = "ba"

# Suffix appended to other column names to represent 90th-percentile values.
SFX_90 = "_90p"

# Name of the generated column containing 90th-percentile live tree height (in
# m).
COL_HEIGHT_90 = f"{COL_HEIGHT}{SFX_90}"

# Name of the generated column containing 90th-percentile live tree diameter (in
# m).
COL_DIAMETER_90 = f"{COL_DIAMETER}{SFX_90}"

# Name of the generated column containing 90th-percentile live tree basal area
# (in m2).
COL_BA_90 = f"{COL_BA}{SFX_90}"

# Possible column names for live biomass in the input data.
COLS_LIVE_BIOMASS = [
    "aboveGroundBiomass_kilograms",
    "aboveGroundLiveBiomass_kilograms",
    "stemBiomass_kilograms"
]

COLS_DEAD_BIOMASS = [
    "standingDeadAboveGroundBiomass_kilograms",
    "aboveGroundDeadBiomass_kilograms"
]

# Dictionary mapping input file names to canonical site names.
_INVENTORY_SITE_NAMES = {
    "Alice_Mulga_diameter_height_biomass_data_lVZI0qg.csv": "AliceSpringsMulga",
    "Boyagin_Wandoo_woodlands_diameter_height_biomass_da_iWYgPOL.csv": "Boyagin",
    "Calperum_Mallee_diameter_height_biomass_data_2be0UM3.csv": "Calperum",
    "Cumberland_Plain_diameter_height_biomass_data.csv": "CumberlandPlain",
    "Gingin_Banksia_stem_diameter_height_biomass_basal_area_data.csv": "Gingin",
    "GWW_diameter_height_biomass_basal_area_sapwood_data.csv": "GreatWesternWoodlands",
    "Litchfield_Savanna_diameter_height_biomass_basal_area_data.csv": "Litchfield",
    "Litchfield_Savanna_stem_diameter_height_biomass_bas_9SCuugH.csv": "Litchfield",
    "Robson_Creek_diameter_height_biomass_data_CuQoBBH.csv": "RobsonCreek",
    "Samford_diameter_height_biomass_data_mVFhel9.csv": "Samford",
    "Tumbarumba_Wet_Eucalypt_diameter_height_biomass_dat_IbRz0nd.csv": "Tumbarumba",
    "Warra_Tall_Eucalypt_diameter_height_biomass_data_teCJffC.csv": "Warra",
    "Whroo_Dry_Eucalypt_diameter_height_biomass_data_15iHyyj.csv": "Whroo",
    "Wombat_Stringybark_Eucalypt_diameter_height_biomass_EYriIKJ.csv": "WombatStateForest",
}

class Options:
    """
    Options for processing inventory data.

    @param input_dir Input directory
    @param output_dir Output directory
    @param live_file File name for live biomass (default: live_biomass.csv.gz)
    @param dead_file File name for dead biomass (default: dead_biomass.csv.gz)
    @param height_file File name for height (default: height.csv.gz)
    @param diameter_file File name for diameter (default: diameter.csv.gz)
    @param ba_file File name for basal area (default: ba.csv.gz)
    @param date_col Date column name in the output file
    @param site_col Site column name in the output file
    @param live_col Live biomass column name in the output file
    @param dead_col Dead biomass column name in the output file
    @param height_col Height column name in the output file
    @param diameter_col Diameter column name in the output file
    @param ba_col Basal area column name in the output file
    @param log_level Log level
    @param show_progress Show progress
    @param date_infmt Format of dates in the input files
    @param date_fmt Date format to be used in the output file
    """
    def __init__(self, input_dir: str, output_dir: str, live_file: str,
                 dead_file: str, height_file: str, diameter_file: str,
                 ba_file: str, date_col: str, site_col: str, live_col: str,
                 dead_col: str, height_col: str, diameter_col: str, ba_col: str,
                 log_level: LogLevel, show_progress: bool, date_infmt: str, date_fmt: str):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.live_file = live_file
        self.dead_file = dead_file
        self.height_file = height_file
        self.diameter_file = diameter_file
        self.ba_file = ba_file
        self.date_col = date_col
        self.site_col = site_col
        self.live_col = live_col
        self.dead_col = dead_col
        self.height_col = height_col
        self.diameter_col = diameter_col
        self.ba_col = ba_col
        self.verbosity = log_level
        self.show_progress = show_progress
        self.date_infmt = date_infmt
        self.date_fmt = date_fmt

def parse_args(opts: list[str]) -> Options:
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description='Process ozflux flux site forest inventory data.')
    parser.add_argument("-i", "--input-dir", required=True, help='Input directory')
    parser.add_argument("-o", "--output-dir", required=True, help='Output directory')
    parser.add_argument("--date-col", default="date", help="Date column name in the output file (default: date)")
    parser.add_argument("--site-col", default="site", help="Site column name in the output file (default: site)")
    parser.add_argument("--live-col", default="live_biomass", help="Live biomass column name in the output file (default: live_biomass)")
    parser.add_argument("--dead-col", default="dead_biomass", help="Dead biomass column name in the output file (default: dead_biomass)")
    parser.add_argument("--height-col", default="height", help="Height column name in the output file (default: height)")
    parser.add_argument("--diameter-col", default="diameter", help="Diameter column name in the output file (default: diameter)")
    parser.add_argument("--ba-col", default="ba", help="Basal area column name in the output file (default: ba)")
    parser.add_argument("--live-file", default = "live_biomass.csv.gz", help="File name for live biomass (default: live_biomass.csv.gz)")
    parser.add_argument("--dead-file", default = "dead_biomass.csv.gz", help="File name for dead biomass (default: dead_biomass.csv.gz)")
    parser.add_argument("--height-file", default = "height.csv.gz", help="File name for height (default: height.csv.gz)")
    parser.add_argument("--diameter-file", default = "diameter.csv.gz", help="File name for diameter (default: diameter.csv.gz)")
    parser.add_argument("--ba-file", default = "ba.csv.gz", help="File name for basal area (default: ba.csv.gz)")
    parser.add_argument("--date-infmt", default = "%d/%m/%Y", help = "Format of dates in the input files (default: '%d/%m/%Y')")
    parser.add_argument("--date-fmt", default = "%Y-%m-%d", help = "Date format to be used in the output file (default: '%Y-%m-%d')")
    parser.add_argument("-v", "--verbosity", default=LogLevel.INFORMATION, type = int, help=f"Log level (0-5, default {LogLevel.INFORMATION})")
    parser.add_argument("-p", "--show-progress", action="store_true", help='Show progress')
    p = parser.parse_args(opts)

    return Options(p.input_dir, p.output_dir, p.live_file, p.dead_file,
                   p.height_file, p.diameter_file, p.ba_file, p.date_col,
                   p.site_col, p.live_col, p.dead_col, p.height_col,
                   p.diameter_col, p.ba_col, p.verbosity, p.show_progress,
                   p.date_infmt, p.date_fmt)

def get_inventory_data_file(in_path: str, site: str) -> str | None:
    """
    Return the path to the inventory data file for the specified site.

    @param inventory_path: Path to a directory containing site-level csv files
                           containing timeseries of inventory data.
    @param site: Site name.
    """
    in_files = [f for f in os.listdir(in_path) if f.endswith(".csv")]
    in_files = [os.path.join(in_path, f) for f in in_files]
    for file in in_files:
        site_name = _INVENTORY_SITE_NAMES[os.path.basename(file)]
        if site_name == site:
            return file
    return None

def get_live_biomass_col(dt: pandas.DataFrame) -> str | None:
    for col in COLS_LIVE_BIOMASS:
        if col in dt.columns:
            return col
    return biom_processing.COL_LIVE

def get_dead_biomass_col(dt: pandas.DataFrame) -> str | None:
    for col in COLS_DEAD_BIOMASS:
        if col in dt.columns:
            return col
    return biom_processing.COL_DEAD

def get_inventory_data(inventory_path: str, site: str, opts: Options) -> pandas.DataFrame:
    """
    Read inventory data from the specified path.

    Returns a tuple containing the height, diameter, live biomass, and dead biomass observations.

    Returns a tuple containing the height, diameter, live biomass, and dead biomass observations.
    Note that any or all of these observations may be empty, depending on the
    input data.

    The outputs are gridcell-level values.

    @param inventory_path: Path to a directory containing site-level csv files containing timeseries of inventory data.
    @param site: Site name.
    """
    file = get_inventory_data_file(inventory_path, site)
    if file is None:
        log_warning(f"No inventory data file name provided for site '{site}'")
        return pandas.DataFrame()

    # Read data from disk.
    df = biom_processing.read_raw_data(file)

    # Add plot area column.
    if biom_processing.COL_PLOT_LENGTH in df and biom_processing.COL_PLOT_WIDTH in df:
        df[COL_PLOT_AREA] = df.apply(lambda row: biom_processing.get_plot_area(df, row.name),
                                     axis=1)
    else:
        missing = []
        if biom_processing.COL_PLOT_LENGTH not in df:
            missing.append(biom_processing.COL_PLOT_LENGTH)
        if biom_processing.COL_PLOT_WIDTH not in df:
            missing.append(biom_processing.COL_PLOT_WIDTH)
        log_warning(f"{site} inventory data is missing columns: [{missing}]")
        df[COL_PLOT_AREA] = numpy.nan

    # Convert start and end visit dates to datetime format.
    if COL_START_DATE in df and not isinstance(df[COL_START_DATE].iloc[0], datetime.datetime):
        df[COL_START_DATE] = pandas.to_datetime(df[COL_START_DATE], format=opts.date_infmt)
    if COL_END_DATE in df and not isinstance(df[COL_END_DATE].iloc[0], datetime.datetime):
        df[COL_END_DATE] = pandas.to_datetime(df[COL_END_DATE], format=opts.date_infmt)

    in_col_diam = biom_processing.get_diameter_col(df)
    if in_col_diam == None:
        log_warning(f"Inventory data file does not contain diameter column")
        return pandas.DataFrame()

    in_col_plot = biom_processing.COL_PATCH
    in_col_biomass_live = get_live_biomass_col(df)
    in_col_biomass_dead = get_dead_biomass_col(df)
    in_col_height = biom_processing.COL_HEIGHT

    # If any data columns are missing from the input file, create them and fill
    # with NaN.
    data_cols = [in_col_biomass_live, in_col_biomass_dead, in_col_height, in_col_diam]
    for col in data_cols:
        if col not in df.columns:
            log_warning(f"{site} inventory data file does not contain column '{col}'")
            df[col] = numpy.nan

    col_height_90 = f"{opts.height_col}{SFX_90}"
    col_diameter_90 = f"{opts.diameter_col}{SFX_90}"

    # Calculate basal area for each row in m2, converting diameter from cm to m.
    df[opts.ba_col] = numpy.pi * (df[in_col_diam] / 200) ** 2

    # Group by the combination of plot ID, start and end date.
    agg_map = {
        COL_PLOT_AREA: (COL_PLOT_AREA, "first"),
        COL_START_DATE: (COL_START_DATE, "first"),
        COL_END_DATE: (COL_END_DATE, "first"),
        # live/dead: sum then divide by plot area to get biomass in kg/m2.
        opts.live_col: (in_col_biomass_live, "sum"),
        opts.dead_col: (in_col_biomass_dead, "sum"),
        # Calculate mean height, diameter, and basal area.
        opts.height_col: (in_col_height, "mean"),
        opts.diameter_col: (in_col_diam, "mean"),
        # Calculate 90th percentile height and diameter, ignoring NaN values.
        col_diameter_90: (in_col_diam, lambda s: numpy.nanpercentile(s, 90)),
        col_height_90: (in_col_height, lambda s: numpy.nanpercentile(s, 90)),
        # Sum basal area.
        opts.ba_col: (opts.ba_col, "sum"),
    }

    group_cols = [in_col_plot, COL_START_DATE, COL_END_DATE]
    group_cols = [c for c in group_cols if c in df.columns and not pandas.isnull(df[c]).any()]
    if len(group_cols) == 0:
        log_warning(f"No valid columns found to group by for {site}. Verify that the site has valid data in the following columns: [{in_col_plot}, {COL_START_DATE}, {COL_END_DATE}]")
        return pandas.DataFrame()
    df = df.groupby(group_cols, as_index=False).agg(**agg_map)

    # Date column is mean of start and end dates.
    # date = start + (end - start) / 2
    if pandas.isnull(df[COL_END_DATE]).any():
        # date = start
        df[opts.date_col] = df[COL_START_DATE]
    elif pandas.isnull(df[COL_START_DATE]).any():
        # date = end
        df[opts.date_col] = df[COL_END_DATE]
    else:
        df[opts.date_col] = df.apply(lambda row: row[COL_START_DATE] + (row[COL_END_DATE] - row[COL_START_DATE]) / 2, axis=1)

    # Convert biomass from kg to kg/m2
    df[opts.live_col] = df[opts.live_col] / df[COL_PLOT_AREA]
    df[opts.dead_col] = df[opts.dead_col] / df[COL_PLOT_AREA]

    # Convert diameters from cm to m.
    df[opts.diameter_col] = df[opts.diameter_col] / 100
    df[col_diameter_90] = df[col_diameter_90] / 100

    # Convert basal area from m2/plot to m2/ha by multiplying by 10000/plotarea.
    df[opts.ba_col] = df[opts.ba_col] * 10_000 / df[COL_PLOT_AREA]

    # Replace 0-values for biomass with NaN.
    df[opts.live_col] = df[opts.live_col].replace(0, numpy.nan)
    df[opts.dead_col] = df[opts.dead_col].replace(0, numpy.nan)

    # If date column is string-typed, we need to convert it to datetime.
    if not df.empty and not isinstance(df[opts.date_col].iloc[0], datetime.datetime):
        df[opts.date_col] = pandas.to_datetime(df[opts.date_col],
                                               format = opts.date_infmt)

    return df

def write(df: pandas.DataFrame, col: str, opts: Options, filename: str):
    """
    Write the per-variable files.
    """
    # Get the columns to be written.
    data_cols = [col]
    col_p90 = f"{col}{SFX_90}"
    if col_p90 in df.columns:
        data_cols.append(col_p90)

    # Filter out rows with NaN in all data columns.
    df_filtered = df.dropna(subset=data_cols, how = "all")

    cols = [opts.date_col, opts.site_col] + data_cols

    # Generate the file path.
    file_path = os.path.join(opts.output_dir, filename)

    log_diagnostic(f"Writing {len(df_filtered)} rows of [{cols}] to {file_path}")
    df_filtered[cols].to_csv(file_path, index=False, date_format=opts.date_fmt)

def main(opts: Options):
    """
    Main entry point function.
    """
    log_diagnostic(f"Processing inventory data from {opts.input_dir} to {opts.output_dir}")
    df_global = pandas.DataFrame()
    for site in _INVENTORY_SITE_NAMES.values():
        log_information(f"Reading {site} inventory data...")
        df = get_inventory_data(opts.input_dir, site, opts)
        if df.empty:
            log_warning(f"No inventory data for site '{site}'")
            continue

        # Add site column to data frame.
        df[opts.site_col] = site
        df_global = pandas.concat([df_global, df], ignore_index=True)

    # Order data frames by site name, then date.
    df_global = df_global.sort_values([opts.site_col, opts.date_col])

    # Write the per-variable files.
    write(df_global, opts.live_col, opts, opts.live_file)
    write(df_global, opts.dead_col, opts, opts.dead_file)
    write(df_global, opts.height_col, opts, opts.height_file)
    write(df_global, opts.diameter_col, opts, opts.diameter_file)
    write(df_global, opts.ba_col, opts, opts.ba_file)

if __name__ == "__main__":
    opts = parse_args(argv[1:])
    set_log_level(opts.verbosity)
    set_show_progress(opts.show_progress)
    main(opts)
