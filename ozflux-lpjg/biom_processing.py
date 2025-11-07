#!/usr/bin/env python3
#
# This script reads site biomass observation files (.csv), aggregates the
# readings across individuals and outputs the observations in convenient units
# and structure for further use by the other benchmarking tools.
#
# The output units are kg/m2.
#

import datetime, math, ozflux_common, pandas, re, time, traceback
import numpy
import chardet
import warnings
from functools import lru_cache
from argparse import ArgumentParser
from enum import Enum
from ozflux_logging import *
from ozflux_netcdf import *
from sys import argv
from typing import Callable

# Number of KiB of data to read from each file to detect its encoding. Note that
# encoding detection only occurs if the file fails to be parsed as UTF-8.
ENCODING_DETECTION_CACHE_KB = 128

# Date format in input data.
_DATE_FMT = "%d/%m/%Y"

# Proportion of time spent reading data.
_read_prop = 1 / 3

# Proportion of time spent procesing data.
_proc_prop = 1 / 3

# Proportion of time spent writing data.
_write_prop = 1 - _read_prop - _proc_prop

################################################################################
# Column names
################################################################################

# Name of the transect ID column.
COL_TRANSECT_ID = "transectId"

# Name of the subplot ID column.
COL_SUBPLOT_ID = "subplotId"

# Name of the subplot dimension column.
COL_SUBPLOT_DIMENSION = "subplotDimension"

# Name of the plot length column.
COL_PLOT_LENGTH = "plotLength_metres"

# Name of the plot width column.
COL_PLOT_WIDTH = "plotWidth_metres"

# Name of the date column.
COL_DATE = "startVisitDate"

# Name of the transect width column.
COL_TRWIDTH = "transectWidth_metres"

# Name of the transect length column.
COL_TRLENGTH = "transectLength_metres"

# Name of the live biomass column.
COL_LIVE = "aboveGroundBiomass_kilograms"

# Name of the dead biomass column.
COL_DEAD = "standingDeadAboveGroundBiomass_kilograms"

# Name of the column which gives the individual ID.
_COL_INDIV = "plantId"

# Name of the column which gives the patch name/ID.
COL_PATCH = "plotId"

# Candidate names of the column which gives the diameter.
COLS_DIAMETER = ["stemDiameter_centimetres", "stemDiameter_centimeters"]

# Name of the column which gives the height.
COL_HEIGHT = "stemHeight_metres"

# Conversion factor from cm to m.
CM_TO_M = 0.01

class Options:
	"""
	Class for storing CLI arguments from the user.

	@param log: Log level.
	@param files: Input files.
	@param odir: Output directory.
	@param prog: True to write progress messages, 0 otherwise.
	"""
	def __init__(self, log : LogLevel, files: list[str], odir: str, prog: bool):
		self.log_level = log
		self.files = files
		self.out_dir = odir
		self.report_progress = prog

class BiomassType(Enum):
	"""
	A biomass type - either alive or dead.
	"""
	Alive = 0,
	Dead = 1

class BiomassReading():
	"""
	Represents a biomass reading on a particular day.
	"""
	def __init__(self, date: datetime.date, live: float, dead: float):
		self.date = date
		self.live = live
		self.dead = dead

class InventoryReading():
	"""
	Represents a biomass reading on a particular day.
	"""
	def __init__(self, date: datetime.date, plot: int, indiv: int, diameter: float, height: float, live_biomass: float, dead_biomass: float):
		"""
		@param date: Date of the reading.
		@param plot: Plot number.
		@param indiv: Individual number.
		@param diameter: Diameter of the tree (m).
		@param height: Height of the tree (m).
		@param live_biomass: Live biomass of the tree (kg).
		@param dead_biomass: Dead biomass of the tree (kg).
		"""
		self.date = date
		self.plot = plot
		self.indiv = indiv
		self.diameter = diameter
		self.height = height
		self.live_biomass = live_biomass
		self.dead_biomass = dead_biomass
	def __str__(self) -> str:
		return f"Patch {self.plot}, Individual {self.indiv}: {self.date} - {self.live_biomass}kg (live), {self.dead_biomass}kg (dead). Height: {self.height}m, Diameter: {self.diameter}m"

class Plot:
	"""
	Represents a plot of land.
	"""
	def __init__(self, id: int, name: str, area: float):
		self.id = id
		self.name = name
		self.area = area # can be NaN if data doesn't contain plot length/width.
	def __str__(self) -> str:
		return f"Plot {self.id}: {self.name} ({self.area}m2)"

class Inventory:
	"""
	Represents the inventory of a gridcell.
	"""
	def __init__(self, plots: list[Plot], readings: list[InventoryReading]):
		self.plots = plots
		self.readings = readings
	def get_area(self, plot: int) -> float:
		"""
		Get the area of a plot.

		@param plot: plot number.
		@return Area of the plot (m2).
		"""
		plots = [p for p in self.plots if p.id == plot]
		if len(plots) != 1:
			raise ValueError(f"Patch {plot} not found")
		return plots[0].area
	def __str__(self) -> str:
		return "\n".join([f"{str(p)}: {len([r for r in self.readings if r.plot == p.id])} readings over {len(numpy.unique([r.date for r in self.readings if r.plot == p.id]))} timesteps" for p in self.plots])

def parse_args(argv: list[str]) -> Options:
	"""
	Parse CLI arguments, return a parsed options object.

	@param argv: Raw CLI arguments.

	@return Parsed Options object.
	"""
	parser = ArgumentParser(prog=argv[0], description = "Formatting ozflux data into a format suitable for consumption by LPJ-Guess")
	parser.add_argument("-v", "--verbosity", type = int, help = "Logging verbosity (1-5, default 3)", nargs = "?", const = LogLevel.INFORMATION, default = LogLevel.INFORMATION)
	parser.add_argument("files", nargs = "+", help = "Input .nc files to be processed")
	parser.add_argument("-o", "--out-dir", required = True, help = "Path to the output directory. Processed files will be saved with the same file name into this directory")
	parser.add_argument("-p", "--show-progress", action = "store_true", help = "Report progress")
	parser.add_argument("--version", action = "version", version = "%(prog)s " + ozflux_common.VERSION)

	p = parser.parse_args(argv[1:])
	return Options(p.verbosity, p.files, p.out_dir, p.show_progress)

def area_from_dim(dim: str) -> float:
	"""
	Read a string in dimension (e.g. 20m x 10m) format and return it as an area
	in m2.

	@param dim: The input dimension string.
	"""
	#                 20                  m                x      20                   m
	pattern = r'[ \t]*(\d+(?:\.\d+)?)[ \t]([A-Za-z]+)[ \t]*x[ \t]*(\d+(?:\.\d+)?)[ \t]*([A-Za-z]+).*'
	match = re.match(pattern, dim)
	if match == None:
		raise ValueError("Cannot parse area from '%s'" % dim)

	# Trim whitespace.
	groups = match.groups()
	length = float(groups[0].strip())
	length_units = groups[1].strip()
	width = float(groups[2].strip())
	width_units = groups[3].strip()

	if math.isnan(length):
		raise ValueError("Cannot parse length from '%s'" % dim)
	if math.isnan(width):
		raise ValueError("Cannot parse width from '%s'" % dim)

	# Units conversion.
	length = to_metres(length, length_units)
	width = to_metres(width, width_units)

	return length * width

def transect_area_from_length_width(data: pandas.DataFrame, row: int) -> float:
	"""
	Get the transect area for a row from the transect length/width column.

	@param data: The data frame.
	@param row: Row index.
	"""
	transect_width = ozflux_common.get_length(data.iloc[row][COL_TRWIDTH])
	transect_length = ozflux_common.get_length(data.iloc[row][COL_TRLENGTH])
	return transect_length * transect_width

def transect_area_from_transect_dim(data: pandas.DataFrame, row: int) -> float:
	"""
	Get the transect area for a row from the transect dimension column.

	@param data: The data frame.
	@param row: Row index.
	"""
	dim = data.iloc[row][COL_SUBPLOT_DIMENSION]
	return area_from_dim(dim)

def get_plot_length(data: pandas.DataFrame, row: int, col: str) -> float:
	"""
	Get the plot length specified on the given row in m.

	@param data: The data frame.
	@param row: The row index.
	"""
	log_debug(f"Attempting to parse {col} {row} from '{data.iloc[row][col]}'")
	raw = data.iloc[row][col]
	if type(raw) == float and math.isnan(raw):
		log_debug(f"Missing {col} in row {row}")
		# Fallback: group data frame by plot ID, get all distinct non-NaN values
		# and if there is one value, use that. If there are multiple, we will
		# need to emit an exception.
		patch = data.iloc[row][COL_PATCH]
		group = data[data[COL_PATCH] == patch]
		lengths = group[COL_PLOT_LENGTH].dropna().unique()
		if len(lengths) != 1:
			raise ValueError(f"Multiple plot lengths found for plot {patch}")
		else:
			return ozflux_common.get_length(lengths[0])
	else:
		return ozflux_common.get_length(raw)

def get_plot_area(data: pandas.DataFrame, row: int) -> float:
	"""
	Get the plot area specified on the given row in m2.

	@param data: The data frame.
	@param row: The row index.
	"""
	if not COL_PLOT_LENGTH in data or not COL_PLOT_WIDTH in data:
		log_warning(f"Inventory data file does not contain plot length/width data")
		return math.nan
	# This will throw if these column don't exist.
	plot_length = get_plot_length(data, row, COL_PLOT_LENGTH)
	plot_width = get_plot_length(data, row, COL_PLOT_WIDTH)
	return plot_length * plot_width

def get_transect_area(data: pandas.DataFrame, row: int) -> float:
	"""
	Get the transect area of the specified row, in m2.

	@param data: The data frame.
	@param row: Row index.
	"""
	# Get transect area in m2.
	transect_area = 0
	if COL_TRWIDTH in data:
		return transect_area_from_length_width(data, row)

	if COL_SUBPLOT_DIMENSION in data:
		return transect_area_from_transect_dim(data, row)

	# Assume measurements are per-plot, rather than per-transect.
	return get_plot_area(data, row)

def parse_float(row: pandas.Series, col_name: str) -> float:
	"""
	Get a specific value for a given row.
	"""
	value = float(row[col_name])

	if math.isnan(value):
		log_debug(f"NaN {col_name} recorded on row {row.name}.")
		return math.nan

	return value

def get_biomass(data: pandas.DataFrame, row: int, biom_type: BiomassType
	, area: float) -> float:
	"""
	Get biomass reading for a given row in kg/ha.

	@param data: Data frame.
	@param row: Index of the row to read.
	@param biom_type: Biomass type to read (alive or dead).
	@param area: Transect area on this row.
	"""
	col_name = COL_LIVE if biom_type == BiomassType.Alive else COL_DEAD
	return parse_float(data, row, col_name) / area

def get_transectid_colname(data: pandas.DataFrame) -> str:
	"""
	Get the name of the transect ID column in the data frame.

	@param data: The data frame.
	"""
	if COL_TRANSECT_ID in data:
		return COL_TRANSECT_ID
	if COL_SUBPLOT_ID in data:
		return COL_SUBPLOT_ID
	# No transect ID column.
	return None

def get_num_transects(data: pandas.DataFrame, rows: list[int]) -> int:
	"""
	Get the number of different transects across the specified rows.

	@param data: The data frame.
	@param rows: Row indices of the rows.
	"""
	col_name = get_transectid_colname(data)
	if col_name == None:
		# Assume no transect - ie measurements across entire plot.
		return 1
	return data.iloc[rows].groupby(col_name).ngroups

def get_diameter_col(df: pandas.DataFrame) -> str:
	"""
	Get the name of the diameter column in the data frame.

	@param df: The data frame.
	"""
	for col in COLS_DIAMETER:
		if col in df.columns:
			return col
	return None

def parse_inventory_reading(row: pandas.Series, patches: list[str], indivs: list[str]) -> InventoryReading:
	"""
	Parse an inventory reading from a row in a data frame.

	@param row: The row to parse.
	@param patches: List of unique patch names.
	@param indivs: List of unique individual IDs (note these may not be ints).
	"""
	date = datetime.datetime.strptime(row[COL_DATE], _DATE_FMT)
	patch = int(numpy.where(patches == row[COL_PATCH])[0][0])
	indiv = int(numpy.where(indivs == row[_COL_INDIV])[0][0])

	diameter = math.nan
	height = math.nan
	live_biomass = math.nan
	dead_biomass = math.nan

	diameter_cols = [c for c in COLS_DIAMETER if c in row.keys()]
	col_diameter = diameter_cols[0] if len(diameter_cols) > 0 else None
	if col_diameter != None:
		diameter = parse_float(row, col_diameter) * CM_TO_M

	if COL_HEIGHT in row.keys():
		height = parse_float(row, COL_HEIGHT) # already in m
	if COL_LIVE in row.keys():
		live_biomass = parse_float(row, COL_LIVE) # kg
	if COL_DEAD in row.keys():
		dead_biomass = parse_float(row, COL_DEAD) # kg

	return InventoryReading(date, patch, indiv, diameter, height, live_biomass,
							dead_biomass)

def parse_patch(data: pandas.DataFrame, patch_name: str, patches: list[str]) -> Plot:
	"""
	Parse a patch from a data frame.

	@param data: The data frame.
	@param patch_name: The name of the patch.
	"""
	patch_rows = data[data[COL_PATCH] == patch_name]
	if len(patch_rows) == 0:
		raise ValueError(f"Patch '{patch_name}' has no data.")
	area = get_plot_area(patch_rows, 0)
	id = int(numpy.where(patches == patch_name)[0][0])
	return Plot(id, patch_name, area)

def read_inventory_data(data: pandas.DataFrame) -> Inventory:
	"""
	Get a list of inventory readings from a steam/diameter/biomass indiv-level
	data file.

	@param data: Data frame containing the data.
	"""
	patch_names = data[COL_PATCH].unique()
	indivs = data[_COL_INDIV].unique()
	readings = data.apply(parse_inventory_reading, axis=1, args=(patch_names, indivs))
	patches = [parse_patch(data, patch_name, patch_names) for patch_name in patch_names]
	return Inventory(patches, readings)

def read_biomass_data(data: pandas.DataFrame, pcb: Callable[[float], None]) \
	-> list[BiomassReading]:
	"""
	Aggregate biomass from per-individual into kg/ha.

	Returns a list of tuples of (date, live biomass, dead biomass).
	"""
	result = []

	groups = data.groupby(COL_DATE)
	n = len(groups.groups)
	i = 0

	for date, rows in groups.groups.items():
		live = 0 # kg/ha
		dead = 0 # kg/ha
		for row in rows:
			area = get_transect_area(data, row)
			live += parse_float(data.iloc[row], COL_LIVE) / area
			dead += parse_float(data.iloc[row], COL_DEAD) / area

		# Each date can contain readings from multiple transects. We need to
		# take the mean of each transect's total biomass.
		n_transects = get_num_transects(data, rows)
		live /= n_transects
		dead /= n_transects

		dt = datetime.datetime.strptime(date, _DATE_FMT)
		result.append(BiomassReading(dt, live, dead))
		i += 1
		pcb(i / n)
	return result

@lru_cache(maxsize=128)
def detect_encoding(file: str) -> str:
	"""
	Detect file encoding by reading a sample of the file.
	Results are cached to avoid re-scanning the same file.
	"""
	# Read at most 128KB to detect encoding
	sample_size = ENCODING_DETECTION_CACHE_KB * 1024
	with open(file, "rb") as f:
		raw = f.read(sample_size)
		result = chardet.detect(raw)
		if result["confidence"] < 1:
			log_warning(f"Encoding for file '{file}' has been guessed as '{result['encoding']}' with confidence {result['confidence']}.")
		return result["encoding"]

def read_biomass_file(file: str, encoding: str) -> pandas.DataFrame:
	"""
	Read raw data from the input file, using the specified encoding.

	@param file: The input file name.
	@param encoding: The encoding to use.
	"""
	# The input files are inconsistent in their representation of missing data.
	nas = ["Na", "NA", "na", "nA"]

	# Using low_memory=False to avoid data type inference issues. The
	# RobsonCreek inventory data file contains plant IDs that are mostly
	# integer, but some look like 604_1. By coincidence, none of these
	# non-numeric IDs appear in the first 16K lines. This causes the type
	# inference engine to give mixed types - int for the first chunk, then str
	# for the rest. Setting low_memory=False causes the entire file to be read
	# in one go, which will cause the type inference engine to run on the entire
	# column, which in this case will give us a string column. Without
	# low_memory=False, we get a mixed column. Actually, we don't even use this
	# column, so it doesn't matter too much, and we have small files and are
	# going to hold the data in memory anyway, so this doesn't really cost us
	# much, and it silences the type inference warning.
	#
	# Note though, that if we ever do want to use this column, we need to
	# remember that it is, in principle string typed, and therefore needs to be
	# handled similarly to the patch column.
	return pandas.read_csv(file, parse_dates = True, na_values = nas,
						   encoding=encoding,low_memory=False)

def read_raw_data(file: str) -> pandas.DataFrame:
	"""
	Read raw data from the input file.
	Attempts UTF-8 first, then tries to detect encoding, falling back to latin1 if needed.
	"""
	try:
		# Try UTF-8 first.
		return read_biomass_file(file, "utf-8")
	except UnicodeDecodeError:
		# If UTF-8 fails, detect encoding. Looking at you, GWW...
		encoding = detect_encoding(file)
		if encoding != "utf-8":
			warn = f"File {file} is not UTF-8 encoded (detected {encoding})"
			log_warning(warn)
		else:
			# Should hopefullly never happen. If it does, the below call will
			# almost certainly fail.
			log_warning(f"File {file} is allegedly UTF-8 encoded. Parsing failed, but the detected encoding is: {encoding}. Something is very wrong here.")

		return read_biomass_file(file, encoding)

def write(data: list[BiomassReading], out_file: str
	, pcb: Callable[[float], None]):
	"""
	Write data to the specified output file in .csv format.
	"""
	log_diagnostic("Writing '%s'..." % out_file)
	with open(out_file, "w") as file:
		i = 0
		file.write("date,live_biomass,dead_biomass\n")
		for (date, live, dead) in data:
			file.write("%s,%.2f,%.2f\n" % (date.strftime("%Y-%m-%d"), live, dead))
			i += 1
			pcb(i / len(data))

def process_file(file: str, out_dir: str, pcb: Callable[[float], None]):
	"""
	Read input data from file, aggregate data, write to a new file in out dir.

	@param file: Input file.
	@param out_dir: Output directory.
	@param pcb: Progress callback function.
	"""
	log_information("Processing '%s'..." % file)

	global _read_prop, _proc_prop, _write_prop

	pcb(0)
	read_start = time.time()
	data = read_raw_data(file)
	read_time = time.time() - read_start
	pcb(_read_prop)

	proc_start = time.time()
	data = read_biomass_data(data, lambda p: pcb(_read_prop + p * _proc_prop))
	proc_time = time.time() - proc_start

	# Determine output file name.
	filename = os.path.basename(file)
	filename = filename.split('_')[0] + '.csv'
	out_file = os.path.join(out_dir, os.path.basename(file))

	start = _read_prop + _proc_prop
	write_start = time.time()
	write(data, out_file, lambda p: pcb(start + _write_prop * p))
	write_time = time.time() - write_start

	total_time = read_time + proc_time + write_time
	_read_prop = read_time / total_time
	_proc_prop = proc_time / total_time
	_write_prop = write_time / total_time
	log_diagnostic("Time distribution: read=%.2f%%, process=%.2f%%, write=%.2f%%" \
		% (_read_prop, _proc_prop, _write_prop))
	log_debug("Successfully processed file '%s'" % file)

def main(opts: Options):
	"""
	Main function.

	@param opts: Parsed CLI options provided by the user..
	"""
	step = 1 / len(opts.files)
	start = 0
	for file in opts.files:
		process_file(file, opts.out_dir, lambda p: log_progress(start + step*p))
		start += step

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
