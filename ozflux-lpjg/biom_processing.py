#!/usr/bin/env python3
#
# This script reads site biomass observation files (.csv), aggregates the
# readings across individuals and outputs the observations in convenient units
# and structure for further use by the other benchmarking tools.
#
# The output units are kg/m2.
#

import datetime, math, ozflux_common, pandas, re, time, traceback

from argparse import ArgumentParser
from enum import Enum
from ozflux_logging import *
from ozflux_netcdf import *
from sys import argv
from typing import Callable

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

def get_plot_area(data: pandas.DataFrame, row: int) -> float:
	"""
	Get the plot area specified on the given row.

	@param data: The data frame.
	@param row: The row index.
	"""
	# This will throw if these column don't exist.
	plot_length = ozflux_common.get_length(data.iloc[row][COL_PLOT_LENGTH])
	plot_width = ozflux_common.get_length(data.iloc[row][COL_PLOT_WIDTH])
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

	biomass = float(data.iloc[row][col_name])

	if math.isnan(biomass):
		log_debug("NaN %s recorded on row %d. Using 0 instead..." % \
			(col_name, row))
		return 0

	return biomass / area

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

def read_biomass_data(data: pandas.DataFrame, pcb: Callable[[float], None]) \
	-> list[BiomassReading]:
	"""
	Aggregate biomass from per-individual into kg/ha.

	Returns a list of tuples of (date, live biomass, dead biomass).
	"""
	# Date format in input data.
	date_fmt = "%d/%m/%y"

	result = []

	groups = data.groupby(COL_DATE)
	n = len(groups.groups)
	i = 0

	for date, rows in groups.groups.items():
		live = 0 # kg/ha
		dead = 0 # kg/ha
		for row in rows:
			area = get_transect_area(data, row)
			live += get_biomass(data, row, BiomassType.Alive, area)
			dead += get_biomass(data, row, BiomassType.Dead, area)

		# Each date can contain readings from multiple transects. We need to
		# take the mean of each transect's total biomass.
		n_transects = get_num_transects(data, rows)
		live /= n_transects
		dead /= n_transects

		dt = datetime.datetime.strptime(date, date_fmt)
		result.append(BiomassReading(dt, live, dead))
		i += 1
		pcb(i / n)
	return result

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

def read_raw_data(file: str) -> pandas.DataFrame:
	"""
	Read raw data from the input file.
	"""
	nas = ["Na", "NA", "na", "nA"]
	return pandas.read_csv(file, parse_dates = True, na_values = nas)

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
