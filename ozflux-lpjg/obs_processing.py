#!/usr/bin/env python3
#
# This script reads the raw .nc files and produces csv files with observations
# used in benchmarking.
#

#%pip install openpyxl numpy pandas

import math
import biom_processing, datetime, glob, numpy, pandas, time, traceback
from multiprocessing import Lock
import ozflux_common
import difflib
from itertools import groupby
from operator import attrgetter
from argparse import ArgumentParser
from ozflux_logging import *
from ozflux_netcdf import *
from ozflux_observations import *
from sys import argv
from netCDF4 import Dataset, date2index, num2date, date2num
from typing import Callable
from ozflux_parallel import JobManager

_ATTR_HEIGHT = "height"

# Name of the latitude variable in output file.
_VAR_LAT = "latitude"

# Name of the longitude variabe in the output file.
_VAR_LON = "longitude"

# Default chunk size of the time chunks in output .nc file. Note that the actual
# chunk size will be less than this if the input data contains less than this
# number of time points.
_TIME_CHUNK_SIZE = 366 * 24

_VAR_DIMS   = (DIM_TIME, DIM_LON, DIM_LAT)

# Total soil profile depth in m.
_SOIL_DEPTH = 1.5

# Name of the soil water column in the output file.
_NAME_SW = "swvol"

# Name of the LAI column in the output file.
_NAME_LAI = "lai"

# Name of the 90cm soil water column in the output file.
_NAME_SMIPS = "swmm_90"

# Name of the 90cm soil moisture index in the output file.
_NAME_SWINDEX = "swindex_90"

# Name of the live biomass column in the output file.
_NAME_LIVE = "live_biomass"

# Name of the dead biomass column in the output file.
_NAME_DEAD = "dead_biomass"

# Name of the height column in the output file.
_NAME_HEIGHT = "height"

# Name of the diameter column in the output file.
_NAME_DIAMETER = "diameter"

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

# Name of the date column in smips input files.
_SMIPS_COL_DATE = "date"

# Name of the SW at 90cm column in smips input files.
_SMIPS_COL_SW = "totalbucket"

# Name of the soil moisture index column in smips input files.
_SMIPS_COL_SWINDEX = "SMindex"

# Minimum allowed LAI value (m2/m2).
_LAI_MIN = 0

# Maximum allowed LAI value (m2/m2).
_LAI_MAX = 20

# Standard variables in the ozflux files which will be used for benchmarking.
_standard_variables = [
	# GPP shouldn't get this high, but this is mainly for error checking anyway.
	ForcingVariable("GPP_LT", "gpp", "gC/m2/day", numpy.mean, -20000, 20000),
	ForcingVariable("ER_LT", "resp", "gC/m2/day", numpy.mean, -10000, 10000),
	ForcingVariable("NEE_LT", "nee", "gC/m2/day", numpy.mean, -10000, 10000),
	# todo: are these bounds sensible?
	ForcingVariable("ET", "transpiration", "mm/day", numpy.mean, -100, 100),
]

class Layer:
	def __init__(self, name: str, depth: int):
		"""
		Create a new Layer object.
		@param name: Name of the variable in the netcdf file.
		@param depth: Depth of this layer in cm.
		"""
		self.name = name
		self.depth = depth

# Map of site names to a list of names of soil water layer variables in the flux
# data file. The site name should be the .nc file name before the first '_'.
_layer_variable_map = {
	"AdelaideRiver": [Layer("Sws", 5)],
	"AliceSpringsMulga": [Layer("Sws_mulga", 0), Layer("Sws_10cm_mulga_2b", 10), Layer("Sws_60cm_mulga_2", 60), Layer("Sws_100cm_mulga_2", 100)],
	# # todo: what is the difference between Sws_05cma_N, Sws_05cma_S, Sws_05cmb_N, Sws_05cmb_S???
	"Boyagin": [Layer("Sws_05cma_N", 5), Layer("Sws_08cm_main", 8), Layer("Sws_10cma_N", 10), Layer("Sws_15cm_main", 15), Layer("Sws_20cma_N", 20), Layer("Sws_40cm_N", 40), Layer("Sws_80cm_N", 80), Layer("Sws_110cm_N", 110)],
	"Calperum": [Layer("Sws", 10), Layer("Sws_25cm", 25), Layer("Sws_50cm", 50), Layer("Sws_100cm", 100)],
	"CapeTribulation": [Layer("Sws", 10), Layer("Sws_75cma", 75)],
	"Collie": [Layer("Sws", 5), Layer("Sws_10cm", 10), Layer("Sws_30cm", 30)],
	"CowBay": [Layer("Sws", 15), Layer("Sws_75cma", 75)],
	"CumberlandPlain": [Layer("Sws", 8), Layer("Sws_20cm", 20)],
	"DalyPasture": [Layer("Sws", 5), Layer("Sws_50", 50)],
	"DalyUncleared": [Layer("Sws", 5), Layer("Sws_50cm", 50)],
	"DryRiver": [Layer("Sws", 1.5), Layer("Sws_50cm", 50)],
	# # Note: Emerald does have Sws, but it has no height attribute!!!!
	"Emerald": [Layer("Sws_5cmd", 5), Layer("Sws_15cmd", 15), Layer("Sws_22p5cmd", 22.5), Layer("Sws_30cmd", 30)],
	# Note: FletcherView has Sws, Sws_bare, Sws_grass, Sws_shrub, and Sws_tree
	"Fletcherview": [Layer("Sws", 0)],
	# # Note FoggDam does have Sws, but it has no height attribute!!!!
	"FoggDam": [Layer("Sws_10", 10)],
	"Gingin": [Layer("Sws", 6), Layer("Sws_20cm", 20), Layer("Sws_40cm", 40), Layer("Sws_80cm", 80)],
	"GreatWesternWoodlands": [Layer("Sws_5cmb", 5), Layer("Sws", 10), Layer("Sws_20cmb", 20), Layer("Sws_30cmd", 30), Layer("Sws_50cmb", 50), Layer("Sws_70cmb", 70)],
	"HowardSprings": [Layer("Sws", 10), Layer("Sws40cm", 40), Layer("Sws100cm", 100), Layer("Sws120cm", 120), Layer("Sws140cm", 140)],
	"Litchfield": [Layer("Sws", 5), Layer("Sws20cm", 20), Layer("Sws50cm", 50), Layer("Sws100cm", 100), Layer("Sws140cm", 140)],
	"Longreach": [Layer("Sws", 10)],
	"Otway": [Layer("Sws", 0)],
	# # Note RedDirt does have Sws, but it has no height attribute!!!!
	"RedDirtMelonFarm": [Layer("Sws_5cm", 5), Layer("Sws_50cm", 50)],
	"Ridgefield": [Layer("Sws", 5), Layer("Sws_10cm", 10), Layer("Sws_20cm", 20), Layer("Sws_40cm", 40), Layer("Sws_80cm", 80)],
	"RiggsCreek": [Layer("Sws", 5), Layer("Sws_50cm", 50)],
	"RobsonCreek": [Layer("Sws", 6), Layer("Sws_75cma", 75)],
	"Samford": [Layer("Sws", 10)],
	"SilverPlains": [Layer("Sws", 10)],
	"SturtPlains": [Layer("Sws", 5), Layer("Sws_50cm", 50)],
	# # Note: TiTreeEast has Sws_grass, Sws_mulga, and Sws_spinifex
	"TiTreeEast": [Layer("Sws", 0)],
	"Tumbarumba": [Layer("Sws", 0)],
	"WallabyCreek": [Layer("Sws", 10), Layer("Sws_15cma", 15), Layer("Sws_20cma", 20), Layer("Sws_40cma", 40), Layer("Sws_80cma", 80)],
	"Warra": [Layer("Sws", 8), Layer("Sws_20cmd", 20), Layer("Sws_40cmd", 40), Layer("Sws_80cmd", 80), Layer("Sws_1p2m", 120)],
	"Whroo": [Layer("Sws", 10)],
	"WombatStateForest": [Layer("Sws", 10), Layer("Sws_35cma", 35), Layer("Sws_65cma", 65), Layer("Sws_95cma", 95)],
	"Yanco": [Layer("Sws_3cm", 3), Layer("Sws", 10), Layer("Sws_15cm", 15), Layer("Sws_45cm", 45), Layer("Sws_75cm", 75)],
	"YarCon": [Layer("Sws_5cm", 5), Layer("Sws", 10)],
	"YarIrr": [Layer("Sws_5cm", 5), Layer("Sws", 10)]
}

_inventory_site_names = {
	"Alice_Mulga_diameter_height_biomass_data_lVZI0qg.csv": "AliceSpringsMulga",
	"Boyagin_Wandoo_woodlands_diameter_height_biomass_da_iWYgPOL.csv": "Boyagin",
	"Calperum_Mallee_diameter_height_biomass_data_2be0UM3.csv": "Calperum",
	"Cumberland_Plain_diameter_height_biomass_data.csv": "CumberlandPlain",
	"Gingin_Banksia_stem_diameter_height_biomass_basal_area_data.csv": "Gingin",
	"GWW_diameter_height_biomass_basal_area_sapwood_data.csv": "GreatWesternWoodlands",
	"Litchfield_Savanna_diameter_height_biomass_basal_area_data.csv": "Litchfield",
	"Robson_Creek_diameter_height_biomass_data_CuQoBBH.csv": "RobsonCreek",
	"Samford_diameter_height_biomass_data_mVFhel9.csv": "Samford",
	"Tumbarumba_Wet_Eucalypt_diameter_height_biomass_dat_IbRz0nd.csv": "Tumbarumba",
	"Warra_Tall_Eucalypt_diameter_height_biomass_data_teCJffC.csv": "Warra",
	"Whroo_Dry_Eucalypt_diameter_height_biomass_data_15iHyyj.csv": "Whroo",
	"Wombat_Stringybark_Eucalypt_diameter_height_biomass_EYriIKJ.csv": "WombatStateForest",
}

class Options:
	"""
	Class for storing CLI arguments from the user.

	@param log: Log level.
	@param files: Input files.
	@param odir: Output directory.
	@param prog: True to write progress messages, 0 otherwise.
	@param parallel: True to process files in parallel.
	@param timestep: Desired output timestep, in hours.
	@param smips_path: Path to file containing SW content at 90cm as obtained from SMIPS.
	@param smips_index_path: Path to file containing SMIPS index.
	@param biomass_file: Optional file containing biomass readings
	@param biomass_searchpath: Search path for biomass readings. If provided, will search under here for a file containing biomass readings for the site
	@param lai_file: Path to file containing LAI observations.
	@param greenness_file: Path to a .csv or .xlsx file containing greenness observations.
	@param inventory_path: Path to directory containing site-level inventory data.
	@param format: Output format (csv or netcdf).
	@param compression_level: Compression quality for output file [0, 9]. 0 = none, 1 = fastest compression, largest filesize, 9 = slowest compression, smallest filesize.
	@param compression_type: Compression algorithm to be used (default 'zlib').
	"""
	def __init__(self, log : LogLevel, files: list[str], odir: str, prog: bool,
		parallel: bool, timestep: int, smips_path: str, smips_index_path: str
		, biomass_file: str, biomass_searchpath: str, lai_file: str
		, greenness_file: str, inventory_path: str, format: str, compression_level: int
		, compression_type: str):
		self.log_level = log
		self.files = files
		self.out_dir = odir
		self.report_progress = prog
		self.parallel = parallel
		self.timestep = timestep
		self.biomass_file = biomass_file
		self.biomass_searchpath = biomass_searchpath
		self.lai_file = lai_file
		self.greenness_file = greenness_file
		self.smips_path = smips_path
		self.smips_index_path = smips_index_path
		self.inventory_path = inventory_path
		self.output_format = format
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
	parser.add_argument("-o", "--out-dir", required = True, help = "Path to the output directory. Processed files will be saved with the same file name into this directory")
	parser.add_argument("-p", "--show-progress", action = "store_true", help = "Report progress")
	parser.add_argument("-P", "--parallel", action = "store_true", help = "Process files in parallel")
	parser.add_argument("-t", "--timestep", type = int, required = True, help = "Output timestep in hours")
	parser.add_argument("-f", "--format", default = "netcdf", help = "Output format (allowed values are 'csv' or 'netcdf'). Defaults to netcdf if not specified.")
	parser.add_argument("--version", action = "version", version = "%(prog)s " + ozflux_common.VERSION)
	parser.add_argument("-l", "--lai-file", required = True, help = "Path to a file containing lai observations")
	parser.add_argument("-g", "--greenness-file", required = True, help = "Path to a .csv or .xlsx file containing greenness data")
	parser.add_argument("-s", "--smips-path", required = True, help = "Path to directory containing site-level SW content at 90cm as obtained from SMIPS. Each file in this directory should be named site.csv")
	parser.add_argument("-i", "--inventory-path", required = True, help = "Path to directory containing site-level inventory data. Each file should be named site.csv")
	parser.add_argument("--compression-level", type = int, nargs = "?", default = 5, help = "Compression quality for output file [0, 9]. 0 = none, 1 = fastest compression, largest filesize, 9 = slowest compression, smallest filesize (default 5)")
	parser.add_argument("--compression-type", default = "zlib", help = "Compression algorithm to be used (default 'zlib')")
	parser.add_argument("--smips-index-path", help = "Optional path to directory containing site-level soil moisture index as obtained from smips. Each file should be named site.csv")

	biomass_group = parser.add_mutually_exclusive_group()
	biomass_group.add_argument("-b", "--biomass-readings", required = False, help = "Optional file containing biomass readings")
	biomass_group.add_argument("--biomass-searchpath", required = False, help = "Search path for biomass readings. If provided, will search under here for a file containing biomass readings for the site")

	p = parser.parse_args(argv[1:])

	return Options(p.verbosity, p.files, p.out_dir, p.show_progress, p.parallel,
		p.timestep, p.smips_path, p.smips_index_path, p.biomass_readings,
		p.biomass_searchpath, p.lai_file, p.greenness_file, p.inventory_path,
		p.format, p.compression_level, p.compression_type)

def multiply_timedelta(delta: datetime.timedelta, n: int) -> datetime.timedelta:
	"""
	Return a timedelta which represents N multiples of the original delta.

	@param delta: A time delta.
	@param n: Any real number.
	"""
	return datetime.timedelta(seconds = n * delta.total_seconds())

def _create_observations(name: str, long_name: str, units: str, data: list[float]
		, start: datetime.datetime, timestep: datetime.timedelta):

	obs = [Observation(start + multiply_timedelta(timestep, i), data[i]) for (i, _) in enumerate(data)]
	return Observations(name, long_name, units, obs)

def write_csv(outfile: str, variables: list[Observations]
	#names: list[str], data: list[list[float]]
	#, biomass_data: list[biom_processing.BiomassReading]
	, start_date: datetime.datetime, end_date: datetime.datetime, timestep: int
	#, lai_data: list[Observation], sw90_data: list[Observation]
	, pcb: Callable[[float], None], delim = ","):
	"""
	Write the data to a csv file.

	@param outfile: Output file path.
	@param names: Output file column names.
	@param data: 2D data array where each element is a single column of data.
	@param start_date: First date in the data.
	@param timestep: Output timestep (in hours).
	@param lai_data: LAI observations. May be empty list if none available.
	"""

	write_hours = timestep % 24 != 0

	delta = datetime.timedelta(hours = timestep)
	nrow = int((end_date - start_date).total_seconds() // delta.total_seconds())

	variables = [v for v in variables if len(v.data) > 0]
	indices = [next(i for i in range(len(var.data)) if var.data[i].date >= start_date) for var in variables]

	with open(outfile, "w") as csv:
		csv.write("year%sdoy%s" % (delim, delim))
		if write_hours:
			csv.write("hour%s" % delim)
		csv.write("%s" % str.join(delim, [v.name for v in variables]))
		csv.write("\n")
		date = start_date
		for i in range(nrow):
			if i % PROGRESS_CHUNK_SIZE == 0:
				pcb(i / nrow)
			# %j is day of year (1-366), but lpj-guess writes day of year as
			# (0-365), so we subtract one here to put it in the same "units".
			doy = str(int(date.strftime("%j")) - 1)
			csv.write("%s%s%s%s" % (date.year, delim, doy, delim))

			# Write the hour column.
			if write_hours:
				csv.write("%s%s" % (date.hour, delim))

			for i in range(len(variables)):
				row = indices[i]
				var = variables[i]
				if row < len(var.data):
					reading = var.data[row]
					if reading.date == date:
						csv.write("%s" % reading.value)
						indices[i] += 1
					if i < (len(variables) - 1):
						csv.write(delim)
			csv.write("\n")
			date += delta
	pcb(1)

def write_netcdf(outfile: str, variables: list[Observations]
	, lat: float, lon: float, compression_level: int, compression_alg: str
	, lock: Lock, timestep: int, pcb: Callable[[float], None]):
	"""
	Write all observations to a netcdf file.

	@param outfile: Output file path.
	@param variables: List of columns of data.
	@param start_date: First date in the data.
	@param end_date: Last date in the data.
	@param timestep: Output timestep (in hours).
	@param lat: Latitude of this site.
	@param lon: Longitude of this site.
	@param lock: Mutex used to serialise access to the output .nc file.
	@param timestep: Output timestep in hours.
	@param pcb: Progress callback function.
	"""
	with lock:
		with open_netcdf(outfile, write = True) as nc:
			time = nc[VAR_TIME]
			time_values = time[:]
			calendar = getattr(time, ATTR_CALENDAR)
			date_values = num2date(time_values, time.units, calendar)
			first_date = min(date_values)
			last_date = max(date_values)

			(ilon, ilat) = get_coord_indices(nc, lon, lat)
			dates = numpy.sort(numpy.unique([d.date for v in variables for d in v.data]))
			nums = date2num(dates, time.units, time.calendar)
			allowed_dates = [n for n in nums if n in time_values]
			if len(allowed_dates) != len(nums):
				log_warning(f"Data contains dates outside of flux data range. Total number of dates: {len(nums)}, number of dates within ozflux range: {len(allowed_dates)}. Data points removed={len(nums) - len(allowed_dates)}")
			# time[:] = nums # dangerous!

			for obs_var in variables:
				if len(obs_var.data) == 0:
					log_warning(f"Observed variable {obs_var.name} has no data")
					continue
				log_diagnostic(f"Writing variable {obs_var.name}...")
				data = [d for d in obs_var.data if d.date >= first_date and d.date <= last_date]
				if len(data) != len(obs_var.data):
					log_warning(f"{obs_var.name}: Discarded {len(obs_var.data) - len(data)} data points which lie outside of ozflux timespan")
				times = [x.date for x in data]
				if timestep % HOURS_PER_DAY == 0:
					# If daily timestep, discard the time part of the times.
					times = [datetime.datetime(x.year, x.month, x.day) for x in times]
				dates = date2index(times, time, calendar)

				name = obs_var.name
				if not name in nc.variables:
					# Variable doesn't exist, create it.
					create_var_if_not_exists(nc, name, FORMAT_FLOAT, _VAR_DIMS
						, compression_level, compression_alg
						, (_TIME_CHUNK_SIZE, 1, 1))
					# Set metadata.
					# todo: std_name?
					nc.variables[name].units = obs_var.units
					setattr(nc.variables[name], ATTR_LONG_NAME, obs_var.long_name)
				var = nc.variables[name]
				values = [x.value for x in data]
				n = len(times)
				if n == 1:
					# date2index returns a scalar when n==1 for some reason.
					# I love this language so much.
					dates = [dates]
				for i in range(0, n, _TIME_CHUNK_SIZE):
					upper = min(n, i + _TIME_CHUNK_SIZE)
					var[dates[i:upper], ilon, ilat] = values[i:upper]
					pcb(upper / n)

def get_start_date(filename: str) -> datetime.datetime:
	"""
	Get the first date in the input file.

	@param filename: Path to a .nc file.
	"""
	with open_netcdf(filename) as nc:
		var = nc.variables[VAR_TIME]
		toffset = var[0]
		units = getattr(var, ATTR_UNITS)
		calendar = get_calendar(var)
		return num2date(toffset, units, calendar)

def get_end_date(filename: str) -> datetime.datetime:
	"""
	Get the first date in the input file.

	@param filename: Path to a .nc file.
	"""
	with open_netcdf(filename) as nc:
		var = nc.variables[VAR_TIME]
		toffset = var[len(var) - 1]
		units = getattr(var, ATTR_UNITS)
		calendar = get_calendar(var)
		return num2date(toffset, units, calendar)

def init_outfile(filename: str, infiles: list[str], timestep: int
		, vars: list[ForcingVariable], compression_level: int, compression_type: str
		, start_date: datetime.datetime, end_date: datetime.datetime):
	"""
	Initialise the output file by creating dimensions and variables as required.

	@param nc_out: The output file.
	@param first_infile: The first input file. All other input files must have
	identical variables and dimensionality to this file.
	@param compression_level: Level of compression (0-9).
	@param compression_type: Compression algorithm.
	"""
	if len(infiles) < 1:
		raise ValueError("No input files")

	# todo: this doesn't consider the possibility of input files containing
	# duplicate grid points, in which case we will create a dimension with a
	# fixed size larger than is strictly necessary.
	num_gridpoints = sum([count_gridpoints(file) for file in infiles])

	flux_start = min([get_start_date(file) for file in infiles])
	flux_end = max([get_end_date(file) for file in infiles])
	start_date = min(start_date, flux_start)
	end_date = max(end_date, flux_end)
	timestep_seconds = timestep * SECONDS_PER_HOUR
	ntime = math.ceil((end_date - start_date).total_seconds() / timestep_seconds) + 1

	# Dimension length must be an even number.
	ntime += ntime % 2
	if ntime < _TIME_CHUNK_SIZE:
		_time_chunk_size = ntime

	with open_netcdf(filename, write = True) as nc:
		create_dim_if_not_exists(nc, DIM_TIME)
		create_dim_if_not_exists(nc, DIM_LON, num_gridpoints)
		create_dim_if_not_exists(nc, DIM_LAT, num_gridpoints)
		create_var_if_not_exists(nc, VAR_TIME, FORMAT_UINT, (DIM_TIME,)
			, compression_level, compression_type, (_time_chunk_size,))
		create_var_if_not_exists(nc, _VAR_LON, FORMAT_FLOAT, (DIM_LON,)
			, compression_level, compression_type, (1,))
		create_var_if_not_exists(nc, _VAR_LAT, FORMAT_FLOAT, (DIM_LAT,)
			, compression_level, compression_type, (1,))
		var_lat = nc.variables[_VAR_LAT]
		var_lon = nc.variables[_VAR_LON]
		var_time = nc.variables[VAR_TIME]

		# Set attributes required by CF spec.

		setattr(var_lat, ATTR_STD_NAME, STD_LAT)
		setattr(var_lon, ATTR_STD_NAME, STD_LON)
		setattr(var_time, ATTR_STD_NAME, STD_TIME)

		setattr(var_lat, ATTR_LONG_NAME, STD_LAT)
		setattr(var_lon, ATTR_LONG_NAME, STD_LON)
		setattr(var_time, ATTR_LONG_NAME, STD_TIME)

		setattr(var_lat, ATTR_UNITS, "degree_north")
		setattr(var_lon, ATTR_UNITS, "degree_east")
		time_unit = "days" if timestep % HOURS_PER_DAY == 0 else "hours"
		setattr(var_time, ATTR_UNITS, f"{time_unit} since 1900-01-01 00:00:00")

		setattr(var_time, ATTR_CALENDAR, CALENDAR_GREGORIAN)

		# Create data variables.
		dims =   (DIM_TIME        , DIM_LON, DIM_LAT)
		chunks = (_time_chunk_size, 1      , 1      )
		for var in vars:
			create_var_if_not_exists(nc, var.out_name, FORMAT_FLOAT, dims
			    , compression_level, compression_type, chunks)

			with open_netcdf(infiles[0]) as nc_in:
				if var.in_name in nc_in.variables:
					var_in = nc_in.variables[var.in_name]
					var_out = nc.variables[var.out_name]
					copy_attr(var_in, var_out, ATTR_STD_NAME)
					copy_attr(var_in, var_out, ATTR_LONG_NAME)
			setattr(var_out, ATTR_UNITS, var.out_units)

		# Populate the time variable.
		dlt = datetime.timedelta(seconds = timestep_seconds)
		dates = [start_date + i * dlt for i in range(ntime)]
		# Don't write time data yet, as the exact time coverage depends on the
		# input data.
		nums = date2num(dates, var_time.units, var_time.calendar)
		var_time[0:ntime] = nums

def get_biomass_data(nc_file: str, biom_file: str, biom_searchpath: str,
	pcb: Callable[[float], None]) \
	-> list[biom_processing.BiomassReading]:
	"""
	Attempt to read biomass data corresponding to the site specified by the .nc
	file name.

	@param nc_file: Input .nc file (may be used to determine site name).
	@param biom_file: -b CLI argument from user (or None if not given).
	@param biom_searchpath: -s CLI argument from user (or None if not given).
	"""
	if biom_file != None:
		raw = biom_processing.read_raw_data(biom_file)
		return biom_processing.read_biomass_data(raw, pcb)

	if biom_searchpath == None:
		log_information(f"No biomass input data file was given. No biomass data will be included.")
		return []

	# Determine site name, convert to lower case.
	site_name = get_site_name(nc_file)
	site_name = site_name.replace(" ", "_").lower()

	if "CumberlandPlain" in nc_file:
		site_name = "cumberland_plains"

	site_dir = "%s/%s" % (biom_searchpath, site_name)
	dir = site_dir if os.path.exists(site_dir) else biom_searchpath

	for file in glob.glob("%s/*.csv" % dir):
		log_diagnostic("Attempting to read biomass data from '%s'..." % file)
		try:
			raw = biom_processing.read_raw_data(file)
			return biom_processing.read_biomass_data(raw, pcb)
		except BaseException as error:
			# log_warning("Failed to read biomass data from file '%s'. This is not necessarily a problem." % file)
			# log_warning(error)
			continue
	# Cannot find a .csv file containing valid biomass readings.
	return []

def get_layer_variables(file: str) -> list[Layer]:
	"""
	Get the names of variables representing layers in this 

	@param nc: The input .nc file name.
	"""
	site = get_site_name_from_filename(file)
	if not site in _layer_variable_map:
		raise ValueError("Unknown layer variables for site '%s'" % site)
	return _layer_variable_map[site]

def get_layer_height(nc: Dataset, layer: str, next_layer: str) -> tuple[float, float]:
	"""
	Get the height of the specified layer in m. (ie distance between top and
	bottom of the layer.) The return value is a tuple containing the top of the
	layer and the layer height in that order.

	@param nc: Input .nc file.
	@param layer: Name of a layered soil water variable in the .nc file. Must
	              have a height attribute.
	@param next_layer: Name of the nc variable containing the next layer, or
	                   None if this is the bottom layer. Must have a height
					   attribute.
	"""
	height = getattr(nc.variables[layer], _ATTR_HEIGHT)
	depths = ozflux_common.get_lengths(height)

	layer_top = depths[0]

	if len(depths) == 3:
		# Some sites oh so helpfully have a height that looks like:
		# "-0.08m vertically,0.08 - 0.38 m vertically"
		if ozflux_common.floats_equal(depths[0], depths[1]):
			depths = [depths[0], depths[2]]
		else:
			raise ValueError("Unable to interpret layer depths from height string '%s'" % height)

	if len(depths) == 2 and depths[0] > depths[1]:
		tmp = depths[0]
		depths[0] = depths[1]
		depths[1] = tmp

	if len(depths) == 2 and ozflux_common.floats_equal(depths[0], depths[1]):
		depths = [depths[0]]

	if len(depths) == 2:
		# Height was of the form 0.1m - 0.2m.
		layer_bottom = depths[1]
		return (layer_top, layer_bottom - layer_top)

	if next_layer == None:
		# Bottom layer - assume this layer extends to bottom of profile.
		if layer_top > _SOIL_DEPTH:
			# If we get to here, the configured layer names in the lookup table
			# are probably incorrect for this site.
			m = "Bottom layer starts at %.2fm, which is below the hardcoded soil profile depth (%.2fm), and has unknown height"
			raise ValueError(m % (layer_top, _SOIL_DEPTH))

		dlt = _SOIL_DEPTH - layer_top
		m = "Site '%s': extending bottom soil layer by %.2fm"
		log_warning(m % (nc.site_name, dlt) )
		return (layer_top, dlt)

	height_next_layer = getattr(nc.variables[next_layer], _ATTR_HEIGHT)
	next_layer_depths = ozflux_common.get_lengths(height_next_layer)

	top_next_layer = next_layer_depths[0]

	if top_next_layer < layer_top:
		m = "Unknown height of layer '%s' in site '%s'; layer starts at %.2fm but the next layer down ('%s') starts at %.2fm. This probably indicates an issue in the layer name lookup table for this site."
		raise ValueError(m % (layer, nc.site_name, layer_top, next_layer, top_next_layer))

	return (layer_top, top_next_layer - layer_top)

def get_sw_layer_name(depth: int):
	"""
	Get a name for a soil water layer at the specified depth, which will be used
	as a column name in the output file.

	@param depth: Depth in cm.
	"""
	return "%s_%d" % (_NAME_SW, depth)

def get_sw_data(file: str, timestep: int, pcb: Callable[[float], None]) \
-> list[Observations]:
	"""
	Read SW timeseries data from the given file name. Return a timeseries for
	each available layer at the specified output timestep.

	@param file: Path to the input .nc file.
	@param timestep: Desired output timestep in hours.
	@param pcb: Progress callback function used for progress reporting.
	"""
	# Get output timestep in minutes.
	timestep_min = timestep * MINUTES_PER_HOUR

	# Create an array containing running total of SW for each timestep.
	layers: list[Observations] = []

	with Dataset(file, "r", format=ozflux_common.NC_FORMAT) as nc:
		first_date = ozflux_common.parse_date(nc.time_coverage_start)
		start_date = get_next_year(first_date)

		# Determine which layer variables we should read for this site.
		layer_variables = get_layer_variables(file)
		step_start = 0
		step_size = 1 / len(layer_variables)
		for i in range(len(layer_variables)):
			layer = layer_variables[i]

			# Read volumetric SW content for this layer.
			var = ForcingVariable(layer.name, layer.name, "m3/m3", numpy.mean
			 	, 0, 1)
			layer_vol = get_data(nc, var, timestep_min, True, True
				, lambda p: pcb(step_start + step_size * p))

			long_name = f"Volumetric soil water content at {layer.depth}cm"

			delta_t = datetime.timedelta(hours = timestep)
			layer = _create_observations(get_sw_layer_name(layer.depth)
				, long_name, var.out_units, layer_vol, start_date, delta_t)
			layers.append(layer)

			step_start += step_size
	pcb(1)

	return layers

def bit_set(x: int, n: int) -> bool:
	"""
	Check if the N-th bit is set in x.

	@param x
	"""
	y = 1 << n
	return x & y == y

def _lai_qc_filter(qc: int) -> bool:
	"""
	Check if a QC flag indicates usable data.

	@param qc: A QC flag value.
	"""
	# See here for a complete technical guide to the product, which includes
	# documentation of the QC flag:
	# https://lpdaac.usgs.gov/documents/624/MOD15_User_Guide_V6.pdf

	# If bit 2 is set, dead detectors caused >50% adjacent detector retrieval.
	if bit_set(qc, 2):
		return False

	# Meaning of bit 3 depends on bit 4.
	if bit_set(qc, 3):
		# If bit 4 is also set, significant clouds were present. Otherwise,
		# cloud state not defined (ie unknown). Either way, we don't trust this
		# value.
		return False

	# If bit 4 set but bit 3 not set, mixed clouds were present. Treat this as
	# acceptable for now.

	# If bit 6 set, RT method failed.
	if bit_set(qc, 6):
		return False

	# If bit 7 set, pixel not produced at all, value couldn't be retrieved.
	if bit_set(qc, 7):
		return False

	return True

def lai_qc_filter(qcflags: pandas.Series) -> list[bool]:
	"""
	Filter out invalid QC flags from the input data.

	@param qcflags: Series of QC flags.
	"""
	return [_lai_qc_filter(int(qc)) for qc in qcflags]

def get_lai_data(data: pandas.DataFrame, site: str, timestep: int, pcb: Callable[[float], None]) -> list[Observation]:
	"""
	Read LAI observations for the specified site from the given input file.

	@param data: Pre-parsed LAI data. (This is parsed once because the input file contains data for all sites.)
	@param site: Name of the site.
	@param timestep: Desired output timestep in hours.
	@param pcb: Progress reporting function.
	"""
	data = data[data[_LAI_COL_ID] == site]
	data = data.groupby([_LAI_COL_DATE]).mean(numeric_only = True)
	data = data.iloc[lai_qc_filter(data[_LAI_COL_QCFLAG])]

	lai = [Observation(datetime.datetime.combine(date.date(), datetime.datetime.min.time()), float(row[_LAI_COL_LAI])) for (date, row) in data.iterrows()]
	lai = [l for l in lai if l.value >= _LAI_MIN and l.value <= _LAI_MAX]
	return lai

def get_greenness_data(data: pandas.DataFrame, lon: float, lat: float
		, pcb: Callable[[float], None]) -> list[Observation]:
	"""
	Read greenness observations for the specified site from the given input file.

	@param data: Pre-parsed greenness data. (This is parsed once because the input file contains data for all sites.)
	@param lon: Longitude of the site.
	@param lat: Latitude of the site.
	@param pcb: Progress callback function.
	"""
	data = data[data[_GREENNESS_COL_LON].apply(lambda x: floats_equal(x, lon))]
	data = data[data[_GREENNESS_COL_LAT].apply(lambda x: floats_equal(x, lat))]
	if len(data) < 1:
		return []
	data = data.groupby([_GREENNESS_COL_DATE]).mean(numeric_only = True)
	greenness = [
		Observation(date, float(row[_GREENNESS_COL_GREENNESS]))
		for (date, row) in data.iterrows()
	]
	return greenness

def get_inventory_data_file(inventory_path: str, site: str) -> str | None:
	"""
	Return the path to the inventory data file for the specified site.

	@param inventory_path: Path to a directory containing site-level csv files containing timeseries of inventory data.
	@param site: Site name.
	"""
	inventory_files = [f for f in os.listdir(inventory_path) if f.endswith(".csv")]
	inventory_files = [os.path.join(inventory_path, f) for f in inventory_files]
	for file in inventory_files:
		site_name = _inventory_site_names[os.path.basename(file)]
		if site_name == site:
			return file
	return None

def get_inventory_data(inventory_path: str, site: str) -> tuple[list[Observation], list[Observation], list[Observation], list[Observation]]:
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
		log_warning(f"No inventory data for site '{site}'")
		return ([], [], [], [])

	# Read data from disk.
	df = biom_processing.read_raw_data(file)
	inventory = biom_processing.read_inventory_data(df)

	# Initialise output lists.
	heights: list[Observation] = []
	diameters: list[Observation] = []
	live_biomass: list[Observation] = []
	dead_biomass: list[Observation] = []

	# Probably not very pythonic but it works...
	sorted_readings = sorted(inventory.readings, key=attrgetter("date"))
	grouped = groupby(sorted_readings, key=attrgetter("date"))
	for date, group in grouped:
		sorted_patches = sorted(group, key=attrgetter("patch"))
		patch_groups = groupby(sorted_patches, key=attrgetter("patch"))
		date_heights = []
		date_diameters = []
		date_live_biomass = []
		date_dead_biomass = []
		for patch, patch_group in patch_groups:
			# Convert iterator to list to avoid exhaustion
			patch_group_list = list(patch_group)
			# Get all values for this patch on this timestep.
			patch_heights = [r.height for r in patch_group_list if not math.isnan(r.height)]
			patch_diameters = [r.diameter for r in patch_group_list if not math.isnan(r.diameter)]
			patch_live_biomass = [r.live_biomass for r in patch_group_list if not math.isnan(r.live_biomass)]
			patch_dead_biomass = [r.dead_biomass for r in patch_group_list if not math.isnan(r.dead_biomass)]

			# Store mean height/diameter (m), and total biomass (in kg/m2) for
			# this patch.
			mean_height = math.nan if len(patch_heights) == 0 else numpy.mean(patch_heights) # m
			mean_diameter = math.nan if len(patch_diameters) == 0 else numpy.mean(patch_diameters) # m
			total_live_biomass = math.nan if len(patch_live_biomass) == 0 else numpy.sum(patch_live_biomass) / inventory.get_area(patch) # kg/m2
			total_dead_biomass = math.nan if len(patch_dead_biomass) == 0 else numpy.sum(patch_dead_biomass) / inventory.get_area(patch) # kg/m2
			if not math.isnan(mean_height):
				date_heights.append(mean_height)
			if not math.isnan(mean_diameter):
				date_diameters.append(mean_diameter)
			if not math.isnan(total_live_biomass) and total_live_biomass > 0:
				date_live_biomass.append(total_live_biomass)
			if not math.isnan(total_dead_biomass) and total_dead_biomass > 0:
				date_dead_biomass.append(total_dead_biomass)

		# Mean value over all patches on this timestep. This should be
		# comparable to gridcell-level model outputs (e.g. cmass.out).
		mean_height = math.nan if len(date_heights) == 0 else numpy.mean(date_heights)
		mean_diameter = math.nan if len(date_diameters) == 0 else numpy.mean(date_diameters)
		mean_live_biomass = math.nan if len(date_live_biomass) == 0 else numpy.mean(date_live_biomass)
		mean_dead_biomass = math.nan if len(date_dead_biomass) == 0 else numpy.mean(date_dead_biomass)

		if not math.isnan(mean_height):
			heights.append(Observation(date, mean_height))
		if not math.isnan(mean_diameter):
			diameters.append(Observation(date, mean_diameter))
		if not math.isnan(mean_live_biomass):
			live_biomass.append(Observation(date, mean_live_biomass))
		if not math.isnan(mean_dead_biomass):
			dead_biomass.append(Observation(date, mean_dead_biomass))

	return (heights, diameters, live_biomass, dead_biomass)

def get_smips_data_from_file(file: str, timestep: int, col: str, pcb: Callable[[float], None]) \
		-> list[Observation]:
	"""
	Read smips SW content at 90cm from the site file in the specified path.

	@param file: Input csv file.
	@param timestep: Output timestep in hours.
	@param pcb: Progress reporting function.
	"""
	data = pandas.read_csv(file, parse_dates = True)
	data.iloc[:,0] = [datetime.datetime.strptime(x, "%Y-%m-%d") for x in data.iloc[:,0]]
	parsed = [Observation(row[_SMIPS_COL_DATE], row[col]) for row in data.iloc]
	return parsed

def get_smips_data(site: str, path: str, timestep: int, col: str, pcb: Callable[[float], None]) \
		-> list[Observation]:
	"""
	Read smips SW content at 90cm from the site file in the specified path.

	@param site: Name of the site.
	@param path: Path to a directory containing site-level csv files containing timeseries of SW content at 90cm.
	@param timestep: Output timestep in hours.
	@param col: Column name of data to be read.
	@param pcb: Progress reporting function.
	"""
	file = os.path.join(path, "%s.csv" % site)
	return get_smips_data_from_file(file, timestep, col, pcb)

def process_file(file: str, out_dir: str, timestep: int, biomass_file: str
	, biomass_searchpath: str, lai_data: pandas.DataFrame
	, greenness_data: pandas.DataFrame, smips_path: str, smips_index_path: str
	, inventory_path: str
	, lock: Lock, pcb: Callable[[float], None]):
	"""
	Extract the standard variables from the input file, write them to a .csv
	file in the specified output directory in the specified timestep.

	@param file: Input file path.
	@param out_dir: Output directory.
	@param timestep: Output timestep in hours.
	@param smips_path: Path to smips input files.
	@param smips_path: Path to smips SW index files.
	@param pcb: Progress callback function.
	@param lai_file: Path to a file containing LAI observations.
	@param lock: Mutex used to serialise access to the output .nc file.
	"""
	log_information("Processing '%s'..." % file)

	read_prop = 0.76
	write_prop = 0.03
	biom_prop = 0.0024
	sw_prop = 0.2
	smips_prop = 0.002
	lai_prop = 1 - write_prop - read_prop - biom_prop - sw_prop - smips_prop
	if lai_prop < 0:
		raise ValueError("Invalid default time proportions")

	flux_vars = _standard_variables

	out_variables: list[Observations] = []

	step_start = 0

	lai_start = time.time()
	site = get_site_name(file)
	lai_data = get_lai_data(lai_data, site, timestep
			, lambda p: pcb(step_start + lai_prop * p))
	out_variables.append(Observations(_NAME_LAI, "Leaf Area Index", "m2/m2", lai_data))
	lai_time = time.time() - lai_start

	log_debug("Opening input file '%s' for reading..." % file)

	step_start += lai_prop

	read_start = time.time()
	timestep_minutes = timestep * MINUTES_PER_HOUR
	start_date: datetime.datetime
	end_date: datetime.datetime
	delta = datetime.timedelta(hours = timestep)
	timeseries: list[datetime.datetime] = []
	latitude: float = 0
	longitude: float = 0

	site_canonical = get_site_name_from_filename(file)
	(heights, diameters, live_biomass, dead_biomass) = get_inventory_data(inventory_path, site_canonical)
	log_diagnostic(f"{site_canonical}: Retrieved {len(heights)} height observations, {len(diameters)} diameter observations, {len(live_biomass)} live biomass observations, and {len(dead_biomass)} dead biomass observations")
	if len(heights) > 0:
		height_data = Observations(_NAME_HEIGHT, "Stem height", "m", heights)
		out_variables.append(height_data)

	if len(diameters) > 0:
		diameter_data = Observations(_NAME_DIAMETER, "Stem diameter", "m", diameters)
		out_variables.append(diameter_data)

	if len(live_biomass) > 0:
		live_data = Observations(_NAME_LIVE, "Above-Ground Live biomass", "kg/m2", live_biomass)
		out_variables.append(live_data)

	if len(dead_biomass) > 0:
		dead_data = Observations(_NAME_DEAD, "Above-Ground Dead biomass", "kg/m2", dead_biomass)
		out_variables.append(dead_data)

	with open_netcdf(file) as nc:
		latitude = nc.variables["latitude"][0]
		longitude = nc.variables["longitude"][0]
		step_size = read_prop / len(flux_vars)
		start_date = get_next_year(ozflux_common.parse_date(nc.time_coverage_start))
		for var in flux_vars:
			if not var.in_name in nc.variables:
				raise ValueError(f"Variable {var.in_name} does not exist in file {file}")
			variable = nc.variables[var.in_name]
			if not hasattr(variable, ATTR_UNITS):
				raise ValueError(f"Variable {var.in_name} in {file} has no units")
			units = getattr(variable, ATTR_UNITS)
			# in_file: Dataset
			# var: ForcingVariable
			# output_timestep: int
			# qc_filter: bool
			# nan_filter: bool
			# progress_cb: Callable[[float], None]
			#
			# in_file, var.in_name, output_timestep, var.out_units, var.aggregator, var.lbound, var.ubound, var.invert,
			# qc_filter, nan_filter, progress_cb
			data = get_data_timeseries(nc, var.in_name, timestep_minutes,
							  		   var.out_units, var.aggregator,
									   var.lbound, var.ubound, var.invert, True,
									   False, lambda p: \
									   pcb(step_start + step_size * p))
			if len(timeseries) == 0:
				# timeseries = [start_date + multiply_timedelta(delta, i) for i in range(len(raw))]
				# end_date = start_date + multiply_timedelta(delta, len(raw))
				timeseries = data.index
				end_date = timeseries[len(timeseries) - 1]
			parsed = [Observation(x, y) for x, y in data.itertuples()]
			long_name = var.out_name
			if hasattr(variable, ATTR_LONG_NAME):
				long_name = getattr(variable, ATTR_LONG_NAME)
			out_variables.append(Observations(var.out_name, long_name, units, parsed))
			step_start += step_size
	read_time = time.time() - read_start

	greenness_data = get_greenness_data(greenness_data, longitude, latitude
		, lambda _: ...)
	greenness = Observations(_VAR_GREENNESS, "NDVI", _GREENNESS_UNITS, greenness_data)
	out_variables.append(greenness)

	sw_start = time.time()
	sw_layers = get_sw_data(file, timestep
		    , lambda p: pcb(step_start + sw_prop * p) )
	for sw_layer in sw_layers:
		out_variables.append(sw_layer)
	sw_time = time.time() - sw_start
	step_start += sw_prop

	smips_start = time.time()
	site_name = get_site_name_from_filename(file)
	sw90_data = get_smips_data(site_name, smips_path, timestep, _SMIPS_COL_SW
			    , lambda p: pcb(step_start + 0.5 * smips_prop * p))
	out_variables.append(Observations(_NAME_SMIPS, _SMIPS_COL_SW, "?mm?", sw90_data))
	if os.path.exists(smips_index_path):
		swindex_data = get_smips_data(site_name, smips_index_path, timestep
				, _SMIPS_COL_SWINDEX
				, lambda p: pcb(step_start + smips_prop * (0.5 * p + 0.5)))
		out_variables.append(Observations(_NAME_SWINDEX, _SMIPS_COL_SWINDEX, "0-1", swindex_data))
	smips_time = time.time() - smips_start
	step_start += smips_prop

	biom_start = time.time()
	biomass_data: list[biom_processing.BiomassReading] = []
	if biomass_file is not None or biomass_searchpath is not None:
		biomass_data = get_biomass_data(file, biomass_file, biomass_searchpath
			, lambda p: pcb(step_start + biom_prop * p))
		biomass_data = [b for b in biomass_data if b.date >= start_date and b.date <= end_date]
		biomass_data = sorted(biomass_data, key = lambda b: b.date)
	if len(biomass_data) > 0:
		live_data = Observations(_NAME_LIVE, "Live Biomass", "kg/ha", [Observation(d.date, d.live) for d in biomass_data])
		dead_data = Observations(_NAME_DEAD, "Dead Biomass", "kg/ha", [Observation(d.date, d.dead) for d in biomass_data])
		out_variables.append(live_data)
		out_variables.append(dead_data)

	biom_time = time.time() - biom_start
	step_start += biom_prop

	# step_start += read_prop

	file = os.path.basename(file)
	ext = "csv" if opts.output_format.lower() == "csv" else "nc"
	site = get_site_name_from_filename(file)
	filename = f"{site}.{ext}"
	outfile = os.path.join(out_dir, filename)

	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	write_start = time.time()

	timestep_hours = timestep
	if opts.output_format.lower() == "csv":
		write_csv(outfile, out_variables, start_date, end_date, timestep_hours
		   , lambda p: pcb(step_start + write_prop * p))
	else:
		outfile = get_netcdf_output_filename(out_dir)
		write_netcdf(outfile, out_variables, latitude, longitude
	       , opts.compression_level, opts.compression_type, lock, timestep
		   , lambda p: pcb(step_start + write_prop * p))
	write_time = time.time() - write_start

	total_time = read_time + write_time + biom_time + sw_time + lai_time + smips_time

	read_prop = read_time / total_time
	write_prop = write_time / total_time
	biom_prop = biom_time / total_time
	sw_prop = sw_time / total_time
	lai_prop = lai_time / total_time
	smips_prop = smips_time / total_time

	log_diagnostic("Read time: %.2f%%; Write time: %.2f%%; Biomass time: %.2f%%; SW time: %.2f%%; lai time: %.2f%%; smips time: %.2f%%"
	% (100 * read_prop, 100 * write_prop, 100 * biom_prop, 100 * sw_prop, 100 * lai_prop, 100 * smips_prop))

class Processor:
	def __init__(self, file: str, out_dir: str, timestep: int, biomass_file: str, biomass_searchpath: str
	      , lai_data: pandas.DataFrame, greenness_data: pandas.DataFrame
		  , smips_path: str, smips_index_path: str, inventory_path: str, lock: Lock):
		self.file = file
		self.out_dir = out_dir
		self.timestep = timestep
		self.biomass_file = biomass_file
		self.biomass_searchpath = biomass_searchpath
		self.lai_data = lai_data
		self.smips_path = smips_path
		self.smips_index_path = smips_index_path
		self.inventory_path = inventory_path
		self.lock = lock
		self.greenness_data = greenness_data
	def exec(self, pcb: Callable[[float], None]):
		process_file(self.file, self.out_dir, self.timestep, self.biomass_file
	       , self.biomass_searchpath, self.lai_data, self.greenness_data
		   , self.smips_path, self.smips_index_path, self.inventory_path
		   , self.lock, pcb)

def get_netcdf_output_filename(outpath: str) -> str:
	return os.path.join(outpath, "ozflux-obs.nc")

def read_csv_or_excel_data(file_path: str, cols: list[str], date_col: str
						   , date_fmt: str) -> pandas.DataFrame:
	"""
	Read an excel (.xlsx) or .csv file.
	"""
	if not os.path.exists(file_path):
		log_error(f"File does not exist: {file_path}")
	ext = os.path.splitext(file_path)[1]
	if ext.lower() == ".csv":
		data = pandas.read_csv(file_path, usecols = cols, parse_dates = False)
	else:
		data = pandas.read_excel(file_path, _SHEET_LAI, usecols = cols
								 , parse_dates = False)
	data[date_col] = pandas.to_datetime(data[date_col], format = date_fmt)
	return data

def read_lai_data(lai_file: str) -> pandas.DataFrame:

	log_information("Reading LAI data...")
	cols = [_LAI_COL_ID, _LAI_COL_DATE, _LAI_COL_LAI, _LAI_COL_QCFLAG]
	# _DATE_FMT = "%d/%m/%Y" # 28/07/2002
	_DATE_FMT = "%m/%d/%Y" # 07/28/2002
	return read_csv_or_excel_data(lai_file, cols, _LAI_COL_DATE, _DATE_FMT)

def read_greenness_data(greenness_file: str) -> pandas.DataFrame:
	log_information("Reading greenness data...")
	cols = [
		_GREENNESS_COL_DATE,
		_GREENNESS_COL_LON,
		_GREENNESS_COL_LAT,
		_GREENNESS_COL_GREENNESS
	]
	return read_csv_or_excel_data(greenness_file, cols, _GREENNESS_COL_DATE
		, _GREENNESS_DATE_FMT)
	

def main(opts: Options):
	"""
	Main function.

	@param opts: Parsed CLI options provided by the user..
	"""

	if not os.path.exists(opts.out_dir):
		os.makedirs(opts.out_dir)

	lai_data = read_lai_data(opts.lai_file)
	greenness_data = read_greenness_data(opts.greenness_file)

	if opts.output_format.lower() == "netcdf":
		outfile = get_netcdf_output_filename(opts.out_dir)
		if os.path.exists(outfile):
			os.remove(outfile)
		lai_start = min(lai_data[_LAI_COL_DATE])
		lai_end = max(lai_data[_LAI_COL_DATE])
		greenness_start = min(greenness_data[_GREENNESS_COL_DATE])
		greenness_end = max(greenness_data[_GREENNESS_COL_DATE])
		start_date = min(lai_start, greenness_start)
		end_date = max(lai_end, greenness_end)
		vars = [v for v in _standard_variables]
		vars.append(ForcingVariable("lai", "lai", "m2/m2", lambda x: x, _LAI_MIN, _LAI_MAX))
		init_outfile(outfile, opts.files, opts.timestep, vars
	      , opts.compression_level, opts.compression_type, start_date, end_date)

	job_manager = JobManager()

	lock = Lock()
	for file in opts.files:
		p = Processor(file, opts.out_dir, opts.timestep, opts.biomass_file
				, opts.biomass_searchpath, lai_data, greenness_data
				, opts.smips_path, opts.smips_index_path, opts.inventory_path
				, lock)
		job_manager.add_job(p, os.path.getsize(file))

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
