import datetime, numpy, os, re, threading, time
import multiprocessing, multiprocessing.connection
import traceback, sys
from ozflux_common import *
from ozflux_logging import *
from netCDF4 import Dataset, Variable, Dimension, default_fillvals, num2date
from enum import IntEnum
from sys import argv
from typing import Callable

# Name of the time dimension in the output file.
DIM_TIME = "time"

# Name of the longitude dimension in the output file.
DIM_LON = "longitude"

# Name of the latitude dimension in the output file.
DIM_LAT = "latitude"

# Inputs/outputs will be read/written in chunks of this size.
# Increasing this will improve performance but increase memory usage.
CHUNK_SIZE = 16384

# Progress updates inside long loops are written every N iterations.
PROGRESS_CHUNK_SIZE = 1024

# Whether to aggregate timesteps in paralllel. This doesn't save a huge amount
# of time in the files that I've tested (which are fairly small), but it doesn't
# hurt performance, and should really help with larger files. I don't really see
# any reason to set this to false, which is why this isn't a CLI option.
PARALLEL_AGGREGATION = True

# If this is true, the latitude/longitude dimensions will be created unlimited
# in size. This is not recommended, as it slows down reading of the NetCDF. When
# this is false, the dimensions will have exactly the length required, which is
# the number of input files. (This assumes 1 grid point per input file.)
UNLIMITED_DIMS = False

# Data format in output file.
FORMAT_FLOAT = "f8"

# Data format for single-precision floating point numbers.
FORMAT_SINGLE = "f4"

# Data format of unsigned long (uint64_t) in the output file.
FORMAT_UINT = "u8"

# Name of the time variable in the input files.
VAR_TIME = "time"

# Standard names defined by the CF spec. See here for all defined values:
# https://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html

# Name of the 'standard name' attribute (from the CF spec).
ATTR_STD_NAME = "standard_name"

# Name of the 'long name' attribute (from the CF spec).
ATTR_LONG_NAME = "long_name"

# Name of the 'units' attribute (from the CF spec).
ATTR_UNITS = "units"

# Name of the 'calendar' attribute (from the CF spec).
ATTR_CALENDAR = "calendar"

# Name of the 'missing value' attribute (from the CF spec).
_ATTR_MISSING_VAL = "missing_value"

# Name of the gregorian value of the calendar attribute (from the CF spec).
CALENDAR_GREGORIAN = "gregorian"

# Standard name of the latitude variable (from the CF spec).
STD_LAT = "latitude"

# Standard name of the longitude variable (from the CF spec).
STD_LON = "longitude"

# Standard name of the time variable (from the CF spec).
STD_TIME = "time"

# Global estimates of how much time it takes to perform various tasks (as a
# proportion of total time spent in get_data() function [0, 1]). These get
# updated during program execution and are used for progress reporting.
get_data_read_prop = 0.5
get_data_fixnan_prop = 0.05
get_data_t_agg_prop = 0.35
get_data_unit_prop = 0.05
get_data_bounds_prop = 0.05

# Global time estimates for reading vs writing data, as a proportion of time
# spent reading (which includes aggregation, unit conversion, etc) vs writing.
global_read_prop = 0.95
global_write_prop = 0.05

# Global mutex
mutex = multiprocessing.BoundedSemaphore(1)

# get_data() time proportion lock. Used to synchronise access to the
# get_data_*_prop variables.
gd_tp_lock = multiprocessing.BoundedSemaphore(1)

# copy_data_*() time proportion lock. Used to synchronise access to the
# global_*_prop variables.
cd_tp_lock = multiprocessing.BoundedSemaphore(1)

# Common alternative names for various units which we may reasonably encounter.
_units_synonyms = [
	["W/m2", "W/m^2", "Wm^-2"],
	["kg/m2/s", "kg/m^2/s", "kgm^-2s^-1", "kg m-2 s-1"],
	["K", "k"],
	["m/s", "ms^-1"],
	["Pa", "pa"],
	["kg/kg", "", "mm/mm", "m/m"], # debatable
	["ppm", "umol/mol"],
	["degC", "°C", "degrees C"],
	["umol/m2/s", "umol/m^2/s"],
	["m3/m3", "m^3/m^3"]
]

#
# Recipes for unit conversions.
#
# The keys are a tuple of two strings:
#
# 1. The input units
# 2. The output units.
#
# The values are functions which convert a value from the input units into the
# output units. These functions take two arguments:
#
# x: The input value in input units.
# t: The timestep length (in seconds).
#
# Note that units are changed after timestep aggregation occurs. So the input
# variable will already be in the output timestep at this point in time.
units_conversions = {
	("g", "kg"): lambda x, _: x * KG_PER_G, # (as an example)
	("mm", "kg/m2/s"): lambda x, t: x / t, # Assuming mm means mm per timestep
	("kg/m2/s", "mm"): lambda x, t: t * x,
	("kg/m2/s", "mm/day"): lambda x, t: SECONDS_PER_DAY * x,
	("degC", "K"): lambda x, _: x + DEG_C_TO_K,
	("K", "degC"): lambda x, _: x - DEG_C_TO_K,
	("kPa", "Pa"): lambda x, _: x * PA_PER_KPA,
	("umol/m2/s", "kgC/m2/day"): lambda x, _: x * MOL_PER_UMOL * G_C_PER_MOL * KG_PER_G * SECONDS_PER_DAY,
	("umol/m2/s", "gC/m2/day"): lambda x, _: x * MOL_PER_UMOL * G_C_PER_MOL * SECONDS_PER_DAY
}

_qc_flag_definitions = {
	0: (True, "Data has passed all QC checks"),
	1: (False, "Data missing from L1 Excel spreadsheet"),
	2: (False, "Failed range check"),
	3: (False, "Failed CSAT check, Diag_CSAT flag!=0 (do_CSATcheck)"),
	4: (False, "Failed 7500 check, Diag_7500 flag!=0 plus dependencies on AGC_7500, Ah_7500_Sd, Cc_7500_Sd, AhAh, CcCc (do_7500check)"),
	5: (False, "Failed diurnal check"),
	6: (False, "Date/time range excluded"),
	7: (False, "Hour range excluded"),
	8: (False, "Value of -9999 found with QC flag value of 0"),
	9: (False, "Fre set to missing when (Fsd>threshold) or (ustar<threshold)"),
	10: (True, "Linear correction or linear drift correction applied to data"),
	11: (False, "Dependent data rejected during 2D coordinate rotation"),
	12: (False, "Data rejected due to dependencies when calculating Massman frequency corrections (MassmanStandard)"),
	13: (False, "Fh rejected during conversion of Fhv to Fh (FhvtoFh)"),
	14: (False, "Fc rejected during WPL correction due to dependencies (Fc_WPLcov/Fe_WPL/Fe_WPLcov)"),
	15: (True, "Ta calculated from CSAT Tv rejected due to dependencies (TaFromTv)"),
	16: (False, "Data rejected at L3 due to failed range check (do_qcchecks)"),
	17: (False, "Data rejected at L3 due to failed diurnal check (do_qcchecks)"),
	18: (False, "Ustar below threshold (FilterUstar)"),
	19: (False, "Data rejected during coordination of gaps in flux series"),
	20: (True, "GapFilling: Driver gap filled using alternate data CONFLICT: same code used by ReplaceRotatedCovariance"),
	21: (True, "Missing rotated covariance replaced with non-rotated value CONFLICT: same code used in gfalternate_main"),
	22: (True, "Soil moisture set to default value in CorrectFgForStorage"),
	30: (True, "GapFilling: Flux Gap Filled by ANN (SOLO)"),
	31: (True, "GapFilling: Flux Gap not Filled by ANN"),
	35: (False, "Data replaced by alternate value when difference between data and alternate value exceeds threshold (ReplaceWhenDiffExceedsRange)"),
	38: (False, "Data rejected at L4 due to failed range check (do_qcchecks)"),
	39: (False, "Data rejected at L4 due to failed diurnal check (do_qcchecks)"),
	40: (False, "Gap filled from climatology"),
	41: (False, "No value available from GapFillUsingClimatology"),
	50: (False, "Gap filled by interpolation"),
	51: (False, "Fsd below threshold when calculating albedo"),
	52: (False, "Hour outside range 1000 to 1400 when calculating albedo"),
	60: (False, "Flux data generated by GapFillFluxFromDayRatio"),
	61: (False, "Stomatal resistance less than 0 (get_stomatalresistance)"),
	62: (False, "Fe less than threshold when calculating stomatal resistance (get_stomatalresistance)"),
	63: (False, "Fsd less than threshold when calculating stomatal resistance (get_stomatalresistance)"),
	64: (False, "Wind speed = 0 when calculating stomatal resistance (get_stomatalresistance)"),
	70: (False, "Partitioning Night: Re computed from exponential temperature response curves"),
	80: (False, "Partitioning Day: GPP/Re computed from light-response curves, GPP = Re - Fc"),
	81: (False, "Partitioning Day: GPP night mask"),
	82: (False, "Partitioning Day: Fc > Re, GPP = 0, Re = Fc"),
}

class ForcingVariable():
	# Create a new ForcingVariable instance.
	# @param invert: If true, the variable will be multiplied by -1.
	# @param lbound: Lower bound of the variable (in output units). Any data less than this value will be set to this value.
	# @param lbound: Upper bound of the variable (in output units). Any data greater than this value will be set to this value.
	def __init__(self, in_name: str, out_name: str, out_units: str
		, aggregator: Callable[[list[float]], float], lbound: float
		, ubound: float, invert: bool = False):
		self.in_name = in_name
		self.out_name = out_name
		self.out_units = out_units
		self.aggregator = aggregator
		self.lbound = lbound
		self.ubound = ubound
		self.invert = invert
	def __repr__(self):
		return f"{self.in_name} -> {self.out_name} ({self.out_units})"

def zeroes(n: int) -> list[float]:
	"""
	Create an array of zeroes of the specified length.
	"""
	return [0.0] * n

def get_units(var_id: Variable) -> str:
	"""
	Get the units for the specified variable in the .nc file.
	"""
	return var_id.units

def get_data(in_file: Dataset \
	, var: ForcingVariable \
	, output_timestep: int \
	, qc_filter: bool \
	, nan_filter: bool \
	, progress_cb: Callable[[float], None]) -> list[float]:
	"""
	Get all data from the input file and convert it to a format suitable for the
	output file.

	@param in_file: Input .nc file.
	@param var: Describes how variable is to be read/aggregated.
	@param progress_cb: Progress callback function.
	@param qc_filter: Iff true, data with invalid QC flags will be treated as NaN.
	@param nan_filter: Iff true, NaN data will be interpolated from neighbouring values.
	"""
	return _get_data(in_file, var.in_name, output_timestep, var.out_units
	, var.aggregator, var.lbound, var.ubound, var.invert, qc_filter, nan_filter
	, progress_cb)

def _get_qc_var_name(name: str):
	"""
	Return the name of the QF flag variable corresponding to the given variable.
	@param name: Name of a NetCDF file variable.
	"""
	return "%s_QCFlag" % name

def _is_good_qc(flag: float) -> bool:
	"""
	Return true iff the given QC flag indicates good data.
	@param flag: A QC flag.
	"""
	if flag in _qc_flag_definitions:
		(result, msg) = _qc_flag_definitions[flag]
		if not result:
			log_diagnostic("Filtering QC flag '%s'; reason: %s" % msg)
		return result
	return True

def _filter_qc(data: list[float], qc: list[float]
	, pcb: Callable[[float], None]) -> list[float]:
	"""
	Mask all data which have a corresponding QC flag which indicates bad/invalid
	data.

	@param data: Data to be filtered.
	@param qc: List of same length as data containing QC flags.
	@param pcb: Progress reporting function.
	"""
	if len(data) != len(qc):
		raise ValueError("Unable to filter from QC flag: data length mismatch")
	for i in range(len(data)):
		if i % PROGRESS_CHUNK_SIZE == 0:
			pcb(i / len(data))
		qc_flag = int(qc[i])
		if qc_flag in _qc_flag_definitions:
			(result, msg) = _qc_flag_definitions[qc_flag]
			if not result:
				data[i].mask = True
				log_diagnostic("Filtering QC flag %d: %s" % (i, msg))
	pcb(1)
	return data

def _get_data(in_file: Dataset \
	, in_name: str
	, output_timestep: int \
	, out_units: str \
	, aggregator: Callable[[list[float]], float] \
	, lower_bound: float \
	, upper_bound: float \
	, invert: bool \
	, qc_filter: bool \
	, nan_filter: bool \
	, progress_cb: Callable[[float], None]) -> list[float]:
		"""
		Get all data from the input file and convert it to a format
		suitable for the output file.

		@param in_name: Name of the variable in the input file.
		@param out_name: Name of the variable in the output file.
		@param output_timestep: Output timestep width in minutes.
		@param qc_filter: Iff true, data with invalid QC flags will be treated as NaN.
		@param nan_filter: Iff true, NaN data will be interpolated from neighbouring values.
		"""
		# Some variables are not translated directly into the output
		# file and can be ignored (ie filled with zeroes).
		if in_name == "zero":
			log_diagnostic("Using zeroes for %s" % in_name)
			input_timestep = int(in_file.time_step)
			timestep_ratio = output_timestep // input_timestep
			n = in_file.variables[VAR_TIME].size // timestep_ratio
			data = zeroes(n)
			return trim_to_start_year(in_file, output_timestep, data)
		if not in_name in in_file.variables:
			m = "Variable %s does not exist in input file"
			raise ValueError(m % in_name)

		var_id = in_file.variables[in_name]

		# Get units of the input variable.
		in_units = get_units(var_id)

		# Check if a unit conversion is required.
		matching_units = units_match(in_units, out_units)

		global get_data_read_prop, get_data_fixnan_prop, get_data_t_agg_prop
		global get_data_unit_prop, get_data_bounds_prop

		units_time_proportion = 0 if matching_units else get_data_unit_prop

		prop_tot = get_data_read_prop + get_data_fixnan_prop + \
			get_data_t_agg_prop + units_time_proportion + get_data_bounds_prop

		read_time_proportion = get_data_read_prop / prop_tot

		fixnan_time_proportion = get_data_fixnan_prop / prop_tot if nan_filter else 0

		aggregation_time_proportion = get_data_t_agg_prop / prop_tot

		bounds_time_prop = get_data_bounds_prop / prop_tot

		units_time_proportion /= prop_tot

		step_start = 0

		step_size = read_time_proportion

		qcvar_name = _get_qc_var_name(in_name)
		if not qcvar_name in in_file.variables:
			if qc_filter:
				msg = "Unable to QC filter variable '%s': QC variable '%s' does not exist in file"
				log_warning(msg % (in_name, qcvar_name))
			qc_filter = False

		if qc_filter:
			step_size /= 2

		# Read data from netcdf file.
		read_start = time.time()
		data = read_data(var_id, lambda p:progress_cb(step_size * p))
		read_tot = time.time() - read_start
		log_debug("Successfully read %d values from netcdf" % len(data))

		step_start += step_size

		if qc_filter:
			qcvar = in_file.variables[qcvar_name]
			step_size /= 2
			qc_data = read_data(qcvar, lambda p: progress_cb(step_start + step_size * p))
			step_start += step_size
			data = _filter_qc(data, qc_data, lambda p: progress_cb(step_start + step_size * p))
			step_start += step_size

		# Replace NaN with a mean of nearby values.
		log_diagnostic("Removing NaNs")
		step_start = read_time_proportion
		step_size = fixnan_time_proportion
		fixnan_start = time.time()
		if nan_filter:
			data = remove_nans(data
			, lambda p: progress_cb(step_start + step_size * p))
			step_start += fixnan_time_proportion
		fixnan_tot = time.time() - fixnan_start

		# Convert to flat array.
		data = numpy.array(data).flatten()

		# Change timestep to something suitable for lpj-guess.
		log_diagnostic("Aggregating %s to output timestep" % in_name)
		t_agg_start = time.time()
		data = temporal_aggregation(in_file, data, output_timestep, aggregator
		, lambda p: progress_cb(step_start + aggregation_time_proportion * p))
		step_start += aggregation_time_proportion
		t_agg_tot = time.time() - t_agg_start

		# Trim data to the start of a year.
		data = trim_to_start_year(in_file, output_timestep, data)

		# Ensure units are correct.
		unit_start = time.time()
		if not matching_units:
			log_diagnostic("Converting %s from %s to %s" % \
				(in_name, in_units, out_units))
			data = fix_units(data, in_units, out_units, output_timestep, invert,
				lambda p: progress_cb( \
					step_start + units_time_proportion * p))
		else:
			log_diagnostic("No units conversion required for %s" % in_name)
		unit_tot = time.time() - unit_start

		step_start += units_time_proportion

		bounds_start = time.time()
		data = bounds_checks(data, lower_bound, upper_bound, lambda p: \
			progress_cb(step_start + bounds_time_prop * p))
		bounds_tot = time.time() - bounds_start

		time_tot = read_tot + fixnan_tot + t_agg_tot + bounds_tot + unit_tot
		read_prop = read_tot / time_tot
		fixnan_prop = fixnan_tot / time_tot
		t_agg_prop = t_agg_tot / time_tot
		bounds_prop = bounds_tot / time_tot
		unit_prop = unit_tot / time_tot

		# Update global time estimates of all activities.

		# First acquire the global time proportion mutex. This isn't really
		# necessary unless running in parallel mode, but it shouldn't really
		# cost much time if not running in parallel mode. And it's hard to know
		# if parallel mode is enabled from in here.
		global gd_tp_lock
		gd_tp_lock.acquire()

		get_data_read_prop = read_prop
		get_data_fixnan_prop = fixnan_prop
		get_data_t_agg_prop = t_agg_prop
		get_data_bounds_prop = bounds_prop
		if not matching_units:
			get_data_unit_prop = unit_prop

		# Release global time proportion mutex.
		gd_tp_lock.release()

		# Done!
		return data

def read_data(variable: Variable, progress_callback: Callable[[float], None]) \
	-> list[float]:
	"""
	Read all data for a variable from the .nc input file.
	"""
	arr = [float] * variable.size
	n = len(arr)
	for i in range(0, n, CHUNK_SIZE):
		lower = i
		upper = i + CHUNK_SIZE
		arr[lower:upper] = variable[lower:upper]
		progress_callback(i / n)
	return arr

def get_steps_per_year(timestep: int) -> int:
	"""
	Calculate the number of timesteps in a year.
	"""
	step_width_hour = timestep / MINUTES_PER_HOUR
	return int(HOURS_PER_DAY / step_width_hour * DAYS_PER_YEAR)

def get_coord_index(dim: Dimension, var: Variable, value: float) -> int:
	index = index_of(var, value)
	if index == -1:
		if dim.isunlimited():
			index = len(var)
		else:
			edges = numpy.ma.notmasked_edges(var)
			if edges is None:
				index = 0
			else:
				index = edges[1] + 1
				if index >= len(var):
					m = "Unable to insert coordinate %.2f: dimension %s is not long enough"
					raise ValueError(m % (value, dim.name))
		var[index] = value
	return index

def get_coord_indices(nc_out: Dataset, lon: float, lat: float) \
	-> tuple[float, float]:
	"""
	Get the indices of the longitude and latitude, adding them to the file if
	they do not already exist. The return value is a tuple of
	(longitud_index, latitude_index).

	@param nc_out: The output NetCDF file.
	@param lon: Longitude.
	@param lat: Latitude.
	"""
	var_lon = nc_out.variables[DIM_LON]
	dim_lon = get_dim_from_std_name(nc_out, DIM_LON)
	index_lon = get_coord_index(dim_lon, var_lon, lon)

	var_lat = nc_out.variables[DIM_LAT]
	dim_lat = get_dim_from_std_name(nc_out, DIM_LAT)
	index_lat = get_coord_index(dim_lat, var_lat, lat)
	return (index_lon, index_lat)

def trim_to_start_year(in_file: Dataset, timestep: int
	, data: list[float]) -> list[float]:
	"""
	Trim the given data to the start of the next year.

	@param in_file: Input file.
	@param timestep: Output timestep in minutes.
	@param data: The data to be trimmed.
	"""
	start_date = parse_date(in_file.time_coverage_start)
	if (start_date.minute == 30):
		# temporal_aggregation() trims the first value if it lies on the
		# half-hour boundary. We need to do the same here.
		log_debug("Data starts on 30m boundary: 1st value is ignored")
		start_date = start_date + datetime.timedelta(minutes = 30)
	next_year = get_next_year(start_date)

	delta = next_year - start_date
	hours_per_timestep = timestep / MINUTES_PER_HOUR
	seconds_per_timestep = SECONDS_PER_HOUR * hours_per_timestep
	n_trim = int(delta.total_seconds() / seconds_per_timestep)
	m = "Trimming %d values to get to start of %d"
	log_debug(m % (n_trim, next_year.year))
	return data[n_trim:]

def units_match(unit0: str, unit1: str) -> str:
	"""
	Check if the two units are equivalent.
	E.g. m/s and ms^-1 would return true, but not m and kg.
	"""
	if unit0 == unit1:
		return True

	for case in _units_synonyms:
		if unit0 in case and unit1 in case:
			return True
	return False

def remove_nans(data: list[float]
	, progress_cb: Callable[[float], None]) -> list[float]:
	"""
	Replace any NaN in the list with a mean of nearby values.
	"""
	n = len(data)
	N_NEIGHBOUR = 5
	for i in range(n):
		if data[i].mask:
			# Use mean of 5 closest values if value is nan.
			x = neighbouring_mean(data, i, N_NEIGHBOUR)
			m = "Replacing NaN at %d with mean %.2f (n = %d)"
			log_debug(m % (i, x, N_NEIGHBOUR))
			data[i] = x
		if i % PROGRESS_CHUNK_SIZE == 0:
			progress_cb(i / n)
	return data

def aggregate(data: list[float], start: int, stop: int, chunk_size: int \
	, aggregator: Callable[[list[float]], float] \
	, progress_cb: Callable[[float], None]) -> list[float]:
	"""
	Aggregate the values in the specified dataset between the start and
	stop indices, using the specified aggregator function.

	@param data: The data.
	@param start: Aggregate only values after this index.
	@param stop: Aggregate only values before this index.
	@param chunk_size: Aggregate values in chunks of this size.
	@param aggregator: A function to aggregate a chunk of values.
	@param progress_cb: A function which handles progress reporting.
	"""
	n = len(data)
	if start < 0:
		raise ValueError("start cannot be less than 0 (is %d)" % start)
	if stop > n:
		m = "stop cannot be past end of list (is %d and n = %d)"
		raise ValueError(m % (stop, n))
	if start >= stop:
		m = "start cannot be greater than stop (start = %d, stop = %d)"
		raise ValueError(m % (start, stop))

	niter = ((stop - start) // chunk_size)
	result = [0.0] * niter
	for i in range(niter):
		lower = start + i * chunk_size
		upper = min(n, start + (i + 1) * chunk_size)
		result[i] = aggregator(data[lower:upper])
		if i % PROGRESS_CHUNK_SIZE == 0:
			progress_cb(i / niter)
	return result

def temporal_aggregation(in_file: Dataset, data: list[float] \
	, output_timestep: float, aggregator: Callable[[list[float]], float]
	, progress_callback: Callable[[float], None]) -> list[float]:
	"""
	Aggregate the input data to destination timestep.

	@param in_file: Input netcdf file.
	@param data: Data to be aggregated.
	@param forcing: The variable represented by the data.
	@param output_timestep: The desired output timestep.
	@param aggregator: The function used to aggregate values (typically mean or sum).
	@param progress_callback: Function used for progress reporting.
	"""
	# todo: support other timesteps (both source and target)?
	input_timestep = int(in_file.time_step)
	if input_timestep == output_timestep:
		return data
	if output_timestep < input_timestep:
		m = "Invalid output timestep (%d); must be >= input timestep (%d)"
		raise ValueError(m % (output_timestep, input_timestep))
	if output_timestep % input_timestep != 0:
		m = "Invalid timestep; output (%d) must be an integer multiple of input\
timestep (%d)"
		raise ValueError(m % (output_timestep, input_timestep))

	# Record the amount of data at the start of the function. This is
	# written as a diagnostic later.
	n_start = len(data)

	start_minute = get_start_minute(in_file)
	if start_minute == 30 and input_timestep == 30:
		# Don't include the first timestep value if it lies on the half-
		# hour boundary. trim_to_start_year() relies on this assumption,
		# so it will need to be updated if we change this and vice versa.
		m = "Input data starts on 30m boundary. First value will be removed."
		log_debug(m)
		data = data[1:]
	elif start_minute != 0:
		m = "Strange start time: '%s'. Data must start on the hour or half-hour"
		raise ValueError(m % in_file.time_coverage_start)

	# Guaranteed to be an integer division, given the above error checks.
	timestep_ratio = output_timestep // input_timestep

	out: list[float]
	out = []

	if PARALLEL_AGGREGATION:
		# This can be quite slow for large datasets, so I've parallelised it.

		# Get number of CPUs.
		ncpu = multiprocessing.cpu_count()

		# Number of timesteps to be processed by each thread.
		chunk_size = len(data) // ncpu

		# As we're aggregating every 2 timesteps together, we need to ensure
		# that each thread is processing an even number of values.
		if chunk_size % 2 == 1:
			chunk_size += 1

		global thread_progress
		thread_progress = [0.0] * ncpu
		def update_progress(progress: float, tid: int):
			global thread_progress
			global mutex
			with mutex:
				thread_progress[tid] = progress
				progress_callback(numpy.mean(thread_progress))

		class Aggregator(threading.Thread):
			def __init__(self, data: list[float], start: int, stop: int
				, chunk_size: int, tid: int
				, aggregation_func: Callable[[list[float]], float]):
				threading.Thread.__init__(self)
				self.data = data
				self.istart = start
				self.stop = stop
				self.chunk_size = chunk_size
				self.thread_id = tid
				self.aggregation_func = aggregation_func
			def run(self):
				self.result = aggregate(self.data, self.istart, self.stop
				, self.chunk_size, self.aggregation_func
				, lambda p: update_progress(p, self.thread_id))

		# Create an array to hold the thread objects.
		threads = []
		for i in range(ncpu):
			# Calculate start and stop indices for this thread. Each thread
			# is operating on the same input data array, but they operate
			# on different parts (slices) of the array.
			start = i * chunk_size
			stop = len(data) if i == ncpu - 1 else start + chunk_size

			# Create the thread object.
			t = Aggregator(data, start, stop, timestep_ratio, i, aggregator)

			# Start the thread and store the object reference for later.
			t.start()
			threads.append(t)

		# Now we wait for all of the threads to finish, and retrieve their
		# results.
		for t in threads:
			t.join()
			out.extend(t.result)
	else:
		out = aggregate(data, 0, len(data), timestep_ratio, aggregator
		, progress_callback)

	n_end = len(out)
	m = "Temporal aggregation reduced array length by %d values"
	log_debug(m % (n_start - n_end))

	return out

def fix_units(data: list[float], current_units: str, desired_units: str, \
	timestep: int, invert: bool, progress_callback: Callable[[float], None]) \
		-> list[float]:
	"""
	Convert data to the units required for the output file.
	This will modify the existing array.

	@param data: Input data.
	@param current_units: The input units.
	@param desired_units: The output units.
	@param timestep: The input timestep length in minutes.
	@param progress_callback: Function for progress reporting.
	"""
	conversion = find_units_conversion(current_units, desired_units)
	n = len(data)
	timestep *= SECONDS_PER_MINUTE
	scalar = -1 if invert else 1
	for i in range(n):
		data[i] = conversion(data[i], timestep) * scalar
		if i % PROGRESS_CHUNK_SIZE == 0:
			progress_callback(i / n)
	return data

def bounds_checks(data: list[float], xmin: float, xmax: float
	, progress_cb: Callable[[float], None]) -> list[float]:
	"""
	Bounds checking. Any values in the list which exceed theses bounds will be
	set to the boundary value.

	@param data: The data to be checked.
	@param xmin: Lower boundary.
	@param xmax: Upper boundary.
	@param progress_cb: Progress reporting function.
	"""
	n = len(data)
	for i in range(n):
		if data[i] < xmin:
			m = "Value %.2f in row %d exceeds lower bound of %.2f"
			log_debug(m % (data[i], i, xmin))
			data[i] = xmin
		elif data[i] > xmax:
			m = "Value %.2f in row %d exceeds upper bound of %.2f"
			log_debug(m % (data[i], i, xmax))
			data[i] = xmax
		if i % PROGRESS_CHUNK_SIZE == 0:
			progress_cb(i / n)
	return data

def index_of(xarr: list[float], x: float) -> int:
	"""
	Return the index of x in the given list, or -1 if not found.

	@param xarr: A list of floats.
	@param x: The value for which to search.
	"""
	for i in range(0, len(xarr)):
		if floats_equal(x, xarr[i]):
			return i
	return -1

def index_of_throw(needle, haystack, get_error_msg: Callable[[None], str]) -> int:
	"""
	Find the index of a specific element in the list, or throw if not found.

	@param needle: The item to search for.
	@param haystack: The list to be searched.
	@param get_error_msg: A function which returns an error message if needle is not found.
	"""
	for i in range(len(haystack)):
		if haystack[i] == needle:
			return i
	raise ValueError(get_error_msg())

def is_start_of_year(date: datetime.datetime) -> bool:
	"""
	Check if a given datetime object represents the start of a year.
	"""
	return date.second == 0 \
	and date.minute == 0 \
	and date.hour == 0 \
	and date.day == 1 \
	and date.month == 1

def get_next_year(start_date: datetime.datetime) -> datetime.datetime:
	"""
	Get a date representing the first valid date time in the first year
	on or after the given date. So if start_date lies on the exact start
	of a year, start_date will be returned. Otherwise the first day of
	the next year will be returned.
	"""
	if is_start_of_year(start_date):
		m = "get_next_year(): Start date %s lies on start of year"
		log_debug(m % start_date)
		return start_date
	return datetime.datetime(start_date.year + 1, 1, 1, 0, 0, 0)

def find_units_conversion(current_units: str, desired_units: str) \
	-> Callable[[float], float]:
	"""
	Find a conversion between two different units. Throw if not found.
	The return value is a function which takes and returns a float.
	"""
	# units_conversions is a dict mapping unit combos to a conversion.
	# units_conversions: dict[tuple[str, str], Callable[[float],float]]
	combination = (current_units, desired_units)
	if combination in units_conversions:
		return units_conversions[combination]
	for units_type in _units_synonyms:
		if current_units in units_type:
			for synonym in units_type:
				combination = (synonym, desired_units)
				if (combination in units_conversions):
					return units_conversions[combination]

	m = "No unit conversion exists from '%s' to '%s'"
	raise ValueError(m % (current_units, desired_units))

def units_match(x: str, y: str) -> bool:
	"""
	Check if the two units are identical or equivalent.
	"""
	if x == y:
		return True

	for units_type in _units_synonyms:
		if x in units_type and y in units_type:
			return True

	return False

def find_units_conversion_opt(current: str, desired: str) -> Callable[[float], float]:
	"""
	Find a conversion between two different units, if one is required. This will
	be returned in the form of a function which may be called, taking two
	arguments:

	1. The data (n-dimensional list of floats) to be converted
	2. The timestep width in seconds

	If no units conversion is required, this will return a function which
	returns the original data.

	An exception will be raised if the input units cannot be converted to the
	specified output units.
	"""
	if units_match(current, desired):
		return lambda x, _: x
	return find_units_conversion(current, desired)

def get_start_minute(in_file: Dataset):
	"""
	Determine the minute at which the dataset starts (0-59).
	This is done by checking the time_coverage_start attribute of the dataset.
	"""
	start_time = in_file.time_coverage_start
	# yyyy-MM-dd hh:mm:ss
	pattern = r"\d{4}-\d{1,2}-\d{1,2} \d{1,2}:(\d{1,2}):\d{1,2}"
	matches = re.findall(pattern, start_time)
	if len(matches) < 1:
		m = "Unable to parse start time; expected yyyy-MM-dd hh:mm:ss format, but was '%s'"
		raise ValueError(m % in_file.time_coverage_start)
	return int(matches[0])

def write_common_metadata(nc: Dataset, timestep: int):
	"""
	Write standard metadata to the specified output file.

	@param nc: Output nc file.
	@param timestep: Output timestep in minutes.
	"""
	setattr(nc, "processor_script_version", VERSION)
	setattr(nc, "time_step",str(timestep)) # in minutes.

def get_site_name_from_filename(file: str) -> str:
	"""
	Get the abbreviated site name given an input file name, without opening the
	.nc file. E.g.

	"AdelaideRiver_L6_20071017_20090524.nc" -> "AdelaideRiver".

	@param file: input file name.
	"""
	return os.path.basename(file).split("_")[0]

def get_site_name(nc_file: str) -> str:
	"""
	Get the site name attribute in the specified .nc file.

	@param nc_file: netcdf file containing met forcing data for the site.
	"""
	with open_netcdf(nc_file) as nc:
		return nc.site_name

def create_var_if_not_exists(nc: Dataset, name: str, format: str
	, dims: tuple[str], compression_level: int = 0, compression_type: str = None
	, chunksizes: tuple[int] = None):
	"""
	Create a variable in the NetCDF file if it does not already exist.

	@param nc: The NetCDF file.
	@param name: Name of the variable.
	@param format: Format of the variable.
	@param dims: Dimensions of the variable.
	@param compression_level: Compression level [1-9], 9 = slowest/smallest.
	@param compression_type: Compression algorithm to be used (or None).
	"""
	compression = None if compression_level == 0 else compression_type
	if not name in nc.variables:
		log_diagnostic(f"Creating variable {name} with format {format}, dimensions {dims}, compression {compression} L{compression_level}, and chunk sizes {chunksizes}")
		nc.createVariable(name, format, dims, compression = compression
		, complevel = compression_level, chunksizes = chunksizes)
		log_debug(f"Variable {name} created successfully")

		log_debug(f"Setting attribute {_ATTR_MISSING_VAL} on variable {name}")
		var = nc.variables[name]
		setattr(var, _ATTR_MISSING_VAL, default_fillvals[format])
		log_debug(f"Attribute {_ATTR_MISSING_VAL} created successfully on variable {name}")
	else:
		log_diagnostic(f"Variable {name} will not be created because it already exists")

def create_dim_if_not_exists(nc: Dataset, name: str, size: int = 0):
	"""
	Create an unlimited dimension in the NetCDF file if it does not already
	exist.

	@param nc: The NetCDF file.
	@param name: Name of the dimension.
	"""
	if not name in nc.dimensions:
		log_diagnostic(f"Creating dimension {name} with size {size}")
		nc.createDimension(name, size)
		log_debug(f"Dimension {name} created successfully")
	else:
		log_diagnostic(f"Dimension {name} not created because it already exists")

def open_netcdf(file: str, write: bool = False, format: str = NC_FORMAT) -> Dataset:
	"""
	Open a netcdf file.

	@param file: The file to be opened.
	@param mode: The opening mode. 'r' for readonly (default) or 'w' for write.
	"""
	mode = "r+" if write else "r"
	log_debug(f"Opening netcdf file '{file}' in mode '{mode}'...")
	dataset = Dataset(file, mode, format = format)
	log_diagnostic(f"Successfully opened netcdf file '{file}' in mode '{mode}'...")
	return dataset

def find_dimension(nc: Dataset, predicate: Callable[[Dimension], bool]
	, err_txt: str = "Failed to find dimension") -> Dimension:
	"""
	Return the first dimension in the NetCDF file which satisfies the given
	condition. If no matching dimension is found, an exception will be raised.

	@param nc: The NetCDF file to be searched.
	@param predicate: A function used to filter the dimension. Should return
					  true if the dimension is a match.
	"""
	for name in nc.dimensions:
		dim = nc.dimensions[name]
		if predicate(dim):
			return dim
	raise ValueError(err_txt)

def get_dim_with_attr(nc: Dataset, attr_name: str, attr_value: str
					  , err_txt: str) -> Dimension:
	"""
	Return the first dimension in a NetCDF whose corresponding variable has an
	attribute matching the specified value.

	@param nc: The input NetCDF file.
	@param attr_name: Name of the attribute to check in the variable.
	@param attr_value: Desired value of the attribute.
	@param err_txt: Text of the error message.s
	"""
	return find_dimension(nc, lambda dim: \
		getattr(nc.variables[dim.name], attr_name) == attr_value, err_txt)


def get_dim_from_std_name(nc, standard_name) -> Dimension:
	"""
	Return the dimension with the specified standard name. Throw if not found.

	@param nc: The netcdf file.
	"""
	return get_dim_with_attr(nc, ATTR_STD_NAME, standard_name,
		f"Input file contains no dimension with standard name {standard_name}")

def count_gridpoints(file: str) -> int:
	"""
	Count the number of gridpoints in the specified netcdf file.
	"""
	with open_netcdf(file) as nc:
		lat = get_dim_from_std_name(nc, DIM_LAT)
		lon = get_dim_from_std_name(nc, DIM_LON)
		return lat.size * lon.size

def find_dim(nc: Dataset, dim_name: str) -> Dimension:
	"""
	Find a dimension with the specified name or standard name.

	@param nc: The opened input NetCDF file.
	@param dim_name: The dimension name to search for. Will return any dimension
	whose name or standard_name matches this value.
	"""
	if dim_name in nc.dimensions:
		return nc.dimensions[dim_name]
	return get_dim_from_std_name(nc, dim_name)

def get_dimension_variable(nc: Dataset, dim: Dimension) -> Variable:
	"""
	Find the variable in a NetCDF file which holds the variable behind the
	specified dimension.

	@param nc: The input NetCDF file.
	@param dim: The dimension.
	"""
	# is this correct? Can we have variable lat behind latitude??
	return nc.variables[dim.name]

def copy_attr(var_in: Variable, var_out: Variable, attr: str):
	"""
	Copy the specified attribute value from one variable to another.

	@param var_in: Input variable.
	@param var_out: Output variable.
	@param attr: Name of the attribute.
	"""
	if hasattr(var_in, attr):
		value = getattr(var_in, attr)
		setattr(var_out, attr, value)

def copy_attributes(nc_in, nc_out):
	"""
	Copy all attributes from the input object to the output object. These may
	be netcdf datasets (ie files), variables, dimensions, etc.
	"""
	for attr in nc_in.ncattrs():
		if not attr[0] == "_":
			setattr(nc_out, attr, getattr(nc_in, attr))

def get_compression_type(filters: dict) -> str:
	"""
	Determine the compression algorithm used for a variable, given that
	variable's filters.

	@param filters: The variable filters, obtained by variable.filters().
	"""
	types = ["zlib", "szip", "zstd", "bzip2", "blosc", "fletcher32"]
	for compression_type in types:
		if filters[compression_type]:
			return compression_type
	return ""

def get_dimension_indices(nc_in: Dataset, nc_out: Dataset, name: str
	, low: int, high: int) -> list[int]:
	"""
	Get the dimension indices in the output file which correspond to the
	dimension indices in the input file within the specified range.

	@param nc_in: Input .nc file.
	@param nc_out: Output .nc file.
	@param name: Dimension name.
	@param low: Start of the desired dimension index range.
	@param high: End of the desired dimension index range.
	"""

	##### fixme: hi drew: dimension can have different name to the variable (e.g. time_bnds has dimension bnds. glhf! :)
	values = nc_in.variables[name][low:high]
	dim = nc_out.dimensions[name]
	var = nc_out.variables[name]
	return [get_coord_index(dim, var, value) for value in values]

def copy_3d(nc_in: Dataset, nc_out: Dataset, name: str, min_chunk_size: int
		, pcb: Callable[[float], None]):
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

	for i in range(nd):
		if dim_names_in[i] != dim_names_out[i]:
			m = "Inconsistent dimension order in variable '%s': dim[%d] in input file is '%s' but in output file it is '%s'"
			raise ValueError(m % (name, i, dim_names_in[i], dim_names_out[i]))

	dims = [nc_in.dimensions[name] for name in dim_names_in]

	chunk_sizes = var_in.chunking()
	chunk_sizes = [max(min_chunk_size, c) for c in chunk_sizes]
	shape = var_in.shape
	niter = [math.ceil(s / c) for (s, c) in zip(shape, chunk_sizes)]

	it_max = niter[0] * niter[1] * niter[2]
	it = 0
	samei = False
	samej = False
	samek = False
	for i in range(niter[0]):
		ilow = i * chunk_sizes[0]
		ihigh = min(shape[0], ilow + chunk_sizes[0])
		ir = range(ilow, ihigh) if samei else get_dimension_indices(nc_in, nc_out, dims[0].name, ilow, ihigh)
		if not samei and ir[0] == ilow and ir[len(ir) - 1] == ihigh - 1:
			samei = True
		for j in range(niter[1]):
			jlow = j * chunk_sizes[1]
			jhigh = min(shape[1], jlow + chunk_sizes[1])
			jr = range(jlow, jhigh) if samej else get_dimension_indices(nc_in, nc_out, dims[1].name, jlow, jhigh)
			if not samej and jr[0] == jlow and jr[len(jr) - 1] == jhigh - 1:
				samej = True
			for k in range(niter[2]):
				klow = k * chunk_sizes[2]
				khigh = min(shape[2], klow + chunk_sizes[2])
				kr = range(klow, khigh) if samek else get_dimension_indices(nc_in, nc_out, dims[2].name, klow, khigh)
				if not samek and kr[0] == klow and kr[len(kr) - 1] == khigh - 1:
					samek = True
				var_out[ir, jr, kr] = var_in[ilow:ihigh, jlow:jhigh, klow:khigh]
				it += 1
				pcb(it / it_max)

def copy_2d(nc_in: Dataset, nc_out: Dataset, name: str, min_chunk_size: int
		, pcb: Callable[[float], None]):
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

	for i in range(nd):
		if dim_names_in[i] != dim_names_out[i]:
			m = "Inconsistent dimension order in variable '%s': dim[%d] in input file is '%s' but in output file it is '%s'"
			raise ValueError(m % (name, i, dim_names_in[i], dim_names_out[i]))

	dims = [nc_in.dimensions[name] for name in dim_names_in]

	chunk_sizes = var_in.chunking()
	chunk_sizes = [max(min_chunk_size, c) for c in chunk_sizes]
	shape = var_in.shape
	niter = [math.ceil(s / c) for (s, c) in zip(shape, chunk_sizes)]

	it_max = niter[0] * niter[1]
	it = 0
	samei = False
	samej = False
	for i in range(niter[0]):
		ilow = i * chunk_sizes[0]
		ihigh = min(shape[0], ilow + chunk_sizes[0])
		ir = range(ilow, ihigh) if samei else get_dimension_indices(nc_in, nc_out, dims[0].name, ilow, ihigh)
		if not samei and ir[0] == ilow and ir[len(ir) - 1] == ihigh - 1:
			samei = True
		for j in range(niter[1]):
			jlow = j * chunk_sizes[1]
			jhigh = min(shape[1], jlow + chunk_sizes[1])
			jr = range(jlow, jhigh) if samej or not dims[1].name in nc_in.variables else get_dimension_indices(nc_in, nc_out, dims[1].name, jlow, jhigh)
			if not samej and jr[0] == jlow and jr[len(jr) - 1] == jhigh - 1:
				samej = True
			var_out[ir, jr] = var_in[ilow:ihigh, jlow:jhigh]
			it += 1
			pcb(it / it_max)


def copy_1d(nc_in: Dataset, nc_out: Dataset, name: str, min_chunk_size: int
		, pcb: Callable[[float], None]):
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
	chunk_size = var_in.chunking()[0]
	if var_in.chunking() == "contiguous":
		if min_chunk_size == 0:
			chunk_size = 64
		else:
			chunk_size = min_chunk_size
	else:
		chunk_size = max(min_chunk_size, chunk_size)
	n = len(var_in)
	for i in range(0, n, chunk_size):
		upper = min(n, i + chunk_size)
		var_out[i:upper] = var_in[i:upper]
		pcb(upper / n)

def copy_transpose_3d(nc_in: Dataset, nc_out: Dataset, name: str
	, min_chunk_size: int, pcb: Callable[[float], None]):
	"""
	Copy a 3d variable from the input file to the output file, reordering the
	dimensions as specified. If the dimensionality is identical in the two
	files, it would be faster to call copy_variable().

	@param nc_in: Input netcdf file.
	@param nc_out: Output netcdf file.
	@param name: Name of the variable to be copied.
	@param opts: Reshape options.
	@param pcb: Progress callback function.
	"""
	var_in = nc_in.variables[name]
	var_out = nc_out.variables[name]

	dims_in = var_in.dimensions
	dims_out = var_out.dimensions

	# This is currently also checked by the caller.
	if len(dims_in) != len(dims_out):
		raise ValueError(f"Inconsistent number of dimensions for variable {name}")

	# Indices of input dimensions in the order of the output file.
	dimension_indices = [index_of_throw(dim, dims_in, lambda: f"Dimension {dim} not found in input file") for dim in dims_out]

	chunk_sizes = var_in.chunking()
	chunk_sizes = [max(min_chunk_size, c) for c in chunk_sizes]

	# Length of each dimension.
	shape = var_in.shape

	# Number of iterations required for each dimension with given chunk sizes.
	niter = [math.ceil(s / c) for (s, c) in zip(shape, chunk_sizes)]

	# Total number of iterations required.
	it_max = niter[0] * niter[1] * niter[2]

	# Index of the current iteration.
	it = 0

	for i in range(niter[0]):
		# Index of the first element to be read on the i-th dimension.
		ilow = i * chunk_sizes[0]

		# Index of the last element to be read on the i-th dimension.
		ihigh = min(shape[0], ilow + chunk_sizes[0])

		# Index range of the i-th dimension of the input file.
		ir = range(ilow, ihigh)

		for j in range(niter[1]):
			# Index of the first element to be read on the j-th dimension.
			jlow = j * chunk_sizes[1]

			# Index of the last element to be read on the j-th dimension.
			jhigh = min(shape[1], jlow + chunk_sizes[1])

			# Index range of the j-th dimension of the input file.
			jr = range(jlow, jhigh)

			for k in range(niter[2]):
				# Index of the first element to be read on the k-th dimension.
				klow = k * chunk_sizes[2]

				# Index of the last element to be read on the k-th dimension.
				khigh = min(shape[2], klow + chunk_sizes[2])

				# Index range of the k-th dimension of the input file.
				kr = range(klow, khigh)

				# Read a chunk of data from the input file.
				chunk = var_in[ilow:ihigh, jlow:jhigh, klow:khigh]

				# Transpose this chunk into the shape in the output file.
				chunk = numpy.transpose(chunk, dimension_indices)

				# Get the index ranges for the output file.
				in_ranges = [ir, jr, kr]
				out_ranges = [in_ranges[i] for i in dimension_indices]

				# Write data to the output file.
				var_out[out_ranges[0], out_ranges[1], out_ranges[2]] = chunk

				# Progress reporting.
				it += 1
				pcb(it / it_max)

def copy_transpose_2d(nc_in: Dataset, nc_out: Dataset, name: str
	, min_chunk_size: int, pcb: Callable[[float], None]):
	"""
	Copy a 2d variable from the input file to the output file, reordering the
	dimensions as specified. If the dimensionality is identical in the two
	files, it would be faster to call copy_variable().

	@param nc_in: Input netcdf file.
	@param nc_out: Output netcdf file.
	@param name: Name of the variable to be copied.
	@param opts: Reshape options.
	@param pcb: Progress callback function.
	"""
	raise ValueError("TBI")

def copy_variable_transpose(nc_in: Dataset, nc_out: Dataset, name: str
	, min_chunk_size: int, pcb: Callable[[float], None]):
	"""
	Copy a variable from the input file to the output file, reordering the
	dimensions as specified. If the dimensionality is identical in the two
	files, it would be faster to call copy_variable().

	@param nc_in: Input netcdf file.
	@param nc_out: Output netcdf file.
	@param name: Name of the variable to be copied.
	@param opts: Reshape options.
	@param pcb: Progress callback function.
	"""
	ndim = len(nc_in.variables[name].dimensions)
	if ndim > 3:
		raise ValueError(f"Unable to copy variable {name}: >3 dimensions not implemented (variable {name} has {ndim} dimensions)")
	elif ndim == 3:
		copy_transpose_3d(nc_in, nc_out, name, min_chunk_size, pcb)
	elif ndim == 2:
		copy_transpose_2d(nc_in, nc_out, name, min_chunk_size, pcb)
	elif ndim == 1:
		raise ValueError(f"Transposing 1-dimensional variable {name} doesn't make sense.")

def consistent_dimensionality(var_in: Variable, var_out: Variable) -> bool:
	"""
	Check if the two variables have consistent dimensionality (ie the same
	dimensions in the same order).

	This will throw if the two variables have a different number of dimensions.

	@param var_in: The first variable.
	@param var_out: The second variable.
	"""
	name = var_in.name

	dims_in = var_in.dimensions
	dims_out = var_out.dimensions

	if len(dims_in) != len(dims_out):
		raise ValueError(f"Unable to copy variable {name}: input file has {len(dims_in)} dimensions but output file has {len(dims_out)} dimensions")

	for i in range(len(dims_in)):
		if dims_in[i] != dims_out[i]:
			return False

	return True

def _copy_variable(nc_in: Dataset, nc_out: Dataset, name: str
	, min_chunk_size: int, pcb: Callable[[float], None]):

	# Call the appropriate function based on number of dimensions.
	ndim = len(nc_in.variables[name].dimensions)

	if ndim > 3:
		raise ValueError(f"Unable to copy variable {name}: >3 dimensions not implemented (variable {name} has {ndim} dimensions)")
	elif ndim == 3:
		copy_3d(nc_in, nc_out, name, min_chunk_size, pcb)
	elif ndim == 2:
		copy_2d(nc_in, nc_out, name, min_chunk_size, pcb)
	elif ndim == 1:
		copy_1d(nc_in, nc_out, name, min_chunk_size, pcb)

def copy_variable(nc_in: Dataset, nc_out: Dataset, name: str
	, min_chunk_size: int, pcb: Callable[[float], None]):
	"""
	Copy the contents of the specified variable in the input file into the
	output file.

	This will transpose data as necessary if the output file has different
	dimensionality for this variable.

	@param nc_in: The input .nc file.
	@param nc_out: The output .nc file.
	@param anme: The name of the variable to be copied.
	@param min_chunk_size: Minimum chunk size used when copying data.
	@param pcb: Progress callback function.
	"""
	var_in = nc_in.variables[name]
	var_out = nc_out.variables[name]

	if consistent_dimensionality(var_in, var_out):
		_copy_variable(nc_in, nc_out, name, min_chunk_size, pcb)
	else:
		copy_variable_transpose(nc_in, nc_out, name, min_chunk_size, pcb)

def check_ndims(var, expected):
	actual = len(var.dimensions)
	if actual != expected:
		raise ValueError(f"Unable to copy variable {var.name} as {expected}-dimensional variable: variable has {actual} dimensions")

def _append_1d(nc_in: Dataset, nc_out: Dataset, name: str, min_chunk_size: int, pcb: Callable[[float], None]):
	"""
	Append (along the time axis) the contents of the specified 1-dimensional
	variable in the input file to the variable in the output file.

	This assumes that the input data is temporally adjacent to the data in the
	output file (if indeed there is any data in the output file).

	@param nc_in: Input netcdf file.
	@param nc_out: Output netcdf file.
	@param name: Name of the variable to be copied.
	@param min_chunk_size: Minimum chunk size used when copying data.
	@param pcb: Progress reporting function.
	"""
	var_in = nc_in.variables[name]
	var_out = nc_out.variables[name]

	check_ndims(var_in, 1)
	check_ndims(var_out, 1)

	chunk_size = var_in.chunking()[0]
	chunk_size = max(chunk_size, min_chunk_size)

	n = len(var_in)
	nexist = len(var_out)
	for i in range(0, n, chunk_size):
		upper = min(n, i + chunk_size)
		var_out[i + nexist:upper + nexist] = var_in[i:upper]
		pcb(upper / n)

def _append_2d(nc_in: Dataset, nc_out: Dataset, name:str, min_chunk_size: int, pcb: Callable[[float], None]):
	"""
	Append (along the time axis) the contents of the specified 2-dimensional
	variable in the input file to the variable in the output file.

	This assumes that the input data is temporally adjacent to the data in the
	output file (if indeed there is any data in the output file).

	@param nc_in: Input netcdf file.
	@param nc_out: Output netcdf file.
	@param name: Name of the variable to be copied.
	@param min_chunk_size: Minimum chunk size used when copying data.
	@param pcb: Progress reporting function.
	"""
	var_in = nc_in.variables[name]
	var_out = nc_out.variables[name]

	check_ndims(var_in, 2)
	check_ndims(var_out, 2)

	dim_names_in = var_in.dimensions
	dim_names_out = var_out.dimensions
	nd = len(dim_names_in)

	if nd != len(dim_names_out):
		m = "Inconsistent dimensionality in variable '%s': input file has %d dimensions (%s), output file has %d dimensions (%s)"
		raise ValueError(m % (name, nd, ", ".join(dim_names_in), len(dim_names_out), ", ".join(dim_names_out)))

	for i in range(nd):
		if dim_names_in[i] != dim_names_out[i]:
			m = "Inconsistent dimension order in variable '%s': dim[%d] in input file is '%s' but in output file it is '%s'"
			raise ValueError(m % (name, i, dim_names_in[i], dim_names_out[i]))

	dims = [nc_in.dimensions[name] for name in dim_names_in]

	# The chunk size of each dimension in the input file.
	chunk_sizes = var_in.chunking()
	chunk_sizes = [max(min_chunk_size, cs) for cs in chunk_sizes]

	# The length of each dimension in the input file.
	shape = var_in.shape

	# The number of iterations required for each dimension, given the above
	# chunk sizes.
	niter = [math.ceil(s / c) for (s, c) in zip(shape, chunk_sizes)]

	# Total number of iterations required.
	it_max = niter[0] * niter[1]

	# The current iteration.
	it = 0

	# The offset into the output file's time dimension occupied by the data
	# being copied from the current input file.
	time_offset = nc_out.dimensions[DIM_TIME].size - nc_in.dimensions[DIM_TIME].size

	# These tell us whether we can use the same indices in the input and output
	# files for a particular dimension. This will be true for spatial dimensions
	# (which are assumed to be identical between files) and it will be false for
	# temporal dimensions (which is the dimension being appended).
	offseti = time_offset if dims[0].name == DIM_TIME else 0
	offsetj = time_offset if dims[1].name == DIM_TIME else 0

	for i in range(niter[0]):
		# ilow/ihigh are the indices used for reading data.
		# ir are the indices used for writing data.

		# Start index for this iteration on the i-th dimension.
		ilow = i * chunk_sizes[0]

		# End index for this iteration on the i-th dimension.
		ihigh = min(shape[0], ilow + chunk_sizes[0])

		# Index range for the i-th dimension.
		ir = range(ilow + offseti, ihigh + offseti)

		for j in range(niter[1]):
			# jlow/jhigh are the indices used for reading data.
			# jr are the indices used for writing data.

			# Start index for this iteration of the j-th dimension.
			jlow = j * chunk_sizes[1]

			#  End index for this iteration on the j-th dimension.
			jhigh = min(shape[1], jlow + chunk_sizes[1])

			# Index range for the j-th dimension.
			jr = range(jlow + offsetj, jhigh + offsetj)

			# Copy the data.
			var_out[ir, jr] = var_in[ilow:ihigh, jlow:jhigh]

			# Progress tracking/reporting.
			it += 1
			pcb(it / it_max)

def get_var_from_std_name(nc: Dataset, name: str) -> Variable:
	"""
	Find a variable in the input file with the specified standard name. Raise
	an exception if no matching variable is found.

	@param nc: The input netcdf file.
	@param name: The standard name for which we search.
	"""
	for var in nc.variables:
		var = nc.variables[var]
		if hasattr(var, ATTR_STD_NAME):
			if getattr(var, ATTR_STD_NAME) == name:
				return var
	raise ValueError(f"No variable found with standard_name '{name}'")

def get_timestep(nc: Dataset) -> int:
	"""
	Get the timestep of the netcdf file, in seconds.
	"""
	var_time = get_var_from_std_name(nc, STD_TIME)
	units = getattr(var_time, ATTR_UNITS)
	calendar = getattr(var_time, ATTR_CALENDAR)

	chunk_size = var_time.chunking()
	chunk_size = 64 if chunk_size == "contiguous" else chunk_size[0]
	n = len(var_time)
	for i in range(0, n, chunk_size):
		upper = min(n, i + chunk_size)
		raw = var_time[i:upper]
		times = num2date(raw, units, calendar, only_use_cftime_datetimes = False, only_use_python_datetimes = True)
		first_date = times[0].date()
		for j in range(len(times)):
			if times[j].date() != first_date:
				timestep = SECONDS_PER_DAY // j
				log_diagnostic(f"Input timestep is {timestep} seconds ({j} steps per day)")
				return timestep

	raise ValueError(f"Unable to determine timestep width: all data appears to be on the same day")

def get_nexist(var: Variable) -> int:
	"""
	Get the number of existing elements in the variable.

	@param var: Variable in a NetCDF file.
	"""
	# return len(var_out)
	edges = numpy.ma.notmasked_edges(var)
	return 0 if edges is None else edges[1] + 1 #numpy.ma.flatnotmasked_edges(var_out)[0]

def _append_3d(nc_in: Dataset, nc_out: Dataset, name: str, min_chunk_size: int
			   , pcb: Callable[[float], None]):
	"""
	Append (along the time axis) the contents of the specified 3-dimensional
	variable in the input file to the variable in the output file.

	This assumes that the input data is temporally adjacent to the data in the
	output file (if indeed there is any data in the output file).

	@param nc_in: Input netcdf file.
	@param nc_out: Output netcdf file.
	@param name: Name of the variable to be copied.
	@param min_chunk_size: Minimum chunk size used when copying data.
	@param pcb: Progress reporting function.
	"""
	var_in = nc_in.variables[name]
	var_out = nc_out.variables[name]

	check_ndims(var_in, 3)
	check_ndims(var_out, 3)

	dim_names_in = var_in.dimensions
	dim_names_out = var_out.dimensions
	nd = len(dim_names_in)

	if nd != len(dim_names_out):
		m = "Inconsistent dimensionality in variable '%s': input file has %d dimensions (%s), output file has %d dimensions (%s)"
		raise ValueError(m % (name, nd, ", ".join(dim_names_in), len(dim_names_out), ", ".join(dim_names_out)))

	for i in range(nd):
		if dim_names_in[i] != dim_names_out[i]:
			m = "Inconsistent dimension order in variable '%s': dim[%d] in input file is '%s' but in output file it is '%s'"
			raise ValueError(m % (name, i, dim_names_in[i], dim_names_out[i]))

	dims = [nc_in.dimensions[name] for name in dim_names_in]

	# The chunk size of each dimension in the input file.
	chunk_sizes = var_in.chunking()
	chunk_sizes = [max(min_chunk_size, cs) for cs in chunk_sizes]

	# The length of each dimension in the input file.
	shape = var_in.shape

	# The number of iterations required for each dimension, given the above
	# chunk sizes.
	niter = [math.ceil(s / c) for (s, c) in zip(shape, chunk_sizes)]

	# Total number of iterations required.
	it_max = niter[0] * niter[1] * niter[2]

	# The current iteration.
	it = 0

	var_time = nc_out.variables[VAR_TIME]
	time_offset = get_nexist(var_time) - nc_in.dimensions[DIM_TIME].size
	log_diagnostic(f"time_offset = {time_offset}")

	# These tell us whether we can use the same indices in the input and output
	# files for a particular dimension. This will be true for spatial dimensions
	# (which are assumed to be identical between files) and it will be false for
	# temporal dimensions (which is the dimension being appended).
	offseti = time_offset if dims[0].name == DIM_TIME else 0
	offsetj = time_offset if dims[1].name == DIM_TIME else 0
	offsetk = time_offset if dims[2].name == DIM_TIME else 0

	units_conversion = None
	if hasattr(var_in, ATTR_UNITS) and hasattr(var_out, ATTR_UNITS):
		in_units = getattr(var_in, ATTR_UNITS)
		out_units = getattr(var_out, ATTR_UNITS)
		units_conversion = find_units_conversion_opt(in_units, out_units)
	else:
		units_conversion = lambda x, _: x
	timestep = get_timestep(nc_in)

	for i in range(niter[0]):
		# ilow/ihigh are the indices used for reading data.
		# ir are the indices used for writing data.

		# Start index for this iteration on the i-th dimension.
		ilow = i * chunk_sizes[0]

		# End index for this iteration on the i-th dimension.
		ihigh = min(shape[0], ilow + chunk_sizes[0])

		# Index range for the i-th dimension.
		ir = range(ilow + offseti, ihigh + offseti)

		for j in range(niter[1]):
			# jlow/jhigh are the indices used for reading data.
			# jr are the indices used for writing data.

			# Start index for this iteration on the j-th dimension.
			jlow = j * chunk_sizes[1]

			# End index for this iteration on the j-th dimension.
			jhigh = min(shape[1], jlow + chunk_sizes[1])

			# Index range for the j-th dimension.
			jr = range(jlow + offsetj, jhigh + offsetj)

			for k in range(niter[2]):
				# klow/khigh are the indices used for reading data.
				# kr are the indices used for writing data.

				# Start index for this iteration on the k-th dimension.
				klow = k * chunk_sizes[2]

				# End index for this iteration on the k-th dimension.
				khigh = min(shape[2], klow + chunk_sizes[2])

				# Index range for the k-th dimension.
				kr = range(klow + offsetk, khigh + offsetk)

				# Copy the data.
				chunk = var_in[ilow:ihigh, jlow:jhigh, klow:khigh]
				var_out[ir, jr, kr] = units_conversion(chunk, timestep)

				# Progress tracking/reporting.
				it += 1
				pcb(it / it_max)

def append_time(nc_in: Dataset, nc_out: Dataset, name: str, min_chunk_size: int
				, pcb: Callable[[float], None]):
	"""
	Append (along the time axis) the contents of the specified variable in the
	input file to the variable in the output file.

	This assumes that the input data is temporally adjacent to the data in the
	output file (if indeed there is any data in the output file).

	@param nc_in: Input netcdf file.
	@param nc_out: Output netcdf file.
	@param name: Name of the variable to be copied.
	@param min_chunk_size: Minimum chunk size used when copying data.
	@param pcb: Progress reporting function.
	"""
	ndim = len(nc_in.variables[name].dimensions)
	if ndim > 3:
		raise ValueError(f"{ndim}-dimensional variable encountered: {name}. >3 dimensions is not supported")
	elif ndim == 3:
		# Typical spatial variables.
		_append_3d(nc_in, nc_out, name, min_chunk_size, pcb)
	elif ndim == 2:
		# Boundary variable (time_bnds)
		_append_2d(nc_in, nc_out, name, min_chunk_size, pcb)
	elif ndim == 1:
		# Coordinate variable (e.g. time).
		_append_1d(nc_in, nc_out, name, min_chunk_size, pcb)

def get_nc_datatype(fmt: str) -> str:
	if fmt == "float64":
		return FORMAT_FLOAT
	if fmt == "float32":
		return FORMAT_SINGLE
	return fmt
