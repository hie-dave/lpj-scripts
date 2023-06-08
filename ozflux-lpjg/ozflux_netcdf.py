import datetime, numpy, os, re, threading, time
import multiprocessing, multiprocessing.connection
import traceback, sys
from ozflux_common import *
from ozflux_logging import *
from netCDF4 import Dataset, Variable
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

# Name of the time variable in the input files.
VAR_TIME = "time"

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
	["kg/m2/s", "kg/m^2/s", "kgm^-2s^-1"],
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
	var_lat = nc_out.variables[DIM_LAT]
	index_lon = index_of(var_lon, lon)
	if index_lon == -1:
		if nc_out.dimensions[DIM_LON].isunlimited():
			index_lon = len(var_lon)
		else:
			edges = numpy.ma.notmasked_edges(var_lon)
			if edges is None:
				index_lon = 0
			else:
				index_lon = edges[1] + 1
				if index_lon >= len(var_lon):
					m = "Unable to insert coordinate (%.2f, %.2f): dimension is not long enough"
					raise ValueError(m % (lon, lat))
		var_lon[index_lon] = lon
	index_lat = index_of(var_lat, lat)
	if index_lat == -1:
		if nc_out.dimensions[DIM_LAT].isunlimited():
			index_lat = len(var_lat)
		else:
			edges = numpy.ma.notmasked_edges(var_lat)
			if edges is None:
				index_lat = 0
			else:
				index_lat = edges[1] + 1
				if index_lon >= len(var_lon):
					m = "Unable to insert coordinate (%.2f, %.2f): dimension is not long enough"
					raise ValueError(m % (lon, lat))
		var_lat[index_lat] = lat
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
	with Dataset(nc_file, "r", format=NC_FORMAT) as nc:
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
	log_diagnostic("Compression = %s L%d" % (compression, compression_level))
	if not name in nc.variables:
		nc.createVariable(name, format, dims, compression = compression
		, complevel = compression_level, chunksizes = chunksizes)

def create_dim_if_not_exists(nc: Dataset, name: str, size: int = 0):
	"""
	Create an unlimited dimension in the NetCDF file if it does not already
	exist.

	@param nc: The NetCDF file.
	@param name: Name of the dimension.
	"""
	if not name in nc.dimensions:
		nc.createDimension(name, size)
