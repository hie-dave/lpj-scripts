import datetime, numpy, os, ozflux_common, re, threading, time
import multiprocessing, multiprocessing.connection
import traceback, sys
from ozflux_logging import *
from netCDF4 import Dataset, Variable
from enum import IntEnum
from ozflux_common import Forcing
from sys import argv
from typing import Callable

# Name of the time dimension in the output file.
OUT_DIM_NAME_TIME = "time"

# Name of the longitude dimension in the output file.
OUT_DIM_NAME_LON = "longitude"

# Name of the latitude dimension in the output file.
OUT_DIM_NAME_LAT = "latitude"

# g to kg. This is not really used except as an example.
G_TO_KG = 1e-3

# KPa to Pa
KPA_TO_PA = 1e3

# °C to K
DEG_C_TO_K = 273.15

# Number of seconds per minute.
SECONDS_PER_MINUTE = 60

# Number of minutes per hour.
MINUTES_PER_HOUR = 60

# Number of seconds per hour.
SECONDS_PER_HOUR = SECONDS_PER_MINUTE * MINUTES_PER_HOUR

# Number of hours per day.
HOURS_PER_DAY = 24

# Number of days per year.
DAYS_PER_YEAR = 365 # ha

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
units_synonyms = [
	["W/m2", "W/m^2", "Wm^-2"],
	["kg/m2/s", "kg/m^2/s", "kgm^-2s^-1"],
	["K", "k"],
	["m/s", "ms^-1"],
	["Pa", "pa"],
	["kg/kg", "", "mm/mm", "m/m"], # debatable
	["ppm", "umol/mol"],
	["degC", "°C", "degrees C"]
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
	("g", "kg"): lambda x, _: x * G_TO_KG, # (as an example)
	("mm", "kg/m2/s"): lambda x, t: x / t, # Assuming mm means mm per timestep
	("kg/m2/s", "mm"): lambda x, t: t * x,
	("degC", "K"): lambda x, _: x + DEG_C_TO_K,
	("K", "degC"): lambda x, _: x - DEG_C_TO_K,
	("kPa", "Pa"): lambda x, _: x * KPA_TO_PA
}

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
	, in_name: str
	, out_name: str \
	, output_timestep: int \
	, in_units: str \
	, out_units: str \
	, aggregator: Callable[[list[float]], float] \
	, lower_bound: float \
	, upper_bound: float \
	, parallel: bool \
	, progress_cb: Callable[[float], None]) -> list[float]:
		"""
		Get all data from the input file and convert it to a format
		suitable for the output file.

		@param in_name: Name of the variable in the input file.
		@param out_name: Name of the variable in the output file.
		@param output_timestep: Output timestep width in minutes.
		"""
		# Some variables are not translated directly into the output
		# file and can be ignored (ie filled with zeroes).
		if in_name == "zero":
			log_diagnostic("Using zeroes for %s" % out_name)
			input_timestep = int(in_file.time_step)
			timestep_ratio = output_timestep // input_timestep
			n = in_file.variables["time"].size // timestep_ratio
			data = zeroes(n)
			return trim_to_start_year(in_file, output_timestep, data)
		if not in_name in in_file.variables:
			m = "Variable %s does not exist in input file"
			raise ValueError(m % in_name)

		var_id = in_file.variables[in_name]
		# in_units = get_units(var_id)
		# out_units = forcing_units[forcing]
		matching_units = units_match(in_units, out_units)

		global get_data_read_prop, get_data_fixnan_prop, get_data_t_agg_prop
		global get_data_unit_prop, get_data_bounds_prop

		# About 1/3 of time to fix units (if units need fixing).
		units_time_proportion = 0 if matching_units else get_data_unit_prop

		prop_tot = get_data_read_prop + get_data_fixnan_prop + \
			get_data_t_agg_prop + units_time_proportion + get_data_bounds_prop

		# About 2/3 of remaining time to read the data from input file.
		read_time_proportion = get_data_read_prop / prop_tot

		# 5% of time to remove NaNs. Need to check this.
		fixnan_time_proportion = get_data_fixnan_prop / prop_tot

		# Rest of time is spent aggregating over the timestep.
		aggregation_time_proportion = get_data_t_agg_prop / prop_tot

		bounds_time_prop = get_data_bounds_prop / prop_tot

		# Read data from netcdf file.
		read_start = time.time()
		data = read_data(var_id, lambda p:progress_cb(read_time_proportion * p))
		read_tot = time.time() - read_start
		log_debug("Successfully read %d values from netcdf" % len(data))

		# Replace NaN with a mean of nearby values.
		log_diagnostic("Removing NaNs")
		step_start = read_time_proportion
		step_size = fixnan_time_proportion
		fixnan_start = time.time()
		data = remove_nans(data
		, lambda p: progress_cb(step_start + step_size * p))
		step_start += fixnan_time_proportion
		fixnan_tot = time.time() - fixnan_start

		# Change timestep to something suitable for lpj-guess.
		log_diagnostic("Aggregating %s to hourly timestep" % out_name)
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
				(out_name, in_units, out_units))
			data = fix_units(data, in_units, out_units, output_timestep, \
				lambda p: progress_cb( \
					step_start + units_time_proportion * p))
		else:
			log_diagnostic("No units conversion required for %s" % out_name)
		unit_tot = time.time() - unit_start

		step_start += units_time_proportion

		bounds_start = time.time()
		data = bounds_checks(data, lower_bound, upper_bound, lambda p: \
			progress_cb(step_start + bounds_time_prop * p))
		bounds_tot = time.time() - bounds_start

		time_tot = read_tot + fixnan_tot + t_agg_tot + bounds_tot + unit_tot
		read_prop = 100.0 * read_tot / time_tot
		fixnan_prop = 100.0 * fixnan_tot / time_tot
		t_agg_prop = 100.0 * t_agg_tot / time_tot
		bounds_prop = 100.0 * bounds_tot / time_tot
		unit_prop = 100.0 * unit_tot / time_tot

		# Update global time estimates of all activities.
		if parallel:
			global gd_tp_lock
			gd_tp_lock.acquire()
		get_data_read_prop = read_prop
		get_data_fixnan_prop = fixnan_prop
		get_data_t_agg_prop = t_agg_prop
		get_data_bounds_prop = bounds_prop
		if not matching_units:
			get_data_unit_prop = unit_prop
		if parallel:
			# global gd_tp_lock
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
	var_lon = nc_out.variables[OUT_DIM_NAME_LON]
	var_lat = nc_out.variables[OUT_DIM_NAME_LAT]
	index_lon = index_of(var_lon, lon)
	if index_lon == -1:
		if nc_out.dimensions[OUT_DIM_NAME_LON].isunlimited():
			index_lon = len(var_lon)
		else:
			index_lon = numpy.ma.flatnotmasked_edges(var_lon)[0]
		var_lon[index_lon] = lon
	index_lat = index_of(var_lat, lat)
	if index_lat == -1:
		if nc_out.dimensions[OUT_DIM_NAME_LAT].isunlimited():
			index_lat = len(var_lat)
		else:
			index_lat = numpy.ma.flatnotmasked_edges(var_lat)[0]
		var_lat[index_lat] = lat
	return (index_lon, index_lat)

def trim_to_start_year(in_file: Dataset, timestep: int
	, data: list[float]) -> list[float]:
	"""
	Trim the given data to the start of the next year.
	"""
	start_date = ozflux_common.parse_date(in_file.time_coverage_start)
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

	for case in units_synonyms:
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
			x = ozflux_common.neighbouring_mean(data, i, N_NEIGHBOUR)
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
	timestep: int, progress_callback: Callable[[float], None]) -> list[float]:
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
	for i in range(n):
		data[i] = conversion(data[i], timestep)
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

def floats_equal(x: float, y: float) -> bool:
	"""
	Check if two floating point numbers are equal.

	@param x: The first number.
	@param y: The second number.
	"""
	EPSILON = 1e-10
	return abs(x - y) < EPSILON

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
