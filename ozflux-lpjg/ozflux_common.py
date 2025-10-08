import datetime, math, re
from enum import IntEnum
from typing import Callable

## Constants

# Semantic version number.
VERSION = "1.0"

# Number of seconds in half an hour (leap seconds are a myth).
SECONDS_PER_MINUTE = 60

# NetCDF file format used for input and output.
NC_FORMAT = "NETCDF4"

# Date format used in the ozflux .nc files.
DATE_FORMAT = r"%Y-%m-%d %H:%M:%S"

# Times are represented in the .nc file as hours since this date.
DATE_BASELINE = datetime.datetime(1800, 1, 1)

################################################################################
# Upper and lower bounds for error checking of forcing data.
################################################################################

# Temperature lower bound (°C).
MIN_TEMP = -100

# Temperature upper bound (°C).
MAX_TEMP = 100

# Atmospheric pressure lower bound (Pa). todo: is this too conservative?
MIN_PS = 50000

# Atmospheric pressure upper bound (Pa). todo: is this too conservative?
MAX_PS = 150000

# SWDOWN lower bound (units).
MIN_SWDOWN = 0

# SWDOWN upper bound (units).
MAX_SWDOWN = 1e5

# PARDF lower bound (units).
MIN_PARDF = 0

# PARDF upper bound (units).
MAX_PARDF = 1e5

# PARDR lower bound (units).
MIN_PARDR = 0

# PARDR upper bound (units).
MAX_PARDR = 1e5

# LWDOWN lower bound (units).
MIN_LWDOWN = 0

# LWDOWN upper bound (units).
MAX_LWDOWN = 1e5

# PRECLS lower bound (units).
MIN_PRECLS = 0

# PRECLS upper bound (units).
MAX_PRECLS = 1e4

# PRECCO lower bound (units).
MIN_PRECCO = 0

# PRECCO upper bound (units).
MAX_PRECCO = 1e3

# TAIR lower bound (units).
MIN_TAIR = 200

# TAIR upper bound (units).
MAX_TAIR = 373

# UAIR lower bound (units).
MIN_UAIR = 0

# UAIR upper bound (units).
MAX_UAIR = 100

# VAIR lower bound (units).
MIN_VAIR = 0

# VAIR upper bound (units).
MAX_VAIR = 1

# QAIR lower bound (units).
MIN_QAIR = 0

# QAIR upper bound (units).
MAX_QAIR = 0.2

# PSURF lower bound (units).
MIN_PSURF = 3e4

# PSURF upper bound (units).
MAX_PSURF = 1.5e5

################################################################################
# Variable names in the input files.
################################################################################

# Name of the temperature variable in the input file.
IN_TEMP = "Ta"

# Name of the shortwave radiation variable in the input file.
IN_SWDOWN = "Fsd"

# Name of the longwave radiation variable in the input file.
IN_LWDOWN = "Fld"

# Name of the precipitation variable in the input file.
IN_PRECLS = "Precip"

# Name of the co2 concentration variable in the input file.
IN_CO2 = "CO2"

# Name of the VPD variable in the input file.
IN_VPD = "VPD"

# Name of the atmospheric pressure variable in the input file.
IN_PS = "ps"

# Name of the UAIR variable in the input file.
IN_UAIR = "Ws"

# Name of the QAIR variable in the input file.
IN_QAIR = "SH"

################################################################################
# Units conversions.
################################################################################

# g to kg. This is not really used except as an example.
KG_PER_G = 1e-3

# KPa to Pa
PA_PER_KPA = 1e3

# Number of pascals in a hectopascal.
PA_PER_HPA = 100

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

# Moles per micromole.
MOL_PER_UMOL = 1e-6

# Atomic mass of carbon (g/mol).
G_C_PER_MOL = 12.011

# Number of square metres in a hectare.
M2_PER_HA = 10_000

# Number of metres in a kilometre.
M_PER_KM = 1_000

# Number of metres in an astronomical unit (just in case :).
M_PER_AU = 1.495979e+11

# Number of centimetres in a metre.
CM_PER_M = 100

# Number of millimetres in a metre.
MM_PER_M = 1_000

# Conversion ratios from metres to various other units.
_LENGTH_CONVERSIONS = {
	"m": 1,
	"km": M_PER_KM,
	"au": M_PER_AU, # Just in case :)
	"cm": 1 / CM_PER_M,
	"mm": 1 / MM_PER_M,
}

# Common alternative names for various units which we may reasonably encounter.
_units_synonyms = [
	["W/m2", "W/m^2", "Wm^-2"],
	["kg/m2/s", "kg/m^2/s", "kgm^-2s^-1", "kg m-2 s-1"],
	["kg/m2/day", "kg/m^2/day", "mm/day"],
	["K", "k"],
	["m/s", "ms^-1"],
	["Pa", "pa"],
	["kg/kg", "kg kg-1", "", "mm/mm", "m/m", "1"], # debatable
	["ppm", "umol/mol"],
	["degC", "°C", "degrees C"],
	["umol/m2/s", "umol/m^2/s"],
	["m3/m3", "m^3/m^3"],
	["gC/m^2/day", "gC/m2/day"]
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
_units_conversions = {
	("g", "kg"): lambda x, _: x * KG_PER_G, # (as an example)
	("mm", "kg/m2/s"): lambda x, t: x / t, # Assuming mm means mm per timestep
	("kg/m2/s", "mm"): lambda x, t: t * x,
	("kg/m2/s", "mm/day"): lambda x, t: SECONDS_PER_DAY * x,
	("degC", "K"): lambda x, _: x + DEG_C_TO_K,
	("K", "degC"): lambda x, _: x - DEG_C_TO_K,
	("kPa", "Pa"): lambda x, _: x * PA_PER_KPA,
	("umol/m2/s", "kgC/m2/day"): lambda x, _: x * MOL_PER_UMOL * G_C_PER_MOL * KG_PER_G * SECONDS_PER_DAY,
	("umol/m2/s", "gC/m2/day"): lambda x, _: x * MOL_PER_UMOL * G_C_PER_MOL * SECONDS_PER_DAY,
	("umol/m2/s", "gC/m2"): lambda x, t: x * MOL_PER_UMOL * G_C_PER_MOL * t,
	("Pa", "kPa"): lambda x, _: x / PA_PER_KPA,
	("hPa", "Pa"): lambda x, _: x * PA_PER_HPA,
	("hPa", "kPa"): lambda x, _: x * PA_PER_HPA / PA_PER_KPA,
	("Mg/ha", "kg/m2"): lambda x, _: x * 0.1,
	("/ha", "/m2"): lambda x, _: x / 10000
}

# Number of seconds per day.
SECONDS_PER_DAY = SECONDS_PER_MINUTE * MINUTES_PER_HOUR * HOURS_PER_DAY

EPS = 1e-16

def _six_digit_string(x: float) -> str:
	"""
	Return a number, encoded as a 6 digit string.
	"""
	# Round x to 2 digits.
	res = str(abs(int(round(x, 2) * 1000)))
	pad = "0" * (6 - len(res))
	return "%s%s" % (pad, res)

def get_met_filename(lon: float, lat: float) -> str:
	"""
	The LSM input module expects a very specific file name for its met
	data. This is a function of the latitude and longitude of the site.
	"""
	lonstr = _six_digit_string(lon)
	latstr = _six_digit_string(lat)

	londir = "W" if lon < 0 else "E"
	latdir = "S" if lat < 0 else "N"

	return "%s%s_%s%s.nc" % (lonstr, londir, latstr, latdir)

def seconds_per_n_minutes(n: int) -> int:
	"""
	Number of seconds in N minutes.
	"""
	return n * SECONDS_PER_MINUTE

def num_timesteps_between_dates(n: int, d0: datetime.datetime, d1: datetime.datetime) -> int:
	"""
	Determine the number of n-minute timesteps between two dates.

	d1 must be later than d1.

	@param n: Timestep size in minutes.
	@param d0: The first date (must be earlier than d1).
	@param d1: The second date (must be later than d0).
	"""
	return int((d1 - d0).total_seconds() / seconds_per_n_minutes(n))

def parse_date(datestr: str) -> datetime.datetime:
	"""
	Parse a date object from a string in yyyy-MM-dd hh:mm:ss format.
	"""
	return datetime.datetime.strptime(datestr, DATE_FORMAT)

def neighbouring_indices(data: list[float], index: int, nindex
	, predicate: Callable[[float], bool]) -> list[int]:
	"""
	Return N closest neighbouring indices for the given index.
	@param index: The start point.
	@param ndata: Length of the list.
	@param nindex: Number of indices to return.
	"""
	ndata = len(data)
	if nindex > ndata:
		raise IndexError("Attempted to interpolate with more values than exist in list")

	indices = [0] * nindex
	attempt = 0
	idx = index
	for i in range(nindex):
		# Given the error check at start of function, we should
		# eventualy find another neighbour. There is probably a
		# more efficient search algorithm than this, but this will
		# be called very rarely and I have bigger fish to fry.
		def can_use(idx: int) -> bool:
			if idx < 0 or idx >= ndata or idx in indices:
				return False
			if predicate != None and not predicate(idx):
				return False
			return True

		while (not can_use(idx)):
			sign = 1 if attempt % 2 == 0 else -1
			attempt += 1
			idx = (index + (attempt + 2) // 2 * sign) % ndata

		indices[i] = idx
	return indices

def neighbouring_mean(data: list[float], index: int, n: int) -> float:
	"""
	Return the mean of the N closest neighbours in the list.
	"""
	indices = neighbouring_indices(data, index, n, lambda x: not data.mask[x])
	mean = 0
	for i in indices:
		mean += data[i]
	return mean / n

def to_metres(length: float, units: str):
	"""
	Convert the value to metres.
	@param length: The numeric value.
	@param units: Units of length.
	"""
	units = units.lower().strip()
	if units in _LENGTH_CONVERSIONS:
		return length * _LENGTH_CONVERSIONS[units]

	raise ValueError("Unknown units: '%s'" % units)

def get_length(length: str) -> float:
	"""
	Read a string with a units suffix and attempt to return the length in
	metres. E.g.

	- 5m -> 5
	- 1.024 km -> 1024

	@param length: The length with a units suffix.
	"""
	return get_lengths(length)[0]

def get_lengths(length: str) -> list[float]:
	"""
	Read a string with a units suffix and attempt to return the length in
	metres. This will throw if no lengths are be parsed, so result length is
	always >= 1. This will also only return the absolute value of the length.
	E.g.

	- 5m -> 5
	- 1.024 km -> 1024
	- -100 to -500 cm -> [1, 5]
	- 30 - 50m -> [30, 50]

	@param length: The length with a units suffix.
	"""
	conversions = "|".join(["(?:%s)" % x for x in _LENGTH_CONVERSIONS.keys()])
	pattern = r'[ \t]*-?([0-9]*\.?[0-9]+)[ \t]*'
	pattern += "(%s)?" % conversions
	pattern += "[ \t]*" # allow trailing whitespace
	matches = re.findall(pattern, length)
	if matches == None:
		raise ValueError("Cannot parse length from '%s'" % length)

	if len(matches) == 0:
		raise ValueError("Cannot parse any lengths from '%s'" % length)

	results = []
	units_prev = ""
	matches.reverse()
	for match in matches:
		value = float(match[0])
		suffix = match[1]

		if suffix == "":
			suffix = units_prev

		length_metres = to_metres(value, suffix)

		if math.isnan(length_metres):
			raise ValueError("Unable to parse length from '%s'" % length)
		results.insert(0, length_metres)
		units_prev = suffix

	return results

def less_than(x: float, y: float) -> bool:
	"""
	Return true iff x is less than y.

	@param x: Any real number.
	@param y: Any real number.
	"""
	return y - x > EPS

def greater_than(x: float, y: float) -> bool:
	"""
	Return true iff x is greater than y.

	@param x: Any real number.
	@param y: Any real number.
	"""
	return x - y > EPS

def floats_equal(x: float, y: float) -> bool:
	"""
	Return true iff x is equal to y.

	@param x: Any real number.
	@param y: Any real number.
	"""
	return abs(x - y) < EPS

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

def find_units_conversion(current_units: str, desired_units: str) \
	-> Callable[[float], float]:
	"""
	Find a conversion between two different units. Throw if not found.
	The return value is a function which takes and returns a float.
	"""
	# units_conversions is a dict mapping unit combos to a conversion.
	# units_conversions: dict[tuple[str, str], Callable[[float],float]]
	combination = (current_units, desired_units)
	if combination in _units_conversions:
		return _units_conversions[combination]
	for units_type in _units_synonyms:
		if current_units in units_type:
			for synonym in units_type:
				combination = (synonym, desired_units)
				if (combination in _units_conversions):
					return _units_conversions[combination]

	m = "No unit conversion exists from '%s' to '%s'"
	raise ValueError(m % (current_units, desired_units))

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
