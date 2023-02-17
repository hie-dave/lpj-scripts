import datetime, math
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

# Number of seconds per day.
SECONDS_PER_DAY = SECONDS_PER_MINUTE * MINUTES_PER_HOUR * HOURS_PER_DAY

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
	for i in range(nindex):
		sign = 1 if i % 2 == 0 else -1
		dlt = 1 + i // 2 * sign # 1, -1, 2, -2, ...
		idx = index + dlt

		# Given the error check at start of function, we should
		# eventualy find another neighbour. There is probably a
		# more efficient search algorithm than this, but this will
		# be called very rarely and I have bigger fish to fry.
		def can_use(idx: int) -> bool:
			if idx < 0 or idx >= ndata or idx in indices:
				return False
			if predicate != None and not predicate(data[idx]):
				return False
			return True

		while (not can_use(idx)):
			if idx < 0:
				sign = 1
			elif idx >= ndata:
				sign = -1
			idx += sign

		indices[i] = idx
	return indices

def neighbouring_mean(data: list[float], index: int, n: int) -> float:
	"""
	Return the mean of the N closest neighbours in the list.
	"""
	indices = neighbouring_indices(data, index, n, lambda x: not x.mask)
	mean = 0
	for i in indices:
		mean += data[i]
	return mean / n
