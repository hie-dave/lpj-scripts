import datetime, math
from enum import IntEnum
from typing import Callable

## Constants

# Number of seconds in half an hour (leap seconds are a myth).
SECONDS_PER_MINUTE = 60

# NetCDF file format used for input and output.
NC_FORMAT = "NETCDF4"

# Date format used in the ozflux .nc files.
DATE_FORMAT = r"%Y-%m-%d %H:%M:%S"

# Forcing variable IDs (and column IDs too - this defines the order of
# the columns in the output file).
class Forcing(IntEnum):
	"""
	Indices of variables in the output .nc file.
	"""
	SWDOWN = 0,
	PARDF = 1,
	PARDR = 2,
	LWDOWN = 3,
	PRECLS = 4,
	PRECCO = 5,
	TAIR = 6,
	UAIR = 7,
	VAIR = 8,
	QAIR = 9,
	PSURF = 10
	# Insert new values here if required.
	NFORCINGS = 11 # This must always be the last/highest enum value

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
