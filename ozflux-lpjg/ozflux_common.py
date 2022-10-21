from enum import IntEnum

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
