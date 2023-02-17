from enum import IntEnum
from ozflux_common import *
from ozflux_netcdf import *

# Forcing variable IDs (and column IDs too - this defines the order of
# the columns in the output file).
class LsmVariables(IntEnum):
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

class LsmVariable(ForcingVariable):
	def __init__(self, col: int, in_name: str, out_name: str, out_units: str
		, aggregator: Callable[[list[float]], float], lbound: float
		, ubound: float, invert: bool = False):
		super().__init__()
		self.col = col

IN_PARDF = "zero"
IN_PARDR = "zero"
IN_PRECCO = "zero"
IN_VAIR = "zero"

# Variable names for the ozflux input file.
_in_names = {
	LsmVariables.SWDOWN: IN_SWDOWN,
	LsmVariables.PARDF: IN_PARDF,
	LsmVariables.PARDR: IN_PARDR,
	LsmVariables.LWDOWN: IN_LWDOWN,
	LsmVariables.PRECLS: IN_PRECLS,
	LsmVariables.PRECCO: IN_PRECCO,
	LsmVariables.TAIR: IN_TEMP,
	LsmVariables.UAIR: IN_UAIR,
	LsmVariables.VAIR: IN_VAIR,
	LsmVariables.QAIR: IN_QAIR,
	LsmVariables.PSURF: IN_PS,
}

# Name of the single variable created in the output file.
OUT_VARIABLE_LSM = "forcing_data"

# Name of the forcings dimension in the output file.
_DIM_FORCINGS = "forcings"

# Variable names for the output file.
lsm_var_names_out = {
	LsmVariables.SWDOWN: "swdown",
	LsmVariables.PARDF: "pardf",
	LsmVariables.PARDR: "pardr",
	LsmVariables.LWDOWN: "lwdown",
	LsmVariables.PRECLS: "precls",
	LsmVariables.PRECCO: "precco",
	LsmVariables.TAIR: "tair",
	LsmVariables.UAIR: "uair",
	LsmVariables.VAIR: "vair",
	LsmVariables.QAIR: "qair",
	LsmVariables.PSURF: "psurf"
}

# From DMB:
# > Radiation (SWDOWN, PARDF, PARDR, LWDOWN)-> W/m^2
# > Precipitation (PRECLS, PRECCO) -> kg/m^2/s
# > Air temperature (TAIR) -> K
# > Wind speed (UAIR, VAIR) -> m/s (u = eastward, v = northward)
# > Air specific humidity (QAIR) -> kg/kg (unitless)
# > Pressure (PSURF) -> Pa
forcing_units = {
	LsmVariables.SWDOWN: "W/m2",
	LsmVariables.PARDF: "W/m2",
	LsmVariables.PARDR: "W/m2",
	LsmVariables.LWDOWN: "W/m2",
	LsmVariables.PRECLS: "kg/m2/s",
	LsmVariables.PRECCO: "kg/m2/s",
	LsmVariables.TAIR: "K",
	LsmVariables.UAIR: "m/s",
	LsmVariables.VAIR: "m/s",
	LsmVariables.QAIR: "",
	LsmVariables.PSURF: "Pa"
}

# Min and max values for each forcing variable, in output units.
forcing_bounds = {
	LsmVariables.SWDOWN: (MIN_SWDOWN, MAX_SWDOWN),
	LsmVariables.PARDF: (MIN_PARDF, MAX_PARDF),
	LsmVariables.PARDR: (MIN_PARDR, MAX_PARDR),
	LsmVariables.LWDOWN: (MIN_LWDOWN, MAX_LWDOWN),
	LsmVariables.PRECLS: (MIN_PRECLS, MAX_PRECLS),
	# Always zero, not used.
	LsmVariables.PRECCO: (MIN_PRECCO, MAX_PRECCO),
	# Kelvins
	LsmVariables.TAIR: (MIN_TAIR, MAX_TAIR),
	LsmVariables.UAIR: (MIN_UAIR, MAX_UAIR),
	# Always zero, not used.
	LsmVariables.VAIR: (MIN_VAIR, MAX_VAIR),
	LsmVariables.QAIR: (MIN_QAIR, MAX_QAIR),
	LsmVariables.PSURF: (MIN_PSURF, MAX_PSURF)
}

# These are the functions used to aggregate data temporally, when converting
# from one timestep to another.
indata_aggregators = {
	LsmVariables.SWDOWN: numpy.mean,
	LsmVariables.PARDF: numpy.mean,
	LsmVariables.PARDR: numpy.mean,
	LsmVariables.LWDOWN: numpy.mean,
	LsmVariables.PRECLS: numpy.sum, # precip in mm
	LsmVariables.PRECCO: numpy.sum, # precip in mm
	LsmVariables.TAIR: numpy.mean,
	LsmVariables.UAIR: numpy.mean,
	LsmVariables.VAIR: numpy.mean,
	LsmVariables.QAIR: numpy.mean,
	LsmVariables.PSURF: numpy.mean,
}

def lsm_var(forcing: LsmVariables, out_name: str, out_units: str) -> LsmVariable:
	"""
	Create an LSM variable instance.
	"""
	return LsmVariable(forcing
		, _in_names[forcing]
		, out_name
		, out_units
		, indata_aggregators[forcing]
		, forcing_bounds[forcing][0]
		, forcing_bounds[forcing][1])

def get_lsm_vars() -> list[LsmVariable]:
	"""
	Get the LSM variables as a list of forcing variable instances.
	"""
	vars = []
	for forcing in LsmVariables:
		if forcing == LsmVariables.NFORCINGS:
			break
		out_name = lsm_var_names_out[forcing]
		var = lsm_var(forcing, out_name, forcing_units[forcing])
		vars.append(var)

def create_dimensions(nc: Dataset):
	"""
	Create the required dimensions in the output file.
	"""
	nc.createDimension(_DIM_FORCINGS, LsmVariables.NFORCINGS)

def create_variables(nc: Dataset):
	"""
	Create the required variables in the output file.

	## todo: call create_variable_if_not_exists() to unify compression code.
	"""
	# LSM mode only requires a single variable.
	nc.createVariable(OUT_VARIABLE_LSM, FORMAT_FLOAT \
	, (DIM_TIME, _DIM_FORCINGS), compression = 'zlib')

def write_data(file: Dataset, data: list[float], forcing: LsmVariable, \
	pcb: Callable[[float], None]):
	"""
	Write data to the output file (LSM mode).

	@param variable: Variable in the output file, to which data will be written.
	@param data: A column of data to be written.
	@param forcing: Column of the variable to which the data will be written.
	@pcb: Called to report progress.
	"""
	n = len(data)

	# Always use at least this many chunks.
	MIN_NUM_CHUNKS = 10

	variable = file.variables[OUT_VARIABLE_LSM]

	chunk_size = min(CHUNK_SIZE, n / MIN_NUM_CHUNKS)
	for i in range(0, n, chunk_size):
		row = i
		upper = min(n, i + chunk_size)
		col = forcing.col
		variable[row:upper, col] = data[i:upper]
		pcb(i / n)

	log_debug("Successfully wrote %d items of %s" % (n, forcing.out_name))

def write_metadata(nc: Dataset):
	"""
	Write LSM-specific metadata to the output file.
	"""
	vars = get_lsm_vars()
	i = 0
	for i in range(len(vars)):
		oname = vars[i].out_name
		units = vars[i].out_units
		attr_name = "col_%d_%s_units" % (i, oname)
		setattr(nc, attr_name, units)
