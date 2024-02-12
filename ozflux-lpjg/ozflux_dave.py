from ozflux_netcdf import *
from ozflux_common import *
import numpy

# Name of the mean temperature variable created in the dailygrass output file.
_OUT_TEMP = "tav"

# Name of the insol variable created in the dailygrass output file.
_OUT_INSOL = "insol"

# Name of the prec variable created in the dailygrass output file.
_OUT_PREC = "prec"

# Name of the co2 variable created in the dailygrass output file.
_OUT_CO2 = "co2"

# Name of the VPD variable created in the dailygrass output file.
_OUT_VPD = "vpd"

# Name of the tmax variable created in the dailygrass output file.
_OUT_TMAX = "tmax"

# Name of the tmin variable created in the dailygrass output file.
_OUT_TMIN = "tmin"

# Name of the pressure variable created in the dailygrass output file.
_OUT_PS = "ps"

# Output units for temperature variables.
_TEMP_UNITS = "degC"

# Output units for shortwave radiation.
_SWDOWN_UNITS = "W/m2"

# Output units for precipitation.
_PREC_UNITS = "mm"

# Output units for co2 concentration.
_CO2_UNITS = "ppm"

# Output units for VPD.
_VPD_UNITS = "kPa"

# Output units for atmospheric pressure.
_PS_UNITS = "Pa"

# Name of the wind speed input variable.
IN_WIND = "Ws"

# Name of the wind speed output variable.
_OUT_WIND = "wind"

# Output units for wind speed.
_WIND_UNITS = "m/s"

def temp_var(i: str, o: str, aggregator: Callable[[list[float]], float]) \
	-> ForcingVariable:
	"""
	Create a 
	"""
	return ForcingVariable(i, o, _TEMP_UNITS, aggregator, MIN_TEMP, MAX_TEMP)

def get_dailygrass_vars(timestep: int) -> list[ForcingVariable]:
	"""
	Get the list of variables required for the daily grass version.

	@param timestep: Output timestep in minutes.
	"""
	tair = temp_var(IN_TEMP, _OUT_TEMP, numpy.mean)
	swdown = ForcingVariable(IN_SWDOWN, _OUT_INSOL, _SWDOWN_UNITS, numpy.mean \
				, MIN_SWDOWN, MAX_SWDOWN)
	precip = ForcingVariable(IN_PRECLS, _OUT_PREC, _PREC_UNITS, numpy.sum \
				, MIN_PRECLS, MAX_PRECLS)
	co2 = ForcingVariable(IN_CO2, _OUT_CO2, _CO2_UNITS, numpy.mean, 250, 900)
	vpd = ForcingVariable(IN_VPD, _OUT_VPD, _VPD_UNITS, numpy.mean, 0, 100)
	pressure = ForcingVariable(IN_PS, _OUT_PS, _PS_UNITS, numpy.mean, MIN_PS \
				, MAX_PS)
	wind = ForcingVariable(IN_WIND, _OUT_WIND, _WIND_UNITS, numpy.mean, 0, 100)

	vars = [tair, swdown, precip, co2, vpd, pressure, wind]

	# If generating a daily file, need to include tmax/tmin variables in output.
	if timestep / MINUTES_PER_HOUR == 24:
		vars.append(temp_var(IN_TEMP, _OUT_TMAX, numpy.amax))
		vars.append(temp_var(IN_TEMP, _OUT_TMIN, numpy.amin))

	return vars

def create_dimensions(nc: Dataset, ngridcell: int):
	"""
	Create dimensions required by the dave-daily-grass version.

	@param nc: Output .nc file.
	@param ngridcell: The number of grid points being written to this file.
	"""
	# Dailygrass mode needs lon/lat dimensions (as it's gridded).
	# Chunk size for geographic (ie lat/lon) variables.
	geo_chunk_size = (1,)
	dim_size = 0 if UNLIMITED_DIMS else ngridcell
	create_dim_if_not_exists(nc, DIM_LON, dim_size)
	create_dim_if_not_exists(nc, DIM_LAT, dim_size)
	create_var_if_not_exists(nc, DIM_LON, FORMAT_FLOAT, DIM_LON, chunksizes = geo_chunk_size)
	create_var_if_not_exists(nc, DIM_LAT, FORMAT_FLOAT, DIM_LAT, chunksizes = geo_chunk_size)
	var_lon = nc.variables[DIM_LON]
	var_lat = nc.variables[DIM_LAT]
	var_lon.long_name = "longitude"
	var_lat.long_name = "latitude"
	var_lon.standard_name = "longitude"
	var_lat.standard_name = "latitude"
	var_lon.units = "degree_east"
	var_lat.units = "degree_north"

def create_dailygrass_variables(nc: Dataset, timestep: int
		, compression_level: int, compression_type: str):
	"""
	Create the variables required for dailygrass mode.

	@param nc: The output netcdf file.
	@param timestep: The output timestep length in minutes.
	@param compression_level: Compression level [1-9], 9 = slowest/smallest.
	@param compression_type: Compression algorithm to be used (or None).
	"""
	# TODO: configurable dimension order (achieved via constants)
	dims = (DIM_LAT, DIM_LON, DIM_TIME)
	time_chunksize = get_steps_per_year(timestep)
	chunksizes = (1, 1, time_chunksize)
	format = FORMAT_FLOAT

	for var in get_dailygrass_vars(timestep):
		create_var_if_not_exists(nc, var.out_name, format, dims
			, compression_level, compression_type, chunksizes)

	# We also need a time dimension.
	time_chunks = (time_chunksize,)
	create_var_if_not_exists(nc, DIM_TIME, FORMAT_FLOAT, DIM_TIME
		, compression_level, compression_type, time_chunks)

def write_metadata(nc: Dataset):
	"""
	Write metadata specific to the dave-daily-grass version.

	@param nc: Output .nc file.
	"""
	var_time = nc.variables[DIM_TIME]
	var_time.calendar = "gregorian"
	var_time.long_name = DIM_TIME
	var_time.standard_name = DIM_TIME
	time_start = DATE_BASELINE.strftime(DATE_FORMAT)
	var_time.units = "hours since %s" % time_start
