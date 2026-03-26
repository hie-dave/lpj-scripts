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
_TEMP_UNITS_DAVE = "degC"

# Output units for temperature variables.
_TEMP_UNITS_TRUNK = "K"

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

# Input name for relative humidity.
_IN_RH = "RH"

# Output name for relative humidity.
_OUT_RH = "relhum"

# Output units for relative humidity.
_RH_OUT_UNITS = "1"

# fixme: CLI option
_TRUNK = False

################################################################################
# Ameriflux input variable names.
################################################################################

# Name of the air temperature variable in the ameriflux input files.
_AMERIFLUX_TAIR = "TA_F"

# Name of the downwelling shortwave radiation variable in the ameriflux input files.
_AMERIFLUX_SWDOWN = "SW_IN_F"

# Name of the precipitation variable in the ameriflux input files.
_AMERIFLUX_PRECIP = "P_F"

# Name of the co2 concentration variable in the ameriflux input files.
_AMERIFLUX_CO2 = "CO2_F_MDS"

# Name of the vpd variable in the ameriflux input files.
_AMERIFLUX_VPD = "VPD_F"

# Name of the air pressure variable in the ameriflux input files.
_AMERIFLUX_PRESSURE = "PA_F"

# Name of the wind speed variable in the ameriflux input files.
_AMERIFLUX_WIND = "WS_F"

def temp_var(i: str, o: str,
             aggregator: Callable[[list[float]], float]) -> ForcingVariable:
    """
    Create a 
    """
    units = _TEMP_UNITS_TRUNK if _TRUNK else _TEMP_UNITS_DAVE
    return ForcingVariable(i, o, units, aggregator, MIN_TEMP, MAX_TEMP)

def get_ozflux_vars(timestep: int) -> list[ForcingVariable]:
    """
    Get the list of variables required for the daily grass version.

    @param timestep: Output timestep in minutes.
    """
    tair = temp_var(IN_TEMP, _OUT_TEMP, numpy.mean)
    swdown = ForcingVariable(IN_SWDOWN, _OUT_INSOL, _SWDOWN_UNITS, numpy.mean,
                             MIN_SWDOWN, MAX_SWDOWN)
    precip = ForcingVariable(IN_PRECLS, _OUT_PREC, _PREC_UNITS, numpy.sum,
                             MIN_PRECLS, MAX_PRECLS)
    co2 = ForcingVariable(IN_CO2, _OUT_CO2, _CO2_UNITS, numpy.mean, 250, 900)
    vpd = ForcingVariable(IN_VPD, _OUT_VPD, _VPD_UNITS, numpy.mean, 0, 100)
    pressure = ForcingVariable(IN_PS, _OUT_PS, _PS_UNITS, numpy.mean, MIN_PS,
                               MAX_PS)
    wind = ForcingVariable(IN_WIND, _OUT_WIND, _WIND_UNITS, numpy.mean, 0, 100)
    rh = ForcingVariable("RH", "rh", "1", numpy.mean, 0, 1)

    vars = [tair, swdown, precip, co2, vpd, pressure, wind, rh]

    # If generating a daily file, need to include tmax/tmin variables in output.
    if timestep / MINUTES_PER_HOUR == 24:
        vars.append(temp_var(IN_TEMP, _OUT_TMAX, numpy.amax))
        vars.append(temp_var(IN_TEMP, _OUT_TMIN, numpy.amin))
        vars.append(ForcingVariable(_IN_RH, _OUT_RH, _RH_OUT_UNITS,
                                    numpy.mean, 0, 100))

    return vars

def get_ameriflux_vars(timestep: int) -> list[ForcingVariable]:
    """
    Get the list of variables required for ameriflux processing.

    @param timestep: Output timestep in minutes.
    """
    tair = temp_var(_AMERIFLUX_TAIR, _OUT_TEMP, numpy.mean)
    swdown = ForcingVariable(_AMERIFLUX_SWDOWN, _OUT_INSOL, _SWDOWN_UNITS,
                             numpy.mean , MIN_SWDOWN, MAX_SWDOWN)
    precip = ForcingVariable(_AMERIFLUX_PRECIP, _OUT_PREC, _PREC_UNITS,
                             numpy.sum, MIN_PRECLS, MAX_PRECLS)
    co2 = ForcingVariable(_AMERIFLUX_CO2, _OUT_CO2, _CO2_UNITS, numpy.mean, 250,
                          900)
    vpd = ForcingVariable(_AMERIFLUX_VPD, _OUT_VPD, _VPD_UNITS, numpy.mean, 0,
                          100)
    pressure = ForcingVariable(_AMERIFLUX_PRESSURE, _OUT_PS, _PS_UNITS,
                               numpy.mean, MIN_PS, MAX_PS)
    wind = ForcingVariable(_AMERIFLUX_WIND, _OUT_WIND, _WIND_UNITS, numpy.mean,
                           0, 100)

    vars = [tair, swdown, precip, co2, vpd, pressure, wind]

    # If generating a daily file, need to include tmax/tmin variables in output.
    if timestep / MINUTES_PER_HOUR == 24:
        vars.append(temp_var(_AMERIFLUX_TAIR, _OUT_TMAX, numpy.amax))
        vars.append(temp_var(_AMERIFLUX_TAIR, _OUT_TMIN, numpy.amin))

    return vars

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
