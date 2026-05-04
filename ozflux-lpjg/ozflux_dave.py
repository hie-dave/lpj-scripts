from ozflux_netcdf import *
from ozflux_common import *
from ozflux_time import *
from ozflux_cf import *

import numpy

# TODO: CLI (ie runtime) option
_TRUNK = False

################################################################################
# Variable names in the output files.
################################################################################

_OUTPUT_NAMES: dict[VariableType, str] = {
    VariableType.TEMP: "tav",
    VariableType.SWR: "insol",
    VariableType.PR: "prec",
    VariableType.CO2: "co2",
    VariableType.VPD: "vpd",
    VariableType.PS: "ps",
    VariableType.WS: "wind",
    VariableType.RH: "relhum",
    VariableType.TMIN: "tmin",
    VariableType.TMAX: "tmax",
    VariableType.SH: "sh",
    VariableType.NDEP_NHX: "ndep_nhx",
    VariableType.NDEP_NOY: "ndep_noy",
    VariableType.LAT: "lat",
    VariableType.LON: "lon",
    VariableType.TIME: "time"
}

################################################################################
# Output units.
################################################################################

_OUTPUT_UNITS_DAVE: dict[VariableType, str] = {
    VariableType.TEMP: "degC",
    VariableType.SWR: "W/m2",
    VariableType.PR: "mm",
    VariableType.CO2: "ppm",
    VariableType.VPD: "kPa",
    VariableType.PS: "Pa",
    VariableType.WS: "m/s",
    VariableType.RH: "1",
    VariableType.TMIN: "degC",
    VariableType.TMAX: "degC",
    VariableType.SH: "kg/kg"
}

_OUTPUT_UNITS_CFINPUT: dict[VariableType, str] = {
    VariableType.TEMP: "K",
    VariableType.SWR: "W m-2",
    VariableType.PR: "kg m-2",
    # VariableType.CO2: "ppm", # not supported by cfinput
    # VariableType.VPD: "kPa"  # not supported by cfinput
    VariableType.PS: "Pa",
    VariableType.WS: "m s-1",
    VariableType.RH: "1",
    VariableType.TMIN: "K",
    VariableType.TMAX: "K",
    VariableType.SH: "1"
}

################################################################################
# Variable names in the input files.
################################################################################

_INPUT_NAMES: dict[VariableType, str] = {
    VariableType.TEMP: "Ta",
    VariableType.SWR: "Fsd",
    VariableType.PR: "Precip",
    VariableType.CO2: "CO2",
    VariableType.VPD: "VPD",
    VariableType.PS: "ps",
    VariableType.WS: "Ws",
    VariableType.RH: "RH",
    VariableType.SH: "SH"
}

################################################################################
# Ameriflux input variable names.
################################################################################

_INPUT_NAMES_AMERIFLUX: dict[VariableType, str] = {
    VariableType.TEMP: "TA_F",
    VariableType.SWR: "SW_IN_F",
    VariableType.PR: "P_F",
    VariableType.CO2: "CO2_F_MDS",
    VariableType.VPD: "VPD_F",
    VariableType.PS: "PA_F",
    VariableType.WS: "WS_F"
}

################################################################################
# Functions to be called to temporally aggregate the data (e.g. 1hr -> 3hr). The
# data will be in its output units when these are called.
################################################################################

_AGGREGATORS: dict[VariableType, Callable[[list[float]], float]] = {
    VariableType.TEMP: numpy.mean,
    VariableType.SWR: numpy.mean,
    VariableType.PR: numpy.sum, # Precip is always an amount in output units
    VariableType.CO2: numpy.mean,
    VariableType.VPD: numpy.mean,
    VariableType.PS: numpy.mean,
    VariableType.WS: numpy.mean,
    VariableType.RH: numpy.mean,
    VariableType.TMIN: numpy.amin,
    VariableType.TMAX: numpy.amax,
    VariableType.SH: numpy.mean,
    VariableType.NDEP_NHX: numpy.mean,
    VariableType.NDEP_NOY: numpy.mean
}

################################################################################
# Bounds for each variable type, expressed in output units. Any data outside these
# bounds will be set to the bound value.
################################################################################

# Air temperature lower bounds (degC).
_MIN_TEMP = -100

# Air temperature upper bounds (degC).
_MAX_TEMP = 100

# Upper and lower bounds for each variable type, expressed in dave units.
_BOUNDS: dict[VariableType, tuple[float, float]] = {
    VariableType.TEMP: (_MIN_TEMP, _MAX_TEMP),
    VariableType.SWR: (0, 1e5),
    VariableType.PR: (0, 1e4),
    VariableType.CO2: (250, 900),
    VariableType.VPD: (0, 100),
    VariableType.PS: (50_000, 150_000),
    VariableType.WS: (0, 100),
    VariableType.RH: (0, 1),
    VariableType.TMIN: (_MIN_TEMP, _MAX_TEMP),
    VariableType.TMAX: (_MIN_TEMP, _MAX_TEMP),
    VariableType.SH: (0, 1),
    VariableType.NDEP_NHX: (0, 10),
    VariableType.NDEP_NOY: (0, 10)
}

"""
Check if a variable type corresponds to an extensive variable (i.e. one that
should be summed when temporally aggregating, rather than averaged).
"""
def is_extensive(type: VariableType) -> bool:
    if type == VariableType.PR:
        return True
    return False

"""
Get the data padder used to pad out missing values.
"""
def get_padder(type: VariableType) -> DataPadder:
    if is_extensive(type):
        return ValuePadder(0)
    return NearestValuePadder()

def ozflux_var(type: VariableType, timestep: int) -> ForcingVariable:
    """
    Create a ForcingVariable with the specified parameters.

    @param type: Variable type.
    @param timestep: Output timestep in seconds.
    """
    in_name = _INPUT_NAMES[type]
    out_name = _OUTPUT_NAMES[type]
    out_units = _OUTPUT_UNITS_CFINPUT[type] if _TRUNK else _OUTPUT_UNITS_DAVE[type]
    aggregator = _AGGREGATORS[type]
    padder = get_padder(type)
    lower_bound = _BOUNDS[type][0]
    upper_bound = _BOUNDS[type][1]
    std_name = get_std_name(type, out_units)

    # Sometimes the bounds may need conversion to the actual output units.
    out_units_dave = _OUTPUT_UNITS_DAVE[type]
    if not units_match(out_units, out_units_dave):
        conversion = find_units_conversion(out_units_dave, out_units)
        lower_bound = conversion(lower_bound, timestep)
        upper_bound = conversion(upper_bound, timestep)

    return ForcingVariable(in_name, out_name, out_units, aggregator, padder,
                           lower_bound, upper_bound, std_name)

def ameriflux_var(type: VariableType, timestep: int) -> ForcingVariable:
    """
    Create a ForcingVariable for ameriflux processing with the specified parameters.

    @param type: Variable type.
    @param timestep: Output timestep in seconds.
    """
    in_name = _INPUT_NAMES_AMERIFLUX[type]
    out_name = _OUTPUT_NAMES[type]
    out_units = _OUTPUT_UNITS_CFINPUT[type] if _TRUNK else _OUTPUT_UNITS_DAVE[type]
    aggregator = _AGGREGATORS[type]
    padder = get_padder(type)
    lower_bound = _BOUNDS[type][0]
    upper_bound = _BOUNDS[type][1]
    std_name = get_std_name(type, out_units)

    # Sometimes the bounds may need conversion to the actual output units.
    out_units_dave = _OUTPUT_UNITS_DAVE[type]
    if not units_match(out_units, out_units_dave):
        conversion = find_units_conversion(out_units_dave, out_units)
        lower_bound = conversion(lower_bound, timestep)
        upper_bound = conversion(upper_bound, timestep)

    return ForcingVariable(in_name, out_name, out_units, aggregator, padder,
                           lower_bound, upper_bound, std_name)

def get_ozflux_vars(timestep: int) -> list[ForcingVariable]:
    """
    Get the list of variables required for the daily grass version.

    @param timestep: Output timestep in minutes.
    """
    types = [VariableType.TEMP, VariableType.SWR, VariableType.PR,
             VariableType.CO2, VariableType.VPD, VariableType.PS,
             VariableType.WS, VariableType.RH]

    # If generating a daily file, need to include tmax/tmin variables in output.
    # I'm including RH for good measure, in case cf input is being used.
    if timestep / MINUTES_PER_HOUR == 24:
        types.append(VariableType.TMAX)
        types.append(VariableType.TMIN)
        types.append(VariableType.RH)

    ts = timestep * SECONDS_PER_MINUTE
    vars = [ozflux_var(v, ts) for v in types]

    return vars

def get_ameriflux_vars(timestep: int) -> list[ForcingVariable]:
    """
    Get the list of variables required for ameriflux processing.

    @param timestep: Output timestep in minutes.
    """
    types = [VariableType.TEMP, VariableType.SWR, VariableType.PR,
             VariableType.CO2, VariableType.VPD, VariableType.PS,
             VariableType.WS]

    # If generating a daily file, need to include tmax/tmin variables in output.
    if timestep / MINUTES_PER_HOUR == 24:
        types.append(VariableType.TMAX)
        types.append(VariableType.TMIN)
        # FIXME: probably need RH in this case.

    ts = timestep * SECONDS_PER_MINUTE
    vars = [ameriflux_var(v, ts) for v in types]

    return vars
