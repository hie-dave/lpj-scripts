from netCDF4 import *
from enum import Enum

# Standard names defined by the CF spec. See here for all defined values:
# https://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html

# Name of the 'standard name' attribute (from the CF spec).
ATTR_STD_NAME = "standard_name"

# Name of the 'long name' attribute (from the CF spec).
ATTR_LONG_NAME = "long_name"

# Name of the 'units' attribute (from the CF spec).
ATTR_UNITS = "units"

# Name of the 'calendar' attribute (from the CF spec).
ATTR_CALENDAR = "calendar"

# Name of the 'missing value' attribute (from the CF spec).
ATTR_MISSING_VAL = "missing_value"

# Name of the attribute which specifies the valid range of a variable.
ATTR_VALID_RANGE = "valid_range"

# Name of the attribute which specifies the minimum valid value of a variable.
ATTR_VALID_MIN = "valid_min"

# Name of the attribute which specifies the maximum valid value of a variable.
ATTR_VALID_MAX = "valid_max"

# Name of the gregorian value of the calendar attribute (from the CF spec).
CALENDAR_GREGORIAN = "gregorian"

# Standard name of the latitude variable (from the CF spec).
STD_LAT = "latitude"

# Standard name of the longitude variable (from the CF spec).
STD_LON = "longitude"

# Standard name of the time variable (from the CF spec).
STD_TIME = "time"

# Standard name of the air temperature variable (from the CF spec).
STD_TEMP = "air_temperature"

# Standard name of the shortwave radiation variable (from the CF spec).
STD_SWR_RAD = "surface_downwelling_shortwave_flux_in_air"

# Standard name of the precipitation variable (from the CF spec) when
# precipitation is given as an amount (e.g. mm).
STD_PR_AMOUNT = "precipitation_amount"

# Standard name of the precipitation variable (from the CF spec) when
# precipitation is given as a rate (e.g. mm/s).
STD_PR_RATE = "precipitation_flux"

# Standard name of the co2 concentration variable (from the CF spec).
STD_WIND = "wind_speed"

# Standard name of the co2 concentration variable (from the CF spec).
STD_CO2 = "mole_fraction_of_carbon_dioxide_in_air"

# Standard name of the vpd variable (from the CF spec).
STD_VPD = "water_vapor_saturation_deficit_in_air"

# Standard name of the atmospheric pressure variable (from the CF spec).
STD_PS = "air_pressure"

# Standard name of the maximum air temperature variable. This is not part of the
# CF spec, but it's what the dave input module requires.
STD_TMAX = "air_temperature_maximum"

# Standard name of the minimum air temperature variable. This is not part of the
# CF spec, but it's what the dave input module requires.
STD_TMIN = "air_temperature_minimum"

# Standard name of the relative humidity variable (from the CF spec).
STD_RH = "relative_humidity"

# Standard name of the specific humidity variable (from the CF spec).
STD_SH = "specific_humidity"

# Standard name of NHX nitrogen deposition variable. This is not part of the CF
# spec, but it's what the dave input module requires.
STD_NDEP_NHX = "ndep_nhx"

# Standard name of NOY nitrogen deposition variable. This is not part of the CF
# spec, but it's what the dave input module requires.
STD_NDEP_NOY = "ndep_noy"

# Long name of the latitude variable (from the CF spec).
LONG_LAT = "latitude"

# Long name of the longitude variable (from the CF spec).
LONG_LON = "longitude"

# Long name of the time variable (from the CF spec).
LONG_TIME = "time"

# Units used in the latitude variable.
UNITS_LAT = "degrees_north"

# Units used in the longitude variable.
UNITS_LON = "degrees_east"

class VariableType(Enum):
    """
    Enum for variable types. This is used to identify the type of a variable
    based on its standard name, and to determine which variable(s) are required
    for a given variable type (e.g. relative humidity requires both temperature
    and specific humidity).
    """
    TEMP = 1
    SWR = 2
    PR = 3
    CO2 = 5
    VPD = 6
    PS = 7
    TMAX = 8
    TMIN = 9
    RH = 10
    SH = 11
    NDEP_NHX = 12
    NDEP_NOY = 13
    WS = 14
    LAT = 15
    LON = 16
    TIME = 17

_STD_NAME_LOOKUP: dict[VariableType, str] = {
    VariableType.TEMP: STD_TEMP,
    VariableType.SWR: STD_SWR_RAD,
    VariableType.CO2: STD_CO2,
    VariableType.VPD: STD_VPD,
    VariableType.PS: STD_PS,
    VariableType.TMAX: STD_TMAX,
    VariableType.TMIN: STD_TMIN,
    VariableType.RH: STD_RH,
    VariableType.SH: STD_SH,
    VariableType.NDEP_NHX: STD_NDEP_NHX,
    VariableType.NDEP_NOY: STD_NDEP_NOY,
    VariableType.WS: STD_WIND,
    VariableType.LAT: STD_LAT,
    VariableType.LON: STD_LON,
    VariableType.TIME: STD_TIME
}

def _is_precip_amount(units: str) -> bool:
    """
    Check if the given units correspond to a precipitation amount variable (e.g. mm).
    """
    return units in ["mm", "kg m-2"]

def _is_precip_rate(units: str) -> bool:
    """
    Check if the given units correspond to a precipitation rate variable (e.g. mm/s).
    """
    return units in ["mm/s", "mm s-1", "kg m-2 s-1", "mm/d", "mm d-1", "mm/day", "mm day-1", "kg m-2 d-1", "kg m-2 day-1"]

def get_std_name(var_type: VariableType, out_units: str) -> str:
    """
    Get the standard name corresponding to the given variable type.
    """
    if var_type == VariableType.PR:
        if _is_precip_amount(out_units):
            return STD_PR_AMOUNT
        elif _is_precip_rate(out_units):
            return STD_PR_RATE
        else:
            raise ValueError(f"Invalid units for precipitation variable: {out_units}")
    if var_type not in _STD_NAME_LOOKUP:
        raise RuntimeError(f"No standard name defined for variable type {var_type}")
    return _STD_NAME_LOOKUP[var_type]

def get_var_from_std_name(nc: Dataset, name: str) -> Variable:
    """
    Find a variable in the input file with the specified standard name. Raise
    an exception if no matching variable is found.

    @param nc: The input netcdf file.
    @param name: The standard name for which we search.
    """
    for var in nc.variables:
        var = nc.variables[var]
        if hasattr(var, ATTR_STD_NAME):
            if getattr(var, ATTR_STD_NAME) == name:
                return var
    raise ValueError(f"No variable found in file '{nc.filepath()}' with standard_name '{name}'")

def get_units(var_id: Variable) -> str:
    """
    Get the units for the specified variable in the .nc file.
    """
    if not hasattr(var_id, ATTR_UNITS):
        raise RuntimeError(f"Variable {var_id.name} has no units attribute")
    return var_id.units

def get_calendar(var: Variable) -> str:
    if hasattr(var, ATTR_STD_NAME) and getattr(var, ATTR_STD_NAME) != STD_TIME:
        raise RuntimeError(f"Cannot get calendar attribute for variable {var.name}: not a time variable")

    if hasattr(var, ATTR_CALENDAR):
        return getattr(var, ATTR_CALENDAR)

    raise RuntimeError(f"Variable {var.name} has no calendar attribute")
