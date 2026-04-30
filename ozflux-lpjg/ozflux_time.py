from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
import datetime, numpy
from ozflux_netcdf import *
from netCDF4 import Dataset, num2date, date2num

# Abstract base class for padding options. Implementations may pad by repeating
# nearest value, extrapolating, etc.
class DataPadder:
    def pad(self, values: numpy.ndarray, npad_start: int,
            npad_end: int) -> numpy.ndarray:
        raise NotImplementedError()

class ExtrapolationPadder(DataPadder):
    def pad(self, values: numpy.ndarray, npad_start: int,
            npad_end: int) -> numpy.ndarray:
        if npad_start > 0:
            # Extrapolate linearly from the first two values.
            delta = values[1] - values[0]
            padding = values[0] - delta * numpy.arange(npad_start, 0, -1)
            values = numpy.concatenate([padding, values])

        if npad_end > 0:
            # Extrapolate linearly from the last two values.
            delta = values[-1] - values[-2]
            padding = values[-1] + delta * numpy.arange(1, npad_end + 1)
            values = numpy.concatenate([values, padding])

        return values

class NearestValuePadder(DataPadder):
    def pad(self, values: numpy.ndarray, npad_start: int,
            npad_end: int) -> numpy.ndarray:
        if npad_start > 0:
            # Repeat the first value npad_start times at the start of the array.
            padding = numpy.repeat(values[0], npad_start)
            values = numpy.concatenate([padding, values])

        if npad_end > 0:
            # Repeat the last value npad_end times at the end of the array.
            padding = numpy.repeat(values[-1], npad_end)
            values = numpy.concatenate([values, padding])

        return values

class ValuePadder(DataPadder):
    def __init__(self, value):
        self.value = value

    def pad(self, values: numpy.ndarray, npad_start: int,
            npad_end: int) -> numpy.ndarray:
        if npad_start > 0:
            padding = numpy.repeat(self.value, npad_start)
            values = numpy.concatenate([padding, values])

        if npad_end > 0:
            padding = numpy.repeat(self.value, npad_end)
            values = numpy.concatenate([values, padding])

        return values

class BoundaryMode(str, Enum):
    # Ignore missing timesteps between first/last timestep and a year boundary.
    IGNORE = "ignore"
    # Pad missing timesteps between first/last timestep and a year boundary with
    # the nearest available value.
    PAD = "pad"
    # Trim missing timesteps between first/last timestep and a year boundary.
    TRIM = "trim"
    # Raise an error if there are missing timesteps between first/last timestep
    # and a year boundary.
    ERROR = "error"

class AlignmentMode(str, Enum):
    # Use the alignment of the input data.
    NONE = "none"
    # Align to the nearest half hour.
    HALF_HOUR = "half_hour"
    # Align to the nearest hour (only valid with timestep >= hourly).
    STRICT_HOUR = "strict_hour"

@dataclass(frozen=True)
class TimePolicy:
    # Desired output timestep in minutes. Must be a multiple of the input
    # timestep.
    output_timestep_minutes: int
    # Behavior when the input time variable is missing timesteps at the start of
    # the year.
    start_boundary_mode: BoundaryMode = BoundaryMode.PAD
    # Behavior when the input time variable is missing timesteps at the end of
    # the year.
    end_boundary_mode: BoundaryMode = BoundaryMode.IGNORE
    # Alignment mode for the output timeseries. This determines how the output
    # time values are aligned to calendar boundaries (e.g. half hour, hour).
    alignment_mode: AlignmentMode = AlignmentMode.NONE
    # Maximum number of timesteps to be prepended when missing first timestep(s)
    # of year at the start of the timeseries.
    max_pad: int = 1

@dataclass(frozen=True)
class InputTimeInfo:
    # Units of the input time variable.
    units: str
    # Calendar type of the input time variable.
    calendar: str
    # Timestep of the input time variable, in seconds.
    input_timestep_seconds: int
    # First timestep in the input time variable.
    first_raw: datetime.datetime
    # Last timestep in the input time variable.
    last_raw: datetime.datetime
    # Total number of timesteps in the input time variable.
    n_raw: int
    # Number of timesteps between the first timestep and a year boundary.
    # 0 iff the first timestep is on a year boundary.
    nmissing_start: int
    # Number of timesteps between the last timestep and a year boundary.
    # 0 iff the last timestep is on a year boundary.
    nmissing_end: int
    # The raw datetime values of the input time variable, in order.
    times: numpy.ndarray[datetime.datetime]

@dataclass(frozen=True)
class AggregationInfo:
    # Number of values to be trimmed from the start of the input data.
    ntrim_start: int
    # Number of values to be trimmed from the end of the input data.
    ntrim_end: int
    # Number of values to be padded at the start of the input data.
    npad_start: int
    # Number of values to be padded at the end of the input data.
    npad_end: int
    # Ratio of output timestep to input timestep.
    aggregate_ratio: int
    # Timestep of the output data in seconds.
    timestep: int

@dataclass(frozen=True)
class TimePlan:
    # Time policy used to generate this plan.
    policy: TimePolicy
    # Information about the input time variable.
    input: InputTimeInfo
    # Information about aggregation to be applied to data variables.
    aggregation: AggregationInfo
    # Start of the input timeseries.
    source_start: datetime.datetime
    # End of the input timeseries.
    source_end: datetime.datetime
    # Start of the output timeseries (after trim/align/pad/aggregate).
    output_start: datetime.datetime
    # End of the output timeseries (after trim/align/pad/aggregate).
    output_end: datetime.datetime
    # Number of values in the output timeseries.
    output_n: int
    # Date time values for the output time variable, in order.
    output_time_values: numpy.ndarray
    # Units for the output time variable.
    output_time_units: str
    # Calendar type for the output time variable.
    output_time_calendar: str

def align_datetime(dt: datetime.datetime,
                   quantum: datetime.timedelta) -> datetime.datetime:
    # Align the datetime to the nearest multiple of the quantum.
    day_start = datetime.datetime(dt.year, dt.month, dt.day)
    offset = dt - day_start
    q = quantum.total_seconds()
    n = round(offset.total_seconds() / q)
    return day_start + datetime.timedelta(seconds=n * q)

def align_times(times: numpy.ndarray, mode: AlignmentMode) -> numpy.ndarray:
    if mode == AlignmentMode.NONE:
        return times
    elif mode == AlignmentMode.HALF_HOUR:
        quantum = datetime.timedelta(minutes = 30)
    elif mode == AlignmentMode.STRICT_HOUR:
        quantum = datetime.timedelta(hours = 1)
    else:
        raise ValueError(f"Invalid alignment mode: {mode}")

    return numpy.array([align_datetime(t, quantum) for t in times])

def inspect_input_time(nc: Dataset) -> InputTimeInfo:
    time = get_var_from_std_name(nc, STD_TIME)
    if not hasattr(time, ATTR_UNITS):
        raise RuntimeError("Time variable is missing 'units' attribute")
    if not hasattr(time, ATTR_CALENDAR):
        raise RuntimeError("Time variable is missing 'calendar' attribute")

    units = getattr(time, ATTR_UNITS)
    calendar = getattr(time, ATTR_CALENDAR)
    times = num2date(time[:], units, calendar,
                     only_use_cftime_datetimes = False,
                     only_use_python_datetimes = True)
    if len(times) == 0:
        raise RuntimeError("Input time variable contains no timesteps")
    if len(times) < 2:
        raise RuntimeError("Input time variable contains only one timestep; at least two are required to determine the timestep length")

    timestep_seconds = int((times[1] - times[0]).total_seconds())
    if timestep_seconds <= 0:
        raise RuntimeError(f"Input time variable has non-positive timestep of {timestep_seconds} seconds between first two timesteps")
    if SECONDS_PER_DAY % timestep_seconds != 0:
        raise RuntimeError(f"Input time variable has timestep of {timestep_seconds} seconds, which does not divide evenly into a day")

    # Verify that all timesteps are the same.
    mismatch = numpy.diff(times) != datetime.timedelta(seconds=timestep_seconds)
    if numpy.any(mismatch):
        non_uniform_indices = numpy.where(mismatch)[0]
        msg = f"Input time variable contains non-uniform timesteps"
        msg += f"\n  Expected timestep: {timestep_seconds} seconds"
        msg += f"\n  Non-uniform timestep indices: {non_uniform_indices}"
        for i in non_uniform_indices:
            msg += f"\n    Index {i}: {times[i]} to {times[i+1]}"
            delta = int((times[i+1] - times[i]).total_seconds())
            msg += f" ({delta} seconds)"
        raise RuntimeError(msg)

    # Get number of missing timesteps.
    # Could probably optimise this by adding nsteps_per_day for all doy days
    # except for the date of times[0].
    delta = datetime.timedelta(seconds = timestep_seconds)
    dt = times[0] - delta
    nmissing_start = 0
    while dt.year == times[0].year:
        nmissing_start += 1
        dt -= delta
    assert(nmissing_start >= 0)

    dt = times[-1] + delta
    nmissing_end = 0
    while dt.year == times[-1].year:
        nmissing_end += 1
        dt += delta
    assert(nmissing_end >= 0)

    return InputTimeInfo(
        units = units,
        calendar = calendar,
        input_timestep_seconds = timestep_seconds,
        first_raw = times[0],
        last_raw = times[-1],
        n_raw = len(times),
        nmissing_start = int(nmissing_start),
        nmissing_end = int(nmissing_end),
        times = times
    )

def make_time_plan(nc: Dataset, policy: TimePolicy) -> TimePlan:
    tinfo = inspect_input_time(nc)

    ntrim_start = 0
    ntrim_end = 0
    npad_start = 0
    npad_end = 0
    times = tinfo.times

    if len(times) < 2:
        raise RuntimeError(f"Input time variable contains only {len(times)} timestep(s); at least two are required to determine the timestep length")

    if tinfo.nmissing_start > 0:
        log_diagnostic(f"Input time variable is missing {tinfo.nmissing_start} timestep(s) at the start of the year (first time is {tinfo.first_raw})")

    if tinfo.nmissing_start == 0:
        log_diagnostic(f"Input time variable starts on year boundary")

    elif policy.start_boundary_mode == BoundaryMode.ERROR:
        raise RuntimeError(f"Input data contains missing timesteps at start")

    elif policy.start_boundary_mode == BoundaryMode.PAD:
        if tinfo.nmissing_start > policy.max_pad:
            raise RuntimeError(f"Number of missing timesteps exceeds the maximum pad length of {policy.max_pad}")

        # Total number of seconds to pad.
        npad_start = tinfo.nmissing_start
        pad_s = float(tinfo.input_timestep_seconds * tinfo.nmissing_start)
        output_start = tinfo.first_raw - datetime.timedelta(seconds = pad_s)

        log_diagnostic(f"Input data will be padded at start with {tinfo.nmissing_start} timestep(s) to {output_start}")

    elif policy.start_boundary_mode == BoundaryMode.TRIM:

        # Get start of next year.
        dt = tinfo.first_raw
        delta = datetime.timedelta(seconds = tinfo.input_timestep_seconds)
        while dt.year == tinfo.first_raw.year:
            dt += delta
            ntrim_start += 1
        if ntrim_start > len(times):
            raise RuntimeError(f"Cannot trim {ntrim_start} timesteps from start of timeseries with only {len(times)} timesteps")
        output_start = dt
        log_diagnostic(f"First {ntrim_start} timestep(s) of input data will be trimmed; data will start at {output_start}")
    elif policy.start_boundary_mode == BoundaryMode.IGNORE:
        log_diagnostic(f"Ignoring missing steps and starting at {tinfo.first_raw}")

    if tinfo.nmissing_end > 0:
        log_diagnostic(f"Input time variable is missing {tinfo.nmissing_end} timestep(s) at the end of the year (last time is {tinfo.last_raw})")

    if tinfo.nmissing_end == 0:
        log_diagnostic(f"Input time variable ends on year boundary")

    elif policy.end_boundary_mode == BoundaryMode.ERROR:
        raise RuntimeError(f"Input data contains missing timesteps at end")

    elif policy.end_boundary_mode == BoundaryMode.PAD:
        # Total number of seconds to pad.
        npad_end = tinfo.nmissing_end
        pad_s = float(tinfo.input_timestep_seconds * tinfo.nmissing_end)
        output_end = tinfo.last_raw + datetime.timedelta(seconds = pad_s)

        log_diagnostic(f"Input data will be padded at end with {tinfo.nmissing_end} timestep(s) to {output_end}")

    elif policy.end_boundary_mode == BoundaryMode.TRIM:
        # Get end of next year.
        dt = tinfo.last_raw
        delta = datetime.timedelta(seconds = tinfo.input_timestep_seconds)
        while dt.year == tinfo.last_raw.year:
            dt -= delta
            ntrim_end += 1
        if ntrim_end > len(times):
            raise RuntimeError(f"Cannot trim {ntrim_end} timesteps from end of timeseries with only {len(times)} timesteps")

        output_end = dt
        log_diagnostic(f"Last {ntrim_end} timestep(s) of input data will be trimmed; data will end at {output_end}")

    elif policy.end_boundary_mode == BoundaryMode.IGNORE:
        log_diagnostic(f"Ignoring missing timesteps at end of timeseries and ending at {tinfo.last_raw}")

    output_timestep_seconds = policy.output_timestep_minutes * SECONDS_PER_MINUTE
    if output_timestep_seconds % tinfo.input_timestep_seconds != 0:
        raise RuntimeError(f"Output timestep of {policy.output_timestep_minutes} minutes is not a multiple of input timestep of {tinfo.input_timestep_seconds} seconds")
    if policy.alignment_mode == AlignmentMode.STRICT_HOUR and output_timestep_seconds < SECONDS_PER_HOUR:
        raise RuntimeError(f"STRICT_HOUR alignment mode requires output timestep >= 1 hour, but output timestep is only {policy.output_timestep_minutes} minutes")

    aggregate_ratio = int(output_timestep_seconds // tinfo.input_timestep_seconds)

    # Aggregate to new timestep if needed.
    if output_timestep_seconds != tinfo.input_timestep_seconds:
        log_diagnostic(f"Output timestep of {policy.output_timestep_minutes} minutes is a multiple of input timestep of {tinfo.input_timestep_seconds} seconds; data will be aggregated with ratio {aggregate_ratio}")

    aggregation = AggregationInfo(ntrim_start, ntrim_end, npad_start, npad_end,
                                  aggregate_ratio, output_timestep_seconds)

    # First time in each window becomes the new timestamp.
    aggregator = lambda values, axis: numpy.take(values, 0, axis = axis)
    padder = ExtrapolationPadder()
    times = apply_time_plan(tinfo.times, aggregation, padder, aggregator)

    # Align times as per user-specified alignment mode.
    times = align_times(times, policy.alignment_mode)

    # Validate alignment.
    if len(set(times)) != len(times):
        raise RuntimeError(f"Alignment mode {policy.alignment_mode} causes duplicate timestamps in output time variable")
    expected = datetime.timedelta(seconds = output_timestep_seconds)
    if numpy.any(numpy.diff(times) != expected):
        raise RuntimeError(f"Alignment mode {policy.alignment_mode} causes non-uniform timesteps in output time variable")

    return TimePlan(
        policy = policy,
        input = tinfo,
        aggregation = aggregation,
        source_start = tinfo.first_raw,
        source_end = tinfo.last_raw,
        output_start = times[0],
        output_end = times[-1],
        output_n = len(times),
        output_time_values = times,
        output_time_units = tinfo.units,
        output_time_calendar = tinfo.calendar
    )

def apply_time_plan(values: numpy.ndarray,
                    options: AggregationInfo,
                    padder: DataPadder,
                    aggregator: Callable[[numpy.ndarray, int], float]) -> numpy.ndarray:
    """
    Pads, trims, and aggregates a single variable timeseries according to plan.
    No unit conversion/bounds checks here.
    """
    # 1. Pad.
    # TODO: ensure that rainfall uses NumericPadder with 0.
    values = padder.pad(values, options.npad_start, options.npad_end)

    # 2. Trim.
    if options.ntrim_start > 0:
        values = values[options.ntrim_start:]
    if options.ntrim_end > 0:
        values = values[:len(values) - options.ntrim_end]

    # 3. Aggregate.
    if options.aggregate_ratio > 1:
        # Trim any extra values at the end that don't fit into a full group.
        values = values[:len(values) - (len(values) % options.aggregate_ratio)]

        # Reshape to (n_output_timesteps, aggregate_ratio).
        values = values.reshape(-1, options.aggregate_ratio)
        values = aggregator(values, axis = 1)

    return values

def validate_time_plan(plan: TimePlan):
    # We cannot trim and pad. This should never happen because they require
    # mutually exclusive modes, but check just in case.
    if plan.aggregation.ntrim_start > 0 and plan.aggregation.npad_start > 0:
        raise RuntimeError(f"Cannot trim {plan.aggregation.ntrim_start} values and pad {plan.aggregation.npad_start} values to the start of the timeseries at the same time")

    if plan.aggregation.ntrim_end > 0 and plan.aggregation.npad_end > 0:
        raise RuntimeError(f"Cannot trim {plan.aggregation.ntrim_end} values and pad {plan.aggregation.npad_end} values to the end of the timeseries at the same time")

def index_to_datetime(plan: TimePlan, idx: int) -> datetime.datetime:
    return plan.output_time_values[idx]

def datetime_to_index(plan: TimePlan, when: datetime.datetime) -> int:
    idx = numpy.searchsorted(plan.output_time_values, when)
    if idx == len(plan.output_time_values):
        raise ValueError(f"Datetime {when} is after end of output time range")
    if plan.output_time_values[idx] != when:
        raise ValueError(f"Datetime {when} is not on output time grid")
    return idx
