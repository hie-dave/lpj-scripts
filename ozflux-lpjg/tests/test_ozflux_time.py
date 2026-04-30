import datetime
import sys
from dataclasses import replace
from pathlib import Path

import numpy
import pytest
from netCDF4 import Dataset, date2num

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ozflux_time import (  # noqa: E402
    AggregationInfo,
    AlignmentMode,
    BoundaryMode,
    ExtrapolationPadder,
    NearestValuePadder,
    TimePolicy,
    ValuePadder,
    align_datetime,
    align_times,
    apply_time_plan,
    datetime_to_index,
    index_to_datetime,
    inspect_input_time,
    make_time_plan,
    validate_time_plan,
)

TIME_UNITS = "hours since 2000-01-01 00:00:00"
TIME_CALENDAR = "gregorian"

def dt(year, month, day, hour=0, minute=0, second=0):
    return datetime.datetime(year, month, day, hour, minute, second)

def date_range(start, count, step):
    return numpy.array([start + i * step for i in range(count)])

def write_time_nc(path, dates, units=TIME_UNITS, calendar=TIME_CALENDAR,
                  with_units=True, with_calendar=True, standard_name="time"):
    with Dataset(path, "w") as nc:
        nc.createDimension("time", len(dates))
        time = nc.createVariable("time", "f8", ("time",))
        if standard_name is not None:
            time.standard_name = standard_name
        if with_units:
            time.units = units
        if with_calendar:
            time.calendar = calendar
        if len(dates) > 0:
            time[:] = date2num(dates, units, calendar)

def inspect_dates(tmp_path, dates, **kwargs):
    path = tmp_path / "input.nc"
    write_time_nc(path, dates, **kwargs)
    with Dataset(path) as nc:
        return inspect_input_time(nc)

def make_plan(tmp_path, dates, policy):
    path = tmp_path / "input.nc"
    write_time_nc(path, dates)
    with Dataset(path) as nc:
        return make_time_plan(nc, policy)

def assert_datetimes_equal(actual, expected):
    assert list(actual) == list(expected)

def test_inspect_accepts_phase_offset_on_first_timestep(tmp_path):
    dates = date_range(dt(2000, 1, 1, 0, 5), 3, datetime.timedelta(hours=1))

    info = inspect_dates(tmp_path, dates)

    assert info.input_timestep_seconds == 3600
    assert info.first_raw == dt(2000, 1, 1, 0, 5)
    assert info.last_raw == dt(2000, 1, 1, 2, 5)
    assert info.n_raw == 3
    assert info.nmissing_start == 0

def test_inspect_counts_missing_start_on_input_phase_grid(tmp_path):
    dates = date_range(dt(2000, 1, 1, 1, 5), 3, datetime.timedelta(hours=1))

    info = inspect_dates(tmp_path, dates)

    assert info.nmissing_start == 1

def test_inspect_counts_missing_end_on_input_phase_grid(tmp_path):
    dates = date_range(dt(2000, 12, 31, 21, 5), 2, datetime.timedelta(hours=1))

    info = inspect_dates(tmp_path, dates)

    assert info.nmissing_end == 1

def test_inspect_requires_units_attribute(tmp_path):
    dates = date_range(dt(2000, 1, 1), 2, datetime.timedelta(hours=1))

    with pytest.raises(RuntimeError, match="units"):
        inspect_dates(tmp_path, dates, with_units=False)

def test_inspect_requires_calendar_attribute(tmp_path):
    dates = date_range(dt(2000, 1, 1), 2, datetime.timedelta(hours=1))

    with pytest.raises(RuntimeError, match="calendar"):
        inspect_dates(tmp_path, dates, with_calendar=False)

def test_inspect_rejects_empty_time_axis(tmp_path):
    with pytest.raises(RuntimeError, match="no timesteps"):
        inspect_dates(tmp_path, [])

def test_inspect_rejects_single_timestep(tmp_path):
    with pytest.raises(RuntimeError, match="only one timestep"):
        inspect_dates(tmp_path, [dt(2000, 1, 1)])

def test_inspect_rejects_duplicate_first_two_timesteps(tmp_path):
    dates = numpy.array([dt(2000, 1, 1), dt(2000, 1, 1)])

    with pytest.raises(RuntimeError, match="non-positive timestep"):
        inspect_dates(tmp_path, dates)

def test_inspect_rejects_descending_first_two_timesteps(tmp_path):
    dates = numpy.array([dt(2000, 1, 1, 1), dt(2000, 1, 1)])

    with pytest.raises(RuntimeError, match="non-positive timestep"):
        inspect_dates(tmp_path, dates)

def test_inspect_rejects_timestep_that_does_not_divide_day(tmp_path):
    dates = date_range(dt(2000, 1, 1), 2, datetime.timedelta(minutes=7))

    with pytest.raises(RuntimeError, match="does not divide evenly into a day"):
        inspect_dates(tmp_path, dates)

def test_inspect_rejects_non_uniform_timestep(tmp_path):
    dates = numpy.array([
        dt(2000, 1, 1, 0),
        dt(2000, 1, 1, 1),
        dt(2000, 1, 1, 3),
    ])

    with pytest.raises(RuntimeError, match="non-uniform timesteps"):
        inspect_dates(tmp_path, dates)

def test_align_datetime_to_nearest_half_hour():
    assert align_datetime(dt(2000, 1, 1, 0, 5),
                          datetime.timedelta(minutes=30)) == dt(2000, 1, 1)
    assert align_datetime(dt(2000, 1, 1, 0, 20),
                          datetime.timedelta(minutes=30)) == dt(2000, 1, 1, 0, 30)

def test_align_datetime_can_round_to_next_day():
    assert align_datetime(dt(2000, 1, 1, 23, 50),
                          datetime.timedelta(hours=1)) == dt(2000, 1, 2)

def test_align_times_none_returns_input_values_unchanged():
    times = date_range(dt(2000, 1, 1, 0, 5), 2, datetime.timedelta(hours=1))

    aligned = align_times(times, AlignmentMode.NONE)

    assert aligned is times

def test_align_times_half_hour():
    times = numpy.array([dt(2000, 1, 1, 0, 5), dt(2000, 1, 1, 0, 35)])

    aligned = align_times(times, AlignmentMode.HALF_HOUR)

    assert_datetimes_equal(aligned, [dt(2000, 1, 1), dt(2000, 1, 1, 0, 30)])

def test_align_times_strict_hour():
    times = numpy.array([dt(2000, 1, 1, 0, 5), dt(2000, 1, 1, 1, 5)])

    aligned = align_times(times, AlignmentMode.STRICT_HOUR)

    assert_datetimes_equal(aligned, [dt(2000, 1, 1), dt(2000, 1, 1, 1)])

def test_align_times_rejects_invalid_mode():
    with pytest.raises(ValueError, match="Invalid alignment mode"):
        align_times(numpy.array([dt(2000, 1, 1)]), "bad")

def test_make_time_plan_without_alignment_preserves_input_phase(tmp_path):
    dates = date_range(dt(2000, 1, 1, 0, 5), 4, datetime.timedelta(minutes=30))
    policy = TimePolicy(output_timestep_minutes=60)

    plan = make_plan(tmp_path, dates, policy)

    assert plan.aggregation.aggregate_ratio == 2
    assert plan.output_n == 2
    assert_datetimes_equal(plan.output_time_values,
                           [dt(2000, 1, 1, 0, 5), dt(2000, 1, 1, 1, 5)])

def test_make_time_plan_half_hour_alignment_labels_output_grid(tmp_path):
    dates = date_range(dt(2000, 1, 1, 0, 5), 4, datetime.timedelta(minutes=15))
    policy = TimePolicy(output_timestep_minutes=30,
                        alignment_mode=AlignmentMode.HALF_HOUR)

    plan = make_plan(tmp_path, dates, policy)

    assert_datetimes_equal(plan.output_time_values,
                           [dt(2000, 1, 1), dt(2000, 1, 1, 0, 30)])

def test_make_time_plan_strict_hour_alignment_labels_output_grid(tmp_path):
    dates = date_range(dt(2000, 1, 1, 0, 5), 3, datetime.timedelta(hours=1))
    policy = TimePolicy(output_timestep_minutes=60,
                        alignment_mode=AlignmentMode.STRICT_HOUR)

    plan = make_plan(tmp_path, dates, policy)

    assert_datetimes_equal(plan.output_time_values,
                           [dt(2000, 1, 1), dt(2000, 1, 1, 1), dt(2000, 1, 1, 2)])

def test_make_time_plan_strict_hour_requires_at_least_hourly_output(tmp_path):
    dates = date_range(dt(2000, 1, 1, 0, 5), 3, datetime.timedelta(minutes=30))
    policy = TimePolicy(output_timestep_minutes=30,
                        alignment_mode=AlignmentMode.STRICT_HOUR)

    with pytest.raises(RuntimeError, match="STRICT_HOUR"):
        make_plan(tmp_path, dates, policy)

def test_make_time_plan_rejects_alignment_that_duplicates_timestamps(tmp_path):
    dates = date_range(dt(2000, 1, 1, 0, 5), 2, datetime.timedelta(minutes=10))
    policy = TimePolicy(output_timestep_minutes=10,
                        alignment_mode=AlignmentMode.HALF_HOUR)

    with pytest.raises(RuntimeError, match="duplicate timestamps"):
        make_plan(tmp_path, dates, policy)

def test_make_time_plan_rejects_alignment_that_changes_timestep_spacing(tmp_path):
    dates = date_range(dt(2000, 1, 1, 0, 5), 3, datetime.timedelta(minutes=20))
    policy = TimePolicy(output_timestep_minutes=20,
                        alignment_mode=AlignmentMode.HALF_HOUR)

    with pytest.raises(RuntimeError, match="non-uniform timesteps"):
        make_plan(tmp_path, dates, policy)

def test_make_time_plan_pads_one_missing_start_timestep(tmp_path):
    dates = date_range(dt(2000, 1, 1, 1, 5), 2, datetime.timedelta(hours=1))
    policy = TimePolicy(output_timestep_minutes=60)

    plan = make_plan(tmp_path, dates, policy)

    assert plan.aggregation.npad_start == 1
    assert_datetimes_equal(plan.output_time_values,
                           [dt(2000, 1, 1, 0, 5), dt(2000, 1, 1, 1, 5),
                            dt(2000, 1, 1, 2, 5)])

def test_make_time_plan_rejects_start_padding_beyond_max_pad(tmp_path):
    dates = date_range(dt(2000, 1, 1, 2, 5), 2, datetime.timedelta(hours=1))
    policy = TimePolicy(output_timestep_minutes=60)

    with pytest.raises(RuntimeError, match="maximum pad length"):
        make_plan(tmp_path, dates, policy)

def test_make_time_plan_start_boundary_error_mode_raises(tmp_path):
    dates = date_range(dt(2000, 1, 1, 1, 5), 2, datetime.timedelta(hours=1))
    policy = TimePolicy(output_timestep_minutes=60,
                        start_boundary_mode=BoundaryMode.ERROR)

    with pytest.raises(RuntimeError, match="missing timesteps at start"):
        make_plan(tmp_path, dates, policy)

def test_make_time_plan_trims_to_next_year_on_input_phase_grid(tmp_path):
    dates = date_range(dt(2000, 12, 31, 22, 5), 3, datetime.timedelta(hours=1))
    policy = TimePolicy(output_timestep_minutes=60,
                        start_boundary_mode=BoundaryMode.TRIM)

    plan = make_plan(tmp_path, dates, policy)

    assert plan.aggregation.ntrim_start == 2
    assert_datetimes_equal(plan.output_time_values, [dt(2001, 1, 1, 0, 5)])

def test_make_time_plan_pads_missing_end_timestep(tmp_path):
    dates = date_range(dt(2000, 12, 31, 21, 5), 2, datetime.timedelta(hours=1))
    policy = TimePolicy(output_timestep_minutes=60,
                        start_boundary_mode=BoundaryMode.IGNORE,
                        end_boundary_mode=BoundaryMode.PAD)

    plan = make_plan(tmp_path, dates, policy)

    assert plan.aggregation.npad_end == 1
    assert_datetimes_equal(plan.output_time_values,
                           [dt(2000, 12, 31, 21, 5),
                            dt(2000, 12, 31, 22, 5),
                            dt(2000, 12, 31, 23, 5)])

def test_make_time_plan_end_boundary_error_mode_raises(tmp_path):
    dates = date_range(dt(2000, 12, 31, 21, 5), 2, datetime.timedelta(hours=1))
    policy = TimePolicy(output_timestep_minutes=60,
                        start_boundary_mode=BoundaryMode.IGNORE,
                        end_boundary_mode=BoundaryMode.ERROR)

    with pytest.raises(RuntimeError, match="missing timesteps at end"):
        make_plan(tmp_path, dates, policy)

def test_make_time_plan_rejects_output_timestep_that_is_not_multiple_of_input(tmp_path):
    dates = date_range(dt(2000, 1, 1), 2, datetime.timedelta(hours=1))
    policy = TimePolicy(output_timestep_minutes=90)

    with pytest.raises(RuntimeError, match="not a multiple"):
        make_plan(tmp_path, dates, policy)

def test_make_time_plan_output_length_matches_applied_data_length(tmp_path):
    dates = date_range(dt(2000, 1, 1), 49, datetime.timedelta(minutes=30))
    policy = TimePolicy(output_timestep_minutes=60)
    plan = make_plan(tmp_path, dates, policy)
    values = numpy.arange(len(dates), dtype=float)

    out = apply_time_plan(values, plan.aggregation, NearestValuePadder(), numpy.mean)

    assert plan.output_n == 24
    assert len(out) == plan.output_n

def test_extrapolation_padder_extends_numeric_sequence():
    values = numpy.array([10, 20])

    out = ExtrapolationPadder().pad(values, npad_start=2, npad_end=1)

    assert list(out) == [-10, 0, 10, 20, 30]

def test_nearest_value_padder_repeats_edge_values():
    values = numpy.array([1, 2])

    out = NearestValuePadder().pad(values, npad_start=2, npad_end=1)

    assert list(out) == [1, 1, 1, 2, 2]

def test_value_padder_uses_fixed_value():
    values = numpy.array([1, 2])

    out = ValuePadder(0).pad(values, npad_start=1, npad_end=2)

    assert list(out) == [0, 1, 2, 0, 0]

def test_apply_time_plan_trims_pads_and_aggregates_values():
    values = numpy.array([10, 20, 30, 40, 50, 60], dtype=float)
    options = AggregationInfo(
        ntrim_start=1,
        ntrim_end=1,
        npad_start=1,
        npad_end=1,
        aggregate_ratio=2,
        timestep=7200,
    )

    out = apply_time_plan(values, options, NearestValuePadder(), numpy.mean)

    assert list(out) == [15, 35, 55]

def test_apply_time_plan_drops_incomplete_final_aggregate_window():
    values = numpy.array([1, 2, 3, 4, 999], dtype=float)
    options = AggregationInfo(
        ntrim_start=0,
        ntrim_end=0,
        npad_start=0,
        npad_end=0,
        aggregate_ratio=2,
        timestep=7200,
    )

    out = apply_time_plan(values, options, NearestValuePadder(), numpy.mean)

    assert list(out) == [1.5, 3.5]

def test_index_to_datetime_and_datetime_to_index_round_trip(tmp_path):
    dates = date_range(dt(2000, 1, 1), 3, datetime.timedelta(hours=1))
    policy = TimePolicy(output_timestep_minutes=60)
    plan = make_plan(tmp_path, dates, policy)

    when = index_to_datetime(plan, 1)

    assert when == dt(2000, 1, 1, 1)
    assert datetime_to_index(plan, when) == 1

def test_datetime_to_index_rejects_time_after_plan_end(tmp_path):
    dates = date_range(dt(2000, 1, 1), 3, datetime.timedelta(hours=1))
    policy = TimePolicy(output_timestep_minutes=60)
    plan = make_plan(tmp_path, dates, policy)

    with pytest.raises(ValueError, match="after end"):
        datetime_to_index(plan, dt(2000, 1, 1, 3))

def test_datetime_to_index_rejects_time_not_on_grid(tmp_path):
    dates = date_range(dt(2000, 1, 1), 3, datetime.timedelta(hours=1))
    policy = TimePolicy(output_timestep_minutes=60)
    plan = make_plan(tmp_path, dates, policy)

    with pytest.raises(ValueError, match="not on output time grid"):
        datetime_to_index(plan, dt(2000, 1, 1, 0, 30))

def test_validate_time_plan_accepts_valid_plan(tmp_path):
    dates = date_range(dt(2000, 1, 1), 3, datetime.timedelta(hours=1))
    policy = TimePolicy(output_timestep_minutes=60)
    plan = make_plan(tmp_path, dates, policy)

    validate_time_plan(plan)

def test_validate_time_plan_rejects_start_trim_and_pad(tmp_path):
    dates = date_range(dt(2000, 1, 1), 3, datetime.timedelta(hours=1))
    policy = TimePolicy(output_timestep_minutes=60)
    plan = make_plan(tmp_path, dates, policy)
    invalid = replace(
        plan,
        aggregation=replace(plan.aggregation, ntrim_start=1, npad_start=1),
    )

    with pytest.raises(RuntimeError, match="Cannot trim .* pad .* start"):
        validate_time_plan(invalid)

def test_validate_time_plan_rejects_end_trim_and_pad(tmp_path):
    dates = date_range(dt(2000, 1, 1), 3, datetime.timedelta(hours=1))
    policy = TimePolicy(output_timestep_minutes=60)
    plan = make_plan(tmp_path, dates, policy)
    invalid = replace(
        plan,
        aggregation=replace(plan.aggregation, ntrim_end=1, npad_end=1),
    )

    with pytest.raises(RuntimeError, match="Cannot trim .* pad .* end"):
        validate_time_plan(invalid)
