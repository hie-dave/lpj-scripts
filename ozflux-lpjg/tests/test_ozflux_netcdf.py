import sys
from pathlib import Path

import numpy
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import ozflux_netcdf  # noqa: E402
from ozflux_netcdf import _filter_qc, bounds_checks, remove_nans  # noqa: E402

def masked_array(values, mask):
	return numpy.ma.array(values, mask=mask, dtype=float)

def assert_array_equal_with_nan(actual, expected):
	numpy.testing.assert_allclose(actual, expected, equal_nan=True)

def assert_mask_equal(actual, expected):
	numpy.testing.assert_array_equal(numpy.ma.getmaskarray(actual), expected)

def test_remove_nans_returns_float_ndarray_for_plain_array():
	data = numpy.array([1, 2, 3], dtype=int)

	result = remove_nans(data)

	assert isinstance(result, numpy.ndarray)
	assert result.dtype == float
	assert_array_equal_with_nan(result, [1.0, 2.0, 3.0])

def test_remove_nans_returns_unmasked_values_unchanged():
	data = masked_array([1.5, 2.5, 3.5], [False, False, False])

	result = remove_nans(data)

	assert_array_equal_with_nan(result, [1.5, 2.5, 3.5])

def test_remove_nans_interpolates_single_masked_value_between_valid_values():
	data = masked_array([1.0, 999.0, 5.0], [False, True, False])

	result = remove_nans(data)

	assert_array_equal_with_nan(result, [1.0, 3.0, 5.0])

def test_remove_nans_interpolates_consecutive_masked_values():
	data = masked_array([0.0, 999.0, 999.0, 6.0], [False, True, True, False])

	result = remove_nans(data)

	assert_array_equal_with_nan(result, [0.0, 2.0, 4.0, 6.0])

def test_remove_nans_fills_leading_masked_values_with_first_valid_value():
	data = masked_array([999.0, 999.0, 3.0, 6.0], [True, True, False, False])

	result = remove_nans(data)

	assert_array_equal_with_nan(result, [3.0, 3.0, 3.0, 6.0])

def test_remove_nans_fills_trailing_masked_values_with_last_valid_value():
	data = masked_array([1.0, 4.0, 999.0, 999.0], [False, False, True, True])

	result = remove_nans(data)

	assert_array_equal_with_nan(result, [1.0, 4.0, 4.0, 4.0])

def test_remove_nans_handles_masked_values_at_both_edges():
	data = masked_array(
		[999.0, 2.0, 999.0, 8.0, 999.0],
		[True, False, True, False, True],
	)

	result = remove_nans(data)

	assert_array_equal_with_nan(result, [2.0, 2.0, 5.0, 8.0, 8.0])

def test_remove_nans_does_not_replace_unmasked_nan_values():
	data = masked_array([0.0, numpy.nan, 999.0, 6.0], [False, False, True, False])

	result = remove_nans(data)

	assert_array_equal_with_nan(result, [0.0, numpy.nan, 4.0, 6.0])

def test_remove_nans_does_not_use_unmasked_nan_as_interpolation_support():
	data = masked_array([0.0, numpy.nan, 999.0, 12.0], [False, False, True, False])

	result = remove_nans(data)

	assert_array_equal_with_nan(result, [0.0, numpy.nan, 8.0, 12.0])

def test_remove_nans_does_not_mutate_input_masked_array():
	data = masked_array([1.0, 999.0, 5.0], [False, True, False])
	original_values = data.data.copy()
	original_mask = data.mask.copy()

	remove_nans(data)

	assert_array_equal_with_nan(data.data, original_values)
	numpy.testing.assert_array_equal(data.mask, original_mask)

def test_remove_nans_raises_when_all_values_are_masked():
	data = masked_array([999.0, 999.0], [True, True])

	with pytest.raises(RuntimeError, match="no valid values"):
		remove_nans(data)

def test_remove_nans_raises_when_masked_values_have_no_finite_support():
	data = masked_array([numpy.nan, 999.0], [False, True])

	with pytest.raises(RuntimeError, match="no valid values"):
		remove_nans(data)

def test_filter_qc_masks_values_with_failing_qc_flags():
	data = numpy.ma.array([1.0, 2.0, 3.0, 4.0])
	qc = numpy.array([0, 1, 10, 50])

	result = _filter_qc(data, qc)

	numpy.testing.assert_array_equal(result.data, [1.0, 2.0, 3.0, 4.0])
	assert_mask_equal(result, [False, True, False, True])

def test_filter_qc_keeps_values_with_passing_qc_flags():
	data = numpy.ma.array([1.0, 2.0, 3.0, 4.0])
	qc = numpy.array([0, 10, 15, 30])

	result = _filter_qc(data, qc)

	assert_mask_equal(result, [False, False, False, False])

def test_filter_qc_keeps_values_with_unknown_qc_flags():
	data = numpy.ma.array([1.0, 2.0, 3.0])
	qc = numpy.array([0, 999, -1])

	result = _filter_qc(data, qc)

	assert_mask_equal(result, [False, False, False])

def test_filter_qc_preserves_existing_mask():
	data = masked_array([1.0, 2.0, 3.0, 4.0], [False, True, False, False])
	qc = numpy.array([0, 0, 1, 0])

	result = _filter_qc(data, qc)

	assert_mask_equal(result, [False, True, True, False])

def test_filter_qc_accepts_plain_ndarray_data():
	data = numpy.array([1.0, 2.0, 3.0])
	qc = numpy.array([0, 2, 0])

	result = _filter_qc(data, qc)

	assert isinstance(result, numpy.ma.MaskedArray)
	numpy.testing.assert_array_equal(result.data, [1.0, 2.0, 3.0])
	assert_mask_equal(result, [False, True, False])

def test_filter_qc_casts_float_qc_flags_to_int():
	data = numpy.ma.array([1.0, 2.0, 3.0])
	qc = numpy.array([0.0, 1.0, 10.0])

	result = _filter_qc(data, qc)

	assert_mask_equal(result, [False, True, False])

def test_filter_qc_raises_when_lengths_differ():
	data = numpy.ma.array([1.0, 2.0, 3.0])
	qc = numpy.array([0, 1])

	with pytest.raises(RuntimeError, match="data array has length 3 but QC array has length 2"):
		_filter_qc(data, qc)

def test_filter_qc_logs_one_summary_per_failing_qc_flag(monkeypatch):
	messages = []
	monkeypatch.setattr(ozflux_netcdf, "log_debug", messages.append)

	data = numpy.ma.array([1.0, 2.0, 3.0, 4.0, 5.0])
	qc = numpy.array([1, 1, 2, 10, 999])

	_filter_qc(data, qc)

	assert len(messages) == 2
	assert any("Filtering 2 values with QC flag 1" in msg for msg in messages)
	assert any("Filtering 1 values with QC flag 2" in msg for msg in messages)

def test_filter_qc_does_not_log_when_no_values_are_filtered(monkeypatch):
	messages = []
	monkeypatch.setattr(ozflux_netcdf, "log_debug", messages.append)

	data = numpy.ma.array([1.0, 2.0, 3.0])
	qc = numpy.array([0, 10, 999])

	_filter_qc(data, qc)

	assert messages == []

def test_bounds_checks_returns_passed_true_when_all_values_are_in_bounds(monkeypatch):
	messages = []
	monkeypatch.setattr(ozflux_netcdf, "log_diagnostic", messages.append)
	data = numpy.array([1.0, 2.0, 3.0])

	result, passed = bounds_checks(data, 0.0, 4.0)

	assert passed is True
	assert result is data
	numpy.testing.assert_array_equal(result, [1.0, 2.0, 3.0])
	assert messages == []

def test_bounds_checks_treats_values_equal_to_bounds_as_passing(monkeypatch):
	messages = []
	monkeypatch.setattr(ozflux_netcdf, "log_diagnostic", messages.append)
	data = numpy.array([0.0, 1.0, 4.0])

	result, passed = bounds_checks(data, 0.0, 4.0)

	assert passed is True
	numpy.testing.assert_array_equal(result, [0.0, 1.0, 4.0])
	assert messages == []

def test_bounds_checks_clips_lower_bound_exceedances_and_logs_summary(monkeypatch):
	messages = []
	monkeypatch.setattr(ozflux_netcdf, "log_diagnostic", messages.append)
	data = numpy.array([-2.0, 0.0, -1.5, 3.0])

	result, passed = bounds_checks(data, -1.0, 5.0)

	assert passed is False
	assert result is data
	numpy.testing.assert_array_equal(result, [-1.0, 0.0, -1.0, 3.0])
	assert messages == [
		"2 values exceed lower bound of -1.00; "
		"first 2 row/value pair(s): 0: -2.00, 2: -1.50"
	]

def test_bounds_checks_clips_upper_bound_exceedances_and_logs_summary(monkeypatch):
	messages = []
	monkeypatch.setattr(ozflux_netcdf, "log_diagnostic", messages.append)
	data = numpy.array([1.0, 8.0, 3.0, 9.5])

	result, passed = bounds_checks(data, 0.0, 4.0)

	assert passed is False
	assert result is data
	numpy.testing.assert_array_equal(result, [1.0, 4.0, 3.0, 4.0])
	assert messages == [
		"2 values exceed upper bound of 4.00; "
		"first 2 row/value pair(s): 1: 8.00, 3: 9.50"
	]

def test_bounds_checks_clips_and_logs_lower_and_upper_bound_exceedances(monkeypatch):
	messages = []
	monkeypatch.setattr(ozflux_netcdf, "log_diagnostic", messages.append)
	data = numpy.array([-5.0, 0.0, 10.0, 2.0, -2.0, 7.0])

	result, passed = bounds_checks(data, -1.0, 5.0)

	assert passed is False
	numpy.testing.assert_array_equal(result, [-1.0, 0.0, 5.0, 2.0, -1.0, 5.0])
	assert messages == [
		"2 values exceed lower bound of -1.00; "
		"first 2 row/value pair(s): 0: -5.00, 4: -2.00",
		"2 values exceed upper bound of 5.00; "
		"first 2 row/value pair(s): 2: 10.00, 5: 7.00",
	]

def test_bounds_checks_logs_only_first_configured_number_of_indices(monkeypatch):
	messages = []
	monkeypatch.setattr(ozflux_netcdf, "log_diagnostic", messages.append)
	monkeypatch.setattr(ozflux_netcdf, "MAX_BOUNDS_LOG_INDICES", 3)
	data = numpy.array([-1.0, -2.0, -3.0, -4.0, 10.0, 11.0, 12.0, 13.0])

	result, passed = bounds_checks(data, 0.0, 5.0)

	assert passed is False
	numpy.testing.assert_array_equal(result, [0.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0, 5.0])
	assert messages == [
		"4 values exceed lower bound of 0.00; "
		"first 3 row/value pair(s): 0: -1.00, 1: -2.00, 2: -3.00, ... (1 more)",
		"4 values exceed upper bound of 5.00; "
		"first 3 row/value pair(s): 4: 10.00, 5: 11.00, 6: 12.00, ... (1 more)",
	]

def test_bounds_checks_accepts_plain_list_and_returns_new_array(monkeypatch):
	messages = []
	monkeypatch.setattr(ozflux_netcdf, "log_diagnostic", messages.append)
	data = [-2.0, 1.0, 6.0]

	result, passed = bounds_checks(data, 0.0, 5.0)

	assert passed is False
	assert isinstance(result, numpy.ndarray)
	numpy.testing.assert_array_equal(result, [0.0, 1.0, 5.0])
	assert data == [-2.0, 1.0, 6.0]

def test_bounds_checks_preserves_nan_values_as_passing(monkeypatch):
	messages = []
	monkeypatch.setattr(ozflux_netcdf, "log_diagnostic", messages.append)
	data = numpy.array([1.0, numpy.nan, 3.0])

	result, passed = bounds_checks(data, 0.0, 4.0)

	assert passed is True
	assert_array_equal_with_nan(result, [1.0, numpy.nan, 3.0])
	assert messages == []

def test_bounds_checks_handles_empty_array(monkeypatch):
	messages = []
	monkeypatch.setattr(ozflux_netcdf, "log_diagnostic", messages.append)
	data = numpy.array([])

	result, passed = bounds_checks(data, 0.0, 4.0)

	assert passed is True
	assert result is data
	numpy.testing.assert_array_equal(result, [])
	assert messages == []
