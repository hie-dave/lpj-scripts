# Logging utility functions.

import enum, sys
from io import TextIOWrapper

MAX_PROGRESS_MSG_LEN = 17

class LogLevel(enum.IntEnum):
	NONE = 0,
	ERROR = 1,
	WARNING = 2,
	INFORMATION = 3,
	DIAGNOSTIC = 4,
	DEBUG = 5

# Configuration options, which can be set via one of the below
# functions. Best not to set these directly from user code.
_log_level: LogLevel
_show_progress: bool
_warnings_as_errors: bool

def set_warnings_as_errors(warnings_as_errors: bool):
	"""
	Treat all warnings as errors (true to enable, false to disable).
	"""
	global _warnings_as_errors
	_warnings_as_errors = warnings_as_errors

def set_log_level(log_level: LogLevel):
	"""
	Set the log level, which controls which messages are written to the
	log.
	"""
	global _log_level
	_log_level = log_level

def set_show_progress(show_progress: bool):
	"""
	Set whether progress messages will be displayed.
	"""
	global _show_progress
	_show_progress = show_progress

def clear_line(file: TextIOWrapper):
	"""
	Clear a line by writing empty spaces followed by a carriage return,
	in order to remove any progress report written to this line.
	"""
	print(" " * MAX_PROGRESS_MSG_LEN, file = file, end = "\r")

def log(msg: str, log_level: LogLevel):
	"""
	Write a log message.
	"""
	# I'm putting this check here, rather than in log_warning(), in case
	# this function is called directly from user code.
	if log_level == LogLevel.WARNING and _warnings_as_errors:
		raise ValueError(msg)

	# todo: custom log file as CLI arg?
	if log_level <= _log_level:
		file = sys.stderr if log_level == LogLevel.ERROR else sys.stdout
		clear_line(file = file)
		print(msg, file = file)

def log_error(msg: str):
	"""
	Write an error message.
	"""
	log(msg, LogLevel.ERROR)

def log_warning(msg: str):
	"""
	Write a warning message.
	"""
	log(msg, LogLevel.WARNING)

def log_information(msg: str):
	"""
	Write an information message.
	"""
	log(msg, LogLevel.INFORMATION)

def log_diagnostic(msg: str):
	"""
	Write a diagnostic message.
	"""
	log(msg, LogLevel.DIAGNOSTIC)

def log_debug(msg: str):
	"""
	Write a debug message.
	"""
	log(msg, LogLevel.DEBUG)

def log_progress(progress: float):
	"""
	Write a progress message.
	"""
	global _show_progress
	if _show_progress:
		print("Working: %.2f%%\r" % (100 * progress), end = "")
