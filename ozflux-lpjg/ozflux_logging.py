# Logging functions
import enum, sys

class LogLevel(enum.IntEnum):
	NONE = 0,
	ERROR = 1,
	WARNING = 2,
	INFORMATION = 3,
	DIAGNOSTIC = 4,
	DEBUG = 5

global_log_level: LogLevel

def log(msg: str, log_level: LogLevel):
	"""
	Write a log message.
	"""
	# todo: custom log file as CLI arg?
	if log_level <= global_log_level:
		file = sys.stderr if log_level == LogLevel.ERROR else sys.stdout
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

def dump_netcdf_file(filename):
	"""
	Print all variable names/descriptions in the .nc file.
	"""
	with open(filename) as nc_file:
		for (key, value) in nc_file.variables.items():
			print("%s: %s" % (key, value.long_name))
