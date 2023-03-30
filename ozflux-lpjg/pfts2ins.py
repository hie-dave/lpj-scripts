#!/usr/bin/env python
#
# This script translates an excel spreadsheet of PFTs into an LPJ-Guess
# instruciton file.
from argparse import ArgumentParser
from ozflux_common import *
from ozflux_logging import *
import ozflux_logging
import numpy, os, pandas, sys, traceback

COL_PARAMETER = "Parameter"
COL_DESCRIPTION = "Description"
COL_UNITS = "Unit"
COL_GLOBAL = "Global"

GROUP_GLOBAL = "global"

class _Options:
	def __init__(self, verbosity: int, infile: str, outfile: str, sheet: str
	      , global_sheet: str):
		self.verbosity = verbosity
		self.infile = infile
		self.outfile = outfile
		self.sheet_name = sheet
		self.global_sheet = global_sheet

def parse_args(argv: list[str]) -> _Options:
	"""
	Parse CLI arguments, return a parsed options object.

	@param argv: Raw CLI arguments.

	@return Parsed Options object.
	"""
	parser = ArgumentParser(prog=argv[0], description = "Formatting ozflux data into a format suitable for consumption by LPJ-Guess")
	parser.add_argument("-v", "--verbosity", type = int, help = "Logging verbosity (1-5, default 3)", nargs = "?", const = LogLevel.INFORMATION, default = LogLevel.INFORMATION)
	parser.add_argument("-i", "--in-file", required = True, help = "Path to input excel file")
	parser.add_argument("-o", "--out-file", required = True, help = "Name/path of output file")
	parser.add_argument("-s", "--sheet-name", required = True, help = "Name of the sheet to read in the input file")
	parser.add_argument("-g", "--global-sheet", required = True, help = "Name of the sheet in the input file containing global parameters")
	parser.add_argument("--version", action = "version", version = "%(prog)s " + VERSION)
	parsed = parser.parse_args(argv[1:])

	return _Options(parsed.verbosity, parsed.in_file, parsed.out_file
		, parsed.sheet_name, parsed.global_sheet)

def is_float(s: str) -> bool:
	try:
		f = float(s)
		return True
	except:
		return False

def main(opts: _Options):
	"""
	Main CLI entrypoint function.

	@param opts: Object containing parsed CLI arguments.
	"""
	# Create output directory if it doesn't exist.
	out_dir = os.path.dirname(opts.outfile)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	# Read input data.
	data = pandas.read_excel(opts.infile, opts.sheet_name)
	global_data = pandas.read_excel(opts.infile, opts.global_sheet)

	metadata_cols = [COL_PARAMETER, COL_DESCRIPTION, COL_UNITS]
	# group_names = ["leafphysiognomy", "phenology", ]
	row_names: list[str] = [x for x in data[COL_PARAMETER].iloc]
	with open(opts.outfile, "w") as file:
		file.writelines('group "%s" (\n' % GROUP_GLOBAL)
		global_param_names = [x for x in global_data[COL_PARAMETER].iloc]
		global_param_values = [x for x in global_data[COL_GLOBAL].iloc]
		for (name, value) in zip(global_param_names, global_param_values):
			file.write("  %s %s\n" % (name, value))
		file.writelines(")\n\n")

		for col in data:
			if col in metadata_cols:
				continue
			pft = col

			file.writelines('pft "%s" (\n' % pft)
			file.writelines("  %s\n" % GROUP_GLOBAL)
			file.writelines("  include 1\n")
			file.writelines('  landcover "natural"\n')
			i = 0
			for row in data[col].iloc:
				# if row_names[i] in group_names:
				# 	file.writelines("  %s\n" % row)
				# else:
				if row == "nan" or (is_float(row) and math.isnan(float(row))):
					print("WARNING: NaN value for PFT '%s' parameter '%s'" % (pft, row_names[i]))
				elif not is_float(row) and not all(is_float(x) for x in str(row).split(" ")) and str(row) != "nan":
					file.writelines('  %s "%s"\n' % (row_names[i], row))
				else:
					if row == "-":
						row = "0"
					if all(c.isnumeric() or c == "e" or c == "E" or c == "-" or c == "." for c in str(row)):
						file.writelines("  %s %.8f\n" % (row_names[i], float(row)))
					else:
						file.writelines("  %s %s\n" % (row_names[i], row))
				i += 1
			file.writelines(")\n\n")

if __name__ == "__main__":
	# Parse CLI args
	opts = parse_args(sys.argv)

	set_log_level(opts.verbosity)
	# set_show_progress(opts.report_progress)

	try:
		# Actual logic is in main().
		main(opts)
	except BaseException as error:
		# Basic error handling.
		log_error(traceback.format_exc())
		exit(1)
