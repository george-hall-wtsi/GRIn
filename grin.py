#! /usr/bin/env python

################################################################################
# Copyright (c) 2016 Genome Research Ltd. 
# 
# Author: George Hall <gh10@sanger.ac.uk> 
# 
# This program is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation; either version 3 of the License, or (at your option) any later 
# version. 
# 
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
# details. 
# 
# You should have received a copy of the GNU General Public License along with 
# this program. If not, see <http://www.gnu.org/licenses/>. 
################################################################################


import sys
import subprocess

import scipy.signal
import numpy as np

import custom_argument_parser


def generate_min_list(hist_dict):
	
	min_list = scipy.signal.argrelextrema(np.array(hist_dict.values()), np.less_equal, 
		order = 1)[0].tolist()

	return min_list
	

def generate_max_list(hist_dict):

	max_list = scipy.signal.argrelextrema(np.array(hist_dict.values()), np.greater_equal, 
		order = 1)[0].tolist()

	return max_list


def find_first_peak_max(hist_dict, min_list = None):

	if not min_list:
		min_list = generate_min_list(hist_dict)

	max_list = generate_max_list(hist_dict)
	min_list_minimum = min(min_list)

	for maximum in sorted(max_list):
		if (maximum > min_list_minimum) and (maximum > 10):

			# Return the point 'x' such that the first peak is equidistant 
			# between the first minimum and x
			return maximum

	print "ERROR: Could not find the maximum of the first peak"

	sys.exit(1)


def find_start_first_peak(hist_dict):

	min_list = generate_min_list(hist_dict)
	first_peak_max = find_first_peak_max(hist_dict, min_list)

	for minimum in sorted(min_list)[::-1]:
		if minimum < first_peak_max:
			return minimum

	print "ERROR: Could not find the start of the first peak"

	sys.exit(1)


def find_start_repeat_kmers(hist_dict):

	# Return the point 'x' such that the first peak is equidistant between the first minimum
	# and x

	start_first_peak = find_start_first_peak(hist_dict)

	return ((2 * find_first_peak_max(hist_dict)) - start_first_peak)


def create_hist_dict(in_file):

	hist_dict = {}

	for line in in_file.readlines():
		splat = [float(x) for x in line.strip().split()]
		hist_dict[splat[0]] = splat[1]

	return hist_dict


def calculate_gri(hist_dict, verbose, error_cutoff, start_repetitive_kmers = 0):
	
	# Negative returns signify an error:
	# -1 => error cutoff greater than start of repetitive k-mers

	if not start_repetitive_kmers:
		if verbose:
			print "Estimating start of repetitive k-mers"
		start_repetitive_kmers = find_start_repeat_kmers(hist_dict)
	else:
		if verbose:
			print "User specified start of reptitive k-mers =" , start_repetitive_kmers

	if verbose:
		print "Start of repetitive k-mers" , start_repetitive_kmers

	# error_cutoff: 0 => use entire k-mer spectrum
	#				-1 => Auto error checking
	#				>= 1 => error cutoff manually specified
	if error_cutoff:
		if error_cutoff == -1:
			min_val_cutoff = find_start_first_peak(hist_dict)
		else:
			min_val_cutoff = error_cutoff

		if min_val_cutoff > start_repetitive_kmers:
			return -1

		if verbose:
			print "Using minimum k-mer occurrence of" , min_val_cutoff
	else:
		min_val_cutoff = 0

	total_number_kmers = sum((a * b) for (a, b) in hist_dict.items() if \
			((not error_cutoff) or (a > min_val_cutoff)))

	if verbose:
		print "Total number of k-mers" , total_number_kmers

	number_repetitive_kmers = 0
	for (a, b) in hist_dict.items():
		if (a >= start_repetitive_kmers):
			number_repetitive_kmers += (a * b)

	if verbose:
		print "Number of repetitive k-mers" , number_repetitive_kmers

	return ((1.0 * number_repetitive_kmers) / total_number_kmers)


def create_parser():
	parser = custom_argument_parser.CustomParser()
	parser.add_argument("-v", "--verbose", action = "store_true")
	parser.add_argument("-c", "--repeat-cutoffs", type = int, nargs = '+')
	parser.add_argument("-a", "--analyzer", action = "store_true")
	parser.add_argument("-e", "--manual-error-cutoffs", type = int, nargs = '+')
	parser.add_argument("-i", "--ignore-error", action = "store_true")
	parser.add_argument("-f", "--file", type = str, nargs = '+', required = True)

	return parser


def parser_main():
	parser = create_parser()
	args = parser.parse_args()

	if args.repeat_cutoffs:
		if len(args.file) != len(args.repeat_cutoffs):
			print "ERROR: Need to have the same number of manual repeat cutoffs as files"
			sys.exit(1)
	else:
		args.repeat_cutoffs = [0 for x in args.file]

	if args.manual_error_cutoffs and args.ignore_error:
		print "ERROR: Cannot specify both --manual-error-cutoffs and --ignore-error"
		sys.exit(1)

	if args.manual_error_cutoffs:
		if len(args.file) != len(args.manual_error_cutoffs):
			print "ERROR: Need to have the same number of manual error cutoffs as files"
			sys.exit(1)

		if any(cutoff <= 0 for cutoff in args.manual_error_cutoffs):
			print "ERROR: --manual-error-cuttoffs must be positive"
			sys.exit(1)

	return args


def set_error_cutoffs(ignore_error, manual_error_cutoffs, file_list):
	error_cutoffs = []
	if not ignore_error and not manual_error_cutoffs:
		# User doesn't want to do anything about errors:
		error_cutoffs = [0 for x in file_list]
	elif ignore_error:
		# User wants errror cutoff to be automatically determined
		error_cutoffs = [-1 for x  in file_list]
	elif manual_error_cutoffs:
		# User has manually specified error cutoffs
		error_cutoffs = manual_error_cutoffs
	else:
		# ERROR: Should have been caught in parser_main()
		print "ERROR: Option parsing error checking let through some mutually exclusive options"
		sys.exit(1)

	return error_cutoffs


def main():
	args = parser_main()

	manual_repeat_cutoffs = args.repeat_cutoffs
	file_paths = args.file
	verbose = args.verbose
	analyzer = args.analyzer
	
	error_cutoffs = set_error_cutoffs(args.ignore_error, args.manual_error_cutoffs, args.file)

	for (file_name, repeat_cutoff, error_cutoff) in \
	zip(file_paths, manual_repeat_cutoffs, error_cutoffs):

		if analyzer:
			subprocess.call(["kmerspectrumanalyzer", file_name])
			with open(file_name + ".fit.detail.csv", 'r') as f, \
			open(file_name + ".fit.detail.csv.hist", 'w') as g:
					for line in f.readlines():
						splat = line.strip().split()
						g.write(splat[0] + " " + splat[2] + "\n")
				
			file_name += ".fit.detail.csv.hist"

		with open(file_name, 'r') as f:
			print "Started processing" , file_name
			hist_dict = create_hist_dict(f)
			gri = calculate_gri(hist_dict, verbose, error_cutoff, repeat_cutoff)
			if gri == -1:
				print "ERROR: Error cutoff greater than start of repetitive k-mers. Skipping this file..."
			else:
				print "GRI = %0.4f" %(gri)
			print "Finished processing" , file_name


if __name__ == "__main__":
	main()

