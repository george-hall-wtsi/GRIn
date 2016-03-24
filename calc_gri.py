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
import scipy.signal
import numpy as np


def find_start_repeat_kmers(hist_dict):

	order_num = 1

	min_list = scipy.signal.argrelextrema(np.array(hist_dict.values()), np.less_equal, 
		order = order_num)[0].tolist()
	max_list = scipy.signal.argrelextrema(np.array(hist_dict.values()), np.greater_equal, 
		order = order_num)[0].tolist()

	min_list_minimum = min(min_list)

	for maximum in sorted(max_list):
		if (maximum > min_list_minimum) and (maximum > 10):
			first_peak = maximum
			break

	# Return the point 'x' such that the first peak is equidistant 
	# between the first minimum and x
	return ((2 * first_peak) - min_list_minimum)


def create_hist_dict(in_file):
	hist_dict = {}
	for line in in_file.readlines():
		splat = [int(x) for x in line.strip().split()]
		hist_dict[splat[0]] = splat[1]

	return hist_dict


def calculate_gri(hist_dict, verbosity):
	start_repetitive_kmers = find_start_repeat_kmers(hist_dict)
	if verbosity != 0:
		print "Start of repetitive k-mers" , start_repetitive_kmers

	total_number_kmers = sum((a * b) for (a, b) in hist_dict.items())
	if verbosity != 0:
		print "Total number of k-mers" , total_number_kmers

	number_repetitive_kmers = 0
	for (a, b) in hist_dict.items():
		if (a >= start_repetitive_kmers):
			number_repetitive_kmers += (a * b)

	if verbosity != 0:
		print "Number of repetitive k-mers" , number_repetitive_kmers

	return ((1.0 * number_repetitive_kmers) / total_number_kmers)


def print_usage(argv):
	print "usage: %s [-v] <file>" %(argv[0])


if __name__ == "__main__":

	if not (2 <= len(sys.argv) <= 3):
		print "ERROR: Incorrect number of options"
		print_usage(sys.argv)
		sys.exit()

	if len(sys.argv) == 3:
		if "-v" == sys.argv[1]:
			verbosity = 1
		else:
			print "ERROR: Unrecognised option: %s" %(sys.argv[1])
			print_usage(sys.argv)
			sys.exit()
	else:
		verbosity = 0

	with open(sys.argv[-1], 'r') as f:
		hist_dict = create_hist_dict(f)
		gri = calculate_gri(hist_dict, verbosity)
		print "GRI =" , gri
