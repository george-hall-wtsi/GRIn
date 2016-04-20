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


from __future__ import print_function, division
import sys
import subprocess

import scipy.signal
import numpy as np

import custom_argument_parser


def generate_min_list(hist_dict):

    """
    Returns a list of the x values corresponding to the minima in hist_dict.
    """
    
    min_list = scipy.signal.argrelextrema(np.array(hist_dict.values()), 
        np.less_equal, order = 3)[0].tolist()

    return min_list
	

def generate_max_list(hist_dict):

    """
    Returns a list of the x values corresponding to the maxima in hist_dict.
    """

    max_list = scipy.signal.argrelextrema(np.array(hist_dict.values()), 
        np.greater_equal, order = 3)[0].tolist()

    return max_list


def find_main_peak(hist_dict, min_list = None):

    """
    Returns the x value of the main peak in hist_dict. The main peak is the 
    peak which should correspond to homozygous k-mers, not that which occurs 
    right at the beginning of the k-mer spectrum, which is due to base errors 
    in the reads.
    """

    min_list = generate_min_list(hist_dict)
    max_list = generate_max_list(hist_dict)
    min_list_minimum = min(min_list)

    max_y_vals = [hist_dict[x] for x  in max_list if x in hist_dict.keys()]

    for maximum in sorted(max_y_vals)[::-1]:
        # Although this looks really inefficient, we hopefully
        # shouldn't have to do it more than twice (hopefully just once)
        # and I can't quickly think of a better way to do it
        for (k, v) in hist_dict.iteritems():
            if v == maximum and k > min_list_minimum and k > 10:
                return k

    print("ERROR: Could not find the maximum of the main peak", 
        file = sys.stderr)
    sys.exit(1)


def find_start_main_peak(hist_dict):

    """
    Returns the smallest minimum value. 
    """

    return min(generate_min_list(hist_dict))


def find_start_repeat_kmers(hist_dict, verbose = False):

    """
    Returns a point 'a' such that the x co-ordinate of the main peak is
    equidistant between the start of the main peak and 'a'. 
    """

    start_first_peak = find_start_main_peak(hist_dict)
    if verbose:
        print("Start of first peak =", start_first_peak)

    return int((2 * find_main_peak(hist_dict)) - start_first_peak)


def create_hist_dict(in_file):

    """
    Returns a dictionary with number of occurrences as the keys and the
    frequency corresponding to each occurence as its value. 
    """

    hist_dict = {}

    for line in in_file.readlines():
        splat = [float(x) for x in line.strip().split()]
        hist_dict[splat[0]] = splat[1]

    return hist_dict


def calculate_gri(hist_dict, verbose, error_cutoff, upper_bound, 
    start_repetitive_kmers = 0):

    """
    Returns the GRI, which we have defined to be the percentage of
    repetitive k-mers in a k-mer spectrum. The user can choose to ignore
    k-mers caused by base errors when counting the total number of k-mers,
    thus increasing the GRI. 
    """
    
    # Negative returns signify an error:
    # -1 => error cutoff greater than start of repetitive k-mers

    if not start_repetitive_kmers:
        if verbose:
            print("Estimating start of repetitive k-mers")
        start_repetitive_kmers = find_start_repeat_kmers(hist_dict, verbose)
    else:
        if verbose:
            print("User specified start of reptitive k-mers =", 
				start_repetitive_kmers)

    if verbose:
        print("Start of repetitive k-mers", start_repetitive_kmers)

    # error_cutoff: 0 => use entire k-mer spectrum
    #				-1 => Auto error checking
    #				>= 1 => error cutoff manually specified
    if error_cutoff:
        if error_cutoff == -1:
            min_val_cutoff = find_start_main_peak(hist_dict)
        else:
            min_val_cutoff = error_cutoff

        if min_val_cutoff > start_repetitive_kmers:
            return -1

        if verbose:
            print("Using minimum k-mer occurrence of", min_val_cutoff)
    else:
        min_val_cutoff = 0

    if upper_bound is None:
        upper_bound = max(hist_dict.keys())
    else:
        if verbose:
            print("Using upper bound of", upper_bound)

    total_number_kmers = sum((a * b) for (a, b) in hist_dict.items() if \
        ((a <= upper_bound) and ((not error_cutoff) or (a > min_val_cutoff))))

    if verbose:
        print("Total number of k-mers", total_number_kmers)

    number_repetitive_kmers = 0
    for (a, b) in hist_dict.items():
        if ((a <= upper_bound) and (a >= start_repetitive_kmers)):
            number_repetitive_kmers += (a * b)

    if verbose:
        print("Number of repetitive k-mers", number_repetitive_kmers)

    return (number_repetitive_kmers / total_number_kmers)


def create_parser():
	
    """
    Returns a parser based on my custom parser, which is stored in 
    custom_argument_parser.py. I have done it in this way because I wanted to 
    be able to write my own help and usage messages. 
    """

    parser = custom_argument_parser.CustomParser()
    parser.add_argument("-v", "--verbose", action = "store_true")
    parser.add_argument("-c", "--repeat-cutoffs", type = int, nargs = '+')
    parser.add_argument("-a", "--analyzer", action = "store_true")
    parser.add_argument("-e", "--manual-error-cutoffs", type = int,
        nargs = '+')
    parser.add_argument("-E", "--single-error-cutoff", type = int, nargs = '?')
    parser.add_argument("-i", "--ignore-error", action = "store_true")
    parser.add_argument("-u", "--upper-bound", type = int)
    parser.add_argument("-f", "--file", type = str, nargs = '+',
        required = True)

    return parser


def parser_main():

    """
    Parses command line arguments and performs extra error checking.
    Returns these arguments as a Namespace object. 
    """

    parser = create_parser()
    args = parser.parse_args()

    if args.repeat_cutoffs:
        if len(args.file) != len(args.repeat_cutoffs):
            print("ERROR: Need to have the same number of manual repeat",
                "cutoffs as files", file = sys.stderr)
            sys.exit(1)
    else:
        args.repeat_cutoffs = [0 for _ in args.file]

    if args.manual_error_cutoffs and args.ignore_error:
        print("ERROR: Cannot specify both --manual-error-cutoffs and",
            "--ignore-error", file = sys.stderr)
        sys.exit(1)

    if args.manual_error_cutoffs and args.single_error_cutoff:
        print("ERROR: Cannot specify both --manual-error cutoffs and",
            "--single-error-cutoff", file = sys.stderr)
        sys.exit(1)

    if args.ignore_error and args.single_error_cutoff:
        print("ERROR: Cannot specify both --ignore-error and",
            "--single-error-cutoff", file = sys.stderr)
        sys.exit(1)

    if args.manual_error_cutoffs:
        if len(args.file) != len(args.manual_error_cutoffs):
            print("ERROR: Need to have the same number of manual error",
                "cutoffs as files", file = sys.stderr)
            sys.exit(1)

        if any(cutoff <= 0 for cutoff in args.manual_error_cutoffs):
            print("ERROR: --manual-error-cuttoffs must be a positive integer", 
                file = sys.stderr)
            sys.exit(1)

    if args.single_error_cutoff:
        if args.single_error_cutoff <= 0:
            print("ERROR: --single-error-cutoff must be a positive integer", 
                file = sys.stderr)
            sys.exit(1)

    if args.upper_bound:
        if args.upper_bound <= 0:
            print("ERROR: --upper-bound must be a positive integer",
                file = sys.stderr)
            sys.exit(1)

    return args


def set_error_cutoffs(ignore_error, manual_error_cutoffs, single_error_cutoff, 
        file_list):

    """
    Determines if the user wants to ignore the k-mers attributed to base errors
    or not, and then sets the error cutoffs accordingly, such that
    calculate_gri() behaves correctly. 
    """ # Check that only one of the three options has been set:
    assert 0 <= sum([bool(x) for x in 
        [ignore_error, manual_error_cutoffs, single_error_cutoff]]) <= 1, \
        "Can only set one of ignore_error, manual_error_cutoffs, " + \
        "single_error_cutoff"

    error_cutoffs = []
    if not single_error_cutoff and \
        not ignore_error and \
        not manual_error_cutoffs:

        # User doesn't want to do anything about errors:
        error_cutoffs = [0 for x in file_list]

    elif ignore_error:
        # User wants errror cutoff to be automatically determined
        error_cutoffs = [-1 for x  in file_list]
    elif single_error_cutoff or manual_error_cutoffs:
        # User has manually specified error cutoffs
        if manual_error_cutoffs:
            error_cutoffs = manual_error_cutoffs
        else:
            error_cutoffs = [single_error_cutoff for x in file_list]

    return error_cutoffs


def run_analyzer(file_name):

    """
    First run KMERSPECTRUMANALYZER (see README for citation) and then calculate 
    GRI based on the histogram output by this program. This idea is to give a 
    k-mer spectrum with more noticable and accurate repeat peaks, but this is 
    often not the case.
    """

    subprocess.call(["kmerspectrumanalyzer", file_name])
    ksa_output_file = file_name + ".fit.detail.csv"
    with open(ksa_output_file, 'r') as f, \
        open(ksa_output_file + ".hist", 'w') as g:
        
        for line in f.readlines():
            splat = line.strip().split()
            g.write(splat[0] + " " + splat[2] + "\n")


def process_histogram_file(file_name, verbose, error_cutoff, upper_bound, 
        repeat_cutoff):

    """
    Compute GRI for a histogram file.
    """

    with open(file_name, 'r') as f:
        print("Processing", file_name)
        hist_dict = create_hist_dict(f)
        gri = calculate_gri(hist_dict, verbose, error_cutoff, upper_bound, 
            repeat_cutoff)
        if gri == -1:
            print("ERROR: Error cutoff greater than start of",
                "repetitive k-mers. Skipping this file...", file = sys.stderr)
        else:
            print("GRI = %0.4f" %(gri))


def main():

    """
    Main driver function. 
    """

    args = parser_main()

    manual_repeat_cutoffs = args.repeat_cutoffs
    file_paths = args.file
    verbose = args.verbose
    analyzer = args.analyzer
    upper_bound = args.upper_bound
    
    error_cutoffs = set_error_cutoffs(args.ignore_error, 
        args.manual_error_cutoffs, args.single_error_cutoff, args.file)

    for (file_name, repeat_cutoff, error_cutoff) in \
    zip(file_paths, manual_repeat_cutoffs, error_cutoffs):

        if analyzer:
            run_analyzer(file_name)	
            file_name += ".fit.detail.csv.hist"

        try:
            process_histogram_file(file_name, verbose, error_cutoff, 
                upper_bound, repeat_cutoff)
        except IOError:
            print("ERROR: Could not open file \"" + file_name + "\".", 
                "Skipping...", file = sys.stderr)


if __name__ == "__main__":
    main()

