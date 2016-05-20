#! /usr/bin/env python

###############################################################################
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
###############################################################################


"""
GRIn is a program which computed a set of read's Genome Repeat Index (GRI). The
GRI is defined as being the proportion of k-mers in a set of read's k-mer
spectrum which are repetitive. All of this is explained in much more detail in
the paper accompanying this program (still being written as of 16/5/16). For
anything to do with this program, contact George Hall (gh10@sanger.ac.uk).
"""


from __future__ import print_function, division
import sys

import scipy.signal as sig
import numpy as np

import custom_argument_parser


# For use in some more serious error messages
MY_EMAIL = "gh10@sanger.ac.uk"


def generate_min_list(hist_dict):

    """
    Returns a list of the x values corresponding to the minima in hist_dict.
    """

    min_list = sig.argrelextrema(np.array(hist_dict.values()), np.less_equal,
                                 order=3)[0].tolist()

    return min_list


def generate_max_list(hist_dict):

    """
    Returns a list of the x values corresponding to the maxima in hist_dict.
    """

    max_list = sig.argrelextrema(np.array(hist_dict.values()),
                                 np.greater_equal, order=3)[0].tolist()

    return max_list


def find_kmer_depth(hist_dict, min_list=None):

    """
    Returns the x value of the maximum of the main peak in hist_dict. The main
    peak is the peak which should correspond to homozygous k-mers, it is not
    the downwards curve which normally appears at the beginning of the k-mer
    spectrum, which is due to base errors in the reads.
    """

    min_list = generate_min_list(hist_dict)
    max_list = generate_max_list(hist_dict)
    min_list_minimum = min(min_list)

    max_y_vals = [hist_dict[x] for x  in max_list if x in hist_dict.keys()]

    for maximum in sorted(max_y_vals)[::-1]:
        # Although this looks really inefficient, we hopefully
        # shouldn't have to do it more than twice (hopefully just once)
        # and I can't quickly think of a better way to do it
        for (key, value) in hist_dict.iteritems():
            if value == maximum and key > min_list_minimum and key > 10:
                return key

    print("ERROR: Could not find the maximum of the main peak",
          file=sys.stderr)
    sys.exit(1)


def find_start_main_peak(hist_dict):

    """
    Returns the smallest value in the list of minima. This normally constitutes
    the beginning of the main peak.
    """

    return min(generate_min_list(hist_dict))


def find_start_repeat_kmers(hist_dict, verbose=False):

    """
    Returns a point 'a' such that the x co-ordinate of the k-mer depth is
    equidistant between the start of the main peak and 'a'.
    """

    start_first_peak = find_start_main_peak(hist_dict)
    if verbose:
        print("Start of first peak =", start_first_peak)

    return (2 * find_kmer_depth(hist_dict)) - start_first_peak


def create_hist_dict(in_file):

    """
    Returns a dictionary containing the frequency distribution of k-mers in the
    set of reads. Each key is a specific number of occurrences, and the value
    paired with this key is the number of distinct k-mers which occur this many
    times in the reads.
    """

    hist_dict = {}

    for line in in_file.readlines():
        splat = [int(x) for x in line.strip().split()]
        hist_dict[splat[0]] = splat[1]

    return hist_dict


def count_num_kmers(hist_dict, lower_bound, upper_bound):

    """
    Count the total number of k-mers contained in the hist dict between
    lower_bound and upper_bound.
    """

    kmer_count = 0
    for (occ, freq) in hist_dict.items():
        if lower_bound <= occ <= upper_bound:
            kmer_count += (occ * freq)

    return kmer_count


def calculate_gri(number_repetitive_kmers, total_number_kmers):

    """
    Returns the percentage of repetitive k-mers.
    """

    gri = number_repetitive_kmers / total_number_kmers

    return gri


def create_parser():

    """
    Returns a parser based on my custom parser, which is stored in
    custom_argument_parser.py. I have done it in this way because I wanted to
    be able to write my own help and usage messages. These messages obviously
    need to be kept up-to-date as features are added/removed.
    """

    parser = custom_argument_parser.CustomParser()
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-c", "--indiv-repeat-cutoffs", type=int, nargs='+')
    parser.add_argument("-C", "--single-repeat-cutoff", type=int, nargs='?')
    parser.add_argument("-e", "--indiv-error-cutoffs", type=int, nargs='+')
    parser.add_argument("-E", "--single-error-cutoff", type=int, nargs='?')
    parser.add_argument("-u", "--indiv-upper-cutoffs", type=int, nargs='+')
    parser.add_argument("-U", "--single-upper-cutoff", type=int, nargs='?')
    parser.add_argument("-f", "--file", type=str, nargs='+', required=True)

    return parser


def parser_main():

    """
    Creates the parser and returns the parsed arguments.
    """

    parser = create_parser()
    args = parser.parse_args()

    return args


def error_check_user_cutoffs(indiv_cutoffs, single_cutoff, num_files,
                             cutoff_name):

    """
    Error checking for user specified cutoffs. This checks that mutually
    exclusive options have not both be specified, and also that the correct
    number of values have been specified (if appropriate).
    """

    if indiv_cutoffs and single_cutoff:
        print("ERROR: Cannot specify both --indiv-", cutoff_name,
              "-cutoffs and --single-", cutoff_name, "-cutoff", sep='',
              file=sys.stderr)
        sys.exit(1)

    if indiv_cutoffs is not None:
        if num_files != len(indiv_cutoffs):
            print("ERROR: Need to have the same number of individual",
                  cutoff_name, "cutoffs as files", file=sys.stderr)
            sys.exit(1)

        if any(cutoff <= 0 for cutoff in indiv_cutoffs):
            print("ERROR: --indiv-", cutoff_name, "-cutoffs must all be ",
                  "positive integers", sep='', file=sys.stderr)
            sys.exit(1)

    if single_cutoff is not None:
        if single_cutoff <= 0:
            print("ERROR: --single-", cutoff_name, "-cutoff must be a ",
                  "positive integer", sep='', file=sys.stderr)
            sys.exit(1)

    return


def construct_cutoff_list(indiv_cutoffs, single_cutoff, num_files):

    """
    This function takes the list of infividual cutoffs specified by the user
    (will be None if it was not specified); a single cutoff if specified by the
    user (will be None if it was not specified); and the number of input files.

    If individual cutoffs were specified then that list is returned; If a
    single cutoff was specified then a list is returned with that cutoff
    repeated once for each file; Otherwise, it is assumed that the user wants
    the value of this cutoff to be automatically estimated. As this is
    impossible to do at this stage in the program, a list is returned
    containing a 0 for each file. This is used later in the program as a
    signal that the cutoff still needs to be estimated.
    """

    # Check that only one of these has been set (user can only specify a single
    # cutoff for the entire set of files, or an individual cutoff for each
    # file, but not both). This is just here to make sure that I didn't mess up
    # my error checking earlier.
    assert not (indiv_cutoffs and single_cutoff), \
        "For error, repeat, and upper cutoffs, setting individual cutoffs " + \
        "and a single cutoff to be applied to the whole set of files is " + \
        "mutually exclusive. This should have been checked in the " + \
        "argument parsing error checking, but was clearly missed, and " + \
        "appeared here. Please email me about seeing this message (" + \
        MY_EMAIL + ")"

    if indiv_cutoffs and not single_cutoff:
        return indiv_cutoffs

    elif single_cutoff and not indiv_cutoffs:
        return [single_cutoff for _ in xrange(num_files)]

    else:
        return [0 for _ in xrange(num_files)]


def est_error_cutoff_if_required(hist_dict, verbose, initial_error_cutoff):

    """
    Carry out error cutoff estimation if required.
    """

    if initial_error_cutoff == 0:
        error_cutoff = find_start_main_peak(hist_dict)
        if verbose:
            print("Estimated error cutoff as", error_cutoff)

        return error_cutoff

    else:
        if verbose:
            print("Using user specified error cutoff of", initial_error_cutoff)

        return initial_error_cutoff


def est_repeat_cutoff_if_required(hist_dict, verbose, initial_repeat_cutoff):

    """
    Carry out repeat cutoff estimation if required.
    """

    if initial_repeat_cutoff == 0:
        repeat_cutoff = find_start_repeat_kmers(hist_dict, verbose)
        if verbose:
            print("Estimated start of repetitive k-mers as", repeat_cutoff)

        return repeat_cutoff

    else:
        if verbose:
            print("Using user specified start of reptitive k-mers as",
                  initial_repeat_cutoff)

        return initial_repeat_cutoff


def est_upper_cutoff_if_required(hist_dict, verbose, initial_upper_cutoff):

    """
    Carry out upper cutoff estimation if required.
    """

    if initial_upper_cutoff == 0:
        upper_cutoff = 20 * find_kmer_depth(hist_dict)
        if verbose:
            print("Estimated upper cutoff as", upper_cutoff)

        return upper_cutoff

    else:
        if verbose:
            print("Using user specified upper bound of", initial_upper_cutoff)

        return initial_upper_cutoff


def check_cutoff_consistency(error_cutoff, repeat_cutoff, upper_cutoff):

    """
    Check we aren't going to have problems with these cutoffs later. This is
    simple stuff, i.e. ensure that error cutoff < repeat cutoff < upper cutoff.
    etc.
    """

    if not (error_cutoff < repeat_cutoff < upper_cutoff):
        print("ERROR: The error cutoff must be smaller than the repeat",
              "cutoff, which in turn must be smaller than the upper cutoff.",
              file=sys.stderr)
        print("That is: error cutoff < repeat cutoff < upper cutoff.",
              file=sys.stderr)
        print("Skipping this file...", file=sys.stderr)

        return -1

    return 0


def process_histogram_file(file_name, verbose, in_error_cutoff,
                           in_upper_cutoff, in_repeat_cutoff):

    """
    Main function for interacting with an individual histogram file. This
    function creates a hist dict for the file, sets any cutoffs which still
    need to be set, and computes and prints the GRI for the file.
    """

    with open(file_name, 'r') as hist_file:

        print("Processing", file_name)

        hist_dict = create_hist_dict(hist_file)

        error_cutoff = est_error_cutoff_if_required(hist_dict, verbose,
                                                    in_error_cutoff)
        repeat_cutoff = est_repeat_cutoff_if_required(hist_dict, verbose,
                                                      in_repeat_cutoff)
        upper_cutoff = est_upper_cutoff_if_required(hist_dict, verbose,
                                                    in_upper_cutoff)

        if check_cutoff_consistency(error_cutoff, repeat_cutoff,
                                    upper_cutoff) == -1:
            return

        total_number_kmers = count_num_kmers(hist_dict, error_cutoff,
                                             upper_cutoff)
        number_repetitive_kmers = count_num_kmers(hist_dict, repeat_cutoff,
                                                  upper_cutoff)

        if verbose:
            print("K-mer depth =", find_kmer_depth(hist_dict))
            print("Total number of k-mers", total_number_kmers)
            print("Number of repetitive k-mers", number_repetitive_kmers)

        gri = calculate_gri(number_repetitive_kmers, total_number_kmers)

        if not (0 <= gri <= 1):
            print("ERROR: GRI is not between 0 and 1. This should never ",
                  "happen!!! Please email me and let me know about this (",
                  MY_EMAIL, ")", file=sys.stderr, sep='')

        print("GRI = %0.4f" %(gri))

    return


def process_user_cutoffs(indiv_cutoffs, single_cutoff, num_files, name):

    """
    This function takes the cutoff command line arguments for a specific type
    of cutoff (i.e. error, repeat, or upper) and constructs the correct cutoff
    list for this cutoff.
    """

    error_check_user_cutoffs(indiv_cutoffs, single_cutoff, num_files, name)

    return construct_cutoff_list(indiv_cutoffs, single_cutoff, num_files)


def main():

    """
    This is the main function for the program. It just calls the argument
    parser, calls the cutoff list constructors, and then iterates over the
    input files, calculating and printing their GRIs.
    """

    args = parser_main()
    file_paths = args.file
    num_files = len(file_paths)
    verbose = args.verbose

    error_cutoffs = process_user_cutoffs(args.indiv_error_cutoffs,
                                         args.single_error_cutoff,
                                         num_files, "error")
    repeat_cutoffs = process_user_cutoffs(args.indiv_repeat_cutoffs,
                                          args.single_repeat_cutoff,
                                          num_files, "repeat")
    upper_cutoffs = process_user_cutoffs(args.indiv_upper_cutoffs,
                                         args.single_upper_cutoff,
                                         num_files, "upper")

    for (file_name, repeat_cutoff, error_cutoff, upper_cutoff) in \
    zip(file_paths, repeat_cutoffs, error_cutoffs, upper_cutoffs):

        try:
            process_histogram_file(file_name, verbose, error_cutoff,
                                   upper_cutoff, repeat_cutoff)
        except IOError:
            print("ERROR: Could not open file \"" + file_name + "\".",
                  "Skipping...", file=sys.stderr)


if __name__ == "__main__":
    main()

