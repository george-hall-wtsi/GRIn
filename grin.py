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
    Returns a list of the occurrence values corresponding to the minima in the
    k-mer spectra represented by hist_dict.
    """

    # hist_keys is a list containing the keys of hist_dict
    # hist_vals is a list containing the values of hist_dict
    (hist_keys, hist_vals) = (list(hist_dict.keys()), list(hist_dict.values()))

    min_list = sig.argrelmin(np.array(hist_vals), order=3)[0].tolist()

    return [hist_keys[x] for x in min_list]


def generate_max_list(hist_dict):

    """
    Returns a list of the occurrence values corresponding to the maxima in the
    k-mer spectra represented by hist_dict.
    """

    # hist_keys is a list containing the keys of hist_dict
    # hist_vals is a list containing the values of hist_dict
    (hist_keys, hist_vals) = (list(hist_dict.keys()), list(hist_dict.values()))

    max_list = sig.argrelmax(np.array(hist_vals), order=3)[0].tolist()

    return [hist_keys[x] for x in max_list]


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

    max_y_vals = [hist_dict[x] for x in max_list if x in hist_dict.keys()]

    for maximum in sorted(max_y_vals)[::-1]:
        # Although this looks really inefficient, we hopefully
        # shouldn't have to do it more than twice (hopefully just once)
        # and I can't quickly think of a better way to do it
        for (key, value) in hist_dict.items():
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


def find_start_repeat_kmers(hist_dict, error_cutoff, verbose):

    """
    Returns a point 'a' such that the x co-ordinate of the k-mer depth is
    equidistant between the start of the main peak and 'a'.
    """

    if verbose:
        print("Using error cutoff of", error_cutoff, "in repeat cutoff",
              "in repeat cutoff estimation")

    return (2 * find_kmer_depth(hist_dict)) - error_cutoff


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
    Return the GRI (that is, the percentage of repetitive k-mers).
    """

    if total_number_kmers == 0:
        print("ERROR: It seems that there were zero total k-mers. It's weird ",
              "this has happened - it's probably an issue with your cutoffs. ",
              "Please email me and let me know about this happening (",
              MY_EMAIL, "). Skipping this file...", file=sys.stderr, sep='')
        return -1

    gri = number_repetitive_kmers / total_number_kmers

    if not (0 <= gri <= 1):
        print("ERROR: GRI is not between 0 and 1. This should never ",
              "happen!!! Please email me and let me know about this (",
              MY_EMAIL, "). Skipping this file...", file=sys.stderr, sep='')
        return -1

    return gri


def error_check_user_cutoffs(args):

    """
    Error checking for user specified cutoffs. This checks that mutually
    exclusive options have not both be specified, and also that the correct
    number of values have been specified (if appropriate).
    """

    num_files = len(args.file)

    for (cutoff_name, indiv_cutoffs, single_cutoff) in \
            [\
            ("error", args.indiv_error_cutoffs, args.single_error_cutoff),\
            ("repeat", args.indiv_repeat_cutoffs, args.single_repeat_cutoff),\
            ("upper", args.indiv_upper_cutoffs, args.single_upper_cutoff)\
            ]:

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
        return [single_cutoff for _ in range(num_files)]

    else:
        return [0 for _ in range(num_files)]


def set_error_cutoff(hist_dict, initial_error_cutoff, verbose):

    """
    Determine whether it is necessary to estimate the error cutoff, or if the
    user has already specified one. If it is necessary to estimate one, then
    carry out this estimation, and return the resulting cutoff. Otherwise,
    simply return the cutoff which the user has specified.
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


def set_repeat_cutoff(hist_dict, initial_repeat_cutoff, error_cutoff, verbose):

    """
    Determine whether it is necessary to estimate the repeat cutoff, or if the
    user has already specified one. If it is necessary to estimate one, then
    carry out this estimation, and return the resulting cutoff. Otherwise,
    simply return the cutoff which the user has specified.
    """

    if initial_repeat_cutoff == 0:
        repeat_cutoff = find_start_repeat_kmers(hist_dict, error_cutoff,
                                                verbose)
        if verbose:
            print("Estimated start of repetitive k-mers as", repeat_cutoff)

        return repeat_cutoff

    else:
        if verbose:
            print("Using user specified start of reptitive k-mers as",
                  initial_repeat_cutoff)

        return initial_repeat_cutoff


def set_upper_cutoff(hist_dict, initial_upper_cutoff, verbose):

    """
    Determine whether it is necessary to estimate the upper cutoff, or if the
    user has already specified one. If it is necessary to estimate one, then
    carry out this estimation, and return the resulting cutoff. Otherwise,
    simply return the cutoff which the user has specified.
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


def process_histogram_file(file_name, initial_error_cutoff,
                           initial_repeat_cutoff, initial_upper_cutoff,
                           verbose):

    """
    Main function for interacting with an individual histogram file. This
    function creates a hist dict for the file, sets the cutoffs for the file,
    and computes and prints the GRI for the file.
    """

    with open(file_name, 'r') as hist_file:

        print("Processing", file_name)

        hist_dict = create_hist_dict(hist_file)

        error_cutoff = set_error_cutoff(hist_dict, initial_error_cutoff,
                                        verbose)
        repeat_cutoff = set_repeat_cutoff(hist_dict,
                                          initial_repeat_cutoff, error_cutoff,
                                          verbose)
        upper_cutoff = set_upper_cutoff(hist_dict, initial_upper_cutoff,
                                        verbose)

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

        if gri != -1:
            print("GRI = %0.4f" %(gri))

    return


def main():

    """
    This is the main function for the program. It just calls the argument
    parser, calls the cutoff list constructors, and then iterates over the
    input files, calculating and printing their GRIs.
    """

    args = custom_argument_parser.parser_main()
    file_paths = args.file
    num_files = len(file_paths)
    verbose = args.verbose

    if verbose:
        print("Command ran:", " ".join(sys.argv))

    # Check user has not set illegal cutoffs
    error_check_user_cutoffs(args)

    # Construct cutoff lists
    error_cutoffs = construct_cutoff_list(args.indiv_error_cutoffs,
                                          args.single_error_cutoff, num_files)
    repeat_cutoffs = construct_cutoff_list(args.indiv_repeat_cutoffs,
                                           args.single_repeat_cutoff,
                                           num_files)
    upper_cutoffs = construct_cutoff_list(args.indiv_upper_cutoffs,
                                          args.single_upper_cutoff, num_files)

    for (file_name, repeat_cutoff, error_cutoff, upper_cutoff) in \
    zip(file_paths, repeat_cutoffs, error_cutoffs, upper_cutoffs):

        try:
            process_histogram_file(file_name, error_cutoff, repeat_cutoff,
                                   upper_cutoff, verbose)
        except IOError:
            print("ERROR: Could not open file \"" + file_name + "\".",
                  "Skipping...", file=sys.stderr)


if __name__ == "__main__":
    main()

