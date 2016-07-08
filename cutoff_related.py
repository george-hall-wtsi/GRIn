
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
Collection of functions related to dealing with cutoffs in GRIn.
"""

from __future__ import print_function, division

import sys

import spectrum_related as spec

# Duplicated - need to fix...
MY_EMAIL = "gh10@sanger.ac.uk"

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


def any_cutoff_set(args):

    """Return True if any cutoff has been set by the user."""

    return any([args.indiv_error_cutoffs, args.single_error_cutoff,
                args.indiv_repeat_cutoffs, args.single_repeat_cutoff,
                args.indiv_upper_cutoffs, args.single_upper_cutoff])


def construct_all_cutoff_lists(args):

    """
    Return a tuple containing a cutoff list for each of the three cutoffs
    """

    if not args.full_auto:

        # User has passed in histogram

        # Check user has not set illegal cutoffs
        error_check_user_cutoffs(args)

        num_files = len(args.file)

        # Construct cutoff lists
        error_cutoffs = construct_cutoff_list(args.indiv_error_cutoffs,
                                              args.single_error_cutoff,
                                              num_files)
        repeat_cutoffs = construct_cutoff_list(args.indiv_repeat_cutoffs,
                                               args.single_repeat_cutoff,
                                               num_files)
        upper_cutoffs = construct_cutoff_list(args.indiv_upper_cutoffs,
                                              args.single_upper_cutoff,
                                              num_files)

    else:
        # User has passed in fast{a,q} file
        error_cutoffs = [0]
        repeat_cutoffs = [0]
        upper_cutoffs = [0]

    return (error_cutoffs, repeat_cutoffs, upper_cutoffs)


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


def set_error_cutoff(hist_dict, initial_error_cutoff, verbosity):

    """
    Determine whether it is necessary to estimate the error cutoff, or if the
    user has already specified one. If it is necessary to estimate one, then
    carry out this estimation, and return the resulting cutoff. Otherwise,
    simply return the cutoff which the user has specified.
    """

    if initial_error_cutoff == 0:
        error_cutoff = spec.find_start_main_peak(hist_dict)

        if error_cutoff == -1:
            return -1

        if verbosity > 0:
            print("Estimated error cutoff as", error_cutoff)

        return error_cutoff

    else:
        if verbosity > 0:
            print("Using user specified error cutoff of", initial_error_cutoff)

        return initial_error_cutoff


def set_repeat_cutoff(hist_dict, initial_repeat_cutoff, error_cutoff,
                      verbosity):

    """
    Determine whether it is necessary to estimate the repeat cutoff, or if the
    user has already specified one. If it is necessary to estimate one, then
    carry out this estimation, and return the resulting cutoff. Otherwise,
    simply return the cutoff which the user has specified.
    """

    if initial_repeat_cutoff == 0:
        repeat_cutoff = spec.find_start_repeat_kmers(hist_dict, error_cutoff,
                                                     verbosity)

        if repeat_cutoff == -1:
            return -1

        if verbosity > 0:
            print("Estimated start of repetitive k-mers as", repeat_cutoff)

        return repeat_cutoff

    else:
        if verbosity > 0:
            print("Using user specified start of repetitive k-mers as",
                  initial_repeat_cutoff)

        return initial_repeat_cutoff


def create_window_generator(hist_dict, window_size):

    """
    Create a generator which generates all windows of size window_size in
    hist_dict. These are generated in a sorted manner: the smallest number of
    occurrences in an (occ, freq) pair in window 'n' will be less than the
    smallest number of occurrences in an (occ, freq) pair in window 'n+1'.
    """

    sorted_items = sorted(hist_dict.items(), key=lambda item: item[0])

    for i in range(0, len(sorted_items)):
        next_window = sorted_items[i:i + window_size]
        if len(next_window) == window_size:
            yield next_window


def mean_diff(window):

    """
    Return the mean difference between each number in the window and its
    predecessor.
    """

    diff_list = []

    for (_, freq1), (_, freq2) in zip(window, window[1:]):
        diff_list.append(abs(freq1 - freq2))

    return (sum(diff_list) / len(diff_list))


def pad_hist_dict(hist_dict):

    """
    'Pad' hist_dict with zeroes when an occurrence value is missing. For
    example, if there is no frequency value recorded for k-mers occurring 10
    times in the set of reads, then add the pair (occ = 10, freq = 0) to
    hist_dict.
    """

    for occ in range(1, sorted(list(hist_dict.keys()))[-1]):
        hist_dict.setdefault(occ, 0)

    return hist_dict


def fluctuation_method(hist_dict, kmer_depth, verbosity):

    """
    Return the number of occurrences which corresponds to the midpoint of the
    first window of size window_size in which the mean difference between a
    point's frequency and its neighbouring point's frequency is less than
    difference_cutoff. Increase the difference cutoff if no Upper Cutoff could
    be calculated using the current difference cutoff A window is a list of
    length window_size which consists of (occ, freq) pairs from hist_dict.
    """

    window_size = 6
    plateau_cutoff = 2
    padded_hist_dict = pad_hist_dict(hist_dict)

    if verbosity > 0:
        print("Number of plateaus required for successful estimation:",
              plateau_cutoff)

    for difference_cutoff in [0, 0.1, 0.5, 1, 2, 5]:

        if verbosity > 0:
            print("Looking for plateaus using difference cutoff of",
                  difference_cutoff)

        # Keep track of how many plateaus we have found thus far. Once we
        # have found a number of plateaus equal to plateau_cutoff, return the
        # most recently found one as the Upper Cutoff.
        plateaus_found = 0

        window_generator = create_window_generator(padded_hist_dict,
                                                   window_size)

        for window in window_generator:
            midpoint_occ_num = window[int(window_size/2)][0]

            if midpoint_occ_num <= kmer_depth:
                continue

            if mean_diff(window) <= difference_cutoff:
                if plateaus_found < plateau_cutoff:
                    plateaus_found += 1
                    if verbosity > 0:
                        print("Found plateau",
                              str(plateaus_found) + "/" + str(plateau_cutoff))
                    continue

                else:
                    upper_cutoff = midpoint_occ_num
                    if verbosity > 0:
                        print("Estimated upper cutoff as", upper_cutoff,
                              "using difference cutoff of",
                              difference_cutoff)

                    return upper_cutoff

    return -1


def set_upper_cutoff(hist_dict, initial_upper_cutoff, verbosity):

    """
    Determine whether it is necessary to estimate the upper cutoff, or if the
    user has already specified one. If it is necessary to estimate one, then
    carry out this estimation, and return the resulting cutoff. Otherwise,
    simply return the cutoff which the user has specified.
    """

    if initial_upper_cutoff == 0:

        # Upper Cutoff needs to be estimated

        kmer_depth = spec.find_kmer_depth(hist_dict)

        if kmer_depth == -1:
            return -1

        upper_cutoff_est = fluctuation_method(hist_dict, kmer_depth, verbosity)

        if upper_cutoff_est != -1:
            # Fluctuation estimation method was successful
            return upper_cutoff_est

        # Fluctuation estimation method was unsuccessful
        print("WARNING: Failed to estimate the upper cutoff using the",
              "fluctuation method.")
        print("Will use an upper cutoff of 20 * the k-mer depth.")

        fallback_estimate = 20 * kmer_depth
        print("Estimated upper cutoff as", fallback_estimate)

        return fallback_estimate

    else:
        if verbosity > 0:
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
