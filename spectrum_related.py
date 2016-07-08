
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
Collection of functions related to dealing with k-mer spectra.
"""

from __future__ import print_function, division

import sys

import import_checks

# Non-breaking imports

try:
    import numpy as np
except ImportError:
    NUMPY_PRESENT = False
else:
    NUMPY_PRESENT = True

try:
    import scipy.signal as sig
except ImportError:
    SCIPY_PRESENT = False
else:
    SCIPY_PRESENT = True


def generate_min_list(hist_dict):

    """
    Returns a list of the occurrence values corresponding to the minima in the
    k-mer spectra represented by hist_dict.
    """

    import_checks.check_scipy_present(SCIPY_PRESENT)
    import_checks.check_numpy_present(NUMPY_PRESENT)

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

    import_checks.check_scipy_present(SCIPY_PRESENT)
    import_checks.check_numpy_present(NUMPY_PRESENT)

    # hist_keys is a list containing the keys of hist_dict
    # hist_vals is a list containing the values of hist_dict
    (hist_keys, hist_vals) = (list(hist_dict.keys()), list(hist_dict.values()))

    max_list = sig.argrelmax(np.array(hist_vals), order=3)[0].tolist()

    return [hist_keys[x] for x in max_list]


def find_kmer_depth(hist_dict):

    """
    Returns the x value of the maximum of the main peak in hist_dict. The main
    peak is the peak which should correspond to homozygous k-mers, it is not
    the downwards curve which normally appears at the beginning of the k-mer
    spectrum, which is due to base errors in the reads.
    """

    max_list = generate_max_list(hist_dict)
    start_main_peak = find_start_main_peak(hist_dict)

    if start_main_peak == -1:
        return -1

    max_y_vals = [hist_dict[x] for x in max_list if x in hist_dict.keys()]

    for maximum in sorted(max_y_vals)[::-1]:
        # Although this looks really inefficient, we hopefully
        # shouldn't have to do it more than twice (hopefully just once)
        # and I can't quickly think of a better way to do it
        for (key, value) in hist_dict.items():
            if value == maximum and key > start_main_peak and key > 10:
                return key

    print("ERROR: Could not find the maximum of the main peak",
          file=sys.stderr)
    sys.exit(1)


def find_start_main_peak(hist_dict):

    """
    Returns the smallest value in the list of minima. This normally constitutes
    the beginning of the main peak.
    """

    min_list = generate_min_list(hist_dict)

    if not min_list:
        print("ERROR: Could not generate a list of minima for this file.",
              "Skipping...")
        return -1

    return min(min_list)


def find_start_repeat_kmers(hist_dict, error_cutoff, verbosity):

    """
    Returns a point 'a' such that the x co-ordinate of the k-mer depth is
    equidistant between the start of the main peak and 'a'.
    """

    if verbosity > 0:
        print("Using error cutoff of", error_cutoff, "in repeat cutoff",
              "in repeat cutoff estimation")

    kmer_depth = find_kmer_depth(hist_dict)

    if kmer_depth == -1:
        return -1

    return (2 * kmer_depth) - error_cutoff


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

