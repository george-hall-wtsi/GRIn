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
the paper accompanying this program (still being written as of this commit).
For anything to do with this program, contact George Hall (gh10@sanger.ac.uk).
"""


from __future__ import print_function, division

import math
import sys
import subprocess as sp

import spectrum_related as spec
import custom_argument_parser
import cutoff_related as cutoffs

# Non-breaking imports:

try:
    import matplotlib.pyplot as plt
except ImportError:
    MATPLOTLIB_PRESENT = False
else:
    MATPLOTLIB_PRESENT = True

try:
    import numpy as np
except ImportError:
    NUMPY_PRESENT = False
else:
    NUMPY_PRESENT = True

# For use in some more serious error messages
MY_EMAIL = "gh10@sanger.ac.uk"


def check_matplotlib_present():

    """
    Print error message and exit if Matplotlib was not successfully imported.
    """

    if not MATPLOTLIB_PRESENT:
        print("ERROR: Could not find Matplotlib installation. Exiting.",
              file=sys.stderr)
        sys.exit(1)
    else:
        return


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


def run_jellyfish(file_paths, verbosity):

    """Generate histogram using Jellyfish"""

    # Options used for Jellyfish. Change them here if you want:
    JELLYFISH_BIN = "jellyfish"
    K_MER_SIZE = "31"
    HASH_TABLE_SIZE = "100M" # Can use S.I. units M & G,
    NUM_THREADS = "25"

    if verbosity > 0:
        print("Counting k-mers with Jellyfish...")

    sp.call([JELLYFISH_BIN, "count", "-m", K_MER_SIZE, "-s", HASH_TABLE_SIZE,
             "-t", NUM_THREADS, "-C"] + file_paths)

    hist_name = generate_hist_file_name(file_paths)

    if verbosity > 0:
        print("Storing histogram in file '", hist_name, "'", sep='')

    with open(hist_name, 'w') as hist_file:
        sp.call([JELLYFISH_BIN, "histo", "mer_counts.jf"], stdout=hist_file)

    return


def generate_hist_file_name(file_names):

    """
    Return the name of the histogram file, computed by concatenating all input
    file names together using underscores
    """

    return "_".join(file_names) + ".hist"


def plot_histogram(hist_dict, error_cutoff, repeat_cutoff, upper_cutoff):

    """Plot histogram using hist_dict"""

    data = [[], []]
    data[0] = list(hist_dict.keys())
    data[1] = list(hist_dict.values())

    plt.plot(data[0], data[1])

    plt.axvline(error_cutoff)
    plt.axvline(repeat_cutoff)
    plt.axvline(upper_cutoff)

    x_vals_array = np.array(data[0])
    y_vals_array = np.array(data[1])

    boundary_func = lambda lower, upper: np.logical_and(x_vals_array >= lower,
                                                        x_vals_array <= upper)

    # Shade error curve
    plt.fill_between(x_vals_array, y_vals_array,
                     where=(boundary_func(1, error_cutoff)), interpolate=True,
                     alpha=0.5, color="red")

    # Shade main peak
    plt.fill_between(x_vals_array, y_vals_array,
                     where=(boundary_func(error_cutoff, repeat_cutoff)),
                     interpolate=True, alpha=0.5, color="blue")

    # Shade repetitive k-mers
    plt.fill_between(x_vals_array, y_vals_array,
                     where=(boundary_func(repeat_cutoff, upper_cutoff)),
                     interpolate=True, alpha=0.5, color="green")

    # Shade abundant k-mers
    plt.fill_between(x_vals_array, y_vals_array,
                     where=(x_vals_array > upper_cutoff), interpolate=True,
                     alpha=0.5, color="black")

    plt.xlim(1, upper_cutoff * 1.05)
    plt.ylim(1, max(hist_dict.values()) * 1.1)

    return


def convert_bp_to_SI(num_bp):

    """
    Returns a string with the number of base pairs num_bp expressed in 'SI'
    units (i.e. Mbp) with the correct suffix attached.
    """

    thousands_power = math.log(num_bp, 1000)

    if thousands_power < 1:
        return str(num_bp) + "bp"

    elif 1 <= thousands_power < 2:
        num_kbp = num_bp / 1000
        return "{0:.1f}Kbp".format(num_kbp)

    elif 2 <= thousands_power < 3:
        num_mbp = num_bp / 1000000
        return "{0:.1f}Mbp".format(num_mbp)

    else:
        num_gbp = num_bp / 1000000000
        return "{0:.1f}Gbp".format(num_gbp)


def process_histogram_file(file_name, initial_error_cutoff,
                           initial_repeat_cutoff, initial_upper_cutoff,
                           verbosity):

    """
    Main function for interacting with an individual histogram file. This
    function creates a hist dict for the file, sets the cutoffs for the file,
    and computes and prints the GRI for the file.
    """

    with open(file_name, 'r') as hist_file:

        print("Processing", file_name)

        hist_dict = spec.create_hist_dict(hist_file)

        error_cutoff = cutoffs.set_error_cutoff(hist_dict,
                                                initial_error_cutoff,
                                                verbosity)
        if error_cutoff == -1:
            return -1

        repeat_cutoff = cutoffs.set_repeat_cutoff(hist_dict,
                                                  initial_repeat_cutoff,
                                                  error_cutoff,
                                                  verbosity)
        if repeat_cutoff == -1:
            return -1

        upper_cutoff = cutoffs.set_upper_cutoff(hist_dict,
                                                initial_upper_cutoff,
                                                verbosity)
        if upper_cutoff == -1:
            return -1

        if cutoffs.check_cutoff_consistency(error_cutoff, repeat_cutoff,
                                            upper_cutoff) == -1:
            return -1

        total_num_kmers_used = spec.count_num_kmers(hist_dict, error_cutoff,
                                                    upper_cutoff)
        number_repetitive_kmers = spec.count_num_kmers(hist_dict,
                                                       repeat_cutoff,
                                                       upper_cutoff)

        if verbosity > 0:
            print("Total number of k-mers", total_num_kmers_used)
            print("Number of repetitive k-mers", number_repetitive_kmers)

            # Don't exit if not able to import Scipy and Numpy
            if spec.SCIPY_PRESENT and NUMPY_PRESENT:
                kmer_depth = spec.find_kmer_depth(hist_dict)
                print("K-mer depth =", kmer_depth)

                genome_size_est = int(total_num_kmers_used / kmer_depth)
                print("Genome size estmination =", genome_size_est,
                      "(" + convert_bp_to_SI(genome_size_est) + ")")

        if verbosity == 2:
            plot_histogram(hist_dict, error_cutoff, repeat_cutoff,
                           upper_cutoff)
            plt.title(file_name)

        gri = calculate_gri(number_repetitive_kmers, total_num_kmers_used)

        if gri != -1:
            print("GRI = %0.4f" %(gri))

    return


def error_check_user_input(args):

    """
    Perform the following checks on the command line options given by the user:

        * Check that --full-auto has not been set if any manual cutoff has also
          been set
        * Check that --verbosity is no greater than 2
    """

    if args.full_auto and cutoffs.any_cutoff_set(args):
        print("ERROR: Cannot set both --full-auto and any manual cutoff",
              file=sys.stderr)
        sys.exit(1)

    if args.verbosity > 2:
        print("ERROR: --verbosity cannot be greater than 2 (either -v or -vv",
              "must be used)", file=sys.stdout)
        sys.exit(1)
    return


def generate_subplot_func(num_subplots):

    """
    Returns a function which can be called to generate the next subplot. I've
    made it work like this because the same function needs to be called every
    time except the value x needs to be incremented by 1.
    """

    return lambda x: plt.subplot(num_subplots, num_subplots, x,
                                 xlabel="Number of Occurrences",
                                 ylabel="Distinct k-mers with Occurence",
                                 yscale="log")


def main():

    """
    This is the main function for the program. It just calls the argument
    parser, calls the cutoff list constructors, and then iterates over the
    input files, calculating and printing their GRIs.
    """

    args = custom_argument_parser.parser_main()
    file_paths = args.file
    verbosity = args.verbosity

    if verbosity > 0:
        print("Command ran:", " ".join(sys.argv))

    if verbosity == 2:

        check_matplotlib_present()

        # num_subplots is the required number of subplots per row and column
        # in the window in order to accomodate all files
        num_subplots = math.ceil(math.sqrt(len(args.file)))
        # file_counter keeps track of which file we are on
        file_counter = 1
        subplot_func = generate_subplot_func(num_subplots)

    error_check_user_input(args)

    if args.full_auto:
        # Jellyfish needs to be run first in order to generate histogram file
        run_jellyfish(file_paths, verbosity)
        file_paths = [generate_hist_file_name(file_paths)]

    (error_cutoffs,
     repeat_cutoffs,
     upper_cutoffs) = cutoffs.construct_all_cutoff_lists(args)

    for (file_name, repeat_cutoff, error_cutoff, upper_cutoff) in \
    zip(file_paths, repeat_cutoffs, error_cutoffs, upper_cutoffs):

        try:
            if verbosity == 2:
                subplot_func(file_counter)
                file_counter += 1
            process_histogram_file(file_name, error_cutoff, repeat_cutoff,
                                   upper_cutoff, verbosity)

        except IOError:
            print("ERROR: Could not open file \"" + file_name + "\".",
                  "Skipping...", file=sys.stderr)

    if verbosity == 2:
        plt.gcf().canvas.set_window_title("GRIn Histogram plot")
        plt.show()


if __name__ == "__main__":
    main()

