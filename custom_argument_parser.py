
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
Custom error messages for the parser used in GRIn, as I couldn't otherwise get
the usage messages to do what I wanted
"""


import sys
import argparse


class CustomParser(argparse.ArgumentParser):

    """Subclass to get custom error messages"""

    @staticmethod
    def generate_usage_message():

        """Custom usage message"""

        return "usage: grin [-h] [options] -f file [file ...]\n"

    def generate_help_message(self):

        """Custom help message"""

        store_str = "\nGRIn version 1.1.0\n\n"
        store_str += self.generate_usage_message()
        store_str += "\nrequired arguments:\n"
        store_str += "\t-f, --file                   "
        store_str += "one or more input histogram files\n"
        store_str += "\t                             (as output by "
        store_str += "Jellyfish)\n\n"

        store_str += "optional arguments:\n"
        store_str += "\t-v, --verbose                "
        store_str += "print more output\n\n"
        store_str += "\t-e, --indiv-error-cutoffs    "
        store_str += "list of error cutoffs for specifying the end of\n"
        store_str += "\t                             the error curve\n"
        store_str += "\t-E, --single-error-cutoff    "
        store_str += "single error cutoff to be applied to all files\n\n"
        store_str += "\t-r, --indiv-repeat-cutoffs   "
        store_str += "list of repeat cutoffs for start of repetitive\n"
        store_str += "\t                             k-mers\n"
        store_str += "\t-R, --single-repeat-cutoff   "
        store_str += "single repeat cutoff to be applied to all files\n\n"
        store_str += "\t-u, --indiv-upper-cutoffs    "
        store_str += "list of upper cutoffs after which k-mers do not\n"
        store_str += "\t                             contribute to total or "
        store_str += "repetitive k-mer counts\n"
        store_str += "\t-U, --single-upper-cutoff    "
        store_str += "single upper cutoff to be applied to all files\n\n"

        store_str += "other arguments:\n"
        store_str += "\t-h, --help                   print this message\n\n"

        store_str += "notes:\n"
        store_str += "\t* If specified, the number of individual "
        store_str += "cutoffs must equal the number of\n"
        store_str += "\t  input files.\n"

        store_str += "\t* All cutoffs must be positive integers.\n"

        store_str += "\t* For each type of cutoff (error, repeat, upper), the "
        store_str += "single cutoff and\n"
        store_str += "\t  individual cutoffs options are mutual exclusive.\n"

        store_str += "\t* If any cutoff is not set by the user then it will "
        store_str += "be estimated and set by\n"
        store_str += "\t  GRIn: The error cutoff will be set to be the "
        store_str += "minimum following the\n"
        store_str += "\t  error curve; The repeat cutoff will be set to be "
        store_str += "the point R such that the \n"
        store_str += "\t  k-mer depth is equidistant between R and the error "
        store_str += "cutoff; the upper\n"
        store_str += "\t  cutoff will be set to be 20 * the k-mer depth.\n\n"

        store_str += "Copyright Genome Research Limited 2016. Licenced under "
        store_str += "the GNU GPL (see COPYING).\n\n"

        store_str += "See https://github.com/george-hall/GRIn for more "
        store_str += "information.\n\n"

        return store_str

    def print_usage(self, file=None):

        """Print custom usage message"""

        if file is None:
            file = sys.stdout
        self._print_message(self.generate_usage_message(), file)

    def print_help(self, file=None):

        """Print custom help message"""

        if file is None:
            file = sys.stdout
        self._print_message(self.generate_help_message(), file)


def create_parser():

    """
    Returns my custom parser with the arguments added.
    """

    parser = CustomParser()
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-r", "--indiv-repeat-cutoffs", type=int, nargs='+')
    parser.add_argument("-R", "--single-repeat-cutoff", type=int, nargs='?')
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

