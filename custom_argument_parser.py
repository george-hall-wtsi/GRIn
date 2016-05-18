
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
        store_str = "usage: grin [-h] [-v] [-a] [-u upper-bound] [-i] "
        store_str += "[-E single-error-cutoff] "
        store_str += "[-e error-cutoff [error-cutoff ...]] "
        store_str += "[-c cutoff [cutoff ...]] -f file [file ...]\n"

        return store_str

    def generate_help_message(self):
        """Custom help message"""
        store_str = "\n"
        store_str += self.generate_usage_message()
        store_str += "\nrequired arguments:\n"
        store_str += "\t-f, --file                 "
        store_str += "one or more input histogram files\n\n"

        store_str += "optional arguments:\n"
        store_str += "\t-v, --verbose              "
        store_str += "print more output\n"
        store_str += "\t-a, --analyzer             "
        store_str += "run kmerspectrumanalyzer and compute GRI from that "
        store_str += "histogram instead\n"
        store_str += "\t-u, --upper-bound          "
        store_str += "upper bound after which k-mers do not contribute to "
        store_str += "total or repetitive k-mer counts\n"
        store_str += "\t-i, --ignore-error         "
        store_str += "estimate the erroneous k-mers and do not include them "
        store_str += "in the total k-mer count\n"
        store_str += "\t-e, --error-cutoffs        "
        store_str += "manual list of the end of the error peak (one per "
        store_str += "file)\n"
        store_str += "\t-E, --single-error-cutoff  manual cutoff of error "
        store_str += "peak, applied to all files\n"
        store_str += "\t-c, --cutoffs              "
        store_str += "list of manual cutoffs for start of repetitive "
        store_str += "k-mers\n\n"

        store_str += "other arguments:\n"
        store_str += "\t-h, --help                 print this message\n\n"

        store_str += "notes:\n"
        store_str += "\t* if specified, the number of manually specified "
        store_str += "cutoffs must equal the number of input files\n"
        store_str += "\t* if specified, the number of manually specified "
        store_str += "error cutoffs must equal the number of input files\n"
        store_str += "\t* if specified, --upper-bound, --single-error-cutoff, "
        store_str += "and --error-cutoffs must be integers\n"
        store_str += "\t* --single-error-cutoff, --error-cutoffs, and "
        store_str += "--ignore-error are mutually exclusive\n\n"

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
