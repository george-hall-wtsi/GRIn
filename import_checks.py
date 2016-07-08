
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
Functions to check whether or not specific modules have been imported.
"""

from __future__ import print_function, division

import sys

def check_scipy_present(SCIPY_PRESENT):

    """
    Print error message and exit if Scipy was not successfully imported.
    """

    if not SCIPY_PRESENT:
        print("ERROR: Could not find Scipy installation. Exiting.",
              file=sys.stderr)
        print("If you were trying to estimate cutoffs, you need this package",
              "installed.", file=sys.stderr)
        print("You can run GRIn with manually specified cutoffs without Scipy",
              "being installed.", file=sys.stderr)
        sys.exit(1)
    else:
        return


def check_numpy_present(NUMPY_PRESENT):

    """
    Print error message and exit if Numpy was not successfully imported.
    """

    if not NUMPY_PRESENT:
        print("ERROR: Could not find Numpy installation. Exiting.",
              file=sys.stderr)
        print("If you were trying to estimate cutoffs, you need this package",
              "installed.", file=sys.stderr)
        print("You can run GRIn with manually specified cutoffs without Numpy",
              "being installed", file=sys.stderr)
        sys.exit(1)
    else:
        return

