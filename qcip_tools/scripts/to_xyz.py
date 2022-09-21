#!/usr/bin/env python3
"""
Convert any chemistry file that contain a geometry to a xyz file
"""

import argparse
import sys

import qcip_tools.scripts
from qcip_tools.chemistry_files import helpers, PropertyNotPresent, xyz

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

__longdoc__ = """
*Nothing yet*
"""


def get_arguments_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    parser.add_argument(
        'infile',
        nargs='?',
        default=sys.stdin,
        action=helpers.create_open_chemistry_file_action(),
        help='source of the derivatives')

    return parser


def main():
    parser = get_arguments_parser()
    args = parser.parse_args()

    try:
        mol = args.infile.property('molecule')
    except PropertyNotPresent:
        return qcip_tools.scripts.exit_failure('cannot find molecule ({})'.format(args.infile.file_type))

    xf = xyz.File.from_molecule(mol)
    print(xf.to_string())


if __name__ == '__main__':
    main()
