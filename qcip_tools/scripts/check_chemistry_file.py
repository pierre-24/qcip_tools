#!/usr/bin/env python3
"""
checks whether qcip_tools identify a given file. If so, it can also check if it have a given property.
"""

from qcip_tools.chemistry_files import helpers
import argparse
import sys

from qcip_tools.scripts import commons  # noqa

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

__longdoc__ = """
Prints 1 if it is the case, 0 if not.

.. note::

    Having this property does not mean that the property in question is indeed available in this given file,
    just that it *could* present this property.
"""


# program options
def get_arguments_parser():
    arguments_parser = argparse.ArgumentParser(description=__doc__)
    arguments_parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    arguments_parser.add_argument(
        'infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='file to check')
    arguments_parser.add_argument(
        '-p', '--property', action='store', help='List of atoms')
    arguments_parser.add_argument(
        '-T', '--not-trust-extension', action='store_true', help='Do not trust extension of the file (slower)')

    return arguments_parser


# main
def main():
    args = get_arguments_parser().parse_args()

    try:
        f = helpers.open_chemistry_file(args.infile, trust_extension=not args.not_trust_extension)

        if args.property:
            print(1 if f.has_property(args.property) else 0)
        else:
            print(1)
    except helpers.ProbablyNotAChemistryFile:
        print(0)


if __name__ == '__main__':
    main()
