#!/usr/bin/env python3
"""
Extract the excited states information
"""

import argparse
import sys
import math

import qcip_tools.scripts
from qcip_tools import derivatives_exci, derivatives_e
from qcip_tools.chemistry_files import helpers, PropertyNotPresent

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

__longdoc__ = """
You can

+ Change the unit of the energy
+ Limit the number of state that appears
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

    parser.add_argument('-l', '--limit', type=float, help='Only print excitation above a given limit')
    parser.add_argument('-L', '--limit-CSFs', type=float, help='Only print CSF above a limit', default=1)

    parser.add_argument(
        '-E', '--unit', type=str, help='Unit for the energy', choices=['au', 'eV', 'nm', 'cm-1'], default='au')

    return parser


def main():
    parser = get_arguments_parser()
    args = parser.parse_args()

    if not args.infile.has_property('excitations'):
        qcip_tools.scripts.exit_failure('cannot find excitations ({})'.format(args.infile.file_type))

    try:
        excitations = args.infile.property('excitations')
    except PropertyNotPresent as e:
        return qcip_tools.scripts.exit_failure('cannot find {} ({})'.format(e, args.infile.file_type))

    if '!F' not in excitations:
        return qcip_tools.scripts.exit_failure('no transition dipoles (!F)')

    if '!' not in excitations and '#' not in excitations:
        return qcip_tools.scripts.exit_failure('no excitation energies (# or !)')

    has_S2 = '_<S2>' in excitations
    if has_S2:
        s2_0 = excitations['_<S2>'][0]

    has_CSFs = '_CSFs' in excitations

    handler = derivatives_exci.Excitations(transition_energies=excitations['!'], transition_dipoles=excitations['!F'])

    if handler.nstates > 1:
        print('-' * (68 + (14 if has_S2 else 0)))
        print('#    {:10}   X            Y            Z            f'.format('E ({})'.format(args.unit)), end='')
        if has_S2:
            print('            <Ŝ²>   Δ<Ŝ²>', end='')
        print()

        print('-' * (68 + (14 if has_S2 else 0)))

        for i in range(1, handler.nstates):
            f_ge = handler.oscillator_strength(i)
            if args.limit and f_ge < args.limit:
                continue

            dipole = handler.transition_dipole(i)
            print('{:<3} {: .5e} {: .5e} {: .5e} {: .5e} {: .5e}'.format(
                i,
                derivatives_e.convert_energy_to(handler.transition_energy(i), args.unit),
                dipole[0], dipole[1], dipole[2],
                f_ge), end='')

            if has_S2:
                s2_i = excitations['_<S2>'][i]
                print(' {: .3f} {: .3f} {}'.format(
                    s2_i, s2_i - s2_0, '!' if math.fabs(s2_i - s2_0) > .2 else ''), end='')
            if has_CSFs:
                print('|', excitations['_CSFs'][i].to_string(limit=args.limit_CSFs), end='')
            print()

        print('-' * (68 + (14 if has_S2 else 0)))


if __name__ == '__main__':
    main()
