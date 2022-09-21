#!/usr/bin/env python3
"""
Compute the Boltzmann population of molecules
"""

import argparse
import sys
import math

from qcip_tools.chemistry_files import helpers, ChemistryFile, PropertyNotPresent
from qcip_tools import derivatives_g


__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

__longdoc__ = """
Compute the Boltzmann population of molecules, based on:

+ their total energy (``-c E``, default),
+ their internal energy (``-c U``),
+ their enthalpy (``-c H``),
+ their free Gibbs energy (``-c G``).

Checks if the molecular formula of all inputs matches
(but that's it, e.g. not that they have been obtained with the same method).
"""


def get_arguments_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    parser.add_argument('-t', '--temperature', type=float, default=298.15, help='Temperature (in K)')
    parser.add_argument('-p', '--pressure', type=float, default=1.01325e5, help='Pressure (in Pa)')

    parser.add_argument(
        '-c', '--criterion', choices=['E', 'H', 'G', 'U'], default='E', help='Criterion for the population')
    parser.add_argument('-f', '--factor', type=float, default=1.0, help='multiply energies by a given factor')

    parser.add_argument(
        'infiles',
        nargs='*',
        type=argparse.FileType('r'),
        default=sys.stdin,
        help='sources')

    return parser


def get_energies(f: ChemistryFile, T: float, p: float, n: int = 1, req: str = 'E') -> float:
    total_energy = f.property('computed_energies')['total']
    energies = {'E': total_energy}

    if req != 'E':
        try:
            geoms = f.property('geometrical_derivatives')
            if 'GG' in geoms:
                mwh = derivatives_g.MassWeightedHessian(f.molecule, cartesian_hessian=geoms['GG'])

                energies.update(**{
                    'U': total_energy + sum(mwh.compute_internal_energy(T)),
                    'H': total_energy + sum(mwh.compute_enthalpy(T)),
                    'G': total_energy + sum(mwh.compute_gibbs_free_energy(n, T, p))
                })
        except PropertyNotPresent:
            pass

    return energies[req]


def main():
    parser = get_arguments_parser()
    args = parser.parse_args()

    # read files and extract energies
    infiles = []
    names = []
    energies = []

    formula = ''

    for i in args.infiles:
        fp = helpers.open_chemistry_file(i)

        if not formula:
            formula = fp.molecule.formula()
        elif formula != fp.molecule.formula():
            raise Exception('molecular formulas do not match! ({}!={})'.format(fp.molecule.formula(), formula))

        infiles.append(fp)
        names.append(i.name)
        energies.append(get_energies(fp, args.temperature, args.pressure, req=args.criterion))

    # compute and output
    lowest = min(energies)
    deltas_energy = [e - lowest for e in energies]
    R = 8.31446261815324 * derivatives_g.MassWeightedHessian.ENERGY_IN_AU_CONVERSION
    denominator = sum(math.exp(-d / (R * args.temperature)) for d in deltas_energy)

    for i, delta in enumerate(deltas_energy):
        print('{}\t{:.5e}\t{:.2f} %'.format(
            names[i],
            delta * args.factor,
            math.exp(-delta / (R * args.temperature)) / denominator * 100))


if __name__ == '__main__':
    main()
