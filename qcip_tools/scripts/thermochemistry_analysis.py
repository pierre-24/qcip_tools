#!/usr/bin/env python3
"""
Perform a thermochemistry analysis
"""

import argparse
import sys

import qcip_tools.scripts
from qcip_tools import derivatives_g, molecule as qcip_molecule, symmetry
from qcip_tools.chemistry_files import helpers, PropertyNotPresent

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

__longdoc__ = """
Try to fetch the hessian, and compute thermochemisty data out of that
Rely on the availability of ``geometrical_derivatives`` and ``computed_energies``.

Currently implemented for:

+ Gaussian FCHK
+ Dalton archive output (``DALTON.HES``)
+ Dalton LOG (only Hessian)
+ GAMESS output (only HF)

To do:
+ Gaussian LOG (only Hessian?)
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

    parser.add_argument('-n', '--symmetry-number', type=int, default=1, help='Symmetry number (for the rotation)')
    parser.add_argument('-t', '--temperature', type=float, default=298.15, help='Temperature (in K)')
    parser.add_argument('-p', '--pressure', type=float, default=1.01325e5, help='Pressure (in Pa)')
    parser.add_argument(
        '-s', '--scale', type=float, default=1.0, help='Scaling factor for the vibrational frequencies')
    parser.add_argument('-x', '--exclude', action='store', help='Exclude frequencies')
    parser.add_argument(
        '-V', '--verbose', help='Gives the detail of the different contributions', action='store_true')
    parser.add_argument('-f', '--factor', type=float, default=1.0, help='multiply energies by a given factor')

    parser.add_argument(
        '-g', '--guess-symmetry', help='Guess the symmetry number', action='store_true')

    return parser


def main():
    parser = get_arguments_parser()
    args = parser.parse_args()

    if not args.infile.has_property('geometrical_derivatives'):
        qcip_tools.scripts.exit_failure('cannot find geometrical derivatives ({})'.format(args.infile.file_type))

    if not args.infile.has_property('computed_energies'):
        qcip_tools.scripts.exit_failure('cannot find energies ({})'.format(args.infile.file_type))

    try:
        geometrical_derivatives = args.infile.property('geometrical_derivatives')
        energies = args.infile.property('computed_energies')
    except PropertyNotPresent as e:
        return qcip_tools.scripts.exit_failure('cannot find {} ({})'.format(e, args.infile.file_type))

    if 'GG' not in geometrical_derivatives:
        return qcip_tools.scripts.exit_failure('no Hessian in this file')

    if 'total' not in energies:
        return qcip_tools.scripts.exit_failure('cannot find total energy!')

    molecule = args.infile.molecule
    mwh = derivatives_g.MassWeightedHessian(
        molecule, cartesian_hessian=geometrical_derivatives['GG'], scale=args.scale)

    # modify list of frequencies
    if args.exclude:
        frequencies_to_exlude = args.exclude.split(',')
        for f in frequencies_to_exlude:
            try:
                n = int(f)
            except ValueError:
                print('{} is not a valid number'.format(f))
                sys.exit(1)

            if not (1 <= n <= mwh.dof):
                return qcip_tools.scripts.exit_failure(
                    '{} is not in the range of allowed frequencies (1-{})'.format(n, mwh.dof))

            try:
                index = mwh.included_modes.index(n - 1)
                if index >= 0:
                    mwh.included_modes.pop(index)
            except ValueError:
                return qcip_tools.scripts.exit_failure('{} is not in the possibilities'.format(n))

    symmetry_number = args.symmetry_number
    if args.guess_symmetry:
        try:
            mol = args.infile.property('molecule')
        except PropertyNotPresent:
            return qcip_tools.scripts.exit_failure('cannot find molecule ({})'.format(args.infile.file_type))

        sf = qcip_molecule.MolecularSymmetryFinder(mol, tol=1e-3)
        g, _, _ = sf.find_symmetry()

        if g.symbol in [
                symmetry.PointGroupType.cyclic, symmetry.PointGroupType.pyramidal, symmetry.PointGroupType.reflexion]:
            symmetry_number = g.order
            if g.order == -1:
                symmetry_number = 1
        elif g.symbol in [
                symmetry.PointGroupType.dihedral,
                symmetry.PointGroupType.prismatic,
                symmetry.PointGroupType.antiprismatic
        ]:
            symmetry_number = g.order
            if g.order == -1:
                symmetry_number = 2
        elif g.symbol == symmetry.PointGroupType.improper_rotation:
            symmetry_number = g.order / 2
        elif g.symbol in [symmetry.PointGroupType.tetrahedral_achiral, symmetry.PointGroupType.tetrahedral_chiral]:
            symmetry_number = 12
        elif g.symbol in [symmetry.PointGroupType.octahedral_achiral, symmetry.PointGroupType.octahedral_chiral]:
            symmetry_number = 24

        print('(found {}, using symmetry number={})'.format(g, symmetry_number))

    print('With T={} K and P={} Pa'.format(args.temperature, args.pressure))

    # compute thermo
    U = list(mwh.compute_internal_energy(args.temperature))
    H = list(mwh.compute_enthalpy(args.temperature))
    S = list(mwh.compute_entropy(args.symmetry_number, args.temperature, args.pressure))
    G = list(mwh.compute_gibbs_free_energy(symmetry_number, args.temperature, args.pressure))

    # prepare data
    zpva = mwh.compute_zpva()
    electronic_energy = energies['total']

    def treat(x, f=False):
        if not f:
            x.append(x[-1] + zpva)
            x.append(sum(x[:3]) + zpva)
            x.append(x[-1] + electronic_energy)
        else:
            x.append(x[-1])
            x.append(sum(x[:3]))
            x.append(x[-1])

        if not args.verbose:
            x = x[-1:]

        x = list(i * args.factor for i in x)

        return x

    S = treat(S, True)
    U = treat(U)
    H = treat(H)
    G = treat(G)

    # print data
    def format_stuffs(name, lst):
        s = '{:5}'.format(name)
        for e in lst:
            s += ' {: .8e}'.format(e)

        return s

    dashes = '-' * (5 + 16 * (6 if args.verbose else 1))

    print(dashes)

    print('     ', end='')
    if args.verbose:
        print('  trans.          rot.            vib.            vib. + ZPVA     total w/o scf ', end='')
    print('  total w. scf')

    print(dashes)
    print(format_stuffs('S', S))
    print(format_stuffs('U(T)', U))
    print(format_stuffs('H(T)', H))
    print(format_stuffs('G(T)', G))
    print(dashes)


if __name__ == '__main__':
    main()
