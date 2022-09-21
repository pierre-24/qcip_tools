#!/usr/bin/env python3
"""
Fetch the different derivatives of the energy with respect to geometrical modifications
"""

import argparse
import sys
from qcip_tools import chemistry_files, derivatives_g
from qcip_tools.chemistry_files import helpers, PropertyNotPresent

from qcip_tools.scripts import commons

__version__ = '0.2'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

__longdoc__ = """
Try to fetch the gradient and hessian. If the second one is found, compute vibrational frequencies.
Rely on the availability of ``geometrical_derivatives``.

Currently implemented for:

+ Gaussian FCHK
+ Dalton archive output (``DALTON.HES`` if ``.HESPUN``!)
+ Dalton LOG (CC and ``RESPONSE`` gradient and Hesssian, if available)

To do:

+ GAMESS output (only Hessian?)
+ Gaussian LOG (only Hessian?)
"""


def print_tensor(geometrical_derivatives, representation, molecule):
    if representation in geometrical_derivatives:
        name = derivatives_g.NAMES[representation]
        print(name)
        if molecule is not None:
            print(geometrical_derivatives[representation].to_string(molecule=molecule))
        else:
            print(geometrical_derivatives[representation].to_string())


def make_NN(geometrical_derivatives, molecule):
    """

    :rtype: derivatives_g.MassWeightedHessian
    """
    return derivatives_g.MassWeightedHessian(molecule, geometrical_derivatives['GG'].components)


def get_arguments_parser():
    arguments_parser = argparse.ArgumentParser(description=__doc__)
    arguments_parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    arguments_parser.add_argument(
        'infile',
        nargs='?',
        default=sys.stdin,
        action=helpers.create_open_chemistry_file_action(),
        help='source of the derivatives')

    arguments_parser.add_argument(
        '-N', '--no-vibs', action='store_true', help='Do not perform the vibrational analysis')

    return arguments_parser


def main():
    args = get_arguments_parser().parse_args()

    if not args.infile.has_property('geometrical_derivatives'):
        return commons.exit_failure('cannot find geometrical derivatives ({})'.format(args.infile.file_type))

    try:
        geometrical_derivatives = args.infile.property('geometrical_derivatives')
    except PropertyNotPresent:
        return commons.exit_failure('cannot find geometrical derivatives ({})'.format(args.infile.file_type))

    molecule = None
    if issubclass(type(args.infile), chemistry_files.WithMoleculeMixin):
        molecule = args.infile.molecule

    print_tensor(geometrical_derivatives, 'G', molecule)
    print_tensor(geometrical_derivatives, 'GG', molecule)

    if 'GG' in geometrical_derivatives and molecule is not None and not args.no_vibs:
        hessian = geometrical_derivatives['GG']

        mwh = make_NN(geometrical_derivatives, molecule)
        if molecule.linear():
            print('!! This molecule is linear')

        mwh.included_modes = list(range(0, mwh.dof))  # include all modes

        print('Low frequencies (cm⁻¹):')
        for i in range(hessian.trans_plus_rot):
            print('{: 10.4f}'.format(
                mwh.frequencies[i] * derivatives_g.MassWeightedHessian.HARTREE_TO_WAVENUMBER_CONVERSION), end=' ')
        print('\n')

        print('Frequencies (cm⁻¹): ', end='')
        for i in range(hessian.trans_plus_rot, hessian.spacial_dof):
            if (i - hessian.trans_plus_rot) % 6 == 0:
                print('')

            print('{: 10.4f}'.format(
                mwh.frequencies[i] * derivatives_g.MassWeightedHessian.HARTREE_TO_WAVENUMBER_CONVERSION), end=' ')

        print('\n')

        print('Total number of normal modes:', hessian.spacial_dof)
        print('Number of vibrational modes: ', hessian.spacial_dof - hessian.trans_plus_rot)
        print('')

        print(mwh.output_displacements())


if __name__ == '__main__':
    main()
