#!/usr/bin/env python3
"""
Fetch the different derivatives of the energy with respect to electric field, and compute related quantities
"""

import argparse
import sys
from scipy import constants

import qcip_tools.scripts
from qcip_tools import derivatives_e, quantities
from qcip_tools.chemistry_files import helpers, PropertyNotPresent

__version__ = '0.3'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

__longdoc__ = """
Try to fetch the dipole moment and the frequency dependant (hyper)polarizabilities, then print the tensors, and
related quantities. Rely on the availability of ``electrical_derivatives``.

Currently implemented for:

+ Dalton archive output (CC and normal responses).
  Probably not working for explicit operator settings (please use ``.FREQUE`` along with things like ``.SHG`` or so),
  and other type of operators (different from ``DIPLEN``).
+ Gaussian FCHK, with sign correction (but **not the second hyperpolarizability tensors**).

To do:

+ GAMESS output (not easy, and don't forget the DFT one)
+ Gaussian LOG (all located in the same place)
+ Dalton LOG (not located in the same place)

.. warning::

    By default, second hyperpolarizability with HF or DFT does not compute all components of the gamma tensor, but only
    the one that contribute to :math:`\\gamma_{||}`.
"""


def to_nanometer(val):
    """Convert frequency to nanometer

    :param val:
    :return:
    """

    if type(val) is str:
        return val

    converted = \
        constants.h * constants.c / (val * quantities.convert(quantities.ureg.hartree, quantities.ureg.joule)) * 1e9

    return '{:.1f}nm'.format(converted)


def print_tensors(electrical_derivatives, representation):
    if representation in electrical_derivatives:
        freqs = [x for x in electrical_derivatives[representation].keys()]
        freqs.sort(key=lambda x: derivatives_e.convert_frequency_from_string(x))

        name = derivatives_e.NAMES[representation]

        # include dipole if found
        kw = {}
        if len(representation) == 3:
            if 'F' in electrical_derivatives:
                kw['dipole'] = electrical_derivatives['F']['static'].components

        for freq in freqs:
            print('{}, w={} ({:.6f} a.u.)'.format(
                name, to_nanometer(freq), derivatives_e.convert_frequency_from_string(freq)))

            print(electrical_derivatives[representation][freq].to_string(**kw))


def get_arguments_parser():
    arguments_parser = argparse.ArgumentParser(description=__doc__)
    arguments_parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    arguments_parser.add_argument(
        'infile',
        nargs='?',
        default=sys.stdin,
        action=helpers.create_open_chemistry_file_action(),
        help='source of the derivatives')

    return arguments_parser


def main():
    args = get_arguments_parser().parse_args()

    if not args.infile.has_property('electrical_derivatives'):
        return qcip_tools.scripts.exit_failure('cannot find electrical derivatives ({})'.format(args.infile.file_type))

    try:
        electrical_derivatives = args.infile.property('electrical_derivatives')
    except PropertyNotPresent:
        return qcip_tools.scripts.exit_failure('cannot find electrical derivatives ({})'.format(args.infile.file_type))

    # mu
    if 'F' in electrical_derivatives:
        print('dipole moment:')
        print(electrical_derivatives['F']['static'].to_string())
        print('')

    # alpha
    print_tensors(electrical_derivatives, 'FF')
    print_tensors(electrical_derivatives, 'dD')

    # beta:
    print_tensors(electrical_derivatives, 'FFF')
    print_tensors(electrical_derivatives, 'dDF')
    print_tensors(electrical_derivatives, 'FDd')
    print_tensors(electrical_derivatives, 'XDD')

    # gamma
    print_tensors(electrical_derivatives, 'FFFF')
    print_tensors(electrical_derivatives, 'dDFF')
    print_tensors(electrical_derivatives, 'dFFD')
    print_tensors(electrical_derivatives, 'XDDF')
    print_tensors(electrical_derivatives, 'dDDd')
    print_tensors(electrical_derivatives, 'XDDD')


if __name__ == '__main__':
    main()
