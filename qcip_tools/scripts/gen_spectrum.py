#!/usr/bin/env python3
"""
Generate a spectrum
"""

import argparse
import sys
import math
import numpy

from scipy import constants

from qcip_tools import derivatives_exci, derivatives_e, quantities
from qcip_tools.chemistry_files import helpers, PropertyNotPresent

from qcip_tools.scripts import commons

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

__longdoc__ = """
Generate spectrum the point to plot a spectrum in a given window of energy (or wavelength).

The `y` values are in arbitrary units.

Currently, **only fetch UV intensities** and plot them, for

+ Gaussian (FCHK)
+ Dalton (archive)

To do, with the Hessian:

+ IR spectra (I need the derivatives of the dipole moment) ;
+ Raman spectra (I need the derivatives of the polarizability).

"""

UNITS = ['nm', 'au', 'ev', 'cm-1']


def is_interval(s_interval):
    """get interval

    :param s_interval: string of the form `min:max:unit`
    :type s_interval: str
    :rtype: tuple
    """

    inf = s_interval.split(':')
    if len(inf) != 3:
        raise argparse.ArgumentTypeError('limits must be `min:max:unit`')

    try:
        min_, max_ = float(inf[0]), float(inf[1])
    except ValueError:
        raise argparse.ArgumentTypeError('min and max must be float')

    if min_ < .0 or max_ < .0:
        raise argparse.ArgumentTypeError('minimum should be larger or equal to zero')

    unit = inf[2].lower()
    if unit not in UNITS:
        raise argparse.ArgumentTypeError('unit must be {}'.format(', '.join(UNITS)))

    if unit == 'nm' and .0 in [min_, max_]:
        raise argparse.ArgumentTypeError('0nm requested, which is impossible')

    return min_, max_, unit


def is_list(lst):
    final_list = []
    elements = lst.split(',')
    for element in elements:
        if ':' in element:  # range
            inf = element.split(':')
            if len(inf) != 2:
                raise argparse.ArgumentTypeError('range {} should contain two elements'.format(element))
            try:
                min_, max_ = int(inf[0]), int(inf[1])
            except ValueError:
                raise argparse.ArgumentTypeError('{} are not numbers'.format(element))

            final_list.extend(range(min_, max_ + 1))

        else:
            try:
                value = int(element)
            except ValueError:
                raise argparse.ArgumentTypeError('{} is not a number'.format(element))

            final_list.append(value)

    if len(final_list) == 0:
        return None

    return list(set(final_list))


def is_fwhm(fwhm):
    """Convert FWHM into an energy in Hartree

    :param fwhm: fwhm, as string
    :type fwhm: str
    :return: energy, in Hartree
    :rtype: float
    """

    return derivatives_e.convert_frequency_from_string(fwhm)


def get_arguments_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    parser.add_argument(
        'infile',
        nargs='?',
        default=sys.stdin,
        action=helpers.create_open_chemistry_file_action(),
        help='source of the derivatives')

    parser.add_argument(
        '-l', '--limits', help='Limit and units of the graph: `min:max:unit`', type=is_interval, default='0:9:eV')

    parser.add_argument(
        '-e', '--each', help='Interval (in unit of limits)', default=.1, type=float)

    parser.add_argument(
        '-I', '--impulses', help='Add impulses (at the end)', action='store_true')

    parser.add_argument(
        '-L', '--limit-impulses', help='Only output impulses above that threshold', default=.0, type=float)

    parser.add_argument(
        '-m', '--maximums', help='Add maximums (at the end, after the impulses if any)', action='store_true')

    # scaling:
    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        '-n', '--normalize', help='Set maximum intensity to 1', action='store_true')

    group.add_argument(
        '-s', '--scale', help='Scale the intensities', type=float)

    # exclude/include
    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        '-E', '--exclude', help='Exclude some peaks (list starts at 0)', type=is_list)

    group.add_argument(
        '-i', '--include', help='Include only some peaks (list starts at 0)', type=is_list)

    group.add_argument(
        '-D', '--decontaminate',
        help='Exclude peaks for which Δ<Ŝ²> is larger than this threshold (only relevant if <S²> is available)',
        type=float,
        default=1e3)

    subparsers = parser.add_subparsers(dest='subparser_name', help='source of the spectrum (`uv` or `ir`)')

    # ---- UV
    parser_uv = subparsers.add_parser('uv')
    parser_uv.add_argument(
        '-f', '--fwhm', help='Full Width at Half Minimum (FWHM)', default='.3eV', type=is_fwhm)

    parser_uv.add_argument('-L', '--lorentzian', help='Use a Lorentzian instead of a Gaussian', action='store_true')

    # ---- IR
    parser_ir = subparsers.add_parser('ir')
    parser_ir.add_argument(
        '-f', '--fwhm', help='Full Width at Half Minimum (FWHM)', default='10cm-1', type=is_fwhm)

    parser_ir.add_argument('-G', '--gaussian', help='Use a Gaussian instead of a Lorentzian', action='store_true')

    return parser


def gen_gaussian(x, mu, fwhm, area):
    """Compute a gaussian.

    .. math::

        f(x) = a\\,e^{-\\frac{4\\,ln(2)\\,(x-\\mu)^2}{\\omega^2}},

    with :math:`\\omega` the FWHM (Full Width at Half Maximum), :math:`\\mu` the center of the gaussian, and

    .. math::

        a = \\frac{2\\,A}{\\omega}\\,\\sqrt{\\frac{ln(2)}{2\\pi}},

    with :math:`A` the total area under the curve of the gaussian.

    :param x: table of values in eV
    :type x: numpy.ndarray|float
    :param mu: center of the gaussian (mean value)
    :type mu: float
    :param fwhm: FWHM
    :type fwhm: float
    :param area: total area (under the curve) of the gaussian
    :type area: float
    :rtype: numpy.ndarray
    """

    a = area * 2 / fwhm * math.sqrt(math.log(2) / (2 * math.pi))
    return a * numpy.exp(-4 * math.log(2) * (x - mu) ** 2 / fwhm**2)


def gen_lorentzian(x, mu, fwhm, area):
    """Compute a Laurentzian.

    .. math::

        f(x) = a\\,\\frac{1}{1+\\frac{4}{\\omega}\\,(x-\\mu)^2}

    with :math:`\\omega` the FWHM (Full Width at Half Maximum), :math:`\\mu` the center of the gaussian, and

    .. math::

        a = \\frac{2*A}{\\omega}\\,\\sqrt{\\frac{ln(2)}{\\pi}},

    with :math:`A` the total area under the curve of the gaussian.

    :param x: table of values in eV
    :type x: numpy.ndarray
    :param mu: center of the gaussian (mean value)
    :type mu: float
    :param fwhm: FWHM
    :type fwhm: float
    :param area: total area (under the curve) of the gaussian
    :type area: float
    :rtype: numpy.ndarray
    """
    a = area * 2 / (fwhm * math.pi)
    return a / (1 + 4 * (x - mu)**2 / fwhm**2)


def main():
    parser = get_arguments_parser()
    args = parser.parse_args()

    peaks = []

    if args.subparser_name == 'uv':
        # fetch uv transitions
        try:
            excitations = args.infile.property('excitations')
        except PropertyNotPresent as e:
            return commons.exit_failure('cannot find {} ({})'.format(e, args.infile.file_type))

        if '!F' not in excitations:
            return commons.exit_failure('no transition dipoles (!F)')

        if '!' not in excitations and '#' not in excitations:
            return commons.exit_failure('no excitation energies (# or !)')

        handler = derivatives_exci.Excitations(
            transition_energies=excitations['!'], transition_dipoles=excitations['!F'])

        s2_0 = 0
        has_S2 = '_<S2>' in excitations
        if has_S2:
            s2_0 = excitations['_<S2>'][0]

        for i in range(1, handler.nstates):
            # remove "contaminated" peak
            if has_S2:
                s2_i = excitations['_<S2>'][i]
                if math.fabs(s2_i - s2_0) > args.decontaminate:
                    continue

            peaks.append((handler.transition_energy(i), handler.oscillator_strength(i)))

        # broadening:
        f_broadening = gen_gaussian
        if args.lorentzian:
            f_broadening = gen_lorentzian

    else:
        raise NotImplementedError(args.subparser_name)

    # min and max
    min_, max_, unit = args.limits
    x = numpy.arange(min_, max_, args.each if min_ < max_ else -args.each, dtype=float)
    y = numpy.zeros(x.shape[0])

    if unit == 'nm':
        e = constants.c * constants.h / 1e-9 * quantities.convert(quantities.ureg.joule, quantities.ureg.hartree)
        xp = e / x
    elif unit == 'cm-1':
        xp = x * quantities.convert(quantities.ureg.wavenumber, quantities.ureg.hartree)
    elif unit == 'au':
        xp = x
    else:
        xp = x * quantities.convert(quantities.ureg.eV, quantities.ureg.hartree)

    impulses = []

    for i, (peak_energy, peak_area) in enumerate(peaks):
        if args.exclude is not None and i in args.exclude:
            continue
        if args.include is not None and i not in args.include:
            continue

        if peak_area > .0:
            y += f_broadening(xp, peak_energy, args.fwhm, peak_area)

        impulses.append((derivatives_e.convert_energy_to(peak_energy, unit), peak_area))

    scaling_factor = 1.0
    if args.normalize:
        scaling_factor = numpy.max(y)
    elif args.scale:
        scaling_factor = args.scale

    print('# intensities wrt energy (in {}), scaling={}'.format(unit, scaling_factor))
    for i in zip(x, y):
        print('{:.6e}\t{:.6e}'.format(i[0], i[1] / scaling_factor))

    if args.impulses:
        print('\n')
        print('# impulses (fourth column → height of the broadening function)')
        if args.limit_impulses > .0:
            print('# NOTE: only impulses with f >= {}'.format(args.limit_impulses))
        for i, impulse in enumerate(impulses):
            if min_ < impulse[0] < max_ and impulse[1] >= args.limit_impulses:
                print('{:<4}\t{:.6e}\t{:.6e}\t{:.6e}'.format(
                    i + 1, *impulse, f_broadening(0, 0, args.fwhm, impulse[1]) / scaling_factor))

    if args.maximums:
        # Due to ben741 (https://gist.github.com/ben741/d8c70b608d96d9f7ed231086b237ba6b)
        print('\n')
        print('# maximums')
        positions = numpy.where((y[1:-1] > y[0:-2]) * (y[1:-1] > y[2:]))[0] + 1
        for p in positions:
            print('{:.6e}\t{:.6e}'.format(x[p], y[p] / scaling_factor))


if __name__ == '__main__':
    main()
