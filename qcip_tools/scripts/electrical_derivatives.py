#!/usr/bin/env python3
"""
Fetch the different derivatives of the energy with respect to electric field, and compute related quantities
"""

import argparse
import sys
from scipy import constants
import pandas as pd
import platform

import qcip_tools.scripts
from qcip_tools import derivatives_e, quantities
from qcip_tools.chemistry_files import helpers, PropertyNotPresent
from qcip_tools.transformations import TransformationMatrix as TM

__version__ = '0.3'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

__longdoc__ = """
Try to fetch the dipole moment and the frequency dependant (hyper)polarizabilities, then print the tensors, and
related quantities. Rely on the availability of ``electrical_derivatives``.

.. warning::

    In Dalton, by default, second hyperpolarizability with HF or DFT does not compute
    all components of the gamma tensor, but only the one that contribute to :math:`\\gamma_{||}`.
    Use ``.GAMALL`` in ``*CUBIC`` to do so.
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


def print_tensors(electrical_derivatives, representation, rm=None):
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

            if rm is not None:
                a = TM.rotate_matrix(electrical_derivatives[representation][freq].components, rm)
                electrical_derivatives[representation][freq].components = a

            print(electrical_derivatives[representation][freq].to_string(**kw))


def get_arguments_parser():
    arguments_parser = argparse.ArgumentParser(description=__doc__)
    arguments_parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    arguments_parser.add_argument('-csv', default=False, action='store_true')
    arguments_parser.add_argument('--save-components', default=False, action='store_true')
    arguments_parser.add_argument('-q', '--quiet', default=False, action='store_true')
    arguments_parser.add_argument('--vector', default=None,
                                  help='rotate tensors to match the (a b c) direction.', nargs=3, type=float)

    arguments_parser.add_argument(
        'infile',
        nargs='+',
        default=sys.stdin,
        help='Files containing the source of the derivatives')

    return arguments_parser


def append_data(electrical_derivatives, to_data_frames, file_id, save_components=False):
    representations = [
        'F',
        'FF',
        'dD',
        'FFF',
        'dDF',
        'FDd',
        'XDD',
        'FFFF',
        'dDFF',
        'dFFD',
        'XDDF',
        'dDDd',
        'xDDD',
    ]
    components = {
        'mu': [f'mu_{i}' for i in 'x y z'.split()],
        'alpha': [f'alpha_{i}{j}' for i in 'x y z'.split()
                  for j in 'x y z'.split()],
        'beta': [f'beta_{i}{j}{k}' for i in 'x y z'.split()
                 for j in 'x y z'.split()
                 for k in 'x y z'.split()],
        'gamma': [f'gamma_{i}{j}{k}{m}' for i in 'x y z'.split()
                  for j in 'x y z'.split()
                  for k in 'x y z'.split()
                  for m in 'x y z'.split()],
    }

    for representation in representations:
        if representation in electrical_derivatives:
            freqs = [x for x in electrical_derivatives[representation].keys()]
            freqs.sort(key=lambda x: derivatives_e.convert_frequency_from_string(x))

            name = derivatives_e.NAMES[representation].replace('(', '_').replace(')', '')

            # include dipole if found
            kw = {}
            if len(representation) == 3:
                if 'F' in electrical_derivatives:
                    kw['dipole'] = electrical_derivatives['F']['static'].components

            for freq in freqs:
                electrical_derivatives[representation][freq].compute_properties(**kw)
                df_name = '{}_w={}'.format(name, to_nanometer(freq))
                if to_data_frames.get(df_name) is None:
                    to_data_frames[df_name] = {}
                to_data_frames[df_name][file_id] = electrical_derivatives[representation][freq].properties
                if save_components:
                    prop = name.split('_')[0]
                    components_dict = {}
                    for c, v in zip(components[prop],
                                    electrical_derivatives[representation][freq].components.flatten()):
                        components_dict[c] = v
                    to_data_frames[df_name][file_id].update(components_dict)


def save_data(to_data_frames):
    for df_name, data in to_data_frames.items():
        df = pd.DataFrame(data).transpose()
        df.index.names = ['names']
        df.to_csv(df_name + '.csv')


def main():
    args = get_arguments_parser().parse_args()
    to_data_frames = {}

    if 'Windows' in platform.system():
        from glob import glob
        files = []
        for f in args.infile:
            files.extend(glob(f))
    else:
        files = args.infile

    for file_name in files:
        f = open(file_name)
        infile = helpers.open_chemistry_file(f)
        f.close()
        infile.file_name = file_name

        if not infile.has_property('electrical_derivatives'):
            return qcip_tools.scripts.exit_failure('cannot find electrical derivatives ({})'.format(infile.file_type))

        try:
            electrical_derivatives = infile.property('electrical_derivatives')
        except PropertyNotPresent:
            return qcip_tools.scripts.exit_failure('cannot find electrical derivatives ({})'.format(infile.file_type))

        rotation_matrix = None
        if args.vector is not None:
            if 'F' in electrical_derivatives:
                dipole = electrical_derivatives['F']['static'].components
                rotation_matrix = TM.rotation_matrix_from_vectors(dipole, args.vector)
                electrical_derivatives['F']['static'].components = TM.rotate_matrix(dipole, rotation_matrix)
            else:
                print('No dipole to define a rotation matrix')

        if not args.quiet:
            # mu
            if 'F' in electrical_derivatives:
                print('dipole moment:')
                print(electrical_derivatives['F']['static'].to_string())
                print('')

            # alpha
            print_tensors(electrical_derivatives, 'FF', rm=rotation_matrix)
            print_tensors(electrical_derivatives, 'dD', rm=rotation_matrix)

            # beta:
            print_tensors(electrical_derivatives, 'FFF', rm=rotation_matrix)
            print_tensors(electrical_derivatives, 'dDF', rm=rotation_matrix)
            print_tensors(electrical_derivatives, 'FDd', rm=rotation_matrix)
            print_tensors(electrical_derivatives, 'XDD', rm=rotation_matrix)

            # gamma
            print_tensors(electrical_derivatives, 'FFFF', rm=rotation_matrix)
            print_tensors(electrical_derivatives, 'dDFF', rm=rotation_matrix)
            print_tensors(electrical_derivatives, 'dFFD', rm=rotation_matrix)
            print_tensors(electrical_derivatives, 'XDDF', rm=rotation_matrix)
            print_tensors(electrical_derivatives, 'dDDd', rm=rotation_matrix)
            print_tensors(electrical_derivatives, 'XDDD', rm=rotation_matrix)

        if args.csv:
            append_data(electrical_derivatives, to_data_frames, infile.file_name, args.save_components)

    if args.csv:
        save_data(to_data_frames)


if __name__ == '__main__':
    main()
