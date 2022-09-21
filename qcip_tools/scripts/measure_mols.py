#!/usr/bin/env python3
"""
Measure distances/angles/dihedrals in molecules
"""

import argparse
from typing import List, Tuple, Union
import sys
import numpy

from qcip_tools.chemistry_files import helpers, ChemistryFile
from qcip_tools import math as qm


__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

__longdoc__ = """
Measure distances/angles/dihedrals in molecules.
"""


def input_criteria(value) -> List[Tuple[int, ...]]:
    if value[0] == '@':
        try:
            with open(value[1:]) as f:
                value = f.read()
        except FileNotFoundError:
            raise argparse.ArgumentTypeError('unable to open {}'.format(value[1:]))

    criteria_ = value.split(',')
    final = []
    for criterion in criteria_:
        info = criterion.split('-')
        if not (2 <= len(info) <= 4):
            raise argparse.ArgumentTypeError('{} should contain 2-4 indices'.format(criterion))

        try:
            info = [int(i) - 1 for i in info]
        except ValueError:
            raise argparse.ArgumentTypeError('{} should contains only integers'.format(criterion))

        final.append(tuple(info))

    return final


def get_observations(
        structures: List[ChemistryFile],
        requests: List[Union[Tuple[int, ...], str]],
) -> Tuple[numpy.ndarray, List[str]]:

    labels = []

    for req in requests:
        labels.append({
            2: 'l({},{})',
            3: 'a({},{},{})',
            4: 'd({},{},{},{})'
        }[len(req)].format(*[a + 1 for a in req]))

    data = numpy.zeros((len(structures), len(requests)))

    for i, structure in enumerate(structures):
        for j, request in enumerate(requests):
            coos = [structure.molecule[e].position for e in request]

            data[i, j] = {
                2: qm.distance,
                3: qm.angle,
                4: qm.torsion_angle
            }[len(request)](*coos)

    return data, labels


def get_arguments_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    parser.add_argument(
        'infiles',
        nargs='*',
        type=argparse.FileType('r'),
        default=sys.stdin,
        help='sources')

    parser.add_argument('-c', '--criteria', help='list of criteria to measure', type=input_criteria, required=True)

    return parser


def main():
    parser = get_arguments_parser()
    args = parser.parse_args()

    infiles = []
    names = []
    for i in args.infiles:
        infiles.append(helpers.open_chemistry_file(i))
        names.append(i.name)

    values, labels = get_observations(infiles, args.criteria)

    # reminder of input
    print('#', ','.join('-'.join(str(c + 1) for c in criterion) for criterion in args.criteria))

    # results
    max_name_length = max(len(n) for n in names)
    print(' ' * max_name_length, ' '.join('{:>15}'.format(label) for label in labels))
    N = 16 * len(args.criteria) + max_name_length
    print('-' * N)
    for i in range(values.shape[0]):
        print(
            ('{{:<{}}}'.format(max_name_length)).format(names[i]),
            ' '.join('{:>15.5f}'.format(d) for d in values[i]))


if __name__ == '__main__':
    main()
