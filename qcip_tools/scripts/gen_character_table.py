#!/usr/bin/env python3
"""
Generate character table for a given group
"""

import argparse
import re

from qcip_tools import symmetry

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

__longdoc__ = """
Generate the character table for a group.

Input is, for example, ``C2v`` (or ``C_2v``), ``D3d`` (or ``D_3d``), ``Td`` (or ``T_d``), ...

.. warning::

    + No warranty is given. This has been tested, but not on every group ;
    + The order of the class is (more than probably) not correct ;
    + Complex conjugate representations are not merged together, but are marked with a star.
    + Large groups may take time to generate.
    + Does not work for :math:`I_h` !
"""


GROUP = re.compile('(?P<main>[CDSTOI])_?(?P<order>[0-9]?)(?P<flavor>[vhd]?)')


def is_group(group):
    m = GROUP.match(group)
    if not m:
        raise argparse.ArgumentTypeError('Unknown point group {}'.format(group))

    m, o, f = m.groups()
    order = 0

    if m == 'C':
        try:
            order = int(o)
        except ValueError:
            raise argparse.ArgumentTypeError('point group {}: missing order'.format(group))
        if f == '':
            symbol = symmetry.PointGroupType.cyclic
        elif f == 'v':
            symbol = symmetry.PointGroupType.pyramidal
        elif f == 'h':
            symbol = symmetry.PointGroupType.reflexion
        else:
            raise argparse.ArgumentTypeError('point group {}: wrong type'.format(group))
    elif m == 'S':
        try:
            order = int(o)
        except ValueError:
            raise argparse.ArgumentTypeError('point group {}: missing order'.format(group))

        if f != '':
            raise argparse.ArgumentTypeError('point group {}: wrong type'.format(group))

        symbol = symmetry.PointGroupType.improper_rotation

    elif m == 'D':
        try:
            order = int(o)
        except ValueError:
            raise argparse.ArgumentTypeError('point group {}: missing order'.format(group))
        if f == '':
            symbol = symmetry.PointGroupType.dihedral
        elif f == 'd':
            symbol = symmetry.PointGroupType.antiprismatic
        elif f == 'h':
            symbol = symmetry.PointGroupType.prismatic
        else:
            raise argparse.ArgumentTypeError('point group {}: wrong type'.format(group))
    elif m == 'T':
        if o != '':
            raise argparse.ArgumentTypeError('point group {}: order not needed'.format(group))
        if f == '':
            symbol = symmetry.PointGroupType.tetrahedral_chiral
        elif f == 'h':
            symbol = symmetry.PointGroupType.pyritohedral
        elif f == 'd':
            symbol = symmetry.PointGroupType.tetrahedral_achiral
        else:
            raise argparse.ArgumentTypeError('point group {}: wrong type'.format(group))
    elif m == 'O':
        if o != '':
            raise argparse.ArgumentTypeError('point group {}: order not needed'.format(group))
        if f == '':
            symbol = symmetry.PointGroupType.octahedral_chiral
        elif f == 'h':
            symbol = symmetry.PointGroupType.octahedral_achiral
        else:
            raise argparse.ArgumentTypeError('point group {}: wrong type'.format(group))
    elif m == 'I':
        if o != '':
            raise argparse.ArgumentTypeError('point group {}: order not needed'.format(group))
        if f == '':
            symbol = symmetry.PointGroupType.icosahedral_chiral
        elif f == 'h':
            symbol = symmetry.PointGroupType.icosahedral_achiral
        else:
            raise argparse.ArgumentTypeError('point group {}: wrong type'.format(group))
    else:
        raise argparse.ArgumentTypeError('point group {}: unknown type'.format(group))

    return symmetry.PointGroupDescription(symbol, order)


def get_arguments_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    parser.add_argument(
        'group',
        help='group for which character table must be printed',
        type=is_group)

    return parser


def main():
    parser = get_arguments_parser()
    args = parser.parse_args()

    group = args.group.gen_point_group(lowers_infinite=8)
    print(group.character_table())


if __name__ == '__main__':
    main()
