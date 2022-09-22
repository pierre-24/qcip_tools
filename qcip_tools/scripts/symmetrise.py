#!/usr/bin/env python3
"""
Detect point group, symmetrise molecule and outputs it as XYZ
"""

import argparse
import sys

import qcip_tools.scripts
from qcip_tools import molecule
from qcip_tools.chemistry_files import helpers, PropertyNotPresent, xyz


__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

__longdoc__ = """
Detect point group and symmetrise molecule.

The tolerance option (``--tolerance``) depends of the quality of an eventual optimisation: ``1e-3`` (the default) is ok
in most cases, but may results in a larger group than it should in a few cases, so lower this threshold if any.

By default, the script only orient the molecule w.r.t. symmetry elements, so that it matches its group.
But ``--symmetrise`` and/or ``--uniques`` implies to generate the whole group, which may take time for
large ones (thus very symmetrical molecules, with :math:`\\#\\mathcal{G} \\geq 60`). Then, it generates uniques atoms,
and eventually the whole molecule back (if only ``--symmetrise`` is set).

.. warning::

    ``--symmetrise`` and ``--uniques`` do not work for icosahedral molecules!!

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

    parser.add_argument(
        '-t', '--tolerance', help='Tolerance threshold', default=1e-3, type=float)

    parser.add_argument(
        '-s', '--symmetrise',
        help='Fully symmetrise the molecule (involve small displacements and the full generation of the group)',
        action='store_true')

    parser.add_argument(
        '-u', '--uniques', help='Get unique atoms (implies `-s`)', action='store_true')

    parser.add_argument(
        '-S', '--same-output', help='Output in the same format as the input', action='store_true')

    return parser


def main():
    parser = get_arguments_parser()
    args = parser.parse_args()

    try:
        mol = args.infile.property('molecule')
    except PropertyNotPresent:
        return qcip_tools.scripts.exit_failure('cannot find molecule ({})'.format(args.infile.file_type))

    sf = molecule.MolecularSymmetryFinder(mol, tol=args.tolerance)

    title = mol.formula()
    if args.symmetrise or args.uniques:
        g, uniques = sf.symmetrise_molecule(lowers_infinite=8)
        title += ' (point group {}'.format(g.description)
        if args.uniques:
            title += ', unique atoms'
            lst = [mol[i] for i in uniques]
            mol = molecule.Molecule(atom_list=lst)

        title += ')'
    else:
        description = sf.orient_molecule()
        title += ' (point group {})'.format(description)

    if not args.same_output:
        xf = xyz.File.from_molecule(mol, title=title)
        print(xf.to_string())
    else:
        args.infile.molecule = mol
        if hasattr(args.infile, 'title'):
            args.infile.title = title
        print(args.infile.to_string())


if __name__ == '__main__':
    main()
