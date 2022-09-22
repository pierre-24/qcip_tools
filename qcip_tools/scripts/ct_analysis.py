#!/usr/bin/env python3
"""
Compute charge transfer (CT) quantities
"""

import argparse
import numpy
from qcip_tools import atom, derivatives_e
from qcip_tools.chemistry_files import gaussian, helpers, xyz

__version__ = '0.2'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

__longdoc__ = """
Theory
******

Based on an implementation of this theory by D. Jacquemin (see supporting informations of the corresponding paper).

The charge in a given region of the space, located by :math:`\\mathbf{r}` and in an element of volume
:math:`d\\mathbf{r}`, is given by

.. math::

    q(\\mathbf{r}) = \\int  \\rho(\\mathbf{r})\\,d\\mathbf{r}.

Charges at different point of the space (a "density") can be obtained by using the
``cubegen`` utility `provided by Gaussian <http://gaussian.com/cubegen/>`_. In particular,

+ ``cubegen 0 density=scf x.fchk out.cub`` permit to obtain the "density" of the ground state, and
+ ``cubegen 0 density=ci x.fchk out.cub`` permit to obtain the "density" of the excited state.

Note that you need to launch the Gaussian calculations with ``Density=(all)`` so that it stores the CI density
(if you use solvatation, please **make sure** to use ``TD=(NonEqSolv)``!).

Let :math:`\\delta q(r)` be the difference between the charge of the excited and the ground state.
This quantity can be splitted into increasing areas and decreasing ones, :math:`\\delta q_+(\\mathbf{r})` and
:math:`\\delta q_-(\\mathbf{r})`, where

.. math::

    \\delta q_+(\\mathbf{r}) = \\left\\{
    \\begin{array}{ll}
        \\delta q(\\mathbf{r}) & \\text{if }\\delta q(\\mathbf{r}) > 0, \\\\
        0 & \\text{otherwise.}
    \\end{array}
    \\right.

and,

.. math::

    \\delta q_-(\\mathbf{r}) = \\left\\{
    \\begin{array}{ll}
        \\delta q(\\mathbf{r}) & \\text{if }\\delta q(\\mathbf{r}) < 0, \\\\
        0 & \\text{otherwise.}
    \\end{array}
    \\right.

One can therefore compute:

- The transferred charge between ground and excited state:

  .. math::

    q_{CT} = \\frac{1}{2}\\,\\sum_{\\mathbf{r}_i} \\delta q_+(\\mathbf{r}_i) - \\delta q_-(\\mathbf{r}_i).

  Note that the original papers states that

  .. math::

    q_{CT} =\\sum_{\\mathbf{r}_i} \\delta q_+(\\mathbf{r}_i) =-\\sum_{\\mathbf{r}_i} \\delta q_-(\\mathbf{r}_i),

  but the implementation of D. Jacquemin reports and uses the average.

- The starting and ending point of the charge transfer, :math:`\\mathbf{r}_+` and :math:`\\mathbf{r}_-`:

  .. math::

    \\mathbf{r}_+ = \\sum_{\\mathbf{r}_i} \\frac{\\mathbf{r}_i\\,q_+(\\mathbf{r}_i)}{q_{CT}},

    \\mathbf{r}_- = \\sum_{\\mathbf{r}_i} \\frac{\\mathbf{r}_i\\,q_-(\\mathbf{r}_i)}{q_{CT}}.

  Those are the barycenters of the positive and the negative densities.
  The vector between those two barycenter is the charge transfer vector, defined as

  .. math::

    \\mathbf{v}_{CT} = \\mathbf{r}_--\\mathbf{r}_+.

  In particular, the charge transfer distance is the norm of this vector, :math:`d_{CT} = |\\mathbf{v}_{CT}|`.
  Notice the usage of the so called *chemist convention*, where the dipole is defined from positive to negative
  positions.

+ The norm of variation of dipole moment between the ground and excited state:

  .. math::

    |\\mu_{CT}| = q_{CT}\\,d_{CT}.


Implementation
**************

Only works with gaussian cube.
Note that external programs may be abble to generate those as well
(if this is not a density but a probability, like with MO, squaring it gives the density, so use the ``-S`` option).

The program reports :math:`\\mathbf{v}_{CT}` as well as :math:`q_{CT}`, :math:`d_{CT}`, and
:math:`|\\mu_{CT}|`.

It allows to save the difference cube (for visualization) and an xyz file containing two dummy atoms (one for each
barycenter, first :math:`\\mathbf{r}_-` and then :math:`\\mathbf{r}_+`).

.. warning::

    Vector and :math:`d_{CT}` are given in Angstrom, :math:`q_{CT}` is in \\|e\\| (electron charge),
    and :math:`|\\mu_{CT}|` is therefore in Angstrom \\|e\\|.

Source
******

+ T. Le Bahers *et al.* *J. Chem. Theory. Comput.* **7**, 2498 (2011)
  `10.1021/ct200308m <http://dx.doi.org/10.1021/ct200308m>`_.
+ D. Jacquemin *et al.* *Phys Chem Chem Phys.* **28**, 5383 (2012)
  `10.1039/c2cp40261k <http://dx.doi.org/10.1039/c2cp40261k>`_.
"""


def get_arguments_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    parser.add_argument(
        '-g', '--ground', required=True, action=helpers.create_open_chemistry_file_action(must_be=[gaussian.Cube]),
        help='Ground state density')
    parser.add_argument(
        '-e', '--excited', required=True, action=helpers.create_open_chemistry_file_action(must_be=[gaussian.Cube]),
        help='Excited state density')

    parser.add_argument('-S', '--square', action='store_true', help='square cube before make the difference')

    parser.add_argument(
        '-o', '--output', type=argparse.FileType('w'), help='difference cube')

    parser.add_argument(
        '-D', '--output-with-diff',
        action='store_true',
        help='Store all the densities (total, positive, negative)')

    parser.add_argument(
        '-d', '--dummy', type=argparse.FileType('w'),
        help='Create an xyz file with dummy atoms in positions of both barycenters (negative then positive)')

    return parser


def main():
    parser = get_arguments_parser()
    args = parser.parse_args()

    if args.square:
        args.excited = args.excited ** 2
        args.ground = args.ground ** 2

    diff: gaussian.Cube = args.excited - args.ground
    ct = diff.compute_charge_transfer()

    dipole = derivatives_e.ElectricDipole(dipole=ct.vector)

    print('Barycenters and charge transfer vector (angstrom):')
    print('            x            y             z')
    print('c(m) {: .6e} {: .6e} {: .6e}'.format(*ct.barycenter_m))
    print('c(p) {: .6e} {: .6e} {: .6e}'.format(*ct.barycenter_p))
    print('v    {: .6e} {: .6e} {: .6e}'.format(*dipole.components))
    print()
    print('q_CT (|e|)     {: .5e}'.format(ct.charge))
    print('d_CT (ang)     {: .5e}'.format(dipole.norm()))
    print('µ_CT (|e|*ang) {: .5e}'.format(ct.charge * ct.distance))

    if args.output:
        if args.output_with_diff:
            records = diff.records[:, :, :, 0].copy()
            rho_p = records.copy()
            rho_p[numpy.where(rho_p < 0.0)] = 0.0
            rho_m = records.copy()
            diff.records = numpy.zeros((*diff.records_per_direction, 3))
            rho_m[numpy.where(rho_p > 0.0)] = 0.0

            diff.records[:, :, :, 0] = records
            diff.records[:, :, :, 1] = rho_p
            diff.records[:, :, :, 2] = rho_m

            diff.data_per_record = 3

        diff.write(args.output)

    if args.dummy:
        fx = xyz.File.from_molecule(
            diff.molecule, title='generated from CT (q={:.3f}, d={:.3f})'.format(ct.charge, ct.distance))
        fx.molecule.insert(atom.Atom(symbol='Mo', position=ct.barycenter_m))
        fx.molecule.insert(atom.Atom(symbol='Po', position=ct.barycenter_p))

        fx.write(args.dummy)


if __name__ == '__main__':
    main()
