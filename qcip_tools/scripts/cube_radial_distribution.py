#!/usr/bin/env python3
"""
Radial distribution
"""

import argparse
import math
import numpy
import itertools
from typing import Tuple
from qcip_tools import quantities
from qcip_tools.chemistry_files import gaussian, helpers

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

__longdoc__ = """
Report the radial distribution of a cube around a given center [by default :math:`(0,0,0)`].

The charge in a given region of the space, located by :math:`\\mathbf{r}` and in an element of volume
:math:`d\\mathbf{r}`, is given by

.. math::

    q(\\mathbf{r}) = \\rho(\\mathbf{r})\\,d\\mathbf{r}.

Integration over whole space gives the number of particles, :math:`Q`.
In spherical coordinates, :math:`d\\mathbf{r} = r^2\\sin{\\theta}\\,dr\\,d\\theta\\,d\\phi`,
this integral becomes

.. math::

    Q =
    \\int_0^{2\\pi}\\int_0^{\\pi}\\int_0^{\\infty}
    \\rho(r,\\theta,\\phi)\\,r^2\\,\\sin{\\theta}\\,dr\\,d\\theta\\,d\\phi.

Thus, the radial distribution is given by

.. math::
    :label: tr

    \\frac{dQ(r)}{dr} =
    r^2\\,\\int_0^{2\\pi}\\int_0^{\\pi}
    \\rho(r,\\theta,\\phi)\\sin{\\theta}\\,d\\theta\\,d\\phi,

Equation :eq:`tr` is obtained numerically by interpolation over the cube.
"""


def trilinear_interpolation(cube: gaussian.Cube, c: Tuple[float, float, float], data=0):
    """Interpolate, assuming the sides have a length of 1.
    See https://en.wikipedia.org/wiki/Trilinear_interpolation.

    ``c`` is in cube coordinates.

    Note:

    ```
       (011)        (111)
       o----------o
      /|         /|
     / |        / |
    o-(010)-----o (101)
    |  |       |  |
    |  o-------|--o (110)
    | / (010)  | /
    |/ (000)   |/
    o----------o (100)

    with,

    y z
    |/
    o-- x

    ```
    """

    def li(a, b, x):
        return a * (1 - x) + b * x

    ce = list(int(ci) for ci in c)

    rc = numpy.zeros((2, 2, 2))
    for i in itertools.product((0, 1), repeat=3):
        coo = tuple(ce[j] + i[j] for j in range(3))
        if all(0 <= coo[j] < cube.records_per_direction[j] for j in range(3)):
            rc[i] = cube.records[coo][data]

    # interpolate along x:
    xd = c[0] - int(c[0])
    c00 = li(rc[0, 0, 0], rc[1, 0, 0], xd)
    c01 = li(rc[0, 0, 1], rc[1, 0, 1], xd)
    c10 = li(rc[0, 1, 0], rc[1, 1, 0], xd)
    c11 = li(rc[0, 1, 1], rc[1, 1, 1], xd)

    # interpolate along y:
    yd = c[1] - int(c[1])
    c0 = li(c00, c10, yd)
    c1 = li(c01, c11, yd)

    # interpolate along z:
    return li(c0, c1, c[2] - int(c[2]))


def to_cartesian(r, theta, phi):
    st = math.sin(theta)
    ct = math.cos(theta)
    sp = math.sin(phi)
    cp = math.cos(phi)
    return r * st * cp, r * st * sp, r * ct


def radial_distribution(cube: gaussian.Cube, center, max_radius, dr, data=0, npolar=20, nazim=40):
    """Note: everything in cube dimension (C.D.)
    """

    dtheta = math.pi / npolar
    dphi = 2 * math.pi / nazim

    radial_distrib = numpy.zeros((math.ceil(max_radius / dr), 3))

    r = .0
    ri = 0
    while r < max_radius:
        sum_dR = 0

        for t in range(npolar + 1):
            for p in range(nazim + 1):
                coo = to_cartesian(r, t * dtheta, p * dphi)
                val = trilinear_interpolation(cube, tuple(center[i] + coo[i] for i in range(3)), data=data)
                sum_dR += val * r ** 2 * math.sin(t * dtheta) * dtheta * dphi

        radial_distrib[ri] = [r, sum_dR * dr, sum_dR]
        ri += 1
        r = ri * dr

    return radial_distrib


def is_coordinate(coo):
    if coo[0] == ':':
        coo = coo[1:]

    e = coo.split(',')
    if len(e) != 3:
        raise argparse.ArgumentError(message='Coordinates must be a triplet of floats')

    try:
        return tuple(float(x) for x in e)
    except ValueError:
        raise argparse.ArgumentError(message='coordinate must contain float')


def get_arguments_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    parser.add_argument(
        'infile', action=helpers.create_open_chemistry_file_action(must_be=[gaussian.Cube]),
        help='Density cube')

    parser.add_argument('-S', '--square', action='store_true', help='square cube before')

    parser.add_argument(
        '-c', '--center', help='Center', type=is_coordinate, default='0,0,0')

    parser.add_argument('--dr', help='Increment of radius', type=float, default=.05)

    parser.add_argument('--n-polar', help='Number of subdivision for the integration over theta', type=int, default=20)
    parser.add_argument(
        '--n-azimuthal', help='Number of subdivision for the integration over phi', type=int, default=40)

    parser.add_argument(
        '-d', '--data', help='Center', type=int, default=0)

    return parser


def main():
    parser = get_arguments_parser()
    args = parser.parse_args()

    cube: gaussian.Cube = args.infile
    incr = cube.increments[0]

    if args.square:
        cube = cube ** 2

    if cube.data_per_record <= args.data:
        raise Exception('only {} data per record in cube'.format(cube.data_per_record))

    # convert center in bohr, then get the result in cube dimensions (C.D.)
    center = tuple(
        (x / quantities.AuToAngstrom - cube.origin[i]) / cube.increments[i] for i, x in enumerate(args.center))

    print('# center (in C.D.): {:.3f}, {:.3f}, {:.3f}'.format(*center))

    # radius
    radius = math.sqrt(max(
        sum(center[i] ** 2 for i in range(3)),
        sum(cube.records_per_direction[i] - center[i] for i in range(3))))

    print('# maximal radius: {:.3f} Å  ({:.3f} C.D.)'.format(radius * incr * quantities.AuToAngstrom, radius))

    # dr
    dr = args.dr / quantities.AuToAngstrom / incr
    print('# dr: {} Å ({:.3f} C.D.)'.format(args.dr, dr))

    # data
    if args.data:
        data = [args.data]
    else:
        data = range(cube.data_per_record)

    # compute, then outputs
    for index, d in enumerate(data):
        if index > 0:
            print('\n\n')  # change index

        print('# --- working on data #{}'.format(d))
        print('# sum = {:.3f}'.format(cube.records[:, :, :, d].sum() * numpy.prod(cube.increments)))

        distribution = radial_distribution(
            cube, center, radius, dr=dr, data=d, npolar=args.n_polar, nazim=args.n_azimuthal)

        print('# r\tQ(r)\tdQ(r)/dr')
        for i in range(distribution.shape[0]):
            print('{:.3f}\t{:.5e}\t{:.5e}'.format(
                distribution[i][0] * incr * quantities.AuToAngstrom,
                distribution[i][1] * incr ** 3,
                distribution[i][2] * incr ** 2 / quantities.AuToAngstrom
            ))

        qct = distribution[:, 2].sum() * dr * incr ** 3
        print('# sums\t{:.4f}\t{:.4f}'.format(distribution[:, 1].sum() * incr ** 3, qct))
        print('# mean radii\t\t{:.4f}'.format(
            (distribution[:, 0] * distribution[:, 2]).sum() * dr * incr ** 4 * quantities.AuToAngstrom / qct
        ))


if __name__ == '__main__':
    main()
