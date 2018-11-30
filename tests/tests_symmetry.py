import sympy
import numpy
import random

from tests import QcipToolsTestCase
from qcip_tools import symmetry


class SymmetryTestCase(QcipToolsTestCase):

    def test_group(self):
        f = symmetry.BinaryOperation(symmetry.Set(range(4)), lambda e: (e[0] + e[1]) % 4)
        self.assertTrue(f.check_surjectivity())
        self.assertTrue(f.check_associativity())

        g = symmetry.Group(f)  # just a cyclic group of order 4, then
        self.assertEqual(g.identity(), 0)
        self.assertEqual(g.inverse(0), 0)
        self.assertEqual(g.inverse(1), 3)
        self.assertEqual(g.inverse(2), 2)

        self.assertTrue(g.abelian)

        for i in range(4):
            self.assertEqual(g.inverse(i) ** -1, i)

        n3 = g.G[3]
        self.assertEqual(n3, 3)
        self.assertEqual(n3 ** 1, 3)
        self.assertEqual(n3 ** 2, 2)
        self.assertEqual(n3 ** 3, 1)
        self.assertEqual(n3 ** -1, 1)
        self.assertEqual(n3 ** -2, 2)

    def test_symmetry_element(self):
        p = numpy.array([random.randint(-3, 3), random.randint(-3, 3), random.randint(-3, 3)])

        # basic elements
        E = symmetry.SymElement.E()
        self.assertArraysAlmostEqual(E.apply(p), p)

        i = symmetry.SymElement.i()
        self.assertArraysAlmostEqual(i.apply(p), -p)

        Oz = symmetry.SymElement.sigma()
        self.assertArraysAlmostEqual(Oz.apply(p), numpy.array([*p[:2], -p[2]]))

        # test identity
        t = Oz * Oz
        self.assertArraysAlmostEqual(t.apply(p), p)
        t = i * i
        self.assertArraysAlmostEqual(t.apply(p), p)

        # Check that C_3 is cyclic
        C3z = symmetry.SymElement.C(3)
        C23z = symmetry.SymElement.C(3, 2)

        t = C3z * C3z
        self.assertArraysAlmostEqual(t.apply(p), C23z.apply(p))

        t = t * C3z
        self.assertArraysAlmostEqual(t.apply(p), p)

        # Check that S_3 is cyclic
        S3z = symmetry.SymElement.S(3)
        S53z = symmetry.SymElement.S(3, 5)

        t = S3z * S3z
        self.assertArraysAlmostEqual(t.apply(p), C23z.apply(p))

        t = t * S3z
        self.assertArraysAlmostEqual(t.apply(p), Oz.apply(p))

        t = t * S3z
        self.assertArraysAlmostEqual(t.apply(p), C3z.apply(p))

        t = t * S3z
        self.assertArraysAlmostEqual(t.apply(p), S53z.apply(p))

        t = t * S3z
        self.assertArraysAlmostEqual(t.apply(p), p)

        # Check composition: C_n(z) * O(z) = S_n(z)
        n = 5
        Cnz = symmetry.SymElement.C(n)
        Snz = symmetry.SymElement.S(n)
        t = Cnz * Oz
        self.assertArraysAlmostEqual(t.apply(p), Snz.apply(p))

        # Check composition: C_2(y) * C_2(x) = C_2(z)
        C2z = symmetry.SymElement.C(2, axis=numpy.array([0, 0, 1.]))
        C2y = symmetry.SymElement.C(2, axis=numpy.array([0, 1., 0]))
        C2x = symmetry.SymElement.C(2, axis=numpy.array([1., 0, 0]))
        t = C2y * C2x
        self.assertArraysAlmostEqual(t.apply(p), C2z.apply(p))

        # Check composition: O(y) * O(x) = C_2(z)
        Ox = symmetry.SymElement.sigma(axis=numpy.array([1, 0., 0]))
        Oy = symmetry.SymElement.sigma(axis=numpy.array([0, 1., 0]))
        t = Oy * Ox
        self.assertArraysAlmostEqual(t.apply(p), C2z.apply(p))

        # Check composition: O(y) * O(xy') = C_4(z)
        Oxy = symmetry.SymElement.sigma(axis=numpy.array([1., -1., 0]))
        C4z = symmetry.SymElement.C(4)
        t = Oy * Oxy
        self.assertArraysAlmostEqual(t.apply(p), C4z.apply(p))

    def test_random(self):
        """That test will disappear in the future
        """

        class Quat:
            """
            NT: according to https://www.mdpi.com/2073-8994/2/3/1423/pdf,
            improper rotation are defined as S_m = R*_n(pi + 2pi/m),
            where R* is the rotation followed by an inversion
            (so the negation of the whole product → simply the conjugation of the point quaternion).
            """

            @classmethod
            def from_axangle(cls, axe, angle=2 * sympy.pi, improper=False):
                q = sympy.Quaternion.from_axis_angle(axe, angle)
                return cls(q, improper)

            def __init__(self, q, improper=False):
                self.q = q
                self.improper = improper

            def apply(self, pin):
                qt = sympy.Quaternion(0, *pin)
                if self.improper:
                    qt = qt.conjugate()

                qe = self.q.mul(qt.mul(self.q.conjugate()))
                return qe.b, qe.c, qe.d

            def __mul__(self, other):
                if not isinstance(other, Quat):
                    raise TypeError('other must be Quat')

                if other.q != self.q:
                    qa = self.q
                    qb = other.q

                    q = qb.mul(qa)
                else:
                    q = self.q ** 2

                return Quat(
                    q,
                    improper=(self.improper and not other.improper) or (other.improper and not self.improper))

        p = [3, -2, 1]

        E = Quat.from_axangle([0, 0, 1])
        self.assertEqual(E.apply(p), tuple(p))

        Oz = Quat.from_axangle([0, 0, 1], angle=sympy.pi + 2 * sympy.pi, improper=True)
        self.assertEqual(Oz.apply(p)[:2], tuple(p[:2]))
        self.assertEqual(Oz.apply(p)[2], -p[2])

        i = Quat.from_axangle([0, 0, 1], angle=sympy.pi + sympy.pi, improper=True)
        self.assertEqual(i.apply(p), (-p[0], -p[1], -p[2]))
