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

        n3 = symmetry.GroupElement(3, g)
        self.assertEqual(n3, 3)
        self.assertEqual(n3 ** 1, 3)
        self.assertEqual(n3 ** 2, 2)
        self.assertEqual(n3 ** 3, 1)
        self.assertEqual(n3 ** -1, 1)
        self.assertEqual(n3 ** -2, 2)

    def test_symmetry_element(self):
        p = numpy.array([random.randint(-3, 3), random.randint(-3, 3), random.randint(-3, 3)])

        # basic elements
        E = symmetry.Operation.E()
        self.assertArraysAlmostEqual(E.apply(p), p)

        i = symmetry.Operation.i()
        self.assertArraysAlmostEqual(i.apply(p), -p)

        Oz = symmetry.Operation.sigma()
        self.assertArraysAlmostEqual(Oz.apply(p), numpy.array([*p[:2], -p[2]]))

        # test identity
        t = Oz * Oz
        self.assertEqual(t, E)
        self.assertArraysAlmostEqual(t.apply(p), p)
        t = i * i
        self.assertEqual(t, E)
        self.assertArraysAlmostEqual(t.apply(p), p)

        # Check that C_3 is cyclic
        C3z = symmetry.Operation.C(3)
        C23z = symmetry.Operation.C(3, 2)

        t = C3z * C3z
        self.assertEqual(t, C23z)
        self.assertEqual(t.get_description(), C23z.get_description())
        self.assertArraysAlmostEqual(t.apply(p), C23z.apply(p))

        t = t * C3z
        self.assertEqual(t, E)
        self.assertArraysAlmostEqual(t.apply(p), p)

        # Check that S_3 is cyclic
        S3z = symmetry.Operation.S(3)
        S53z = symmetry.Operation.S(3, 5)

        t = S3z * S3z
        self.assertEqual(t, C23z)
        self.assertEqual(t.get_description(), C23z.get_description())
        self.assertArraysAlmostEqual(t.apply(p), C23z.apply(p))

        t = t * S3z
        self.assertEqual(t, Oz)
        # TODO: self.assertEqual(t.get_description(), Oz.get_description()) → k = 2
        self.assertArraysAlmostEqual(t.apply(p), Oz.apply(p))

        t = t * S3z
        self.assertEqual(t, C3z)
        self.assertEqual(t.get_description(), C3z.get_description())
        self.assertArraysAlmostEqual(t.apply(p), C3z.apply(p))

        t = t * S3z
        self.assertEqual(t, S53z)
        self.assertEqual(t.get_description(), S53z.get_description())
        self.assertArraysAlmostEqual(t.apply(p), S53z.apply(p))

        t = t * S3z
        self.assertEqual(t, E)
        # TODO: self.assertEqual(t.get_description(), E.get_description()) → all that matter is the symbol
        self.assertArraysAlmostEqual(t.apply(p), p)

        # Check composition: C_n(z) * O(z) = S_n(z)
        n = 5
        Cnz = symmetry.Operation.C(n)
        Snz = symmetry.Operation.S(n)
        t = Cnz * Oz
        self.assertEqual(t, Snz)
        self.assertEqual(t.get_description(), Snz.get_description())
        self.assertArraysAlmostEqual(t.apply(p), Snz.apply(p))

        # Check composition: C_2(y) * C_2(x) = C_2(z)
        C2z = symmetry.Operation.C(2, axis=numpy.array([0, 0, 1.]))
        C2y = symmetry.Operation.C(2, axis=numpy.array([0, 1., 0]))
        C2x = symmetry.Operation.C(2, axis=numpy.array([1., 0, 0]))
        t = C2y * C2x
        self.assertEqual(t, C2z)
        self.assertEqual(t.get_description(), C2z.get_description())
        self.assertArraysAlmostEqual(t.apply(p), C2z.apply(p))

        # Check composition: O(y) * O(x) = C_2(z)
        Ox = symmetry.Operation.sigma(axis=numpy.array([1, 0., 0]))
        Oy = symmetry.Operation.sigma(axis=numpy.array([0, 1., 0]))
        t = Oy * Ox
        self.assertEqual(t, C2z)
        self.assertEqual(t.get_description(), C2z.get_description())
        self.assertArraysAlmostEqual(t.apply(p), C2z.apply(p))

        # Check composition: O(y) * O(xy') = C_4(z)
        Oxy = symmetry.Operation.sigma(axis=numpy.array([1., -1., 0]))
        C4z = symmetry.Operation.C(4)
        t = Oy * Oxy
        self.assertEqual(t, C4z)
        self.assertEqual(t.get_description(), C4z.get_description())
        self.assertArraysAlmostEqual(t.apply(p), C4z.apply(p))

    def test_point_groups(self):

        # checked against http://gernot-katzers-spice-pages.com/character_tables/index.html
        groups = [
            (symmetry.PointGroup.C_n(4), 4),
            (symmetry.PointGroup.C_nv(2), 4),
            (symmetry.PointGroup.C_nv(3), 6),
            (symmetry.PointGroup.C_nh(1), 2),  # = C_s
            (symmetry.PointGroup.C_nh(2), 4),  # != S_2
            (symmetry.PointGroup.C_nh(3), 6),
            (symmetry.PointGroup.C_nh(4), 8),  # != S_4
            (symmetry.PointGroup.S_n(2), 2),  # = C_i
            (symmetry.PointGroup.S_n(4), 4),
        ]

        for g, n_elements in groups:
            print(g)
            self.assertEqual(len(g.G), n_elements)
            self.assertTrue(g.binary_operation.check_surjectivity())
            self.assertTrue(g.binary_operation.check_associativity())
