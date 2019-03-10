import numpy
import random

import qcip_tools.math
from tests import QcipToolsTestCase
from qcip_tools import symmetry


class SymmetryTestCase(QcipToolsTestCase):

    def test_closest_fraction(self):
        self.assertEqual(qcip_tools.math.closest_fraction(-.499999, max_denominator=1000), (-1, 2))
        self.assertEqual(qcip_tools.math.closest_fraction(.499999, max_denominator=1000), (1, 2))

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
        self.assertArraysAlmostEqual(numpy.dot(E.matrix_representation(), p), p)

        i = symmetry.Operation.i()
        self.assertArraysAlmostEqual(i.apply(p), -p)
        self.assertArraysAlmostEqual(numpy.dot(i.matrix_representation(), p), -p)

        Oz = symmetry.Operation.sigma()
        self.assertArraysAlmostEqual(Oz.apply(p), numpy.array([*p[:2], -p[2]]))
        self.assertArraysAlmostEqual(numpy.dot(Oz.matrix_representation(), p), numpy.array([*p[:2], -p[2]]))

        S52z = symmetry.Operation.S(5, 2, axis=numpy.array([1, 1., 0.]))
        self.assertArraysAlmostEqual(S52z.apply(p), numpy.dot(S52z.matrix_representation(), p))

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
            (symmetry.PointGroup.C_n(4), 4, True, 4),
            (symmetry.PointGroup.S_n(2), 2, True, 2),  # = C_i
            (symmetry.PointGroup.S_n(4), 4, True, 4),
            (symmetry.PointGroup.C_nv(2), 4, True, 4),
            (symmetry.PointGroup.C_nv(3), 6, False, 3),
            (symmetry.PointGroup.C_nv(4), 8, False, 5),
            (symmetry.PointGroup.C_nh(1), 2, True, 2),  # = C_s
            (symmetry.PointGroup.C_nh(2), 4, True, 4),  # != S_2
            (symmetry.PointGroup.C_nh(3), 6, True, 6),
            (symmetry.PointGroup.C_nh(4), 8, True, 8),  # != S_4
            (symmetry.PointGroup.D_n(4), 8, False, 5),
            (symmetry.PointGroup.D_nh(3), 12, False, 6),
            (symmetry.PointGroup.D_nh(4), 16, False, 10),
            (symmetry.PointGroup.D_nd(2), 8, False, 5),
            (symmetry.PointGroup.D_nd(3), 12, False, 6),
            (symmetry.PointGroup.T(), 12, False, 4),
            (symmetry.PointGroup.T_h(), 24, False, 8),
            (symmetry.PointGroup.T_d(), 24, False, 5),
            (symmetry.PointGroup.O(), 24, False, 5),
            (symmetry.PointGroup.O_h(), 48, False, 10),
            (symmetry.PointGroup.I(), 60, False, 5),  # ~2 seconds
            (symmetry.PointGroup.I_h(), 120, False, 10),  # ~ 5 seconds
        ]

        for g, n_elements, is_abelian, number_of_class in groups:
            self.assertEqual(len(g.G), n_elements)
            self.assertEqual(g.abelian, is_abelian)
            self.assertEqual(g.number_of_class, number_of_class)
            if len(g.G) <= 12:
                self.assertTrue(g.binary_operation.check_surjectivity())  # O(n²)
                self.assertTrue(g.binary_operation.check_associativity())  # O(n³)

    def test_conjugacy_classes(self):
        """Tests if the conjugacy classes comes in order"""

        C_3v = symmetry.PointGroup.C_nv(3)
        self.assertTrue(
            list(len(x) for x in C_3v.conjugacy_classes) == [1, 2, 3])  # ... Which gives 6 elements in total

        self.assertEqual(C_3v.conjugacy_classes[0], {C_3v.e})  # first class is always identity
        self.assertEqual(
            {e.element for e in C_3v.conjugacy_classes[1]}, {symmetry.Operation.C(3), symmetry.Operation.C(3, 2)})

    def test_character_table(self):
        """Generate the character table"""

        # a few groups to check
        groups = [
            symmetry.PointGroup.C_nv(2),
            symmetry.PointGroup.C_nh(3),
            symmetry.PointGroup.C_nh(4),
            symmetry.PointGroup.D_nd(3),
            symmetry.PointGroup.T_h(),
            symmetry.PointGroup.I()
        ]

        for g in groups:
            t = g.character_table()
            print(t)
            # print('----')

            # check orthogonality of the lines
            for i in range(g.number_of_class):
                # print(t.irreducible_representations[i])
                for j in range(i + 1, g.number_of_class):
                    self.assertAlmostEqual(t.irreducible_representations[i].dot(t.irreducible_representations[j]), .0)

    def test_representations(self):
        """Check the direct representation sum and product

        Checked against http://gernot-katzers-spice-pages.com/character_tables/C6v.html
        """

        g = symmetry.PointGroup.C_nv(6)
        t = g.character_table()

        # linear characters
        A1 = t.irreducible_representations[0]  # A1
        A2 = t.irreducible_representations[1]  # A2
        B1 = t.irreducible_representations[2]  # B1
        B2 = t.irreducible_representations[3]  # B2
        # second order ones:
        E1 = t.irreducible_representations[4]
        E2 = t.irreducible_representations[5]

        # test direct sum:
        r = A1 + A2 + E2
        self.assertEqual(r.size, 3)

        self.assertEqual(r.irreducible_representations[0], 1)
        self.assertIn(A1, r)
        self.assertEqual(r.irreducible_representations[1], 1)
        self.assertIn(A2, r)
        self.assertEqual(r.irreducible_representations[2], 0)
        self.assertNotIn(B1, r)
        self.assertEqual(r.irreducible_representations[3], 0)
        self.assertNotIn(B2, r)
        self.assertEqual(r.irreducible_representations[4], 0)
        self.assertNotIn(E1, r)
        self.assertEqual(r.irreducible_representations[5], 1)
        self.assertIn(E2, r)

        # test direct product
        for ri in [A1, A2, B1, B2, E1, E2]:
            r = A1 * ri  # trivial product gives the same irrep
            self.assertEqual(r.size, 1)
            self.assertIn(ri, r)

        r = E1 * B2  # = E2
        self.assertEqual(r.size, 1)
        self.assertIn(E2, r)

        r = E1 * E2  # = B1 + B2 + E1
        self.assertEqual(r.size, 3)
        self.assertNotIn(A1, r)
        self.assertNotIn(A2, r)
        self.assertIn(B1, r)
        self.assertIn(B2, r)
        self.assertIn(E1, r)
        self.assertNotIn(E2, r)

        r = E1 * E1  # = A1 + A2 + E2
        self.assertEqual(r.size, 3)
        self.assertIn(A1, r)
        self.assertIn(A2, r)
        self.assertNotIn(B1, r)
        self.assertNotIn(B2, r)
        self.assertNotIn(E1, r)
        self.assertIn(E2, r)

    def test_symmetry_finder(self):
        """Test if one is able to detect symmetry"""

        def s(points, tol=1e-5):
            d = symmetry.SymmetryFinder(points, tol).find_symmetry()[0]
            return d.symbol, d.order

        """O_h geometry:

           1 1
           |/
        1--3--1
          /|
         1 1
        """

        p = numpy.array([
            (3., 0, 0, 0),
            (1., 1, 0, 0),
            (1., -1, 0, 0),
            (1., 0, 1, 0),
            (1., 0, -1, 0),
            (1., 0, 0, 1),
            (1., 0, 0, -1),
        ])

        self.assertEqual(s(p), (symmetry.PointGroupType.octahedral_achiral, 0))
        self.assertEqual(s(p[:-1]), (symmetry.PointGroupType.pyramidal, 4))
        self.assertEqual(s(p[:-2]), (symmetry.PointGroupType.prismatic, 4))
        self.assertEqual(s(p[:-3]), (symmetry.PointGroupType.pyramidal, 2))
        self.assertEqual(s(p[:-4]), (symmetry.PointGroupType.prismatic, -1))
        self.assertEqual(s(p[:-5]), (symmetry.PointGroupType.pyramidal, -1))
