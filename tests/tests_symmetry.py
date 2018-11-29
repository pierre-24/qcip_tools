import sympy

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
        s = symmetry.SymmetryElement.S(4, 7)
        print(s)

    def test_random(self):
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

                qa = self.q
                if self.improper:
                    qa = qa.conjugate()
                qb = other.q
                if other.improper:
                    qb = qb.conjugate()

                return Quat(
                    qb.mul(qa),
                    improper=self.improper and not other.improper or other.improper and not self.improper)

        p = [4, 2, 3]

        # Check C_n(z) * C_n(z) = C^2_n(z)
        C3z = Quat.from_axangle([0, 0, 1], 2 * sympy.pi / 3)
        C23z = Quat.from_axangle([0, 0, 1], 2 * sympy.pi / 3 * 2)
        t = C3z * C3z
        self.assertEqual(t.q, C23z.q)
        self.assertEqual(t.apply(p), C23z.apply(p))

        # S3z = Quat.from_axangle([0, 0, 1], sympy.pi + 2 * sympy.pi / 3, improper=True)
        # S23z = Quat.from_axangle([0, 0, 1], 2 * sympy.pi / 3 * 2)
        # t = S3z * S3z
        # self.assertEqual(t.q, S23z.q)
        # self.assertEqual(t.apply(p), S23z.apply(p))

        # check composition of C_n(z) * O(z) → S_n(z)
        C4z = Quat.from_axangle([0, 0, 1], sympy.pi / 2)
        Oz = Quat.from_axangle([0, 0, 1], angle=sympy.pi + 2 * sympy.pi, improper=True)
        S4z = Quat.from_axangle([0, 0, 1], angle=sympy.pi + 2 * sympy.pi / 4, improper=True)
        t = C4z * Oz
        self.assertEqual(t.q, S4z.q)
        self.assertEqual(t.apply(p), S4z.apply(p))

        C5z = Quat.from_axangle([0, 0, 1], 2 * sympy.pi / 5)
        S5z = Quat.from_axangle([0, 0, 1], sympy.pi + 2 * sympy.pi / 5, improper=True)
        t = C5z * Oz
        self.assertEqual(t.q, S5z.q)
        self.assertEqual(t.apply(p), S5z.apply(p))

        # check C_2(z) * C_2(x) → C2z
        C2z = Quat.from_axangle([0, 0, 1], angle=sympy.pi)
        C2y = Quat.from_axangle([0, 1, 0], angle=sympy.pi)
        C2x = Quat.from_axangle([1, 0, 0], angle=sympy.pi)
        t = C2y * C2x
        self.assertEqual(t.q, C2z.q)
        self.assertEqual(t.apply(p), C2z.apply(p))

        # check Ox * Oy → C2z
        Ox = Quat.from_axangle([1, 0, 0], angle=sympy.pi + 2 * sympy.pi, improper=True)
        Oy = Quat.from_axangle([0, 1, 0], angle=sympy.pi + 2 * sympy.pi, improper=True)
        t = Oy * Ox
        self.assertEqual(t.q, C2z.q)
        self.assertEqual(t.apply(p), C2z.apply(p))
