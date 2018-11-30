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

                if other.q != self.q:
                    qa = self.q
                    qb = other.q

                    q = qb.mul(qa)
                else:
                    q = self.q ** 2

                return Quat(
                    q,
                    improper=(self.improper and not other.improper) or (other.improper and not self.improper))

        p = [-1, -2, -1]

        E = Quat.from_axangle([0, 0, 1])
        self.assertEqual(E.apply(p), tuple(p))

        Oz = Quat.from_axangle([0, 0, 1], angle=sympy.pi + 2 * sympy.pi, improper=True)
        self.assertEqual(Oz.apply(p)[:2], tuple(p[:2]))
        self.assertEqual(Oz.apply(p)[2], -p[2])

        # Check identity for Oz * Oz
        t = Oz * Oz
        self.assertEqual(t.apply(p), E.apply(p))

        # Check C_n(z) cyclicity
        C3z = Quat.from_axangle([0, 0, 1], 2 * sympy.pi / 3)
        C23z = Quat.from_axangle([0, 0, 1], 2 * sympy.pi / 3 * 2)
        t = C3z * C3z
        self.assertEqual(t.q, C23z.q)
        self.assertEqual(t.apply(p), C23z.apply(p))

        t = t * C3z
        self.assertEqual(E.q, t.q)
        self.assertEqual(t.apply(p), E.apply(p))
        self.assertEqual(E.apply(p), tuple(p))

        # Check S_n(z) cyclicity
        S3z = Quat.from_axangle([0, 0, 1], sympy.pi + 2 * sympy.pi / 3, improper=True)
        S53z = Quat.from_axangle([0, 0, 1], sympy.pi + 2 * sympy.pi / 3 * 5, improper=True)

        t = S3z * S3z
        self.assertEqual(sympy.simplify(t.apply(p)), sympy.simplify(C23z.apply(p)))

        t = t * S3z
        self.assertEqual(sympy.simplify(t.apply(p)), sympy.simplify(Oz.apply(p)))

        t = t * S3z
        self.assertEqual(sympy.simplify(t.apply(p)), sympy.simplify(C3z.apply(p)))

        t = t * S3z
        self.assertEqual(sympy.simplify(t.apply(p)), sympy.simplify(S53z.apply(p)))

        t = t * S3z
        self.assertEqual(sympy.simplify(t.apply(p)), sympy.simplify(E.apply(p)))

        # check composition of C_n(z) * O(z) → S_n(z)
        n = 7
        Snz = Quat.from_axangle([0, 0, 1], sympy.pi + 2 * sympy.pi / n, improper=True)
        Cnz = Quat.from_axangle([0, 0, 1], 2 * sympy.pi / n)
        t = Cnz * Oz
        self.assertEqual(sympy.simplify(t.apply(p)), sympy.simplify(Snz.apply(p)))

        # check C_2(z) * C_2(x) → C2z
        C2z = Quat.from_axangle([0, 0, 1], angle=sympy.pi)
        C2y = Quat.from_axangle([0, 1, 0], angle=sympy.pi)
        C2x = Quat.from_axangle([1, 0, 0], angle=sympy.pi)
        t = C2y * C2x
        self.assertEqual(t.q, C2z.q)
        self.assertEqual(t.apply(p), C2z.apply(p))

        # check O(x) * O(y) → C2z
        Ox = Quat.from_axangle([1, 0, 0], angle=sympy.pi + 2 * sympy.pi, improper=True)
        Oy = Quat.from_axangle([0, 1, 0], angle=sympy.pi + 2 * sympy.pi, improper=True)
        t = Oy * Ox
        self.assertEqual(t.q, C2z.q)
        self.assertEqual(t.apply(p), C2z.apply(p))

        # Check O(y) * O(xy') → C4z
        Od = Quat.from_axangle([1, -1, 0], angle=sympy.pi + 2 * sympy.pi, improper=True)
        t = Oy * Od
        C4z = Quat.from_axangle([0, 0, 1], angle=2 * sympy.pi / 4)
        self.assertEqual(t.apply(p), C4z.apply(p))
