from tests import QcipToolsTestCase
from qcip_tools import basis_set, basis_set_esml


class BasisSetTestCase(QcipToolsTestCase):

    def setUp(self):
        # STO-3G for hydrogen:
        self.bf_H_S = basis_set.Function()
        self.bf_H_S.add_primitive(basis_set.Primitive(3.42525091, 0.15432897))
        self.bf_H_S.add_primitive(basis_set.Primitive(0.62391373, 0.53532814))
        self.bf_H_S.add_primitive(basis_set.Primitive(0.16885540, 0.44463454))

        # STO-3G for oxygen (with SP):
        self.bf_O_S = basis_set.Function()
        self.bf_O_S.add_primitive(basis_set.Primitive(130.7093200, 0.15432897))
        self.bf_O_S.add_primitive(basis_set.Primitive(23.8088610, 0.53532814))
        self.bf_O_S.add_primitive(basis_set.Primitive(6.4436083, 0.44463454))

        self.bf_O_SP = basis_set.Function()
        self.bf_O_SP.add_primitive(basis_set.Primitive(5.0331513, -0.09996723, 0.15591627))
        self.bf_O_SP.add_primitive(basis_set.Primitive(1.1695961, 0.39951283, 0.60768372))
        self.bf_O_SP.add_primitive(basis_set.Primitive(0.3803890, 0.70011547, 0.39195739))

        # Oyxgen d-type polarization function [from 6-31G(d)]
        self.bf_O_D = basis_set.Function()
        self.bf_O_D.add_primitive(basis_set.Primitive(1.2920000, 1.))

    def test_basis_function(self):
        """Test basis functions"""

        # test ordering for primitive (largest exponent first)
        f1 = basis_set.Function()
        f1.add_primitive(basis_set.Primitive(1.))
        f1.add_primitive(basis_set.Primitive(5.))
        f1.add_primitive(basis_set.Primitive(3.))

        self.assertEqual(f1.primitives[0].exponent, 5.)
        self.assertEqual(f1.primitives[1].exponent, 3.)

        # test is_diffuse()
        fe = basis_set.Function()
        self.assertFalse(fe.is_diffuse())

        fe.add_primitive(basis_set.Primitive(.1, .1))
        self.assertFalse(fe.is_diffuse())  # contraction coefficient != 1

        fe.primitives[0] = basis_set.Primitive(5.5, 1.)
        self.assertFalse(fe.is_diffuse())  # exponent is too large

        fe.primitives[0] = basis_set.Primitive(.01, 1., 4.)
        self.assertFalse(fe.is_diffuse())  # p coefficient is != 1

        fe.primitives[0] = basis_set.Primitive(.01, 1.)
        self.assertFalse(fe.is_diffuse(exponent_threshold=.001))  # threshold

        fe.primitives[0] = basis_set.Primitive(.01, 1)
        self.assertTrue(fe.is_diffuse())  # OK!

    def test_atomic_basis_set(self):
        """Test the behavior of atomic basis set"""

        b_H = basis_set.AtomicBasisSet('H')
        self.assertEqual(b_H.atom.atomic_number, 1)
        b_H.add_basis_function('s', self.bf_H_S)

        self.assertTrue('S' in b_H)
        self.assertFalse('P' in b_H)
        self.assertEqual(len(b_H.basis_functions_per_shell['S']), 1)

        self.assertEqual(str(b_H), 'H [3s|1s]')

        b_O = basis_set.AtomicBasisSet('O')
        self.assertEqual(b_O.atom.atomic_number, 8)

        b_O.add_basis_function('s', self.bf_O_S)
        b_O.add_basis_function('sp', self.bf_O_SP)

        self.assertTrue('S' in b_O)
        self.assertTrue('SP' in b_O)
        self.assertFalse('D' in b_O)
        self.assertEqual(len(b_O.basis_functions_per_shell['S']), 1)
        self.assertEqual(len(b_O.basis_functions_per_shell['SP']), 1)

        self.assertEqual(str(b_O), 'O [6s3p|2s1p]')

        b_O.add_basis_function('d', self.bf_O_D)
        self.assertTrue('D' in b_O)
        self.assertEqual(len(b_O.basis_functions_per_shell['D']), 1)
        self.assertEqual(str(b_O), 'O [6s3p1d|2s1p1d]')

        # fire exceptions
        with self.assertRaises(Exception):
            b_O.add_basis_function('D', basis_set.Function())  # empty basis function

        with self.assertRaises(Exception):
            f = basis_set.Function()
            f.add_primitive(basis_set.Primitive(1., 1.))
            b_O.add_basis_function('P', f)  # P while SP

        with self.assertRaises(Exception):
            f = basis_set.Function()
            f.add_primitive(basis_set.Primitive(1., 1.))
            b_O.add_basis_function('X', f)  # non-existing shell

        with self.assertRaises(Exception):
            'x' in b_O  # non-existing shell

        # test ordering for basis function (ordered by size, then largest exponent)
        f1 = basis_set.Function()
        f1.add_primitive(basis_set.Primitive(1.))
        f1.add_primitive(basis_set.Primitive(5.))
        f1.add_primitive(basis_set.Primitive(3.))

        f2 = basis_set.Function()
        f2.add_primitive(basis_set.Primitive(10.))
        f2.add_primitive(basis_set.Primitive(20.))

        f3 = basis_set.Function()
        f3.add_primitive(basis_set.Primitive(40.))
        f3.add_primitive(basis_set.Primitive(60.))

        bf_X = basis_set.AtomicBasisSet(atomic_number=2)  # whatever
        bf_X.add_basis_function('s', f2)
        bf_X.add_basis_function('s', f3)
        bf_X.add_basis_function('s', f1)

        self.assertEqual(bf_X.basis_functions_per_shell['S'][0], f1)
        self.assertEqual(bf_X.basis_functions_per_shell['S'][1], f3)

    def test_basis_set(self):
        """Test the behavior of the basis set"""

        b_H = basis_set.AtomicBasisSet('H')
        b_H.add_basis_function('s', self.bf_H_S)

        b_O = basis_set.AtomicBasisSet('O')
        self.assertEqual(b_O.atom.atomic_number, 8)

        b_O.add_basis_function('s', self.bf_O_S)
        b_O.add_basis_function('sp', self.bf_O_SP)

        sto3g = basis_set.BasisSet('STO-3G')

        sto3g.add_atomic_basis_set(b_O)
        sto3g.add_atomic_basis_set(b_H)

        self.assertTrue('O' in sto3g)
        self.assertTrue(8 in sto3g)  # can also use atomic numbers
        self.assertTrue('H' in sto3g)
        self.assertFalse('Ne' in sto3g)

        # exceptions
        with self.assertRaises(Exception):
            'x' in sto3g  # non-existing atom

        with self.assertRaises(Exception):
            44444 in sto3g  # non-existing atom

    def test_ESML(self):
        """Test scrapping data from ESML basis set exchange"""

        with self.assertRaises(basis_set_esml.ESMLBasisSetError):
            basis_set_esml.get_atomic_basis_set('xxx', ['C', 'H'])  # basis set does not exists
        with self.assertRaises(basis_set_esml.ESMLBasisSetError):
            basis_set_esml.get_atomic_basis_set('STO-3G', ['C', 'H', 'D'])  # unknown atom
        with self.assertRaises(basis_set_esml.ESMLBasisSetError):
            basis_set_esml.get_atomic_basis_set('STO-3G', [1, 162])  # unknown atom
        with self.assertRaises(basis_set_esml.ESMLBasisSetError):
            basis_set_esml.get_atomic_basis_set('STO-3G', [])  # empty list of atom
        with self.assertRaises(basis_set_esml.ESMLBasisSetError):
            basis_set_esml.get_atomic_basis_set('STO-3G', ['C', 'H'], basis_set_format='x')  # unknown format

        b = basis_set_esml.get_atomic_basis_set('STO-3G', ['C', 'H'])
        self.assertIn('C', b)
        self.assertIn('H', b)
