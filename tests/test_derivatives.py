import unittest
import numpy

from qcip_tools import derivatives, derivatives_e
from tests import float_almost_equals, array_almost_equals


class DerivativesTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_derivatives(self):
        """Test the behavior of the Derivate object"""

        # create energy
        e = derivatives.Derivative()
        self.assertEqual(e.diff_representation, '')
        self.assertEqual(e.representation(), '')
        self.assertEqual(e.dimension(), 1)
        self.assertEqual(e.shape(), [1])
        self.assertIsNone(e.basis)
        self.assertEqual(e.order(), 0)

        # test smart iterator (yes, also with energy which is basically a scalar)
        num_smart_iterator_call = 0

        r = numpy.zeros(e.shape()).flatten()
        for i in e.smart_iterator():
            num_smart_iterator_call += 1
            self.assertTrue(i < e.dimension(), i)
            for j in e.inverse_smart_iterator(i):
                r[j] += 1

        self.assertTrue(numpy.all(r == 1))
        self.assertEqual(num_smart_iterator_call, 1)

        # create derivative
        d0 = derivatives.Derivative(from_representation='FF')
        self.assertEqual(d0.diff_representation, 'FF')
        self.assertEqual(d0.representation(), 'FF')
        self.assertEqual(d0.dimension(), 3 * 3)
        self.assertEqual(d0.shape(), [3, 3])
        self.assertIsNone(d0.basis)
        self.assertEqual(d0.order(), 2)

        # test smart iterator
        num_smart_iterator_call = 0

        r = numpy.zeros(d0.shape()).flatten()
        for i in d0.smart_iterator():
            num_smart_iterator_call += 1
            self.assertTrue(i < d0.dimension(), i)
            for j in d0.inverse_smart_iterator(i):
                r[j] += 1

        self.assertTrue(numpy.all(r == 1))
        self.assertEqual(num_smart_iterator_call, 6)  # Note: 6 = 3 * (3+1) / 2

        # differentiate
        d1 = derivatives.Derivative(
            from_representation='G', basis=d0, spacial_dof=6)

        self.assertEqual(d1.diff_representation, 'G')
        self.assertEqual(d1.basis.diff_representation, 'FF')
        self.assertEqual(d1.representation(), 'GFF')

        self.assertEqual(d1.dimension(), 6 * 3 * 3)
        self.assertEqual(d1.shape(), [6, 3, 3])
        self.assertEqual(d1.order(), 3)

        # test smart iterator
        num_smart_iterator_call = 0

        r = numpy.zeros(d1.shape()).flatten()
        for i in d1.smart_iterator():
            num_smart_iterator_call += 1
            self.assertTrue(i < d1.dimension(), i)
            for j in d1.inverse_smart_iterator(i):
                r[j] += 1

        self.assertTrue(numpy.all(r == 1))
        self.assertEqual(num_smart_iterator_call, 6 * 6)

        # differentiate again:
        d2 = d1.differentiate('G')
        self.assertEqual(d2.representation(), 'GGFF')
        self.assertEqual(d2.basis.representation(), d1.representation())
        self.assertEqual(d2.order(), 4)

        # test smart iterator
        num_smart_iterator_call = 0

        r = numpy.zeros(d2.shape()).flatten()
        for i in d2.smart_iterator():
            num_smart_iterator_call += 1
            self.assertTrue(i < d2.dimension(), i)
            for j in d2.inverse_smart_iterator(i):
                r[j] += 1

        self.assertTrue(numpy.all(r == 1))
        self.assertEqual(num_smart_iterator_call, 6 * 21)  # Note: 21 = 6 * (6+1) / 2

        # make the code cry:
        with self.assertRaises(derivatives.RepresentationError):
            derivatives.Derivative(from_representation='NE', spacial_dof=6)  # "E"
        with self.assertRaises(Exception):
            derivatives.Derivative(from_representation='N')  # missing dof
        with self.assertRaises(Exception):
            derivatives.Derivative(from_representation='G')  # missing dof
        with self.assertRaises(ValueError):
            d0.differentiate('')  # nothing
        with self.assertRaises(derivatives.RepresentationError):
            d0.differentiate('E')  # "E"
        with self.assertRaises(Exception):
            d0.differentiate('N')  # missing dof
        with self.assertRaises(Exception):
            d0.differentiate('G')  # missing dof

        # test comparison
        self.assertTrue(d1 == 'GFF')
        self.assertTrue(d1 == 'FGF')  # no matter the order
        self.assertTrue(d1 != 'FGD')  # ... But the type matter
        self.assertFalse(d1 == d2)
        self.assertTrue(d1 != d2)

    def test_tensor(self):
        """Test the behavior of the Tensor object"""

        beta_tensor = numpy.array([a * (-1) ** a for a in range(27)]).reshape((3, 3, 3))
        # Note: a tensor is suppose to be symmetric, which is not the case here
        t = derivatives.Tensor('FDD', components=beta_tensor, frequency='static')

        self.assertEqual(t.representation.representation(), 'FDD')
        self.assertEqual(t.representation.dimension(), 27)
        self.assertEqual(t.frequency, 'static')
        self.assertIsNone(t.spacial_dof)
        self.assertTrue(numpy.array_equal(beta_tensor, t.components))

        # make the code cry:
        with self.assertRaises(ValueError):
            derivatives.Tensor('FDD')  # no frequency
        with self.assertRaises(Exception):
            derivatives.Tensor('N')  # no dof

    def test_electrical_derivatives_tensors(self):
        """Test the objects in derivatives_e.py

        Note: no gamma for the moment
        """

        # static water, HF/d-aug-cc-pVDZ (Gaussian)
        dipole = numpy.array([.0, .0, -0.767392])

        d = derivatives_e.ElectricDipole(dipole=dipole)
        self.assertTrue(float_almost_equals(d.norm(), 0.767392))

        alpha = numpy.array(
            [[0.777020e1, .0, .0],
             [.0, 0.895381e1, .0],
             [.0, .0, 0.832281e1]]
        )

        a = derivatives_e.PolarisabilityTensor(tensor=alpha)

        self.assertTrue(float_almost_equals(a.isotropic_value(), 0.834894e1))
        self.assertTrue(float_almost_equals(a.anisotropic_value(), 0.102578e1))

        beta = numpy.array(
            [[[.0, .0, -0.489173],
              [.0, .0, .0],
              [-0.489173, .0, .0]],
             [[.0, .0, .0],
              [.0, .0, 9.080568],
              [.0, 9.080568, .0]],
             [[-0.489173, .0, .0],
              [.0, 9.080568, .0],
              [.0, .0, 4.276656]]]
        )

        b = derivatives_e.FirstHyperpolarisabilityTensor(tensor=beta)

        self.assertTrue(float_almost_equals(b.beta_squared_zzz(), 29.4147))
        self.assertTrue(float_almost_equals(b.beta_squared_zxx(), 8.5707))
        self.assertTrue(float_almost_equals(b.beta_hrs(), 6.1632))
        self.assertTrue(float_almost_equals(b.depolarization_ratio(), 3.4320))
        self.assertTrue(float_almost_equals(b.octupolar_contribution_squared(), 167.0257))
        self.assertTrue(float_almost_equals(b.dipolar_contribution_squared(), 99.3520))
        self.assertTrue(float_almost_equals(b.nonlinear_anisotropy(), 1.2966))

        self.assertTrue(float_almost_equals(b.beta_parallel(dipole), -7.7208))
        self.assertTrue(float_almost_equals(b.beta_perpendicular(dipole), -2.5736))
        self.assertTrue(float_almost_equals(b.beta_kerr(dipole), -7.7208))

        self.assertTrue(array_almost_equals(b.beta_vector(), [.0, .0, 12.8681]))

        # static NH3, HF/d-aug-cc-pVDZ (Gaussian)
        dipole = numpy.array([.0, .0, 0.625899])

        d = derivatives_e.ElectricDipole(dipole=dipole)
        self.assertTrue(float_almost_equals(d.norm(), 0.625899))

        alpha = numpy.array(
            [[0.125681e2, .0, -0.485486e-4],
             [.0, 0.125681e2, .0],
             [-0.485486e-4, .0, 0.132024e2]]
        )

        a = derivatives_e.PolarisabilityTensor(tensor=alpha)

        self.assertTrue(float_almost_equals(a.isotropic_value(), 0.127795e2))
        self.assertTrue(float_almost_equals(a.anisotropic_value(), 0.634350))

        beta = numpy.array(
            [[[9.258607, -0.012368, -6.097955],
              [-0.012368, -9.257993, .0],
              [-6.097955, .0, -0.000073]],
             [[-0.012368, -9.257993, .0],
              [-9.257993, 0.012368, -6.097633],
              [.0, -6.097633, .0]],
             [[-6.097955, .0, -0.000073],
              [.0, -6.097633, .0],
              [-0.000073, .0, -6.483421]]]
        )

        b = derivatives_e.FirstHyperpolarisabilityTensor(tensor=beta)

        self.assertTrue(float_almost_equals(b.beta_squared_zzz(), 64.6483))
        self.assertTrue(float_almost_equals(b.beta_squared_zxx(), 19.8385))
        self.assertTrue(float_almost_equals(b.beta_hrs(), 9.1917))
        self.assertTrue(float_almost_equals(b.depolarization_ratio(), 3.2587))
        self.assertTrue(float_almost_equals(b.octupolar_contribution_squared(), 398.6438))
        self.assertTrue(float_almost_equals(b.dipolar_contribution_squared(), 209.3432))
        self.assertTrue(float_almost_equals(b.nonlinear_anisotropy(), 1.3799))

        self.assertTrue(float_almost_equals(b.beta_parallel(dipole), -11.2074))
        self.assertTrue(float_almost_equals(b.beta_perpendicular(dipole), -3.7358))
        self.assertTrue(float_almost_equals(b.beta_kerr(dipole), -11.2074))

        self.assertTrue(array_almost_equals(b.beta_vector(), [.0, .0, -18.6790]))

        # static CH4, HF/d-aug-cc-pVDZ (Gaussian)
        alpha = numpy.array(
            [[0.159960e2, .0, .0],
             [.0, 0.159960e2, .0],
             [.0, .0, 0.159960e2]]
        )

        a = derivatives_e.PolarisabilityTensor(tensor=alpha)

        self.assertTrue(float_almost_equals(a.isotropic_value(), 0.159960e2))
        self.assertTrue(float_almost_equals(a.anisotropic_value(), .0))

        beta = numpy.array(
            [[[.0, .0, .0],
              [.0, .0, -11.757505],
              [.0, -11.757505, .0]],
             [[.0, .0, -11.757505],
              [.0, .0, .0],
              [-11.757505, .0, .0]],
             [[.0, -11.757505, .0],
              [-11.757505, .0, .0],
              [.0, .0, .0]]]
        )

        b = derivatives_e.FirstHyperpolarisabilityTensor(tensor=beta)

        self.assertTrue(float_almost_equals(b.beta_squared_zzz(), 47.3962))
        self.assertTrue(float_almost_equals(b.beta_squared_zxx(), 31.5975))
        self.assertTrue(float_almost_equals(b.beta_hrs(), 8.8878))
        self.assertTrue(float_almost_equals(b.depolarization_ratio(), 1.5))
        self.assertTrue(float_almost_equals(b.octupolar_contribution_squared(), 829.4336))
        self.assertTrue(float_almost_equals(b.dipolar_contribution_squared(), .0))

        # since CH4 has no dipole moment, the rest of the calculations failed ;)
