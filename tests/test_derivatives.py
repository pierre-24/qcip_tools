import random
import math

import numpy
from qcip_tools import derivatives, derivatives_e, math as qcip_math, derivatives_g, atom as qcip_atom, \
    molecule as qcip_molecule, numerical_differentiation
from tests import QcipToolsTestCase, factories


class DerivativesTestCase(QcipToolsTestCase):

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
        for i in e.smart_iterator(as_flatten=True):
            num_smart_iterator_call += 1
            self.assertTrue(i < e.dimension(), i)
            for j in e.inverse_smart_iterator(i, as_flatten=True):
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
        for i in d0.smart_iterator(as_flatten=True):
            num_smart_iterator_call += 1
            self.assertTrue(i < d0.dimension(), i)
            for j in d0.inverse_smart_iterator(i, as_flatten=True):
                r[j] += 1

        self.assertTrue(numpy.all(r == 1))
        self.assertEqual(num_smart_iterator_call, 6)  # Note: 6 = 3 * (3+1) / 2

        # differentiate
        d1 = derivatives.Derivative(
            from_representation='G', basis=d0, spacial_dof=9)

        self.assertEqual(d1.diff_representation, 'G')
        self.assertEqual(d1.basis.diff_representation, 'FF')
        self.assertEqual(d1.representation(), 'GFF')

        self.assertEqual(d1.dimension(), 9 * 3 * 3)
        self.assertEqual(d1.shape(), [9, 3, 3])
        self.assertEqual(d1.order(), 3)

        # test smart iterator
        num_smart_iterator_call = 0

        r = numpy.zeros(d1.shape()).flatten()
        for i in d1.smart_iterator(as_flatten=True):
            num_smart_iterator_call += 1
            self.assertTrue(i < d1.dimension(), i)
            for j in d1.inverse_smart_iterator(i, as_flatten=True):
                r[j] += 1

        self.assertTrue(numpy.all(r == 1))
        self.assertEqual(num_smart_iterator_call, 9 * 6)

        # differentiate again:
        d2 = d1.differentiate('G')
        self.assertEqual(d2.representation(), 'GGFF')
        self.assertEqual(d2.basis.representation(), d1.representation())
        self.assertEqual(d2.order(), 4)

        # test smart iterator
        num_smart_iterator_call = 0

        r = numpy.zeros(d2.shape()).flatten()
        for i in d2.smart_iterator(as_flatten=True):
            num_smart_iterator_call += 1
            self.assertTrue(i < d2.dimension(), i)
            for j in d2.inverse_smart_iterator(i, as_flatten=True):
                r[j] += 1

        self.assertTrue(numpy.all(r == 1))
        self.assertEqual(num_smart_iterator_call, 6 * 45)  # Note: 45 = 9 * (9+1) / 2

        # tricky one:
        d4 = derivatives.Derivative(from_representation='FDF')
        self.assertEqual(d4.diff_representation, 'FDF')
        self.assertEqual(d4.representation(), 'FDF')
        self.assertEqual(d4.dimension(), 3 * 3 * 3)
        self.assertEqual(d4.shape(), [3, 3, 3])
        self.assertIsNone(d4.basis)
        self.assertEqual(d4.order(), 3)

        # test smart iterator
        num_smart_iterator_call = 0

        r = numpy.zeros(d4.shape()).flatten()
        for i in d4.smart_iterator(as_flatten=True):
            num_smart_iterator_call += 1
            self.assertTrue(i < d4.dimension(), i)
            for j in d4.inverse_smart_iterator(i, as_flatten=True):
                r[j] += 1

        self.assertEqual(num_smart_iterator_call, 18)  # 3 * 6
        self.assertTrue(numpy.all(r == 1))

        # another tricky one:
        d5 = derivatives.Derivative(from_representation='FDDF')
        self.assertEqual(d5.diff_representation, 'FDDF')
        self.assertEqual(d5.representation(), 'FDDF')
        self.assertEqual(d5.dimension(), 3 * 3 * 3 * 3)
        self.assertEqual(d5.shape(), [3, 3, 3, 3])
        self.assertIsNone(d5.basis)
        self.assertEqual(d5.order(), 4)

        # test smart iterator
        num_smart_iterator_call = 0

        r = numpy.zeros(d5.shape()).flatten()
        for i in d5.smart_iterator(as_flatten=True):
            num_smart_iterator_call += 1
            self.assertTrue(i < d5.dimension(), i)
            for j in d5.inverse_smart_iterator(i, as_flatten=True):
                r[j] += 1

        self.assertEqual(num_smart_iterator_call, 36)  # 6 * 6
        self.assertTrue(numpy.all(r == 1))

        # once again, but with full components (not flatten indices)
        num_smart_iterator_call = 0

        r = numpy.zeros(d5.shape())
        for i in d5.smart_iterator(as_flatten=False):
            num_smart_iterator_call += 1
            self.assertEqual(len(i), d5.order())
            self.assertTrue(all(x < y for x, y in zip(i, d5.shape())))
            for j in d5.inverse_smart_iterator(i, as_flatten=False):
                r[j] += 1

        self.assertEqual(num_smart_iterator_call, 36)  # 6 * 6
        self.assertTrue(numpy.all(r.flatten() == 1))

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

        gamma_tensor = numpy.array([a * (-1) ** a for a in range(81)]).reshape((3, 3, 3, 3))
        # Note: a tensor is suppose to be symmetric, which is not the case here
        t = derivatives.Tensor('FdDD', components=gamma_tensor, frequency='static')

        self.assertEqual(t.representation.representation(), 'FdDD')
        self.assertEqual(t.representation.dimension(), 81)
        self.assertEqual(t.frequency, 'static')
        self.assertIsNone(t.spacial_dof)
        self.assertTrue(numpy.array_equal(gamma_tensor, t.components))

        # make the code cry:
        with self.assertRaises(ValueError):
            derivatives.Tensor('FDD')  # no frequency
        with self.assertRaises(Exception):
            derivatives.Tensor('N')  # no dof

    def test_electrical_derivatives_tensors(self):
        """Test the objects in derivatives_e.py

        Also test that the properties that are invariant under rotation remain invariant !
        """

        angles_set = [
            (0, 0, 0),  # first test without any rotation
            (180, 60, -60),
            (15, -15, 25),
            (450, -75, 75),
            # and the last one, to be certain
            (random.randrange(-360, 360), random.randrange(-360, 360), random.randrange(-360, 360))
        ]

        # static water, HF/d-aug-cc-pVDZ (Gaussian)
        dipole = numpy.array([.0, .0, -0.767392])

        d = derivatives_e.ElectricDipole(dipole=dipole)
        self.assertAlmostEqual(d.norm(), 0.767392)

        alpha = numpy.array(
            [[0.777020e1, .0, .0],
             [.0, 0.895381e1, .0],
             [.0, .0, 0.832281e1]]
        )

        for angles in angles_set:  # invariants remain invariants under rotation
            new_alpha = qcip_math.tensor_rotate(alpha, *angles)
            na = derivatives_e.PolarisabilityTensor(tensor=new_alpha)

            self.assertAlmostEqual(na.isotropic_value(), 0.834894e1)
            self.assertAlmostEqual(na.anisotropic_value(), 0.102578e1, places=3)

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

        self.assertAlmostEqual(b.beta_parallel(dipole), -7.7208, places=3)
        self.assertAlmostEqual(b.beta_perpendicular(dipole), -2.5736, places=3)
        self.assertAlmostEqual(b.beta_kerr(dipole), -7.7208, places=3)
        self.assertArraysAlmostEqual(b.beta_vector(), [.0, .0, 12.8681])

        # NOTE: above properties are also invariant to rotation, but in a less funny way.

        for angles in angles_set:
            new_beta = qcip_math.tensor_rotate(beta, *angles)
            nb = derivatives_e.FirstHyperpolarisabilityTensor(tensor=new_beta)

            self.assertAlmostEqual(nb.beta_squared_zzz(), 29.4147, places=3)
            self.assertAlmostEqual(nb.beta_squared_zxx(), 8.5707, places=3)
            self.assertAlmostEqual(nb.beta_hrs(), 6.1632, places=3)
            self.assertAlmostEqual(nb.depolarization_ratio(), 3.4320, places=3)
            self.assertAlmostEqual(nb.octupolar_contribution_squared(), 167.0257, places=3)
            self.assertAlmostEqual(nb.dipolar_contribution_squared(), 99.3520, places=3)
            self.assertAlmostEqual(nb.quadrupolar_contribution_squared(), .0, places=3)
            self.assertAlmostEqual(nb.nonlinear_anisotropy(), 1.2966, places=3)

        # static NH3, HF/d-aug-cc-pVDZ (Gaussian)
        dipole = derivatives_e.ElectricDipole(dipole=[.0, .0, 0.625899])

        d = derivatives_e.ElectricDipole(dipole=dipole.components)
        self.assertAlmostEqual(d.norm(), 0.625899)

        alpha = numpy.array(
            [[0.125681e2, .0, -0.485486e-4],
             [.0, 0.125681e2, .0],
             [-0.485486e-4, .0, 0.132024e2]]
        )

        for angles in angles_set:
            new_alpha = qcip_math.tensor_rotate(alpha, *angles)
            na = derivatives_e.PolarisabilityTensor(tensor=new_alpha)

            self.assertAlmostEqual(na.isotropic_value(), 0.127795e2, places=3)
            self.assertAlmostEqual(na.anisotropic_value(), 0.634350, places=3)

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

        self.assertAlmostEqual(b.beta_parallel(dipole), -11.2074, places=3)
        self.assertAlmostEqual(b.beta_perpendicular(dipole), -3.7358, places=3)
        self.assertAlmostEqual(b.beta_kerr(dipole), -11.2074, places=3)
        self.assertArraysAlmostEqual(b.beta_vector(), [.0, .0, -18.6790])

        for angles in angles_set:
            new_beta = qcip_math.tensor_rotate(beta, *angles)
            nb = derivatives_e.FirstHyperpolarisabilityTensor(tensor=new_beta)

            self.assertAlmostEqual(nb.beta_squared_zzz(), 64.6483, places=3)
            self.assertAlmostEqual(nb.beta_squared_zxx(), 19.8385, places=3)
            self.assertAlmostEqual(nb.beta_hrs(), 9.1917, places=3)
            self.assertAlmostEqual(nb.depolarization_ratio(), 3.2587, places=3)
            self.assertAlmostEqual(nb.octupolar_contribution_squared(), 398.6438, places=3)
            self.assertAlmostEqual(nb.dipolar_contribution_squared(), 209.3432, places=3)
            self.assertAlmostEqual(nb.quadrupolar_contribution_squared(), .0, places=3)
            self.assertAlmostEqual(nb.nonlinear_anisotropy(), 1.3799, places=3)

        # static CH4, HF/d-aug-cc-pVDZ (Gaussian)
        alpha = numpy.array(
            [[0.159960e2, .0, .0],
             [.0, 0.159960e2, .0],
             [.0, .0, 0.159960e2]]
        )

        for angles in angles_set:
            new_alpha = qcip_math.tensor_rotate(alpha, *angles)
            na = derivatives_e.PolarisabilityTensor(tensor=new_alpha)

            self.assertAlmostEqual(na.isotropic_value(), 0.159960e2, places=3)
            self.assertAlmostEqual(na.anisotropic_value(), .0, places=3)

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

        for angles in angles_set:
            new_beta = qcip_math.tensor_rotate(beta, *angles)
            nb = derivatives_e.FirstHyperpolarisabilityTensor(tensor=new_beta)

            self.assertAlmostEqual(nb.beta_squared_zzz(), 47.3962, places=3)
            self.assertAlmostEqual(nb.beta_squared_zxx(), 31.5975, places=3)
            self.assertAlmostEqual(nb.beta_hrs(), 8.8878, places=3)
            self.assertAlmostEqual(nb.depolarization_ratio(), 1.5, places=3)
            self.assertAlmostEqual(nb.octupolar_contribution_squared(), 829.4336, places=3)
            self.assertAlmostEqual(nb.dipolar_contribution_squared(), .0, places=3)
            self.assertAlmostEqual(nb.quadrupolar_contribution_squared(), .0, places=3)

        # ... since CH4 has no dipole moment, the rest of the properties failed ;)

        # dynamic (911.3nm) water BLYP/d-aug-cc-pVTZ (Dalton)
        beta = numpy.array(
            [[[0.000000e+00, 0.000000e+00, -1.164262e+01],
              [0.000000e+00, 0.000000e+00, 0.000000e+00],
              [-1.164262e+01, 0.000000e+00, 0.000000e+00]],
             [[0.000000e+00, 0.000000e+00, 0.000000e+00],
              [0.000000e+00, 0.000000e+00, -1.451377e+01],
              [0.000000e+00, -1.451377e+01, 0.000000e+00]],
             [[-9.491719e+00, 0.000000e+00, 0.000000e+00],
              [0.000000e+00, -1.475714e+01, 0.000000e+00],
              [0.000000e+00, 0.000000e+00, -2.222017e+01]]]
        )

        for angles in angles_set:
            new_beta = qcip_math.tensor_rotate(beta, *angles)
            nb = derivatives_e.FirstHyperpolarisabilityTensor(tensor=new_beta)

            # Values obtained directly from the contribution matrices
            self.assertAlmostEqual(nb.octupolar_contribution_squared(), 123.3727, places=3)
            # self.assertAlmostEqual(nb.dipolar_contribution_squared(), 1367.5055 + 1.2128, places=3)
            self.assertAlmostEqual(nb.quadrupolar_contribution_squared(), 1.9108, places=3)

            self.assertAlmostEqual(nb.beta_squared_zzz(), 280.551, places=3)
            self.assertAlmostEqual(nb.beta_squared_zxx(), 31.304, places=3)
            self.assertAlmostEqual(nb.beta_hrs(), 17.659, places=3)
            self.assertAlmostEqual(nb.depolarization_ratio(), 8.962, places=3)

        # static CH2Cl2, CCS/d-aug-cc-pVDZ (dalton)
        gamma = numpy.array(
            [[[[-4.959904e+03, 4.730379e-03, 0.000000e+00],
               [4.730379e-03, -1.790611e+03, -1.524347e-03],
               [0.000000e+00, -1.524347e-03, -2.052200e+03]],
              [[4.730379e-03, -1.790611e+03, -1.524347e-03],
               [-1.790611e+03, 2.387447e-03, 0.000000e+00],
               [-1.524347e-03, 0.000000e+00, 3.412659e-03]],
              [[0.000000e+00, -1.524347e-03, -2.052200e+03],
               [-1.524347e-03, 0.000000e+00, 3.412659e-03],
               [-2.052200e+03, 3.412659e-03, 0.000000e+00]]],
             [[[4.730379e-03, -1.790611e+03, -1.524347e-03],
               [-1.790611e+03, 2.387447e-03, 0.000000e+00],
               [-1.524347e-03, 0.000000e+00, 3.412659e-03]],
              [[-1.790611e+03, 2.387447e-03, 0.000000e+00],
               [2.387447e-03, -5.193209e+03, -6.751678e-04],
               [0.000000e+00, -6.751678e-04, -2.207921e+03]],
              [[-1.524347e-03, 0.000000e+00, 3.412659e-03],
               [0.000000e+00, -6.751678e-04, -2.207921e+03],
               [3.412659e-03, -2.207921e+03, -2.799551e-03]]],
             [[[0.000000e+00, -1.524347e-03, -2.052200e+03],
               [-1.524347e-03, 0.000000e+00, 3.412659e-03],
               [-2.052200e+03, 3.412659e-03, 0.000000e+00]],
              [[-1.524347e-03, 0.000000e+00, 3.412659e-03],
               [0.000000e+00, -6.751678e-04, -2.207921e+03],
               [3.412659e-03, -2.207921e+03, -2.799551e-03]],
              [[-2.052200e+03, 3.412659e-03, 0.000000e+00],
               [3.412659e-03, -2.207921e+03, -2.799551e-03],
               [0.000000e+00, -2.799551e-03, -9.412690e+03]]]]
        )

        orig_g = derivatives_e.SecondHyperpolarizabilityTensor(tensor=gamma)

        for angles in angles_set:
            new_gamma = qcip_math.tensor_rotate(gamma, *angles)
            ng = derivatives_e.SecondHyperpolarizabilityTensor(tensor=new_gamma)

            self.assertAlmostEqual(ng.gamma_parallel(), orig_g.gamma_parallel(), places=3)
            self.assertAlmostEqual(ng.gamma_perpendicular(), orig_g.gamma_perpendicular(), places=3)
            self.assertAlmostEqual(ng.gamma_kerr(), orig_g.gamma_kerr(), places=3)

            self.assertAlmostEqual(ng.gamma_squared_zzzz(), orig_g.gamma_squared_zzzz(), places=3)
            self.assertAlmostEqual(ng.gamma_squared_zxxx(), orig_g.gamma_squared_zxxx(), places=3)
            self.assertAlmostEqual(ng.gamma_ths(), orig_g.gamma_ths())
            self.assertAlmostEqual(ng.depolarization_ratio(), orig_g.depolarization_ratio(), places=3)

            GJ4 = ng.hexacadecapolar_contribution_squared(first_version=True)
            GJ2 = ng.quadrupolar_contribution_squared(first_version=True)
            GJ0 = ng.isotropic_contribution_squared()

            rho_2 = GJ2 / GJ0
            rho_4 = GJ4 / GJ0
            rho_0m = GJ0 / GJ2
            rho_4m = GJ4 / GJ2

            self.assertAlmostEqual(ng.gamma_squared_zzzz(), 8 / 315 * GJ4 + 4 / 35 * GJ2 + 1 / 5 * GJ0, places=3)
            self.assertAlmostEqual(ng.gamma_squared_zxxx(), 1 / 63 * GJ4 + 3 / 140 * GJ2, places=3)
            self.assertAlmostEqual(
                ng.depolarization_ratio(), (32 * rho_4 + 144 * rho_2 + 252) / (20 * rho_4 + 27 * rho_2), places=3)
            self.assertAlmostEqual(
                ng.depolarization_ratio(), (32 * rho_4m + 252 * rho_0m + 144) / (20 * rho_4m + 27), places=3)

        # static CH4, CCS/d-aug-cc-pVDZ (dalton)
        gamma = numpy.array([
            [[[-2.229953e+03, -3.839926e-05, -3.112595e-05],
              [-3.839926e-05, -8.607108e+02, 1.515257e-05],
              [-3.112595e-05, 1.515257e-05, -8.607108e+02]],
             [[-3.839926e-05, -8.607108e+02, 1.515257e-05],
              [-8.607108e+02, -2.006040e-05, -1.040868e-05],
              [1.515257e-05, -1.040868e-05, -1.487595e-05]],
             [[-3.112595e-05, 1.515257e-05, -8.607108e+02],
              [1.515257e-05, -1.040868e-05, -1.487595e-05],
              [-8.607108e+02, -1.487595e-05, -5.166972e-05]]],
            [[[-3.839926e-05, -8.607108e+02, 1.515257e-05],
              [-8.607108e+02, -2.006040e-05, -1.040868e-05],
              [1.515257e-05, -1.040868e-05, -1.487595e-05]],
             [[-8.607108e+02, -2.006040e-05, -1.040868e-05],
              [-2.006040e-05, -2.229954e+03, 1.370702e-05],
              [-1.040868e-05, 1.370702e-05, -8.607109e+02]],
             [[1.515257e-05, -1.040868e-05, -1.487595e-05],
              [-1.040868e-05, 1.370702e-05, -8.607109e+02],
              [-1.487595e-05, -8.607109e+02, 4.011427e-05]]],
            [[[-3.112595e-05, 1.515257e-05, -8.607108e+02],
              [1.515257e-05, -1.040868e-05, -1.487595e-05],
              [-8.607108e+02, -1.487595e-05, -5.166972e-05]],
             [[1.515257e-05, -1.040868e-05, -1.487595e-05],
              [-1.040868e-05, 1.370702e-05, -8.607109e+02],
              [-1.487595e-05, -8.607109e+02, 4.011427e-05]],
             [[-8.607108e+02, -1.487595e-05, -5.166972e-05],
              [-1.487595e-05, -8.607109e+02, 4.011427e-05],
              [-5.166972e-05, 4.011427e-05, -2.229954e+03]]]
        ])

        orig_g = derivatives_e.SecondHyperpolarizabilityTensor(tensor=gamma)

        for angles in angles_set:
            new_gamma = qcip_math.tensor_rotate(gamma, *angles)
            ng = derivatives_e.SecondHyperpolarizabilityTensor(tensor=new_gamma)

            self.assertAlmostEqual(ng.gamma_parallel(), orig_g.gamma_parallel())
            self.assertAlmostEqual(ng.gamma_perpendicular(), orig_g.gamma_perpendicular())
            self.assertAlmostEqual(ng.gamma_kerr(), orig_g.gamma_kerr())

            self.assertAlmostEqual(ng.gamma_squared_zzzz(), orig_g.gamma_squared_zzzz())
            self.assertAlmostEqual(ng.gamma_squared_zxxx(), orig_g.gamma_squared_zxxx())
            self.assertAlmostEqual(ng.gamma_ths(), orig_g.gamma_ths())
            self.assertAlmostEqual(ng.depolarization_ratio(), orig_g.depolarization_ratio())

            GJ4 = ng.hexacadecapolar_contribution_squared(first_version=True)
            GJ2 = ng.quadrupolar_contribution_squared(first_version=True)
            GJ0 = ng.isotropic_contribution_squared()

            rho_2 = GJ2 / GJ0
            rho_4 = GJ4 / GJ0
            rho_0m = GJ0 / GJ2
            rho_4m = GJ4 / GJ2

            self.assertAlmostEqual(ng.gamma_squared_zzzz(), 8 / 315 * GJ4 + 4 / 35 * GJ2 + 1 / 5 * GJ0, places=3)
            self.assertAlmostEqual(ng.gamma_squared_zxxx(), 1 / 63 * GJ4 + 3 / 140 * GJ2, places=3)
            self.assertAlmostEqual(
                ng.depolarization_ratio(), (32 * rho_4 + 144 * rho_2 + 252) / (20 * rho_4 + 27 * rho_2), places=3)
            self.assertAlmostEqual(
                ng.depolarization_ratio(), (32 * rho_4m + 252 * rho_0m + 144) / (20 * rho_4m + 27), places=3)

        # test conversion
        self.assertAlmostEqual(derivatives_e.convert_frequency_from_string('1064nm'), 0.0428, places=3)
        self.assertAlmostEqual(derivatives_e.convert_frequency_from_string('2eV'), 0.073, places=3)
        self.assertAlmostEqual(derivatives_e.convert_frequency_from_string('1500cm-1'), 0.0068, places=4)

        # order and name
        g = derivatives_e.BaseElectricalDerivativeTensor(input_fields=(0, 1))
        self.assertEqual(g.representation.representation(), 'FDF')
        self.assertEqual(g.name, 'beta(-w;w,0)')
        self.assertEqual(g.rank(), 3)

        # just check that DFWM is now possible:
        g = derivatives_e.BaseElectricalDerivativeTensor(input_fields=(1, 1, -1))
        self.assertEqual(g.representation.representation(), 'FDDd')
        self.assertEqual(g.name, 'gamma(-w;w,w,-w)')
        self.assertEqual(g.rank(), 4)

        g = derivatives_e.BaseElectricalDerivativeTensor(input_fields=(-1, 1, 1))
        self.assertEqual(g.representation.representation(), 'FDDd')  # reordering
        self.assertEqual(g.name, 'gamma(-w;w,w,-w)')

    def test_geometrical_derivatives(self):
        """Test geometrical ones.

        Note: only tests the Hessian
        """

        # H2O molecule:
        atom_list = [
            qcip_atom.Atom(symbol='O', position=[0, .0, 0.11393182]),
            qcip_atom.Atom(symbol='H', position=[0, -0.75394266, -0.45572727]),
            qcip_atom.Atom(symbol='H', position=[0, 0.75394266, -0.45572727])
        ]

        water_molecule = qcip_molecule.Molecule(atom_list=atom_list)

        # Water (geometry above) HF/Sadlej-POL (cartesian set)
        hessian = numpy.array(
            [[-0.00011, .0, .0, 0.00005, .0, .0, 0.00005, .0, .0],
             [.0, 0.80982, .0, .0, -0.40491, -0.30598, .0, -0.40491, 0.30598],
             [.0, .0, 0.53532, .0, -0.23905, -0.26766, .0, 0.23905, -0.26766],
             [0.00005, .0, .0, .0, .0, .0, -0.00005, .0, .0],
             [.0, -0.40491, -0.23905, .0, 0.43721, 0.27252, .0, -0.03230, -0.03346],
             [.0, -0.30598, -0.26766, .0, 0.27252, 0.24945, .0, 0.03346, 0.01821],
             [0.00005, .0, .0, -0.00005, .0, .0, .0, .0, .0],
             [.0, -0.40491, 0.23905, .0, -0.03230, 0.03346, .0, 0.43721, -0.27252],
             [.0, 0.30598, -0.26766, .0, -0.03346, 0.01821, .0, -0.27252, 0.24945]]
        )

        h = derivatives_g.BaseGeometricalDerivativeTensor(
            3 * len(water_molecule), 5 if water_molecule.linear() else 6, 'GG', components=hessian)

        mwh = derivatives_g.MassWeightedHessian(water_molecule, hessian)
        self.assertEqual(mwh.vibrational_dof, 3)
        self.assertFalse(mwh.linear)

        self.assertEqual(len(mwh.frequencies), 9)

        self.assertArraysAlmostEqual(
            [-38.44, -3.8, -0.03, 0.0, 18.1, 36.2, 1760.4, 4136.5, 4244.3],
            [a * 219474.63 for a in mwh.frequencies],
            places=1)

        # projected hessian must contain square of frequency on the diagonal
        projected_h = h.project_over_normal_modes(mwh)
        self.assertEqual(projected_h.representation.representation(), 'NN')  # change type

        for i in range(h.spacial_dof):
            self.assertAlmostEqual(
                math.fabs(projected_h.components[i, i]), mwh.frequencies[i] ** 2, places=10)

        # check output
        self.assertIn('dQ(1)', projected_h.to_string())
        self.assertNotIn('dQ(1)', projected_h.to_string(skip_trans_plus_rot_dof=6))


class TensorNumDiff(QcipToolsTestCase):

    def test_numerical_differentiation(self):
        energy = 150.
        mu = factories.FakeElectricDipole()
        alpha = factories.FakePolarizabilityTensor(input_fields=(0,))
        beta = factories.FakeFirstHyperpolarizabilityTensor(input_fields=(0, 0))

        min_field = 0.004
        k_max = 5
        ratio = 2.

        def energy_exp(fields, h0, basis, inverse, component, **kwargs):
            """Taylor series of the energy"""

            r_field = numerical_differentiation.real_fields(fields, h0, ratio)

            x = energy
            x += numpy.tensordot(mu.components, r_field, axes=1)
            x += 1 / 2 * numpy.tensordot(numpy.tensordot(alpha.components, r_field, axes=1), r_field, axes=1)
            x += 1 / 6 * numpy.tensordot(
                numpy.tensordot(numpy.tensordot(beta.components, r_field, axes=1), r_field, axes=1), r_field, axes=1)

            return x

        def dipole_exp(fields, h0, basis, inverse, component, **kwargs):
            """Taylor series of the dipole moment"""

            r_field = numerical_differentiation.real_fields(fields, h0, ratio)

            x = mu.components.copy()
            x += numpy.tensordot(alpha.components, r_field, axes=1)
            x += 1 / 2 * numpy.tensordot(numpy.tensordot(beta.components, r_field, axes=1), r_field, axes=1)

            return x[component]

        # compute polarizability
        t, triangles = derivatives.compute_numerical_derivative_of_tensor(
            derivatives.Derivative(from_representation='F'),
            derivatives.Derivative('F'),
            dipole_exp, k_max, min_field, ratio)

        self.assertArraysAlmostEqual(alpha.components, t.components, places=3)

        t, triangles = derivatives.compute_numerical_derivative_of_tensor(
            derivatives.Derivative(from_representation=''),
            derivatives.Derivative('FF'),
            energy_exp, k_max, min_field, ratio)

        self.assertArraysAlmostEqual(alpha.components, t.components, places=3)

        # compute first polarizability
        t, triangles = derivatives.compute_numerical_derivative_of_tensor(
            derivatives.Derivative(from_representation='F'),
            derivatives.Derivative('FF'),
            dipole_exp, k_max, min_field, ratio)

        self.assertArraysAlmostEqual(beta.components, t.components, delta=.001)

        t, triangles = derivatives.compute_numerical_derivative_of_tensor(
            derivatives.Derivative(from_representation=''),
            derivatives.Derivative('FFF'),
            energy_exp, k_max, min_field, ratio)

        self.assertArraysAlmostEqual(beta.components, t.components, delta=.01)
