import math
import os
import numpy

from tests import QcipToolsTestCase
from qcip_tools import derivatives_g, derivatives, derivatives_e
from qcip_tools.chemistry_files import gaussian, dalton, PropertyNotPresent


class PropertiesTestCase(QcipToolsTestCase):
    """test the property() function"""

    def get_property(self, klass, path, property):
        """Get a given property in a file represented by a class.

        :param klass: the object
        :type klass: object
        :param path: path to the file
        :type path: str
        :param property: property to get
        :type property: str
        :rtype: str|dict
        """

        copied_path = self.copy_to_temporary_directory(os.path.join('properties', path))
        fx = klass()

        with open(copied_path, 'r' if not klass.requires_binary_mode else 'rb') as f:
            fx.read(f)
        p = fx.property(property)

        return fx, copied_path, p

    def test_geometry(self):
        """Test the geometry"""

        # in FCHK:
        fchk_file, path, geometry = self.get_property(
            gaussian.FCHK, 'electrical_derivatives/gaussian_output.fchk', 'molecule')

        self.assertEqual(fchk_file.molecule, geometry)

        # ... And one is enough (as long as it derives from WithMoleculeMixin)

    def test_energies(self):
        """Test energy"""

        # 1. In FCHK:
        fchk_file, path, energies = self.get_property(
            gaussian.FCHK, 'computed_energies/gaussian_output.fchk', 'computed_energies')

        found_energies = ['SCF/DFT', 'total', 'HF', 'CCSD', 'CCSD(T)', 'MP2', 'MP3', 'MP4D', 'MP4DQ', 'MP4SDQ']
        for e in found_energies:
            self.assertIn(e, energies)

        fchk_file, path, energies = self.get_property(
            gaussian.FCHK, 'computed_energies/gaussian_output_BLYP.fchk', 'computed_energies')

        found_energies = ['SCF/DFT', 'total', 'BLYP']
        for e in found_energies:
            self.assertIn(e, energies)

        # 2. In dalton archive
        archive_file, path, energies = self.get_property(
            dalton.ArchiveOutput, 'computed_energies/dalton_output.tar.gz', 'computed_energies')

        found_energies = ['total', 'SCF/DFT']
        for e in found_energies:
            self.assertIn(e, energies)
        archive_file, path, energies = self.get_property(
            dalton.ArchiveOutput, 'computed_energies/dalton_output_CC.tar.gz', 'computed_energies')

        found_energies = ['total', 'CCSD']
        for e in found_energies:
            self.assertIn(e, energies)

        # 3. In dalton output
        archive_file, path, energies = self.get_property(
            dalton.Output, 'computed_energies/dalton_output_CC.out', 'computed_energies')

        found_energies = ['total', 'SCF/DFT', 'MP2', 'CCSD']
        for e in found_energies:
            self.assertIn(e, energies)

    def test_electrical_derivatives(self):
        """Test electrical properties"""

        # 1. In FCHK:
        fchk_file, path, electrical_derivatives = self.get_property(
            gaussian.FCHK, 'electrical_derivatives/gaussian_output.fchk', 'electrical_derivatives')

        self.assertIn('F', electrical_derivatives)
        self.assertIn('FF', electrical_derivatives)
        self.assertIn('FFF', electrical_derivatives)
        self.assertIn('dD', electrical_derivatives)
        self.assertIn('dDF', electrical_derivatives)
        self.assertIn('XDD', electrical_derivatives)

        self.assertIn(0.02, electrical_derivatives['dDF'])
        self.assertIn(0.04, electrical_derivatives['dDF'])
        self.assertIn(0.06, electrical_derivatives['dDF'])
        self.assertIn(0.08, electrical_derivatives['dDF'])
        self.assertIn(0.10, electrical_derivatives['dDF'])

        self.assertAlmostEqual(electrical_derivatives['dD'][0.02].isotropic_value(), 0.835791e1, places=5)
        self.assertAlmostEqual(electrical_derivatives['XDD'][0.02].beta_hrs(), 6.2390, places=3)

        # 2. In dalton archive:
        f = 0.0428226504

        # Coupled Cluster
        archive_file, path, electrical_derivatives = self.get_property(
            dalton.ArchiveOutput,
            'electrical_derivatives/dalton_output_CC.tar.gz',
            'electrical_derivatives')

        self.assertIn('F', electrical_derivatives)
        self.assertIn('FF', electrical_derivatives)
        self.assertIn('dD', electrical_derivatives)
        self.assertIn('dDF', electrical_derivatives)
        self.assertIn('FFF', electrical_derivatives)
        self.assertIn('XDD', electrical_derivatives)
        self.assertIn('FFFF', electrical_derivatives)
        self.assertIn('dFFD', electrical_derivatives)
        self.assertIn('XDDF', electrical_derivatives)
        self.assertIn('XDDD', electrical_derivatives)
        self.assertIn('dDDd', electrical_derivatives)

        self.assertIn(f, electrical_derivatives['dD'])

        tests_in_tensor = [
            ('F', 'static', (2,), 0.64738),
            ('FF', 'static', (0, 0), 0.04757),
            ('dD', f, (0, 0), 0.047968),
            ('FFF', 'static', (0, 0, 2), 0.073362),
            ('dDF', f, (0, 0, 2), 0.072835),
            ('XDD', f, (0, 0, 2), 0.073591),
            ('FFFF', 'static', (0, 0, 0, 0), -0.351403),
            ('dFFD', f, (0, 0, 0, 0), -0.35619),
            ('XDDF', f, (0, 0, 0, 0), -0.36612),
            ('dDDd', f, (0, 0, 0, 0), -0.36104),
            ('XDDD', f, (0, 0, 0, 0), -0.382162)
        ]

        for(tensor, freq, coo, value) in tests_in_tensor:
            self.assertAlmostEqual(electrical_derivatives[tensor][freq].components[coo], value, places=4)

            d = derivatives.Derivative(tensor)
            for i in d.smart_iterator():
                v = electrical_derivatives[tensor][freq].components[i]
                for j in d.inverse_smart_iterator(i):
                    self.assertArraysAlmostEqual(electrical_derivatives[tensor][freq].components[j], v)

        # Response alpha:
        archive_file, path, electrical_derivatives = self.get_property(
            dalton.ArchiveOutput,
            'electrical_derivatives/dalton_output_RSP_alpha.tar.gz',
            'electrical_derivatives')

        self.assertIn('F', electrical_derivatives)
        self.assertIn('FF', electrical_derivatives)
        self.assertIn('dD', electrical_derivatives)

        self.assertNotIn('FFF', electrical_derivatives)
        self.assertNotIn('FFFF', electrical_derivatives)

        self.assertIn(f, electrical_derivatives['dD'])

        tests_in_tensor = [
            ('F', 'static', (2,), 0.67340),
            ('FF', 'static', (2, 2), 2.0408),
            ('dD', f, (2, 2), 2.0459)
        ]

        for(tensor, freq, coo, value) in tests_in_tensor:
            self.assertAlmostEqual(electrical_derivatives[tensor][freq].components[coo], value, places=4, msg=tensor)

            d = derivatives.Derivative(tensor)
            for i in d.smart_iterator():
                v = electrical_derivatives[tensor][freq].components[i]
                for j in d.inverse_smart_iterator(i):
                    self.assertArraysAlmostEqual(electrical_derivatives[tensor][freq].components[j], v)

        # Response beta:
        archive_file, path, electrical_derivatives = self.get_property(
            dalton.ArchiveOutput,
            'electrical_derivatives/dalton_output_RSP_beta.tar.gz',
            'electrical_derivatives')

        self.assertIn('F', electrical_derivatives)
        self.assertIn('FFF', electrical_derivatives)
        self.assertIn('dDF', electrical_derivatives)
        self.assertIn('XDD', electrical_derivatives)

        self.assertNotIn('FF', electrical_derivatives)
        self.assertNotIn('FFFF', electrical_derivatives)

        old_f = f
        f = 0.042823
        self.assertIn(f, electrical_derivatives['XDD'])

        tests_in_tensor = [
            ('F', 'static', (2,), 0.67340),
            ('FFF', 'static', (2, 2, 2), -3.28102875),
            ('dDF', f, (2, 2, 2), -3.30365796),
            ('XDD', f, (2, 2, 2), -3.34979625)
        ]

        for(tensor, freq, coo, value) in tests_in_tensor:
            self.assertAlmostEqual(electrical_derivatives[tensor][freq].components[coo], value, places=4, msg=tensor)

            d = derivatives.Derivative(tensor)
            for i in d.smart_iterator():
                v = electrical_derivatives[tensor][freq].components[i]
                for j in d.inverse_smart_iterator(i):
                    self.assertArraysAlmostEqual(electrical_derivatives[tensor][freq].components[j], v)

        # Response beta (patched version, in DALTON.PROP):
        archive_file, path, electrical_derivatives = self.get_property(
            dalton.ArchiveOutput,
            'electrical_derivatives/dalton_output_RSP_beta_patched.tar.gz',
            'electrical_derivatives')

        self.assertIn('F', electrical_derivatives)
        self.assertIn('FFF', electrical_derivatives)
        self.assertIn('dDF', electrical_derivatives)
        self.assertIn('XDD', electrical_derivatives)

        self.assertNotIn('FF', electrical_derivatives)
        self.assertNotIn('FFFF', electrical_derivatives)

        f = old_f
        self.assertIn(f, electrical_derivatives['XDD'])

        tests_in_tensor = [
            ('F', 'static', (2,), 0.44317),
            ('FFF', 'static', (2, 2, 2), -9.76387326),
            ('dDF', f, (2, 2, 2), -9.82193674),
            ('XDD', f, (2, 2, 2), -9.94003880)
        ]

        for (tensor, freq, coo, value) in tests_in_tensor:
            self.assertAlmostEqual(electrical_derivatives[tensor][freq].components[coo], value, places=4, msg=tensor)

            d = derivatives.Derivative(tensor)
            for i in d.smart_iterator():
                v = electrical_derivatives[tensor][freq].components[i]
                for j in d.inverse_smart_iterator(i):
                    self.assertArraysAlmostEqual(electrical_derivatives[tensor][freq].components[j], v)

        # Response gamma:
        archive_file, path, electrical_derivatives = self.get_property(
            dalton.ArchiveOutput,
            'electrical_derivatives/dalton_output_RSP_gamma.tar.gz',
            'electrical_derivatives')

        self.assertIn('F', electrical_derivatives)
        self.assertIn('FFFF', electrical_derivatives)
        self.assertIn('dDFF', electrical_derivatives)
        self.assertIn('XDDF', electrical_derivatives)
        self.assertIn('dDDd', electrical_derivatives)
        self.assertIn('XDDD', electrical_derivatives)

        self.assertNotIn('FF', electrical_derivatives)
        self.assertNotIn('FFF', electrical_derivatives)

        self.assertIn(f, electrical_derivatives['XDDD'])

        tests_in_tensor = [
            ('F', 'static', (2,), 0.44317),
            ('FFFF', 'static', (1, 0, 0, 0), -12.59711626),
            ('dDFF', f, (1, 0, 0, 0), -12.72601993),
            ('XDDF', f, (1, 0, 0, 0), -13.00065914),
            ('dDDd', f, (1, 0, 0, 0), -12.85136396),
            ('XDDD', f, (1, 0, 0, 0), -13.44430492),
        ]

        for(tensor, freq, coo, value) in tests_in_tensor:
            self.assertAlmostEqual(electrical_derivatives[tensor][freq].components[coo], value, places=4, msg=tensor)

            d = derivatives.Derivative(tensor)
            for i in d.smart_iterator():
                v = electrical_derivatives[tensor][freq].components[i]
                for j in d.inverse_smart_iterator(i):
                    self.assertArraysAlmostEqual(electrical_derivatives[tensor][freq].components[j], v)

        # Response gamma (patched version):
        archive_file, path, electrical_derivatives = self.get_property(
            dalton.ArchiveOutput,
            'electrical_derivatives/dalton_output_RSP_gamma_patched.tar.gz',
            'electrical_derivatives')

        self.assertIn('F', electrical_derivatives)
        self.assertIn('FFFF', electrical_derivatives)
        self.assertIn('dDFF', electrical_derivatives)
        self.assertIn('XDDF', electrical_derivatives)
        self.assertIn('dDDd', electrical_derivatives)
        self.assertIn('XDDD', electrical_derivatives)

        self.assertNotIn('FF', electrical_derivatives)
        self.assertNotIn('FFF', electrical_derivatives)

        self.assertIn(f, electrical_derivatives['XDDD'])

        tests_in_tensor = [
            ('F', 'static', (2,), 0.44317),
            ('FFFF', 'static', (1, 0, 0, 0), -12.59711626),
            ('dDFF', f, (1, 0, 0, 0), -12.72601993),
            ('XDDF', f, (1, 0, 0, 0), -13.00065914),
            ('dDDd', f, (1, 0, 0, 0), -12.85136396),
            ('XDDD', f, (1, 0, 0, 0), -13.44430492),
        ]

        for(tensor, freq, coo, value) in tests_in_tensor:
            self.assertAlmostEqual(electrical_derivatives[tensor][freq].components[coo], value, places=4, msg=tensor)

            d = derivatives.Derivative(tensor)
            for i in d.smart_iterator():
                v = electrical_derivatives[tensor][freq].components[i]
                for j in d.inverse_smart_iterator(i):
                    self.assertArraysAlmostEqual(electrical_derivatives[tensor][freq].components[j], v)

    def test_geometrical_derivatives(self):
        """Test geometrical properties"""

        # 1. In FCHK:
        fchk_file, path, geometrical_derivatives = self.get_property(
            gaussian.FCHK, 'geometrical_derivatives/gaussian_output.fchk', 'geometrical_derivatives')

        self.assertIn('G', geometrical_derivatives)
        self.assertIn('GG', geometrical_derivatives)

        mwh = derivatives_g.MassWeightedHessian(fchk_file.molecule, geometrical_derivatives['GG'])
        frequencies_in_wavenumber = [
            a * derivatives_g.MassWeightedHessian.HARTREE_TO_WAVENUMBER_CONVERSION for a in mwh.frequencies]
        frequencies_and_occurrences = [
            (1430, 3),
            (1657, 2),
            (3147, 1),
            (3271, 3)]  # from log

        for f, occurrence in frequencies_and_occurrences:
            close_to_value = [math.fabs(fx - f) < 1 for fx in frequencies_in_wavenumber]
            self.assertIn(True, close_to_value, msg='f={}'.format(f))
            self.assertEqual(close_to_value.count(True), occurrence)

        omega_e, omega_t, omega_r, omega_v = mwh.compute_partition_functions(12)
        self.assertEqual(omega_e, 1)
        self.assertAlmostEqual(omega_t, 0.252295e7, delta=1e2)  # from log
        self.assertAlmostEqual(omega_r, 0.363253e2, places=3)  # from log
        self.assertAlmostEqual(omega_v, 0.100370e1, places=3)  # from log, "Vib (V=0)"

        self.assertAlmostEqual(mwh.compute_zpva(), 0.046857, places=5)
        self.assertAlmostEqual(mwh.compute_zpva() + sum(mwh.compute_internal_energy()), 0.049715, places=5)
        self.assertAlmostEqual(mwh.compute_zpva() + sum(mwh.compute_enthalpy()), 0.050659, places=5)
        self.assertAlmostEqual(mwh.compute_zpva() + sum(mwh.compute_gibbs_free_energy(12)), 0.029544, places=5)

        # Test on the frequency calculation of a transition state
        # Gaussian simply choose to ignore the imaginary frequency mode
        fchk_file, path, geometrical_derivatives = self.get_property(
            gaussian.FCHK, 'geometrical_derivatives/gaussian_output_TS.fchk', 'geometrical_derivatives')

        mwh = derivatives_g.MassWeightedHessian(fchk_file.molecule, geometrical_derivatives['GG'])
        self.assertTrue(mwh.frequencies[0] * derivatives_g.MassWeightedHessian.HARTREE_TO_WAVENUMBER_CONVERSION < -400)

        mwh.included_modes.pop(0)
        # Note: actually the imaginary mode is the first one, since it have the lowest frequency, so this line above
        # remove a mode which is actually rotation (but whatever).

        self.assertAlmostEqual(mwh.compute_zpva(), 0.040758, places=5)
        self.assertAlmostEqual(mwh.compute_zpva() + sum(mwh.compute_internal_energy()), 0.045152, places=5)
        self.assertAlmostEqual(mwh.compute_zpva() + sum(mwh.compute_enthalpy()), 0.046096, places=5)
        self.assertAlmostEqual(mwh.compute_zpva() + sum(mwh.compute_gibbs_free_energy(1)), 0.013912, places=5)

        # test with an atom alone
        fchk_file, path, geometrical_derivatives = self.get_property(
            gaussian.FCHK, 'geometrical_derivatives/gaussian_atom_alone.fchk', 'geometrical_derivatives')

        mwh = derivatives_g.MassWeightedHessian(fchk_file.molecule, geometrical_derivatives['GG'])
        self.assertEqual(mwh.compute_zpva(), .0)
        self.assertAlmostEqual(mwh.compute_zpva() + sum(mwh.compute_internal_energy()), 0.001416, places=5)
        self.assertAlmostEqual(mwh.compute_zpva() + sum(mwh.compute_enthalpy()), 0.002360, places=5)
        self.assertAlmostEqual(mwh.compute_zpva() + sum(mwh.compute_gibbs_free_energy(1)), -0.015023, places=5)

        # test force
        fchk_file, path, geometrical_derivatives = self.get_property(
            gaussian.FCHK, 'geometrical_derivatives/gaussian_output_force.fchk', 'geometrical_derivatives')

        self.assertIn('G', geometrical_derivatives)
        self.assertNotEqual(numpy.linalg.norm(geometrical_derivatives['G'].components), .0)

        # test with a single point (cartesian gradient is given, but null)
        fchk_file = gaussian.FCHK()
        with open(os.path.join(self.tests_files_directory, 'properties/computed_energies/gaussian_output.fchk')) as f:
            fchk_file.read(f)

        with self.assertRaises(PropertyNotPresent):
            fchk_file.property('geometrical_derivatives')

        gradient = fchk_file.get('Cartesian Gradient')
        self.assertEqual(numpy.linalg.norm(gradient), .0)

        # 2. test dalton
        out_file, path, geometrical_derivatives = self.get_property(
            dalton.Output, 'geometrical_derivatives/dalton_output_CC.out', 'geometrical_derivatives')

        self.assertIn('G', geometrical_derivatives)

        tar_file, path, geometrical_derivatives = self.get_property(
            dalton.ArchiveOutput, 'geometrical_derivatives/dalton_output.tar.gz', 'geometrical_derivatives')

        self.assertIn('GG', geometrical_derivatives)

        out_file, path, geometrical_derivatives_2 = self.get_property(
            dalton.Output, 'geometrical_derivatives/dalton_output.out', 'geometrical_derivatives')

        self.assertIn('G', geometrical_derivatives_2)
        self.assertIn('GG', geometrical_derivatives_2)

        self.assertArraysAlmostEqual(
            geometrical_derivatives_2['GG'].components, geometrical_derivatives['GG'].components, places=4)

        tar_file, path, geometrical_derivatives = self.get_property(
            dalton.ArchiveOutput, 'geometrical_derivatives/dalton_output_2.tar.gz', 'geometrical_derivatives')

        self.assertIn('GG', geometrical_derivatives)

        out_file, path, geometrical_derivatives_2 = self.get_property(
            dalton.Output, 'geometrical_derivatives/dalton_output_2.out', 'geometrical_derivatives')

        self.assertIn('G', geometrical_derivatives_2)
        self.assertIn('GG', geometrical_derivatives_2)

        self.assertArraysAlmostEqual(
            geometrical_derivatives_2['GG'].components, geometrical_derivatives['GG'].components, places=4)

    def test_excitations(self):
        """Test excitations"""

        # 1. Test gaussian
        fchk_file, path, excitations = self.get_property(
            gaussian.FCHK, 'excitations/gaussian_output.fchk', 'excitations')

        self.assertIn('!', excitations)
        self.assertIn('!F', excitations)

        self.assertEqual(excitations['!'].components[0], .0)  # ground state has no energy wrt ground state
        self.assertAlmostEqual(derivatives_e.convert_energy_to(excitations['!'].components[1], 'eV'), 4.6202, places=4)
        self.assertAlmostEqual(
            derivatives_e.ElectricDipole(dipole=excitations['!F'].components[1]).norm(), 0.4886, places=4)
