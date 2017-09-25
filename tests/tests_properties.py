import math

from tests import QcipToolsTestCase
from qcip_tools import derivatives_g, quantities
from qcip_tools.chemistry_files import gaussian


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

        copied_path = self.copy_to_temporary_directory(path)
        fx = klass()

        with open(copied_path) as f:
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

    def test_electrical_derivatives(self):
        """Test electrical properties"""

        # in FCHK:
        fchk_file, path, electrical_derivatives = self.get_property(
            gaussian.FCHK, 'electrical_derivatives/gaussian_output.fchk', 'electrical_derivatives')

        self.assertIn('F', electrical_derivatives)
        self.assertIn('FF', electrical_derivatives)
        self.assertIn('FFF', electrical_derivatives)
        self.assertIn('FD', electrical_derivatives)
        self.assertIn('FDF', electrical_derivatives)
        self.assertIn('FDD', electrical_derivatives)

        self.assertIn(0.02, electrical_derivatives['FDF'])
        self.assertIn(0.04, electrical_derivatives['FDF'])
        self.assertIn(0.06, electrical_derivatives['FDF'])
        self.assertIn(0.08, electrical_derivatives['FDF'])
        self.assertIn(0.10, electrical_derivatives['FDF'])

        self.assertAlmostEqual(electrical_derivatives['FD'][0.02].isotropic_value(), 0.835791e1, places=5)
        self.assertAlmostEqual(electrical_derivatives['FDD'][0.02].beta_hrs(), 6.2390, places=3)

    def test_geometrical_derivatives(self):
        """Test geometrical properties"""

        HartreeToWavenumber = quantities.convert(quantities.ureg.hartree, quantities.ureg.wavenumber)

        # in FCHK:
        fchk_file, path, geometrical_derivatives = self.get_property(
            gaussian.FCHK, 'geometrical_derivatives/gaussian_output.fchk', 'geometrical_derivatives')

        self.assertIn('G', geometrical_derivatives)
        self.assertIn('GG', geometrical_derivatives)

        mwh = derivatives_g.MassWeightedHessian(fchk_file.molecule, geometrical_derivatives['GG'])
        frequencies_in_wavenumber = [a * HartreeToWavenumber for a in mwh.frequencies]
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
