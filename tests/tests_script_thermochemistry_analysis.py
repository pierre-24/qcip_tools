import subprocess
from tests import QcipScriptsTestCase


class ThermochemistryAnalysisTestCase(QcipScriptsTestCase):
    def setUp(self):
        # gaussian
        self.gaussian_fchk_no_GG = self.copy_to_temporary_directory(
            'properties/geometrical_derivatives/gaussian_output_no_GG.fchk')
        self.gaussian_fchk = self.copy_to_temporary_directory(
            'properties/geometrical_derivatives/gaussian_output.fchk')
        self.gaussian_fchk_TS = self.copy_to_temporary_directory(
            'properties/thermochemistry_analysis/gaussian_output_TS.fchk')

    def test_command(self):
        """Test the command on gaussian stuffs"""

        # TODO: to check if the values are good, compare with outputs.

        # without hessian
        process = self.run_python_script(
            'qcip_tools/scripts/thermochemistry_analysis.py',
            [self.gaussian_fchk_no_GG],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertIsNotNone(stderr_t)
        self.assertNotEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertIsNotNone(stdout_t)
        self.assertEqual(len(stdout_t), 0)

        stdout_decoded = stdout_t.decode()

        self.assertNotIn('G(T)', stdout_decoded)

        # with vibrational analysis
        process = self.run_python_script(
            'qcip_tools/scripts/thermochemistry_analysis.py',
            [self.gaussian_fchk, '-n', '12'],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertIsNotNone(stderr_t)
        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertIsNotNone(stdout_t)
        self.assertNotEqual(len(stdout_t), 0)

        stdout_decoded = stdout_t.decode()

        self.assertIn('G(T)', stdout_decoded)
        G_line = stdout_decoded.splitlines()[-2].split()
        self.assertAlmostEqual(float(G_line[1]), -40.173859, places=5)

        # with vibrational analysis and guessing symmetry
        process = self.run_python_script(
            'qcip_tools/scripts/thermochemistry_analysis.py',
            [self.gaussian_fchk, '-g'],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertIsNotNone(stderr_t)
        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertIsNotNone(stdout_t)
        self.assertNotEqual(len(stdout_t), 0)

        stdout_decoded = stdout_t.decode()

        self.assertIn('G(T)', stdout_decoded)
        G_line = stdout_decoded.splitlines()[-2].split()
        self.assertAlmostEqual(float(G_line[1]), -40.173859, places=5)

        # TS: exclude first frequency
        process = self.run_python_script(
            'qcip_tools/scripts/thermochemistry_analysis.py',
            [self.gaussian_fchk_TS, '-x', '7'],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertIsNotNone(stderr_t)
        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertIsNotNone(stdout_t)
        self.assertNotEqual(len(stdout_t), 0)

        stdout_decoded = stdout_t.decode()

        self.assertIn('G(T)', stdout_decoded)
        G_line = stdout_decoded.splitlines()[-2].split()
        self.assertAlmostEqual(float(G_line[1]), -598.517192, places=5)
