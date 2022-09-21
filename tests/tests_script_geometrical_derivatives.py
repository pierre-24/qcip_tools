import subprocess
from tests import QcipScriptsTestCase


class GeometricalDerivativesTestCase(QcipScriptsTestCase):
    def setUp(self):
        # gaussian
        self.gaussian_fchk_no_GG = self.copy_to_temporary_directory(
            'properties/geometrical_derivatives/gaussian_output_no_GG.fchk')
        self.gaussian_fchk = self.copy_to_temporary_directory(
            'properties/geometrical_derivatives/gaussian_output.fchk')

    def test_command(self):
        """Test the command on dalton archive"""

        # TODO: to check if the values are good, compare with outputs.

        # without hessian
        process = self.run_python_script(
            'qcip_tools/scripts/geometrical_derivatives.py',
            [self.gaussian_fchk_no_GG],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertIsNotNone(stderr_t)
        self.assertNotEqual(len(stderr_t), 0)
        self.assertIsNotNone(stdout_t)
        self.assertEqual(len(stdout_t), 0, msg=stdout_t.decode())

        stdout_decoded = stdout_t.decode()

        self.assertNotIn('gradient', stdout_decoded)
        self.assertNotIn('hessian', stdout_decoded)

        # with vibrational analysis
        process = self.run_python_script(
            'qcip_tools/scripts/geometrical_derivatives.py',
            [self.gaussian_fchk],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertIsNotNone(stderr_t)
        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertIsNotNone(stdout_t)
        self.assertNotEqual(len(stdout_t), 0)

        stdout_decoded = stdout_t.decode()

        self.assertIn('gradient', stdout_decoded)
        self.assertIn('hessian', stdout_decoded)
        self.assertIn('Total number of normal modes: 15', stdout_decoded)
        self.assertIn('Number of vibrational modes:  9', stdout_decoded)

        # without vibrational analysis
        process = self.run_python_script(
            'qcip_tools/scripts/geometrical_derivatives.py',
            [self.gaussian_fchk, '-N'],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertIsNotNone(stderr_t)
        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertIsNotNone(stdout_t)
        self.assertNotEqual(len(stdout_t), 0)

        stdout_decoded = stdout_t.decode()

        self.assertIn('gradient', stdout_decoded)
        self.assertIn('hessian', stdout_decoded)
        self.assertNotIn('Total number of normal modes: 15', stdout_decoded)
        self.assertNotIn('Number of vibrational modes:  9', stdout_decoded)
