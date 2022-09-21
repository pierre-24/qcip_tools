import subprocess
from tests import QcipScriptsTestCase


class ElectricalDerivativesTestCase(QcipScriptsTestCase):
    def setUp(self):
        # Dalton:
        self.dalton_archive_output_no_RSP = self.copy_to_temporary_directory(
            'properties/electrical_derivatives/dalton_output_no_RSP.tar.gz')
        self.dalton_archive_output_cc = self.copy_to_temporary_directory(
            'properties/electrical_derivatives/dalton_output_CC.tar.gz')

    def test_command(self):
        """Test the command on dalton archive"""

        # TODO: to check if the values are good, compare with outputs.

        # file without derivatives
        process = self.run_python_script(
            'qcip_tools/scripts/electrical_derivatives.py',
            [self.dalton_archive_output_no_RSP],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertIsNotNone(stderr_t)
        self.assertNotEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertIsNotNone(stdout_t)
        self.assertEqual(len(stdout_t), 0)

        stderr_decoded = stderr_t.decode()
        self.assertIn('cannot find electrical derivatives', stderr_decoded)

        # files that works
        process = self.run_python_script(
            'qcip_tools/scripts/electrical_derivatives.py',
            [self.dalton_archive_output_cc],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertIsNotNone(stderr_t)
        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertIsNotNone(stdout_t)
        self.assertNotEqual(len(stdout_t), 0)

        stdout_decoded = stdout_t.decode()

        self.assertIn('dipole', stdout_decoded)
        self.assertIn('alpha(0;0), w=static', stdout_decoded)
        self.assertIn('alpha(-w;w), w=1064.0nm', stdout_decoded)
        self.assertIn('beta(0;0,0), w=static', stdout_decoded)
        self.assertIn('beta(-w;w,0), w=1064.0nm', stdout_decoded)
        self.assertIn('beta(0;w,-w), w=1064.0nm', stdout_decoded)
        self.assertIn('beta(-2w;w,w), w=1064.0nm', stdout_decoded)
        self.assertIn('B_||', stdout_decoded)  # it did use the value of the dipole moment
        self.assertIn('gamma(0;0,0,0), w=static', stdout_decoded)
        self.assertIn('gamma(-w;0,0,w), w=1064.0nm', stdout_decoded)
        self.assertIn('gamma(-2w;w,w,0), w=1064.0nm', stdout_decoded)
        self.assertIn('gamma(-w;w,w,-w), w=1064.0nm', stdout_decoded)
        self.assertIn('gamma(-3w;w,w,w), w=1064.0nm', stdout_decoded)
