import subprocess
from tests import QcipScriptsTestCase

from qcip_tools.scripts import commons # noqa
from qcip_tools.chemistry_files import gaussian


class ElectricalDerivativesTestCase(QcipScriptsTestCase):
    def setUp(self):
        # Dalton:
        self.dalton_archive_output_no_RSP = self.copy_to_temporary_directory(
            'properties/electrical_derivatives/dalton_output_no_RSP.tar.gz')
        self.dalton_archive_output_cc = self.copy_to_temporary_directory(
            'properties/electrical_derivatives/dalton_output_CC.tar.gz')

        # gaussian
        self.gaussian_fchk = self.copy_to_temporary_directory('properties/electrical_derivatives/gaussian_output.fchk')
        self.gaussian_log = self.copy_to_temporary_directory('properties/electrical_derivatives/gaussian_output.log')

    def test_gaussian_log(self):

        fchk_file = gaussian.FCHK()
        with open(self.gaussian_fchk) as f:
            fchk_file.read(f)

        log_file = gaussian.Output()
        with open(self.gaussian_log) as f:
            log_file.read(f)

        log_file_ed = log_file.property('electrical_derivatives')
        fchk_file_ed = fchk_file.property('electrical_derivatives')

        test_equals = [
            ('FF', 'static', 'static'),
            ('FFF', 'static', 'static'),
            ('dDF', 0.04, 0.04),
            ('XDD', 0.04, 0.04)
        ]

        for prop, freq1, freq2 in test_equals:
            self.assertArraysAlmostEqual(
                log_file_ed[prop][freq1].components,
                fchk_file_ed[prop][freq2].components)

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
