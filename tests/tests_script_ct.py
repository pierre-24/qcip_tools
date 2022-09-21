import subprocess
from tests import QcipScriptsTestCase


class CTAnalysisTestCase(QcipScriptsTestCase):
    def setUp(self):
        # gaussian
        self.gs = self.copy_to_temporary_directory('ct_analysis/gaussian_gs.cub')
        self.es = self.copy_to_temporary_directory('ct_analysis/gaussian_es.cub')

    def test_command(self):
        """Test the command on gaussian output"""
        process = self.run_python_script(
            'qcip_tools/scripts/ct_analysis.py',
            ['-g', self.gs, '-e', self.es],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertIsNotNone(stderr_t)
        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertIsNotNone(stdout_t)
        self.assertNotEqual(len(stdout_t), 0)

        lines = stdout_t.decode().splitlines()
        self.assertAlmostEqual(float(lines[-3].split()[-1]), 0.742, places=3)  # q_CT
        self.assertAlmostEqual(float(lines[-2].split()[-1]), 0.411, places=3)  # d_CT


class RadialDistributionTestCase(QcipScriptsTestCase):
    def setUp(self):
        self.cube = self.copy_to_temporary_directory('ct_analysis/H_1s.cub')

    def test_command(self):
        """Test the command on gaussian output"""

        process = self.run_python_script(
            'qcip_tools/scripts/cube_radial_distribution.py',
            [self.cube],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertIsNotNone(stderr_t)
        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertIsNotNone(stdout_t)
        self.assertNotEqual(len(stdout_t), 0)

        lines = stdout_t.decode().splitlines()

        # mean radius
        self.assertAlmostEqual(float(lines[-1].split()[-1]), .755, places=3)

        # sums of electrons (1 electron in the density)
        self.assertAlmostEqual(float(lines[-2].split()[-1]), 1., places=2)
