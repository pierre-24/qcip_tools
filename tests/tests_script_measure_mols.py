import subprocess
from tests import QcipScriptsTestCase

from qcip_tools.chemistry_files import helpers
from qcip_tools import math as qm


class MeasureMolsTestCase(QcipScriptsTestCase):
    def setUp(self):
        self.gaussian_fchk = self.copy_to_temporary_directory('properties/geometrical_derivatives/gaussian_output.fchk')

    def test_command(self):
        """Test the command on gaussian output"""

        criteria = [(1, 2), (2, 1, 3), (2, 1, 3, 4)]

        process = self.run_python_script(
            'qcip_tools/scripts/measure_mols.py',
            [self.gaussian_fchk, '-c', ','.join('-'.join(str(i) for i in c) for c in criteria)],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertIsNotNone(stderr_t)
        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertIsNotNone(stdout_t)
        self.assertNotEqual(len(stdout_t), 0)

        lines = stdout_t.decode().splitlines()

        with open(self.gaussian_fchk) as f:
            fp = helpers.open_chemistry_file(f)
        mol = fp.molecule

        values = lines[-1].split()[1:]
        for i, c in enumerate(criteria):
            coos = [mol[j - 1].position for j in c]
            self.assertAlmostEqual(
                float(values[i]),
                {
                    2: qm.distance,
                    3: qm.angle,
                    4: qm.torsion_angle
                }[len(c)](*coos), places=4)
