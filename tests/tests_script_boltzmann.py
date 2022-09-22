import subprocess
from tests import QcipScriptsTestCase


class BoltzmannTestCase(QcipScriptsTestCase):
    def setUp(self):
        self.conformers = [
            self.copy_to_temporary_directory('conformers/confo1.fchk'),
            self.copy_to_temporary_directory('conformers/confo2.fchk'),
            self.copy_to_temporary_directory('conformers/confo3.fchk')
        ]

    def command(self, confos, T: float = 298.15, criterion: str = 'E'):
        process = self.run_python_script(
            'qcip_tools/scripts/boltzmann_population.py',
            confos + ['-t', str(T), '-c', criterion],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertIsNotNone(stderr_t)
        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertIsNotNone(stdout_t)
        self.assertNotEqual(len(stdout_t), 0)

        return stdout_t.decode()

    def test_command(self):
        """Test the command on gaussian outputs"""

        # one conformer: 100%
        out = self.command([self.conformers[0]])
        self.assertAlmostEqual(float(out.splitlines()[-1].split()[-2]), 100.)

        # low temperature: 100% of the most stable
        most_stable = -1
        out = self.command(self.conformers, T=10)
        self.assertAlmostEqual(float(out.splitlines()[most_stable].split()[-2]), 100.)

        # large temperature: 33.33333% for all
        out = self.command(self.conformers, T=1e7)
        lines = out.splitlines()
        for i in [-1, -2, -3]:
            self.assertAlmostEqual(float(lines[i].split()[-2]), 33.3, places=1)

        # use "G":
        out1 = float(self.command(self.conformers, criterion='E').splitlines()[most_stable].split()[-2])
        out2 = float(self.command(self.conformers, criterion='G').splitlines()[most_stable].split()[-2])
        self.assertTrue(out1 > out2)  # takes entropy (always positive) into account, so it lowers the energy diff
