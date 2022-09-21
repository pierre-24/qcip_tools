import subprocess
from tests import QcipScriptsTestCase

from qcip_tools import quantities, derivatives_exci
from qcip_tools.chemistry_files import dalton


class ExcitationsTestCase(QcipScriptsTestCase):
    def setUp(self):
        # dalton
        self.exc = self.copy_to_temporary_directory('properties/excitations/dalton_output.tar.gz')

    def test_command(self):
        """Test the command on dalton output"""
        process = self.run_python_script(
            'qcip_tools/scripts/excitations.py',
            [self.exc, '-E', 'eV'],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertIsNotNone(stderr_t)
        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertIsNotNone(stdout_t)
        self.assertNotEqual(len(stdout_t), 0)

        lines = stdout_t.decode().splitlines()

        # get excitations:
        with open(self.exc, 'rb') as f:
            fx = dalton.ArchiveOutput()
            fx.read(f)

        p = fx.property('excitations')
        self.assertEqual(len(lines), p['!'].representation.nstates - 1 + 4)

        excis = derivatives_exci.Excitations(p['!'], p['!F'])

        # check if output is correct
        for i in range(3):
            exploded = lines[i + 3].split()
            self.assertAlmostEqual(
                float(exploded[1]),
                excis.transition_energy(i + 1) * quantities.convert(quantities.ureg.hartree, quantities.ureg.eV),
                places=4)

        # check with limit
        process = self.run_python_script(
            'qcip_tools/scripts/excitations.py',
            [self.exc, '-E', 'eV', '-l', '0.01'],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertIsNotNone(stderr_t)
        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertIsNotNone(stdout_t)
        self.assertNotEqual(len(stdout_t), 0)

        lines = stdout_t.decode().splitlines()

        line = 0
        for i in range(3):
            exploded = lines[line + 3].split()
            if excis.oscillator_strength(i + 1) > .01:
                self.assertAlmostEqual(
                    float(exploded[1]),
                    excis.transition_energy(i + 1) * quantities.convert(quantities.ureg.hartree, quantities.ureg.eV),
                    places=4)
                line += 1

        self.assertEqual(line, 3)  # skipped excitation #1 due to f_ge < .01


class GenSpectrumTestCase(QcipScriptsTestCase):
    def setUp(self):

        # dalton
        self.exc = self.copy_to_temporary_directory('properties/excitations/dalton_output.tar.gz')

    def test_command(self):
        """Test the command on dalton output"""

        process = self.run_python_script(
            'qcip_tools/scripts/gen_spectrum.py',
            ['-I', '-l', '0:20:eV', self.exc, 'uv'],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertIsNotNone(stderr_t)
        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertIsNotNone(stdout_t)
        self.assertNotEqual(len(stdout_t), 0)

        lines = stdout_t.decode().splitlines()

        # get excitations:
        with open(self.exc, 'rb') as f:
            fx = dalton.ArchiveOutput()
            fx.read(f)

        p = fx.property('excitations')
        excis = derivatives_exci.Excitations(p['!'], p['!F'])

        # check impulses
        impulses = lines[-excis.nstates:]
        self.assertIn('impulses', impulses[0])
        for i in range(1, excis.nstates):
            inf = impulses[i].split()
            self.assertAlmostEqual(float(inf[2]), excis.oscillator_strength(i), places=4)
            self.assertAlmostEqual(
                float(inf[1]),
                excis.transition_energy(i) * quantities.convert(quantities.ureg.hartree, quantities.ureg.eV),
                places=4)
