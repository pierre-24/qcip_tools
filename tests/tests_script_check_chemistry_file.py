import subprocess
from tests import QcipScriptsTestCase


class CheckChemistryFileTestCase(QcipScriptsTestCase):
    def setUp(self):
        self.dalton_archive_output_cc = self.copy_to_temporary_directory(
            'properties/electrical_derivatives/dalton_output_CC.tar.gz')

    def command_test(self, file, property='', expected_return=1):
        args = [file]

        if property:
            args.extend(['-p', property])

        process = self.run_python_script(
            'qcip_tools/scripts/check_chemistry_file.py',
            args,
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertIsNotNone(stderr_t)
        # self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertIsNotNone(stdout_t)
        self.assertNotEqual(len(stdout_t), 0)
        self.assertEqual(stdout_t, '{}\n'.format(expected_return).encode())

    def test_command(self):
        """Test the command on dalton files"""

        # Dalton archive
        self.command_test(self.dalton_archive_output_cc)
        self.command_test(self.dalton_archive_output_cc, property='file_type')
        self.command_test(self.dalton_archive_output_cc, property='xxxx', expected_return=0)
