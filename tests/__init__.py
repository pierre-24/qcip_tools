import unittest
import numpy
import os
import tempfile
import shutil


def array_almost_equals(a, b, places=7, delta=None, msg=''):
    """Check if two arrays containing float number are almost equals"""

    atol = 10 ** (-places)

    if delta is not None:
        atol = delta

    return numpy.testing.assert_allclose(a, b, atol=atol, err_msg=msg)


class QcipToolsTestCase(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.tests_files_directory = os.path.join(os.path.dirname(__file__), 'tests_files')
        self.temporary_directory = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temporary_directory)

    def assertArrayAlmostEqual(self, a, b, places=3, delta=None, msg=''):
        return array_almost_equals(a, b, places=places, delta=delta, msg=msg)

    def copy_to_temporary_directory(self, path, new_name=''):
        """Copy the content of a file to the temporary directory

        :param path: path to the file to copy
        :type path: str
        :param new_name: the new name of the file in the temporary directory (if blank, the one from path is used)
        :type new_name: str
        :rtype: str
        """

        path_in_test = os.path.join(self.tests_files_directory, path)

        if not os.path.exists(path_in_test):
            raise FileNotFoundError(path_in_test)

        if not new_name:
            new_name = os.path.basename(path)

        path_in_temp = os.path.join(self.temporary_directory, new_name)

        if os.path.exists(path_in_temp):
            raise FileExistsError(path_in_temp)

        with open(path_in_temp, 'wb') as f:
            with open(os.path.join(self.tests_files_directory, path), 'rb') as fx:
                f.write(fx.read())

        return path_in_temp
