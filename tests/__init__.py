import unittest
import numpy
import os


def array_almost_equals(a, b, places=7, delta=None, msg=''):
    """Check if two arrays containing float number are almost equals"""

    atol = 10 ** (-places)

    if delta is not None:
        atol = delta

    return numpy.testing.assert_allclose(a, b, atol=atol, err_msg=msg)


class QcipToolsTestCase(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.test_directory = os.path.dirname(__file__)

    def assertArrayAlmostEqual(self, a, b, places=3, delta=None, msg=''):
        return array_almost_equals(a, b, places=places, delta=delta, msg=msg)
