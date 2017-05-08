import unittest
import numpy


def array_almost_equals(a, b, places=7, delta=None, msg=''):
    """Check if two arrays containing float number are almost equals"""

    atol = 10 ** (-places)

    if delta is not None:
        atol = delta

    return numpy.testing.assert_allclose(a, b, atol=atol, err_msg=msg)


class QcipToolsTestCase(unittest.TestCase):

    def assertArrayAlmostEqual(self, a, b, places=3, delta=None, msg=''):
        return array_almost_equals(a, b, places=places, delta=delta, msg=msg)
