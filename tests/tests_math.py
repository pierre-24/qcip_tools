import numpy

from tests import QcipToolsTestCase

from qcip_tools import math as qcip_math


class MathTestCase(QcipToolsTestCase):

    def setUp(self):
        pass

    def test_vector_manipulation(self):
        """Test the behavior of vector functions"""

        # test basic stuffs:
        self.assertTrue(numpy.array_equal(qcip_math.normalize([-3, 0, 0]), [-1, 0, 0.]))
        self.assertEqual(qcip_math.distance([-1, 0, 0], [.5, 0, 0]), 1.5)
        self.assertEqual(qcip_math.angle([2, 0, 0], [0, 0, 0], [0, 1, 0]), 90.0)
        self.assertEqual(qcip_math.torsion_angle([2, 0, 0], [0, 0, 0], [0, 1, 0], [-2, 1, 0]), -180.0)

    def test_bla(self):
        """Test the BLA definition"""

        # test BLA:
        position_list_1 = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (2, 1, 0)]
        self.assertEqual(qcip_math.BLA(position_list_1), 0.0)  # bonds length are equal, BLA=0

        position_list_2 = [(0, 0, 0), (1, 0, 0), (1, .5, 0), (2, .5, 0)]
        self.assertEqual(qcip_math.BLA(position_list_2), .5)  # first bond is long, so BLA > 0

        position_list_3 = [(0, 0, 0), (.5, 0, 0), (.5, 1, 0), (1, 1, 0)]
        self.assertEqual(qcip_math.BLA(position_list_3), -.5)  # first bond is short, so BLA < 0

        with self.assertRaises(Exception):
            qcip_math.BLA([(0, 0, 0), (1, 0, 0)])  # not enough bonds

        with self.assertRaises(Exception):
            qcip_math.BLA([(0, 0, 0), (1, 0, 0), (1, 1, 0)])  # not enough bonds

    def test_permutations(self):
        """Test the permutations function"""

        self.assertEqual(len([a for a in qcip_math.unique_permutations(['a', 'b', 'c'])]), 6)
        self.assertEqual(len([a for a in qcip_math.unique_permutations(['a', 'b', 'b'])]), 3)
        self.assertEqual(len([a for a in qcip_math.unique_permutations(['b', 'b', 'b'])]), 1)

        self.assertEqual(qcip_math.num_of_unique_permutations(['a', 'b', 'c']), 6)
        self.assertEqual(qcip_math.num_of_unique_permutations(['a', 'b', 'b']), 3)
        self.assertEqual(qcip_math.num_of_unique_permutations(['b', 'b', 'b']), 1)

        self.assertEqual([a for a in qcip_math.unique_everseen([1, 2, 3, 1, 4, 1, 2])], [1, 2, 3, 4])

        def tolower_(a):
            return a.lower()

        self.assertEqual(
            [a for a in qcip_math.unique_everseen(['a', 'A', 'b', 'c', 'B'], key=tolower_)], ['a', 'b', 'c'])
