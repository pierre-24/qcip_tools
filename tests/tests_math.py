import unittest
import numpy

from tests import array_almost_equals

from qcip_tools import math as qcip_math


class MathTesCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_vector_manipulation(self):
        """Test the behavior of rodrigues_rotation() normalize(), distance(), angle(), torsion_angle() and BLA()"""

        # test basic stuffs:
        self.assertTrue(array_almost_equals(qcip_math.rodrigues_rotation([1, 0, 0], [0, 0, 1], 90), [0, 1, 0]))
        self.assertTrue(array_almost_equals(qcip_math.rodrigues_rotation([1, 0, 0], [0, 1, 0], 90), [0, 0, -1]))

        self.assertTrue(numpy.array_equal(qcip_math.normalize([-3, 0, 0]), [-1, 0, 0.]))
        self.assertEqual(qcip_math.distance([-1, 0, 0], [.5, 0, 0]), 1.5)
        self.assertEqual(qcip_math.angle([2, 0, 0], [0, 0, 0], [0, 1, 0]), 90.0)
        self.assertEqual(qcip_math.torsion_angle([2, 0, 0], [0, 0, 0], [0, 1, 0], [-2, 1, 0]), -180.0)

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

    def test_bounding_box(self):
        b1 = qcip_math.BoundingBox([0, .5, 0], size=[1, 1, 1])
        b2 = qcip_math.BoundingBox([0, .5, 0], maximum=[1, 1.5, 1])  # same bounding box, but defined different

        self.assertTrue(numpy.array_equal(b2.size, b1.size))

        self.assertTrue(numpy.array_equal(b2.maximum(), [1, 1.5, 1]))
        self.assertTrue(numpy.array_equal(b1.maximum(), b2.maximum()))

        # in and out points
        points_in = [(0, .5, 0), (1, 1.5, 1), (.5, .5, .5), (.5, 1, 0), (.3, .7, 1)]
        for p in points_in:
            self.assertTrue(b1.contains(p))
            self.assertTrue(b2.contains(p))

        points_out = [(0, .5, 2), (2, 1, 1), (.1, .6, 6), (4, 4, 4)]
        for p in points_out:
            self.assertFalse(b1.contains(p))
            self.assertFalse(b2.contains(p))

        # local coordinates:
        self.assertTrue(numpy.array_equal(b1.local_coordinates([0, 0, 0]), [0, -.5, 0]))
        self.assertTrue(numpy.array_equal(b2.local_coordinates([0, 0, 0]), [0, -.5, 0]))

        # extends by changing size
        point_out = (1, 1.75, 1)
        self.assertFalse(b1.contains(point_out))

        b1.update(point_out)
        self.assertTrue(b1.contains(point_out))

        self.assertTrue(numpy.array_equal(b1.origin, [0, .5, 0]))  # origin stays the same
        self.assertTrue(numpy.array_equal(b1.size, [1, 1.25, 1]))  # ... but size change
        self.assertTrue(numpy.array_equal(b1.maximum(), [1, 1.75, 1]))  # ... and maximum has well

        for p in points_in:
            self.assertTrue(b1.contains(p))  # points in are still in

        # extends by changing origin
        point_out = (1, -.5, 1)
        self.assertFalse(b2.contains(point_out))

        b2.update(point_out)
        self.assertTrue(b2.contains(point_out))

        self.assertTrue(numpy.array_equal(b2.origin, [0, -.5, 0]))  # origin has changed
        self.assertTrue(numpy.array_equal(b2.size, [1, 2, 1]))  # ... but size as well !
        self.assertTrue(numpy.array_equal(b2.maximum(), [1, 1.5, 1]))  # ... and maximum has not.

        for p in points_in:
            self.assertTrue(b2.contains(p))  # points in are still in

        # test max[i] < min[i]: those 4 bounding box should have the same origin and size
        bounding_boxes = [
            qcip_math.BoundingBox(origin=[0, 0], size=[1, -2]),
            qcip_math.BoundingBox(origin=[0, 0], maximum=[1, -2]),
            qcip_math.BoundingBox(origin=[0, -2], maximum=[1, 0]),
            qcip_math.BoundingBox(origin=[0, -2], size=[1, 2])
        ]

        true_origin = (0, -2)
        true_size = (1, 2)

        for bb in bounding_boxes:
            self.assertTrue(numpy.array_equal(bb.origin, true_origin))
            self.assertTrue(numpy.array_equal(bb.size, true_size))

        # dimension mismatch in creation
        with self.assertRaises(qcip_math.DimensionAreNotEquals):
            qcip_math.BoundingBox([0, .5, 0], size=[1, 1])
        with self.assertRaises(qcip_math.DimensionAreNotEquals):
            qcip_math.BoundingBox([0, .5], size=[1, 1, 1])
        with self.assertRaises(qcip_math.DimensionAreNotEquals):
            qcip_math.BoundingBox([0, .5, 0], maximum=[1, 1])
        with self.assertRaises(qcip_math.DimensionAreNotEquals):
            qcip_math.BoundingBox([0, .5], maximum=[1, 1, 1])

        # dimension mismatch in local_coordinates()
        with self.assertRaises(qcip_math.DimensionAreNotEquals):
            b1.local_coordinates([2, 2])

        # dimension mismatch in contains()
        with self.assertRaises(qcip_math.DimensionAreNotEquals):
            b1.contains([2, 2])

        # dimension mismatch in update()
        with self.assertRaises(qcip_math.DimensionAreNotEquals):
            b1.update([2, 2])

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
