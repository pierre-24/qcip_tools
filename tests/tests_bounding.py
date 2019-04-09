import numpy

from qcip_tools import bounding
from tests import QcipToolsTestCase


class BoundingTestCase(QcipToolsTestCase):

    def setUp(self):
        pass

    def test_aa_bounding_box(self):
        b1 = bounding.AABoundingBox([0, .5, 0], size=[1, 1, 1])
        b2 = bounding.AABoundingBox([0, .5, 0], maximum=[1, 1.5, 1])  # same bounding box, but defined different

        self.assertTrue(numpy.array_equal(b2.size, b1.size))

        self.assertTrue(numpy.array_equal(b2.maximum(), [1, 1.5, 1]))
        self.assertTrue(numpy.array_equal(b1.maximum(), b2.maximum()))

        # in and out points
        points_in = [(0, .5, 0), (1, 1.5, 1), (.5, .5, .5), (.5, 1, 0), (.3, .7, 1)]
        for p in points_in:
            self.assertTrue(p in b1)
            self.assertTrue(p in b2)

        points_out = [(0, .5, 2), (2, 1, 1), (.1, .6, 6), (4, 4, 4)]
        for p in points_out:
            self.assertFalse(p in b1)
            self.assertFalse(p in b2)

        # local coordinates:
        self.assertTrue(numpy.array_equal(b1.local_coordinates([0, 0, 0]), [0, -.5, 0]))
        self.assertTrue(numpy.array_equal(b2.local_coordinates([0, 0, 0]), [0, -.5, 0]))

        # extends by changing size
        point_originally_out = [1, 1.75, 1]
        self.assertFalse(b1.contains(point_originally_out))

        b1.update(point_originally_out)
        self.assertTrue(b1.contains(point_originally_out))

        self.assertTrue(numpy.array_equal(b1.origin, [0, .5, 0]))  # origin stays the same
        self.assertTrue(numpy.array_equal(b1.size, [1, 1.25, 1]))  # ... but size change
        self.assertTrue(numpy.array_equal(b1.maximum(), [1, 1.75, 1]))  # ... and maximum has well

        for p in points_in:
            self.assertTrue(b1.contains(p))  # points in are still in

        # extends by changing origin
        point_originally_out = [1, -.5, 1]
        self.assertFalse(b2.contains(point_originally_out))

        b2.update(point_originally_out)
        self.assertTrue(b2.contains(point_originally_out))

        self.assertTrue(numpy.array_equal(b2.origin, [0, -.5, 0]))  # origin has changed
        self.assertTrue(numpy.array_equal(b2.size, [1, 2, 1]))  # ... but size as well !
        self.assertTrue(numpy.array_equal(b2.maximum(), [1, 1.5, 1]))  # ... and maximum has not.

        for p in points_in:
            self.assertTrue(b2.contains(p))  # points in are still in

        # test max[i] < min[i]: those 4 bounding box should have the same origin and size
        bounding_boxes = [
            bounding.AABoundingBox(origin=[0, 0], size=[1, -2]),
            bounding.AABoundingBox(origin=[0, 0], maximum=[1, -2]),
            bounding.AABoundingBox(origin=[0, -2], maximum=[1, 0]),
            bounding.AABoundingBox(origin=[0, -2], size=[1, 2])
        ]

        true_origin = (0, -2)
        true_size = (1, 2)

        for bb in bounding_boxes:
            self.assertTrue(numpy.array_equal(bb.origin, true_origin))
            self.assertTrue(numpy.array_equal(bb.size, true_size))

        # dimension mismatch in creation
        with self.assertRaises(bounding.DimensionAreNotEquals):
            bounding.AABoundingBox([0, .5, 0], size=[1, 1])
        with self.assertRaises(bounding.DimensionAreNotEquals):
            bounding.AABoundingBox([0, .5], size=[1, 1, 1])
        with self.assertRaises(bounding.DimensionAreNotEquals):
            bounding.AABoundingBox([0, .5, 0], maximum=[1, 1])
        with self.assertRaises(bounding.DimensionAreNotEquals):
            bounding.AABoundingBox([0, .5], maximum=[1, 1, 1])

        # dimension mismatch in local_coordinates()
        with self.assertRaises(bounding.DimensionAreNotEquals):
            b1.local_coordinates([2, 2])

        # dimension mismatch in contains()
        with self.assertRaises(bounding.DimensionAreNotEquals):
            b1.contains([2, 2])

        # dimension mismatch in update()
        with self.assertRaises(bounding.DimensionAreNotEquals):
            b1.update([2, 2])

    def test_bounding_sphere(self):
        """Test bounding sphere"""

        b = bounding.BoundingSphere([0, .5, 0], radius=2)

        # test in and out
        points_in = [
            b.origin,
            b.origin - [2.0, 0, 0],
            b.origin + [1.0, 1.0, 1.0],
            b.origin + [1.0, -1.0, .5]
        ]

        for p in points_in:
            self.assertTrue(p in b)

        points_out = [
            b.origin - [2.1, 0, 0],
            b.origin + [1.0, 2, 1.0],
            b.origin + [1.0, -1.0, 1.5],
            b.origin + [5, 0, 0]
        ]

        for p in points_out:
            self.assertTrue(p not in b)

        # test update:
        point_originally_out = b.origin + [2.1, 0, 0]
        self.assertTrue(point_originally_out not in b)

        b.update(point_originally_out)
        self.assertTrue(point_originally_out in b)

        self.assertAlmostEqual(b.radius, 2.1, places=1)
        self.assertArraysAlmostEqual(b.origin, [0, .5, 0])  # origin has not changed !
        self.assertArraysAlmostEqual(b.local_coordinates(point_originally_out), [2.1, 0, 0], places=1)

    def test_bounding_set(self):
        """Test the set of bounding objects"""

        bs = bounding.BoundingSphere([.5, 0, 0], radius=1)
        point_in_sphere = bs.origin + [.25, .1, .1]
        point_in_both = [.0, .0, .0]
        self.assertTrue(point_in_sphere in bs)
        self.assertTrue(point_in_both in bs)

        bc = bounding.AABoundingBox([-.5, -.5, -.5], size=[.5, .5, .5])
        point_in_cube = bc.origin + [.1, .1, .25]
        self.assertTrue(point_in_cube in bc)
        self.assertTrue(point_in_both in bc)
        self.assertTrue(point_in_sphere not in bc)
        self.assertTrue(point_in_cube not in bs)

        # create set
        s = bs + bc
        self.assertTrue(all([point_in_sphere in s, point_in_cube in s, point_in_both in s]))

        # create set by excluding sphere:
        s = bc - bs
        self.assertTrue(all([point_in_sphere not in s, point_in_both not in s, point_in_cube in s]))

        # create set by excluding cube:
        s = bs - bc
        self.assertTrue(all([point_in_sphere in s, point_in_both not in s, point_in_cube not in s]))
