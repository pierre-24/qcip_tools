import copy
import math

import numpy


class DimensionAreNotEquals(Exception):

    def __init__(self, var1='a', var2='b'):
        super().__init__('{} and {} have not the same dimension'.format(var1, var2))


class BoundingObject:
    """Define a general bounding object

    :param origin: the origin
    :type origin: list|numpy.ndarray
    """

    def __init__(self, origin):
        self.dimension = len(origin)
        self.origin = numpy.array(origin)

    def contains(self, item):
        """Test if point in bounding box

        :param item: coordinates of the point
        :type item: list|numpy.ndarray
        :rtype: bool
        """

        if type(item) not in [numpy.ndarray, list, tuple]:
            raise TypeError(type(item))

        if len(item) != len(self.origin):
            raise DimensionAreNotEquals('coordinates', 'box')

        raise NotImplementedError()

    def __contains__(self, item):
        return self.contains(item)

    def update(self, coo):
        """Update the size of the bounding box to include a given point

        :param coo: coordinates of the point
        :type coo: list|numpy.ndarray
        """

        raise NotImplementedError()

    def local_coordinates(self, coo):
        """Returns coordinates with respect to the origin of the bounding box

        :param coo: coordinates of the point
        :type coo: list|numpy.ndarray
        :rtype: numpy.ndarray
        """

        if len(coo) != self.dimension:
            raise DimensionAreNotEquals('coordinates', 'box')

        return numpy.array([coo[i] - self.origin[i] for i in range(self.dimension)])

    def __add__(self, other):
        """Create a BoundingSet by the addition of two bounding objects"""

        n = BoundingSet()
        n.include(self)
        n.include(other)

        return n

    def __sub__(self, other):
        """Create a BoundingSet by the subtraction of two bounding objects"""

        n = BoundingSet()
        n.include(self)
        n.exclude(other)

        return n


class BoundingSet:
    """Define a set of bounding objects"""

    def __init__(self):

        self.inclusion_set = []
        self.exclusion_set = []
        self.dimension = -1

    def include(self, bounding_object):
        """Add an object to the inclusion set

        :param bounding_object: object
        :type bounding_object: BoundingObject
        """

        if not isinstance(bounding_object, BoundingObject):
            raise TypeError(type(BoundingObject))

        if self.dimension == -1:
            self.dimension = bounding_object.dimension

        if self.dimension != bounding_object.dimension:
            raise DimensionAreNotEquals('set', 'object')

        self.inclusion_set.append(bounding_object)

    def __add__(self, other):

        n = copy.deepcopy(self)
        n.include(other)

        return n

    def __iadd__(self, other):
        self.include(other)
        return self

    def exclude(self, bounding_object):
        """Add an object to the inclusion set

        :param bounding_object: object
        :type bounding_object: BoundingObject
        """

        if not isinstance(bounding_object, BoundingObject):
            raise TypeError(type(BoundingObject))

        if self.dimension == -1:
            self.dimension = bounding_object.dimension

        if self.dimension != bounding_object.dimension:
            raise DimensionAreNotEquals('set', 'object')

        self.exclusion_set.append(bounding_object)

    def __sub__(self, other):

        n = copy.deepcopy(self)
        n.exclude(other)

        return n

    def __isub__(self, other):
        self.exclude(other)
        return self

    def contains(self, item):
        """Test if point in the set

        :param item: coordinates of the point
        :type item: list|numpy.ndarray|tuple
        :rtype: bool
        """

        if type(item) not in [numpy.ndarray, list, tuple]:
            raise TypeError(type(item))

        if len(item) != self.dimension:
            raise DimensionAreNotEquals('coordinates', 'set')

        coordinates_in = False

        for i in self.inclusion_set:
            if item in i:
                coordinates_in = True
                break

        if coordinates_in:
            for i in self.exclusion_set:
                if item in i:
                    coordinates_in = False
                    break

        return coordinates_in

    def __contains__(self, item):
        return self.contains(item)


class AABoundingBox(BoundingObject):
    """Defines a Axis Aligned (AA) bounding box. You must define either ``maximum`` or ``size``.

    .. note::

        Works for any dimension.

    .. note::

        The ``origin`` is corrected so that is the lowest coordinates possible, and ``dimension`` is set
        accordingly, so that it contains only positive numbers.

    :param origin: the origin
    :type origin: list|numpy.ndarray
    :param maximum: the maximum
    :type maximum: list|numpy.ndarray
    :param size: size of bounding box
    :type size: list|numpy.ndarray
    """

    def __init__(self, origin, maximum=None, size=None):
        super().__init__(origin)

        if maximum is not None:
            if self.dimension != len(maximum):
                raise DimensionAreNotEquals('min', 'max')

            self.size = [maximum[a] - origin[a] for a in range(len(origin))]

        elif size is not None:

            if self.dimension != len(size):
                raise DimensionAreNotEquals('min', 'size')

            self.size = size

        else:
            ValueError('you must give either maximum or size')

        for i, val in enumerate(self.size):  # change origin if maximum was actually the minimum
            if val < 0.0:
                self.origin[i] += val
                self.size[i] = -val

    def contains(self, coo):
        """Test if point in bounding box

        :param coo: coordinates of the point
        :type coo: list|numpy.ndarray
        :rtype: bool
        """

        if type(coo) not in [numpy.ndarray, list, tuple]:
            raise TypeError(type(coo))

        if len(coo) != len(self.origin):
            raise DimensionAreNotEquals('coordinates', 'box')

        max_ = self.maximum()

        for i in range(self.dimension):
            if coo[i] < self.origin[i] or coo[i] > max_[i]:
                return False

        return True

    def update(self, coo):
        """Update the size of the bounding box to include a given point

        :param coo: coordinates of the point
        :type coo: list|numpy.ndarray
        """

        if not self.contains(coo):
            maximum = self.maximum()

            for i in range(self.dimension):

                self.origin[i] = coo[i] if coo[i] < self.origin[i] else self.origin[i]
                new_size = maximum[i] - self.origin[i]
                new_size_from_point = coo[i] - self.origin[i]

                self.size[i] = new_size if new_size > new_size_from_point else new_size_from_point

    def maximum(self):
        """Returns the position of the maximum of the bounding box in "global" coordinates

        :rtype: numpy.ndarray
        """

        return numpy.array([self.origin[i] + self.size[i] for i in range(self.dimension)])


class BoundingSphere(BoundingObject):
    """Define a bounding sphere

    .. note::

        Works for any dimension.

    :param origin: the origin
    :type origin: list|numpy.ndarray
    :param radius: the radius
    :type radius: float
    """

    def __init__(self, origin, radius):

        super().__init__(origin)
        self.radius = radius
        self.squared_radius = radius ** 2

    def contains(self, item, tolerance=1e-3):
        """Test if point in bounding sphere

        :param item: coordinates of the point
        :type item: list|numpy.ndarray
        :param tolerance: threshold of tolerance
        :type tolerance: float
        :rtype: bool
        """

        if type(item) not in [numpy.ndarray, list, tuple]:
            raise TypeError(type(item))

        if len(item) != len(self.origin):
            raise DimensionAreNotEquals('coordinates', 'sphere')

        if numpy.sum((numpy.array(item) - self.origin) ** 2) < self.squared_radius + tolerance:
            return True

        return False

    def update(self, coo):
        """Update the radius of the bounding sphere to include a given point

        :param coo: coordinates of the point
        :type coo: list|numpy.ndarray
        """

        if len(coo) != len(self.origin):
            raise DimensionAreNotEquals('coordinates', 'sphere')

        self.squared_radius = numpy.sum((numpy.array(coo) - self.origin) ** 2)
        self.radius = math.sqrt(self.squared_radius)
