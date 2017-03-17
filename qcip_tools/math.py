import math
import numpy
import collections
from functools import reduce
import operator
import itertools


class DimensionAreNotEquals(Exception):

    def __init__(self, var1='a', var2='b'):
        super().__init__('{} and {} have not the same dimension'.format(var1, var2))


class BoundingBox:

    def __init__(self, origin, maximum=None, size=None):
        """Defines a bounding box

        :param origin:
        :param maximum:
        :param size:
        """

        self.dimension = len(origin)
        self.origin = origin

        if maximum is not None:
            if self.dimension != len(maximum):
                raise DimensionAreNotEquals('min', 'max')

            self.size = [maximum[a] - origin[a] for a in range(len(origin))]

        elif size is not None:

            if self.dimension != len(size):
                raise DimensionAreNotEquals('min', 'size')

            self.size = size

        else:
            RuntimeError('you must give either maximum or size')

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

    def local_coordinates(self, coo):
        """Returns coordinates with respect to the origin of the bounding box

        :param coo: coordinates of the point
        :type coo: list|numpy.ndarray
        :rtype: numpy.ndarray
        """

        if len(coo) != self.dimension:
            raise DimensionAreNotEquals('coordinates', 'box')

        return numpy.array([coo[i] - self.origin[i] for i in range(self.dimension)])

    def maximum(self):
        """Returns the position of the maximum of the bounding box in "global" coordinates

        :rtype: numpy.ndarray
        """

        return numpy.array([self.origin[i] + self.size[i] for i in range(self.dimension)])


def rodrigues_rotation(v, rot_axis, angle):
    """
    Rotate vector `v` around rotation axis `rot_axis`, of an angle `angle`
    Note: `v` and `rot_axis` must be defined via the same origin.

    :param v: rotated vector
    :type v: list|numpy.ndarray
    :param rot_axis: axis of rotation (SHOULD BE NORMALIZED !!)
    :type rot_axis: list|numpy.ndarray
    :param angle: angle (in degree)
    :type angle: float
    :return: the rotated vector
    :rtype: numpy.ndarray
    """

    ang_ = numpy.radians(angle)
    ra_ = numpy.array(rot_axis)
    v_ = numpy.array(v)

    v_rot = v_ * math.cos(ang_) + numpy.cross(ra_, v_) * math.sin(ang_) + \
        ra_ * numpy.dot(ra_, v) * (1 - math.cos(ang_))

    return v_rot


def normalize(c):
    """Normalize and return the vector

    :param c: vector
    :type c: numpy.ndarray|list
    :rtype: numpy.ndarray
    """
    n = numpy.linalg.norm(c)
    return c * 1 / n


def distance(c1, c2):
    """
    Get the distance between two position.

    :param c1: coordinate
    :type c1: list
    :param c2: coordinate
    :type c2: list
    :rtype: float
    """
    if c1 is not c2:
        return numpy.linalg.norm(numpy.array(c2) - numpy.array(c1))

    return 0.0


def angle(c1, c2, c3):
    """
    Get the angle (in degree) between three position. Note that `c2` is the center coordinate !!

    :param c1: coordinate
    :type c1: list
    :param c2: coordinate
    :type c2: list
    :param c3: coordinate
    :type c3: list
    :rtype: float
    """

    return numpy.degrees(angle_vector(numpy.array(c2) - numpy.array(c1), numpy.array(c2) - numpy.array(c3)))


def angle_vector(v1, v2):
    """Get the angle (in degree) between two vectors

    :param v1: vector
    :type v1: numpy.ndarray
    :param v2: vector
    :type v2: numpy.ndarray
    :rtype: float
    """

    return math.acos(numpy.dot(v1, v2) / (numpy.linalg.norm(v1) * numpy.linalg.norm(v2)))


def torsion_angle(c1, c2, c3, c4):
    """
    Get the torsion angle (in degree) between four position.

    :param c1: coordinate
    :type c1: list
    :param c2: coordinate
    :type c2: list
    :param c3: coordinate
    :type c3: list
    :param c4: coordinate
    :type c4: list
    :rtype: float
    """

    b1 = numpy.array(c2) - numpy.array(c1)
    b2 = numpy.array(c3) - numpy.array(c2)
    b3 = numpy.array(c4) - numpy.array(c3)
    v1 = normalize(numpy.cross(b1, b2))
    v2 = normalize(numpy.cross(b2, b3))

    return -numpy.degrees(math.atan2(numpy.dot(numpy.cross(v1, normalize(b2)), v2), numpy.dot(v1, v2)))


def BLA(lst):
    """Return the BLA of the atom in the positions givens.

    :param lst: list of position
    :type lst: list
    :rtype: float
    """

    number_of_bonds = len(lst) - 1
    if number_of_bonds < 3:
        raise Exception('not enought bonds')

    distances = []
    previous = lst[0]
    for i in lst[1:]:
        distances.append(distance(previous, i))
        previous = i
    acc = 0
    for i in range(number_of_bonds - 1):
        acc += (distances[i + 1] - distances[i]) * (-1) ** (i - 1)
    return (1 / (number_of_bonds - 1)) * acc


def conjugate(z):
    """Get the conjugate complex of z

    :param z: complex number
    :type z: complex
    :return: the conjugate
    :rtype: complex
    """
    if z.imag != 0:
        return z.real - z.imag * 1j
    else:
        return z


def prod(iterable):
    """Same as sum(), but multiplying

    :param iterable: an iterable
    :type iterable: iterable
    :return: multiplication
    :rtype: float
    """
    return reduce(operator.mul, iterable, 1)


def unique_permutations(elements):
    """Like itertools.permutation(), but yield UNIQUE elements.
    Iterative. May be not as fast as possible.

    Credit to http://stackoverflow.com/a/30558049.

    :param elements: a set of elements
    :type elements: iterable
    """

    if len(elements) == 1:
        yield (elements[0],)
    else:
        unique_elements = set(elements)
        for first_element in unique_elements:
            remaining_elements = list(elements)
            remaining_elements.remove(first_element)
            for sub_permutation in unique_permutations(remaining_elements):
                yield (first_element,) + sub_permutation


def num_of_unique_permutations(elements):
    """Get the number of unique elements. Compute the
    `multinomial coefficients <https://en.wikipedia.org/wiki/Multinomial_theorem#Multinomial_coefficients>`_:

    .. math::

        \\prod_i^m \\left(\\begin{matrix}\\sum_j^i k_j\\\\k_i\\end{matrix}\\right) = \\frac{n!}{\\prod_i k_i!}.

    where :math:`n` is the number of elements in the set and :math:`k_i` is the number of :math:`i` in the set (their
    multiplicities).

    This is equivalent to ``len(unique_permutation(elements))`` (but faster).

    :param elements: a set
    :type elements: iterable
    :return: length
    :rtype: int
    """

    n = len(elements)
    each = collections.Counter(elements)
    return int(math.factorial(n) / prod(math.factorial(each[i]) for i in each))


def unique_everseen(iterable, key=None):
    """List unique elements, preserving order. Remember all elements ever seen.

    Credit goes to the "recipes" in https://docs.python.org/dev/library/itertools.html.

    :param iterable: any iterable
    :type iterable: iterable
    :param key: apply a function to each element before checking if seen
    :type key: function
    """

    seen = set()
    seen_add = seen.add
    if key is None:
        for element in itertools.filterfalse(seen.__contains__, iterable):
            seen_add(element)
            yield element
    else:
        for element in iterable:
            k = key(element)
            if k not in seen:
                seen_add(k)
                yield element
