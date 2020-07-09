import math
import numpy
import collections
from functools import reduce
import operator
import itertools


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
    if number_of_bonds < 2:
        raise Exception('not enough bonds')

    pos = numpy.vstack(lst)
    distances = numpy.linalg.norm(pos[:-1] - pos[1:], axis=1)
    acc = 0
    for i in range(number_of_bonds - 1):
        acc += (distances[i + 1] - distances[i]) * (-1) ** i
    return (1 / (number_of_bonds - 1)) * acc


def conjugate(z):
    """Get the conjugate complex of z

    :param z: complex number
    :type z: complex
    :return: the conjugate
    :rtype: complex
    """
    if type(z) is not complex:
        raise TypeError(z)

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


def closest_fraction(x, max_denominator=2 * 3 * 4 * 5 * 7 * 9 * 11 * 13):
    """A fairly simple algorithm to get the fraction out of a floating number.
    Try to get as close as possible to ``x``.

    Notice that the maximum denominator is chosen to be the multiplication of the primes (and their even power)

    :param x: the float
    :type x: float
    :param max_denominator: the maximum denominator
    :type max_denominator: int
    :rtype: tuple(int, int)
    """
    e = int(x * max_denominator)
    u = math.fabs(e / max_denominator - x)
    if math.fabs((e + 1) / max_denominator - x) <= u:
        e += 1
    elif math.fabs((e - 1) / max_denominator - x) <= u:
        e -= 1
    g = math.gcd(e, max_denominator)
    return e // g, max_denominator // g


def closest_int(x):
    """Compute the closest integer

    :param x: a float
    :type x: float
    :rtype: int
    """
    c = math.ceil(x)
    f = math.floor(x)
    if math.fabs(x - c) < math.fabs(x - f):
        return int(c)
    else:
        return int(f)
