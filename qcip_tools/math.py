import math
import numpy


def rodrigues_rotation(rot_axis, angle, v):
    """
    Rotate vector `v` around rotation axis `rot_axis`, of an angle `angle`
    Note: `v` and `rot_axis` must be defined via the same origin.

    :param rot_axis: axis of rotation (SHOULD BE NORMALIZED !!)
    :param angle: angle (in radian)
    :param v: rotated vector
    :return: the rotated vector
    """
    v_rot = v * math.cos(angle) + numpy.cross(rot_axis, v) * math.sin(angle) + \
        rot_axis * numpy.dot(rot_axis, v) * (1 - math.cos(angle))
    return v_rot


def normalize(c):
    """Normalize and return the vector

    :param c: vector
    :type c: numpy.ndarray
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
    Get the angle between three position. Note that `c2` is the center coordinate !!

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
    """Get the angle between two vectors

    :param v1: vector
    :type v1: numpy.ndarray
    :param v2: vector
    :type v2: numpy.ndarray
    :rtype: float
    """

    return math.acos(numpy.dot(v1, v2) / (numpy.linalg.norm(v1) * numpy.linalg.norm(v2)))


def torsion_angle(c1, c2, c3, c4):
    """
    Get the torsion angle between four position.

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
        return 0
    distances = []
    previous = lst[0]
    for i in lst[1:]:
        distances.append(distance(previous, i))
        previous = i
    acc = 0
    for i in range(number_of_bonds - 1):
        acc += (distances[i + 1] - distances[i]) * (-1) ** (i - 1)
    return (1 / (number_of_bonds - 1)) * acc
