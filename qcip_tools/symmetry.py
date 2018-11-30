"""
Symmetry handling.

See https://en.wikipedia.org/wiki/Group_(mathematics)

Base implementation inspired from https://github.com/naftaliharris/Abstract-Algebra/

See https://physics.stackexchange.com/questions/351372/generate-all-elements-of-a-point-group-from-generating-set
and https://pdfs.semanticscholar.org/4edd/1ac673ea4cab251bb0b64bf0f5b65f18cc9d.pdf
and http://www.euclideanspace.com/maths/geometry/affine/reflection/quaternion/index.htm
"""

import numpy
import math
from transforms3d import quaternions


class Set(set):
    """Define a (finite) set of (unique) elements
    """

    def __mul__(self, other):
        """Cartesian product:

        .. math::

            A\\times B = \\{(x,y) | x \\in A, y \\in B\\}

        :rtype: Set
        """

        if not isinstance(other, Set):
            raise TypeError('other must be a set')

        return Set((x, y) for x in self for y in other)


class BinaryOperation:
    """Define a binary operation :math:`f` such that,

    .. math::

        f:S\\times S \rightarrow S.

    Thus, we only define a binary operation trough its codomain.

    See https://en.wikipedia.org/wiki/Binary_operation
    """

    def __init__(self, codomain, func):
        self.function = func
        self.codomain = codomain

    def __call__(self, e):
        """Call the function on an element of the domain

        :return: an element of the codomain
        """
        if e[0] not in self.codomain or e[1] not in self.codomain:
            raise TypeError('not an element of the domain: {}'.format(e))

        return self.function(e)

    def image(self):
        """Get the image set (which is part of the codomain)

        :rtype: Set
        """

        return Set(self(e) for e in self.codomain * self.codomain)

    def check_closure(self):
        """Check that the relation gives elements that are in the codomain

        :rtype: bool
        """

        return self.image() <= self.codomain

    def check_surjectivity(self):
        """Check surjectivity (image equals codomain)

        :rtype: bool
        """

        return self.image() == self.codomain

    def check_associativity(self):
        """Check if the binary operation is associative for all the element of the domain

        :rtype: bool
        """

        return all(
            self((a, self((b, c)))) == self((self((a, b)), c))
            for a in self.codomain for b in self.codomain for c in self.codomain)


class GroupElement:
    """Syntaxic sugar
    """

    def __init__(self, element, group):
        if not isinstance(group, Group):
            raise TypeError('Group must be group')
        if element not in group.binary_operation.codomain:
            raise TypeError('{} is not in group element'.format(element))

        self.element = element
        self.group = group

    def __repr__(self):
        return repr(self.element)

    def __mul__(self, other):
        e = other
        if isinstance(other, GroupElement):
            e = other.element
        if e not in self.group:
            raise ValueError('{} is not in group domain')

        return GroupElement(self.group.binary_operation((self.element, e)), self.group)

    def __rmul__(self, other):
        e = other
        if isinstance(other, GroupElement):
            e = other.element

        return GroupElement(self.group.binary_operation((e, self.element)), self.group)

    def __eq__(self, other):
        e = other
        if isinstance(other, GroupElement):
            e = other.element

        return self.element == e

    def __pow__(self, power, modulo=None):
        if not isinstance(power, int):
            raise TypeError('power must be an integer')

        if power == 0:
            return self.group.e
        elif power == 1:
            return self
        elif power < 0:
            return self.group.inverse(self.element) ** -power
        elif power % 2 == 1:
            return self * (self ** (power - 1))
        else:
            return (self * self) ** int(power / 2)

    def __hash__(self):
        return hash(self.element)


class Group:
    """A set :math:`G`, together with an operation :math:`\\star` (group law),
    form a group :math:`(G, \\star)` if it satisfies 4 requirements:

    - Closure (not checked here) ;
    - Associativity (not checked here) ;
    - (unique) Identity ;
    - Inverse.

    A group may also be Abelian.
    """

    def __init__(self, binary_operation):
        if not isinstance(binary_operation, BinaryOperation):
            raise TypeError('binary_operation must be a BinaryOperation')

        self.binary_operation = binary_operation
        self.G = [GroupElement(e, self) for e in binary_operation.codomain]

        # get identity
        self.e = None
        for i in self.G:
            if all(i * a == a for a in self.G):
                self.e = i
                break

        if self.e is None:
            raise RuntimeError('G does not contain an identity element')

        # get inverses (and check if Abelian in the meantime)
        self.inverses = {}
        invs = []
        self.abelian = True
        for i in self.G:
            inverse_found = False
            for j in self.G:
                if i * j == self.e:
                    if j in invs:
                        raise RuntimeError(
                            '{} is already the inverse of an element, so cannot be the one of {} as well'.format(j, i))

                    self.inverses[i] = j
                    invs.append(j)
                    inverse_found = True
                    if j in self.inverses and self.inverses[j] != i:
                        self.abelian = False
                    break
            if not inverse_found:
                raise RuntimeError('{} does not have an inverse!'.format(i))

    def identity(self):
        """Get identity

        :return: an element of G
        """

        return self.e

    def inverse(self, element):
        """Return an inverse of a given element of G

        :param element: the element
        :return: another element of G
        """

        return self.inverses[element]

    def __contains__(self, item):
        """is an element part of the group

        :param item: the element
        :rtype: bool
        """
        if isinstance(item, GroupElement):
            return item in self.G
        else:
            return item in self.binary_operation.codomain


class SymElement:
    """Define a symmetry element (fully based on quaternion).

    According `Mamone et al. <https://www.mdpi.com/2073-8994/2/3/1423/pdf>`_,

    .. math::

        \\hat{S}^k_n = \\hat{R}(\\pi + \\frac{2\\pi}{n}),

    where :math:`\\hat{R}` is the normal rotation operator.
    To use this operator on a point, the opposite of the sandwich product must be taken, so that the conjugate of the
    quaternion that represent the point is actually used.

    TODO: __equal__, __pow__, __hash__ (for usage with Group)
    """

    def __init__(self, q, improper=False):
        self.q = q
        self.improper = improper

    def __mul__(self, other):
        if not isinstance(other, SymElement):
            raise TypeError('other must be SymElement')

        q = quaternions.qmult(other.q, self.q)
        return SymElement(q, (self.improper and not other.improper) or (other.improper and not self.improper))

    def apply(self, pin):
        """Apply the operation on a point

        :param pin: the point
        :type pin: numpy.ndarray
        :rtype: numpy.ndarray
        """

        qp = numpy.array([0., *pin])
        if self.improper:
            qp = quaternions.qconjugate(qp)

        qe = quaternions.qmult(self.q, quaternions.qmult(qp, quaternions.qconjugate(self.q)))
        return qe[1:]

    @classmethod
    def from_axangle(cls, axis, angle=2 * numpy.pi, improper=False, in_degree=False):
        """Create from axe and angle

        :param axis: the rotation axe
        :type axis: numpy.ndarray
        :param angle: the rotation angle
        :type angle: float
        :param improper: is it an improper rotation ?
        :type improper: bool
        :param in_degree: is the angle in degree
        :type in_degree: bool
        :rtype: SymElement
        """

        if in_degree:
            angle = numpy.radians(angle)
        if improper:
            angle += numpy.pi

        axis /= numpy.linalg.norm(axis)

        q = numpy.array([math.cos(angle / 2), *(math.sin(angle / 2) * axis)])
        return cls(q, improper)

    @classmethod
    def E(cls):
        """Generate the identity

        :rtype: SymElement
        """
        return cls(quaternions.qeye())

    @classmethod
    def C(cls, n, k=1, axis=numpy.array([0, 0, 1.])):
        """Generate the representation of a proper rotation (parallel to ẑ)

        :param n: order of the axis
        :type n: int
        :param k: number of step
        :type k: int
        :param axis: to which axis this rotation axis should be parallel ?
        :type axis: numpy.ndarray
        :rtype: SymElement
        """

        if n < 1:
            raise ValueError('n<1')

        k %= n

        if k == 0:
            return cls.E()

        g = math.gcd(n, k)
        if g != 1:
            n = int(n / g)
            k = int(k / g)

        return cls.from_axangle(axis, 2 * numpy.pi / n * k)

    @classmethod
    def S(cls, n, k=1, axis=numpy.array([0, 0, 1.])):
        """Generate the representation of a proper rotation (parallel to ẑ)

        :param n: order of the axis
        :type n: int
        :param k: number of step
        :type k: int
        :param axis: to which axis this rotation axis should be parallel ?
        :type axis: numpy.ndarray
        :rtype: SymElement
        """
        if n < 1:
            raise ValueError('n<1')

        if n % 2 == 1:
            k %= 2 * n
        else:
            k %= n

        if k == 0:
            return cls.E()

        if n % 2 == 0 and k % 2 == 0:
            return cls.C(k, n, axis=axis)
        else:
            return cls.from_axangle(axis, 2 * numpy.pi / n * k, improper=True)

    @classmethod
    def i(cls):
        """Generate an inversion center

        :rtype: SymElement
        """
        return cls.S(2)

    @classmethod
    def sigma(cls, axis=numpy.array([0, 0, 1.])):
        """Generate a symmetry plane

        :param axis: normal axis to the plane
        :type axis: numpy.ndarray
        :rtype: SymElement
        """
        return cls.S(1, axis=axis)
