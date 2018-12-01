"""
Symmetry handling.

See https://en.wikipedia.org/wiki/Group_(mathematics)

Base implementation inspired from https://github.com/naftaliharris/Abstract-Algebra/

See https://physics.stackexchange.com/questions/351372/generate-all-elements-of-a-point-group-from-generating-set
and https://pdfs.semanticscholar.org/4edd/1ac673ea4cab251bb0b64bf0f5b65f18cc9d.pdf
"""

import numpy
import math
from enum import Enum

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
    """
    Define a binary operation :math:`f` such that,

    .. math::

        f:S\\times S \\rightarrow S.

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

    @classmethod
    def generate(cls, generators, func):
        """Create a complete set out of the generators.

        Algorithm from https://physics.stackexchange.com/a/351400.

        :param generators; the generators of the set
        :type generators: list
        :param func: function
        :type func: function
        :rtype: BinaryOperation
        """

        g = g1 = generators[0]
        L = [g]
        while True:
            g = func((g, g1))
            if g in L:
                break
            L.append(g)

        for i in range(1, len(generators)):
            C = [g1]
            L1 = L.copy()
            more = True
            while more:
                more = False
                for g in C:
                    for s in generators[:i + 1]:
                        sg = s * g
                        if sg not in L:
                            C.append(sg)
                            L.extend([func((sg, t)) for t in L1])
                            more = True

        return cls(Set(L), func)


class GroupElement:
    """Syntaxic sugar
    """

    def __init__(self, element, group):
        if not isinstance(group, Group):
            raise TypeError('Group must be group')
        if element not in group.binary_operation.codomain:
            raise TypeError('{} is not in group elements'.format(element))

        self.element = element
        self.group = group

    def __str__(self):
        return str(self.element)

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
        """Compute the power of a group element

        Only the multiplication is required for this to work.

        :param power: the power
        :type power: int
        :param modulo: (**modulo is not used**)
        :rtype: GroupElement
        """
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

    def __iter__(self):
        yield self.e
        for i in self.G:
            if i != self.e:
                yield i


def closest_fraction(x, max_denominator=2 * 3 * 4 * 5 * 7 * 9 * 11 * 13):
    """A fairly simple algorithm to get the fraction out of a floating number.
    Try to get as close as possible to ``x``.

    Notice that the maximum denominator is chosen to be the multiplication of the primes (and their odd power)

    :param x: the float
    :type x: float
    :param max_denominator: the maximum denominator
    :type max_denominator: int
    :rtype: tuple(int, int)
    """
    e = int(x * max_denominator)
    if abs((e + 1) / max_denominator - x) <= abs(e / max_denominator - x):
        e += 1
    g = math.gcd(e, max_denominator)
    return e // g, max_denominator // g


class Operation:
    """Define a symmetry element (fully based on quaternion).

    According `Mamone et al. <https://www.mdpi.com/2073-8994/2/3/1423/pdf>`_,

    .. math::

        \\hat{S}^k_n = \\hat{R}(\\pi + \\frac{2\\pi}{n}),

    where :math:`\\hat{R}` is the normal rotation operator.
    To use this operator on a point, the opposite of the sandwich product must be taken, so that the conjugate of the
    quaternion that represent the point is actually used.
    """

    def __init__(self, q, improper=False, description=None):
        self.q = q
        self.improper = improper
        self.description = description
        self._hash = None

    def __mul__(self, other):
        if not isinstance(other, Operation):
            raise TypeError('other must be Operation')

        q = quaternions.qmult(other.q, self.q)
        return Operation(q, (self.improper and not other.improper) or (other.improper and not self.improper))

    def __eq__(self, other):
        """
        Check equality. Note that a quaternion or its multiplication by -1 gives the same transformation
        (but this is checked by ``transform3d``).

        :param other: the other symmetry element
        :type other: Operation
        :rtype: bool
        """
        if not isinstance(other, Operation):
            return False

        return quaternions.nearly_equivalent(self.q, other.q) and self.improper == other.improper

    def __hash__(self):
        """Hash based on the quaternion

        1. Set the quaternion to be positive ;
        2. Take the fraction corresponding to the 4 components of the quaternion.
        3. Add the improper variable.
        4. Hash the corresponding tuple.
        """

        def _should_negate(q, thresh=.00001):
            for i in range(4):
                if math.fabs(q[i]) < thresh:
                    continue

                return q[i] < .0

        if not self._hash:
            q = self.q
            if _should_negate(q):
                q = -q

            self._hash = hash((*(closest_fraction(x) for x in q), self.improper))

        return self._hash

    def __str__(self):
        return str(self.get_description())

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

    def get_description(self, force=False):
        """Try to recognize the transformation:

        1. Extract axis and angle from quaternion (using ``quat2axangle()``) ;
        2. Treat angle (:math:`\\theta`):

           - if the axis have a negative z, take the opposite and set :math:`\\theta = 2\\pi-\\theta`,
           - if this is an improper rotation, set :math:`\\theta=\\theta-\\pi`,
           - Cast between 0 and :math:`4\\pi`,
           - divide by :math:`2\\pi` (so that the result is between 0 and 2).
        3. Try to extract a simple fraction (thus :math:`k` and :math:`n`) out of this number ;
        4. Further reduce k to be smaller than :math:`n` if proper rotation or odd :math:`n` ;
        5. Recognize the specific :math:`n` cases (1 or 2).

        :param force: use ``self.description`` if ``False``, or try a new recognition otherwise
        :type force: bool
        :rtype: OperationDescription
        """

        def _should_negate(axis, thresh=.001):
            p = sorted(list(axis))
            if math.fabs(p[0] - p[1]) < thresh:
                return False
            else:
                return p[0] < .0

        if self.description is None or force:
            axis, angle = quaternions.quat2axangle(self.q)

            # treatment on angle:
            angle /= 2 * numpy.pi

            if _should_negate(axis):  # return axis so that it always have a positive z
                axis = -axis
                angle = 1 - angle

            if self.improper:
                angle -= .5

            if self.improper:
                angle %= 2
            else:
                angle %= 1

            k, n = closest_fraction(angle)

            if not self.improper or n % 2 == 0:
                k %= n

            # find symbol
            symbol = OperationType.improper_rotation if self.improper else OperationType.proper_rotation

            if n == 1:
                if not self.improper:
                    symbol = OperationType.identity
                else:
                    symbol = OperationType.reflexion_plane

            if n == 2 and self.improper:
                symbol = OperationType.inversion

            self.description = OperationDescription(symbol, axis, n, k)

        return self.description

    @classmethod
    def from_axangle(cls, axis, angle=2 * numpy.pi, improper=False, in_degree=False, description=None):
        """Create from axe and angle

        :param axis: the rotation axe
        :type axis: numpy.ndarray
        :param angle: the rotation angle
        :type angle: float
        :param improper: is it an improper rotation ?
        :type improper: bool
        :param in_degree: is the angle in degree
        :param description: description of the operation
        :type description: OperationDescription
        :type in_degree: bool
        :rtype: Operation
        """

        if in_degree:
            angle = numpy.radians(angle)
        if improper:
            angle += numpy.pi

        axis /= numpy.linalg.norm(axis)

        q = numpy.array([math.cos(angle / 2), *(math.sin(angle / 2) * axis)])
        return cls(q, improper, description=description)

    @classmethod
    def E(cls):
        """Generate the identity

        :rtype: Operation
        """
        return cls(quaternions.qeye(), description=OperationDescription(OperationType.identity))

    @classmethod
    def C(cls, n, k=1, axis=numpy.array([0, 0, 1.])):
        """Generate the representation of a proper rotation

        :param n: order of the axis
        :type n: int
        :param k: number of step
        :type k: int
        :param axis: to which axis this rotation axis should be parallel ?
        :type axis: numpy.ndarray
        :rtype: Operation
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

        return cls.from_axangle(
            axis, 2 * numpy.pi / n * k, description=OperationDescription(OperationType.proper_rotation, axis, n, k))

    @classmethod
    def S(cls, n, k=1, axis=numpy.array([0, 0, 1.])):
        """Generate the representation of a improper rotation

        :param n: order of the axis
        :type n: int
        :param k: number of step
        :type k: int
        :param axis: to which axis this rotation axis should be parallel ?
        :type axis: numpy.ndarray
        :rtype: Operation
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
            symbol = OperationType.improper_rotation
            if n == 1:
                symbol = OperationType.reflexion_plane
            elif n == 2:
                symbol = symbol.inversion

            return cls.from_axangle(
                axis, 2 * numpy.pi / n * k, improper=True, description=OperationDescription(symbol, axis, n, k))

    @classmethod
    def i(cls):
        """Generate an inversion center

        :rtype: Operation
        """
        return cls.S(2)

    @classmethod
    def sigma(cls, axis=numpy.array([0, 0, 1.])):
        """Generate a symmetry plane

        :param axis: normal axis to the plane
        :type axis: numpy.ndarray
        :rtype: Operation
        """
        return cls.S(1, axis=axis)


class OperationType(str, Enum):
    """Symbol of the symmetry element, in `Schoenflies notation <https://en.wikipedia.org/wiki/Schoenflies_notation>`_.

    """
    proper_rotation = 'C'
    improper_rotation = 'S'
    identity = 'E'
    reflexion_plane = 'σ'
    inversion = 'i'


class OperationDescription:
    """Store the description of a symmetry operation
    """

    def __init__(self, symbol, axis=numpy.zeros(3), n=1, k=1):
        if symbol not in OperationType:
            raise ValueError('not an allowed symbol')

        self.symbol = symbol
        self.axis = axis
        self.textual_axis = OperationDescription.find_textual_axis(axis)
        self.n = n
        self.k = k

    def __eq__(self, other):
        if not isinstance(other, OperationDescription):
            return False

        # TODO: being less strict on axis, k and n if E, sigma or i
        return self.symbol == other.symbol and \
            self.n == other.n and self.k == other.k and numpy.allclose(self.axis, other.axis)

    def __str__(self):
        axis = ','.join('{:.3f}'.format(e) for e in self.axis)
        if self.symbol in [OperationType.identity, OperationType.inversion]:
            return self.symbol.value
        if self.symbol == OperationType.reflexion_plane:
            return '{}({})'.format(self.symbol.value, self.textual_axis if self.textual_axis else axis)
        else:
            return '{}({},{},{})'.format(
                self.symbol.value, self.n, self.k, self.textual_axis if self.textual_axis else axis)

    @staticmethod
    def find_textual_axis(axis):
        """Try to find the main axis (does not try a lot, though).

        :param axis: the axis
        :type axis: numpy.ndarray
        :rtype: str|None
        """

        if numpy.allclose(axis, [0, 0, 1]):
            return 'z'
        if numpy.allclose(axis, [0, 1, 0]):
            return 'y'
        if numpy.allclose(axis, [1, 0, 0]):
            return 'x'


class PointGroupType(str, Enum):
    cyclic = 'C_n'
    pyramidal = 'C_nv'
    prismatic = 'C_nh'
    gyro_n_gonal = 'S_n'
    dihedral = 'D_n'
    ortho_n_gonal = 'D_nh'
    antiprismatic = 'D_nd'
    tetrahedral_chiral = 'T'
    pyritohedral = 'T_h'
    tetrahedral_achiral = 'T_d'
    octahedral_chiral = 'O'
    octahedral_achiral = 'O_h'
    icosahedral_chiral = 'I'
    icosahedral_achiral = 'I_h'


class PointGroupDescription:
    def __init__(self, symbol, n=0):
        if symbol not in PointGroupType:
            raise ValueError('unrecognized group type')

        if symbol in [
            PointGroupType.cyclic, PointGroupType.pyramidal, PointGroupType.prismatic,
            PointGroupType.gyro_n_gonal, PointGroupType.dihedral, PointGroupType.ortho_n_gonal,
                PointGroupType.antiprismatic]:

            if n < 1:
                    raise ValueError('order must be equal or superior to 1')
        else:
            if n > 0:
                raise ValueError('order does not matter for point group {}'.format(symbol.value))

        self.symbol = symbol
        self.order = n

    def __str__(self):
        s = self.symbol.value
        if 'n' in s:
            s = s.replace('n', str(self.order))
        return s


class PointGroup(Group):

    def __init__(self, func, description=None):
        super().__init__(func)
        self.description = description

    def __str__(self):
        return '{} <{}>'.format(
            str(self.description) if self.description is not None else 'Point group',
            ','.join(str(i) for i in self)
        )

    @staticmethod
    def product(e):
        return e[0] * e[1]

    @classmethod
    def generate(cls, generators, description=None):
        return cls(BinaryOperation.generate(generators, PointGroup.product), description)

    @classmethod
    def C_n(cls, n=1):
        """Create a cyclic group of order n

        :param n; order of the group
        :type n: int
        :rtype: PointGroup
        """
        return cls.generate([Operation.C(n)], description=PointGroupDescription(PointGroupType.cyclic, n))

    @classmethod
    def C_nv(cls, n=1):
        """Create a pyramidal group of order n

        :param n; order of the group
        :type n: int
        :rtype: PointGroup
        """
        return cls.generate(
            [Operation.C(n), Operation.sigma(axis=numpy.array([1., 0, 0]))],
            description=PointGroupDescription(PointGroupType.pyramidal, n))

    @classmethod
    def S_n(cls, n=1):
        """Create a prismatic (odd n) or gyro-n-gonal (even n) point group of order n

        :param n; order of the group
        :type n: int
        :rtype: PointGroup
        """

        if n % 2 == 0:
            type_ = PointGroupType.gyro_n_gonal
        else:
            type_ = PointGroupType.prismatic

        return cls.generate(
            [Operation.S(n)],
            description=PointGroupDescription(type_, n))

    @classmethod
    def C_nh(cls, n=1):
        """Create a prismatic point group of order n

        :param n; order of the group
        :type n: int
        :rtype: PointGroup
        """

        return cls.generate(
            [Operation.C(n), Operation.sigma()],
            description=PointGroupDescription(PointGroupType.prismatic, n))
