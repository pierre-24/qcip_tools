"""
Symmetry handling.

See https://en.wikipedia.org/wiki/Group_(mathematics)

Base implementation inspired by https://github.com/naftaliharris/Abstract-Algebra/

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

        # print([str(l) for l in L])
        return cls(Set(L), func)


class NotInGroup(Exception):
    pass


class GroupElement:
    """Abstract group element: syntaxic sugar
    """

    def __init__(self, element, group):
        if not isinstance(group, Group):
            raise TypeError('Group must be group')
        if element not in group.binary_operation.codomain:
            raise NotInGroup(element)

        self.element = element
        self.group = group
        self._hash = None

    def __repr__(self):
        return str(self.element)

    def __mul__(self, other):
        e = other
        if isinstance(other, GroupElement):
            e = other.element
        if e not in self.group:
            raise NotInGroup(e)

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

        if not self._hash:
            self._hash = hash(self.element)

        return self._hash


class Group:
    """Abstract group structure.

    A set :math:`G`, together with an operation :math:`\\star` (group law),
    form a group :math:`(G, \\star)` if it satisfies 4 requirements:

    - Closure (not checked here) ;
    - Associativity (not checked here) ;
    - (unique) Identity ;
    - (unique) Inverse.
    """

    def __init__(self, binary_operation):
        if not isinstance(binary_operation, BinaryOperation):
            raise TypeError('binary_operation must be a BinaryOperation')

        self.binary_operation = binary_operation
        self.G = {}

        for g in binary_operation.codomain:
            ge = GroupElement(g, self)
            self.G[ge] = ge

        self.order = len(self.G)

        # get identity
        self.e = None
        for i in self.G:
            if all(i * a == a for a in self.G):
                self.e = i
                break

        if self.e is None:
            raise RuntimeError('G does not contain an identity element')

        # compute Caley table, get inverses (and check if Abelian in the meantime)
        self.abelian = True
        self.cayley_table = {}
        self.inverses = {}
        for gi in self.G:
            self.cayley_table[gi] = {}
            inverse_to_gi = False
            for gj in self.G:
                g = gi * gj
                try:
                    self.cayley_table[gi][gj] = self.G[g]
                except KeyError:
                    raise NotInGroup(g)
                if g == self.e:
                    if not inverse_to_gi:
                        self.inverses[gi] = gj
                        inverse_to_gi = True
                    else:
                        raise ValueError('{} have two inverses'.format(gi))

                if self.abelian and gi != gj:
                    try:
                        self.abelian = self.cayley_table[gi][gj] == self.cayley_table[gj][gi]
                    except KeyError:
                        pass

        if len(self.inverses) != self.order:
            raise ValueError('each element does not have an inverse')

        # find conjugacy classes:
        # for each A and B, computes A⁻¹ * B * A and mix the result in the conjugacy class
        self.conjugacy_classes = []
        self.to_conjugacy_class = {}
        for gi in self.G:
            if gi in self.to_conjugacy_class:
                continue
            else:
                klass = []
                for gj in self.G:
                    p = self.cayley_table[self.cayley_table[self.inverse(gj)][gi]][gj]
                    if p not in klass:
                        self.to_conjugacy_class[p] = len(self.conjugacy_classes)
                        klass.append(p)
                self.conjugacy_classes.append(set(klass))

        self.number_of_class = len(self.conjugacy_classes)

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

        try:
            return self.inverses[element]
        except KeyError:
            return NotInGroup(element)

    def conjugacy_class(self, g):
        """Get the elements of the conjugacy class of ``g``

        :param g: a group element
        :type g: GroupElement
        """

        try:
            index = self.to_conjugacy_class[g]
        except KeyError:
            raise NotInGroup(g)

        if self.abelian:
            return [g]
        else:
            return self.conjugacy_classes[index]

    def __contains__(self, item):
        """is an element part of the group ?

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


class Operation:
    """Define a symmetry element (fully based on quaternion).

    According `Mamone et al. <https://www.mdpi.com/2073-8994/2/3/1423/pdf>`_,

    .. math::

        S^k_n = \\hat{R}\\left(\\pi + \\frac{2\\,k\\,\\pi}{n}\\right),

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
        2. Take the fraction corresponding to the 4 components of the quaternion ;
        3. Add the improper variable ;
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

            # print(self, *(closest_fraction(x) for x in q), self._hash)

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

    def matrix_representation(self):
        """Get a 3x3 matrix representation of the symmetry operation
        """

        mat = quaternions.quat2mat(self.q)
        if self.improper:
            mat = -mat

        return mat

    def get_description(self, force=False):
        """Try to recognize the transformation:

        1. Extract axis and angle from quaternion (using ``quat2axangle()``) ;
        2. Treat angle (:math:`\\theta`):

           - if the axis have a negative z, take the opposite and set :math:`\\theta = 2\\pi-\\theta`,
           - if this is an improper rotation, set :math:`\\theta=\\theta-\\pi`,
           - Cast between 0 and :math:`4\\pi`,
           - divide by :math:`2\\pi` (so that the result is clamped between 0 and 2).

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
        """Try to find the main axis:

        - Find sole axis (x, y, z) ;
        - Find the axis at 45° with respect to two axis ;
        - Find the axis at 45° with respect to three axis (the 8 apexes of a cube).

        .. code-block:: text

            Some example of the different axis:

                (x'yz') +------*-------+ (xyz')
                       /              /|
                (x'y) *      o (y)   * |
                     /              /  * (xz')
             (x'yz) +------*-------+   |
                    |              | o |
                    |              |   + (xy',z')    y
              (x'z) *      o (z)   *  /              |
                    |              | * (y'z)         +--- x
                    |              |/               /
            (x'y'z) +------*-------+ (xy'z)        z

            with

            o : single axis,
            * : double axis,
            + : triple axis.

        :param axis: the axis
        :type axis: numpy.ndarray
        :rtype: str|None
        """

        s2 = math.sqrt(2) / 2
        s3 = math.sqrt(3) / 3

        axes = {
            # single (3)
            'z': [0, 0, 1.],
            'y': [0, 1., 0],
            'x': [1., 0, 0],
            # double (12)
            '[xy]': [s2, s2, 0],
            '[xz]': [s2, 0, s2],
            '[yz]': [0, s2, s2],
            '[x\'y]': [-s2, s2, 0],
            '[x\'z]': [-s2, 0, s2],
            '[y\'z]': [0, -s2, s2],
            '[xy\']': [s2, -s2, 0],
            '[xz\']': [s2, 0, -s2],
            '[yz\']': [0, s2, -s2],
            '[x\'y\']': [-s2, -s2, 0],
            '[x\'z\']': [-s2, 0, -s2],
            '[y\'z\']': [0, -s2, -s2],
            # triple (8)
            '[xyz]': [s3, s3, s3],
            '[xy\'z]': [s3, -s3, s3],
            '[x\'yz]': [-s3, s3, s3],
            '[x\'y\'z]': [-s3, -s3, s3],
            '[xyz\']': [s3, s3, -s3],
            '[xy\'z\']': [s3, -s3, -s3],
            '[x\'yz\']': [-s3, s3, -s3],
            '[x\'y\'z\']': [-s3, -s3, -s3],
        }

        for a, value in axes.items():
            if numpy.allclose(value, axis):
                return a


class PointGroupType(str, Enum):
    # cyclic (https://en.wikipedia.org/wiki/Cyclic_symmetry_in_three_dimensions):
    cyclic = 'C_n'
    pyramidal = 'C_nv'
    reflexion = 'C_nh'
    improper_rotation = 'S_n'
    # dihedral (https://en.wikipedia.org/wiki/Dihedral_symmetry_in_three_dimensions)
    dihedral = 'D_n'
    prismatic = 'D_nh'
    antiprismatic = 'D_nd'
    # cubic:
    # - tetrahedral (https://en.wikipedia.org/wiki/Tetrahedral_symmetry),
    # - octahedral (https://en.wikipedia.org/wiki/Octahedral_symmetry).
    tetrahedral_chiral = 'T'
    pyritohedral = 'T_h'
    tetrahedral_achiral = 'T_d'
    octahedral_chiral = 'O'
    octahedral_achiral = 'O_h'
    # icosahedral (https://en.wikipedia.org/wiki/Icosahedral_symmetry)
    icosahedral_chiral = 'I'
    icosahedral_achiral = 'I_h'


class PointGroupDescription:
    def __init__(self, symbol, n=0):
        if symbol not in PointGroupType:
            raise ValueError('unrecognized group type')

        if symbol in [
            PointGroupType.cyclic, PointGroupType.pyramidal, PointGroupType.reflexion,
                PointGroupType.improper_rotation, PointGroupType.dihedral, PointGroupType.prismatic,
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
    """Concrete implementation of a point group"""

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

        :param n: order of the group
        :type n: int
        :rtype: PointGroup
        """
        return cls.generate([Operation.C(n)], description=PointGroupDescription(PointGroupType.cyclic, n))

    @classmethod
    def S_n(cls, n=1):
        """Create a reflexion (odd n) or improper rotation (even n) point group of order n.

        Note that :math:`S_{2}` is usually noted :math:`C_i`.

        :param n: order of the group
        :type n: int
        :rtype: PointGroup
        """

        if n % 2 == 0:
            type_ = PointGroupType.improper_rotation
        else:
            type_ = PointGroupType.reflexion

        return cls.generate(
            [Operation.S(n)],
            description=PointGroupDescription(type_, n))

    @classmethod
    def C_nv(cls, n=2):
        """Create a pyramidal group of order n:

        .. math::

            C_{nv} = C_n \\times \\sigma_v.

        :param n: order of the group
        :type n: int
        :rtype: PointGroup
        """
        if n < 2:
            raise ValueError('minimum order is 2')

        return cls.generate(
            [Operation.C(n), Operation.sigma(axis=numpy.array([1., 0, 0]))],
            description=PointGroupDescription(PointGroupType.pyramidal, n))

    @classmethod
    def C_nh(cls, n=1):
        """Create a reflexion point group of order n:

        .. math::

            C_{nv} = C_n \\times \\sigma_h.

        Note that :math:`C_{1h}` is usually noted :math:`C_s`.

        :param n: order of the group
        :type n: int
        :rtype: PointGroup
        """

        return cls.generate(
            [Operation.C(n), Operation.sigma()],
            description=PointGroupDescription(PointGroupType.reflexion, n))

    @classmethod
    def D_n(cls, n=2):
        """Create a dihedral point group of order n:

        .. math::

            D_{n} = C_n(z) \\times C_2(x).

        :param n: order of the group
        :type n: int
        :rtype: PointGroup
        """
        if n < 2:
            raise ValueError('minimum order is 2')

        return cls.generate(
            [Operation.C(n), Operation.C(2, axis=numpy.array([1., 0, 0]))],
            description=PointGroupDescription(PointGroupType.dihedral, n))

    @classmethod
    def D_nh(cls, n=2):
        """Create a prismatic point group of order n:

        .. math::

            D_{nh} &= D_n \\times \\sigma_h \\\\
            &= C_n(z) \\times C_2(x) \\times\\sigma_h.

        :param n: order of the group
        :type n: int
        :rtype: PointGroup
        """
        if n < 2:
            raise ValueError('minimum order is 2')

        return cls.generate(
            [Operation.C(n), Operation.C(2, axis=numpy.array([1., 0, 0])), Operation.sigma()],
            description=PointGroupDescription(PointGroupType.prismatic, n))

    @classmethod
    def D_nd(cls, n=2):
        """Create an anti prismatic point group of order n:

        .. math::

            D_{nd} = S_{2n}(z) \\times C_2(x).

        :param n: order of the group
        :type n: int
        :rtype: PointGroup
        """
        if n < 2:
            raise ValueError('minimum order is 2')

        return cls.generate(
            [Operation.S(2 * n), Operation.C(2, axis=numpy.array([1., 0, 0]))],
            description=PointGroupDescription(PointGroupType.antiprismatic, n))

    @classmethod
    def T(cls):
        """Create the chiral tetrahedral group

        .. math::

            T = C_2(z) \\times C_3(1, 1, 1).

        :rtype: PointGroup
        """

        return cls.generate(
            [
                Operation.C(2),
                Operation.C(3, axis=numpy.array([1., 1., 1.]))],
            description=PointGroupDescription(PointGroupType.tetrahedral_chiral))

    @classmethod
    def T_h(cls):
        """Create the pyritohedral group

        .. math::

            T_h = C_2 \\times \\sigma_h \\times C_3(1, 1, 1).

        :rtype: PointGroup
        """

        return cls.generate(
            [
                Operation.C(2),
                Operation.sigma(),
                Operation.C(3, axis=numpy.array([1., 1., 1.]))],
            description=PointGroupDescription(PointGroupType.pyritohedral))

    @classmethod
    def T_d(cls):
        """Create the achiral tetrahedral group

        .. math::

            T_d = S_4 \\times C_3(1, 1, 1).

        :rtype: PointGroup
        """

        return cls.generate(
            [
                Operation.S(4),
                Operation.C(3, axis=numpy.array([1., 1., 1.]))],
            description=PointGroupDescription(PointGroupType.tetrahedral_achiral))

    @classmethod
    def O(cls):  # noqa
        """Create the chiral octahedral group

        .. math::

            O = C_4(z) \\times C_3(1, 1, 1).

        :rtype: PointGroup
        """

        return cls.generate(
            [
                Operation.C(4),
                Operation.C(3, axis=numpy.array([1., 1., 1.]))],
            description=PointGroupDescription(PointGroupType.octahedral_chiral))

    @classmethod
    def O_h(cls):
        """Create the achiral octahedral group

        .. math::

            O_h = C_4 \\times \\sigma_h \\times C_3(1, 1, 1).

        :rtype: PointGroup
        """

        return cls.generate(
            [
                Operation.C(4),
                Operation.sigma(),
                Operation.C(3, axis=numpy.array([1., 1., 1.]))],
            description=PointGroupDescription(PointGroupType.octahedral_achiral))


class SymmetryFinder:
    """Find the symmetry

    Inspired by https://github.com/sunqm/pyscf/blob/master/pyscf/symm/geom.py

    :param points: an Nx4 array of N points with ``(Z, x, y, z)`` for each point. ``Z`` is a label for points
        that are equivalent.
    :type points: numpy.ndarray
    :param tol: tolerance threshold
    :type tol: float
    """

    def __init__(self, points, tol=1e-5):
        self.points = points
        self.tol = tol

        # find center and translate
        tot_Z = .0
        tot_coord = numpy.zeros(3)

        for i in range(len(points)):
            tot_coord += points[i, 1:]
            tot_Z = points[i, 0]

        self.center = tot_coord / tot_Z
        self.points[:, 1:] -= self.center

        # group points
        self.group_points_per_distance = []
        all_indexes = numpy.asarray(range(len(points)))

        decimals = int(-numpy.log10(tol)) - 1

        uniques, indexes = numpy.unique(points[:, 0], return_inverse=True)
        for u, _ in enumerate(uniques):  # group per label
            sub_indexes = all_indexes[indexes == u]
            sub_points = self.points[sub_indexes, 1:]

            uniques_2, indexes_2 = numpy.unique(
                numpy.around(numpy.linalg.norm(sub_points, axis=1), decimals), return_inverse=True)

            for u2, _ in enumerate(uniques_2):  # group per distance
                self.group_points_per_distance.append(sub_indexes[indexes_2 == u2])

    @staticmethod
    def vec_in_vecs(vec, vecs, tol):
        """Check if ``vec`` is in the set of ``vecs``

        :param vec: vector to find
        :param vec: numpy.ndarray
        :param vecs: set of vectors
        :type vecs: numpy.ndarray
        :param tol: tolerence threshold
        :type tol: float
        :rtype: bool
        """
        return min(numpy.sum(numpy.abs(v - vec)) for v in vecs) < tol

    def _symetric_for(self, op):
        """Check if symmetric for a given operation

        :param op: the operation
        :type op: numpy.ndarray
        """

        for lst in self.group_points_per_distance:
            p = self.points[lst, 1:]
            yield all(SymmetryFinder.vec_in_vecs(x, p, self.tol) for x in numpy.dot(p, op))

    def symmetric_for(self, op):
        """Check if symmetric for a given operation

        :param op: the operation
        :type op: Operation
        """

        return all(self._symetric_for(op.matrix_representation()))

    def has_inversion(self):
        return self.symmetric_for(Operation.i())

    def has_rotation(self, axis, n):
        return self.symmetric_for(Operation.C(n, axis=axis))

    def has_mirror(self, axis):
        return self.symmetric_for(Operation.sigma(axis))

    def has_improper_rotation(self, n, axis):
        return self.symmetric_for(Operation.S(n, axis=axis))
