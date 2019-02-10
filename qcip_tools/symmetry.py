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

from qcip_tools import math as qcip_math


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

        Algorithm from https://physics.stackexchange.com/a/351400j, simplified version of
        http://pymatgen.org/_modules/pymatgen/symmetry/analyzer.html

        :param generators; the generators of the set
        :type generators: list
        :param func: function
        :type func: function
        :rtype: BinaryOperation
        """

        full = list(generators)
        for g in full:
            for s in generators:
                op = g * s
                if op not in full:
                    full.append(op)

        return cls(Set(full), func)


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
        meet = {}
        self.to_conjugacy_class = {}
        for gi in self.G:
            if gi in meet:
                continue
            else:
                klass = {}
                for gj in self.G:
                    p = self.cayley_table[self.cayley_table[self.inverse(gj)][gi]][gj]
                    if p not in klass:
                        meet[p] = klass[p] = True
                self.conjugacy_classes.append(set(klass.keys()))

        self._sort_classes()
        self.number_of_class = len(self.conjugacy_classes)

        for i, elmts in enumerate(self.conjugacy_classes):
            for e in elmts:
                self.to_conjugacy_class[e] = i

    def _sort_classes(self):
        """Sort the conjugacy classes, if needed
        """

        pass

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

    def _sort_classes(self):
        """Reimplemented to sort the different class in a given order:

        1. :math:`E` comes first,
        2. Then comes proper rotations,

            a. Proper rotation around z axis comes first, sorted by reverse order ;
            b. Then comes other rotation.

        3. Then comes improper rotations,

            a. :math:`i` comes first ;
            b. Then comes rotoreflexion axis of order :math:`n >2` ;
            c. Finally, reflexion plane, in the order :math:`\\sigma_h`, :math:`\\sigma_v`, :math:`\\sigma_d`.

        Since the order is arbitrary in character tables,
        **this will probably not correspond to the usual "textbook" order**.
        """

        def s(a):
            if self.e in a:  # E is first
                return 0, 0, 0, 0
            else:
                i = iter(a)
                x = next(i).element
                d = x.get_description()
                if not x.improper:
                    return 1, numpy.around(-d.axis[-1], 2), -d.n, d.k
                else:
                    t = 0
                    e = 0
                    if d.n > 2:
                        t = -d.n
                        e = d.k
                        if len(a) > 1 and d.k > d.n - d.k:
                            e = d.n - d.k
                    elif d.n == 1:
                        t = 100
                        e = 0 if d.textual_axis in ['x', 'y'] else 1

                    return 2, t, numpy.around(-d.axis[-1], 2), e

        self.conjugacy_classes.sort(key=lambda a: s(a))

    @staticmethod
    def product(e):
        return e[0] * e[1]

    def gen_character_table(self):
        """Try to generate the character table
        """

        pass

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
        """Create a reflexion (odd n, :math:`C_{(2k+1)h}`)
        or improper rotation (even n, :math:`S_{2k}`) point group of
        order n.

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

    @classmethod
    def I(cls):  # noqa
        """Create the chiral octahedral group

        .. math::

            I = C_5(1, 0, \\phi) \\times C_3(1, 1, 1) \\times C_2(z), \\text{ with }
            \\phi = \\frac{1+\\sqrt 5}{2}.

        :rtype: PointGroup
        """
        golden = (1 + 5 ** 0.5) / 2
        return cls.generate(
            [
                Operation.C(2),
                Operation.C(3, axis=numpy.array([1., 1., 1.])),
                Operation.C(5, axis=numpy.array([1, 0, golden]))
            ],
            description=PointGroupDescription(PointGroupType.icosahedral_chiral))

    @classmethod
    def I_h(cls):
        """Create the chiral octahedral group

        .. math::

            I_h = C_5(1, 0, \\phi) \\times C_3(1, 1, 1) \\times C_2(z) \\times \\sigma_h, \\text{ with }
            \\phi = \\frac{1+\\sqrt 5}{2}.

        :rtype: PointGroup
        """
        golden = (1 + 5 ** 0.5) / 2
        return cls.generate(
            [
                Operation.C(2),
                Operation.C(3, axis=numpy.array([1., 1., 1.])),
                Operation.C(5, axis=numpy.array([1, 0, golden])),
                Operation.sigma(),
            ],
            description=PointGroupDescription(PointGroupType.icosahedral_achiral))


class PointGroupError(Exception):
    pass


class PointGroupDescription:
    def __init__(self, symbol, n=0):
        if symbol not in PointGroupType:
            raise ValueError('unrecognized group type')

        if symbol in [
            PointGroupType.cyclic, PointGroupType.pyramidal, PointGroupType.reflexion,
                PointGroupType.improper_rotation, PointGroupType.dihedral, PointGroupType.prismatic,
                PointGroupType.antiprismatic]:

            if n < 1 and n != -1:
                raise ValueError('order must be equal or superior to 1')
        else:
            if n > 0:
                raise ValueError('order does not matter for point group {}'.format(symbol.value))

        self.symbol = symbol
        self.order = n

    def __str__(self):
        s = self.symbol.value
        if 'n' in s:
            if self.order > 0:
                s = s.replace('n', str(self.order))
            else:
                s = s.replace('n', 'oo')
        return s

    def gen_point_group(self, lowers_infinite=-1):
        """Generate the point group out of the description.

        Not able to generate infinite groups :math:`D_{\\infty h}` or :math:`C_{\\infty v}`, except
        if ``lowers_infinite`` is set.

        :param lowers_infinite: axis order used instead of :math:`\\infty`.
        :type lowers_infinite: int
        :rtype: PointGroup
        """
        mapping = {
            # with order:
            PointGroupType.cyclic: PointGroup.C_n,
            PointGroupType.pyramidal: PointGroup.C_nv,
            PointGroupType.reflexion: PointGroup.C_nh,
            PointGroupType.improper_rotation: PointGroup.S_n,
            PointGroupType.dihedral: PointGroup.D_n,
            PointGroupType.prismatic: PointGroup.D_nh,
            PointGroupType.antiprismatic: PointGroup.D_nd,
            # no order:
            PointGroupType.tetrahedral_chiral: PointGroup.T,
            PointGroupType.pyritohedral: PointGroup.T_h,
            PointGroupType.tetrahedral_achiral: PointGroup.T_d,
            PointGroupType.octahedral_chiral: PointGroup.O,
            PointGroupType.octahedral_achiral: PointGroup.O_h,
            PointGroupType.icosahedral_chiral: PointGroup.I,
            PointGroupType.icosahedral_achiral: PointGroup.I_h,
        }

        if self.symbol not in mapping:
            raise PointGroupError('cannot generate {}'.format(self.symbol))

        if self.order < 0:
            if lowers_infinite < 2:
                raise PointGroupError('cannot generate infinite groups')
            else:
                mapping[self.symbol](lowers_infinite)
        elif self.order == 0:
            return mapping[self.symbol]()
        else:
            return mapping[self.symbol](self.order)


class SymmetryFinderError(Exception):
    pass


class SymmetryFinder:
    """Find the symmetry

    Inspired by https://github.com/sunqm/pyscf/blob/master/pyscf/symm/geom.py
    and http://pymatgen.org/_modules/pymatgen/symmetry/analyzer.html

    :param points: an Nx4 array of N points with ``(Z, x, y, z)`` for each point. ``Z`` is a label for points
        that are equivalent.
    :type points: numpy.ndarray
    :param tol: tolerance threshold
    :type tol: float
    """

    def __init__(self, points, tol=1e-5):
        self.points = points
        self.tol = tol
        self.decimals = int(-numpy.log10(tol)) - 1

        # find center and translate
        self.center = numpy.einsum('i,ij->j', points[:, 0], points[:, 1:]) / points[:, 0].sum()
        self.points[:, 1:] -= self.center

        # group points
        self.group_points_per_distance = []
        all_indexes = numpy.asarray(range(len(points)))

        uniques, indexes = numpy.unique(points[:, 0], return_inverse=True)
        for u, _ in enumerate(uniques):  # group per label
            sub_indexes = all_indexes[indexes == u]
            sub_points = self.points[sub_indexes, 1:]

            uniques_2, indexes_2 = numpy.unique(
                numpy.around(numpy.linalg.norm(sub_points, axis=1), self.decimals), return_inverse=True)

            for u2, _ in enumerate(uniques_2):  # group per distance
                self.group_points_per_distance.append(sub_indexes[indexes_2 == u2])

    @staticmethod
    def cartesian_tensor(tensor, n):
        """Compute multipole expansion of cartesian tensor of order ``n`` (what you find in the right column of a
        character table) and get its principal components.

        https://en.wikipedia.org/wiki/Multipole_expansion

        :param tensor: tensor
        :type tensor: numpy.ndarray
        :param n: order
        :type n: int
        :return: eigenvalue and eigenvectors of the tensor
        """

        z = tensor[:, 0]
        r = tensor[:, 1:]
        ncart = (n + 1) * (n + 2) // 2
        t = numpy.sqrt(numpy.copy(z).reshape(len(z), -1)) / z.sum()
        for i in range(n):
            t = numpy.einsum('zi,zj->zij', t, r).reshape(len(z), -1)

        e, c = numpy.linalg.eigh(numpy.dot(t.T, t))
        return e[-ncart:], c[:, -ncart:]

    @staticmethod
    def inertia_tensor_moments(tensor):
        """Get inertia egeinvalue and egeinvectors, sorted
        (smallest first, so the main axis is the last eigenvector, if any).

        :param tensor: carthesian tensor
        :type tensor: numpy.ndarray
        :rtype: tuple
        """
        w = tensor[:, 0]

        t = numpy.zeros((3, 3))
        for i in range(3):
            t[i, i] = numpy.sum(w * (tensor[:, 1 + (i + 1) % 3] ** 2 + tensor[:, 1 + (i + 2) % 3] ** 2))

        for i, j in [(0, 1), (0, 2), (1, 2)]:
            x = -numpy.sum(w * tensor[:, 1 + i] * tensor[:, 1 + j])
            t[i, j] = t[j, i] = x

        e, c = numpy.linalg.eigh(t)
        indices = numpy.argsort(e)
        e = e[indices]
        c = c.T[indices]
        return e, c

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

    @staticmethod
    def degeneracies(e, decimals=-5):
        e = numpy.around(e, decimals)
        u, idx = numpy.unique(e, return_inverse=True)
        return [(numpy.count_nonzero(idx == i), u[i]) for i in range(len(u))]

    def find_symmetry(self):
        """(Try to) find symmetry.

        Algorithm is inspired by http://pymatgen.org/_modules/pymatgen/symmetry/analyzer.html (and actually any
        book on group theory):

        1. Find inertia tensor, its eigenvalues and eigenvectors (axes)
        2. Based on the eigenvalues:

          a. Linear molecules have one zero egeinvalue, and belongs to :math:`C_{\\infty v}`
             or :math:`D_{\\infty h}` ;
          b. Spherical top molecules have three equal (nonzero) eigenvalues, which corresponds to :math:`T`,
             :math:`I` or :math:`O` groups ;
          c. Asymmetric top molecules have all different (non zero) eigenvalues, which corresponds to groups with
             at most rotation axis of order 2 (:math:`D_{2d}` or subgroups), and here they are treated
             like symmetric top (since it is essentially the same thing) ;
          d. Symmetric top molecules have two equals eigenvalues, thus corresponding to axial point group with (maybe)
             a unique principal rotation axis. To discriminate a little bit more,

             - Either there is perpendicular :math:`C_2` axes, and the group is a dihedral one
               (depending on the presence or not of mirror planes, :math:`D_{nh}`, :math:`D_{nd}` or :math:`D_n`) ;
             - There is a main axis: depending on the presence of mirror plane, one can discriminate between
               the different cyclic groups (:math:`C_{nh}`, :math:`C_{nv}`, :math:`S_{2n}` or :math:`C_n`) ;
             - There is no main axis, and the group can be either :math:`C_s`, :math:`C_i` or :math:`C_1`.

        TODO: deal with orientation (main axis versus others), check against real molecules

        :return: the point group, the center and the rotational matrix
        :rtype: tuple(PointGroupDescription, numpy.ndarray, numpy.ndarray)
        """

        e, v = SymmetryFinder.inertia_tensor_moments(self.points)
        degeneracies = SymmetryFinder.degeneracies(e, self.decimals)

        group = PointGroupDescription(PointGroupType.cyclic, 1)

        if len(degeneracies) == 1 and degeneracies[0][1] < self.tol:
            raise SymmetryFinderError('belongs to SO(3) rotation group')
        elif degeneracies[0][1] < self.tol:  # linear
            if self.has_inversion():
                group = PointGroupDescription(PointGroupType.prismatic, -1)  # Dooh
            else:
                group = PointGroupDescription(PointGroupType.pyramidal, -1)  # Coov

        elif len(degeneracies) == 1:  # spherical top
            probable_cn = self.find_probable_rotations()

            if len(probable_cn) == 0:
                raise SymmetryFinderError('accidentally spherical top, no rotation')

            n, main_axis = self.find_c_highest(probable_cn)
            other_axes = v.T[numpy.where(abs(numpy.dot(v.T, main_axis)) < self.tol)]

            if n < 3:
                raise SymmetryFinderError('accidentally spherical top, n < 3')
            elif n == 3:
                group = PointGroupDescription(PointGroupType.tetrahedral_chiral)  # T
                if self.has_inversion():
                    group = PointGroupDescription(PointGroupType.pyritohedral)  # T_h
                elif self.has_mirror_type(main_axis, other_axes, 'd'):
                    group = PointGroupDescription(PointGroupType.tetrahedral_achiral)  # T_d
            elif n == 4:
                group = (PointGroupType.octahedral_chiral)  # O
                if self.has_inversion():
                    group = PointGroupDescription(PointGroupType.octahedral_achiral)  # O_h
            elif n == 5:
                group = PointGroupDescription(PointGroupType.icosahedral_chiral)  # I
                if self.has_inversion():
                    group = PointGroupDescription(PointGroupType.icosahedral_achiral)  # I_h
            else:
                raise SymmetryFinderError('spherical top molecule with n>5')

        else:  # (a)symmetric top
            probable_cn = self.find_probable_rotations()
            if len(probable_cn) == 0:
                if self.has_inversion():
                    group = PointGroupDescription(PointGroupType.improper_rotation, 2)  # C_i = S_2
                else:
                    for j in range(3):
                        if self.has_mirror(v.T[j]):
                            group = PointGroupDescription(PointGroupType.reflexion, 1)  # C_s = C_1h
                            break
            else:
                n, main_axis = self.find_c_highest(probable_cn)
                other_axes = v.T[numpy.where(abs(numpy.dot(v.T, main_axis)) < self.tol)]
                if len(probable_cn) >= 2 and self.has_perpendicular_C2(main_axis, probable_cn):
                    if self.has_mirror(main_axis):
                        group = PointGroupDescription(PointGroupType.prismatic, n)  # D_nh
                    elif self.has_mirror_type(main_axis, other_axes, 'd'):
                        group = PointGroupDescription(PointGroupType.antiprismatic, n)  # D_nd
                    else:
                        group = PointGroupDescription(PointGroupType.dihedral, n)  # D_n
                else:
                    if self.has_mirror(main_axis):
                        group = PointGroupDescription(PointGroupType.reflexion, n)  # C_nh
                    elif self.has_mirror_type(main_axis, other_axes, 'v'):
                        group = PointGroupDescription(PointGroupType.pyramidal, n)  # C_nv
                    elif self.has_improper_rotation(axis=main_axis, n=2 * n):
                        group = PointGroupDescription(PointGroupType.improper_rotation, 2 * n)  # S_2n
                    else:
                        group = PointGroupDescription(PointGroupType.cyclic, n)  # C_n

        return group, self.center, v

    def find_probable_rotations(self, parallel_to=None):
        """Find rotations

        Inspired by https://github.com/sunqm/pyscf/blob/master/pyscf/symm/geom.py !

        The idea is to find, for each set of equidistant points, the ones that are at equal distance from the
        first point, which means that a rotation axis may pass by those points.

        :param parallel_to: get only the rotation axis parallel to a given axis (should be normalized)
        :type parallel_to: numpy.ndarray|list
        """

        # find possible C_n
        maybe_cn = []
        for lst in self.group_points_per_distance:
            if len(lst) < 2:
                continue

            # find C_2
            coords = self.points[lst, 1:]
            for i in range(1, len(lst)):
                if abs(coords[0] + coords[i]).sum() > self.tol:
                    maybe_cn.append((2, coords[0] + coords[i]))
                else:  # abs(coords[0] - coords[i]).sum() > self.tol
                    maybe_cn.append((2, coords[0] - coords[i]))

            # Find C_n n > 2
            r0 = coords - coords[0]
            distances = numpy.linalg.norm(r0, axis=1)
            d = abs(distances[:, None] - distances) < self.tol
            for i in range(2, len(lst)):
                for j in numpy.where(d[i, :i])[0]:
                    ang = qcip_math.angle_vector(r0[i], r0[j])
                    n = 2 * numpy.pi / ang
                    if abs(int(n) - n) < self.tol:
                        maybe_cn.append((int(n), numpy.cross(r0[i], r0[j])))

        # remove null vectors
        v = numpy.vstack([x[1] for x in maybe_cn])
        ns = numpy.hstack([x[0] for x in maybe_cn])
        indexes = numpy.linalg.norm(v, axis=1) > self.tol  # no more zero length
        v = v[indexes] / numpy.linalg.norm(v[indexes], axis=1).reshape(-1, 1)  # normalize
        ns = ns[indexes]

        # remove non-parallel (dot product with parallel axis is 1 or -1)
        if parallel_to:
            d = numpy.dot(v, parallel_to)
            indexes = (abs(d + 1.) < self.tol) | (abs(d - 1.) < self.tol)
            v = v[indexes]
            ns = ns[indexes]

        # remove duplicates (opposite direction), if any
        probable_cn = []
        seen = numpy.zeros(len(v), dtype=bool)
        for ki, vi in enumerate(v):
            if not seen[ki]:
                w1 = numpy.where(numpy.einsum('ij->i', abs(v[ki:] - vi)) < self.tol)[0] + ki
                w2 = numpy.where(numpy.einsum('ij->i', abs(v[ki:] + vi)) < self.tol)[0] + ki
                seen[w1] = True
                seen[w2] = True
                e = numpy.einsum('ix->x', v[w1]) - numpy.einsum('ix->x', v[w2])  # average axis
                e /= numpy.linalg.norm(e)
                for n in set(ns[w1]) | set(ns[w2]):
                    probable_cn.append((n, e))

        return probable_cn

    def find_c_highest(self, probable_cn):
        """Among all rotation axis, find highest.

        :param probable_cn: list of rotation axis ``(order, axis)``
        :type probable_cn: list(tuple)
        """
        n_max = 1
        axis_max = numpy.array([0, 0, 1.])
        for n, axis in probable_cn:
            if n > n_max and self.has_rotation(axis, n):
                n_max = n
                axis_max = axis

        return n_max, axis_max

    def has_perpendicular_C2(self, perpendicular_to, probable_cn):
        """Check if a perpendicular :math:`C_2` actually exists

        :param probable_cn: list of rotation axis ``(order, axis)``
        :type probable_cn: list(tuple)
        :param perpendicular_to: main axis
        :type perpendicular_to: numpy.ndarray
        """
        cn = numpy.vstack([x[1] for x in probable_cn if x[0] == 2.])  # keep C_2
        cn = cn[numpy.where(abs(numpy.einsum('ij,j->i', cn, perpendicular_to)) < self.tol)]  # keep perpendicular
        for i in range(len(cn)):
            if self.has_rotation(cn[i], 2):
                return True
        return False

    def find_probable_parallel_mirrors(self, parallel_to):
        """Find possible mirrors parallel to a given axis (so that the normal is perpendicular).

        For each combination of points, find if the corresponding normal is perpendicular
        to the axis (dot product is null), then remove duplicates.

        Returns a list of normal.

        :param parallel_to: axis to which the plane should be parallel to
        :type parallel_to: numpy.ndarray|list
        """
        possible_mirrors = []

        for lst in self.group_points_per_distance:
            if len(lst) < 2:
                continue

            r = self.points[lst, 1:]
            for i in range(1, len(lst)):
                re = r - r[i]
                w = abs(numpy.einsum('ij,j->i', re, parallel_to)) < self.tol
                for j in numpy.where(w)[0]:
                    if i == j:
                        continue
                    possible_mirrors.append(re[j])

        v = numpy.vstack(possible_mirrors)
        v /= numpy.linalg.norm(v, axis=1).reshape(-1, 1)  # normalize

        # remove duplicates
        probable_mirrors = []
        seen = numpy.zeros(len(v), dtype=bool)
        for ki, vi in enumerate(v):
            if not seen[ki]:
                w1 = numpy.where(numpy.einsum('ij->i', abs(v[ki:] - vi)) < self.tol)[0] + ki
                w2 = numpy.where(numpy.einsum('ij->i', abs(v[ki:] + vi)) < self.tol)[0] + ki
                seen[w1] = True
                seen[w2] = True
                e = numpy.einsum('ix->x', v[w1]) - numpy.einsum('ix->x', v[w2])  # average axis
                probable_mirrors.append(e / numpy.linalg.norm(e))

        return probable_mirrors

    def has_mirror_type(self, parallel_to, other_axes, mirror_type):
        """Find the different mirror types

        :param parallel_to: principal axis
        :type parallel_to: numpy.ndarray
        :param other_axes: list of other axis (to determine whether the mirror is "v" or "d")
        :type other_axes: numpy.ndarray
        :param mirror_type: mirror_type
        :type mirror_type: str
        """

        if mirror_type == 'h':
            return self.has_mirror(parallel_to)

        # find vertical/diagonal
        probable_mirrors = self.find_probable_parallel_mirrors(parallel_to)
        for m in probable_mirrors:
            if self.has_mirror(m):
                mirror_type_found = 'd'
                for axis in other_axes:
                    if numpy.dot(m, axis) < self.tol:
                        mirror_type_found = 'v'
                        break

                if mirror_type_found == mirror_type:
                    return True

        return False

    def has_inversion(self):
        return self.symmetric_for(Operation.i())

    def has_rotation(self, axis, n):
        return self.symmetric_for(Operation.C(n, axis=axis))

    def has_mirror(self, axis):
        return self.symmetric_for(Operation.sigma(axis))

    def has_improper_rotation(self, axis, n):
        return self.symmetric_for(Operation.S(n, axis=axis))
