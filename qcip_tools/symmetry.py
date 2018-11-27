"""
Symmetry handling.

See https://en.wikipedia.org/wiki/Group_(mathematics)

See https://physics.stackexchange.com/questions/351372/generate-all-elements-of-a-point-group-from-generating-set
and https://github.com/naftaliharris/Abstract-Algebra/blob/master/absalg/Group.py
and https://pdfs.semanticscholar.org/4edd/1ac673ea4cab251bb0b64bf0f5b65f18cc9d.pdf
and http://www.euclideanspace.com/maths/geometry/affine/reflection/quaternion/index.htm
"""

import numpy
import math


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

        return all(self(e) in self.codomain for e in self.codomain * self.codomain)

    def check_associativity(self):
        """Check if the binary operation is associative

        :rtype: bool
        """

        return all(
            self((a, self((b, c)))) == self((self((a, b)), c))
            for a in self.codomain for b in self.codomain for c in self.codomain)


class Group:
    """A set G, together with an operation * (group law), form a group (G, *) if it satisfies 4 requirements:

    - Closure (not checked here);
    - Associativity (not checked here) ;
    - (unique) Identity ;
    - Inverse.

    A group may also be Abelian.
    """

    def __init__(self, binary_operation):
        if not isinstance(binary_operation, BinaryOperation):
            raise TypeError('binary_operation must be a BinaryOperation')

        self.binary_operation = binary_operation
        self.G = binary_operation.codomain

        # get identity
        self.e = None
        for i in self.G:
            if all(self.binary_operation((i, a)) == a for a in self.G):
                self.e = i
                break

        if self.e is None:
            raise RuntimeError('G does not contain an identity element')

        # get inverses (and check if Abelian)
        self.inverses = {}
        self.abelian = True
        for i in self.G:
            inverse_found = False
            for j in self.G:
                if self.binary_operation((i, j)) == self.e:
                    self.inverses[i] = j
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

        :return: another element of G
        """

        return self.inverses[element]


class SymmetryElement:
    """Define a Symmetry element. Use the Schönflies notation (and restrict to rotations)
    """

    def __init__(self, name, matrix, n=1, k=1, axis='z'):
        self.matrix_representation = matrix
        self.name = name
        self.n = n
        self.k = k
        self.axis = axis

    def copy(self):
        return SymmetryElement(
            self.name, self.matrix_representation.copy, self.n, self.k, self.axis)

    def inverse(self):
        """Return the inverse element

        :rtype: SymmetryElement
        """
        if self.name in ['E', 'i', 'O']:  # those are their own inverse
            return self.copy()
        elif self.name == 'C':
            return SymmetryElement.C(self.n, self.n - self.k)
        elif self.name == 'S':
            k = (self.n if self.n % 2 == 0 else 2 * self.n) - self.k
            return SymmetryElement.S(self.n, k)
        else:  # other symmetry element
            return SymmetryElement(self.name, self.matrix_representation.inv(), n=self.n, k=self.k)

    def __str__(self):
        if self.name in ['E', 'i']:
            return self.name
        elif self.name in ['S', 'C']:
            return '{}_{}{}{}'.format(
                self.name,
                self.n,
                '({})'.format(self.axis),
                '^{}'.format(self.k) if self.k != 1 else '',
            )
        elif self.name == 'O':
            return '{}({})'.format(
                self.name,
                'xy' if self.axis == 'z' else ('xz' if self.axis == 'y' else ('yz' if self.axis == 'x' else 'r')))
        else:
            return '{}(n={},k={},axis={})'.format(self.name, self.n, self.k, self.axis)

    def __eq__(self, other):
        """Compare two symmetry elements on the basis of their matrix representation
        """

        if not isinstance(other, SymmetryElement):
            raise TypeError('other must be symmetry element')

        return numpy.allclose(self.matrix_representation, other.matrix_representation)

    def __pow__(self, power, modulo=None):
        if not isinstance(power, int):
            raise TypeError('power must be int')

        if power == 0:
            return SymmetryElement.E()
        if power < 0:
            return self.inverse() ** -power
        else:
            if self.name in ['E', 'i', 'O']:  # they are their own inverse
                if power % 2 == 0:
                    return SymmetryElement.E()
                else:
                    return self.copy()
            elif self.name == 'C':
                return SymmetryElement.C(self.n, self.k * power, self.axis)
            elif self.name == 'S':
                return SymmetryElement.S(self.n, self.k * power, self.axis)
            else:
                x = self.copy()
                x.matrix_representation = x.matrix_representation**power
                return x

    @classmethod
    def E(cls):
        """Generate identity

        :rtype: SymmetryElement
        """
        return cls('E', numpy.identity(3))

    @staticmethod
    def rotation_axis(n, k=1, axis='z', improper=False):
        """Generate the matrix representation of an (im)proper rotation (parallel to ``axis``)

        :param n: order of the axis
        :type n: int
        :param k: number of step
        :type k: int
        :param axis: to which axis this rotation axis should be parallel ?
        :type axis: str
        :param improper: is this an improper rotation?
        :type improper: bool
        :rtype: numpy.ndarray
        """

        a = 2 * numpy.pi / n * k
        sa, ca = math.sin(a), math.cos(a)

        i = 1 if not improper else (-1 if k % 2 == 1 else 1)

        mat = numpy.identity(3)

        if axis == 'z':
            mat[0, 0] = ca
            mat[0, 1] = -sa
            mat[1, 0] = sa
            mat[1, 1] = ca
            mat[2, 2] = i
        elif axis == 'y':
            mat[0, 0] = ca
            mat[0, 2] = sa
            mat[2, 0] = -sa
            mat[2, 2] = ca
            mat[1, 1] = i
        elif axis == 'x':
            mat[1, 1] = ca
            mat[1, 2] = -sa
            mat[2, 1] = sa
            mat[2, 2] = ca
            mat[0, 0] = i
        else:
            raise ValueError('axis must be x, y or z')

        return mat

    @classmethod
    def C(cls, n, k=1, axis='z'):
        """Generate the representation of a proper rotation (parallel to ẑ)

        :param n: order of the axis
        :type n: int
        :param k: number of step
        :type k: int
        :param axis: to which axis this rotation axis should be parallel ?
        :type axis: str
        :rtype: SymmetryElement
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

        return cls('C', SymmetryElement.rotation_axis(n, k, axis=axis), n, k, axis)

    @classmethod
    def S(cls, n, k=1, axis='z'):
        """Generate the representation of an improper rotation (parallel to axis)

        :param n: order of the axis
        :type n: int
        :param k: number of step
        :type k: int
        :param axis: to which axis this rotation axis should be parallel ?
        :type axis: str
        :rtype: SymmetryElement
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
            return cls.C(n, k, axis=axis)
        else:
            return cls('S', SymmetryElement.rotation_axis(n, k, improper=True, axis=axis), n, k, axis)

    @classmethod
    def sigma(cls, axis=None, parallel_to=None):
        """Generate a reflexion plane parallel to ``parallel_to`` or perpendicular to ``axis``

        :param parallel_to: the two axis to which the plane is parallel
        :type parallel_to: str
        :param axis: to which axis this plane should be parallel ?
        :type axis: str
        :rtype: SymmetryElement
        """

        if (axis is None and parallel_to is None) or (axis is not None and parallel_to is not None):
            raise ValueError('you must provide either axis or parallel_to')

        if parallel_to is not None:
            parallel_to = str(sorted(parallel_to))

            if parallel_to == 'xy':
                axis = 'z'
            elif parallel_to == 'xz':
                axis = 'y'
            elif parallel_to == 'yz':
                axis = 'x'
            else:
                ValueError('{} is not a valid definition'.format(parallel_to))

        return cls('O', SymmetryElement.rotation_axis(1, axis=axis, improper=True), axis=axis)

    @classmethod
    def i(cls):
        """Generate an inversion center

        :rtype: SymmetryElement
        """

        return cls('i', SymmetryElement.rotation_axis(2, improper=True))
