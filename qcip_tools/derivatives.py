import math
import collections
import itertools

from qcip_tools import math as qcip_math

#: In term of derivative of the energy
#: ``F`` = static electric field derivative, ``D`` = dynamic electric field derivative (which can be static),
#: ``G`` = geometrical derivative, ``N`` = normal mode derivative
ALLOWED_DERIVATIVES = ('F', 'D', 'G', 'N')

COORDINATES = {0: 'x', 1: 'y', 2: 'z'}  #: spacial 3D coordinates
COORDINATES_LIST = list(COORDINATES)
ORDERS = {1: 'first', 2: 'second', 3: 'third', 4: 'fourth', 5: 'fifth'}  #: number to x*th*.


class RepresentationError(Exception):
    pass


class Derivative:
    """Represent a quantity, which is derivable

        :param from_representation: representation of the derivative
        :type from_representation: str
        :param basis: basis for the representation, if any
        :type basis: Derivative
        :param spacial_dof: spacial degrees of freedom (3N)
        :type spacial_dof: int
        """

    def __init__(self, from_representation=None, basis=None, spacial_dof=None):

        self.basis = None
        self.spacial_dof = spacial_dof
        self.diff_representation = ''

        if basis:
            if isinstance(basis, Derivative):
                self.basis = basis
                if basis.spacial_dof:
                    if spacial_dof and basis.spacial_dof and spacial_dof != basis.spacial_dof:
                        raise ValueError(spacial_dof)
                    if not spacial_dof and basis.spacial_dof:
                        self.spacial_dof = basis.spacial_dof
            else:
                raise TypeError(basis)

        if from_representation is not None:
            for i in from_representation:
                if i not in ALLOWED_DERIVATIVES:
                    raise RepresentationError(from_representation)

                if i in 'GN' and not spacial_dof:
                    raise Exception('geometrical derivative and no spacial_dof !')

            self.diff_representation = from_representation

    def representation(self):
        """Get the full representation (mix basis and current)

        :rtype: str
        """

        each_diff = collections.Counter(self.diff_representation)
        each_basis = collections.Counter(self.basis.diff_representation if self.basis else '')

        r = ''

        for i in 'GNFD':
            n = 0
            if i in each_basis:
                n += each_basis[i]
            if i in each_diff:
                n += each_diff[i]

            r += i * n

        return r

    def differentiate(self, derivatives_representation, spacial_dof=None):
        """Create a new derivative from the differentiation of of the current one.

        :param derivatives_representation: the representation of the derivatives
        :type derivatives_representation: str
        :param spacial_dof: spacial degrees of freedom
        :type spacial_dof: int
        :return: a new derivative
        :rtype: Derivative
        """

        if derivatives_representation == '':
            raise ValueError(derivatives_representation)

        representation = ''
        basis_representation = self.representation()

        for i in derivatives_representation:
            if i not in ALLOWED_DERIVATIVES:
                raise RepresentationError(derivatives_representation)

            representation += i

        sdof = spacial_dof if spacial_dof else self.spacial_dof

        if 'N' in representation and not sdof:
            raise Exception('No DOF')

        return Derivative(
            from_representation=representation,
            basis=Derivative(from_representation=basis_representation, spacial_dof=sdof),
            spacial_dof=sdof
        )

    def dimension(self):
        """Return the dimension of the (full) flatten tensor

        :return: the size
        :rtype: int
        """

        size = self.basis.dimension() if self.basis else 1

        for i in self.diff_representation:
            size *= 3 if i in 'FD' else self.spacial_dof

        return size

    def shape(self):
        """Return the shape of the (full) tensor

        :return: the shape
        :rtype: list
        """

        shape = [1] if self.representation() == '' else []

        for i in self.representation():
            shape.append(3 if i in ['F', 'D'] else self.spacial_dof)

        return shape

    def order(self):
        """Get the order of derivation with respect to energy

        :rtype: int
        """

        return (self.basis.order() if self.basis else 0) + len(self.diff_representation)

    def smart_iterator(self):
        """An iterator to avoid computing all stuffs but only some by yielding a subset of independent coordinates."""

        if self.representation() == '':  # special case of energy
            yield 0
            return

        shape = self.shape()
        each = collections.Counter(self.representation())

        iterable = None

        for c in 'GNFD':

            if c not in each:
                continue

            permutable = [
                a for a in itertools.combinations_with_replacement(
                    range(3 if c in 'FD' else self.spacial_dof), each[c])
            ]

            if not iterable:
                iterable = permutable
            else:
                prev_iterable = iterable.copy()
                iterable = []
                for i in prev_iterable:
                    e = list(i)
                    for a in permutable:
                        el = e.copy()
                        el.extend(list(a))
                        iterable.append(el)

        for component in iterable:

            n = 0
            for i, e in enumerate(component):
                n += e * qcip_math.prod(shape[i + 1:])

            yield n

    def inverse_smart_iterator(self, element):
        """Back-iterate over all the other components :
        from a coordinates, give all the other ones that are equivalents

        :param element: the coordinates
        :type element: tuple|list
        """

        if self.diff_representation == '':  # special case of energy
            yield 0
            return

        each = collections.Counter(self.representation())
        shape = self.shape()

        components = []
        rest = element
        for i, (s, d) in enumerate(zip(self.shape(), self.representation())):

            current = int(math.floor(rest / qcip_math.prod(shape[i + 1:])))
            components.append(current)
            rest -= current * qcip_math.prod(shape[i + 1:])

        iterable = []
        index = 0

        for c in 'GNFD':
            if c not in each:
                continue

            n = each[c]
            permutations = [a for a in qcip_math.unique_everseen(itertools.permutations(components[index:index + n]))]
            # Since it will be called more than once, it must be a list rather than an iterator !

            if not iterable:
                iterable = list(permutations)
            else:
                prev_iterable = iterable.copy()
                iterable = []
                for i in prev_iterable:
                    e = list(i)
                    for a in permutations:
                        el = e.copy()
                        el.extend(list(a))
                        iterable.append(el)

            index += n

        for component in iterable:
            n = 0
            for i, e in enumerate(component):
                n += e * qcip_math.prod(shape[i + 1:])

            yield n
