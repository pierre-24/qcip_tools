import math
import collections
import itertools
import numpy

from qcip_tools import math as qcip_math

#: In term of derivative of the energy
#: ``F`` = static electric field derivative, ``D`` = dynamic electric field derivative (which can be static),
#: ``d`` = inverse of the dynamic electric field (-w)
#: ``G`` = geometrical derivative, ``N`` = normal mode derivative
#: Note: do not change the order, except if you have good reasons to do so.
ALLOWED_DERIVATIVES = ('G', 'N', 'F', 'D', 'd')

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
                        raise ValueError('basis and current does not have the same dof'.format(
                            spacial_dof, basis.spacial_dof))
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

    def __eq__(self, other):
        if type(other) is str:
            return self == Derivative(from_representation=other, spacial_dof=self.spacial_dof)
        elif isinstance(other, Derivative):
            return self.representation() == other.representation()
        else:
            raise TypeError(other)

    def representation(self):
        """Get the full representation (mix basis and current)

        :rtype: str
        """

        each_ = collections.Counter(self.raw_representation(exclude='FDd'))

        r = ''

        for i in 'GN':
            if i not in each_:
                continue

            r += i * each_[i]

        return r + self.raw_representation(exclude=['G', 'N'])

    def raw_representation(self, exclude=None):
        """Get raw representation (simply the parent representation + obj representation).

        .. warning::

            Dangerous to use, prefer ``representation()``.

        :param exclude: exclude some element from the representation
        :type exclude: list|str
        :rtype: str
        """

        raw = (self.basis.diff_representation if self.basis else '') + self.diff_representation
        if exclude is None:
            return raw

        return ''.join(a for a in raw if a not in exclude)

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
            size *= 3 if i in 'FDd' else self.spacial_dof

        return size

    def shape(self):
        """Return the shape of the (full) tensor

        :return: the shape
        :rtype: list
        """

        shape = [1] if self.representation() == '' else []

        for i in self.representation():
            shape.append(3 if i in 'FDd' else self.spacial_dof)

        return shape

    def order(self):
        """Get the order of derivation with respect to energy

        :rtype: int
        """

        return (self.basis.order() if self.basis else 0) + len(self.diff_representation)

    @classmethod
    def expend_list(cls, iterable, sub_iterable):
        """For each element of ``iterable``, create a new list by adding each element of ``sub_iterable``.

        :param iterable: main list
        :type iterable: list
        :param sub_iterable: sub list
        :type sub_iterable: list
        :rtype list
        """

        if not iterable:
            return sub_iterable
        if not sub_iterable:
            return iterable
        else:
            prev_iterable = iterable.copy()
            iterable = []
            for i in prev_iterable:
                e = list(i)
                for a in sub_iterable:
                    el = e.copy()
                    el.extend(list(a))
                    iterable.append(el)

        return iterable

    def smart_iterator(self):
        """Apply the
        `Shwarz's theorem <https://en.wikipedia.org/wiki/Symmetry_of_second_derivatives#Schwarz.27s_theorem>`_
        and only return as subset of independant coordinates. Order normally guaranteed.
        """

        if self.representation() == '':  # special case of energy
            yield 0
            return

        representation = self.representation()
        each = collections.Counter(representation)

        iterable = []

        for c in ALLOWED_DERIVATIVES:
            if c not in each:
                continue

            permutable = [a for a in itertools.combinations_with_replacement(
                range(self.spacial_dof if c in ['G', 'N'] else 3), each[c])]

            iterable = Derivative.expend_list(iterable, permutable)

        for component in iterable:
            yield self.components_to_flatten_component(component)

    def inverse_smart_iterator(self, element):
        """Back-iterate over all the other components :
        from a coordinates, give all the other ones that are equivalents

        :param element: the coordinates
        :type element: int
        """

        if self.representation() == '':  # special case of energy
            yield 0
            return

        representation = self.representation()
        each = collections.Counter(representation)
        components = self.flatten_component_to_components(element)

        iterable = []
        index = 0

        for c in ALLOWED_DERIVATIVES:
            if c not in each:
                continue

            n = each[c]
            permutations = [a for a in qcip_math.unique_everseen(itertools.permutations(components[index:index + n]))]
            iterable = Derivative.expend_list(iterable, permutations)
            index += n

        for component in iterable:
            yield self.components_to_flatten_component(component)

    def flatten_component_to_components(self, i):
        """From the index if the array would be flatten, gives the components

        :param i: flatten index
        :type i: int
        :rtype: list
        """

        components = []
        rest = i
        shape = self.shape()

        for i, _ in enumerate(self.representation()):
            current = int(math.floor(rest / qcip_math.prod(shape[i + 1:])))
            components.append(current)
            rest -= current * qcip_math.prod(shape[i + 1:])

        return components

    def components_to_flatten_component(self, components):
        """Given a component, give the index if the array would be flatten

        :param components: components
        :type components: list|tuple
        :rtype: int
        """

        shape = self.shape()

        n = 0
        for i, e in enumerate(components):
            n += e * qcip_math.prod(shape[i + 1:])

        return n


def representation_to_operator(representation, component=None, molecule=None):
    if representation not in ALLOWED_DERIVATIVES:
        raise RepresentationError(representation)

    if representation == 'N':
        return 'dQ' + ('({})'.format(component + 1) if component is not None else '')

    if representation == 'G':
        return 'dG' + ('({})'.format(
            component + 1 if not molecule else '{}{}'.format(
                molecule[math.floor(component / 3)].symbol, COORDINATES[component % 3]))
            if component is not None else '')

    if representation in 'FDd':
        return 'dF' + ('({})'.format(COORDINATES[component]) if component is not None else '')


class Tensor:
    """Create a tensor to store a given derivative.

    :param representation: representation of the derivative
    :type representation: str|Derivative
    :param spacial_dof: number of spacial degrees of freedom
    :type spacial_dof: int
    :param frequency: frequency if dynamic electric field
    :type frequency: str|float
    """

    def __init__(self, representation, components=None, spacial_dof=None, frequency=None, name=''):
        if isinstance(representation, Derivative):
            self.representation = representation
        else:
            self.representation = Derivative(from_representation=representation, spacial_dof=spacial_dof)

        if components is not None:
            if type(components) is list:
                components = numpy.array(components)

            if components.shape != self.representation.shape():
                components = components.reshape(self.representation.shape())

        self.spacial_dof = spacial_dof

        if 'D' in self.representation.representation() and frequency is None:
            raise ValueError(frequency)

        self.frequency = frequency
        self.name = name
        self.components = numpy.zeros(self.representation.shape()) if components is None else components

    def to_string(self, threshold=1e-5, columns_per_line=6, molecule=None, **kwargs):
        """Print the tensor in a more or less textual version.

        :param threshold: show 0 instead of the value if lower than threshold
        :type threshold: float
        :param columns_per_line: number of columns per "line" of the tensor
        :type columns_per_line: int
        :param molecule: use molecule to gives the atom insteaad of Gxx
        :type molecule: qcip_tools.molecule.Molecule
        :return: representation
        :rtype: str
        """
        s = ''
        order = self.representation.order()
        shape = self.representation.shape()
        representation = self.representation.representation()
        dimension = self.representation.dimension()

        if self.name:
            s += self.name + '\n'

        for offset in range(0, shape[-1], columns_per_line):
            if offset != 0:
                s += '\n'

            s += (' ' * 8) * (order - 1) + ' ' * 2

            for index in range(0, columns_per_line):
                if (offset + index) >= shape[-1]:
                    break
                s += '{:8}'.format(representation_to_operator(representation[-1], offset + index, molecule)) + ' ' * 6

            s += '\n'

            for idx in range(int(dimension / shape[-1])):

                components = []
                rest = idx
                for i, _ in enumerate(representation[:-1]):
                    current = int(math.floor(rest / qcip_math.prod(shape[i + 1:-1])))
                    components.append(current)
                    rest -= current * qcip_math.prod(shape[i + 1:-1])

                for index, c in enumerate(components):
                    if index < len(components) - 1:
                        if numpy.all(numpy.array(components[index + 1:]) == 0):
                            s += '{:8}'.format(representation_to_operator(representation[index], c, molecule))
                        else:
                            s += '        '
                    else:
                        s += '{:8}'.format(representation_to_operator(representation[index], c, molecule))

                s += ' ' * 1

                for k in range(offset, offset + 6):
                    if k >= shape[-1]:
                        break

                    val = self.components[tuple(components)][k]
                    s += '{: .6e} '.format(0.0 if math.fabs(val) < threshold else val)

                s += '\n'
        return s

    def __repr__(self):
        return self.to_string()
