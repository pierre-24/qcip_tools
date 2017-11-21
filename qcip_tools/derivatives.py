import math
import collections
import itertools
import numpy

from qcip_tools import math as qcip_math, numerical_differentiation

#: In term of derivative of the energy
#: ``G`` = geometrical derivative,
#: ``N`` = normal mode derivative,
#: ``F`` = static electric field derivative,
#: ``D`` = dynamic electric field derivative (which can be static),
#: ``d`` = inverse of the dynamic electric field (-w).
#: Note: do not change the order, except if you know what your are doing
ALLOWED_DERIVATIVES = ('G', 'N', 'F', 'D', 'd')

GEOMETRICAL_DERIVATIVES = ('G', 'N')
ELECTRICAL_DERIVATIVES = ('F', 'D', 'd')


def __is_derivative(derivative, typ):
    if type(derivative) is str:
        for i in typ:
            if i in derivative:
                return True
        return False
    elif type(derivative) is Derivative:
        for i in typ:
            if i in derivative.representation():
                return True
        return False
    else:
        raise TypeError(derivative)


def is_electrical(derivative):
    """Return if the derivatives contains an electrical one

    :type derivative: str|qcip_tools.derivatives.Derivative
    :rtype: bool
    """
    return __is_derivative(derivative, ELECTRICAL_DERIVATIVES)


def is_geometrical(derivative):
    """Return if the derivatives contains a geometrical one

    :type derivative: str|qcip_tools.derivatives.Derivative
    :rtype: bool
    """
    return __is_derivative(derivative, GEOMETRICAL_DERIVATIVES)

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

                if i in GEOMETRICAL_DERIVATIVES and not spacial_dof:
                    raise Exception('geometrical derivative and no spacial_dof !')

            self.diff_representation = from_representation

    def __eq__(self, other):
        if type(other) is str:
            return self == Derivative(from_representation=other, spacial_dof=self.spacial_dof)
        elif isinstance(other, Derivative):
            return self.representation() == other.representation()
        else:
            raise TypeError(other)

    def __repr__(self):
        """Get representation

        :rtype: str
        """

        s = self.representation()
        return s if len(s) != 0 else 'energy'

    def representation(self):
        """Get the full representation (mix basis and current)

        .. note::

            First comme the geometrical derivatives (G or N), then the electrical ones (F, D or d).

        :rtype: str
        """

        return self.raw_representation(exclude=ELECTRICAL_DERIVATIVES) + \
            self.raw_representation(exclude=GEOMETRICAL_DERIVATIVES)

    def raw_representation(self, exclude=None):
        """Get raw representation (simply the parent representation + obj representation).

        .. warning::

            Dangerous to use, prefer ``representation()``.

        :param exclude: exclude some element from the representation
        :type exclude: tuple|list|str
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

        if ('N' in representation and not sdof) or ('G' in representation and not sdof):
            raise Exception('No DOF')

        return Derivative(
            from_representation=representation,
            basis=Derivative(from_representation=basis_representation, spacial_dof=sdof),
            spacial_dof=sdof
        )

    def dimension(self, raw=False):
        """Return the dimension of the (full) flatten tensor

        :param raw: return the shape of the raw representation
        :type raw: bool
        :return: the size
        :rtype: int
        """

        return qcip_math.prod(self.shape(raw=raw))

    def shape(self, raw=False):
        """Return the shape of the (full) tensor

        :param raw: return the shape of the raw representation
        :type raw: bool
        :return: the shape
        :rtype: list
        """

        representation = self.raw_representation() if raw else self.representation()

        shape = [1] if representation == '' else []

        for i in representation:
            shape.append(3 if i in ELECTRICAL_DERIVATIVES else self.spacial_dof)

        return shape

    def order(self):
        """Get the order of derivation with respect to energy

        :rtype: int
        """

        return (self.basis.order() if self.basis else 0) + len(self.diff_representation)

    @classmethod
    def expend_list(cls, iterable, sub_iterable):
        """
        For each element of ``iterable``, create a new list by adding each element of ``sub_iterable``.

        :param iterable: main list
        :type iterable: list
        :param sub_iterable: sub list
        :type sub_iterable: list
        :rtype: list
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

    @classmethod
    def apply_permutation(cls, iterable, permutations):
        """Apply the permutations (list of tuple ``(from, to)``) over iterable

        :param iterable: an iterable that you can copy (a list)
        :type iterable: list
        :param permutations: list of permutation
        :type permutations: list|iterator
        :return: new iterable
        :rtype: list
        """

        n_iterable = iterable.copy()
        for a, b in permutations:
            if a == b:
                continue
            n_iterable[b] = iterable[a]

        return n_iterable

    @classmethod
    def correct_components(cls, ideal_representation, representation, list_of_components):
        """Correct the order of the components

        :param ideal_representation: representation when the different type of derivatives nicelly follows each others
        :type ideal_representation: str
        :param representation: real representation
        :type representation: str
        :param list_of_components: list of components to correct
        :type list_of_components: list
        """
        if representation != ideal_representation:
            diff = tuple(index for index, (x, y) in enumerate(zip(representation, ideal_representation)) if x != y)
            changed = None

            r_list = [a for a in representation]
            ideal_list = [a for a in ideal_representation]
            for p in itertools.permutations(diff):
                if p == diff:
                    continue
                if Derivative.apply_permutation(ideal_list, zip(diff, p)) == r_list:
                    changed = p
                    break

            if changed is None:
                raise Exception(
                    'cannot find permutation to go from {} to {} ?!?'.format(ideal_representation, representation))

            # perform permutation
            for index, components in enumerate(list_of_components):
                list_of_components[index] = Derivative.apply_permutation(components, zip(diff, changed))

    def smart_iterator(self, as_flatten=False):
        """Apply the
        `Shwarz's theorem <https://en.wikipedia.org/wiki/Symmetry_of_second_derivatives#Schwarz.27s_theorem>`_
        and only return as subset of independant coordinates. Order normally guaranteed.

        .. note::

            Ideally, the different type derivatives follows each other. This is not always the case for electrical
            ones, so there is a stage of re-ordering.

        :param as_flatten: yield full components, not flatten ones
        :type as_flatten: bool
        """

        representation = self.representation()

        if representation == '':  # special case of energy
            yield 0
            return

        each = collections.Counter(representation)
        list_of_components = []
        ideal_representation = ''

        for c in ALLOWED_DERIVATIVES:
            if c not in each:
                continue

            perms = [
                a for a in itertools.combinations_with_replacement(
                    range(self.spacial_dof if c in GEOMETRICAL_DERIVATIVES else 3), each[c])]
            list_of_components = Derivative.expend_list(list_of_components, perms)

            ideal_representation += each[c] * c

        # correct the order:
        Derivative.correct_components(ideal_representation, representation, list_of_components)

        for components in list_of_components:
            if as_flatten:
                yield self.components_to_flatten_component(components)
            else:
                yield tuple(components)

    def inverse_smart_iterator(self, element, as_flatten=False):
        """Back-iterate over all the other components :
        from a coordinates, give all the other ones that are equivalents

        :param element: the coordinates, either a flatten index (if ``as_flatten=True``) or a tuple
        :type element: int|tuple
        :param as_flatten: yield full components, not flatten ones
        :type as_flatten: bool
        """

        representation = self.representation()

        if representation == '':  # special case of energy
            yield 0
            return

        each = collections.Counter(representation)

        if as_flatten:
            components = self.flatten_component_to_components(element)
        else:
            if len(element) != self.order():
                raise ValueError('the element is not of the right size ({} != {})'.format(len(element), self.order()))
            components = element

        list_of_components = []
        ideal_representation = ''

        for c in ALLOWED_DERIVATIVES:
            if c not in each:
                continue

            components_with_this_derivative = []
            for index, cp in enumerate(representation):
                if cp == c:
                    components_with_this_derivative.append(components[index])

            permutations = [
                a for a in qcip_math.unique_everseen(itertools.permutations(components_with_this_derivative))]
            list_of_components = Derivative.expend_list(list_of_components, permutations)

            ideal_representation += each[c] * c

        # correct the order:
        Derivative.correct_components(ideal_representation, representation, list_of_components)

        for components_ in list_of_components:
            if as_flatten:
                yield self.components_to_flatten_component(components_)
            else:
                yield tuple(components_)

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
    """Use the representation to find a translation for the operator

    :param representation: the representation
    :type representation: str
    :param component: (flatten) component
    :type component: int
    :param molecule: the molecule (print the symbol of the atoms instead of a number if given)
    :type molecule: qcip_tools.molecule.Molecule
    :rtype: str
    """

    if representation not in ALLOWED_DERIVATIVES:
        raise RepresentationError(representation)

    if representation == 'N':
        return 'dQ' + ('({})'.format(component + 1) if component is not None else '')

    if representation == 'G':
        return 'dG' + ('({})'.format(
            component + 1 if not molecule else '{}{}{}'.format(
                molecule[math.floor(component / 3)].symbol, math.floor(component / 3) + 1, COORDINATES[component % 3]))
            if component is not None else '')

    if representation in ELECTRICAL_DERIVATIVES:
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

        if spacial_dof is not None and spacial_dof % 3 != 0:
            raise ValueError('DOF should be a multiple of 3 ?!?')

        self.spacial_dof = spacial_dof

        if 'D' in self.representation.representation() and frequency is None:
            raise ValueError(frequency)

        self.frequency = frequency
        self.name = name
        self.components = numpy.zeros(self.representation.shape()) if components is None else components

    def project_over_normal_modes(self, mwh):
        """Project the tensor over normal modes (get its normal mode equivalent)

        :param mwh: mass weighted hessian
        :type mwh: qcip_tools.derivatives_g.MassWeightedHessian
        :rtype: qcip_tools.derivatives.Tensor
        """

        b_repr = self.representation.representation()

        if 'N' in b_repr:
            raise ValueError('Tensor already projected ?')

        num_g = b_repr.count('G')

        if num_g == 0:
            raise ValueError('No geometrical derivative to project in {}'.format(self.representation.representation()))

        if self.spacial_dof != mwh.dof:
            raise ValueError('DOF does not matches ({} != {})'.format(self.spacial_dof, mwh.dof))

        if is_electrical(self.representation):
            nshape = [self.spacial_dof] * num_g + [3 ** (self.representation.order() - num_g)]
            projected_tensor = self.components.reshape(nshape)  # reshape
        else:
            projected_tensor = numpy.dot(self.components, mwh.displacements.transpose())
            num_g -= 1

        # project
        for i in range(num_g):
            projected_tensor = numpy.dot(mwh.displacements, projected_tensor)

        if is_electrical(self.representation):
            projected_tensor = projected_tensor.reshape(self.representation.shape())

        return Tensor(
            spacial_dof=self.spacial_dof,
            representation=b_repr.replace('G', 'N'),
            frequency=self.frequency,
            components=projected_tensor)

    def to_string(self, threshold=1e-5, columns_per_line=6, molecule=None, skip_trans_plus_rot_dof=0, **kwargs):
        """Print the tensor in a more or less textual version.

        :param threshold: show 0 instead of the value if lower than threshold
        :type threshold: float
        :param columns_per_line: number of columns per "line" of the tensor
        :type columns_per_line: int
        :param molecule: use molecule to gives the atom instead of Gxx
        :type molecule: qcip_tools.molecule.Molecule
        :param skip_trans_plus_rot_dof: skip the *n* first N (normal) modes
        :type skip_trans_plus_rot_dof: int
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

        for offset in range(skip_trans_plus_rot_dof if representation[-1] == 'N' else 0, shape[-1], columns_per_line):
            if offset != 0:
                s += '\n'

            s += (' ' * 8) * (order - 1) + ' ' * 2

            # columns title
            for index in range(0, columns_per_line):
                if (offset + index) >= shape[-1]:
                    break
                s += '{:8}'.format(representation_to_operator(representation[-1], offset + index, molecule)) + ' ' * 6

            s += '\n'

            # rows
            for idx in range(int(dimension / shape[-1])):
                components = []
                rest = idx
                can_be_skipped = False

                # row title
                for i, c in enumerate(representation[:-1]):
                    current = int(math.floor(rest / qcip_math.prod(shape[i + 1:-1])))

                    if c == 'N' and current < skip_trans_plus_rot_dof:
                        can_be_skipped = True

                    components.append(current)
                    rest -= current * qcip_math.prod(shape[i + 1:-1])

                if can_be_skipped:
                    continue

                for index, c in enumerate(components):
                    if index < len(components) - 1:
                        z = [
                            skip_trans_plus_rot_dof if representation[a] == 'N' else 0
                            for a in range(index + 1, order - 1)
                        ]
                        if components[index + 1:] == z:
                            s += '{:8}'.format(representation_to_operator(representation[index], c, molecule))
                        else:
                            s += ' ' * 8
                    else:
                        s += '{:8}'.format(representation_to_operator(representation[index], c, molecule))

                s += ' ' * 1

                # columns
                for k in range(offset, offset + 6):
                    if k >= shape[-1]:
                        break

                    val = self.components[tuple(components)][k]
                    s += '{: .6e} '.format(0.0 if math.fabs(val) < threshold else val)

                s += '\n'
        return s

    def __repr__(self):
        return self.to_string()


def compute_numerical_derivative_of_tensor(
        basis, derivative_repr, tensor_func, k_max, min_field, ratio, frequency=None, dry_run=False, accuracy_level=1,
        **kwargs):
    """
    Compute numerical derivative of tensor (taking advantage of the so-called *smart iterator*).

    .. note::

        The parameter ``accuracy_level`` is provided for retro-compatibility, but should remain to "1".

    :param basis: basis of differentiation (representation for the base tensors)
    :type basis: qcip_tools.derivatives.Derivative
    :param derivative_repr: representation of the further derivatives of the tensor
    :type derivative_repr: qcip_tools.derivatives.Derivative
    :param tensor_func: function to access to the tensor
    :param dry_run: do not fill the tensor or perform romberg analysis
    :type dry_run: bool
    :param k_max: Romberg triangle k maximum
    :type k_max: int
    :param min_field: minimal value
    :type min_field: float
    :param ratio: ratio (a)
    :type ratio: float
    :param frequency: frequency if electrical derivative
    :type frequency: float|str
    :param accuracy_level: level of accuracy. Default value (1) does not change the behavior,
      but "0" and lower use forward derivatives for order of differentiation == 3 (retro-compatibility behavior)
    :param kwargs: args passed to tensor_func
    :type kwargs: dict
    :return: tensor and romberg triangles
    :rtype: qcip_tools.derivatives.Tensor, collection.OrderedDict
    """

    b_repr = basis.representation()
    d_repr = derivative_repr.representation()

    if 'D' in d_repr:
        raise Exception('derivation with respect to D is impossible')
    if 'N' in d_repr:
        raise Exception('derivation with respect to N is impossible (but projection is!)')
    if 'F' in d_repr and 'G' in d_repr:
        raise Exception('no mixed F and G derivatives')
    if 'D' in basis.representation() and not frequency:
        raise Exception('dynamic electric field require frequency')

    derivative_order = derivative_repr.order()

    if basis.spacial_dof is None:
        basis.spacial_dof = derivative_repr.spacial_dof

    final_derivative = basis.differentiate(derivative_repr.representation())
    tensor = Tensor(final_derivative, spacial_dof=basis.spacial_dof, frequency=frequency)
    romberg_triangles = collections.OrderedDict()

    for derivative_coo in derivative_repr.smart_iterator():
        romberg_triangles[derivative_coo] = collections.OrderedDict()
        inv_derivative_coos = list(derivative_repr.inverse_smart_iterator(derivative_coo))
        each = collections.Counter(derivative_coo)
        deriv_coefs = []

        for ce in each:
            order = each[ce]
            deriv_coefs.append((numerical_differentiation.Coefficients(
                order,
                numerical_differentiation.Coefficients.choose_p_for_centered(order),
                ratio=ratio, method='C'), ce))

        if accuracy_level < 1:
            if derivative_order == 3 and derivative_coo == (0, 1, 2):
                # lower the accuracy for the ijk tensor component
                c1f = numerical_differentiation.Coefficients(1, 1, ratio=ratio, method='F')
                c1 = numerical_differentiation.Coefficients(1, 2, ratio=ratio, method='C')

                if accuracy_level == 0:
                    deriv_coefs = [(c1, 0), (c1f, 1), (c1f, 2)]
                else:  # T-REX
                    deriv_coefs = [(c1f, 0), (c1f, 1), (c1f, 2)]

            if derivative_order == 3 and len(each) == 1:  # lower the accuracy for iii tensor components
                k_max -= 1

        inverse = False
        if d_repr[0] == 'F':
            # because E(-F) = E0 - µ.F + 1/2*a.F² + ..
            #         µ(-F) = µ0 - a.F + 1/2*b.F² + ...
            if basis.order() == 0:
                inverse = True if derivative_order % 2 == 0 else False
            else:
                inverse = True if derivative_order % 2 == 1 else False

        for basis_coo in basis.smart_iterator():
            values_per_k = [
                numerical_differentiation.compute_derivative_of_function(
                    deriv_coefs,
                    tensor_func,
                    k,
                    min_field,
                    derivative_repr.shape()[0],
                    # kwargs:
                    basis=basis, inverse=inverse, component=basis_coo, frequency=frequency, **kwargs)
                for k in range(k_max)
            ]

            if not dry_run:
                trig = numerical_differentiation.RombergTriangle(values_per_k, ratio=ratio, r=2)
                romberg_triangles[derivative_coo][basis_coo] = trig
                val = trig.find_best_value()

                for inv_derivative_coo in inv_derivative_coos:
                    for inv_basis_coo in basis.inverse_smart_iterator(basis_coo):
                        if b_repr != '':
                            if 'F' in d_repr:
                                tensor.components[inv_basis_coo][inv_derivative_coo] = val[1]
                            else:
                                tensor.components[inv_derivative_coo][inv_basis_coo] = val[1]
                        else:
                            tensor.components[inv_derivative_coo] = val[1]

    if not dry_run:
        return tensor, romberg_triangles
