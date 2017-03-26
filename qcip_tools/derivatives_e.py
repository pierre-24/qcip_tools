import itertools
import math

import numpy
from qcip_tools import derivatives

field_to_out = {-1: '-w', 0: '0', 1: 'w'}
in_to_field = {'w': 1, '-w': -1, '0': 0}

#: Correspondance between a name and a representation
REPRESENTATIONS = {
    'mu': 'F',
    'alpha': 'FF',
    'beta': 'FFF',
    'alpha(-w;w)': 'FD',
    'beta(-2w;w,w)': 'FDD',
    'beta(-w;w,0)': 'FFD',
    'gamma': 'FFFF',
    'gamma(-w;w,0,0)': 'FFFD',
    'gamma(-w;-w,w,w)': 'FDDD',
    'gamma(-2w;w,w,0)': 'FFDD',
    'gamma(-3w;w,w,w)': 'FDDD',
}

DERIVATIVES = list(REPRESENTATIONS)

#: List of all simplified names
SIMPLIFIED_NAMES = {
    'energy': 'E',
    'mu': 'µ',
    'alpha': 'α(0;0)',
    'beta': 'β(0;0,0)',
    'alpha(-w;w)': 'α(-w;w)',
    'beta(-2w;w,w)': 'β(-2w;w,w)',
    'beta(-w;w,0)': 'β(-w;w,0)',
    'gamma': 'γ(0;0,0,0)',
    'gamma(-w;w,0,0)': 'γ(-w;w,0,0)',
    'gamma(-w;-w,w,w)': 'γ(-w;-w,w,w)',
    'gamma(-2w;w,w,0)': 'γ(-2w;w,w,0)',
    'gamma(-3w;w,w,w)': 'γ(-3w;w,w,w)',
}


def fields_to_str(input_):
    """Get the string representation of a series of fields and the response, such as "-3w;w,w,w".

    :param input_: list of fields
    :type input_: tuple
    :return: representation
    :rtype: str
    """

    if type(input_) is not int:
        response = sum(input_)
    else:
        response = input_

    if response == 0:
        r = '0;'
    elif abs(response) == 1:
        r = '{}w;'.format('-' if response > 0 else '')
    else:
        r = '{}w;'.format(-response)

    r += ','.join(field_to_out[f] for f in input_)

    return r


class BaseElectricalDerivativeTensor(derivatives.Tensor):
    """Base for all derivatives of the energy with respect to the electric field.

    """

    def __init__(self, tensor=None, input_fields=None, frequency='static'):

        self.input_fields = input_fields if input_fields is not None else []

        if input_fields:
            representation = ''.join(['F' if i == 0 else 'D' for i in input_fields])
        else:
            representation = ''

        super().__init__('F' + representation, frequency=frequency, components=tensor)

    def to_string(self, threshold=1e-5, columns_per_line=6):
        """Rewritten to get a better output with this kind of tensor
        """

        s = ''
        order = self.representation.order()
        s += ' ' * (order - 1)
        s += '            x            y             z\n'
        for i in itertools.product(range(3), repeat=order - 1):
            for c in i:
                s += derivatives.COORDINATES[c]

            s += '     '

            for k in derivatives.COORDINATES_LIST:
                val = self.components[i][k]
                s += '{: .6e} '.format(0.0 if math.fabs(val) < threshold else val)
            s += '\n'
        return s

    def to_name(self):
        """Get the name of the tensor from its representation, by using ``REPRESENTATION``.

        :rtype: str
        """

        r = self.representation.representation()
        for i in DERIVATIVES:
            if REPRESENTATIONS[i] == r:
                return i

        raise KeyError(r)

    def rank(self):
        """Alias for the order

        :rtype: int
        """

        return self.representation.order()


class ElectricDipole(BaseElectricalDerivativeTensor):
    """Dipole moment, commonly written :math:`\\mu`."""

    def __init__(self, dipole=None):
        super().__init__(tensor=dipole)

    def norm(self):
        """Norm of the dipole

        :rtype: float
        """
        return numpy.linalg.norm(self.components)


class PolarisabilityTensor(BaseElectricalDerivativeTensor):
    """
    Polarisability  tensor, commonly written :math:`\\alpha(-\\omega_\\sigma;\\omega_1)`.
    """

    def __init__(self, tensor=None, frequency='static'):
        super().__init__(tensor=tensor, input_fields=(1,), frequency=frequency)

    def isotropic_value(self):
        """Isotropic value:

        .. math::

            \\bar{\\alpha}=\\frac{1}{3}\\sum_i \\alpha_{ii}.

        :rtype: float
        """
        iso = self.components.trace() * 1 / 3
        return iso

    def anisotropic_value(self):
        """Ansotropic value:

        .. math::

            \\Delta\\alpha=
            \\left[\\frac{1}{2}\\sum_{ij} 3\\,\\alpha_{ij}^2 - \\alpha_{ii}\\,\\alpha_{jj}\\right]^{1/2}.

        :rtype: float
        """
        tmp = .0
        for i in range(3):
            for j in range(3):
                tmp += 3 * self.components[i, j] ** 2 - self.components[i, i] * self.components[j, j]
        try:
            return math.sqrt(.5 * tmp)
        except:
            return 0.0


class NotSHG(Exception):
    """Raised when trying to compute a quantity from SHG experiment on a tensor which is not"""


class FirstHyperpolarisabilityTensor(BaseElectricalDerivativeTensor):
    """
    First hyperpolarisability  tensor, commonly written :math:`\\beta(-\\omega_\\sigma;\\omega_1,\\omega_2)`.
    """

    def __init__(self, tensor=None, input_fields=(1, 1), frequency='static'):

        if len(input_fields) != 2:
            raise ValueError(input_fields)

        super().__init__(tensor=tensor, input_fields=input_fields, frequency=frequency)

    def beta_squared_zxx(self):
        """Compute :math:`\\langle\\beta^2_{ZXX}\\rangle`:

        .. math::
            \\begin{align}
                \\langle\\beta_{ZXX}^2 \\rangle &= \\frac{1}{35} \\sum\\limits_{i} \\beta_{iii}^2 \\nonumber\\\\
                & + \\frac{4}{105} \\sum\\limits_{i \\neq j} \\beta_{iii} \\beta_{ijj}
                - \\frac{2}{35} \\sum\\limits_{i \\neq j} \\beta_{iii} \\beta_{jji}
                + \\frac{8}{105} \\sum\\limits_{i \\neq j} \\beta_{iij}^2 \\nonumber\\\\
                &+ \\frac{3}{35} \\sum\\limits_{i \\neq j} \\beta_{ijj}^2
                - \\frac{2}{35} \\sum\\limits_{i \\neq j} \\beta_{iij} \\beta_{jii}  \\nonumber\\\\
                &+ \\frac{1}{35} \\sum\\limits_{i \\neq j \\neq k} \\beta_{ijj} \\beta_{ikk}
                - \\frac{2}{105} \\sum\\limits_{i \\neq j \\neq k} \\beta_{iik} \\beta_{jjk}
                - \\frac{2}{105} \\sum\\limits_{i \\neq j \\neq k} \\beta_{iij} \\beta_{jkk}\\nonumber\\\\
                & + \\frac{2}{35} \\sum\\limits_{i \\neq j \\neq k} \\beta_{ijk}^2
                - \\frac{2}{105} \\sum\\limits_{i \\neq j \\neq k} \\beta_{ijk} \\beta_{jik}
            \\end{align}

        :rtype: float
        """

        if self.input_fields != (1, 1) and self.input_fields != (0, 0) \
                and self.frequency != 'static' and self.frequency != .0:
            raise NotSHG(self.input_fields)

        tmp = 0

        for i in derivatives.COORDINATES_LIST:

            tmp += 1 / 35 * self.components[i, i, i] ** 2

            for j in derivatives.COORDINATES_LIST:
                if i != j:
                    tmp += 4 / 105 * self.components[i, i, i] * self.components[i, j, j]
                    tmp -= 2 / 35 * self.components[i, i, i] * self.components[j, j, i]
                    tmp += 8 / 105 * self.components[i, i, j] ** 2
                    tmp += 3 / 35 * self.components[i, j, j] ** 2
                    tmp -= 2 / 35 * self.components[i, i, j] * self.components[j, i, i]

                    for k in derivatives.COORDINATES_LIST:
                        if i != k and j != k:
                            tmp += 1 / 35 * self.components[i, j, j] * self.components[i, k, k]
                            tmp -= 2 / 105 * self.components[i, i, k] * self.components[j, j, k]
                            tmp -= 2 / 105 * self.components[i, i, j] * self.components[j, k, k]
                            tmp += 2 / 35 * self.components[i, j, k] ** 2
                            tmp -= 2 / 105 * self.components[i, j, k] * self.components[j, i, k]

        return tmp

    def beta_squared_zzz(self):
        """Compute :math:`\\langle\\beta^2_{ZZZ}\\rangle`:

        .. math::
            \\begin{align}
                \\langle \\beta_{ZZZ}^2 \\rangle &= \\frac{1}{7} \\sum\\limits_{i} \\beta_{iii}^2 \\nonumber\\\\
                &+ \\frac{4}{35} \\sum\\limits_{i\\neq j} \\beta_{iij}^2
                + \\frac{2}{35} \\sum\\limits_{i\\neq j} \\beta_{iii}\\beta_{ijj}
                + \\frac{4}{35} \\sum\\limits_{i\\neq j} \\beta_{jii} \\beta_{iij} \\nonumber\\\\
                & + \\frac{4}{35} \\sum\\limits_{i\\neq j} \\beta_{iii} \\beta_{jji}
                + \\frac{1}{35} \\sum\\limits_{i\\neq j} \\beta_{jii}^2  \\nonumber\\\\
                &+ \\frac{4}{105} \\sum\\limits_{i\\neq j \\neq k} \\beta_{iij} \\beta_{jkk}
                + \\frac{1}{105} \\sum\\limits_{i\\neq j \\neq k} \\beta_{jii} \\beta_{jkk}
                + \\frac{4}{105} \\sum\\limits_{i\\neq j \\neq k} \\beta_{iij} \\beta_{kkj}	\\nonumber\\\\
                & + \\frac{2}{105} \\sum\\limits_{i\\neq j \\neq k} \\beta_{ijk}^2
                + \\frac{4}{105} \\sum\\limits_{i\\neq j \\neq k} \\beta_{ijk} \\beta_{jik}
            \\end{align}

        :rtype: float
        """

        if self.input_fields != (1, 1) and self.input_fields != (0, 0) \
                and self.frequency != 'static' and self.frequency != .0:
            raise NotSHG(self.input_fields)

        tmp = 0

        for i in derivatives.COORDINATES_LIST:

            tmp += 1 / 7 * self.components[i, i, i] ** 2

            for j in derivatives.COORDINATES_LIST:
                if i != j:
                    tmp += 4 / 35 * self.components[i, i, j] ** 2
                    tmp += 2 / 35 * self.components[i, i, i] * self.components[i, j, j]
                    tmp += 4 / 35 * self.components[j, i, i] * self.components[i, i, j]
                    tmp += 4 / 35 * self.components[i, i, i] * self.components[j, j, i]
                    tmp += 1 / 35 * self.components[j, i, i] ** 2

                    for k in derivatives.COORDINATES_LIST:
                        if i != k and j != k:
                            tmp += 4 / 105 * self.components[i, i, j] * self.components[j, k, k]
                            tmp += 1 / 105 * self.components[j, i, i] * self.components[j, k, k]
                            tmp += 4 / 105 * self.components[i, i, j] * self.components[k, k, j]
                            tmp += 2 / 105 * self.components[i, j, k] ** 2
                            tmp += 4 / 105 * self.components[i, j, k] * self.components[j, i, k]
        return tmp

    def beta_hrs(self, assume_kleinman=False):
        """Hyper-Rayleigh scattering quantity:

        .. math ::

            \\beta_{HRS}=\\left\\{\\begin{array}{ll}
                |\\beta_{J=1}|\,\\sqrt{\\frac{2}{21}\\rho^2+\\frac{2}{9}} & \\text{if Kleinman's conditions} \\\\
                \\sqrt{\\langle\\beta^2_{ZZZ}\\rangle + \\langle\\beta^2_{XZZ}\\rangle} & \\text{otherwise}
            \\end{array}\\right.

        :param assume_kleinman:
            Assume the Kleinman conditions's (full permutation of the components) and use the alternate definition
        :type assume_kleinman: bool
        :rtype: float
        """

        if self.input_fields != (1, 1) and self.input_fields != (0, 0) \
                and self.frequency != 'static' and self.frequency != .0:
            raise NotSHG(self.input_fields)

        if not assume_kleinman:
            return math.sqrt(self.beta_squared_zxx() + self.beta_squared_zzz())
        else:
            _BJ1_2 = self.dipolar_contribution_squared()
            _BJ3_2 = self.octupolar_contribution_squared()
            return math.sqrt(_BJ1_2 * 2 / 3 * (1 / 3 + 1 / 7 * _BJ3_2 / _BJ1_2))

    def depolarization_ratio(self):
        """Hyper-Rayleigh depolarization ratio:

        .. math::

            DR = \\frac{\\langle\\beta^2_{ZZZ}\\rangle}{\\langle\\beta^2_{XZZ}\\rangle}

        """

        if self.input_fields != (1, 1) and self.input_fields != (0, 0) \
                and self.frequency != 'static' and self.frequency != .0:
            raise NotSHG(self.input_fields)

        return self.beta_squared_zzz() / self.beta_squared_zxx()

    def dipolar_contribution_squared(self):
        """Calculate the square of the dipolar contribution

        .. math::

            \\begin{align}
                |\\beta_{J=1}|^2 &= \\frac{3}{5} \\sum\\limits_{i} \\beta_{iii}^2 \\nonumber\\\\
                &+ \\frac{3}{5} \\sum\\limits_{i\\neq j} \\beta_{iij}^2
                + \\frac{6}{5} \\sum\\limits_{i\\neq j} \\beta_{iii}\\beta_{ijj} \\nonumber\\\\
                &+ \\frac{3}{5} \\sum\\limits_{i\\neq j \\neq k} \\beta_{ijj} \\beta_{ikk}
            \\end{align}

        :rtype: float
        """

        if self.input_fields != (1, 1) and self.input_fields != (0, 0) \
                and self.frequency != 'static' and self.frequency != .0:
            raise NotSHG(self.input_fields)

        tmp = 0
        # 1 :
        for i in derivatives.COORDINATES_LIST:
            tmp += 3 / 5 * self.components[i, i, i] ** 2

            for j in derivatives.COORDINATES_LIST:
                if i != j:
                    tmp += 6 / 5 * self.components[i, i, i] * self.components[i, j, j]
                    tmp += 3 / 5 * self.components[i, j, j] ** 2

                    for k in derivatives.COORDINATES_LIST:
                        if i != k and j != k:
                            tmp += 3 / 5 * self.components[i, j, j] * self.components[i, k, k]
        return tmp

    def octupolar_contribution_squared(self):
        """Calculate the square of the octupolar contribution

        .. math::

            \\begin{align}
                |\\beta_{J=3}|^2 &= \\frac{2}{5} \\sum\\limits_{i} \\beta_{iii}^2 \\nonumber\\\\
                &+ \\frac{12}{5} \\sum\\limits_{i\\neq j} \\beta_{iij}^2
                - \\frac{6}{5} \\sum\\limits_{i\\neq j} \\beta_{iii}\\beta_{ijj} \\nonumber\\\\
                &- \\frac{3}{5} \\sum\\limits_{i\\neq j \\neq k} \\beta_{ijj} \\beta_{ikk}
                + \\sum\\limits_{i\\neq j \\neq k} \\beta_{ijk}^2
            \\end{align}

        :rtype: float
        """

        if self.input_fields != (1, 1) and self.input_fields != (0, 0) \
                and self.frequency != 'static' and self.frequency != .0:
            raise NotSHG(self.input_fields)

        tmp = 0
        # 1 :
        for i in derivatives.COORDINATES_LIST:
            tmp += 2 / 5 * self.components[i, i, i] ** 2

            for j in derivatives.COORDINATES_LIST:
                if i != j:
                    tmp -= 6 / 5 * self.components[i, i, i] * self.components[i, j, j]
                    tmp += 12 / 5 * self.components[i, j, j] ** 2

                    for k in derivatives.COORDINATES_LIST:
                        if i != k and j != k:
                            tmp -= 3 / 5 * self.components[i, j, j] * self.components[i, k, k]
                            tmp += self.components[i, j, k] ** 2
        return tmp

    def dipolar_contribution(self):
        """

        :rtype: float
        """

        if self.input_fields != (1, 1) and self.input_fields != (0, 0) \
                and self.frequency != 'static' and self.frequency != .0:
            raise NotSHG(self.input_fields)

        return math.sqrt(self.dipolar_contribution_squared())

    def octupolar_contribution(self):
        """

        :rtype: float
        """

        if self.input_fields != (1, 1) and self.input_fields != (0, 0) \
                and self.frequency != 'static' and self.frequency != .0:
            raise NotSHG(self.input_fields)

        return math.sqrt(self.octupolar_contribution_squared())

    def nonlinear_anisotropy(self):
        """Compute the nonlinear anisotropy:

        .. math::

            \\rho = \\frac{|\\beta_{J=3}|}{|\\beta_{J=1}|}

        :rtype: float
        """

        if self.input_fields != (1, 1) and self.input_fields != (0, 0) \
                and self.frequency != 'static' and self.frequency != .0:
            raise NotSHG(self.input_fields)

        return math.sqrt(self.octupolar_contribution_squared() / self.dipolar_contribution_squared())

    def beta_vector(self):
        """return the hyperpolarizability vector

        .. math::

            \\beta_i = \\frac{1}{3} \\sum_j \\beta_{ijj} + \\beta_{jij} + \\beta_{jii}

        :rtype: numpy.ndarray
        """

        vec = numpy.zeros(3)
        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                vec[i] += 1 / 3 * (self.components[i, j, j] + self.components[j, i, j] + self.components[j, j, i])
        return vec

    def beta_parallel(self, dipole):
        """Fetch beta value under parallel polarization of all fields:

        .. math::

            \\begin{align}
            \\beta_{||} &= \\frac{1}{5}\\,\\sum_i
            \\frac{\\mu_i}{||\\vec{\\mu}||}\\,\\sum_{j} \\beta_{ijj}
            + \\beta_{jij} + \\beta_{jji} \\\\
            &= \\frac{3}{5}\\,\\sum_i\\frac{\\mu_i\\,\\beta_{i}}{||\\vec{\\mu}||}
            \\end{align}

        where :math:`\\beta_{i}` is a component of the beta "tensor".

        :param dipole: dipole moment of the molecule
        :type dipole: numpy.ndarray
        :rtype: float
        """

        norm = numpy.linalg.norm(dipole)
        val = 0.0
        for i in derivatives.COORDINATES_LIST:
            val_ = 0.0
            for j in derivatives.COORDINATES_LIST:
                val_ += self.components[i, j, j] + self.components[j, i, j] + self.components[j, j, i]
            val += dipole[i] * val_

        return val * 1 / 5 / norm

    def beta_perpendicular(self, dipole):
        """Fetch beta value under perpendicular polarization of optical and static field.

        .. math::

            \\beta_{\\perp} = \\frac{1}{5}\\,\\sum_i
            \\frac{\\mu_i}{||\\vec{\\mu}||}\\,\\sum_{j} 2\\beta_{ijj} - 3\\beta_{jij} + 2 \\beta_{jji}.


        :param dipole: the dipole moment
        :type dipole: numpy.ndarray
        :rtype: float
        """

        norm = numpy.linalg.norm(dipole)
        val = 0.0
        for i in derivatives.COORDINATES_LIST:
            val_ = 0.0
            for j in derivatives.COORDINATES_LIST:
                val_ += 2 * self.components[i, j, j] - 3 * self.components[j, i, j] + 2 * self.components[j, j, i]
            val += dipole[i] * val_

        return val * 1 / 5 / norm

    def beta_kerr(self, dipole):
        """Measured quantity in a dc-Kerr :math:`\\beta(-\\omega;\\omega,0)` experiment.

        .. math::

            \\beta^K = \\frac{3}{2}\,(\\beta_{||}-\\beta_{\\perp})

        :param dipole: dipole moment
        :type dipole: numpy.ndarray
        :return: beta dc-Kerr
        :rtype: float
        """

        return 3 / 2 * (self.beta_parallel(dipole) - self.beta_perpendicular(dipole))


class SecondHyperpolarizabilityTensor(BaseElectricalDerivativeTensor):
    """
    Second hyperpolarisability  tensor, commonly written
    :math:`\\gamma(-\\omega_\\sigma;\\omega_1,\\omega_2,\\omega_3)`.
    """

    def __init__(self, tensor=None, input_fields=(1, 1, 1), frequency='static'):

        if len(input_fields) != 3:
            raise ValueError(input_fields)

        super().__init__(tensor=tensor, input_fields=input_fields, frequency=frequency)

    def gamma_parallel(self):
        """Fetch parallel isotropic average of gamma

        .. math::

            \\gamma_{||} = \\frac{1}{15}\\,\\sum_{ij} \\gamma_{iijj} + \\gamma_{ijji} + \\gamma_{ijij}

        :rtype: float
        """

        val = 0.0
        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                val += self.components[i, i, j, j] + self.components[i, j, j, i] + self.components[i, j, i, j]

        return val * 1 / 15

    def gamma_perpendicular(self):
        """Fetch perpendicular isotropic average of gamma

        .. math::

            \\gamma_{\\perp} = \\frac{2}{15}\\,\\sum_{ij} \\gamma_{ijji} - \\gamma_{iijj}

        :rtype: float
        """

        val = 0.0
        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                val += 2 * self.components[i, j, j, i] - self.components[i, i, j, j]

        return val * 1 / 15

    def gamma_kerr(self):
        """Measured quantity in a dc-Kerr (-w;w,0,0) experiment

        .. math::

            \\gamma^K = \\frac{3}{2}\,(\\gamma_{||}-\\gamma_{\\perp})

        :rtype: float
        """
        return 3 / 2 * (self.gamma_parallel() - self.gamma_perpendicular())
