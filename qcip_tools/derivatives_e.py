import itertools
import math
import collections
from scipy import constants
from contextlib import suppress

import numpy
from qcip_tools import derivatives, quantities


field_to_out = {-1: '-w', 0: '0', 1: 'w'}
in_to_field = dict((b, a) for (a, b) in field_to_out.items())
field_to_representation = {-1: 'd', 0: 'F', 1: 'D'}
representation_to_field = dict((b, a) for (a, b) in field_to_representation.items())

#: Correspondence between a name and a representation
REPRESENTATIONS = {
    'mu': 'F',
    'alpha(0;0)': 'FF',
    'beta(0;0,0)': 'FFF',
    'alpha(-w;w)': 'dD',
    'beta(-2w;w,w)': 'XDD',
    'beta(-w;w,0)': 'dDF',
    'beta(0;w,-w)': 'FDd',
    'gamma(0;0,0,0)': 'FFFF',
    'gamma(-w;0,0,w)': 'dFFD',
    'gamma(-w;w,0,0)': 'dDFF',
    'gamma(-2w;w,w,0)': 'XDDF',
    'gamma(-w;w,w,-w)': 'dDDd',
    'gamma(-3w;w,w,w)': 'XDDD',
}

DERIVATIVES = list(REPRESENTATIONS.values())

NAMES = dict((b, a) for (a, b) in REPRESENTATIONS.items())

PHENOMENON = {
    'FFF': 'static first hyperpolarizability',
    'dDF': 'electro-optical Pockels effect (EOP)',
    'FDd': 'optical rectification',
    'XDD': 'second harmonic generation',
    'FFFF': 'static second hyperpolarizability',
    'dFFD': 'dc-Kerr effect',
    'dDFF': 'dc-Kerr effect (from dalton)',
    'XDDF': '(static) electric field-induced second harmonic generation (EFISHG)',
    'dDDd': 'intensity dependent refractive index (IDRI/DFWM)',
    'XDDD': 'third harmonic generation'
}

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
    'gamma(-w;0,0,w)': 'γ(-w;0,0,w)',
    'gamma(-w;w,0,0)': 'γ(-w;0,0,w)',
    'gamma(-2w;w,w,0)': 'γ(-2w;w,w,0)',
    'gamma(-w;w,w,-w)': 'γ(-w;w,w,-w)',
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


def convert_frequency_from_string(f):
    """Convert the string value to the value in atomic units.

    :param f: the string value
    :type f: str|float
    :rtype: float
    """

    if type(f) is float:
        return f
    elif type(f) is not str:
        raise TypeError(f)

    if f == 'static':
        return 0.0

    if f[-2:].lower() == 'nm':
        freq = float(f[:-2])
        if freq != 0:
            val = constants.h * \
                constants.c / (float(f[:-2]) * 1e-9) * \
                quantities.convert(quantities.ureg.joule, quantities.ureg.hartree)
        else:
            raise Exception('A wavelength of 0nm is requested, which is impossible')
    elif f[-2:].lower() == 'ev':
        val = float(f[:-2]) * quantities.convert(quantities.ureg.electron_volt, quantities.ureg.hartree)
    elif f[-4:].lower() == 'cm-1':
        val = float(f[:-4]) * quantities.convert(quantities.ureg.wavenumber, quantities.ureg.hartree)
    else:
        val = float(f)

    return val


def convert_energy_to(e, unit):
    """Convert energy to something else

    :param e: the energy (in atomic unit)
    :type e: float
    :param unit: the unit (nm, ev, cm-1)
    :type unit: str
    """

    unit = unit.lower()

    if unit not in ['au', 'nm', 'ev', 'cm-1']:
        raise ValueError(unit)

    if unit == 'au':
        return e
    elif unit == 'ev':
        return e * quantities.convert(quantities.ureg.hartree, quantities.ureg.electron_volt)
    elif unit == 'cm-1':
        return e * quantities.convert(quantities.ureg.hartree, quantities.ureg.wavenumber)
    else:
        return (constants.h * constants.c) / \
               (e * quantities.convert(quantities.ureg.hartree, quantities.ureg.joule)) * 1e9


class BaseElectricalDerivativeTensor(derivatives.Tensor):
    """Base for all derivatives of the energy with respect to the electric field.

    """

    def __init__(self, tensor=None, input_fields=None, frequency='static'):

        self.input_fields = input_fields if input_fields is not None else []

        if input_fields:
            representation = ''.join([field_to_representation[i] for i in input_fields])
            sum_fields = -sum(input_fields)
            if sum_fields in field_to_representation:
                s_representation = field_to_representation[sum_fields]
            else:
                s_representation = 'X'
        else:
            representation = ''
            s_representation = 'F'

        if (s_representation + representation) not in DERIVATIVES:
            each = collections.Counter(representation)
            n_repr = ''

            for i in 'DdF':
                n_repr += each[i] * i

            if (s_representation + n_repr) not in DERIVATIVES:
                raise derivatives.RepresentationError(s_representation + n_repr)
            else:
                representation = n_repr

        super().__init__(s_representation + representation, frequency=frequency, components=tensor)
        self.name = self.to_name()

    def to_string(self, threshold=1e-5, **kwargs):
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

    def to_name(self, simplified=False):
        """Get the name of the tensor from its representation, by using ``REPRESENTATION``.

        :rtype: str
        """

        r = self.representation.representation()
        if r in NAMES:
            x = NAMES[r]
            if simplified:
                return SIMPLIFIED_NAMES[x]
            else:
                return x
        else:
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

    def to_string(self, threshold=1e-5, disable_extras=False, **kwargs):
        """Rewritten to add information
        """

        r = super().to_string(threshold=threshold, **kwargs)

        if not disable_extras:
            rx = r.splitlines()
            r = rx[0] + '         norm\n' + rx[1] + '{:.5e}'.format(self.norm())

        return r


class PolarisabilityTensor(BaseElectricalDerivativeTensor):
    """
    Polarisability  tensor, commonly written :math:`\\alpha(-\\omega_\\sigma;\\omega_1)`.
    """

    def __init__(self, tensor=None, input_fields=(1,), frequency='static'):
        super().__init__(tensor=tensor, input_fields=input_fields, frequency=frequency)

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
        except ValueError:
            return 0.0

    def to_string(self, threshold=1e-5, disable_extras=False, **kwargs):
        """Rewritten to add information
        """

        r = super().to_string(threshold=threshold, **kwargs)

        if not disable_extras:
            r += '\n'
            r += 'iso    {: .5e}\n'.format(self.isotropic_value())
            r += 'aniso  {: .5e}\n'.format(self.anisotropic_value())

        return r


class NotSHG(Exception):
    """Raised when trying to compute a quantity from SHG experiment on a tensor which is not"""


class FirstHyperpolarisabilityTensor(BaseElectricalDerivativeTensor):
    """
    First hyperpolarisability  tensor, commonly written :math:`\\beta(-\\omega_\\sigma;\\omega_1,\\omega_2)`.
    """

    def __init__(self, tensor=None, input_fields=(1, 1), frequency='static'):

        if len(input_fields) != 2:
            raise ValueError('There should be 2 input fields')

        super().__init__(tensor=tensor, input_fields=input_fields, frequency=frequency)

    def polarization_angle_dependant_intensity(self, angle):
        """Compute the angle (:math:`\\Psi`) dependant intensity in the SHS setup.

        :param angle: angle (in degree)
        :type angle: float
        :rtype: float
        """

        if self.input_fields != (1, 1) and self.input_fields != (0, 0) \
                and self.frequency != 'static' and self.frequency != .0:
            raise NotSHG(self.input_fields)

        tmp = 0.
        ang = math.cos(math.radians(angle))
        ang_2 = ang ** 2
        ang_4 = ang ** 4

        mat_angs = numpy.dot(numpy.array([
            [2, 8, -4],
            [4, -26, 20],
            [4, 2, -8],
            [1, -10, 12],
            [4, 2, -8]
        ]), numpy.array([1, ang_2, ang_4]).transpose())

        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                    for k in derivatives.COORDINATES_LIST:
                        tmp += numpy.dot([
                            self.components[i, j, k] ** 2,
                            self.components[i, i, j] * self.components[j, k, k],
                            self.components[i, i, j] * self.components[k, j, k],
                            self.components[i, j, j] * self.components[i, k, k],
                            self.components[i, j, k] * self.components[j, i, k]
                        ], mat_angs)

        return 1 / 105 * tmp

    def is_shg(self):
        """Check if the tensor correspond to a SHG one

        :rtype: bool
        """
        if self.input_fields != (1, 1) and self.input_fields != (0, 0) \
                and self.frequency != 'static' and self.frequency != .0:
            return False

        return True

    def beta_squared_zxx(self):
        """Compute :math:`\\langle\\beta^2_{ZXX}\\rangle`

        :rtype: float
        """

        if not self.is_shg():
            raise NotSHG(self.input_fields)

        # tmp = 1 / 9 * self.dipolar_ms_contribution_squared() + 1 / 45 * self.dipolar_fs_contribution_squared() + \
        #       1 / 15 * self.quadrupolar_contribution_squared() + 4 / 105 * self.octupolar_contribution_squared()
        # return tmp

        tmp = 0

        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                    for k in derivatives.COORDINATES_LIST:
                        tmp += -2 * self.components[i, i, j] * self.components[j, k, k]
                        tmp += -2 * self.components[i, i, j] * self.components[k, j, k]
                        tmp += 3 * self.components[i, j, j] * self.components[i, k, k]
                        tmp += 6 * self.components[i, j, k] ** 2
                        tmp += -2 * self.components[i, j, k] * self.components[j, i, k]

        return 1 / 105 * tmp

    def beta_squared_zzz(self):
        """Compute :math:`\\langle\\beta^2_{ZZZ}\\rangle`.

        :rtype: float
        """

        if not self.is_shg():
            raise NotSHG(self.input_fields)

        # tmp = 3 / 15 * self.dipolar_fs_contribution_squared() + 2 / 35 * self.octupolar_contribution_squared()
        # return tmp

        tmp = 0

        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                    for k in derivatives.COORDINATES_LIST:
                        tmp += 4 * self.components[i, i, j] * self.components[j, k, k]
                        tmp += 4 * self.components[i, i, j] * self.components[k, j, k]
                        tmp += self.components[i, j, j] * self.components[i, k, k]
                        tmp += 2 * self.components[i, j, k] ** 2
                        tmp += 4 * self.components[i, j, k] * self.components[j, i, k]

        return 1 / 105 * tmp

    def beta_hrs(self):
        """Hyper-Rayleigh scattering quantity:

        .. math ::

            \\beta_{HRS}=\\left\\{\\begin{array}{ll}
                |\\beta_{J=1}|\,\\sqrt{\\frac{2}{21}\\rho^2+\\frac{2}{9}} & \\text{if Kleinman's conditions} \\\\
                \\sqrt{\\langle\\beta^2_{ZZZ}\\rangle + \\langle\\beta^2_{ZXX}\\rangle} & \\text{otherwise}
            \\end{array}\\right.

        :param assume_kleinman:
            Assume the Kleinman conditions's (full permutation of the components) and use the alternate definition
        :type assume_kleinman: bool
        :rtype: float
        """

        return math.sqrt(self.beta_squared_zxx() + self.beta_squared_zzz())

    def depolarization_ratio(self):
        """Hyper-Rayleigh depolarization ratio:

        .. math::

            DR = \\frac{\\langle\\beta^2_{ZZZ}\\rangle}{\\langle\\beta^2_{XZZ}\\rangle}

        """

        return self.beta_squared_zzz() / self.beta_squared_zxx()

    def dipolar_fs_contribution_squared(self):
        """Calculate the square of the dipolar contribution

        :rtype: float

        .. math::

            \\begin{align}
                |\\beta_{J=1}|^2 &=\\frac{3}{5} \\sum_{ijk} \\beta_{ijj} \\beta_{ikk}
            \\end{align}

        :rtype: float
        """

        if not self.is_shg():
            raise NotSHG(self.input_fields)

        tmp = 0

        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                for k in derivatives.COORDINATES_LIST:
                    tmp += 1 / 15 * self.components[i, j, j] * self.components[i, k, k]
                    tmp += 4 / 15 * self.components[i, i, j] * self.components[j, k, k]
                    tmp += 4 / 15 * self.components[i, i, j] * self.components[k, j, k]
        return tmp

    def dipolar_ms_contribution_squared(self):
        """Calculate the square of the dipolar contribution

        :rtype: float

        .. math::

            \\begin{align}
                |\\beta_{J=1}|^2 &=\\frac{3}{5} \\sum_{ijk} \\beta_{ijj} \\beta_{ikk}
            \\end{align}

        :rtype: float
        """

        if not self.is_shg():
            raise NotSHG(self.input_fields)

        tmp = 0

        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                for k in derivatives.COORDINATES_LIST:
                    tmp += 7 / 15 * self.components[i, j, j] * self.components[i, k, k]
                    tmp -= 8 / 15 * self.components[i, i, j] * self.components[j, k, k]
                    tmp += 1 / 15 * self.components[i, i, j] * self.components[k, j, k]
        # print(tmp)
        return tmp

    def dipolar_contribution_squared(self, old_version=True):
        """Calculate the square of the dipolar contribution

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
        :rtype: float

        .. math::

            \\begin{align}
                |\\beta_{J=1}|^2 &=\\frac{3}{5} \\sum_{ijk} \\beta_{ijj} \\beta_{ikk}
            \\end{align}

        :rtype: float
        """

        if not self.is_shg():
            raise NotSHG(self.input_fields)

        tmp = 0

        if old_version:
            for i in derivatives.COORDINATES_LIST:
                tmp += 3 / 5 * self.components[i, i, i] ** 2
                for j in derivatives.COORDINATES_LIST:
                    if i != j:
                        tmp += 6 / 5 * self.components[i, i, i] * self.components[i, j, j]
                        tmp += 3 / 5 * self.components[i, j, j] ** 2
                        for k in derivatives.COORDINATES_LIST:
                            if i != k and j != k:
                                tmp += 3 / 5 * self.components[i, j, j] * self.components[i, k, k]
        else:
            tmp = self.dipolar_fs_contribution_squared() + self.dipolar_ms_contribution_squared()

        return tmp

    def octupolar_contribution_squared(self, old_version=True):
        """Calculate the square of the octupolar contribution

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
        :rtype: float

        .. math::

            \\begin{align}
                |\\beta_{J=3}|^2 &= - \\frac{1}{5} \\sum_{ijk} 3 \\beta_{ijj} \\beta_{ikk} + 5 \\beta_{ijk}^2
            \\end{align}

        :rtype: float
        """

        if not self.is_shg():
            raise NotSHG(self.input_fields)

        tmp = 0

        if old_version:
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
        else:
            for i in derivatives.COORDINATES_LIST:
                for j in derivatives.COORDINATES_LIST:
                    for k in derivatives.COORDINATES_LIST:
                        tmp -= 1 / 15 * self.components[i, j, j] * self.components[i, k, k]
                        tmp -= 4 / 15 * self.components[i, i, j] * self.components[j, k, k]
                        tmp -= 4 / 15 * self.components[i, i, j] * self.components[k, j, k]
                        tmp += 5 / 15 * self.components[i, j, k] ** 2
                        tmp += 10 / 15 * self.components[i, j, k] * self.components[j, i, k]
        return tmp

    def quadrupolar_contribution_squared(self, old_version=True):
        """Calculate the square of the quadrupolar contribution

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
        :rtype: float

        .. math::

            \\begin{align}
                |\\beta_{J=2}|^2 &= \ldots
            \\end{align}

        :rtype: float
        """

        if not self.is_shg():
            raise NotSHG(self.input_fields)

        if old_version:
            return .0

        tmp = 0
        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                for k in derivatives.COORDINATES_LIST:
                    # tmp -= 3 / 5 * self.components[i, j, j] * self.components[i, k, k]
                    # tmp += self.components[i, j, k] ** 2

                    tmp -= 1 / 3 * self.components[i, j, j] * self.components[i, k, k]
                    tmp += 2 / 3 * self.components[i, i, j] * self.components[j, k, k]
                    tmp -= 1 / 3 * self.components[i, i, j] * self.components[k, j, k]
                    tmp += 2 / 3 * self.components[i, j, k] ** 2
                    tmp -= 2 / 3 * self.components[i, j, k] * self.components[j, i, k]
        return tmp

    def dipolar_contribution(self, old_version=True):
        """

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
        :rtype: float
        """

        return math.sqrt(self.dipolar_contribution_squared(old_version=old_version))

    def quadrupolar_contribution(self, old_version=True):
        """

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
        :rtype: float
        """

        return math.sqrt(self.quadrupolar_contribution_squared(old_version=old_version))

    def octupolar_contribution(self, old_version=True):
        """

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
        :rtype: float
        """

        return math.sqrt(self.octupolar_contribution_squared(old_version=old_version))

    def nonlinear_anisotropy(self, old_version=True):
        """Compute the nonlinear anisotropy:

        .. math::

            \\rho = \\frac{|\\beta_{J=3}|}{|\\beta_{J=1}|}

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
        :rtype: float
        """

        return math.sqrt(
            self.octupolar_contribution_squared(old_version=old_version) /
            self.dipolar_contribution_squared(old_version=old_version))

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
        :type dipole: numpy.ndarray|qcip_tools.derivatives_e.ElectricDipole
        :rtype: float
        """

        if type(dipole) is ElectricDipole:
            dipole = dipole.components

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
        :type dipole: numpy.ndarray|qcip_tools.derivatives_e.ElectricDipole
        :rtype: float
        """

        if type(dipole) is ElectricDipole:
            dipole = dipole.components

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

    def to_string(self, threshold=1e-5, disable_extras=False, dipole=None, **kwargs):
        """Rewritten to add information
        """

        r = super().to_string(threshold=threshold, **kwargs)

        if not disable_extras:
            r += '\n'

            if self.representation.representation() in PHENOMENON:
                r += '({})\n\n'.format(PHENOMENON[self.representation.representation()])

            beta_vector = self.beta_vector()
            r += '||B||     {: .5e}\n'.format(numpy.linalg.norm(beta_vector))

            if dipole is not None:
                r += 'B_||      {: .5e}\n'.format(self.beta_parallel(dipole))
                r += 'theta     {: .3f}°\n'.format(numpy.rad2deg(math.acos(
                    numpy.dot(dipole, beta_vector) / (numpy.linalg.norm(dipole) * numpy.linalg.norm(beta_vector)))))

                if sum(self.input_fields) == 1 or self.frequency == .0 or self.frequency == 'static':  # OEP
                    r += 'beta^K    {: .5e}\n'.format(self.beta_kerr(dipole))

            if sum(self.input_fields) == 2 or self.frequency == .0 or self.frequency == 'static':  # SHG
                B2zzz = self.beta_squared_zzz()
                B2zxx = self.beta_squared_zxx()
                BJ1 = self.dipolar_contribution(old_version=True)
                BJ3 = self.octupolar_contribution(old_version=True)

                r += '<B2zzz>   {: .5e}\n'.format(B2zzz)
                r += '<B2zxx>   {: .5e}\n'.format(B2zxx)
                r += 'beta_HRS  {: .5e}\n'.format(math.sqrt(B2zxx + B2zzz))
                r += 'DR        {: .3f}\n'.format(B2zzz / B2zxx)
                r += 'B|J=1|    {: .5e}\n'.format(BJ1)
                r += 'B|J=3|    {: .5e}\n'.format(BJ3)

                with suppress(ValueError):
                    r += 'rho_3/1   {: .3f}\n'.format(BJ3 / BJ1 if BJ1 != .0 else float('inf'))

        return r


class NotTHG(Exception):
    """Raised when trying to compute a quantity from THG experiment on a tensor which is not"""


class SecondHyperpolarizabilityTensor(BaseElectricalDerivativeTensor):
    """
    Second hyperpolarisability  tensor, commonly written
    :math:`\\gamma(-\\omega_\\sigma;\\omega_1,\\omega_2,\\omega_3)`.
    """

    def __init__(self, tensor=None, input_fields=(1, 1, 1), frequency='static'):

        if len(input_fields) != 3:
            raise ValueError('There should be 3 input fields')

        super().__init__(tensor=tensor, input_fields=input_fields, frequency=frequency)

    def is_thg(self):
        """Check if a tensor is a THG one

        :rtype: bool
        """

        if self.input_fields != (1, 1, 1) and self.input_fields != (0, 0, 0) \
                and self.frequency != 'static' and self.frequency != .0:
            return False

        return True

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

            \\gamma_{\\perp} = \\frac{1}{15}\\,\\sum_{ij} 2\\gamma_{iijj} - \\gamma_{ijij}

        :rtype: float
        """

        val = 0.0
        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                val += 2 * self.components[i, i, j, j] - self.components[i, j, i, j]

        return val * 1 / 15

    def gamma_kerr(self):
        """Measured quantity in a dc-Kerr (-w;w,0,0) experiment

        .. math::

            \\gamma^K = \\frac{3}{2}\,(\\gamma_{||}-\\gamma_{\\perp})

        :rtype: float
        """
        return 3 / 2 * (self.gamma_parallel() - self.gamma_perpendicular())

    def polarization_angle_dependant_intensity(self, angle):
        """Compute the angle (:math:`\\Psi`) dependant intensity in the THS setup.

        :param angle: angle (in degree)
        :type angle: float
        :rtype: float
        """

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        tmp = 0.

        ang = math.cos(math.radians(angle))
        ang_2 = ang ** 2
        ang_4 = ang ** 4
        ang_6 = ang ** 6

        mat_angs = numpy.dot(numpy.array([
            [4, 36, -12, -12],
            [6, -81, 198, -126],
            [24, -108, -72, 144],
            [12, 54, -90, 18],
            [6, -54, 36, 36],
            [6, -81, 198, -126],
            [12, 54, -90, 18]
        ]), numpy.array([1, ang_2, ang_4, ang_6]))

        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                    for k in derivatives.COORDINATES_LIST:
                        for l in derivatives.COORDINATES_LIST:
                            tmp += numpy.dot([
                                self.components[i, j, k, l] ** 2,
                                self.components[i, i, j, j] * self.components[k, l, l, k],
                                self.components[i, i, j, k] * self.components[j, l, l, k],
                                self.components[i, i, j, k] * self.components[l, j, l, k],
                                self.components[i, j, j, k] * self.components[i, k, l, l],
                                self.components[i, j, j, k] * self.components[k, i, l, l],
                                self.components[i, j, k, l] * self.components[j, i, k, l]
                            ], mat_angs)

        return 1 / 630 * tmp

    def gamma_squared_zzzz(self):
        """Compute :math:`\\langle\\gamma^2_{ZZZZ}\\rangle`.

        :rtype: float
        """

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        tmp = 0

        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                    for k in derivatives.COORDINATES_LIST:
                        for l in derivatives.COORDINATES_LIST:
                            tmp += 2 * self.components[i, j, k, l] ** 2
                            tmp += 12 * self.components[i, i, j, k] * self.components[j, l, l, k]
                            tmp += 6 * self.components[i, i, j, k] * self.components[l, j, l, k]
                            tmp += 6 * self.components[i, j, k, l] * self.components[j, i, k, l]
                            tmp += 3 * self.components[i, j, j, k] * self.components[i, k, l, l]
                            tmp += 3 * self.components[i, i, j, j] * self.components[k, l, l, k]
                            tmp += 3 * self.components[i, j, j, k] * self.components[k, i, l, l]

        return 1 / 315 * tmp

    def gamma_squared_zxxx(self):
        """Compute :math:`\\langle\\gamma^2_{ZXXX}\\rangle`.

        :rtype: float
        """

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        tmp = 0

        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                    for k in derivatives.COORDINATES_LIST:
                        for l in derivatives.COORDINATES_LIST:
                            tmp += 16 * self.components[i, j, k, l] ** 2
                            tmp += 24 * self.components[i, j, j, k] * self.components[i, k, l, l]
                            tmp += -12 * self.components[i, i, j, k] * self.components[j, l, l, k]
                            tmp += -6 * self.components[i, i, j, k] * self.components[l, j, l, k]
                            tmp += -6 * self.components[i, j, k, l] * self.components[j, i, k, l]
                            tmp += -3 * self.components[i, i, j, j] * self.components[k, l, l, k]
                            tmp += -3 * self.components[i, j, j, k] * self.components[k, i, l, l]

        return 1 / 630 * tmp

    def gamma_ths(self):
        """Compute :math:`\\gamma_{THS}`, the quantity from a third harmonic scattering experiment.

        .. math::

            \\gamma_{THS} = \\sqrt{\\langle\\gamma^2_{ZZZZ}\\rangle + \\langle\\gamma^2_{ZXXX}\\rangle}

        :rtype: float
        """

        return math.sqrt(self.gamma_squared_zzzz() + self.gamma_squared_zxxx())

    def depolarization_ratio(self):
        """Compute the depolarization ratio:

        .. math::

            DR = \\frac{\\langle\\gamma^2_{ZZZZ}\\rangle}{\\langle\\gamma^2_{ZXXX}\\rangle}

        :rtype: float
        """

        return self.gamma_squared_zzzz() / self.gamma_squared_zxxx()

    def isotropic_contribution_squared(self, first_version=True):
        """Compute the square of the isotropic contribution

        :param first_version: version assuming static limit
        :type first_version: bool
        :rtype: float
        """

        if not first_version:
            raise NotImplementedError('full version!')

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        tmp = 0

        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                    for k in derivatives.COORDINATES_LIST:
                        for l in derivatives.COORDINATES_LIST:
                            tmp += 1 / 5 * self.components[i, i, j, j] * self.components[k, k, l, l]

        return tmp

    def quadrupolar_contribution_squared(self, first_version=True):
        """Compute the square of the quadrupolar contribution

        :param first_version: version assuming static limit
        :type first_version: bool
        :rtype: float
        """

        if not first_version:
            raise NotImplementedError('full version!')

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        tmp = 0

        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                    for k in derivatives.COORDINATES_LIST:
                        for l in derivatives.COORDINATES_LIST:
                            tmp += 6 / 7 * self.components[i, i, j, k] * self.components[j, k, l, l]
                            tmp -= 2 / 7 * self.components[i, i, j, j] * self.components[k, k, l, l]

        return tmp

    def hexadecapolar_contribution_squared(self, first_version=True):
        """Compute the square of the hexadecapolar (bless you!) contribution

        :param first_version: version assuming static limit
        :type first_version: bool
        :rtype: float
        """

        if not first_version:
            raise NotImplementedError('full version!')

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        tmp = 0

        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                    for k in derivatives.COORDINATES_LIST:
                        for l in derivatives.COORDINATES_LIST:
                            tmp += self.components[i, j, k, l] ** 2
                            tmp -= 6 / 7 * self.components[i, i, j, k] * self.components[j, k, l, l]
                            tmp += 3 / 35 * self.components[i, i, j, j] * self.components[k, k, l, l]

        return tmp

    def isotropic_contribution(self, first_version=True):
        """Compute the isotropic contribution

        :param first_version: version assuming static limit
        :type first_version: bool
        :rtype: float
        """

        return math.sqrt(self.isotropic_contribution_squared(first_version=first_version))

    def quadrupolar_contribution(self, first_version=True):
        """Compute the quadrupolar contribution

        :param first_version: version assuming static limit
        :type first_version: bool
        :rtype: float
        """
        return math.sqrt(self.quadrupolar_contribution_squared(first_version=first_version))

    def hexadecapolar_contribution(self, first_version=True):
        """Compute the hexadecapolar contribution

        :param first_version: version assuming static limit
        :type first_version: bool
        :rtype: float
        """

        return math.sqrt(self.hexadecapolar_contribution_squared(first_version=first_version))

    def to_string(self, threshold=1e-5, disable_extras=False, dipole=None, **kwargs):
        """Rewritten to add information
        """

        r = super().to_string(threshold=threshold, **kwargs)

        if not disable_extras:

            r += '\n'

            if self.representation.representation() in PHENOMENON:
                r += '({})\n\n'.format(PHENOMENON[self.representation.representation()])

            para = self.gamma_parallel()
            perp = self.gamma_perpendicular()

            if sum(self.input_fields) == 3 or self.frequency == .0 or self.frequency == 'static':  # THG
                G2zzzz = self.gamma_squared_zzzz()
                G2zxxx = self.gamma_squared_zxxx()

                GJ0 = self.isotropic_contribution(first_version=True)
                GJ2 = self.quadrupolar_contribution(first_version=True)
                GJ4 = self.hexadecapolar_contribution(first_version=True)

                r += '<G2zzzz>  {: .5e}\n'.format(G2zzzz)
                r += '<G2zxxx>  {: .5e}\n'.format(G2zxxx)
                r += 'gamma_THS {: .5e}\n'.format(math.sqrt(G2zzzz + G2zxxx))
                r += 'DR        {: .3f}\n'.format(G2zzzz / G2zxxx)
                r += 'DR\'      {: .3f}\n'.format(28 / 3 * (GJ0 / GJ2) ** 2 + 16 / 3)
                r += 'G|J=0|    {: .5e}\n'.format(GJ0)
                r += 'G|J=2|    {: .5e}\n'.format(GJ2)
                r += 'G|J=4|    {: .5e}\n'.format(GJ4)

                with suppress(ValueError):
                    r += 'rho_0/2   {: .5e}\n'.format(GJ0 / GJ2)
                    r += 'rho_4/2   {: .5e}\n'.format(GJ4 / GJ2)

            r += 'gamma_||  {: .5e}\n'.format(para)
            r += 'gamma_per {: .5e}\n'.format(perp)
            r += 'r         {: .3f}\n'.format(para / perp)

        return r
