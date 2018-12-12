import itertools
import math
import collections
from scipy import constants

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


def _sqrt_or_neg_sqrt(i):
    try:
        return math.sqrt(i)
    except ValueError:
        return -math.sqrt(-i)


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

    def response_to_electric_field(self, field):
        """Response to an electric field, located at the (0, 0, 0) position and for which the direction is given
        by ``field``.

        .. math::

            E_i = \\sum_{ij\\ldot}^{x,y,z} \\chi_{ij\\ldots}\\,F_i\\,F_j\\,\\ldots

        :param field: The electric field
        :type field: numpy.array
        :rtype: numpy.array
        """
        if len(field) != 3:
            raise ValueError('your field should have 3 coordinates')

        response = numpy.zeros(3)

        for i in range(3):
            for set in itertools.product(range(3), repeat=self.rank() - 1):
                tmp = self.components[tuple([i] + list(set))] * field[i]
                for j in set:
                    tmp *= field[j]
                response[i] += tmp

        return response


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
        return _sqrt_or_neg_sqrt(.5 * tmp)

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

        .. math::

            \\langle\\beta^2\\rangle = \\frac{1}{105}\\,\\begin{pmatrix}
            4-26\\cos^2\\Psi+20\\cos^4\\Psi\\\\
            4+2\\cos^2\\Psi-8\\cos^4\\Psi \\\\
            1-10\\cos^2\\Psi+12\\cos^4\\Psi \\\\
            2+8\\cos^2\\Psi-4\\cos^4\\Psi \\\\
            4+2\\cos^2\\Psi-8\\cos^4\\Psi
            \\end{pmatrix}^T\\,\\begin{pmatrix}
            [g\\beta^2]_A\\\\
            [g\\beta^2]_B\\\\
            [g\\beta^2]_C\\\\
            [g\\beta^2]_D\\\\
            [g\\beta^2]_E\\\\
            \\end{pmatrix}.

        :param angle: angle (in **degree**)
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
            [4, -26, 20],
            [4, 2, -8],
            [1, -10, 12],
            [2, 8, -4],
            [4, 2, -8]
        ]), numpy.array([1, ang_2, ang_4]).transpose())

        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                    for k in derivatives.COORDINATES_LIST:
                        tmp += numpy.dot(mat_angs.transpose(), [
                            self.components[i, i, j] * self.components[j, k, k],
                            self.components[i, i, j] * self.components[k, j, k],
                            self.components[i, j, j] * self.components[i, k, k],
                            self.components[i, j, k] ** 2,
                            self.components[i, j, k] * self.components[j, i, k]
                        ])

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
        """Compute :math:`\\langle\\beta^2_{ZXX}\\rangle`:


        .. math::

            \\begin{align}
            \\langle\\beta^2_{ZXX}\\rangle
            &= \\frac{1}{105}\\,
            \\sum_{ijk}^{xyz} 6\\,\\beta_{ijk}^2
            + 3\\,\\beta_{ijj} \\beta_{ikk}
            -2\\, (\\beta_{iij} \\beta_{jkk} +\\beta_{iij} \\beta_{kjk}  +\\beta_{ijk} \\beta_{jik})\\\\
            &\\approx \\frac{1}{105}\\,\\sum_{ijk}^{xyz}4\\, \\beta_{ijk}^2-\\beta_{ijj} \\beta_{ikk}
            \\end{align}

        (actually computed using the polarization dependant intensity)

        :rtype: float
        """

        if not self.is_shg():
            raise NotSHG(self.input_fields)

        return self.polarization_angle_dependant_intensity(0)

    def beta_squared_zzz(self):
        """Compute :math:`\\langle\\beta^2_{ZZZ}\\rangle`:


        .. math::

            \\begin{align}
            \\langle\\beta^2_{ZZZ}\\rangle
            &= \\frac{1}{105}\\,\\sum_{ijk}^{xyz} 2\\,
            \\beta_{ijk}^2 + \\beta_{ijj} \\beta_{ikk}
            + 4\\, (\\beta_{iij} \\beta_{jkk}
            + \\beta_{iij} \\beta_{kjk} + \\beta_{ijk} \\beta_{jik})\\\\
            &\\approx \\frac{1}{35}\\,\\sum_{ijk}^{xyz} 2\\,\\beta_{ijk}^2 + 3 \\beta_{ijj} \\beta_{ikk}
            \\end{align}

        (actually computed using the polarization dependant intensity)

        :rtype: float
        """

        if not self.is_shg():
            raise NotSHG(self.input_fields)

        return self.polarization_angle_dependant_intensity(90)

    def beta_hrs(self):
        """Hyper-Rayleigh scattering quantity:

        .. math ::

            \\beta_{HRS}=\\sqrt{\\langle\\beta^2_{ZZZ}\\rangle + \\langle\\beta^2_{ZXX}\\rangle}

        :rtype: float
        """

        return _sqrt_or_neg_sqrt(self.beta_squared_zxx() + self.beta_squared_zzz())

    def depolarization_ratio(self):
        """Hyper-Rayleigh depolarization ratio:

        .. math::

            DR = \\frac{\\langle\\beta^2_{ZZZ}\\rangle}{\\langle\\beta^2_{XZZ}\\rangle}

        """

        return self.beta_squared_zzz() / self.beta_squared_zxx()

    def compute_sum(self, what):
        """Shorthand to simplify the computation of invariant quantities.
        ``what`` is the list of quantities to compute:
        ``beta_tensor.compute_sum([(2/5, 'ijjikk'), (3/5, 'iijkjk'])`` is equivalent to

        .. math::

            \\sum_{ijk}^{xyz} \\frac{2}{5}\\,\\sbb{ijj}{ijj} + \\frac{2}{5}\\,\\sbb{iij}{kjk}.


        :param what: what to compute ? list of tuple of the form (fraction, string of 6 char  [either i, j or k])
        :type what: list
        :rtype: float
        """

        def _sum(i, j, k):
            r = 0

            trans = {'i': i, 'j': j, 'k': k}
            for f, c in what:
                r += f * \
                    self.components[tuple(trans[x] for x in c[0:3])] * \
                    self.components[tuple(trans[x] for x in c[3:])]
            return r

        for _, c in what:
            if len(c) != 6:
                raise Exception('not 6 characters string')

        tmp = 0
        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                for k in derivatives.COORDINATES_LIST:
                    tmp += _sum(i, j, k)
        return tmp

    def spherical_J1a_contribution_squared(self):
        """
        Compute the :math:`|\\beta_{J=1\\alpha}|^2` contribution

       .. math::

           |\\beta_{J=1\\alpha}|^2 = \\frac{1}{5}\\,\\sum_{ijk}^{xyz} 2\\,\\sbb{ijj}{ikk}

       :rtype: float
       """

        if not self.is_shg():
            raise NotSHG(self.input_fields)

        return 1 / 5 * self.compute_sum([(2, 'ijjikk')])

    def spherical_J1b_contribution_squared(self):
        """
        Compute the :math:`|\\beta_{J=1\\beta}|^2` contribution

       .. math::

           |\\beta_{J=1\\beta}|^2 = \\frac{1}{5}\\,\\sum_{ijk}^{xyz} 3\\,\\sbb{iij}{kjk}

       :rtype: float
       """

        if not self.is_shg():
            raise NotSHG(self.input_fields)

        return 1 / 5 * self.compute_sum([(3, 'iijkjk')])

    def spherical_J1ab_contribution_squared(self):
        """
        Compute the :math:`|\\beta_{J=1\\alpha\\beta}|^2` contribution

       .. math::

           |\\beta_{J=1\\alpha\\beta}|^2 = -\\frac{1}{5}\\,\\sum_{ijk}^{xyz}  \\sbb{iij}{jkk}

       :rtype: float
       """

        if not self.is_shg():
            raise NotSHG(self.input_fields)

        return - 1 / 5 * self.compute_sum([(1, 'iijjkk')])

    def spherical_J1_contribution_squared(self):
        """
        Compute the :math:`|\\beta_{J=1}|^2` contribution

       .. math::

           \\begin{align}
           |\\beta_{J=1}|^2 &= |\\beta_{J=1\\alpha}|^2 + |\\beta_{J=1\\beta}|^2 + 2\\,|\\beta_{J=1\\alpha\\beta}|^2
            &= \\frac{1}{5}\\,\\sum_{ijk}^{xyz} 2\\,\\sbb{ijj}{ikk} + 3\\,\\sbb{iij}{kjk} - 2\\,\\sbb{iij}{jkk}
           \\end{align}

       :rtype: float
       """

        if not self.is_shg():
            raise NotSHG(self.input_fields)

        return 1 / 5 * self.compute_sum([(2, 'ijjikk'), (3, 'iijjkk'), (-2, 'iijkjk')])

    def spherical_J2_contribution_squared(self):
        """
        Compute the :math:`|\\beta_{J=2}|^2` contribution

       .. math::

           |\\beta_{J=2}|^2 = \\frac{1}{3}\\,\\sum_{ijk}^{xyz} 2\\beta_{ijk}^2 -2\\,\\sbb{ijk}{jik}-\\sbb{ijj}{ikk}
           -\\,\\sbb{iij}{kjk} +2\\,\\sbb{iij}{jkk}

       :rtype: float
       """

        if not self.is_shg():
            raise NotSHG(self.input_fields)

        return 1 / 3 * self.compute_sum([
            (2, 'ijkijk'), (-2, 'ijkjik'), (-1, 'ijjikk'), (2, 'iijjkk'), (-1, 'iijkjk')])

    def spherical_J3_contribution_squared(self):
        """
        Compute the :math:`|\\beta_{J=3}|^2` contribution

       .. math::

           |\\beta_{J=3}|^2 =
           \\frac{1}{15}\\,\\sum_{ijk}^{xyz} 5\\beta_{ijk}^2
           + 10\\,\\sbb{ijk}{jik}-\\sbb{ijj}{ikk} -4\\,\\sbb{iij}{kjk} -4\\,\\sbb{iij}{jkk}

       :rtype: float
       """

        if not self.is_shg():
            raise NotSHG(self.input_fields)

        return 1 / 15 * self.compute_sum([
            (5, 'ijkijk'), (10, 'ijkjik'), (-1, 'ijjikk'), (-4, 'iijjkk'), (-4, 'iijkjk')])

    def dipolar_contribution_squared(self, old_version=True):
        """Compute the square of the dipolar contribution

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
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
            tmp = self.spherical_J1_contribution_squared()

        return tmp

    def octupolar_contribution_squared(self, old_version=True):
        """Compute the square of the octupolar contribution

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
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
            tmp = self.spherical_J3_contribution_squared()

        return tmp

    def dipolar_contribution(self, old_version=True):
        """

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
        :rtype: float
        """

        return _sqrt_or_neg_sqrt(self.dipolar_contribution_squared(old_version=old_version))

    def octupolar_contribution(self, old_version=True):
        """

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
        :rtype: float
        """

        return _sqrt_or_neg_sqrt(self.octupolar_contribution_squared(old_version=old_version))

    def nonlinear_anisotropy(self, old_version=True):
        """Compute the nonlinear anisotropy:

        .. math::

            \\rho_{3/1}^2 = \\frac{|\\beta_{J=3}|^2}{|\\beta_{J=1}|^2}

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
        :rtype: float
        """

        return _sqrt_or_neg_sqrt(
            self.octupolar_contribution_squared(
                old_version=old_version) / self.dipolar_contribution_squared(old_version=old_version))

    def beta_vector(self):
        """return the hyperpolarizability vector

        .. math::

            \\beta_i = \\frac{1}{3} \\sum_j \\beta_{ijj} + \\beta_{jij} + \\beta_{jji}

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

            \\beta^K = \\frac{3}{2}\\,(\\beta_{||}-\\beta_{\\perp})

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
            r += '||B||      {: .5e}\n'.format(numpy.linalg.norm(beta_vector))

            if dipole is not None:
                r += 'B_||       {: .5e}\n'.format(self.beta_parallel(dipole))
                r += 'theta      {: .3f}°\n'.format(numpy.rad2deg(math.acos(
                    numpy.dot(dipole, beta_vector) / (numpy.linalg.norm(dipole) * numpy.linalg.norm(beta_vector)))))

                if sum(self.input_fields) == 1 or self.frequency == .0 or self.frequency == 'static':  # OEP
                    r += 'beta^K    {: .5e}\n'.format(self.beta_kerr(dipole))

            if self.is_shg():  # SHG
                B2zzz = self.beta_squared_zzz()
                B2zxx = self.beta_squared_zxx()

                r += '<B2zzz>    {: .5e}\n'.format(B2zzz)
                r += '<B2zxx>    {: .5e}\n'.format(B2zxx)
                r += 'beta_HRS   {: .5e}\n'.format(_sqrt_or_neg_sqrt(B2zxx + B2zzz))
                r += 'DR         {: .3f}\n'.format(B2zzz / B2zxx)

                # "old" version
                BJ1 = self.dipolar_contribution(old_version=True)
                BJ3 = self.octupolar_contribution(old_version=True)
                r += 'B|J=1|     {: .5e}\n'.format(BJ1)
                r += 'B|J=3|     {: .5e}\n'.format(BJ3)
                r += 'rho_3/1    {: .3f}\n'.format(BJ3 / BJ1 if BJ1 != .0 else float('inf'))

                # new version:
                r += 'B|J=1a|*   {: .5e}\n'.format(_sqrt_or_neg_sqrt(self.spherical_J1a_contribution_squared()))
                r += 'B|J=1b|*   {: .5e}\n'.format(_sqrt_or_neg_sqrt(self.spherical_J1b_contribution_squared()))
                r += 'B|J=1ab|*  {: .5e}\n'.format(_sqrt_or_neg_sqrt(self.spherical_J1ab_contribution_squared()))
                r += 'B|J=1|*    {: .5e}\n'.format(_sqrt_or_neg_sqrt(self.spherical_J1_contribution_squared()))
                r += 'B|J=2|*    {: .5e}\n'.format(_sqrt_or_neg_sqrt(self.spherical_J2_contribution_squared()))
                r += 'B|J=3|*    {: .5e}\n'.format(_sqrt_or_neg_sqrt(self.spherical_J3_contribution_squared()))

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

            \\gamma^K = \\frac{3}{2}\\,(\\gamma_{||}-\\gamma_{\\perp})

        :rtype: float
        """
        return 3 / 2 * (self.gamma_parallel() - self.gamma_perpendicular())

    def polarization_angle_dependant_intensity(self, angle):
        """Compute the angle (:math:`\\Psi`) dependant intensity in the THS setup.

        .. math::

            \\langle\\gamma^2\\rangle = \\frac{1}{630}\\,
            \\begin{pmatrix}
            6-81\\cos^2\\Psi+198\\cos^4\\Psi-126\\cos^6\\Psi\\\\
            24-108\\cos^2\\Psi-72\\cos^4\\Psi+144\\cos^6\\Psi\\\\
            12+54\\cos^2\\Psi-90\\cos^4\\Psi+18\\cos^6\\Psi\\\\
            6-54\\cos^2\\Psi+36\\cos^4\\Psi+36\\cos^6\\Psi\\\\
            4+36\\cos^2\\Psi-12\\cos^4\\Psi-12\\cos^6\\Psi \\\\
            6-81\\cos^2\\Psi+198\\cos^4\\Psi-126\\cos^6\\Psi\\\\
            12+54\\cos^2\\Psi-90\\cos^4\\Psi+18\\cos^6\\Psi
            \\end{pmatrix}^T\\,\\begin{pmatrix}
            [g\\gamma^2]_A\\\\
            [g\\gamma^2]_B\\\\
            [g\\gamma^2]_C\\\\
            [g\\gamma^2]_D\\\\
            [g\\gamma^2]_E\\\\
            [g\\gamma^2]_F\\\\
            [g\\gamma^2]_G\\\\
            \\end{pmatrix}.

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
            [6, -81, 198, -126],
            [24, -108, -72, 144],
            [12, 54, -90, 18],
            [6, -54, 36, 36],
            [4, 36, -12, -12],
            [6, -81, 198, -126],
            [12, 54, -90, 18]
        ]), numpy.array([1, ang_2, ang_4, ang_6]))

        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                    for k in derivatives.COORDINATES_LIST:
                        for l in derivatives.COORDINATES_LIST:
                            tmp += numpy.dot(mat_angs.transpose(), [
                                self.components[i, i, j, j] * self.components[k, l, l, k],
                                self.components[i, i, j, k] * self.components[j, l, l, k],
                                self.components[i, i, j, k] * self.components[l, j, l, k],
                                self.components[i, j, j, k] * self.components[i, k, l, l],
                                self.components[i, j, k, l] ** 2,
                                self.components[i, j, j, k] * self.components[k, i, l, l],
                                self.components[i, j, k, l] * self.components[j, i, k, l]
                            ])

        return 1 / 630 * tmp

    def gamma_squared_zzzz(self):
        """Compute :math:`\\langle\\gamma^2_{ZZZZ}\\rangle`:

        .. math::

            \\langle\\gamma_{ZZZZ}^2\\rangle &=
            \\frac{1}{315}\\sum_{ijkl}^{xyz} \\left\\{\\begin{array}{l}
            2\\,\\gamma_{ijkl}^2
            + 12\\,\\gamma_{iijk} \\gamma_{jllk}
            + 6\\,(\\gamma_{iijk} \\gamma_{ljlk}
            + \\gamma_{ijkl} \\gamma_{jikl} )\\\\
            + 3\\,(\\gamma_{ijjk} \\gamma_{ikll}
            +\\gamma_{iijj} \\gamma_{kllk}
            +\\gamma_{ijjk} \\gamma_{kill})
            \\end{array}\\right\\}\\\\
            &\\approx \\frac{1}{315}\\,\\sum_{ijkl}^{xyz}
            8\\,\\gamma^2_{ijkl} + 24\\, \\gamma_{iijk}\\gamma_{jkll}+3\\,\\gamma_{iijj}\\gamma_{kkll}

        (actually computed using the polarization dependant intensity)

        :rtype: float
        """

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        return self.polarization_angle_dependant_intensity(90)

    def gamma_squared_zxxx(self):
        """Compute :math:`\\langle\\gamma^2_{ZXXX}\\rangle`:

        .. math::

            \\langle\\gamma_{ZXXX}^2\\rangle &=
            \\frac{1}{630}\\sum_{ijkl}^{xyz} \\left\\{\\begin{array}{l}
            16\\,\\gamma_{ijkl}^2
            +24\\,\\gamma_{ijjk} \\gamma_{ikll}
            - 12\\,\\gamma_{iijk} \\gamma_{jllk}\\\\
             - 6\\,(\\gamma_{iijk} \\gamma_{ljlk}
             + \\gamma_{ijkl} \\gamma_{jikl} )
             -3\\,(\\gamma_{iijj} \\gamma_{kllk}+\\gamma_{ijjk} \\gamma_{kill})
            \\end{array}\\right\\}\\\\
            &\\approx \\frac{1}{630} \\sum_{ijkl}^{xyz}
            10\\,\\gamma^2_{ijkl} +3\\,\\gamma_{iijk}\\gamma_{jkll}  -3\\gamma_{iijj}\\gamma_{kkll}

        (actually computed using the polarization dependant intensity)

        :rtype: float
        """

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        return self.polarization_angle_dependant_intensity(0)

    def gamma_ths(self):
        """Compute :math:`\\gamma_{THS}`, the quantity from a third harmonic scattering experiment.

        .. math::

            \\gamma_{THS} = \\sqrt{\\langle\\gamma^2_{ZZZZ}\\rangle + \\langle\\gamma^2_{ZXXX}\\rangle}

        :rtype: float
        """

        return _sqrt_or_neg_sqrt(self.gamma_squared_zzzz() + self.gamma_squared_zxxx())

    def depolarization_ratio(self):
        """Compute the depolarization ratio:

        .. math::

            DR = \\frac{\\langle\\gamma^2_{ZZZZ}\\rangle}{\\langle\\gamma^2_{ZXXX}\\rangle}

        :rtype: float
        """

        return self.gamma_squared_zzzz() / self.gamma_squared_zxxx()

    def compute_sum(self, what):
        """Shorthand to simplify the computation of invariant quantities.
        ``what`` is the list of quantities to compute:
        ``gamma_tensor.compute_sum([(2/5, 'ijklijkl'), (3/5, 'iijjkkll'])`` is equivalent to

        .. math::

            \\sum_{ijk}^{xyz} \\frac{2}{5}\\,\\sgg{ijkl}{ijkl} + \\frac{2}{5}\\,\\sgg{iijj}{kkll}.


        :param what: what to compute ? list of tuple of the form (fraction, string of 8 char  [either i, j, k or l])
        :type what: list
        :rtype: float
        """

        def _sum(i, j, k, l):
            r = 0

            trans = {'i': i, 'j': j, 'k': k, 'l': l}
            for f, c in what:
                r += f * \
                    self.components[tuple(trans[x] for x in c[0:4])] * \
                    self.components[tuple(trans[x] for x in c[4:])]
            return r

        for _, c in what:
            if len(c) != 8:
                raise Exception('not 8 characters string')

        tmp = 0
        for i in derivatives.COORDINATES_LIST:
            for j in derivatives.COORDINATES_LIST:
                for k in derivatives.COORDINATES_LIST:
                    for l in derivatives.COORDINATES_LIST:
                        tmp += _sum(i, j, k, l)
        return tmp

    def spherical_J0_contribution_squared(self):
        """
        Compute the :math:`|\\gamma_{J=0}|^2` contribution:

        .. math::

           |\\gamma_{J=0}|^2 = \\frac{1}{5} \\sum^{xyz}_{ijkl} \\sgg{iijj}{kllk}


        :rtype: float
        """

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        return 1 / 5 * self.compute_sum([(1, 'iijjkkll')])

    def spherical_J1_contribution_squared(self):
        """
        Compute the :math:`|\\gamma_{J=1}|^2` contribution:

        .. math::

           |\\gamma_{J=1}|^2 = \\frac{1}{10} \\sum^{xyz}_{ijkl} 3\\,\\sgg{ijjk}{ikll} - 3\\,\\sgg{ijjk}{kill}


        :rtype: float
        """

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        return 1 / 10 * self.compute_sum([(3, 'ijjkikll'), (-3, 'ijjkkill')])

    def spherical_J2a_contribution_squared(self):
        """
        Compute the :math:`|\\gamma_{J=2\\alpha}|^2` contribution:

        .. math::

            |\\gamma_{J=2\\alpha}|^2 =
            \\frac{1}{42}\\,\\sum_{ijkl}^{xyz} 15\\,\\sgg{ijjk}{kill}+15\\,\\sgg{ijjk}{ikll}-10\\,\\sgg{iijj}{kllk}


        :rtype: float
        """

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        return 1 / 42 * self.compute_sum([(15, 'ijjkkill'), (15, 'ijjkikll'), (-10, 'iijjkllk')])

    def spherical_J2b_contribution_squared(self):
        """
        Compute the :math:`|\\gamma_{J=2\\beta}|^2` contribution:

        .. math::

            |\\gamma_{J=2\\beta}|^2 = \\frac{1}{21} \\sum_{ijkl}^{xyz} 15\\,\\sgg{iijk}{ljlk}-5\\,\\sgg{iijj}{kllk}


        :rtype: float
        """

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        return 1 / 21 * self.compute_sum([(15, 'iijkljlk'), (-5, 'iijjkllk')])

    def spherical_J2ab_contribution_squared(self):
        """
        Compute the :math:`|\\gamma_{J=2\\alpha\\beta}|^2` contribution:

        .. math::

            |\\gamma_{J=2\\alpha\\beta}|^2
            =\\frac{1}{21}\\,\\sum_{ijkl}^{xyz} 2\\,\\sgg{iijj}{kllk}-6\\,\\sgg{iijk}{jllk}


        :rtype: float
        """

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        return 1 / 21 * self.compute_sum([(-6, 'iijkjllk'), (2, 'iijjkllk')])

    def spherical_J2_contribution_squared(self):
        """
        Compute the :math:`|\\gamma_{J=2}|^2` contribution:

        .. math::

            \\begin{align}
            |\\gamma_{J=2}|^2 &=
            |\\gamma_{J=2\\alpha}|^2 + |\\gamma_{J=2\\beta}|^2 + 2\\,|\\gamma_{J=2\\alpha\\beta}|^2 \\\\
            &= \\frac{1}{14}\\,\\sum_{ijkl}^{xyz}
            10\\,\\sgg{iijk}{ljlk}+5\\,\\sgg{ijjk}{kill}
            +5\\,\\sgg{ijjk}{ikll}-8\\,\\sgg{iijk}{jllk}-4\\,\\sgg{iijj}{kllk}
            \\end{align}


        :rtype: float
        """

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        return 1 / 14 * self.compute_sum([
            (10, 'iijkljlk'), (5, 'ijjkkill'), (5, 'ijjkikll'), (-8, 'iijkjllk'), (-4, 'iijjkllk')])

    def spherical_J3_contribution_squared(self):
        """
        Compute the :math:`|\\gamma_{J=3}|^2` contribution:

        .. math::

            |\\gamma_{J=3}|^2 = \\frac{1}{20} \\sum^{xyz}_{ijkl}\\left\\{
            \\begin{array}{l}
            15\\,\\gamma_{ijkl}^2+20\\,\\sgg{iijk}{jllk}+\\sgg{ijjk}{kill}\\\\
            -10\\,\\sgg{iijk}{ljlk}-15\\,\\sgg{ijkl}{jikl}-11\\sgg{ijjk}{ikll}
            \\end{array}
            \\right\\}


        :rtype: float
        """

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        return 1 / 20 * self.compute_sum([
            (15, 'ijklijkl'), (20, 'iijkjllk'),
            (-10, 'iijkljlk'), (1, 'ijjkkill'), (-15, 'ijkljikl'), (-11, 'ijjkikll')])

    def spherical_J4_contribution_squared(self):
        """
        Compute the :math:`|\\gamma_{J=4}|^2` contribution:

        .. math::

            |\\gamma_{J=4}|^2 = \\frac{1}{140} \\sum^{xyz}_{ijkl}\\left\\{
            \\begin{array}{l}
            35\\,\\gamma_{ijkl}^2+105\\,\\sgg{ijkl}{jikl}+12\\,\\sgg{iijj}{kllk}\\\\
            -60\\,\\sgg{iijk}{jllk}-30\\,\\sgg{iijk}{ljlk}-15\\sgg{ijjk}{ikll}-15\\sgg{ijjk}{kill}
            \\end{array}
            \\right\\}


        :rtype: float
        """

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        return 1 / 140 * self.compute_sum([
            (35, 'ijklijkl'), (105, 'ijkljikl'),
            (12, 'iijjkllk'), (-60, 'iijkjllk'), (-30, 'iijkljlk'), (-15, 'ijjkkill'), (-15, 'ijjkikll')])

    def isotropic_contribution_squared(self, old_version=True):
        """Compute the square of the isotropic contribution

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
        :rtype: float
        """

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        tmp = 0

        if old_version:
            for i in derivatives.COORDINATES_LIST:
                for j in derivatives.COORDINATES_LIST:
                        for k in derivatives.COORDINATES_LIST:
                            for l in derivatives.COORDINATES_LIST:
                                tmp += 1 / 5 * self.components[i, i, j, j] * self.components[k, k, l, l]
        else:
            tmp = self.spherical_J0_contribution_squared()

        return tmp

    def quadrupolar_contribution_squared(self, old_version=True):
        """Compute the square of the quadrupolar contribution

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
        :rtype: float
        """

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        tmp = 0

        if old_version:
            for i in derivatives.COORDINATES_LIST:
                for j in derivatives.COORDINATES_LIST:
                        for k in derivatives.COORDINATES_LIST:
                            for l in derivatives.COORDINATES_LIST:
                                tmp += 6 / 7 * self.components[i, i, j, k] * self.components[j, k, l, l]
                                tmp -= 2 / 7 * self.components[i, i, j, j] * self.components[k, k, l, l]
        else:
            tmp = self.spherical_J2_contribution_squared()

        return tmp

    def hexadecapolar_contribution_squared(self, old_version=True):
        """Compute the square of the hexadecapolar (bless you!) contribution

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
        :rtype: float
        """

        if not self.is_thg():
            raise NotTHG(self.input_fields)

        tmp = 0

        if old_version:
            for i in derivatives.COORDINATES_LIST:
                for j in derivatives.COORDINATES_LIST:
                        for k in derivatives.COORDINATES_LIST:
                            for l in derivatives.COORDINATES_LIST:
                                tmp += self.components[i, j, k, l] ** 2
                                tmp -= 6 / 7 * self.components[i, i, j, k] * self.components[j, k, l, l]
                                tmp += 3 / 35 * self.components[i, i, j, j] * self.components[k, k, l, l]
        else:
            tmp = self.spherical_J4_contribution_squared()

        return tmp

    def isotropic_contribution(self, old_version=True):
        """Compute the isotropic contribution

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
        :rtype: float
        """

        return _sqrt_or_neg_sqrt(self.isotropic_contribution_squared(old_version=old_version))

    def quadrupolar_contribution(self, old_version=True):
        """Compute the quadrupolar contribution

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
        :rtype: float
        """
        return _sqrt_or_neg_sqrt(self.quadrupolar_contribution_squared(old_version=old_version))

    def hexadecapolar_contribution(self, old_version=True):
        """Compute the hexadecapolar contribution

        :param old_version: Use the previous (with Kleinman's conditions) version
        :type old_version: bool
        :rtype: float
        """

        return _sqrt_or_neg_sqrt(self.hexadecapolar_contribution_squared(old_version=old_version))

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

            r += 'gamma_||   {: .5e}\n'.format(para)
            r += 'gamma_per  {: .5e}\n'.format(perp)
            r += 'r (||/per) {: .3f}\n'.format(para / perp)

            if self.is_thg():  # THG
                G2zzzz = self.gamma_squared_zzzz()
                G2zxxx = self.gamma_squared_zxxx()

                r += '<G2zzzz>   {: .5e}\n'.format(G2zzzz)
                r += '<G2zxxx>   {: .5e}\n'.format(G2zxxx)
                r += 'gamma_THS  {: .5e}\n'.format(_sqrt_or_neg_sqrt(G2zzzz + G2zxxx))
                r += 'DR         {: .3f}\n'.format(G2zzzz / G2zxxx)

                # "old" definition
                GJ0 = self.isotropic_contribution(old_version=True)
                GJ2 = self.quadrupolar_contribution(old_version=True)
                GJ4 = self.hexadecapolar_contribution(old_version=True)

                r += 'G|J=0|     {: .5e}\n'.format(GJ0)
                r += 'G|J=2|     {: .5e}\n'.format(GJ2)
                r += 'G|J=4|     {: .5e}\n'.format(GJ4)
                r += 'rho_0/2    {: .5e}\n'.format(GJ0 / GJ2 if GJ2 != .0 else float('inf'))
                r += 'rho_4/2    {: .5e}\n'.format(GJ4 / GJ2 if GJ2 != .0 else float('inf'))

                # new definition
                r += 'G|J=0|*    {: .5e}\n'.format(_sqrt_or_neg_sqrt(self.spherical_J0_contribution_squared()))
                r += 'G|J=1|*    {: .5e}\n'.format(_sqrt_or_neg_sqrt(self.spherical_J1_contribution_squared()))
                r += 'G|J=2a|*   {: .5e}\n'.format(_sqrt_or_neg_sqrt(self.spherical_J2a_contribution_squared()))
                r += 'G|J=2b|*   {: .5e}\n'.format(_sqrt_or_neg_sqrt(self.spherical_J2b_contribution_squared()))
                r += 'G|J=2ab|*  {: .5e}\n'.format(_sqrt_or_neg_sqrt(self.spherical_J2ab_contribution_squared()))
                r += 'G|J=2|*    {: .5e}\n'.format(_sqrt_or_neg_sqrt(self.spherical_J2_contribution_squared()))
                r += 'G|J=3|*    {: .5e}\n'.format(_sqrt_or_neg_sqrt(self.spherical_J3_contribution_squared()))
                r += 'G|J=4|*    {: .5e}\n'.format(_sqrt_or_neg_sqrt(self.spherical_J4_contribution_squared()))

        return r
