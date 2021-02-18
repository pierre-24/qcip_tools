import numpy

from qcip_tools import derivatives


class Excitations:
    """

    :param energies: transition (or excited state) energies (``!`` or ``#``)
    :type energies: qcip_tools.derivatives.Tensor
    :param transition_dipoles: transition dipoles (``!F``)
    :type transition_dipoles: qcip_tools.derivatives.Tensor
    """

    def __init__(self, transition_energies, transition_dipoles):
        if type(transition_energies) is not derivatives.Tensor or type(transition_dipoles) is not derivatives.Tensor:
            raise TypeError('must be Tensors')

        if transition_energies.representation.representation() not in ['!', '#']:
            raise Exception('this is not transition energies')

        if transition_dipoles.representation.representation() != '!F':
            raise Exception('this is not transition dipole')

        if transition_energies.representation.nstates != transition_dipoles.representation.nstates:
            raise Exception('nstates does not match')

        self.energies = transition_energies
        self.dipoles = transition_dipoles
        self.nstates = transition_energies.representation.nstates

    def oscillator_strength(self, i):
        """Get the oscillator strength

        .. math::

            f_{gi} = \\frac{2}{3}\\,E_{gi}\\,|\\mu_{gi}|^2

        where :math:`i` is the index of the excited state.

        :param i: index (``0 < i < nstate``)
        :type i: int
        :rtype: float
        """
        e = self.transition_energy(i)
        ds = numpy.linalg.norm(self.transition_dipole(i))**2

        return 2. / 3. * e * ds

    def transition_energy(self, i):
        """Get the transition energy

        :param i: index (``0 < i < nstate``)
        :type i: int
        :rtype: float
        """
        if i <= 0 or i >= self.nstates:
            raise ValueError(i)

        if self.energies.representation.raw_representation() == '!':
            return self.energies.components[i]
        else:
            return self.energies.components[i] - self.energies.components[0]

    def transition_dipole(self, i):
        """Get the transition dipole

        :param i: index (``0 < i < nstate``)
        :type i: int
        :rtype: numpy.ndarray
        """
        if i <= 0 or i >= self.nstates:
            raise ValueError(i)

        return self.dipoles.components[i]


class Configuration:
    """Define a determinant as an excitation of the ground state determinant
    (with all electrons stored in the lowest energy MOs)

    ``excitation`` is a list of tuple containing two integer and one string:
    source MO wrt HOMO, destination MO wrt LUMO, spin. Thus,

    + ``(0,0,None)`` is a HOMO→LUMO excitation of spin alpha,
    + ``(-1,2,'a')`` is a HOMO-1→LUMO+2 excitation of spin alpha, etc.
    """
    def __init__(self, excitations):
        self.excitations = excitations

    def n(self):
        """Return the number of excitations

        :rtype: int
        """
        return len(self.excitations)

    def __repr__(self):
        r = ','.join('{s}H{f}→{s}L{t}'.format(**{
            's': s if s is not None else '',
            'f': '' if f == 0 else '{:+d}'.format(f),
            't': '' if t == 0 else '{:+d}'.format(t),
        }) for f, t, s in self.excitations)

        return '({})'.format(r) if self.n() > 1 else r


class ConfigurationStateFunction:
    """(Symmetry-adapted) linear combination of Slater determinants (configuration).

    ``configuration`` is a list of tuple of one float (weight) and one ``Configuration``.
    """

    def __init__(self, configurations):
        self.configurations = sorted(configurations, key=lambda x: abs(x[0]), reverse=True)  # important configs first

    def missing(self):
        """Return the percentage of missing configurations
        """

        return 1 - sum(c[0] ** 2 for c in self.configurations)

    def to_string(self, limit=.01):
        """Text representation, only configurations above ``limit``
        """
        if len(self.configurations) == 0:
            return 'GS'

        return ';'.join(
            '{:.1f}% of {}'.format(100 * coef ** 2, config)
            for coef, config in self.configurations if coef ** 2 >= limit)

    def __repr__(self):
        return self.to_string()
