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

        if transition_energies.representation.nstate != transition_dipoles.representation.nstate:
            raise Exception('nstates does not match')

        self.energies = transition_energies
        self.dipoles = transition_dipoles
        self.nstates = transition_energies.representation.nstate

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

        if self.energies.representation.raw_representation() == '#':
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
