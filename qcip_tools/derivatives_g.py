import numpy
import math
from scipy import constants

from qcip_tools import derivatives, quantities, math as qcip_math

#: Correspondence between a name and a representation
REPRESENTATIONS = {
    'gradient': 'G',
    'hessian': 'GG',
    'cubic_FF': 'GGG',
    'projected gradient': 'N',
    'projected hessian': 'NN'
}

DERIVATIVES = list(REPRESENTATIONS.values())

NAMES = dict((b, a) for (a, b) in REPRESENTATIONS.items())

#: List of all simplified names
SIMPLIFIED_NAMES = {
    'energy': 'E',
    'gradient': 'Fi',
    'hessian': 'Hij',
    'cubic_FF': 'Fijk',
    'projected gradient': 'Fi(N)',
    'projected hessian': 'Hij(N)'
}


class BaseGeometricalDerivativeTensor(derivatives.Tensor):
    """Base for all derivatives of the energy with respect to geometrical parameters
    """

    def __init__(self, spacial_dof, trans_plus_rot=0, representation='', components=None):

        if 'F' in representation or 'D' in representation:
            raise derivatives.RepresentationError(representation)

        self.trans_plus_rot = trans_plus_rot

        super().__init__(representation, components=components, spacial_dof=spacial_dof)

    def to_name(self):
        """Get the name of the tensor from its representation, by using ``REPRESENTATION``.

        :rtype: str
        """

        r = self.representation.representation()
        for i in DERIVATIVES:
            if REPRESENTATIONS[i] == r:
                return i

        raise KeyError(r)


class MassWeightedHessian:
    """To be exact, this is not a derivatives of the energy. So it is not a child of ``derivatives.Tensor`` !

    :param molecule: the molecule
    :type molecule: qcip_tools.molecule.Molecule
    :param cartesian_hessian: cartesian hessian (``GG``).
    :type cartesian_hessian: numpy.array|qcip_tools.derivatives_g.BaseGeometricalDerivativeTensor
    :param scale: scaling factor for the vibrational frequencies. It is not applied directly on the frequencies, but
      but used, for example, in thermochemistry
    :type scale: float
    """

    #: Convert inertia from UMA*A² to kg*m²
    INERTIA_CONVERSION = quantities.convert(
        quantities.ureg.atomic_mass_constant * quantities.ureg.angstrom ** 2,
        quantities.ureg.kilogram * quantities.ureg.meter ** 2)
    #: Convert mass from AMU to kg
    MASS_CONVERSION = quantities.convert(
        quantities.ureg.atomic_mass_constant, quantities.ureg.kilogram)
    #: Convert hartree to hertz
    VIB_ENERGY_CONVERSION = quantities.convert(
        quantities.ureg.hartree / quantities.ureg.planck_constant, quantities.ureg.hertz)
    #: Convert energy from J/mol to Hartree (per particle)
    ENERGY_IN_AU_CONVERSION = quantities.convert(quantities.ureg.joule, quantities.ureg.hartree) / constants.Avogadro

    #: Convert Hartree to wavenumber (cm⁻¹)
    HARTREE_TO_WAVENUMBER_CONVERSION = quantities.convert(quantities.ureg.hartree, quantities.ureg.wavenumber)

    def __init__(self, molecule, cartesian_hessian=None, scale=1.0):

        self.molecule = molecule
        self.linear = self.molecule.linear()
        self.scale = scale

        self.dof = len(molecule) * 3
        self.trans_plus_rot_dof = 5 if self.linear else 6
        self.vibrational_dof = self.dof - self.trans_plus_rot_dof

        weights = [a.mass for a in self.molecule]
        self.mass_matrix = numpy.ones((self.dof, self.dof))

        for i in range(len(self.molecule)):
            self.mass_matrix[i * 3:(i + 1) * 3, :] *= 1 / math.sqrt(weights[i])
            self.mass_matrix[:, i * 3:(i + 1) * 3] *= 1 / math.sqrt(weights[i])

        self.cartesian_hessian = numpy.zeros((self.dof, self.dof))
        self.components = numpy.zeros((self.dof, self.dof))
        self.frequencies = numpy.zeros(self.dof)
        self.displacements = numpy.zeros((self.dof, self.dof))
        self.reduced_masses = numpy.zeros(self.dof)
        self.normal_modes = numpy.zeros((self.dof, self.dof))

        self.included_modes = list(range(self.trans_plus_rot_dof, self.dof))

        if cartesian_hessian is not None:
            if issubclass(type(cartesian_hessian), BaseGeometricalDerivativeTensor):
                self.from_cartesian_hessian(cartesian_hessian.components)
            if issubclass(type(cartesian_hessian), derivatives.Tensor) and cartesian_hessian.representation == 'GG':
                self.from_cartesian_hessian(cartesian_hessian.components)
            elif type(cartesian_hessian) is numpy.ndarray:
                self.from_cartesian_hessian(cartesian_hessian)
            else:
                raise TypeError(cartesian_hessian)

    def from_cartesian_hessian(self, cartesian_hessian):
        """From a given carthesian hessian, create mass-weigted ones, and extracts frequencies and displacements.

        :param cartesian_hessian: the cartesian hessian
        :type cartesian_hessian: numpy.ndarray
        """

        if cartesian_hessian.shape != self.components.shape:
            raise ValueError('Carthesian hessian dimension {} != {}'.format(
                cartesian_hessian.shape, self.components.shape))

        self.cartesian_hessian = cartesian_hessian
        self.components = self.mass_matrix * cartesian_hessian

        if not numpy.allclose(self.components.transpose(), self.components):
            for x in range(self.dof):
                for y in range(0, x):
                    self.components[y, x] = self.components[x, y]

        # normal coordinates analysis:
        eigvals, eigvecs = numpy.linalg.eig(self.components)
        weigths = numpy.zeros(self.dof)

        for i in range(len(self.molecule)):
            weigths[i * 3:(i + 1) * 3] = math.sqrt(self.molecule[i].mass * quantities.AMUToElectronMass)

        no_freqs = [0] * self.dof
        no_disps = numpy.zeros((self.dof, self.dof))

        for i in range(self.dof):
            val = eigvals[i] / quantities.AMUToElectronMass
            if val > 0:
                no_freqs[i] = math.sqrt(val)
            else:
                no_freqs[i] = -math.sqrt(-val)

            self.normal_modes[i] = eigvecs[:, i].copy()
            no_disps[i] = eigvecs[:, i] / weigths

        # order freqs and displacements according to freqs
        ordering = [a for a in range(self.dof)]
        self.frequencies, ordering = (list(t) for t in zip(*sorted(zip(no_freqs, ordering))))
        for idx, pos in enumerate(ordering):
            self.displacements[idx] = no_disps[pos]

        self.reduced_masses = numpy.sum(self.displacements**2, axis=1)**-1

    def output_displacements(self):
        """Gaussian-like formatting of displacements. Frequencies in :math:`cm^{-1}`, reduced masses in *AMU*,
        displacement in :math:`\\text{Å}\\,AMU^{-1/2}`.

        :rtype: str
        """

        r = ''
        max_modes = len(self.included_modes)

        for offset in range(0, max_modes, 5):

            line_mode = ' ' * 10
            line_freqs = 'Frequencies   '
            line_rmass = 'Reduced-masses'

            for mode_i in range(offset, offset + 5):
                if mode_i >= max_modes:
                    break

                mode = self.included_modes[mode_i]

                line_mode += '{:15}'.format(mode_i + 1)
                line_freqs += '{: 15.4f}'.format(self.frequencies[mode] * self.HARTREE_TO_WAVENUMBER_CONVERSION)
                line_rmass += '{: 15.5f}'.format(self.reduced_masses[mode] / quantities.AMUToElectronMass)

            r += line_mode
            r += '\n'
            r += line_freqs
            r += '\n'
            r += line_rmass
            r += '\n'
            r += 'Coord  AN Element:\n'

            for index, a in enumerate(self.molecule):
                for c, cx in enumerate(derivatives.COORDINATES_LIST):
                    r += '{:4} {:4} {:^6}'.format(c + 1, a.atomic_number, a.symbol)

                    for mode_i in range(offset, offset + 5):
                        if mode_i >= max_modes:
                            break

                        mode = self.included_modes[mode_i]

                        r += '{:15.9f}'.format(
                            self.displacements[mode, index * 3 + c] * math.sqrt(quantities.AMUToElectronMass))

                    r += '\n'

        return r

    @staticmethod
    def vibrational_temperature(vibrational_energy, scale=1.0):
        """Convert a vibration energy (hartree) to a vibrational temperature (K)

        :param vibrational_energy: the vibrational energy (in Hartree)
        :type vibrational_energy: float
        :param scale: scale of the vibrational energyµ
        :type scale: float
        :rtype: float
        """
        return constants.h * vibrational_energy * scale * \
            MassWeightedHessian.VIB_ENERGY_CONVERSION / constants.Boltzmann

    @staticmethod
    def rotational_temperature(inertia_moment):
        """Convert an inertia moment (in AMU*angstrom**2) to a rotational temperature (K)

        :param inertia_moment: inertia moment
        :type inertia_moment: float
        :rtype: float
        """

        return constants.hbar ** 2 / (2 * inertia_moment * MassWeightedHessian.INERTIA_CONVERSION * constants.Boltzmann)

    def compute_partition_functions(self, symmetry_number, temperature=298.15, pressure=1.01325e5):
        """Compute the value of the different partition functions

        .. note::

            + Assume that the electronic partition is equal to the multiplicity (T<1000K?)
            + Assume that the rotational partition function is equal to the "high" temperature limit (T>100K?)
            + Assume that the first vibrational level is the zero energy (the "V=0" version of Gaussian)

        :param symmetry_number: symmetry number for the rotation
        :type symmetry_number: int
        :param temperature: temperature considered (in Kelvin)
        :type temperature: float
        :param pressure: pressure (in pascal)
        :type pressure: float
        :return: (omega_e, omega_t, omega_r, omega_v) (as a tuple)
        :rtype: tuple
        """

        omega_e = self.molecule.multiplicity

        inertia_moments, _ = numpy.linalg.eigh(self.molecule.moments_of_inertia())
        inertia_moments = sorted(list(inertia_moments))

        if self.linear:
            omega_r = (2 * self.INERTIA_CONVERSION * inertia_moments[2] * constants.Boltzmann * temperature) / \
                      (symmetry_number * constants.hbar ** 2)
        else:
            omega_r = math.sqrt(math.pi) / symmetry_number * \
                (2 * constants.Boltzmann * temperature / (constants.hbar ** 2)) ** (3 / 2) * \
                math.sqrt(self.INERTIA_CONVERSION ** 3 * qcip_math.prod(inertia_moments))

        omega_t = \
            (2 * math.pi * self.molecule.mass() * self.MASS_CONVERSION * constants.Boltzmann * temperature) \
            ** (3 / 2) * constants.R * \
            temperature / (pressure * constants.h ** 3 * constants.Avogadro)

        omega_v = 1.0

        for index, vib_energy in enumerate(self.frequencies):
            if index not in self.included_modes:
                continue
            omega_v *= 1 / \
                (1 - math.exp(-MassWeightedHessian.vibrational_temperature(vib_energy, self.scale) / temperature))

        return omega_e, omega_t, omega_r, omega_v

    def compute_zpva(self):
        """Compute the ZPVA (zero point vibrational averaging) energy from frequencies,
        but WITHOUT the electronic (SCF) energy

        :rtype: float
        """

        zpva = .0

        for index, vib_energy in enumerate(self.frequencies):
            if index not in self.included_modes:
                continue
            zpva += vib_energy * self.scale / 2

        return zpva

    def compute_internal_energy(self, temperature=298.15):
        """Compute the value of the different contribution (translation, rotation, vibration) to the internal energy,
        **IN HARTREE** !

        .. note::

            + No electronic contribution (since it is the energy of the system)
            + Assume that the rotational partition function is equal to the "high" temperature limit (T>100K?)
            + The vibrational contribution does not include ZPVA

        :param temperature: temperature considered (in Kelvin)
        :type temperature: float
        :return: (U_t, U_r, U_v) (as a tuple, in Hartree)
        :rtype: tuple
        """

        U_t = 3 / 2 * constants.R * temperature * self.ENERGY_IN_AU_CONVERSION

        if len(self.molecule) == 1:
            U_r = .0
        elif self.linear:
            U_r = constants.R * temperature * self.ENERGY_IN_AU_CONVERSION
        else:
            U_r = 3 / 2 * constants.R * temperature * self.ENERGY_IN_AU_CONVERSION

        U_v = .0

        for index, vib_energy in enumerate(self.frequencies):
            if index not in self.included_modes:
                continue
            theta_v = MassWeightedHessian.vibrational_temperature(vib_energy, scale=self.scale)
            U_v += theta_v / (math.exp(theta_v / temperature) - 1)

        U_v *= constants.R * self.ENERGY_IN_AU_CONVERSION

        return U_t, U_r, U_v

    def compute_enthalpy(self, temperature=298.15):
        """Compute the value of the different contribution (translation, rotation, vibration) to the enthalpy,
        **IN HARTREE** !

        .. note::

            Since :math:`H=U+RT`, the code just adds :math:`RT` to the translational contribution to internal energy,
            and returns the whole thing (so the same remarks applies, in particular: NO ZPVA)

        :param temperature: temperature considered (in Kelvin)
        :type temperature: float
        :return: (H_t, H_r, H_v) (as a tuple, in Hartree)
        :rtype: tuple
        """

        U_t, U_r, U_v = self.compute_internal_energy(temperature)
        return U_t + constants.R * temperature * self.ENERGY_IN_AU_CONVERSION, U_r, U_v

    def compute_entropy(self, symmetry_number, temperature=298.15, pressure=1.01325e5):
        """Compute the value of the different contribution (translation, rotation, vibration) to the entropy,
        **IN HARTREE / KELVIN** !

        .. note::

            + Assume that the rotational partition function is equal to the "high" temperature limit (T>100K?)
            + The vibrational contribution does include ZPVA (?)

        :param symmetry_number: symmetry number for the rotation
        :type symmetry_number: int
        :param temperature: temperature considered (in Kelvin)
        :type temperature: float
        :param pressure: pressure (in pascal)
        :type pressure: float
        :return: (S_t, S_r, S_v) (as a tuple, in Hartree / Kelvin)
        :rtype: tuple
        """

        omega_e, omega_t, omega_r, omega_v = self.compute_partition_functions(symmetry_number, temperature, pressure)

        S_t = constants.R * (math.log(omega_t) + 5 / 2) * self.ENERGY_IN_AU_CONVERSION

        if len(self.molecule) != 1:
            S_r = constants.R * (math.log(omega_r) + (1 if self.linear else 3 / 2)) * self.ENERGY_IN_AU_CONVERSION
        else:
            S_r = .0

        S_v = .0

        for index, vib_energy in enumerate(self.frequencies):
            if index not in self.included_modes:
                continue
            theta_over_T = MassWeightedHessian.vibrational_temperature(vib_energy, scale=self.scale) / temperature
            S_v += theta_over_T / (math.exp(theta_over_T) - 1) - math.log(1 - math.exp(-theta_over_T))

        S_v *= constants.R * self.ENERGY_IN_AU_CONVERSION

        return S_t, S_r, S_v

    def compute_gibbs_free_energy(self, symmetry_number, temperature=298.15, pressure=1.01325e5):
        """Compute the value of the different contribution (translation, rotation, vibration) to gibbs free energy,
        **IN HARTREE** !

        .. note::

            Since :math:`G=H-TS`, both enthalpy and entropy are computed, and that is it.
            No ZPVA, no electronic energy (as usual)

        :param symmetry_number: symmetry number for the rotation
        :type symmetry_number: int
        :param temperature: temperature considered (in Kelvin)
        :type temperature: float
        :param pressure: pressure (in pascal)
        :type pressure: float
        :return: (G_t, G_r, G_v) (as a tuple, in Hartree)
        :rtype: tuple
        """

        H_t, H_r, H_v = self.compute_enthalpy(temperature)
        S_t, S_r, S_v = self.compute_entropy(symmetry_number, temperature, pressure)

        return H_t - temperature * S_t, H_r - temperature * S_r, H_v - temperature * S_v
