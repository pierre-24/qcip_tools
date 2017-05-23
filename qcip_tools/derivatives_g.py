import numpy
import math

from qcip_tools import derivatives, quantities

#: Correspondence between a name and a representation
REPRESENTATIONS = {
    'gradient': 'G',
    'hessian': 'GG',
    'cubic_FF': 'GGG',
}

DERIVATIVES = list(REPRESENTATIONS)

#: List of all simplified names
SIMPLIFIED_NAMES = {
    'energy': 'E',
    'gradient': 'Fi',
    'hessian': 'Hij',
    'cubic_FF': 'Fijk',
}


class BaseGeometricalDerivativeTensor(derivatives.Tensor):
    """Base for all derivatives of the energy with respect to geometrical parameters
    """

    def __init__(self, spacial_dof, trans_plus_rot=0, representation='', components=None):

        if 'F' in representation or 'D' in representation:
            raise ValueError(representation)

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

    def project_over_normal_modes(self, displacements):
        """Project over normal modes.

        .. note::

            In the absence of a general formula, this is only valid for gradient, hessian and cubic force field
            projections.

        :param displacements: carthesian displacements (eigenvectors of mass-weighted hessian)
        :type displacements: numpy.ndarray
        :rtype: BaseGeometricalDerivativeTensor
        """
        if 'N' is self.representation:
            raise Exception('already projected ?!?')

        if displacements.shape != (self.spacial_dof, self.spacial_dof):
            raise ValueError(displacements)

        if self.representation == 'G':
            projected_tensor = numpy.dot(displacements, self.components)
        elif self.representation == 'GG':
            projected_tensor = numpy.dot(numpy.dot(displacements, self.components), displacements.transpose())
        elif self.representation == 'GGG':
            projected_tensor = numpy.dot(
                displacements, numpy.dot(numpy.dot(displacements, self.components), displacements.transpose()))
        else:
            raise Exception('no projection defined for {}'.format(self.representation.representation()))

        return BaseGeometricalDerivativeTensor(
            self.spacial_dof,
            self.trans_plus_rot,
            self.representation.representation().replace('G', 'N'),
            components=projected_tensor)


class MassWeightedHessian:
    """To be exact, this is not a derivatives of the energy. So it is not a child of ``derivatives.Tensor`` !

    :type molecule: qcip_tools.molecule.Molecule
    :type carthesian_hessian: numpy.array
    """

    def __init__(self, molecule, carthesian_hessian=None):

        self.molecule = molecule
        self.linear = self.molecule.linear()

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

        if carthesian_hessian is not None:
            self.from_cartesian_hessian(carthesian_hessian)

    def from_cartesian_hessian(self, cartesian_hessian):
        """From a given carthesian hessian, create mass-weigted ones, and extracts frequencies and displacements.

        :param cartesian_hessian: the carthesian hessian
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

        HartreeToWavenumber = quantities.convert(quantities.ureg.hartree, quantities.ureg.wavenumber)

        for offset in range(0, self.vibrational_dof, 5):

            line_mode = ' ' * 10
            line_freqs = 'Frequencies   '
            line_rmass = 'Reduced-masses'

            for mode in range(self.trans_plus_rot_dof + offset, self.trans_plus_rot_dof + offset + 5):
                if mode >= self.dof:
                    break
                line_mode += '{:15}'.format(mode + 1 - self.trans_plus_rot_dof)
                line_freqs += '{: 15.4f}'.format(self.frequencies[mode] * HartreeToWavenumber)
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

                    for mode in range(self.trans_plus_rot_dof + offset, self.trans_plus_rot_dof + offset + 5):
                        if mode >= self.dof:
                            break

                        r += '{:15.9f}'.format(
                            self.displacements[mode, index * 3 + c] * math.sqrt(quantities.AMUToElectronMass))

                    r += '\n'

        return r
