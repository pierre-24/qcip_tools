import struct
import numpy

from qcip_tools import molecule, atom, datafile
from qcip_tools.chemistry_files import ChemistryFile, WithOutputMixin, WithMoleculeMixin, WithIdentificationMixin


class BadChemistryDataFile(Exception):
    pass


class ChemistryDataFile(ChemistryFile, WithOutputMixin, WithMoleculeMixin, WithIdentificationMixin):
    """A file to save quantum chemistry data: geometry and derivatives.

    .. container:: class-members

        + ``self.molecule``: the molecule (``qcip_tools.molecule.Molecule``)
        + ``self.spacial_dof``: spacial degree of freedom = 3N (``int``)
        + ``self.trans_plus_rot_dof``: either 5 or 6 (``int``)
        + ``self.title``: the title (``str``)
        + ``self.derivatives``: derivatives of the energy (``dict``): the key is a representation of the derivative,
          the value is either a ``Tensor``, or a dict in case of electrical derivatives.
          This new dict is subdivided into the different frequencies.

          .. code-block:: text

            + "F"
                + "static" : ElectricDipole
            + "FF":
                + "static" : PolarizabilityTensor
            + "FD":
                + "0.042" : PolarizabilityTensor
                + ...
            + "FDD":
                + "static" : FirstHyperpolarizabilityTensor
                + ...
            + "G": BaseGeometricalDerivativeTensor
            + "GG": BaseGeometricalDerivativeTensor
            + ...

    Internally, use a ``BinaryDataFile`` (with a modified magic number, ``0xCEB1DAF1``), with:

    + ``title`` (str) to store the title
    + ``molecule_charge_and_multiplicity`` (2 int),
      ``molecule_atoms`` (5xn array of floats: Z number, mass, position),
    + ``derivatives_available`` (str, comma separated list), and
      ``d:XXXXX`` (array of floats, one for each derivatives,
      may be followed by the frequency if electrical derivative) to store derivatives.

    """

    #: The identifier
    file_type = 'QCIP_CDF'
    magic_number = 0xCEB1DAF1
    requires_binary_mode = True
    version = 1

    def __init__(self):
        self.molecule = molecule.Molecule()
        self.spacial_dof = 0
        self.trans_plus_rot_dof = 0
        self.title = ''
        self.derivatives = {}

    @classmethod
    def possible_file_extensions(cls):
        return ['chdf']

    @classmethod
    def attempt_identification(cls, f):
        """use binary datafile
        """

        # check magic number
        magic = struct.unpack('<I', f.read(4))[0]

        if magic != ChemistryDataFile.magic_number:
            return False

        f.seek(0)
        fb = datafile.BinaryDataFile()
        fb.magic_number = ChemistryDataFile.magic_number

        try:
            fb.read(f)
        except:
            return False

        must_be = ['version', 'title', 'molecule_charge_and_multiplicity', 'molecule_atoms', 'derivatives_available']
        for a in must_be:
            if a not in fb:
                return False

        return True

    @classmethod
    def from_molecule(cls, molecule, title='', *args, **kwargs):
        obj = super().from_molecule(molecule, *args, **kwargs)
        obj.title = title

        return obj

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        if 'b' not in f.mode:
            raise IOError('file {} must be open in binary mode!'.format(f.name))

        fx = datafile.BinaryDataFile()
        fx.magic_number = self.magic_number
        fx.read(f)
        self.from_read = True

        must_be = ['title', 'molecule_charge_and_multiplicity', 'molecule_atoms', 'derivatives_available', 'version']
        for a in must_be:
            if a not in fx:
                raise BadChemistryDataFile('{} not in chemistry data file', a)

        self.title = fx['title']

        # get molecule
        self.molecule = molecule.Molecule()
        for i in range(fx['molecule_atoms'].shape[0]):
            infos = fx['molecule_atoms'][i]
            a = atom.Atom(atomic_number=int(infos[0]), position=infos[2:])
            a.mass = infos[1]
            self.molecule.insert(a)

        self.molecule.charge, self.molecule.multiplicity = fx['molecule_charge_and_multiplicity']
        self.spacial_dof = 3 * len(self.molecule)
        self.trans_plus_rot_dof = 6 if not self.molecule.linear() else 5

        # get derivatives
        derivatives_available = fx['derivatives_available'].split(',')
        for d in derivatives_available:
            key = 'd:' + d
            if key not in fx:
                raise BadChemistryDataFile('{} is not available, but is present in the list! ({})'.format(
                    key, fx['derivatives_available']))

            self.derivatives[d] = fx[key]

    def to_string(self):
        raise NotImplementedError('to_string for Binary Data!')

    def write(self, f):
        if 'b' not in f.mode:
            raise IOError('file {} must be open in binary mode!'.format(f.name))

        positions_list = numpy.zeros((len(self.molecule), 5))

        for i, a in enumerate(self.molecule):
            positions_list[i][0] = a.atomic_number
            positions_list[i][1] = a.mass
            positions_list[i][2:] = a.position

        fx = datafile.BinaryDataFile()
        fx.magic_number = self.magic_number
        fx.set('version', 'I', [self.version])
        fx.set('title', 'S', self.title)
        fx.set('molecule_charge_and_multiplicity', 'I', [self.molecule.charge, self.molecule.multiplicity])
        fx.set('molecule_atoms', 'A', positions_list)

        derivatives_available = ','.join(self.derivatives.keys())
        fx.set('derivatives_available', 'S', derivatives_available)
        for key, derivative in self.derivatives.items():
            fx.set('d:' + key, 'A', derivative)

        fx.write(f)
