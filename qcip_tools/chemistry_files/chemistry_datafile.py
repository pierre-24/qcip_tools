import numpy
import h5py

from qcip_tools import molecule as qcip_molecule, derivatives, atom as qcip_atom
from qcip_tools.chemistry_files import ChemistryFile, WithOutputMixin, WithMoleculeMixin, WithIdentificationMixin


class BadChemistryDataFile(Exception):
    pass

TYPE_UNICODE = h5py.special_dtype(vlen=bytes)


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

    Internally, use a ``h5py.File``, with:

    + ``/title`` (str) to store the title
    + ``/molecule/charge_and_multiplicity`` (2 int),
      ``/molecule/atoms`` (5xn array of floats: Z number, mass, position),
    + ``/derivatives`` with a ``derivatives_available`` attribute (str, comma separated list), and
      ``/derivatives/XXXX`` (array of floats, one for each derivatives),
      with a ``frequencies`` attribute to give the frequencies if required.

    """

    #: The identifier
    file_type = 'QCIP_CDF'
    version = 1

    def __init__(self):
        self.molecule = qcip_molecule.Molecule()
        self.spacial_dof = 0
        self.trans_plus_rot_dof = 0
        self.title = ''
        self.derivatives = {}

    @classmethod
    def possible_file_extensions(cls):
        return ['h5', 'hdf5']

    @classmethod
    def attempt_identification(cls, f):
        """it is an h5 file, and some specific datasets
        """

        must_be = ['title', 'version', '/molecule/charge_and_multiplicity', '/molecule/atoms', '/derivatives']

        try:
            with h5py.File(f.name, 'r') as fx:
                for m in must_be:
                    if m not in fx:
                        return False
                if fx['version'][0] != 1:
                    return False

                if 'type' not in fx['version'].attrs or fx['version'].attrs['type'] != ChemistryDataFile.file_type:
                    return False
        except:
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

        self.from_read = True
        must_be = ['title', 'version', '/molecule/charge_and_multiplicity', '/molecule/atoms', '/derivatives']

        with h5py.File(f.name, 'r') as fx:
            for m in must_be:
                if m not in fx:
                    raise BadChemistryDataFile('missing {}'.format(m))

            if fx['version'][0] != 1:
                raise BadChemistryDataFile('version > 1 (={})'.format(fx['version'][0]))

            if 'type' not in fx['version'].attrs or fx['version'].attrs['type'] != self.file_type:
                raise BadChemistryDataFile('type is incorect')

            atoms_shape = fx['/molecule/atoms'].shape
            if len(atoms_shape) != 2 or atoms_shape[1] != 5:
                raise BadChemistryDataFile('shape of /molecule/atoms is not correct ({})'.format(atoms_shape))

            self.title = fx['title'][0].decode()

            # molecule:
            self.molecule = qcip_molecule.Molecule()
            for i in range(atoms_shape[0]):
                info = fx['/molecule/atoms'][i]
                a = qcip_atom.Atom(atomic_number=int(info[0]))
                a.mass = info[1]
                a.position = info[2:]
                self.molecule.insert(a)

            self.molecule.charge, self.molecule.multiplicity = fx['/molecule/charge_and_multiplicity']
            self.spacial_dof = 3 * len(self.molecule)
            self.trans_plus_rot_dof = 5 if self.molecule.linear() else 6

            # derivatives
            derivatives_group = fx['/derivatives']

            if 'derivatives_available' not in derivatives_group.attrs:
                raise BadChemistryDataFile('no list of derivatives available ?!?')

            derivatives_available = derivatives_group.attrs['derivatives_available'].split(',')
            for derivative in derivatives_available:
                path = '/derivatives/{}'.format(derivative)
                try:
                    d = derivatives.Derivative(derivative, spacial_dof=self.spacial_dof)
                except derivatives.RepresentationError as e:
                    raise BadChemistryDataFile(e)

                if path not in fx:
                    raise BadChemistryDataFile('could not find {}'.format(path))

                mat = fx[path]
                self.derivatives[derivative] = {}

                if derivatives.is_electrical(d):
                    if 'frequencies' not in fx[path].attrs:
                        raise BadChemistryDataFile('no frequencies for {}'.format(path))
                    frequencies = fx[path].attrs['frequencies'].split(',')
                    shape = d.shape()
                    shape.insert(0, len(frequencies))
                    if mat.shape != tuple(shape):
                        raise BadChemistryDataFile('matrix shape ({}) is not right ({})'.format(mat.shape, shape))
                    for index, frequency in enumerate(frequencies):
                        f = None
                        if frequency[0] == 's':
                            f = frequency[2:]
                        elif frequency[0] == 'f':
                            f = float(frequency[2:])
                        else:
                            raise BadChemistryDataFile('frequency {} is not formated correctly for {}'.format(
                                frequency, derivative))

                        self.derivatives[derivative][f] = mat[index][:]
                else:
                    if mat.shape != tuple(d.shape()):
                        raise BadChemistryDataFile('matrix shape ({}) is not right ({})'.format(mat.shape, d.shape()))
                    self.derivatives[derivative] = mat[:]

    def to_string(self):
        raise NotImplementedError('to_string()')

    def write(self, f):

        with h5py.File(f.name, 'w') as fx:
            dset = fx.create_dataset('version', (1,), dtype='i', data=self.version)
            dset.attrs['type'] = self.file_type

            fx.create_dataset('title', (1,), maxshape=(None,), dtype=TYPE_UNICODE, data=self.title.encode())

            fx.create_dataset(
                '/molecule/charge_and_multiplicity',
                (2,),
                dtype='i',
                data=[self.molecule.charge, self.molecule.multiplicity])

            atoms_data_shape = (len(self.molecule), 5)
            atoms_data = numpy.zeros(atoms_data_shape)

            for i, a in enumerate(self.molecule):
                atoms_data[i][0] = a.atomic_number
                atoms_data[i][1] = a.mass
                atoms_data[i][2:] = a.position

            fx.create_dataset('/molecule/atoms', atoms_data_shape, data=atoms_data)

            derivatives_group = fx.create_group('derivatives')
            derivatives_available = ''
            first_derivative = True

            for key in self.derivatives:
                d = derivatives.Derivative(key, spacial_dof=self.spacial_dof)

                derivatives_available += '{}{}'.format(',' if not first_derivative else '', key)
                first_derivative = False

                shape = list(d.shape())
                if derivatives.is_electrical(key):
                    num_freqs = len(self.derivatives[key])
                    freqs = ''
                    first_freq = True
                    shape.insert(0, num_freqs)
                    super_array = numpy.zeros(shape)
                    f = 0
                    for freq in self.derivatives[key]:
                        freqs += '{}{}:{}'.format(
                            ',' if not first_freq else '', 'f' if type(freq) is float else 's', freq)
                        first_freq = False
                        super_array[f] = self.derivatives[key][freq]
                        f += 1

                    dset = derivatives_group.create_dataset(key, shape, data=super_array)
                    dset.attrs['frequencies'] = freqs
                else:
                    derivatives_group.create_dataset(key, shape, data=self.derivatives[key])

            derivatives_group.attrs['derivatives_available'] = derivatives_available
