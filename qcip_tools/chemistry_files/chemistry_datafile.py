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
                raise BadChemistryDataFile('type is incorrect')

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
            self.derivatives = ChemistryDataFile.read_derivatives_from_group(derivatives_group, self.spacial_dof)

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

            # derivatives
            derivatives_group = fx.create_group('derivatives')
            ChemistryDataFile.write_derivatives_in_group(derivatives_group, self.derivatives, self.spacial_dof)

    @staticmethod
    def write_derivatives_in_group(group, derivatives_, spacial_dof):
        """Write all the derivatives in a given group

        :param group: the group
        :type group: h5py.Group
        :param derivatives_: the derivatives to write
        :type derivatives_: dict
        :param spacial_dof: the DOF
        :type spacial_dof: int
        """
        derivatives_available = []

        for key in derivatives_:
            d = derivatives.Derivative(key, spacial_dof=spacial_dof)

            dataset_name = key
            if key == '':
                dataset_name = 'energy'

            derivatives_available.append(dataset_name)
            ChemistryDataFile.write_derivative_in_dataset(group, dataset_name, d, derivatives_[key])

        group.attrs['derivatives_available'] = ','.join(derivatives_available)

    @staticmethod
    def write_derivative_in_dataset(group, dataset_name, derivative, value):
        """Write a given derivative in a dataset (sadly, a dataset cannot be simply resized when created, so the
        dataset name is given)

        :param group: the group
        :type group: h5py.Group
        :param dataset_name: dataset name
        :type dataset_name: str
        :param derivative: the derivative
        :type derivative: derivatives.Derivative
        :param value: tensor
        :type: numpy.ndarray
        """

        shape = list(derivative.shape())

        if derivatives.is_electrical(derivative):
            num_freqs = len(value)
            freqs = ''
            first_freq = True
            shape.insert(0, num_freqs)
            super_array = numpy.zeros(shape)
            f = 0
            for freq in value:
                freqs += '{}{}:{}'.format(
                    ',' if not first_freq else '', 'f' if type(freq) is float else 's', freq)
                first_freq = False
                super_array[f] = value[freq]
                f += 1

            dset = group.create_dataset(dataset_name, shape, dtype='float64', data=super_array)
            dset.attrs['frequencies'] = freqs
        else:
            group.create_dataset(dataset_name, shape, dtype='float64', data=value)

    @staticmethod
    def read_derivatives_from_group(group, spacial_dof):
        """Read all the derivatives from a given group

        :param group: the group
        :type group: h5py.Group
        :param spacial_dof: the DOF
        :type spacial_dof: int
        :rtype: dict
        """
        if 'derivatives_available' not in group.attrs:
            raise BadChemistryDataFile('no list of derivatives available ?!?')

        derivatives_available = group.attrs['derivatives_available'].split(',')
        derivatives_ = {}

        for derivative in derivatives_available:
            dataset_name = derivative

            if derivative == 'energy':
                derivative = ''

            try:
                d = derivatives.Derivative(derivative, spacial_dof=spacial_dof)
            except derivatives.RepresentationError as e:
                raise BadChemistryDataFile(e)

            if dataset_name not in group:
                raise BadChemistryDataFile('could not find {} in {}'.format(dataset_name, group))

            derivatives_[derivative] = ChemistryDataFile.read_derivative_from_dataset(group[dataset_name], d)

        return derivatives_

    @staticmethod
    def read_derivative_from_dataset(dataset, derivative):
        """Read a derivative from a dataset

        :param dataset: the dataset
        :type dataset: h5py.Dataset
        :param derivative: the derivative
        :type derivative: derivatives.Derivative
        :rtype: numpy.ndarray
        """

        mat = dataset[:]

        if derivatives.is_electrical(derivative):
            if 'frequencies' not in dataset.attrs:
                raise BadChemistryDataFile('no frequencies for {}'.format(str(dataset)))
            frequencies = dataset.attrs['frequencies'].split(',')
            shape = derivative.shape()
            shape.insert(0, len(frequencies))
            if mat.shape != tuple(shape):
                raise BadChemistryDataFile('matrix shape ({}) is not right ({})'.format(mat.shape, shape))

            freqd_derivative = {}

            for index, frequency in enumerate(frequencies):
                if frequency[0] == 's':
                    f = frequency[2:]
                elif frequency[0] == 'f':
                    f = float(frequency[2:])
                else:
                    raise BadChemistryDataFile('frequency {} is not formated correctly for {}'.format(
                        frequency, derivative))

                freqd_derivative[f] = mat[index][:]

            return freqd_derivative
        else:
            if mat.shape != tuple(derivative.shape()):
                raise BadChemistryDataFile('matrix shape ({}) is not right ({})'.format(mat.shape, derivative.shape()))
            return mat
