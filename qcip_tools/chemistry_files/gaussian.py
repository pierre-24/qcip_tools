import os
import math
import re
import numpy
import copy

from qcip_tools import molecule, atom, quantities, basis_set, derivatives, derivatives_e, derivatives_g, \
    derivatives_exci
from qcip_tools.chemistry_files import ChemistryFile, WithOutputMixin, WithMoleculeMixin, ChemistryLogFile, \
    FormatError, WithIdentificationMixin, PropertyNotPresent


class InputFormatError(FormatError):
    pass


class Input(ChemistryFile, WithOutputMixin, WithMoleculeMixin, WithIdentificationMixin):
    """Gaussian input file.

    .. container:: class-members

        + ``self.molecule``: the molecule (``qcip_tools.molecule.Molecule``)
        + ``self.options``: the different options (starting with ``%`` in the file, ``dict`` of ``str``
          [if ``%xxx=yyy``, key is ``xxx`` and value is ``yyy``])
        + ``self.input_card``: the input card (starting with ``#`` in the file, ``list`` of ``str``
          [the different lines])
        + ``self.title``: the title (``str``)
        + ``self.other_blocks``: the additionnal blocks after the molecule
          (``list`` of ``list`` of ``str`` [the different lines])
    """

    #: The identifier
    file_type = 'GAUSSIAN_INP'

    def __init__(self):

        self.molecule = molecule.Molecule()

        self.other_blocks = []
        self.title = ''
        self.options = {}
        self.input_card = []

    @classmethod
    def possible_file_extensions(cls):
        return ['com', 'inp', 'gau']

    @classmethod
    def attempt_identification(cls, f):
        """A gaussian input should contains '%' or '#' as the first character of the line of the "first block"
        """

        found_some_lines = 0
        found_input_line = False

        for line in f.readlines():
            if line.strip() == '':
                break
            if line[0] not in ['%', '#']:
                if not found_some_lines:
                    return False
            else:
                if line[0] == '#':
                    found_input_line = True
                found_some_lines += 1

        if not found_some_lines:
            return False

        if not found_input_line:
            return False

        # TODO: if ``found_some_lines`` is equal to 1, maybe an extra test ? (charge and multiplicity, for example)
        return True

    @classmethod
    def from_molecule(cls, molecule, title='', input_card='', options=None, other_blocks=None, *args, **kwargs):
        """Create a file from molecule

        :param molecule: the molecule
        :type molecule: qcip_tools.molecule.Molecule
        :param title: title of the run
        :type title: str
        :param input_card: the input card (starting with ``#P``)
        :type input_card: str|list
        :param options: calculation options
        :type options: dict
        :param other_blocks: the other input blocks
        :type other_blocks: list
        :rtype: qcip_tools.chemistry_files.gaussian.Input
        """

        obj = super().from_molecule(molecule, *args, **kwargs)
        obj.title = title

        if type(input_card) is str:
            obj.input_card = [input_card]
        elif type(input_card) is list:
            obj.input_card = input_card
        else:
            raise TypeError(input_card)

        obj.options = {} if options is None else options
        obj.other_blocks = [] if other_blocks is None else other_blocks
        return obj

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        self.from_read = True
        lines = f.readlines()

        raw_blocks = []

        temp = []
        for line in lines:
            if len(temp) > 0 and line.strip() == '':
                raw_blocks.append(temp)
                temp = []
            else:
                temp.append(line.strip())

        # add last
        if len(temp) > 0:
            raw_blocks.append(temp)
        if len(raw_blocks) < 3:
            raise InputFormatError('The file contains less than 3 raw_blocks')

        # try to understand first block
        for line in raw_blocks[0]:
            if line[0] == '%':
                opt = line[1:].split('=')
                self.options[opt[0].lower()] = opt[1]
            else:
                self.input_card.append(line)

        # read title, charge and multiplicity:
        self.title = ' '.join(raw_blocks[1])
        charge_and_multiplicity = [int(i) for i in raw_blocks[2][0].strip().split()]
        self.molecule.charge = charge_and_multiplicity[0]
        self.molecule.multiplicity = charge_and_multiplicity[1]

        # read coordinates:
        coordinates = raw_blocks[2][1:]
        for atm in coordinates:
            s = atm.split()
            if len(s) != 4:
                raise InputFormatError('Wrong atom definition "' + atm + '"')
            else:
                coords = [float(i) for i in s[1:]]
                atm_ = s[0]
                if atm_.isnumeric():
                    atm_ = atom.AtomicNumberToSymbol[int(atm_)]
                self.molecule.insert(atom.Atom(symbol=atm_, position=coords))

        # then store "other blocks"
        if len(raw_blocks) > 3:
            self.other_blocks = raw_blocks[3:]

    def to_string(self, chk='', original_chk=False):
        """

        :param chk: new chk name if needed
        :type chk: str
        :param original_chk: use original chk name
        :type original_chk: bool
        :rtype: str
        """

        rstr = ''

        # first, options:
        for key in self.options:
            if not original_chk and key == 'chk':
                rstr += '%{}={}\n'.format(key, chk)
            else:
                rstr += '%{}={}\n'.format(key, self.options[key])

        # then, input card and title
        rstr += '\n'.join(self.input_card) + '\n\n'
        rstr += self.title + '\n\n'

        # write charge and multiplicity
        rstr += str(self.molecule.charge) + ' ' + str(self.molecule.multiplicity) + '\n'

        # write coordinates
        rstr += self.molecule.output_atoms()
        rstr += '\n'

        # write other blocks :
        for i, block in enumerate(self.other_blocks):
            rstr += '\n'.join(block)
            rstr += '\n\n'
        return rstr

    def write(self, f, chk=None, original_chk=False):
        """

        :param f: File
        :type f: file
        """

        if chk is None:
            chk = os.path.basename('.'.join(f.name.split('.')[:-1]))

        f.write(self.to_string(chk=chk, original_chk=original_chk))


FCHK_AUTHORIZED_TYPES = ['I', 'R', 'C']
REGEX_FLOAT = re.compile(r'(?P<base>[0-9](\.[0-9]*)?)\-(?P<exp>[0-9]*)')


def transform_string_from_fchk(data, data_type):
    if data_type == 'I':
        return int(data)
    elif data_type == 'R':
        try:
            return float(data)
        except ValueError:
            if REGEX_FLOAT.match(data):
                return 0.0
    elif data_type == 'C':
        return data
    else:
        raise Exception('Unknown data type {}'.format(data_type))


class FCHKFormatError(FormatError):
    pass


class FCHKChunkInformation:
    """Where to find a given keyword in the FCHK file

    :param keyword: keyword related to this piece of information
    :type keyword: str
    :param data_type: type of this piece of information
    :type data_type: str
    :param data_length: length of this data
    :type data_length: int
    :param line_start: start in the source file
    :type line_start: int|object
    :param line_end: end in the source file
    :type line_end: int|object
    :param is_matrix: is it a matrix ?
    :type is_matrix: bool
    """

    def __init__(self, keyword, data_type, data_length, line_start, line_end, is_matrix=True):
        self.keyword = keyword
        self.data_type = data_type
        self.data_length = data_length
        self.line_start = line_start
        self.line_end = line_end
        self.is_matrix = is_matrix


class FCHK(ChemistryFile, WithMoleculeMixin, WithIdentificationMixin):
    """A FCHK file. Based on the same principle as DataFile (split into chunks, interpret and store after).

    Most of the information are found in the "FChk file" part of the
    `"Interfacing" documentation of Gaussian <http://gaussian.com/interfacing/>`_.

    .. container:: class-members

        + ``self.molecule``: the molecule (``qcip_tools.molecule.Molecule``)
        + ``self.title``: the title (``str``)
        + ``self.calculation_type``: the type of calculation (``str``)
        + ``self.calculation_method``: the method used (``str``)
        + ``self.basis_set``: the basis set (``str``)
        + ``self.chunks_information``: information about every chunk (``dict`` of ``FCHKChunkInformation``)
        + ``self.chunks_parsed``: parsed chunks are stored here (``dict`` of ``int|float|str|numpy.ndarray``)
        + ``self.lines``: the different lines of the file, used for interpretation (``list`` of ``str``)

    """

    #: The identifier
    file_type = 'GAUSSIAN_FCHK'

    def __init__(self):
        self.molecule = molecule.Molecule()

        self.title = ''
        self.calculation_type = ''
        self.calculation_method = ''
        self.basis_set = ''

        self.chunks_information = {}
        self.chunks_parsed = {}
        self.lines = []
        self.contains_ECP = False
        self.frozen_electrons = 0

    @classmethod
    def possible_file_extensions(cls):
        return ['fchk']

    @classmethod
    def attempt_identification(cls, f):
        """A gaussian fchk starts with two two line of info, then "Number of atoms"
        """

        i = 0

        for line in f.readlines():
            if i == 2:
                if 'Number of atoms' not in line:
                    return False
            if i == 3:
                if 'Info1-9' in line:
                    return True
                return False
            i += 1

        return False

    @classmethod
    def from_molecule(cls, molecule, *args, **kwargs):
        raise NotImplementedError('from_molecule')

    def get(self, key):
        """Get the content from a given keyword

        :param key: the keyword
        :type key: str
        :rtype: list|int|float|str
        """
        content = []

        if key not in self.chunks_information:
            raise KeyError(key)

        if key not in self.chunks_parsed:
            chunk = self.chunks_information[key]

            if chunk.is_matrix:
                if chunk.data_type in ['R', 'I']:
                    matrix = self.lines[chunk.line_start:chunk.line_end]
                    matrix = (''.join(matrix)).split()
                    if len(matrix) == chunk.data_length:
                        for d in matrix:
                            content.append(transform_string_from_fchk(d, chunk.data_type))
                    else:
                        raise FCHKFormatError('Size is not the same as defined for {}'.format(key))
                else:
                    content = ''.join(a.strip() for a in self.lines[chunk.line_start:chunk.line_end])
            else:
                content = self.lines[chunk.line_start][49:].strip()
                content = transform_string_from_fchk(content, chunk.data_type)

            self.chunks_parsed[key] = content
        else:
            content = self.chunks_parsed[key]

        return content

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        self.from_read = True
        self.lines = f.readlines()

        # extract the first information :
        self.title = self.lines[0].strip()
        second_line = self.lines[1]
        self.calculation_type = second_line[0:10].strip()
        self.calculation_method = second_line[10:40].strip()
        self.basis_set = second_line[40:].strip()

        # now, crawl trough information
        current_line_index = 2
        while True:
            current_line = self.lines[current_line_index]
            keyword = current_line[:40].strip()  # name is contained in the 40 first characters
            size = 1
            data_type = current_line[43]

            if data_type not in FCHK_AUTHORIZED_TYPES:
                raise FCHKFormatError('Type of {} ({}) is unknown'.format(keyword, data_type))

            is_matrix = current_line[47] == 'N'

            if is_matrix:
                size = int(current_line[49:].strip())
                if data_type in ['I', 'R', 'C']:
                    values_per_line = 5
                    if data_type == 'I':
                        values_per_line = 6
                    num_lines = math.ceil(size / values_per_line)
                else:
                    raise FCHKFormatError('Unknown format {}!'.format(data_type))

                line_start, line_end = current_line_index + 1, current_line_index + 1 + num_lines
                current_line_index = line_end

            else:
                line_start = line_end = current_line_index
                current_line_index += 1

            self.chunks_information[keyword] = FCHKChunkInformation(
                keyword, data_type, size, line_start, line_end, is_matrix)

            if current_line_index >= len(self.lines):
                break

        # finally, read the molecule
        nuclear_charges = self.get('Atomic numbers')
        atoms_coordinates = self.get('Current cartesian coordinates')
        weights = self.get('Real atomic weights')

        for index, c in enumerate(nuclear_charges):
            a = atom.Atom(
                atomic_number=c,
                position=[a * quantities.AuToAngstrom for a in atoms_coordinates[index * 3:index * 3 + 3]])
            if weights[index] > .0:  # for a reason that I don't get, fchk files sometimes contains only zeroes
                a.mass = weights[index]
            self.molecule.insert(a)

        self.molecule.charge = self.get('Charge')
        self.molecule.multiplicity = self.get('Multiplicity')

        num_elec = self.get('Number of electrons')
        if 'ECP-RNFroz' in self.chunks_information:
            self.contains_ECP = True
            self.frozen_electrons = sum(self.get('ECP-RNFroz'))
            num_elec += self.frozen_electrons

        if num_elec != self.molecule.number_of_electrons():
            raise ValueError('Number of electron does not match: {} != {}'.format(
                self.get('Number of electrons'), self.molecule.number_of_electrons()))

    def __getitem__(self, item):
        """
        Implement access trough ``[]``
        `"""
        return self.get(item)

    def __contains__(self, item):
        """
        Implement test with ``in`` keyword
        """
        return item in self.chunks_information


@FCHK.define_property('computed_energies')
def gaussian__fchk__property__computed_energies(obj, *args, **kwargs):
    """Get the energies. Returns a dictionary of the energies at different level of approximation.

    :param obj: object
    :type obj: qcip_tools.chemistry_files.gaussian.FCHK
    :rtype: dict
    """

    possible_energies = [
        # either HF or DFT:
        ('SCF Energy', 'SCF/DFT'),
        # MPx:
        ('MP2 Energy', 'MP2'),
        ('MP3 Energy', 'MP3'),
        ('MP4D Energy', 'MP4D'),
        ('MP4DQ Energy', 'MP4DQ'),
        ('MP4SDQ Energy', 'MP4SDQ'),
        ('Cluster Energy', 'CCSD'),
        ('Cluster Energy with triples', 'CCSD(T)'),
        # CIx:
        ('CIS Energy', 'CIS'),
        ('CIS-MP2 Energy', 'CIS(D)'),
        # Total energy:
        ('Total Energy', 'total')
    ]

    found_energies = {}

    for (key, method) in possible_energies:
        if key in obj:
            found_energies[method] = obj.get(key)

    # try to determine wheter it is HF or DFT
    if 'SCF/DFT' not in found_energies:
        raise PropertyNotPresent('computed_energies')

    if len(found_energies) == 2:
        found_energies[obj.calculation_method[1:]] = found_energies['SCF/DFT']
    else:
        found_energies['HF'] = found_energies['SCF/DFT']

    return found_energies


def beta_SHG_from_fchk(compressed_tensor, offset=0):
    """SHG beta tensor in FCHK is a bit messed up

    :param offset: frequency offset
    :type offset: int
    :rtype: numpy.ndarray
    """
    tensor = numpy.zeros((3, 3, 3))
    for i in range(0, 3):
        for j in range(0, 3):
            for k in range(0, j + 1):
                n = j + k
                if j > 1:
                    n += 1
                tensor[i, j, k] = -1 * compressed_tensor[offset * 18 + i * 6 + n]
                if j != k:
                    tensor[i, k, j] = tensor[i, j, k]

    return tensor


def beta_EOP_from_fchk(compressed_tensor, offset=0):
    """EOP beta tensor in FCHK is a bit messed up

    :param offset: frequency offset
    :type offset: int
    :rtype: numpy.ndarray
    """

    to_consider = [
        (0, 0, 0), (0, 1, 0), (1, 1, 0), (2, 0, 0), (2, 1, 0), (2, 2, 0),
        (0, 0, 1), (0, 1, 1), (1, 1, 1), (2, 0, 1), (2, 1, 1), (2, 2, 1),
        (0, 0, 2), (0, 1, 2), (1, 1, 2), (2, 0, 2), (2, 1, 2), (2, 2, 2),
    ]

    tensor = numpy.zeros((3, 3, 3))

    for index, components in enumerate(to_consider):
        v = - compressed_tensor[offset * 18 + index]
        tensor[components] = v
        if components[0] != components[1]:
            tensor[components[1], components[0], components[2]] = v

    return tensor


def gamma_from_gaussian(compressed_tensor, offset=0, exchange_23=False):
    """gamma tensor in FCHK is (also) a bit messed up

    :param offset: frequency offset
    :type offset: int
    :param exchange_23: exchange the third and the fourth component if True (second and third otherwise)
    :type exchange_23: bool
    :rtype: numpy.ndarray
    """

    to_consider = [
        (0, 0, 0, 0), (0, 1, 0, 0), (0, 1, 1, 0), (0, 2, 0, 0), (0, 2, 1, 0), (0, 2, 2, 0),
        (1, 0, 0, 0), (1, 1, 0, 0), (1, 1, 1, 0), (1, 2, 0, 0), (1, 2, 1, 0), (1, 2, 2, 0),
        (2, 0, 0, 0), (2, 1, 0, 0), (2, 1, 1, 0), (2, 2, 0, 0), (2, 2, 1, 0), (2, 2, 2, 0),
        (0, 0, 0, 1), (0, 1, 0, 1), (0, 1, 1, 1), (0, 2, 0, 1), (0, 2, 1, 1), (0, 2, 2, 1),
        (1, 0, 0, 1), (1, 1, 0, 1), (1, 1, 1, 1), (1, 2, 0, 1), (1, 2, 1, 1), (1, 2, 2, 1),
        (2, 0, 0, 1), (2, 1, 0, 1), (2, 1, 1, 1), (2, 2, 0, 1), (2, 2, 1, 1), (2, 2, 2, 1),
        (0, 0, 0, 2), (0, 1, 0, 2), (0, 1, 1, 2), (0, 2, 0, 2), (0, 2, 1, 2), (0, 2, 2, 2),
        (1, 0, 0, 2), (1, 1, 0, 2), (1, 1, 1, 2), (1, 2, 0, 2), (1, 2, 1, 2), (1, 2, 2, 2),
        (2, 0, 0, 2), (2, 1, 0, 2), (2, 1, 1, 2), (2, 2, 0, 2), (2, 2, 1, 2), (2, 2, 2, 2),
    ]

    tensor = numpy.zeros((3, 3, 3, 3))

    for index, components in enumerate(to_consider):
        v = compressed_tensor[offset * 54 + index]
        tensor[components] = v
        if not exchange_23 and components[1] != components[2]:
            tensor[components[0], components[2], components[1], components[3]] = v
        if exchange_23 and components[2] != components[3]:
            tensor[components[0], components[1], components[3], components[2]] = v

    return tensor


@FCHK.define_property('electrical_derivatives')
def gaussian__fchk__property__electrical_derivatives(obj, *args, **kwargs):
    """Get electrical derivatives in FCHK. Returns a dictionary of dictionary:

    .. code-block:: text

        + "F"
            + static : ElectricDipole
        + "FF":
            + static : PolarizabilityTensor
        + "FD":
            + 0.042 : PolarizabilityTensor
            + ...
        + "FDD":
            + static : FirstHyperpolarizabilityTensor
            + ...
        + ...

    Note that the dipole moment is present in more or less every cases.

    Does not fetch (yet?) the second hyperpolarizability, stored as ``Derivative Beta(-w,w,0)``
    and ``Derivative Beta(w,w,-2w)`` (that will be fun for reorganisation).

    :param obj: object
    :type obj: qcip_tools.chemistry_files.gaussian.FCHK
    :rtype: dict
    """

    electrical_derivatives = {}

    if 'Dipole Moment' in obj:
        electrical_derivatives['F'] = {'static': derivatives_e.ElectricDipole(dipole=obj.get('Dipole Moment'))}

    # statics ones
    if 'Polarizability' in obj:
        compressed_tensor = obj.get('Polarizability')
        o = derivatives_e.PolarisabilityTensor(input_fields=(0,))
        to_consider = [(0, 0), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2)]

        for index, components in enumerate(to_consider):
            for uniques in o.representation.inverse_smart_iterator(components):
                o.components[uniques] = compressed_tensor[index]

        electrical_derivatives['FF'] = {'static': o}

    if 'HyperPolarizability' in obj:
        compressed_tensor = obj.get('HyperPolarizability')
        o = derivatives_e.FirstHyperpolarisabilityTensor(input_fields=(0, 0))
        to_consider = [
            (0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1), (2, 0, 0), (2, 0, 1), (2, 1, 1), (2, 2, 0), (2, 2, 1),
            (2, 2, 2)]

        for index, components in enumerate(to_consider):
            for uniques in o.representation.inverse_smart_iterator(components):
                o.components[uniques] = - compressed_tensor[index]

        electrical_derivatives['FFF'] = {'static': o}

    # dynamic ones
    frequencies = []
    if 'Frequencies for FD properties' in obj:
        frequencies = obj.get('Frequencies for FD properties')

    if frequencies and len(frequencies) > 1:
        if 'Alpha(-w,w)' in obj:
            compressed_tensor = obj.get('Alpha(-w,w)')
            num_of_frequencies = int(len(compressed_tensor) / 9)

            if num_of_frequencies != len(frequencies):
                raise Exception('not the right number of frequencies for alpha(-w;w)')

            data = {}

            for i in range(1, num_of_frequencies):
                o = derivatives_e.PolarisabilityTensor(frequency=frequencies[i], input_fields=(1,))
                o.components = numpy.array(compressed_tensor[9 * i:9 * (i + 1)]).reshape(o.representation.shape())
                data[frequencies[i]] = o

            electrical_derivatives['dD'] = data

        if 'Beta(w,w,-2w)' in obj:
            compressed_tensor = obj.get('Beta(w,w,-2w)')
            num_of_frequencies = int(len(compressed_tensor) / 18)

            if num_of_frequencies != len(frequencies):
                raise Exception('not the right number of frequencies for beta(-2w;w,w)')

            data = {}

            for i in range(1, num_of_frequencies):
                o = derivatives_e.FirstHyperpolarisabilityTensor(frequency=frequencies[i], input_fields=(1, 1))
                o.components = beta_SHG_from_fchk(compressed_tensor, i)
                data[frequencies[i]] = o

            electrical_derivatives['XDD'] = data

        if 'Beta(-w,w,0)' in obj:
            compressed_tensor = obj.get('Beta(-w,w,0)')
            num_of_frequencies = int(len(compressed_tensor) / 18)

            if num_of_frequencies != len(frequencies):
                raise Exception('not the right number of frequencies for beta(-w;w,0)')

            data = {}

            for i in range(1, num_of_frequencies):
                o = derivatives_e.FirstHyperpolarisabilityTensor(frequency=frequencies[i], input_fields=(0, 1))
                o.components = beta_EOP_from_fchk(compressed_tensor, i)
                data[frequencies[i]] = o

            electrical_derivatives['dDF'] = data

        if 'Derivative Beta(-w,w,0)' in obj:
            compressed_tensor = obj.get('Derivative Beta(-w,w,0)')
            num_of_frequencies = int(len(compressed_tensor) / 54)

            if num_of_frequencies != len(frequencies):
                raise Exception('not the right number of frequencies for gamma(-w;0,0,w)')

            # catch static:
            electrical_derivatives['FFFF'] = {'static': derivatives_e.SecondHyperpolarizabilityTensor()}
            electrical_derivatives['FFFF']['static'].components = gamma_from_gaussian(compressed_tensor, 0)

            # catch dynamic:
            data = {}

            for i in range(1, num_of_frequencies):
                o = derivatives_e.SecondHyperpolarizabilityTensor(frequency=frequencies[i], input_fields=(1, 0, 0))
                o.components = gamma_from_gaussian(compressed_tensor, i)
                data[frequencies[i]] = o

            electrical_derivatives['dFFD'] = data

        if 'Derivative Beta(w,w,-2w)' in obj:
            compressed_tensor = obj.get('Derivative Beta(w,w,-2w)')
            num_of_frequencies = int(len(compressed_tensor) / 54)

            if num_of_frequencies != len(frequencies):
                raise Exception('not the right number of frequencies for gamma(-2w;w,w,0)')

            data = {}

            for i in range(1, num_of_frequencies):
                o = derivatives_e.SecondHyperpolarizabilityTensor(frequency=frequencies[i], input_fields=(1, 1, 0))
                o.components = gamma_from_gaussian(compressed_tensor, i)
                data[frequencies[i]] = o

            electrical_derivatives['XDDF'] = data

    if not electrical_derivatives:
        raise PropertyNotPresent('qs:electrical_derivatives')

    return electrical_derivatives


def make_cartesian_hessian_from_fchk(half_cartesian_hessian, spacial_dof):
    """Half the hessian is stored in fchk, get full hessian back

    :param half_cartesian_hessian: the half of the hessian, as retrieved from FCHK
    :type half_cartesian_hessian: list
    :param spacial_dof: spacial degree of freedom (DOF)
    :type spacial_dof: int
    :return:
    """

    if spacial_dof % 3 != 0:
        raise ValueError('DOF should be a multiple of 3 ?!?')

    full_cartesian_hessian = numpy.zeros((spacial_dof, spacial_dof))

    for i in range(spacial_dof):
        for j in range(0, i + 1):
            pos = int((i * (i + 1)) / 2 + j)
            full_cartesian_hessian[i, j] = half_cartesian_hessian[pos]
            if i != j:
                full_cartesian_hessian[j, i] = full_cartesian_hessian[i, j]

    return full_cartesian_hessian


@FCHK.define_property('geometrical_derivatives')
def gaussian__fchk__property__geometrical_derivatives(obj, *args, **kwargs):
    """Get geometrical derivatives in FCHK. Returns a dictionary of dictionary:

    .. code-block:: text

        + "G": BaseGeometricalDerivativeTensor
        + "GG": BaseGeometricalDerivativeTensor

    Does not fetch (yet?) the cubic force constants.

    :param obj: object
    :type obj: qcip_tools.chemistry_files.gaussian.FCHK
    :rtype: dict
    """

    geometrical_derivatives = {}
    spacial_dof = 3 * len(obj.molecule)
    trans_plus_rot_dof = 5 if obj.molecule.linear() else 6

    if obj.calculation_type.lower() not in ['freq', 'force']:
        raise PropertyNotPresent('geometrical_derivatives')

    if 'Cartesian Gradient' in obj:
        gradient = obj.get('Cartesian Gradient')

        if len(gradient) != 3 * len(obj.molecule):
            raise Exception('Wrong number of DOF in gradient')

        geometrical_derivatives['G'] = derivatives_g.BaseGeometricalDerivativeTensor(
            spacial_dof=spacial_dof, trans_plus_rot=trans_plus_rot_dof, representation='G', components=gradient)

    if 'Cartesian Force Constants' in obj:
        half_hessian = obj.get('Cartesian Force Constants')
        geometrical_derivatives['GG'] = derivatives_g.BaseGeometricalDerivativeTensor(
            spacial_dof=spacial_dof,
            trans_plus_rot=trans_plus_rot_dof,
            representation='GG',
            components=make_cartesian_hessian_from_fchk(half_hessian, spacial_dof))

    if not geometrical_derivatives:
        raise PropertyNotPresent('geometrical_derivatives')

    return geometrical_derivatives


@FCHK.define_property('excitations')
def gaussian__fchk__property__excitations(obj, *args, **kwargs):
    """Get excitation properties in FCHK. Returns a dictionary:

    .. code-block:: text

        + "!" : Tensor
        + "!F" : Tensor

    :param obj: object
    :type obj: qcip_tools.chemistry_files.gaussian.FCHK
    :rtype: dict
    """

    excitations = {}

    scf_energy = obj.get('SCF Energy')

    if 'ETran scalars' in obj:
        nstate = obj.get('ETran scalars')[0]
        n_per_state = obj.get('ETran scalars')[1]
        if nstate < 0:
            raise Exception('number of excitation < 0 ?')

        if 'ETran state values' in obj:
            c = obj.get('ETran state values')
            excitations['!'] = derivatives.Tensor('!', nstates=nstate + 1)
            excitations['!F'] = derivatives.Tensor('!F', nstates=nstate + 1)

            for i in range(nstate):
                excitations['!'].components[i + 1] = c[i * n_per_state] - scf_energy
                excitations['!F'].components[i + 1] = c[i * n_per_state + 1:i * n_per_state + 4]

    if not excitations:
        raise PropertyNotPresent('excitations')

    return excitations


class LinkCalled:
    """Remind when a link is called

    :param link: number of the link called
    :type link: int
    :type line_start: int
    :type line_end: int
    """

    def __init__(self, link, line_start, line_end):

        self.link = link
        self.line_start = line_start
        self.line_end = line_end

    def __repr__(self):
        return 'Link {}: {}:{}'.format(self.link, self.line_start, self.line_end)


class OutputFormatError(FormatError):
    pass


class Output(ChemistryLogFile, WithMoleculeMixin, WithIdentificationMixin):
    """Log file of Gaussian. Contains a lot of informations, but not well printed or located.
    If possible, rely on the FCHK rather than on the LOG file.

    .. note::

        ``self.geometry`` contains the input geometry, but it may be reoriented according to the symmetry if
        ``nosym`` was not provided. If it was not reoriented, ``self.input_orientation`` is set to ``True``.

    .. warning::

        Results are not guaranteed if the job ran without ``#P``.

    .. container:: class-members

        + ``self.molecule``: the molecule (``qcip_tools.molecule.Molecule``)
        + ``self.input_orientation``: wether the molecule has been reoriented w.r.t. input (``bool``)
        + ``self.chunks``: the different links called (``list`` of
          `LinkCalled <#qcip_tools.chemistry_files.gaussian.LinkCalled>`_)
        + ``self.lines``: the lines of the file (``list`` of ``str``)

    """

    #: The identifier
    file_type = 'GAUSSIAN_LOG'
    chunk_title_variable = 'link'

    def __init__(self):
        self.molecule = molecule.Molecule()
        self.input_orientation = False

    @classmethod
    def attempt_identification(cls, f):
        """A gaussian log ... Contains a lot of "gaussian" in the beginning (limit to the 150 first lines)"
        """

        count = 0
        num_of_gaussian = 0

        for line in f.readlines():
            if count > 100:
                break
            if 'Gaussian' in line or 'gaussian' in line:
                num_of_gaussian += 1
            count += 1

        return num_of_gaussian > 10

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        super().read(f)

        found_enter_link1 = False

        for index, line in enumerate(self.lines):
            if 'Entering Link 1' in line:
                self.chunks.append(LinkCalled(1, index, -1))
                found_enter_link1 = True
                break

        if not found_enter_link1:
            raise OutputFormatError('this does not seems to be a Gaussian output')

        for index, line in enumerate(self.lines):
            # List links (needs a job that ran with `#p`)
            if '(Enter ' in line:
                self.chunks[-1].line_end = index - 1

                link_num = line[-10:-6]
                if link_num[0] == 'l':
                    link_num = link_num[1:]
                current_link = int(link_num)
                self.chunks.append(LinkCalled(current_link, line_start=index, line_end=-1))

        # fetch molecule:
        charge_and_multiplicity_line = self.search('Charge = ', into=101)
        orientation_line = self.search('orientation:', into=202)

        if -1 in [charge_and_multiplicity_line, orientation_line]:
            raise OutputFormatError('Could not find geometry')

        self.input_orientation = self.lines[orientation_line][26] == 'I'

        for line in self.lines[orientation_line + 5:]:
            if '-------------' in line:
                break
            self.molecule.insert(atom.Atom(
                atomic_number=int(line[12:21].strip()),
                position=[float(a) for a in line[36:].split()]
            ))

        self.molecule.charge = int(self.lines[charge_and_multiplicity_line][9:13])
        self.molecule.multiplicity = int(self.lines[charge_and_multiplicity_line][27:])


def _find_in_links(obj, q, links):
    found = filter(lambda n: n > 0, (obj.search(q, into=link) for link in links))
    try:
        return next(found)
    except StopIteration:
        return -1


@Output.define_property('computed_energies')
def gaussian__output__get_computed_energies(obj, *args, **kwargs):
    """Get the energies... Actually, "only"
    + HF and DFT energy (in link 502/503/506/508)
    + MP2 energy (in link 804/903/905/906)
    + MP3, CCSD and CCSD(T) energy (in link 913)

    :param obj: object
    :type obj: qcip_tools.chemistry_files.gaussian.Output
    :rtype: dict
    """
    energies = {}
    modified_gaussian = False  # pristine Gaussian 16

    # fetch HF or DFT energy
    n = _find_in_links(obj, 'SCF Done:', [502, 503, 506, 508])
    if n > 0:
        chunks = obj.lines[n][12:].split()
        e = float(chunks[2])
        if (len(chunks[2]) - chunks[2].find('.')) != 9:
            modified_gaussian = True
        else:  # try to fetch more decimals above: "E = xxx" contains eleven of them
            for line in reversed(obj.lines[:n]):
                if 'Delta-E=' in line:
                    e = float(line.split()[1])
                    break
                elif 'Cycle' in line:
                    break

        if 'HF' in chunks[0]:
            energies['HF'] = e
        else:
            energies['SCF/DFT'] = e

        energies['total'] = e

    # fetch MP2 energies
    n = _find_in_links(obj, 'EUMP2 =', [804, 903, 905, 906])
    if n > 0:
        chunk = obj.lines[n].split()[-1]
        if not modified_gaussian:
            chunk = chunk.replace('D', 'E')

        energies['MP2'] = float(chunk)
        energies['total'] = energies['MP2']

    # fetch MP3 energies
    n = _find_in_links(obj, 'EUMP3=', [913])
    if n > 0:
        chunk = obj.lines[n].split()[-1]
        if not modified_gaussian:
            chunk = chunk.replace('D', 'E')
        energies['MP3'] = float(chunk)

    # fetch CCSD energies
    n = obj.search('Wavefunction amplitudes converged.', into=913)
    if n > 0:
        chunk = obj.lines[n - 3][:-16].split()[-1]
        if not modified_gaussian:
            chunk = chunk.replace('D', 'E')
        energies['CCSD'] = float(chunk)
        energies['total'] = energies['CCSD']

    # fetch CCSD(T) energies
    n = obj.search('CCSD(T)=', into=913)
    if n > 0:
        chunk = obj.lines[n].split()[-1]
        if not modified_gaussian:
            chunk = chunk.replace('D', 'E')
        energies['CCSD(T)'] = float(chunk)
        energies['total'] = energies['CCSD(T)']

    return energies


@Output.define_property('excitations')
def gaussian__output__property__excitations(obj, *args, **kwargs):
    """Get excitation properties in FCHK. Returns a dictionary:

    .. code-block:: text

        + "!" : Tensor
        + "!F" : Tensor

    :param obj: object
    :type obj: qcip_tools.chemistry_files.gaussian.Output
    :rtype: dict
    """

    excitations = {}

    for chunk in obj.chunks:
        if chunk.link == 914:

            # get transition dipoles
            s = obj.search(
                'Ground to excited state transition electric dipole moments',
                line_start=chunk.line_start,
                line_end=chunk.line_end)

            if s < 0:
                continue

            nstates = 0
            dipoles = [[.0, .0, .0]]

            for li in obj.lines[s + 2:]:
                if 'Ground' in li:
                    break
                nstates += 1
                dipoles.append(list(float(i) for i in li.split()[1:4]))

            if nstates == 0:
                continue

            excitations['!F'] = derivatives.Tensor('!F', nstates=nstates + 1, components=numpy.vstack(dipoles))

            # look for excitation energies
            s = obj.search('Excitation energies and oscillator strengths:', line_start=s, line_end=chunk.line_end)
            if s < 0:
                del excitations['!F']
                continue

            excitations['!'] = derivatives.Tensor('!', nstates=nstates + 1)
            excitations['_<S2>'] = [obj.molecule.square_of_spin_angular_moment()]
            excitations['_CSFs'] = [derivatives_exci.ConfigurationStateFunction([])]
            n = 1

            number_of_e = obj.molecule.number_of_electrons()
            homo_level = number_of_e // 2
            is_closed_shell = obj.molecule.multiplicity == 1

            for index, li in enumerate(obj.lines[s + 2:chunk.line_end]):
                if 'Excited State' in li:
                    se = li.split()
                    excitations['_<S2>'].append(float(se[-1][-5:] if se[-1][-6] == '=' else se[-1][-6:]))
                    excitations['!'].components[n] = float(se[4]) * quantities.convert(
                        quantities.ureg.eV, quantities.ureg.hartree)

                    # get CSF
                    configs = []
                    for li_e in obj.lines[s + 2 + index + 1: chunk.line_end]:
                        if '->' not in li_e:
                            break
                        fr = li_e[:9].strip()
                        to = li_e[11:20].strip()
                        coef = float(li_e[20:].strip())
                        spin = None

                        if not is_closed_shell:
                            spin = fr[-1]
                            fr = int(fr[:-1]) - homo_level - (1 if spin == 'A' else 0)
                            to = int(to[:-1]) - homo_level - 1 - (1 if spin == 'A' else 0)
                        else:
                            fr = int(fr) - homo_level
                            to = int(to) - homo_level - 1
                            coef *= math.sqrt(2)

                        configs.append((coef, derivatives_exci.Configuration(
                            [(fr, to, spin)]
                        )))

                    excitations['_CSFs'].append(derivatives_exci.ConfigurationStateFunction(configs))

                    n += 1

    if not excitations:
        raise PropertyNotPresent('excitations')

    return excitations


@Output.define_property('electrical_derivatives')
def gaussian__output__property__electrical_derivatives(obj, *args, **kwargs):
    """

    :param obj: object
    :type obj: gaussian.Output
    :rtype: dict
    """

    line_start = obj.search(' Property number 1', into=801)
    if line_start < 0:
        return {}

    offset = 0
    data = {}
    while True:
        line = obj.lines[line_start + offset]
        if 'Property' not in line:
            break
        elif 'Alpha' in line:
            frequency = float(line[-10:-2])
            tensor = derivatives_e.PolarisabilityTensor(frequency=frequency if frequency != .0 else 'static')
            for i in range(3):
                cpts = obj.lines[line_start + offset + 2 + i].split()
                tensor.components[i, :] = [float(cpts[j + 1].replace('D', 'e')) for j in range(3)]
            if frequency == .0:
                data.update({'FF': {'static': tensor}})
            else:
                if 'dD' not in data:
                    data['dD'] = {}
                data['dD'].update({frequency: tensor})
            offset += 5
        elif 'Beta' in line:
            frequency = float(line[-10:-2])
            is_SHG = 'w,w,-2w' in line
            tensor = derivatives_e.FirstHyperpolarisabilityTensor(
                input_fields=(0, 0) if frequency == .0 else ((1, 1) if is_SHG else (0, 1)),
                frequency=frequency if frequency != .0 else 'static')

            compressed_tensor = []
            for i in range(18):
                compressed_tensor.append(float(obj.lines[line_start + offset + 2 + i][-15:].replace('D', 'e')))
            if is_SHG:
                tensor.components = beta_SHG_from_fchk(compressed_tensor)
                if frequency == .0:
                    data.update({'FFF': {'static': tensor}})
                else:
                    if 'XDD' not in data:
                        data['XDD'] = {}
                    data['XDD'].update({frequency: tensor})
            else:
                tensor.components = beta_EOP_from_fchk(compressed_tensor)
                if frequency != .0:
                    if 'dDF' not in data:
                        data['dDF'] = {}
                    data['dDF'].update({frequency: tensor})

            offset += 20
        else:  # probably Gamma
            break

    return data


class ChargeTransferInformation:
    """Analyzed charge transfer.

    In it, ``charge`` is the amount of charge (in e) transfered from ``barycenter_m`` to ``barycenter_p``,
    linked by ``vector``, which lead to a certain charge transfer ``distance``.

    :param charge: charge
    :type charge: float
    :param barycenter_m: barycenter of the negative charges
    :type barycenter_m: numpy.array
    :param barycenter_p: barycenter of the positive charges
    :type barycenter_p: numpy.array
    """

    def __init__(self, charge, barycenter_m, barycenter_p):
        self.charge = charge
        self.barycenter_m = barycenter_m
        self.barycenter_p = barycenter_p
        self.vector = barycenter_m - barycenter_p
        self.distance = numpy.linalg.norm(self.vector)


class CubeFormatError(FormatError):
    pass


class Cube(ChemistryFile, WithOutputMixin, WithMoleculeMixin, WithIdentificationMixin):
    """Gaussian cube.

    Documentation on that format can be found `here <http://gaussian.com/cubegen/>`_.

    .. note::

        Increments and origin are stored as found, in bohr.

    .. container:: class-members

        + ``self.molecule``: the molecule (``qcip_tools.molecule.Molecule``)
        + ``self.cube_type``: set to ``MO`` if the cube contains MOs (``str``)
        + ``self.records_per_directions``: number of point in the grid (``list`` of ``int``)
        + ``self.origin``: origin of the grid (3x1 ``numpy.ndarray``)
        + ``self.increments``: increments between each point (3x1 ``numpy.ndarray``)
        + ``self.data_per_records``: data per point (``int``)
        + ``self.records``: the records (``numpy.ndarray``)
        + ``self.title``: the title (``str``)
        + ``self.subtitle``: the second line (``str``)
        + ``self.MOs``: the MOs contained in the file, if any (``list`` of ``int``)

    """

    #: The identifier
    file_type = 'GAUSSIAN_CUBE'

    def __init__(self):

        self.molecule = molecule.Molecule()

        self.cube_type = ''
        self.records_per_direction = [0, 0, 0]
        self.origin = numpy.zeros(3)
        self.increments = numpy.zeros(3)
        self.data_per_record = 1
        self.records = None
        self.title = ''
        self.subtitle = ''
        self.MOs = []

    @classmethod
    def possible_file_extensions(cls):
        return ['cub', 'cube']

    @classmethod
    def attempt_identification(cls, f):
        """A cube file contains, in line 3-6, quadruplets of numbers, the first one being an integer, then 3 floats.
        Then, after the different atoms, line should contains numbers in scientific notation.
        """

        count = 0
        number_of_atoms = 0

        for line in f.readlines():
            count += 1
            # line 3-6:
            if 7 > count > 2:
                c = line.split()
                if 6 < len(c) < 4:
                    return False

                if not c[0].isdigit():
                    return False

                for i in c[1:]:
                    try:
                        float(i)
                    except ValueError:
                        return False

                if count == 3:
                    number_of_atoms = abs(int(c[0]))
            # only numbers in scientific notation in records
            if count == 7 + number_of_atoms:
                c = line.split()

                for i in c:
                    if 'E' not in i:
                        return False
                return True

        return False

    @classmethod
    def from_molecule(cls, molecule, *args, **kwargs):
        raise NotImplementedError('from_molecule')

    def read(self, f):
        """

        .. warning::

            The charge/multiplicity may be wrong !

        :param f: File
        :type f: file
        """

        self.from_read = True

        lines = f.readlines()
        additional_line = False
        self.cube_type = 'unknown'
        self.data_per_record = 1

        # titles:
        self.title = lines[0].strip()
        self.subtitle = lines[1].strip()

        self.MOs = []

        # parameters
        params = lines[2].split()
        num_of_atoms = int(params[0])
        if num_of_atoms < 0:
            additional_line = True
            num_of_atoms = - num_of_atoms
            self.cube_type = 'MO'

        self.origin = numpy.array([float(a) for a in params[1:4]])

        if len(params) > 4:
            self.data_per_record = int(params[4])

        # number of records and increments
        for i in range(0, 3):
            self.increments[i] = float(lines[3 + i].split()[i + 1])

        for i in range(0, 3):
            num = int(lines[i + 3][0:8].strip())
            self.records_per_direction[i] = num

        # atoms:
        self.molecule = molecule.Molecule()
        for i in range(0, num_of_atoms):
            atm = lines[i + 6].split()
            symbol = atom.AtomicNumberToSymbol[int(atm[0])]
            atom_ = atom.Atom(
                symbol=symbol, position=numpy.array([float(a) * quantities.AuToAngstrom for a in atm[2:5]]))
            self.molecule.insert(atom_)

        # records
        start_line = 6 + num_of_atoms

        if additional_line:
            # number of stuffs per line
            extra_line = lines[start_line].strip().split()
            self.data_per_record = int(extra_line[0])

            for i in extra_line[1:]:
                self.MOs.append(int(i))

            start_line += 1

        total_num_of_data = self.number_of_data()
        actual_num_of_data = 0

        self.records = numpy.zeros(total_num_of_data)

        for line in lines[start_line:]:
            raw = [float(a) for a in line.strip().split()]
            l_size = len(raw)
            self.records[actual_num_of_data:actual_num_of_data + l_size] = raw
            actual_num_of_data += l_size

        reshape = self.records_per_direction.copy()
        reshape.append(self.data_per_record)

        self.records = self.records.reshape(tuple(reshape))

        if actual_num_of_data != total_num_of_data:
            raise CubeFormatError(
                'num of record is not equal to total num of records ({}!={})'.format(
                    actual_num_of_data, total_num_of_data))

    def number_of_records(self):
        """Number of records in the file

        .. note::

            According to the documentation:

             |  Cube files have one row per record (i.e., N1*N2 records each of length N3*NVal) (where N1 and N2
               are the record per direction in the X and Y direction, while N3 is for the Z direction, the fastest
               running index).

            So it returns the number of rows, and should not be mistaken with ``number_of_data()``.

        :rtype: int
        """

        return numpy.prod(self.records_per_direction[:2])

    def number_of_data(self):
        """Number of points in the file

        :rtype: int
        """
        return self.number_of_records() * self.records_per_direction[2] * self.data_per_record

    def dV(self):
        """Return the size of the smallest unit of volume in the cube, in angstrom**3.

        :rtype: float
        """

        return numpy.prod(self.increments) * quantities.AuToAngstrom ** 3

    def to_string(self):
        """

        :rtype: str
        """

        r = '{}\n{}\n'.format(self.title, self.subtitle)

        num_of_atoms = len(self.molecule)
        if self.cube_type == 'MO':
            num_of_atoms = - num_of_atoms

        # origin and increments
        r += '{:5d} {:11f} {:11f} {:11f} {:4d}\n'.format(
            num_of_atoms, self.origin[0], self.origin[1], self.origin[2], self.data_per_record)

        r += '{:5d} {:11f} {:11f} {:11f}\n'.format(self.records_per_direction[0], self.increments[0], 0.0, 0.0)
        r += '{:5d} {:11f} {:11f} {:11f}\n'.format(self.records_per_direction[1], 0.0, self.increments[1], 0.0)
        r += '{:5d} {:11f} {:11f} {:11f}\n'.format(self.records_per_direction[2], 0.0, 0.0, self.increments[2])

        # molecule
        for a in self.molecule:
            r += '{:5d} {:11f} {:11f} {:11f} {:11f}\n'.format(
                a.atomic_number, a.atomic_number, *(a.position / quantities.AuToAngstrom))

        if self.cube_type == 'MO':
            r += '    {}  {}\n'.format(len(self.MOs), '  '.join(str(m) for m in self.MOs))

        # records:
        for x in range(self.records_per_direction[0]):
            for y in range(self.records_per_direction[1]):
                for z in range(self.records_per_direction[2]):
                    for index_data in range(self.data_per_record):
                        if (z * self.data_per_record + index_data) % 6 == 0 and z != 0:
                            r += '\n'

                        r += ' {: .5E}'.format(self.records[x][y][z][index_data])

                r += '\n'

        return r

    def alike(self, other):
        """Compare two cubes to check if they have the same size

        :param other: the other cube
        :type other: Cube
        :rtype: bool
        """

        if self.number_of_data() != other.number_of_data():
            return False

        if not numpy.array_equal(self.increments, other.increments):
            return False

        if self.records.shape != other.records.shape:
            return False

        return True

    def __add__(self, other):
        """Sum the values of two cubes

        :param other: other cube
        :type other: Cube
        :rtype: Cube
        """

        if not self.alike(other):
            raise ArithmeticError('cannot sum with this cube')

        new_cube = copy.deepcopy(self)
        new_cube.records += other.records

        return new_cube

    def __sub__(self, other):
        """Subtract the values of two cubes

        :param other: other cube
        :type other: Cube
        :rtype: Cube
        """

        if not self.alike(other):
            raise ArithmeticError('cannot sub with this cube')

        new_cube = copy.deepcopy(self)
        new_cube.records -= other.records

        return new_cube

    def __pow__(self, power, modulo=None):
        """Power the records

        :param power: power
        :rtype power: float
        :param modulo: ?
        :rtype modulo: int
        :rtype: Cube
        """

        new_cube = copy.deepcopy(self)
        new_cube.records **= power

        return new_cube

    def slice(self, index=0):
        """Get a specific part of the data

        :param index: index of the record
        :type index: int
        :rtype: Cube
        """

        if index < 0 or index >= self.data_per_record:
            raise ValueError(index)

        new_cube = copy.deepcopy(self)
        new_shape = list(self.records.shape)
        new_shape[-1] = 1
        new_cube.records = self.records[:, :, :, index].reshape(new_shape)
        new_cube.data_per_record = 1

        if self.cube_type == 'MO':
            new_cube.MOs = [self.MOs[index]]

        return new_cube

    def positions(self):
        """Return an array with the position (in angstrom) of each point.

        :rtype: numpy.array
        """

        s = self.records_per_direction[:]
        s.append(3)
        positons = numpy.zeros(tuple(s))

        for x in range(self.records_per_direction[0]):
            for y in range(self.records_per_direction[1]):
                for z in range(self.records_per_direction[2]):
                    current_position = self.origin + numpy.array([x, y, z]) * self.increments
                    positons[x, y, z] = current_position * quantities.AuToAngstrom
        return positons

    def compute_charge_transfer(self, index=0):
        """
        Analyse CT by computing:

        - Transferred charge ;
        - Barycenter positions (the negative one being the starting point of the transfer) ;
        - Distance between the two barycenters.

        Made according to *J. Chem. Theory. Comput.* **7**, 2498 (2011).

        :param index: slice of the record to use
        :type index: int
        :rtype: ChargeTransferInformation
        """

        dV = self.dV()

        # positive and negative densities:
        rho_p = self.records.copy()
        rho_p[numpy.where(rho_p < 0.0)] = 0.0

        rho_m = self.records.copy()
        rho_m[numpy.where(rho_p > 0.0)] = 0.0

        charge_transferred = .5 * (rho_p.sum() - rho_m.sum()) * dV

        # find barycenter positions (in angstrom)
        positions = self.positions()
        bar_p = numpy.zeros(3)
        bar_m = numpy.zeros(3)

        for i in range(3):
            bar_p[i] = numpy.sum(positions[:, :, :, i] * rho_p[:, :, :, index]) / charge_transferred * dV
            bar_m[i] = - numpy.sum(positions[:, :, :, i] * rho_m[:, :, :, index]) / charge_transferred * dV

        return ChargeTransferInformation(
            charge=charge_transferred / (quantities.AuToAngstrom ** 3),  # charge in |e|
            barycenter_m=bar_m,
            barycenter_p=bar_p)

    def sum_density_of_sets(self, bounding_sets, index=0):
        """Sum the density in a given space, defined by the bounding boxes

        :param bounding_sets: sets of bounding boxes
        :type bounding_sets: dict
        :param index: slice of the record to use
        :type index: int
        :return: dictionary with the sum of the density
        :rtype: dict
        """

        if len(bounding_sets) == 0:
            return {}

        sums = {}

        for k in bounding_sets:
            sums[k] = .0

        dV = self.dV() / (quantities.AuToAngstrom ** 3)

        for x in range(self.records_per_direction[0]):
            for y in range(self.records_per_direction[1]):
                for z in range(self.records_per_direction[2]):
                    current_position = (self.origin + [x, y, z] * self.increments) * quantities.AuToAngstrom
                    current_value = self.records[x, y, z, index] * dV

                    for key, bounding_set in bounding_sets.items():
                        if current_position in bounding_set:
                            sums[key] += current_value

        return sums


class BasisSetFormatError(FormatError):
    pass


class BasisSet(ChemistryFile, WithOutputMixin):
    """File containing a basis set in the gaussian format.

    This file format does not exists, but lets say that copy/paste from ESML basis set exchange should work.

    .. container:: class-members

        + ``self.basis_set``: the basis_set (``qcip_tools.basis_set.BasisSet``)

    """

    file_type = 'GAUSSIAN_BS'

    def __init__(self):
        self.basis_set = basis_set.BasisSet()

    def __contains__(self, item):
        return item in self.basis_set

    def __getitem__(self, item):
        return self.basis_set[item]

    def read(self, f, fail_for_same_atom=False):
        """

        :param f: f
        :type f: file
        :param fail_for_same_atom: raise an exception is the same atom is defined more than once
        :type fail_for_same_atom: bool
        """

        def treat_atomic_basis_set(part):
            lines = [
                a.strip() for a in part.splitlines() if len(a) != 0 and a[0] != '!' and a.strip() not in ['', '****']]

            if len(lines) == 0:
                return  # nothing to search here

            if len(lines) < 3:
                raise BasisSetFormatError('atomic basis set too short for {}'.format(lines))

            atom_def = lines[0].split()

            if len(atom_def) != 2:
                raise BasisSetFormatError('wrong atom definition {}'.format(lines[0]))

            try:
                if atom_def[0] in self.basis_set:
                    if fail_for_same_atom:
                        raise BasisSetFormatError('two times {}'.format(atom_def[0]))
                    else:
                        return
            except KeyError as e:
                raise BasisSetFormatError('not a atomic symbol: {}'.format(e))

            a = basis_set.AtomicBasisSet(atom_def[0])
            next_basis_function = 1

            while True:
                if next_basis_function >= len(lines):
                    break

                def_basis_function = lines[next_basis_function].split()
                if len(def_basis_function) != 3:
                    raise BasisSetFormatError(
                        'wrong definition of shell {} for {}'.format(lines[next_basis_function], atom))

                shell = def_basis_function[0]
                number_of_primitives = int(def_basis_function[1])
                normalization = float(def_basis_function[2])

                if shell not in basis_set.SHELL_TO_TAM:
                    raise BasisSetFormatError('unknown shell {} for {}'.format(shell, atom))

                if next_basis_function + number_of_primitives >= len(lines):
                    raise BasisSetFormatError('expected {} primitives for shell {} of {}'.format(
                        number_of_primitives, shell, atom))

                f = basis_set.Function(normalization=normalization)

                for i in range(number_of_primitives):
                    def_primitive = lines[next_basis_function + i + 1].split()
                    if len(def_primitive) not in [2, 3]:
                        raise BasisSetFormatError('wrong primitive {} in shell {} of {}'.format(
                            def_primitive, shell, atom))

                    if shell == 'SP' and len(def_primitive) == 2:
                        raise BasisSetFormatError('for SP shell, 3 number are needed (in {})'.format(atom))
                    elif shell != 'SP' and len(def_primitive) == 3:
                        raise BasisSetFormatError('for non-SP shell {}, only 2 number are needed (in {})'.format(
                            shell, atom))

                    f.add_primitive(basis_set.Primitive(*[float(a) for a in def_primitive]))

                next_basis_function += number_of_primitives + 1
                a.add_basis_function(shell, f)

            self.basis_set.add_atomic_basis_set(a)

        content = f.read()
        prev_pos = 0
        self.from_read = True

        while True:
            next_separator = content.find('****', prev_pos + 4)

            if next_separator < 0:
                break

            treat_atomic_basis_set(content[prev_pos:next_separator])
            prev_pos = next_separator

        if len(self.basis_set) == 0:
            raise BasisSetFormatError('empty basis set')

    def to_string(self, for_molecule=None):
        """Give basis set

        :param for_molecule: if defined, only output for the atoms contained in the molecule
        :type for_molecule: qcip_tools.molecule.Molecule
        :rtype: str
        """

        to_output = self.basis_set.atomic_basis_sets.keys()
        if for_molecule:
            to_output = for_molecule.symbols_contained

        r = ''

        for i in to_output:
            if i not in self.basis_set:
                raise ValueError('No basis set for {} defined'.format(i))

            atomic_basis_set = self.basis_set[i]
            atom = atomic_basis_set.atom.symbol

            r += '****\n{}     0\n'.format(atom)
            for shell in atomic_basis_set.shells():
                for basis_function in atomic_basis_set[shell]:
                    r += '{:3} {:<3} {:.2f}\n'.format(shell, len(basis_function), basis_function.normalization)
                    for primitive in basis_function.primitives:
                        r += ' {:15.7f}            {: .8f}{}\n'.format(
                            primitive.exponent,
                            primitive.contraction_coefficient,
                            '            {: .8f}'.format(primitive.p_coefficient) if shell == 'SP' else '')

        return r + '****'

    def write(self, f, for_molecule=None):
        f.write(self.to_string(for_molecule=for_molecule))
